/******************************************
Copyright (c) 2023, Marc Thurley, Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#include "cryptominisat5/cryptominisat.h"
#include "solver.h"
#include "GitSHA1.h"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include "src/GitSHA1.h"

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

namespace po = boost::program_options;
using std::string;
using std::vector;
po::options_description main_options = po::options_description("Main options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

using namespace std;
int verb = 1;
int seed = 0;
int do_comp_caching = 1;
uint64_t max_cache = 0;
int do_implicit_bcp = 1;
int do_restart = 1;
int do_pcc = 1;
int hashrange = 1;
double delta = 0.05;
uint32_t first_restart = 100000000;

string ganak_version_info()
{
    std::stringstream ss;
    ss << "c GANAK SHA revision " << GANAK::get_version_sha1() << endl;
    ss << "c GANAK compilation env " << GANAK::get_compilation_env() << endl;
    #ifdef __GNUC__
    ss << "c GANAK compiled with gcc version " << __VERSION__ << endl;
    #else
    ss << "c GANAK compiled with non-gcc compiler" << endl;
    #endif
    CMSat::SATSolver sat_solver;
    ss << "c CMS version: " << sat_solver.get_version_sha1();

    return ss.str();
}

void add_appmc_options()
{

    std::ostringstream my_delta;
    my_delta << std::setprecision(8) << delta;

    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&verb)->default_value(1), "verb")
    ("seed,s", po::value(&seed)->default_value(seed), "Seed")
    ("delta", po::value(&delta)->default_value(delta, my_delta.str()), "Delta")
    ("rstfirst", po::value(&first_restart)->default_value(first_restart), "Run restarts")
    ("restart", po::value(&do_restart)->default_value(do_restart), "Run restarts")
    ("cc", po::value(&do_comp_caching)->default_value(do_comp_caching), "Component caching")
    ("maxcache", po::value(&max_cache)->default_value(max_cache), "Max cache size in MB. 0 == use 80% of free mem")
    ("ibpc", po::value(&do_implicit_bcp)->default_value(do_implicit_bcp), "Implicit Boolean Constraint Prop")
    ("version", "Print version info")
    ("pcc", po::value(&do_pcc)->default_value(do_pcc), "Probabilistic Component Caching")
    ("hashrange", po::value(&hashrange)->default_value(hashrange), "Seed")
    ;

    help_options.add(main_options);
}

void add_supported_options(int argc, char** argv)
{
    add_appmc_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Exact counter" << endl;

            cout
            << "Usage: ./ganak [options] inputfile" << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            cout << ganak_version_info();
            std::exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unknown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only the input CNF can be given as a positional option." << endl;
        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }

}

int main(int argc, char *argv[])
{
  #if defined(__GNUC__) && defined(__linux__)
  feenableexcept(FE_INVALID   |
                 FE_DIVBYZERO |
                 FE_OVERFLOW
                );
  #endif

  //Reconstruct the command line so we can emit it later if needed
  string command_line;
  for(int i = 0; i < argc; i++) {
      command_line += string(argv[i]);
      if (i+1 < argc) {
          command_line += " ";
      }
  }
  Solver solver;
  add_supported_options(argc, argv);

#ifndef DOPCC
  solver.config().perform_pcc = false;
#endif

  solver.config().do_comp_caching = do_comp_caching;
  solver.config().do_failed_lit_probe = do_implicit_bcp;
  solver.config().do_restart = do_restart;
  solver.config().verb = verb;
  solver.config().do_pcc = do_pcc;
  solver.config().randomseed = seed;
  solver.config().hashrange = hashrange;
  solver.config().delta = delta;
  solver.config().first_restart = first_restart;
  solver.statistics().maximum_cache_size_bytes_ = max_cache * 1024ULL*1024ULL;

  if (verb) {
    cout << ganak_version_info() << endl;
    cout << "c called with: " << command_line << endl;
  }

  if (vm.count("input") != 0) {
    vector<string> inp = vm["input"].as<vector<string> >();
    if (inp.size() > 1) {
        cout << "[appmc] ERROR: you must only give one CNF as input" << endl;
        exit(-1);
    }
    solver.solve(inp[0]);
  } else {
    // TODO read stdin, once we are a library.
    cout << "ERROR: must give input file to read" << endl;
    exit(-1);
  }
  return 0;
}
