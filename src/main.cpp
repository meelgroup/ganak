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
#include "cryptominisat5/solvertypesmini.h"
#include "solver.h"
#include "GitSHA1.h"

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <iomanip>
#include <time_mem.h>
#include <boost/program_options.hpp>
#include "src/GitSHA1.h"
#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <arjun/arjun.h>

using CMSat::StreamBuffer;
using CMSat::DimacsParser;
using CMSat::SATSolver;
ArjunNS::Arjun* arjun = NULL;

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
int verb = 0;
int seed = 0;
int do_comp_caching = 1;
uint64_t max_cache = 0;
int do_implicit_bcp = 1;
int do_restart = 1;
int do_pcc = 1;
int do_arjun = 1;
int hashrange = 1;
double delta = 0.05;
uint64_t first_restart_start = 40;
uint64_t first_restart;
CMSat::SATSolver* sat_solver;
uint32_t must_mult_exp2 = 0;
bool indep_support_given = false;
set<uint32_t> indep_support;
int do_check = 0;
MTRand mtrand;

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
    ss << "c CMS version: " << sat_solver->get_version_sha1();

    return ss.str();
}

void add_ganak_options()
{
    std::ostringstream my_delta;
    my_delta << std::setprecision(8) << delta;

    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&verb)->default_value(verb), "verb")
    ("seed,s", po::value(&seed)->default_value(seed), "Seed")
    ("delta", po::value(&delta)->default_value(delta, my_delta.str()), "Delta")
    ("rstfirst", po::value(&first_restart_start)->default_value(first_restart_start), "Run restarts")
    ("restart", po::value(&do_restart)->default_value(do_restart), "Run restarts")
    ("cc", po::value(&do_comp_caching)->default_value(do_comp_caching), "Component caching")
    ("maxcache", po::value(&max_cache)->default_value(max_cache), "Max cache size in MB. 0 == use 80% of free mem")
    ("ibpc", po::value(&do_implicit_bcp)->default_value(do_implicit_bcp), "Implicit Boolean Constraint Prop")
    ("version", "Print version info")
    ("pcc", po::value(&do_pcc)->default_value(do_pcc), "Probabilistic Component Caching")
    ("check", po::value(&do_check)->default_value(do_check), "Check count at every step")
    ("hashrange", po::value(&hashrange)->default_value(hashrange), "Seed")
    ("arjun", po::value(&do_arjun)->default_value(do_arjun)
        , "Use arjun to minimize sampling set")
    ;

    help_options.add(main_options);
}

void parse_supported_options(int argc, char** argv)
{
    add_ganak_options();
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

template<class T>
void parse_file(const std::string& filename, T* reader) {
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, NULL, verb);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, NULL, verb);
  #endif
  if (in == NULL) {
      std::cout << "ERROR! Could not open file '" << filename
      << "' for reading: " << strerror(errno) << endl;
      std::exit(-1);
  }
  if (!parser.parse_DIMACS(in, true)) exit(-1);
  #ifndef USE_ZLIB
  fclose(in);
  #else
  gzclose(in);
  #endif

  indep_support_given = parser.sampling_vars_found;
  if (parser.sampling_vars_found) {
    for(const auto& lit: parser.sampling_vars) indep_support.insert(lit+1);
  } else {
    for(uint32_t i = 1; i < sat_solver->nVars()+1; i++) indep_support.insert(i);
  }
  must_mult_exp2 = parser.must_mult_exp2;
}


void set_up_solver(Solver& solver) {
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
  solver.config().maximum_cache_size_bytes_ = max_cache * 1024ULL*1024ULL;
}

vector<Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<Lit> ganak_cl;
  for(const auto& l: cl) ganak_cl.push_back(Lit(l.var()+1, !l.sign()));
  return ganak_cl;
}

vector<CMSat::Lit> ganak_to_cms_cl(const vector<Lit>& cl) {
  vector<CMSat::Lit> cms_cl;
  for(const auto& l: cl) cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

vector<CMSat::Lit> ganak_to_cms_cl(const Lit& l) {
  vector<CMSat::Lit> cms_cl;
  cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

bool take_solution(vector<CMSat::lbool>& model) {
  /* sat_solver->set_polarity_mode(CMSat::PolarityMode::polarmode_rnd); */
  //solver.set_up_for_sample_counter(100);
  CMSat::lbool ret = sat_solver->solve();
  assert(ret != CMSat::l_Undef);
  if (ret == CMSat::l_False) {
    return false;
  }
  model = sat_solver->get_model();
#if 0
  cout <<"c Model: ";
  for(int i = 0; i < model.size(); i ++)
    cout << (i+1) * (model[i] == CMSat::l_True ? 1 : -1) << " ";
  cout << "0" << endl;;
#endif
  return true;
}

void create_from_sat_solver(Solver& solver, SATSolver& ss) {
  solver.new_vars(sat_solver->nVars());
  ss.start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);
  vector<CMSat::Lit> cms_cl;
  while(ss.get_next_small_clause(cms_cl)) {
    const auto cl = cms_to_ganak_cl(cms_cl);
    solver.add_irred_cl(cl);
  }
  ss.end_getting_small_clauses();
  solver.end_irred_cls();

  uint32_t num_bins = 0;
  ss.start_getting_small_clauses(
      2,
      std::numeric_limits<uint32_t>::max(),
      true);
  while(ss.get_next_small_clause(cms_cl)) {
    const auto cl = cms_to_ganak_cl(cms_cl);
    if (cl.size() == 2) {
      solver.add_red_cl(cl);
      num_bins++;
    }
  }
  ss.end_getting_small_clauses();
  cout << "Num bins from CMS: " << num_bins << endl;
}

mpz_class check_count_independently_no_restart(const vector<CMSat::Lit>& cube) {
  SATSolver sat_solver2;
  sat_solver2.new_vars(sat_solver->nVars());
  sat_solver->start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);

  vector<CMSat::Lit> cl;
  while(sat_solver->get_next_small_clause(cl)) sat_solver2.add_clause(cl);
  sat_solver->end_getting_small_clauses();
  for(const auto& l: cube) {
    cl.clear();
    cl.push_back(~l);
    sat_solver2.add_clause(cl);
  }
  const auto res = sat_solver2.solve();
  if (res == CMSat::l_False) {
    cout << "Check count got UNSAT" << endl;
    return 0;
  }

  Solver solver;
  set_up_solver(solver);
  solver.config().do_restart = false;

  create_from_sat_solver(solver, sat_solver2);
  solver.set_indep_support(indep_support);

  vector<Lit> largest_cube;
  const auto count = solver.solve(largest_cube);
  assert(largest_cube.empty());
  return count;
}


void transfer_bins(Solver& solver, const vector<Lit>& bins)
{
  vector<Lit> cl;
  size_t at = 0;
  while(true) {
    if (at >= bins.size()) break;
    const auto& l = bins[at];
    if (l == SENTINEL_LIT) {
      assert(cl.size() == 2);
      solver.add_red_cl(cl);
      cl.clear();
    } else {
      cl.push_back(l);
    }
    at++;
  }
  cout << "c Transferred " << bins.size()/3 << " bins" << endl;
}

int main(int argc, char *argv[])
{
  double myTime = cpuTime();
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
  if (verb) {
    cout << ganak_version_info() << endl;
    cout << "c called with: " << command_line << endl;
  }
  sat_solver = new SATSolver;
  parse_supported_options(argc, argv);
  string fname;
  if (vm.count("input") != 0) {
    vector<string> inp = vm["input"].as<vector<string> >();
    if (inp.size() > 1) {
        cout << "[appmc] ERROR: you must only give one CNF as input" << endl;
        exit(-1);
    }
    fname = inp[0];
  } else {
    // TODO read stdin, once we are a library.
    cout << "ERROR: must give input file to read" << endl;
    exit(-1);
  }
  parse_file(fname, sat_solver);
  mpz_class count = 0;

  vector<double> act;
  vector<uint8_t> polars;
  double act_inc;
  uint32_t num_cubes = 0;
  vector<Lit> units;
  vector<Lit> bins;
  first_restart = first_restart_start;
  // TODO: add hyper-binary BIN clauses to GANAK
  while (sat_solver->okay()) {
    double call_time = cpuTime();
    Solver solver;
    set_up_solver(solver);
    create_from_sat_solver(solver, *sat_solver);
    if (num_cubes == 0) solver.init_activity_scores();
    solver.set_indep_support(indep_support);

    vector<CMSat::lbool> model;
    if (!take_solution(model)) break;
    solver.set_target_polar(model);
    transfer_bins(solver, bins);
    vector<Lit> largest_cube;
    if (!act.empty()) solver.set_activities(act, polars, act_inc);
    mpz_class this_count = solver.solve(largest_cube);
    solver.get_activities(act, polars, act_inc);
    units.clear();
    solver.get_unit_cls(units);
    for(const auto& l: units) sat_solver->add_clause(ganak_to_cms_cl(l));
    cout << "Transferred " << units.size() << " units -- out of vars: " << sat_solver->nVars() << endl;
    bins.clear();
    solver.get_bin_red_cls(bins);
    count += this_count;
    const auto cms_cl = ganak_to_cms_cl(largest_cube);
    cout << "c cnt for this cube: " << std::setw(15) << std::left << this_count
      << " cube sz: " << std::setw(6) << cms_cl.size()
      << " cube num: " << std::setw(3) << num_cubes
      << " cnt so far: " << std::setw(15) << count
      << " T: " << std::setprecision(2) << std::fixed << (cpuTime() - call_time)
      << endl;
    cout << "c ---> cube: ";
    for(const auto& l: cms_cl) cout << l << " ";
    cout << "0" << endl;

    if (do_check) {
      auto check_count = check_count_independently_no_restart(cms_cl);
      if (check_count != this_count) {
        cout << "Check count says: " << check_count << " ooops." << endl << endl;
      }
      assert(check_count == this_count);
    }
    sat_solver->add_clause(cms_cl);
    num_cubes++;
    first_restart*=2;
    if (first_restart > 50*first_restart_start) first_restart = first_restart_start;
  }
  mpz_mul_2exp(count.get_mpz_t(), count.get_mpz_t(), must_mult_exp2);
  cout << "c Time: " << std::setprecision(2) << std::fixed << (cpuTime() - myTime) << endl;
  if (indep_support_given) cout << "s pmc ";
  else cout << "s mc ";
  cout << count << endl;

  delete sat_solver;
  return 0;
}
