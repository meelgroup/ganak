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

using CMSat::StreamBuffer;
using CMSat::DimacsParser;
using CMSat::SATSolver;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

namespace po = boost::program_options;
using std::string;
using std::vector;
po::options_description main_options = po::options_description("Main options");
po::options_description probe_options = po::options_description("Probe options");
po::options_description lookahead_options = po::options_description("Lookeahead options");
po::options_description restart_options = po::options_description("Restart options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

using namespace std;
CMSat::SATSolver* sat_solver;
uint32_t must_mult_exp2 = 0;
bool indep_support_given = false;
set<uint32_t> indep_support;
int do_check = 0;
int exact = 1;
MTRand mtrand;
CounterConfiguration conf;
int do_hyperbin = 1;
int red_cls_also = 0;

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
    my_delta << std::setprecision(8) << conf.delta;

    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&conf.verb)->default_value(conf.verb), "verb")
    ("seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")
    ("delta", po::value(&conf.delta)->default_value(conf.delta, my_delta.str()), "Delta")
    ("singlebump", po::value(&conf.do_single_bump)->default_value(conf.do_single_bump), "Do single bumping, no double (or triple, etc) bumping of activities. Non-single bump is old ganak")
    ("branch", po::value(&conf.branch_type)->default_value(conf.branch_type), "Branching type. 0 == default, 1 == gpmc (with TD stuff)")

    ("cscore", po::value(&conf.do_cache_score)->default_value(conf.do_cache_score), "Do cache scores")
    ("hyper", po::value(&do_hyperbin)->default_value(do_hyperbin), "Do hyperbinary resolution via intree")
    ("maxcache", po::value(&conf.maximum_cache_size_bytes_)->default_value(conf.maximum_cache_size_bytes_), "Max cache size in BYTES. 0 == use 80% of free mem")
    ("exp", po::value(&conf.exp)->default_value(conf.exp), "Probabilistic Component Caching")
    ("version", "Print version info")
    ("check", po::value(&do_check)->default_value(do_check), "Check count at every step")
    ("red", po::value(&red_cls_also)->default_value(red_cls_also), "Also add redundant clauses from CNF")
    ("alluipincact", po::value(&conf.alluip_inc_act)->default_value(conf.alluip_inc_act), "All UIP should increase activities")
    ;

    restart_options.add_options()
    ("rstfirst", po::value(&conf.first_restart)->default_value(conf.first_restart), "Run restarts")

    ("restart", po::value(&conf.do_restart)->default_value(conf.do_restart), "Run restarts")
    ("rsttype", po::value(&conf.restart_type)->default_value(conf.restart_type), "Check count at every step")
    ("rstcutoff", po::value(&conf.restart_cutoff_mult)->default_value(conf.restart_cutoff_mult), "Multiply cutoff with this")

    ("onpathprint", po::value(&conf.do_on_path_print)->default_value(conf.do_on_path_print), "Print ON-PATH during restart")
    ("exact", po::value(&exact)->default_value(exact), "Exact counting")
    ;

    lookahead_options.add_options()
    ("lookahead", po::value(&conf.do_lookahead)->default_value(conf.do_lookahead), "Do lookahead?")
    ("lookaheaddepth", po::value(&conf.lookahead_depth)->default_value(conf.lookahead_depth), "Lookahead depth")
    ("looknum", po::value(&conf.lookahead_num)->default_value(conf.lookahead_num), "How many to check for lookahead")
    ;

    probe_options.add_options()
    ("probe", po::value(&conf.failed_lit_probe_type)->default_value(conf.failed_lit_probe_type), "Failed Lit Probe Type. 0 == none, 1 == full, 2 == only bottom RATIO, where ratio is given by --probeonlyafter")
    ("probeonlyafter", po::value(&conf.probe_only_after_ratio)->default_value(conf.probe_only_after_ratio), "What ratio of failed lit probe in terms of decision. Only active if '--failed 2'")
    ("probemulti", po::value(&conf.num_probe_multi)->default_value(conf.num_probe_multi), "Multiply by this amount how many variables to probe.")
    ("bprop", po::value(&conf.bprop)->default_value(conf.bprop), "Do bothprop")
    ;

    help_options.add(main_options);
    help_options.add(probe_options);
    help_options.add(lookahead_options);
    help_options.add(restart_options);
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
            << "Approx/Exact counter" << endl;

            cout
            << "Usage: ./ganak [options] inputfile" << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            cout << ganak_version_info() << endl;
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
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, NULL, conf.verb);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, NULL, conf.verb);
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

bool get_target(vector<CMSat::lbool>& model) {
  vector<double> act;
  sat_solver->set_polarity_mode(CMSat::PolarityMode::polarmode_rnd);
  /* sat_solver->set_up_for_sample_counter(100); */
  /* sat_solver->set_activities(act); */
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

void create_from_sat_solver(Counter& counter, SATSolver& ss) {
  counter.new_vars(ss.nVars());

  // Clean the clauses before we add them
  vector<CMSat::Lit> assumps;
  for(const auto& v: indep_support) assumps.push_back(CMSat::Lit(v-1, false));
  string s ("clean-cls");
  ss.simplify(&assumps, &s);

  // Irred cls
  ss.start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);
  vector<CMSat::Lit> cms_cl;
  while(ss.get_next_small_clause(cms_cl)) {
    const auto cl = cms_to_ganak_cl(cms_cl);
    counter.add_irred_cl(cl);
  }
  ss.end_getting_small_clauses();
  counter.end_irred_cls();

  if (red_cls_also) {
    uint32_t red_cl = 0;
    ss.start_getting_small_clauses(
        std::numeric_limits<uint32_t>::max(),
        std::numeric_limits<uint32_t>::max(),
        true);
    while(ss.get_next_small_clause(cms_cl)) {
      const auto cl = cms_to_ganak_cl(cms_cl);
      counter.add_red_cl(cl);
      red_cl++;
    }
    ss.end_getting_small_clauses();
    cout << "c red cl added: " << red_cl << endl;
  }
}

mpz_class check_count_independently_no_restart(const vector<vector<CMSat::Lit>>& cubes) {
  SATSolver sat_solver2;
  sat_solver2.new_vars(sat_solver->nVars());
  sat_solver->start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);

  vector<CMSat::Lit> cl;
  while(sat_solver->get_next_small_clause(cl)) sat_solver2.add_clause(cl);
  sat_solver->end_getting_small_clauses();

  if (cubes.size() > 1) {
    for(uint32_t i = 0; i < cubes.size()-1; i++) {
      sat_solver2.add_clause(cubes[i]);
    }
  }
  for(const auto& l: cubes.back()) {
    cl.clear();
    cl.push_back(~l);
    sat_solver2.add_clause(cl);
  }
  const auto res = sat_solver2.solve();
  if (res == CMSat::l_False) {
    cout << "Check count got UNSAT" << endl;
    return 0;
  }

  CounterConfiguration conf2;
  conf2.verb = 0;
  conf2.do_restart = false;
  Counter counter(conf2);

  create_from_sat_solver(counter, sat_solver2);
  counter.set_indep_support(indep_support);

  vector<Lit> largest_cube;
  const auto count = counter.count(largest_cube);
  assert(largest_cube.empty());
  return count;
}

void add_hyperbins()
{
  string s = "intree-probe";
  vector<CMSat::Lit> dont_elim;
  for(const auto& v: indep_support) dont_elim.push_back(CMSat::Lit(v-1, false));
  sat_solver->simplify(&dont_elim, &s);
}

int main(int argc, char *argv[])
{
  const double start_time = cpuTime();
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
  sat_solver = new SATSolver;
  sat_solver->set_renumber(false);
  parse_supported_options(argc, argv);
  if (conf.verb) {
    cout << ganak_version_info() << endl;
    cout << "c called with: " << command_line << endl;
  }
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
  if (!sat_solver->okay() || sat_solver->solve() == CMSat::l_False) {
    if (indep_support_given) cout << "s pmc ";
    else cout << "s mc ";
    cout << "0" << endl;
    exit(0);
  }
  if (do_hyperbin) add_hyperbins();
  mpz_class count = 0;

  vector<double> act;
  vector<uint32_t> comp_acts;
  vector<uint8_t> polars;
  double act_inc;
  uint32_t num_cubes = 0;
  // TODO: add hyper-binary BIN clauses to GANAK
  conf.next_restart = conf.first_restart;
  Counter* counter = new Counter(conf);
  create_from_sat_solver(*counter, *sat_solver);
  counter->set_indep_support(indep_support);
  counter->init_activity_scores();
  vector<vector<CMSat::Lit>> cubes;
  mpz_class total_check_count = 0;
  double total_check_time = 0;
  act.resize(2*(sat_solver->nVars()+1), 0);
  vector<CMSat::lbool> model;
  while (sat_solver->okay()) {
    double call_time = cpuTime();
    if (!get_target(model)) break;
    counter->set_target_polar(model);
    vector<Lit> largest_cube;
    mpz_class this_count = counter->count(largest_cube);
    count += this_count;
    auto cms_cl = ganak_to_cms_cl(largest_cube);
    cubes.push_back(cms_cl);
    cout << "c cnt for this cube: " << std::setw(15) << std::left << this_count
      << " cube sz: " << std::setw(6) << cms_cl.size()
      << " cube num: " << std::setw(3) << num_cubes
      << " cnt so far: " << std::setw(15) << count
      << " T: " << std::setprecision(2) << std::fixed << (cpuTime() - call_time)
      << endl;
    if (conf.verb >= 2) {
      cout << "c ---> cube: ";
      for(const auto& l: cms_cl) cout << l << " ";
      cout << "0" << endl;
    }
    cout << "c Total time until now: "
      << std::fixed << (cpuTime() - start_time) - total_check_time<< endl;

    if (do_check) {
      double this_check_time = cpuTime();
      auto check_count = check_count_independently_no_restart(cubes);
      total_check_count += check_count;
      if (check_count != this_count && conf.verb >= 2) {
        cout << "Check count says: " << check_count << endl;
      }
      if (exact) release_assert(check_count == this_count);
      else {
        if (conf.verb >=2)
          cout << "Difference rel this cube: " << this_count.get_d()/check_count.get_d()  << endl;
        cout << "Difference rel: " << count.get_d()/total_check_count.get_d() << endl;
      }
      total_check_time += cpuTime() - this_check_time;
    }
    /* cout << "Before miniim: " << cms_cl << endl; */
    /* if (!sat_solver->minimize_clause(cms_cl)) { */
    /*   cout << "After minim: " << cms_cl << endl; */
      sat_solver->add_clause(cms_cl);
    /* } */
    if (exact) counter->get_activities(act, polars, act_inc, comp_acts);
    sat_solver->set_verbosity(0);
    num_cubes++;
    conf.next_restart*=2;
    if (conf.next_restart > 20*conf.first_restart) conf.next_restart = conf.first_restart;
    counter->set_next_restart(conf.next_restart);

    //Exact
    if (exact && sat_solver->okay()) {
      delete counter;
      counter = new Counter(conf);
      create_from_sat_solver(*counter, *sat_solver);
      counter->set_indep_support(indep_support);
      counter->set_activities(act, polars, act_inc, comp_acts);
    }
  }
  mpz_mul_2exp(count.get_mpz_t(), count.get_mpz_t(), must_mult_exp2);
  cout << "c Time: " << std::setprecision(2) << std::fixed << (cpuTime() - start_time) - total_check_time << endl;
  if (indep_support_given) cout << "s pmc ";
  else cout << "s mc ";
  cout << count << endl;

  delete counter;
  delete sat_solver;
  return 0;
}
