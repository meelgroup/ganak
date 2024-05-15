/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

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

#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <cryptominisat5/solvertypesmini.h>

#include "counter.hpp"
#include "GitSHA1.hpp"
#include "common.hpp"
#include "time_mem.hpp"

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <iomanip>
#include <boost/program_options.hpp>
#include "src/GitSHA1.hpp"
#include "breakid.hpp"
#include <arjun/arjun.h>

using CMSat::StreamBuffer;
using CMSat::DimacsParser;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

namespace po = boost::program_options;
using std::string;
using std::vector;
po::options_description main_options = po::options_description("Main options");
po::options_description restart_options = po::options_description("Restart options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

using namespace std;
bool indep_support_given = false;
CounterConfiguration conf;
int arjun_verb = 2;
int do_arjun = 1;
int arjun_gates = 1;
int ignore_indep = 0;
int sbva_steps = 1000;
int sbva_cls_cutoff = 4;
int sbva_lits_cutoff = 5;
int sbva_tiebreak = 1;
int do_bce = 1;
int do_breakid = 0;
ArjunNS::SimpConf simp_conf;

string ganak_version_info()
{
    std::stringstream ss;
    ss << "c o GANAK SHA revision " << GANAK::get_version_sha1() << endl;
    ss << "c o GANAK compilation env " << GANAK::get_compilation_env() << endl;
    #ifdef __GNUC__
    ss << "c o GANAK compiled with gcc version " << __VERSION__;
    #else
    ss << "c o GANAK compiled with non-gcc compiler";
    #endif
    cout << "c o CMS revision: " << CMSat::SATSolver::get_version_sha1() << endl;
    cout << "c o Arjun SHA revision: " << ArjunNS::Arjun ::get_version_info() << endl;
    cout << "c o Arjun SBVA SHA revision: " << ArjunNS::Arjun::get_sbva_version_info() << endl;

    return ss.str();
}

void add_ganak_options()
{
    // see out-ganak-6749880.pbs101-10
    simp_conf.bve_too_large_resolvent = 120;

    std::ostringstream my_delta;
    my_delta << std::setprecision(8) << conf.delta;

    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&conf.verb)->default_value(conf.verb), "verb")
    ("seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")
    ("delta", po::value(&conf.delta)->default_value(conf.delta, my_delta.str()), "Delta")
    ("arjun", po::value(&do_arjun)->default_value(do_arjun), "Use arjun")
    ("arjunverb", po::value(&arjun_verb)->default_value(arjun_verb), "Arjun verb")
    ("arjungates", po::value(&arjun_gates)->default_value(arjun_gates), "Use arjun's gate detection")
    ("ignore", po::value(&ignore_indep)->default_value(ignore_indep), "Ignore indep support given")
    ("extraclbump", po::value(&conf.do_extra_cl_bump)->default_value(conf.do_extra_cl_bump), "Also bump clauses when they propagate. By bump, we mean: set 'used' flag, and update LBD")
    ("td", po::value(&conf.do_td)->default_value(conf.do_td), "Run TD decompose")
    ("tdmaxw", po::value(&conf.td_maxweight)->default_value(conf.td_maxweight), "TD max weight")
    ("tdminw", po::value(&conf.td_minweight)->default_value(conf.td_minweight), "TD min weight")
    ("tddiv", po::value(&conf.td_divider)->default_value(conf.td_divider), "TD divider")
    ("tdweighted", po::value(&conf.do_td_weight)->default_value(conf.do_td_weight), "TD weight enabled")
    ("tdexpmult", po::value(&conf.td_exp_mult)->default_value(conf.td_exp_mult), "TD exponential multiplier")
    ("tdcheckagainstind", po::value(&conf.do_check_td_vs_ind)->default_value(conf.do_check_td_vs_ind), "Check TD against indep size")
    ("tditers", po::value(&conf.td_iters)->default_value(conf.td_iters), "TD flowcutter iterations (restarts)")
    ("tdsteps", po::value(&conf.td_steps)->default_value(conf.td_steps), "TD flowcutter number of steps at most")

    ("rdbclstarget", po::value(&conf.rdb_cls_target)->default_value(conf.rdb_cls_target), "RDB clauses target size (added to this are LBD 3 or lower)")
    ("rdbeveryn", po::value(&conf.reduce_db_everyN)->default_value(conf.reduce_db_everyN), "Reduce the clause DB every N conflicts")
    ("rdbkeepused", po::value(&conf.rdb_keep_used)->default_value(conf.rdb_keep_used), "RDB keeps clauses that are used")
    ("consolidateeveryn", po::value(&conf.consolidate_every_n)->default_value(conf.consolidate_every_n), "Consolidate every N learnt clause")
    ("lbd", po::value(&conf.base_lbd_cutoff)->default_value(conf.base_lbd_cutoff), "Initial LBD cutoff")

    ("bce", po::value(&do_bce)->default_value(do_bce), "Do BCE")
    ("cache", po::value(&conf.do_use_cache)->default_value(conf.do_use_cache), "Use (i.e. store and retrieve) cache")
    ("maxcache", po::value(&conf.maximum_cache_size_MB)->default_value(conf.maximum_cache_size_MB), "Max cache size in MB. 0 == use 80% of free mem")
    ("actexp", po::value(&conf.act_exp)->default_value(conf.act_exp), "Probabilistic Comp Caching")
    ("version", "Print version info")
    ("alluipincact", po::value(&conf.alluip_inc_act)->default_value(conf.alluip_inc_act), "All UIP should increase activities")
    ("polar", po::value(&conf.polar_type)->default_value(conf.polar_type),
     "Use polarity cache. Otherwise, false default polar.")
    ("vivif", po::value(&conf.do_vivify)->default_value(conf.do_vivify), "Vivify clauses")
    ("vivifevery", po::value(&conf.vivif_every)->default_value(conf.vivif_every), "Vivify every N conflicts")
    ("vivifmult", po::value(&conf.vivif_mult)->default_value(conf.vivif_mult), "How much to multiply timeout for vivif")
    ("vivifoutern", po::value(&conf.vivif_outer_every_n)->default_value(conf.vivif_outer_every_n), "How many restarts between outer vivif")

    ("buddy", po::value(&conf.do_buddy)->default_value(conf.do_buddy), "Run BuDDy")
    ("decide", po::value(&conf.decide)->default_value(conf.decide), "Decision type. 0 = sstd-inspired, 1 = gpmc-inspired")
    ("cachetime", po::value(&conf.cache_time_update)->default_value(conf.cache_time_update), "Cache score update type. 0 = standard, 1 = set tp new point, 2 = set to mid-point")
    ("cacherevsort", po::value(&conf.do_cache_reverse_sort)->default_value(conf.do_cache_reverse_sort), "Cache score reverse sort")
    ("sbva", po::value(&sbva_steps)->default_value(sbva_steps), "SBVA steps. 0 = no SBVA")
    ("sbvaclcut", po::value(&sbva_cls_cutoff)->default_value(sbva_cls_cutoff), "SBVA cls cutoff")
    ("sbvalitcut", po::value(&sbva_lits_cutoff)->default_value(sbva_lits_cutoff), "SBVA lits cutoff")
    ("sbvabreak", po::value(&sbva_tiebreak)->default_value(sbva_tiebreak), "SBVA tie breaking. 0 = old, 1 = sbva")
    ("bveresolvmaxsz", po::value(&simp_conf.bve_too_large_resolvent)->default_value(simp_conf.bve_too_large_resolvent), "Puura BVE max resolvent size in literals. -1 == no limit")
    ("bvegrowiter1", po::value(&simp_conf.bve_grow_iter1)->default_value(simp_conf.bve_grow_iter1), "Puura BVE growth allowance iter1")
    ("bvegrowiter2", po::value(&simp_conf.bve_grow_iter2)->default_value(simp_conf.bve_grow_iter2), "Puura BVE growth allowance iter2")

    ("buddymaxcls", po::value(&conf.buddy_max_cls)->default_value(conf.buddy_max_cls), "Run BuDDy")
    ("compsort", po::value(&conf.do_comp_sort)->default_value(conf.do_comp_sort), "Sort components in different order")
    ("varfreqdiv", po::value(&conf.var_freq_divider)->default_value(conf.var_freq_divider), "Var freq divider in score_of")
    ;

    restart_options.add_options()
    ("rstfirst", po::value(&conf.first_restart)->default_value(conf.first_restart), "Run restarts")

    ("restart", po::value(&conf.do_restart)->default_value(conf.do_restart), "Run restarts")
    ("rsttype", po::value(&conf.restart_type)->default_value(conf.restart_type), "Check count at every step")
    ("rstcutoff", po::value(&conf.restart_cutoff_mult)->default_value(conf.restart_cutoff_mult), "Multiply cutoff with this")
    ("rstcheckcnt", po::value(&conf.do_cube_check_count)->default_value(conf.do_cube_check_count), "Check the count of each cube")
    ("rstreadjust", po::value(&conf.do_readjust_for_restart)->default_value(conf.do_readjust_for_restart), "Readjust params for restart")
    ("breakid", po::value(&do_breakid)->default_value(do_breakid), "Enable BreakID")
    ;

    help_options.add(main_options);
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

template<class T> void parse_file(const std::string& filename, T* reader) {
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, nullptr, 0);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, nullptr, 0);
  #endif
  if (in == nullptr) {
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

  if (reader->get_sampl_vars_set() && !ignore_indep) {
    indep_support_given = true;
  } else {
    indep_support_given = false;
    vector<uint32_t> tmp;
    for(uint32_t i = 0; i < reader->nVars(); i++) tmp.push_back(i);
    reader->set_sampl_vars(tmp);
    reader->set_opt_sampl_vars(tmp);
  }
}

vector<Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<Lit> ganak_cl;
  for(const auto& l: cl) ganak_cl.push_back(Lit(l.var()+1, !l.sign()));
  return ganak_cl;
}

double biginteger_log_modified(const mpz_class& x) {
  signed long int ex;
  const double di = mpz_get_d_2exp(&ex, x.get_mpz_t());
  return log10(di) + log10(2) * (double) ex;
}

void print_vars(const vector<uint32_t>& vars) {
  auto tmp = vars;
  std::sort(tmp.begin(), tmp.end());
  for(const auto& v: tmp) cout << v+1 << " ";
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
  parse_supported_options(argc, argv);
  if (conf.verb) {
    cout << ganak_version_info() << endl;
    cout << "c o called with: " << command_line << endl;
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
  ArjunNS::SimplifiedCNF cnf;
  if (!do_arjun) {
    parse_file(fname, &cnf);
    cout << "c o sampl_vars: "; print_vars(cnf.sampl_vars); cout << endl;
    if (cnf.opt_sampl_vars_given) {
      cout << "c o opt sampl_vars: "; print_vars(cnf.opt_sampl_vars); cout << endl;
    }
  } else {
    parse_file(fname, &cnf);
    double my_time = cpuTime();
    ArjunNS::Arjun arjun;
    arjun.set_verb(arjun_verb);
    arjun.set_or_gate_based(arjun_gates);
    arjun.set_xor_gates_based(arjun_gates);
    arjun.set_ite_gate_based(arjun_gates);
    arjun.set_irreg_gate_based(arjun_gates);
    arjun.only_backbone(cnf);
    arjun.only_run_minimize_indep(cnf);
    bool do_extend_indep = true;
    bool do_unate = false;
    arjun.elim_to_file(cnf, indep_support_given, do_extend_indep, do_bce, do_unate, simp_conf, sbva_steps, sbva_cls_cutoff, sbva_lits_cutoff, sbva_tiebreak);
    verb_print(1, "Arjun T: " << (cpuTime()-my_time));
    cout << "c o sampl_vars: "; print_vars(cnf.sampl_vars); cout << endl;
    cout << "c o opt sampl_vars: "; print_vars(cnf.opt_sampl_vars); cout << endl;
  }

  vector<map<Lit, Lit>> generators;
  if (conf.do_restart && do_breakid && cnf.clauses.size() > 1) {
    double my_time = cpuTime();
    BID::BreakID breakid;
    /* breakid.set_useMatrixDetection(conf.useMatrixDetection); */
    /* breakid.set_useFullTranslation(conf.useFullTranslation); */
    breakid.set_verbosity(0);
    breakid.start_dynamic_cnf(cnf.nVars());
    for(const auto& cl: cnf.clauses) {
      breakid.add_clause((BID::BLit*)cl.data(), cl.size());
    }
    breakid.set_steps_lim(4000);
    breakid.end_dynamic_cnf();
    verb_print(1, "[breakid] Num generators: " << breakid.get_num_generators());
    breakid.detect_subgroups();
    if (conf.verb >= 1) breakid.print_generators(std::cout, "c o ");
    vector<unordered_map<BID::BLit, BID::BLit> > orig_gen;
    breakid.get_perms(&orig_gen);
    for(const auto& m: orig_gen) {
      map<Lit, Lit> gen;
      for(const auto& gp: m) {
        gen[Lit(gp.first.var()+1, gp.first.sign())] = Lit(gp.second.var()+1, gp.second.sign());
        if (conf.verb >= 2) cout << "c o " << gp.first << " -> " << gp.second << endl;
      }
      generators.push_back(gen);
    }
    verb_print(1, "[breakid] T: " << (cpuTime()-my_time));
    /* cnf.write_simpcnf("tmp.cnf"); */
  }
  Counter<mpz_class>* counter = new Counter<mpz_class>(conf);
  CMSat::SATSolver* sat_solver = new CMSat::SATSolver;
  counter->new_vars(cnf.nVars());
  counter->set_generators(generators);
  sat_solver->new_vars(cnf.nVars());

  mpz_class cnt = 0;
  for(const auto& cl: cnf.clauses) sat_solver->add_clause(cl);
  for(const auto& cl: cnf.red_clauses) sat_solver->add_clause(cl);
  auto ret = sat_solver->solve();
  if (ret == CMSat::l_True) {
    for(const auto& cl: cnf.clauses) {
      auto cl2 = cms_to_ganak_cl(cl);
      counter->add_irred_cl(cl2);
    }
    counter->end_irred_cls();
    for(const auto& cl: cnf.red_clauses) {
      auto cl2 = cms_to_ganak_cl(cl);
      counter->add_red_cl(cl2);
    }
    set<uint32_t> tmp;
    for(auto const& s: cnf.sampl_vars) tmp.insert(s+1);
    counter->set_indep_support(tmp);
    if (cnf.opt_sampl_vars_given) {
      tmp.clear();
      for(auto const& s: cnf.opt_sampl_vars) tmp.insert(s+1);
      counter->set_optional_indep_support(tmp);
    }

    counter->init_activity_scores();
    cnt = counter->outer_count(sat_solver);
  }
  cout << "c o Total time [Arjun+GANAK]: " << std::setprecision(2) << std::fixed << (cpuTime() - start_time) << endl;

  if (cnt > 0) cout << "s SATISFIABLE" << endl;
  else cout << "s UNSATISFIABLE" << endl;
  if (indep_support_given) cout << "c s type pmc " << endl;
  else cout << "c s type mc" << endl;
  cnt *= cnf.multiplier_weight;
  cout << "c s log10-estimate ";
  if (cnt == 0) {
    cout << "-inf" << endl;
  } else {
    cout << std::setprecision(6) << std::fixed << biginteger_log_modified(cnt) << endl;
  }
  cout << "c s exact arb int " << std::fixed << cnt << endl;


  delete counter;
  delete sat_solver;
  return 0;
}
