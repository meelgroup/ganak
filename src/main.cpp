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
#include <arjun/arjun.h>

using CMSat::StreamBuffer;
using CMSat::DimacsParser;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

/* #define SIMPLE */

namespace po = boost::program_options;
using std::string;
using std::vector;
po::options_description main_options = po::options_description("Main options");
#ifndef SIMPLE
po::options_description restart_options = po::options_description("Restart options");
#endif
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

using namespace std;
bool indep_support_given = false;
int do_check = 0;
CounterConfiguration conf;
int optional_indep = 1;
int arjun_verb = 0;
int do_arjun = 1;
int ignore_indep = 0;
int sbva_steps = 200;

struct CNFHolder {
  vector<vector<CMSat::Lit>> clauses;
  vector<vector<CMSat::Lit>> red_clauses;
  vector<uint32_t> sampling_vars;
  vector<uint32_t> optional_sampling_vars;
  uint32_t nvars = 0;

  uint32_t nVars() const { return nvars; }
  uint32_t new_vars(uint32_t vars) { nvars+=vars; return nvars; }
  uint32_t new_var() { nvars++; return nvars;}
  void add_xor_clause(vector<uint32_t>&, bool) { exit(-1); }
  void add_xor_clause(vector<CMSat::Lit>&, bool) { exit(-1); }
  void add_clause(vector<CMSat::Lit>& cl) { clauses.push_back(cl); }
  void add_red_clause(vector<CMSat::Lit>& cl) { red_clauses.push_back(cl); }
  uint32_t must_mult_exp2 = 0;
};
CNFHolder cnfholder;

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
    ("arjun", po::value(&do_arjun)->default_value(do_arjun), "Use arjun")
    ("arjunverb", po::value(&arjun_verb)->default_value(arjun_verb), "Arjun verb")
    ("ignore", po::value(&ignore_indep)->default_value(ignore_indep), "Ignore indep support given")
    ("optind", po::value(&optional_indep)->default_value(optional_indep), "Ignore indep support given")
    ("singlebump", po::value(&conf.do_single_bump)->default_value(conf.do_single_bump), "Do single bumping, no double (or triple, etc) bumping of activities. Non-single bump is old ganak")
    ("forcebranch", po::value(&conf.force_branch)->default_value(conf.force_branch), "Force branch. 0 = no force, 1 = TD priority, 2 = conflict priority")
#ifndef SIMPLE
    ("extraclbump", po::value(&conf.do_extra_cl_bump)->default_value(conf.do_extra_cl_bump), "Also bump clauses when they propagate. By bump, we mean: set 'used' flag, and update LBD")
    ("td", po::value(&conf.do_td)->default_value(conf.do_td), "Run TD decompose")
    ("tdwithredbins", po::value(&conf.td_with_red_bins)->default_value(conf.td_with_red_bins), "TD computation with redundant binaries")

    ("rdbclstarget", po::value(&conf.rdb_cls_target)->default_value(conf.rdb_cls_target), "RDB clauses target size (added to this are LBD 3 or lower)")
    ("rdbkeepused", po::value(&conf.rdb_keep_used)->default_value(conf.rdb_keep_used), "RDB keeps clauses that are used")
    ("cscore", po::value(&conf.do_cache_hit_scores)->default_value(conf.do_cache_hit_scores), "Do cache scores")
    ("cache", po::value(&conf.do_use_cache)->default_value(conf.do_use_cache), "Use (i.e. store and retrieve) cache")
    ("maxcache", po::value(&conf.maximum_cache_size_MB)->default_value(conf.maximum_cache_size_MB), "Max cache size in MB. 0 == use 80% of free mem")
    ("actexp", po::value(&conf.act_exp)->default_value(conf.act_exp), "Probabilistic Comp Caching")
    ("version", "Print version info")
    ("check", po::value(&do_check)->default_value(do_check), "Check count at every step")
    ("alluipincact", po::value(&conf.alluip_inc_act)->default_value(conf.alluip_inc_act), "All UIP should increase activities")
    ("polar", po::value(&conf.polar_type)->default_value(conf.polar_type),
     "Use polarity cache. Otherwise, false default polar.")
    ("vivif", po::value(&conf.do_vivify)->default_value(conf.do_vivify), "Vivify clauses")
    ("vivifevery", po::value(&conf.vivif_every)->default_value(conf.vivif_every), "Vivify every N conflicts")
    ("vivifmult", po::value(&conf.vivif_mult)->default_value(conf.vivif_mult), "How much to multiply timeout for vivif")
    ("buddy", po::value(&conf.do_buddy)->default_value(conf.do_buddy), "Run BuDDy")
    ("decide", po::value(&conf.decide)->default_value(conf.decide), "Decision type. 0 = sstd-inspired, 1 = gpmc-inspired")
    ("sbva", po::value(&sbva_steps)->default_value(sbva_steps), "SBVA steps. 0 = no SBVA")
    ("buddymaxcls", po::value(&conf.buddy_max_cls)->default_value(conf.buddy_max_cls), "Run BuDDy")
    ("comprevsort", po::value(&conf.do_comp_reverse_sort)->default_value(conf.do_comp_reverse_sort), "Sort components in reverse order")
#endif
    ;

#ifndef SIMPLE
    restart_options.add_options()
    ("rstfirst", po::value(&conf.first_restart)->default_value(conf.first_restart), "Run restarts")

    ("restart", po::value(&conf.do_restart)->default_value(conf.do_restart), "Run restarts")
    ("rsttype", po::value(&conf.restart_type)->default_value(conf.restart_type), "Check count at every step")
    ("rstnext", po::value(&conf.next_restart)->default_value(conf.next_restart), "Restart interval (for type 3)")
    ("rstcutoff", po::value(&conf.restart_cutoff_mult)->default_value(conf.restart_cutoff_mult), "Multiply cutoff with this")
    ;
#endif

    help_options.add(main_options);
#ifndef SIMPLE
    help_options.add(restart_options);
#endif
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

void write_simpcnf(const ArjunNS::SimplifiedCNF& simpcnf,
        const std::string& elimtofile, const uint32_t orig_cnf_must_mult_exp2,
        bool red = true)
{
    uint32_t num_cls = simpcnf.cnf.size();
    std::ofstream outf;
    outf.open(elimtofile.c_str(), std::ios::out);
    outf << "p cnf " << simpcnf.nvars << " " << num_cls << endl;

    //Add projection
    outf << "c p show ";
    for(const auto& v: simpcnf.sampling_vars) {
        assert(v < simpcnf.nvars);
        outf << v+1  << " ";
    }
    outf << "0\n";

    for(const auto& cl: simpcnf.cnf) outf << cl << " 0\n";
    if (red) for(const auto& cl: simpcnf.red_cnf) outf << "c red " << cl << " 0\n";
    outf << "c MUST MULTIPLY BY 2**" << simpcnf.empty_occs+orig_cnf_must_mult_exp2 << endl;
}

template<class T>
void parse_file(const std::string& filename, T* reader) {
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

  indep_support_given = parser.sampling_vars_found && !ignore_indep;
  if (parser.sampling_vars_found && !ignore_indep) {
    cnfholder.sampling_vars = parser.sampling_vars;
    cnfholder.optional_sampling_vars = parser.sampling_vars;
  } else {
    for(uint32_t i = 0; i < reader->nVars(); i++) cnfholder.sampling_vars.push_back(i);
    if (optional_indep && parser.sampling_vars_found)
      cnfholder.optional_sampling_vars = parser.sampling_vars;
    else
      cnfholder.optional_sampling_vars = cnfholder.sampling_vars;
  }
  cnfholder.must_mult_exp2 = parser.must_mult_exp2;
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
#ifdef SIMPLE
  conf.verb = 0;
#endif

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
  if (!do_arjun) {
    parse_file(fname, &cnfholder);
  } else {
    double my_time = cpuTime();
    ArjunNS::Arjun* arjun = new ArjunNS::Arjun;
    arjun->set_seed(conf.seed);
    arjun->set_verbosity(arjun_verb);
    parse_file(fname, arjun);
    arjun->set_starting_sampling_set(cnfholder.sampling_vars);
    cnfholder.sampling_vars = arjun->get_indep_set();
    ArjunNS::SimpConf simp_conf;
    auto ret = arjun->get_fully_simplified_renumbered_cnf(
            cnfholder.sampling_vars, simp_conf, true, false);

    arjun->set_verbosity(1);
    arjun->run_sbva(ret, sbva_steps, 3, 3);
    cnfholder = CNFHolder();
    delete arjun;
    // Extend
    {
      ArjunNS::Arjun arj2;
      arj2.new_vars(ret.nvars);
      arj2.set_verbosity(arjun_verb);
      for(const auto& cl: ret.cnf) arj2.add_clause(cl);
      arj2.set_starting_sampling_set(ret.sampling_vars);
      ret.optional_sampling_vars = arj2.extend_indep_set();
    }
    if (ret.sampling_vars.size() != ret.nvars) ret.renumber_sampling_vars_for_ganak();
    if (!indep_support_given) {
      for(uint32_t i = 0; i < ret.nvars; i++) ret.optional_sampling_vars.push_back(i);
    }
    verb_print(1, "Arjun T: " << (cpuTime()-my_time));
    SLOW_DEBUG_DO(write_simpcnf(ret, fname+"-simplified.cnf", ret.empty_occs, true));

    // set up cnfholder
    cnfholder = CNFHolder();
    cnfholder.clauses = ret.cnf;
    cnfholder.red_clauses = ret.red_cnf;
    cnfholder.nvars = ret.nvars;
    cnfholder.optional_sampling_vars = ret.optional_sampling_vars;
    cnfholder.sampling_vars = ret.sampling_vars;
    cnfholder.must_mult_exp2 += ret.empty_occs;
  }
  Counter* counter = new Counter(conf);
  CMSat::SATSolver* sat_solver = new CMSat::SATSolver;
  counter->new_vars(cnfholder.nVars());
  sat_solver->new_vars(cnfholder.nVars());

  mpz_class cnt = 0;
  for(const auto& cl: cnfholder.clauses) sat_solver->add_clause(cl);
  auto ret = sat_solver->solve();
  if (ret == CMSat::l_True) {
    for(const auto& cl: cnfholder.clauses) {
      auto cl2 = cms_to_ganak_cl(cl);
      counter->add_irred_cl(cl2);
    }
    counter->end_irred_cls();
    for(const auto& cl: cnfholder.red_clauses) {
      auto cl2 = cms_to_ganak_cl(cl);
      counter->add_red_cl(cl2);
    }
    set<uint32_t> tmp;
    for(auto const& s: cnfholder.sampling_vars) tmp.insert(s+1);
    counter->set_indep_support(tmp);
    if (optional_indep || !indep_support_given) {
      tmp.clear();
      for(auto const& s: cnfholder.optional_sampling_vars) tmp.insert(s+1);
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
  for(uint32_t i = 0; i < cnfholder.must_mult_exp2; i++) cnt *= 2;
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
