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
#include <gmpxx.h>
#include "src/GitSHA1.hpp"
#include "breakid.hpp"
#include <arjun/arjun.h>
#include "mpreal.h"
#include "src/argparse.hpp"

using CMSat::StreamBuffer;
using CMSat::DimacsParser;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

#define myopt(name, var, fun, hhelp) \
    program.add_argument(name) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)
#define myopt2(name1, name2, var, fun, hhelp) \
    program.add_argument(name1, name2) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)

using std::string;
using std::vector;
argparse::ArgumentParser program = argparse::ArgumentParser("ganak");

using namespace std;
CounterConfiguration conf;
int arjun_verb = 2;
int do_arjun = 1;
int arjun_gates = 1;
int do_breakid = 0;
int arjun_extend_max_confl = 1000;
int do_pre_backbone = 0;
int do_probe_based = 1;
int arjun_simp_level = 2;
int arjun_backw_maxc = 20000;
ArjunNS::Arjun::ElimToFileConf etof_conf;
ArjunNS::SimpConf simp_conf;
string debug_arjun_cnf;
int do_precise = 1;
int do_backbone_only_optindep = 0;
int arjun_oracle_find_bins = 6;
double arjun_cms_mult = -1.0;
int do_puura = 1;

string print_version()
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
    std::ostringstream my_delta;
    my_delta << std::setprecision(8) << conf.delta;

    myopt2("-v", "--verb", conf.verb, atoi, "Verbosity");
    myopt2("-s", "--seed", conf.seed, atoi, "Seed");
    program.add_argument("-v", "--version") \
        .action([&](const auto&) {print_version(); exit(0);}) \
        .flag()
        .help("Print version and exit");
    myopt("--delta", conf.delta, atof, "Delta");
    myopt("--arjun", do_arjun, atoi, "Use arjun");
    myopt("--arjunverb", arjun_verb, atoi, "Arjun verb");
    myopt("--arjungates", arjun_gates, atoi, "Use arjun's gate detection");
    myopt("--arjunextend", etof_conf.do_extend_indep, atoi, "Extend indep via Arjun's extend system");
    myopt("--arjunextendmaxconfl", arjun_extend_max_confl, atoi, "Max number of conflicts per extend operation in Arjun");
    myopt("--prebackbone", do_pre_backbone, atoi, "Perform backbone before other things");
    myopt("--puura", do_puura, atoi, "Run Puura");
    myopt("--backbonepuura", simp_conf.do_backbone_puura, atoi, "Perform backbone in Puura");
    myopt("--arjunprobe", do_probe_based, atoi, "Probe based arjun");
    myopt("--backboneonlyoptind", do_backbone_only_optindep, atoi, "Backbone only over the opt indep set");
    myopt("--arjunsimplev", arjun_simp_level, atoi, "Arjun simp level");
    myopt("--arjunbackwmaxc", arjun_backw_maxc, atoi, "Arjun backw max confl");
    myopt("--arjunoraclefindbins", arjun_oracle_find_bins, atoi, "Arjun's oracle should find bins or not");
    myopt("--allindep", etof_conf.all_indep, atoi, "All variables can be made part of the indepedent support. Indep support is given ONLY to help the solver.");
    myopt("--td", conf.do_td, atoi, "Run TD decompose");
    myopt("--tdmaxw", conf.td_maxweight, atof, "TD max weight");
    myopt("--tdminw", conf.td_minweight, atof, "TD min weight");
    myopt("--tddiv", conf.td_divider, atof, "TD divider");
    myopt("--tdweighted", conf.do_td_weight, atof, "TD weight enabled");
    myopt("--tdexpmult", conf.td_exp_mult, atof, "TD exponential multiplier");
    myopt("--tdcheckagainstind", conf.do_check_td_vs_ind, atoi, "Check TD against indep size");
    myopt("--tditers", conf.td_iters, atoi, "TD flowcutter iterations (restarts)");
    myopt("--tdsteps", conf.td_steps, atoll, "TD flowcutter number of steps at most");
    myopt("--tdlook", conf.td_lookahead, atoi, "-1 means never");
    myopt("--tdlooktwcut", conf.td_lookahead_tw_cutoff, atoi, "TD lookahead only when TW of current comp is larger than this value");
    myopt("--tdlookiters", conf.td_lookahead_iters, atoi, "TD lookahead iterations");
    myopt("--tdlookonlyweight", conf.td_look_only_weight, atoi, "TD lookahead ONLY update weights");
//
    myopt("--rdbclstarget", conf.rdb_cls_target, atoi, "RDB clauses target size (added to this are LBD 3 or lower)");
    myopt("--rdbeveryn", conf.reduce_db_everyN, atoi, "Reduce the clause DB every N conflicts");
    myopt("--rdbkeepused", conf.rdb_keep_used, atoi, "RDB keeps clauses that are used");
    myopt("--consolidateeveryn", conf.consolidate_every_n, atoi, "Consolidate every N learnt clause");
    myopt("--lbd", conf.base_lbd_cutoff, atoi, "Initial LBD cutoff");
    myopt("--updatelbdcutoff", conf.update_lbd_cutoff, atoi, "Update lbd cutoff");
    myopt("--totusedcutoffvivif", conf.tot_used_cutoff_vivif, atoi, "Total used vivif cutoff");
//
    myopt("--bce", etof_conf.do_bce, atoi, "Do static BCE");
    myopt("--cache", conf.do_use_cache, atoi, "Use (i.e. store and retrieve) cache");
    myopt("--maxcache", conf.maximum_cache_size_MB, atoll, "Max cache size in MB. 0 == use 80% of free mem");
    myopt("--polar", conf.polar_type, atoi, "Use polarity cache. Otherwise, false default polar");
    myopt("--vivif", conf.do_vivify, atoi, "Vivify clauses");
    myopt("--vivifevery", conf.vivif_every, atoi, "Vivify every N conflicts");
    myopt("--vivifmult", conf.vivif_mult, atof, "How much to multiply timeout for vivif");
    myopt("--vivifoutern", conf.vivif_outer_every_n, atoi, "How many restarts between outer vivif");
//
    myopt("--buddy", conf.do_buddy, atoi, "Run BuDDy");
    myopt("--decide", conf.decide, atoi, "1 = gpmc-inspired");
    myopt("--cachetime", conf.cache_time_update, atoi, "2 = set to mid-point");
    myopt("--sbvasteps", etof_conf.num_sbva_steps, atoi, "SBVA steps. 0 = no SBVA");
    myopt("--sbvaclcut", etof_conf.sbva_cls_cutoff, atoi, "SBVA cls cutoff");
    myopt("--sbvalitcut", etof_conf.sbva_lits_cutoff, atoi, "SBVA lits cutoff");
    myopt("--sbvabreak", etof_conf.sbva_tiebreak, atoi, "1 = sbva");
    myopt("--bveresolvmaxsz", simp_conf.bve_too_large_resolvent, atoi, "Puura BVE max resolvent size in literals. -1 == no limit");
    myopt("--bvegrowiter1", simp_conf.bve_grow_iter1, atoi, "Puura BVE growth allowance iter1");
    myopt("--bvegrowiter2", simp_conf.bve_grow_iter2, atoi, "Puura BVE growth allowance iter2");
    myopt("--extraoracle", simp_conf.oracle_extra, atoi, "Extra oracle at the end of puura");
    myopt("--resolvsub", simp_conf.do_subs_with_resolvent_clauses, atoi, "Extra oracle at the end of puura");
    myopt("--arjunoraclegetlearnt", simp_conf.oracle_vivify_get_learnts, atoi, "Arjun's oracle should get learnts");
//
    myopt("--buddymaxcls", conf.buddy_max_cls, atoi, "Run BuDDy");
    myopt("--initact", conf.do_init_activity_scores, atoi, "Init activity scores to var freq");
    myopt("--vsadsadjust", conf.vsads_readjust_every, atoi, "VSADS ajust activity every N");
    myopt("--satrstmult", conf.sat_restart_mult, atoi, "SAT restart multiplier");
    myopt("--precise", do_precise, atoi, "If set to 0, we use so-called 'infinite' precision foats that are far from infinite precision and can give nonsense results. Apparently acceptable by the model counting competition (why? how? what?)");
//
    myopt("--rstfirst", conf.first_restart, atoll, "Run restarts");
    myopt("--restart", conf.do_restart, atoi, "Run restarts");
    myopt("--rsttype", conf.restart_type, atoi, "Check count at every step");
    myopt("--rstcutoff", conf.restart_cutoff_mult, atof, "Multiply cutoff with this");
    myopt("--rstcheckcnt", conf.do_cube_check_count, atoi, "Check the count of each cube");
    myopt("--rstreadjust", conf.do_readjust_for_restart, atoi, "Readjust params for restart");
    myopt("--breakid", do_breakid, atoi, "Enable BreakID");
    myopt("--appmct", conf.appmc_timeout, atof, "after K seconds");
    myopt("--epsilon", conf.appmc_epsilon, atof, "AppMC epsilon");
    myopt("--maxrst", conf.max_num_rst, atoi, "Max number of restarts");
    myopt("--debugarjuncnf", debug_arjun_cnf, string, "Write debug arjun CNF into this file");
    myopt("--arjuncmsmult", arjun_cms_mult, atof,  "Pass this multiplier to CMSat through Arjun");
    myopt("--maxcubesperrst", conf.max_num_cubes_per_restart, atoi,  "Max number of cubes per restart");
    program.add_argument("inputfile").remaining().help("input CNF");
}

void parse_supported_options(int argc, char** argv)
{
    add_ganak_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout << "Probilistic Approcimate Counter" << endl << endl
            << "approxmc [options] inputfile" << endl;
            cout << program << endl;
            exit(0);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        exit(-1);
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

  if (!reader->get_sampl_vars_set()) {
    etof_conf.all_indep = true;
    vector<uint32_t> tmp;
    for(uint32_t i = 0; i < reader->nVars(); i++) tmp.push_back(i);
    reader->set_sampl_vars(tmp); // will automatically set the opt_sampl_vars
  } else {
    // Check if CNF has all vars as indep. Then its's all_indep
    set<uint32_t> tmp;
    for(auto const& s: reader->get_sampl_vars()) {
      if (s >= reader->nVars()) {
        cout << "ERROR: Sampling var " << s+1 << " is larger than number of vars in formula: "
          << reader->nVars() << endl;
        exit(-1);
      }
      tmp.insert(s);
    }
    if (tmp.size() == reader->nVars()) etof_conf.all_indep = true;
    if (!reader->get_opt_sampl_vars_set()) {
      reader->set_opt_sampl_vars(reader->get_sampl_vars());
    }
  }
}

vector<Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<Lit> ganak_cl; ganak_cl.reserve(cl.size());
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


void setup_ganak(const ArjunNS::SimplifiedCNF& cnf, vector<map<Lit, Lit>>& generators,
    OuterCounter& counter) {
  counter.new_vars(cnf.nVars());
  counter.set_generators(generators);

  set<uint32_t> tmp;
  for(auto const& s: cnf.sampl_vars) tmp.insert(s+1);
  counter.set_indep_support(tmp);
  if (cnf.get_opt_sampl_vars_set()) {
    tmp.clear();
    for(auto const& s: cnf.opt_sampl_vars) tmp.insert(s+1);
    counter.set_optional_indep_support(tmp);
  }

  if (cnf.weighted) {
    for(const auto& t: cnf.weights) {
      counter.set_lit_weight(Lit(t.first+1, true), t.second.pos);
      counter.set_lit_weight(Lit(t.first+1, false), t.second.neg);
    }
  }

  for(const auto& cl: cnf.clauses) counter.add_irred_cl(cms_to_ganak_cl(cl));
  for(const auto& cl: cnf.red_clauses) counter.add_red_cl(cms_to_ganak_cl(cl));
  counter.end_irred_cls();
}

auto run_breakid(const ArjunNS::SimplifiedCNF& cnf) {
  double my_time = cpu_time();
  vector<map<Lit, Lit>> generators;
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
  verb_print(1, "[breakid] T: " << (cpu_time()-my_time));
  return generators;
}

void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  double my_time = cpu_time();
  ArjunNS::Arjun arjun;
  if (conf.verb == 0) arjun_verb = 0;
  arjun.set_verb(arjun_verb);
  arjun.set_or_gate_based(arjun_gates);
  arjun.set_xor_gates_based(arjun_gates);
  arjun.set_ite_gate_based(arjun_gates);
  arjun.set_irreg_gate_based(arjun_gates);
  arjun.set_extend_max_confl(arjun_extend_max_confl);
  arjun.set_probe_based(do_probe_based);
  arjun.set_simp(arjun_simp_level);
  arjun.set_backw_max_confl(arjun_backw_maxc);
  arjun.set_backbone_only_optindep(do_backbone_only_optindep);
  arjun.set_oracle_find_bins(arjun_oracle_find_bins);
  arjun.set_cms_mult(arjun_cms_mult);
  if (do_pre_backbone) arjun.standalone_backbone(cnf);
  arjun.standalone_minimize_indep(cnf);
  if (do_puura) arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);
  else cnf.renumber_sampling_vars_for_ganak();
  if (etof_conf.all_indep) {
    cnf.opt_sampl_vars.clear();
    for(uint32_t i = 0; i < cnf.nVars(); i++) cnf.opt_sampl_vars.push_back(i);
  }
  verb_print(1, "Arjun T: " << (cpu_time()-my_time));
}

template<typename T>
void run_weighted_counter(OuterCounter& counter, const ArjunNS::SimplifiedCNF& cnf, const double start_time) {
    T cnt;
    static constexpr bool precise = std::is_same<T, mpq_class>::value;
    if (cnf.multiplier_weight == 0) cnt = 0;
    else {
      if constexpr (!precise) {
          cnt = counter.w_outer_count();
      } else {
          cnt = counter.wq_outer_count();
      }
    }
    cout << "c o Total time [Arjun+GANAK]: " << std::setprecision(2)
      << std::fixed << (cpu_time() - start_time) << endl;
    if (!cnf.get_projected()) cout << "c s type wmc" << endl;
    else cout << "c s type pwmc " << endl;

    if constexpr(!precise) {
      cnt *= cnf.multiplier_weight.get_mpq_t();
    } else {
      cnt *= cnf.multiplier_weight;
    }

    bool neglog = false;
    if (cnt != 0) cout << "s SATISFIABLE" << endl;
    else cout << "s UNSATISFIABLE" << endl;
    if (cnt == 0) cout << "c s log10-estimate -inf" << endl;
    else {
      if (cnt < 0) {
        cout << "c s neglog10-estimate ";
        cnt *= -1;
        neglog = true;
      } else {
        cout << "c s log10-estimate ";
      }
      if constexpr (!precise) {
        cout << std::setprecision(12) << std::fixed << mpfr::log10(cnt) << endl;
        if (neglog) cnt *= -1;
        cout << "c s exact arb float " << std::scientific << std::setprecision(40) << cnt << endl;
      } else {
        cout << std::setprecision(12) << std::fixed << mpfr::log10(cnt.get_mpq_t()) << endl;
        if (neglog) cnt *= -1;
        cout << "c s exact arb float " << std::scientific << std::setprecision(40) << std::flush;
        mpf_set_default_prec(1024); // Set default precision in bits
        mpf_t f;
        mpf_init(f);
        mpf_set_q(f, cnt.get_mpq_t());
        uint32_t n = 50;
        gmp_printf("%.*FE", n, f);
        std::flush(std::cout);
        mpf_clear(f);
        cout << endl;
        cout << "c o exact arb rational " << std::scientific << std::setprecision(40) << cnt << endl;
      }
    }
}

int main(int argc, char *argv[])
{
  const double start_time = cpu_time();
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
    cout << print_version() << endl;
    cout << "c o called with: " << command_line << endl;
  }

  string fname;

  auto files = program.get<std::vector<std::string>>("inputfile");
  if (files.empty()) {
    // TODO read stdin, once we are a library.
    cout << "ERROR: must give input file to read" << endl;
    exit(-1);
  }
  if (files.size() > 1) {
      cout << "[appmc] ERROR: you must only give one CNF as input" << endl;
      exit(-1);
  }
  fname = files[0];
  ArjunNS::SimplifiedCNF cnf;
  parse_file(fname, &cnf);
  if (cnf.get_weighted() && conf.do_buddy) {
    cout << "ERROR: Cannot run BuDDy with weighted CNF" << endl;
    exit(-1);
  }
  if (!do_arjun) cnf.renumber_sampling_vars_for_ganak();
  else run_arjun(cnf);
  if (conf.verb) {
    cout << "c o sampl_vars: "; print_vars(cnf.sampl_vars); cout << endl;
    if (cnf.get_opt_sampl_vars_set()) {
      cout << "c o opt sampl_vars: "; print_vars(cnf.opt_sampl_vars); cout << endl;
    }
  }
  vector<map<Lit, Lit>> generators;
  if (conf.do_restart && do_breakid && cnf.clauses.size() > 1)
    generators = run_breakid(cnf);
  if (!debug_arjun_cnf.empty()) cnf.write_simpcnf(debug_arjun_cnf, true, true);

  OuterCounter counter(conf, cnf.weighted, do_precise);
  setup_ganak(cnf, generators, counter);

  if (cnf.weighted && !do_precise) {
    run_weighted_counter<mpfr::mpreal>(counter, cnf, start_time);
  } else if (cnf.weighted && do_precise) {
    run_weighted_counter<mpq_class>(counter, cnf, start_time);
  } else {
    mpz_class cnt;
    if (cnf.multiplier_weight == 0) cnt = 0;
    else cnt = counter.unw_outer_count();
    cout << "c o Total time [Arjun+GANAK]: " << std::setprecision(2)
      << std::fixed << (cpu_time() - start_time) << endl;
    bool is_appx = counter.get_is_approximate();

    if (cnt != 0) cout << "s SATISFIABLE" << endl;
    else cout << "s UNSATISFIABLE" << endl;
    if (!cnf.get_projected()) cout << "c s type mc" << endl;
    else cout << "c s type pmc " << endl;
    cnt *= cnf.multiplier_weight;
    cout << "c s log10-estimate ";
    if (cnt == 0) cout << "-inf" << endl;
    else cout << std::setprecision(12) << std::fixed << biginteger_log_modified(cnt) << endl;
    if (is_appx) {
      cout << "c s pac guarantees epsilon: " << conf.appmc_epsilon << " delta: " << conf.delta << endl;
      cout << "c s approx arb int " << std::fixed << cnt << endl;
    } else
      cout << "c s exact arb int " << std::fixed << cnt << endl;
  }
  return 0;
}
