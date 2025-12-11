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

#include "ganak.hpp"
#include "GitSHA1.hpp"
#include "common.hpp"
#include "time_mem.hpp"

#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <set>
#include <iomanip>
#include <gmpxx.h>
#include <mpfr.h>
/* #include <breakid/breakid.hpp> */
#include <arjun/arjun.h>
#include "src/argparse.hpp"
#include "mpoly.hpp"
#include "mparity.hpp"
#include "mcomplex.hpp"
#include "mcomplex-mpfr.hpp"
#include <approxmc/approxmc.h>

using CMSat::StreamBuffer;
using CMSat::DimacsParser;
using std::set;
using namespace GanakInt;
using std::setprecision;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

#define add_arg(name, var, fun, hhelp) \
    program.add_argument(name) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)
#define add_arg2(name1, name2, var, fun, hhelp) \
    program.add_argument(name1, name2) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)

using std::string;
using std::vector;
argparse::ArgumentParser program = argparse::ArgumentParser("ganak",
        GANAK::get_version_sha1(),
        argparse::default_arguments::help);
CounterConfiguration conf;
int arjun_verb = 1;
int do_arjun = 1;
int arjun_gates = 1;
/* int do_breakid = 0; */
int arjun_extend_max_confl = 1000;
int do_pre_backbone = 0;
int do_probe_based = 1;
int arjun_simp_level = 2;
int arjun_backw_maxc = 20000;
ArjunNS::Arjun::ElimToFileConf etof_conf;
ArjunNS::SimpConf simp_conf;
string debug_arjun_cnf;
int arjun_oracle_find_bins = 6;
double arjun_cms_glob_mult = -1.0;
int do_puura = 1;
uint32_t arjun_further_min_cutoff = 10;
int arjun_extend_ccnr = 0;
int arjun_autarkies = 0;
int mode = 0;
int poly_nvars = -1;
int prime_field = -1;
int strip_opt_indep = 0;
FG fg = nullptr;
int num_threads = 1;
int bits_jobs = 10;

string print_version()
{
    std::stringstream ss;
    ss << "c o Ganak SHA1: " << GANAK::get_version_sha1() << endl;
    ss << "c o Arjun SHA1: " << ArjunNS::Arjun::get_version_sha1() << endl;
    ss << "c o SBVA SHA1: " << ArjunNS::Arjun::get_sbva_version_sha1() << endl;
    ss << "c o CMS SHA1: " << CMSat::SATSolver::get_version_sha1() << endl;
    ss << "c o ApproxMC SHA1: " << ApproxMC::AppMC::get_version_sha1() << endl;
    /* ss << "c o BreakID SHA1: " << BID::BreakID::get_version_sha1() << endl; */
    ss << ArjunNS::Arjun::get_thanks_info("c o ") << endl;
    ss << CMSat::SATSolver::get_thanks_info("c o ") << endl;
    ss << "c o Using Graph library by Tuukka Korhonen and Matti Jarvisalo" << endl;
    ss << "c o Using Flowcutter by Ben Strasser" << endl;
    ss << "c o Ganak compilation env " << GANAK::get_compilation_env() << endl;
    return ss.str();
}

void add_ganak_options()
{
    std::ostringstream my_delta;
    my_delta << setprecision(8) << conf.delta;

    add_arg2("-v", "--verb", conf.verb, atoi, "Verbosity");
    add_arg2("-s", "--seed", conf.seed, atoi, "Seed");
    program.add_argument("-v", "--version") \
        .action([&](const auto&) {cout << print_version(); exit(0);}) \
        .flag()
        .help("Print version and exit");
    add_arg("--mode", mode , atoi, R"delimiter(0=counting,
1=weighted counting,
2=complex numbers,
3=multivariate polynomials over the rational field,
4=parity counting, 5=counting over prime field,
6=mpfr complex numbers, 7=mpfr normal numbers)delimiter");
    add_arg("--prime", prime_field, atoi, "Prime for prime field counting");
    add_arg("--npolyvars", poly_nvars, atoi, "Number of variables in the polynomial field");
    add_arg("--delta", conf.delta, atof, "Delta");
    /* add_arg("--breakid", do_breakid, atoi, "Enable BreakID"); */
    add_arg("--appmct", conf.appmc_timeout, atof, "after K seconds");
    add_arg("--epsilon", conf.appmc_epsilon, atof, "AppMC epsilon");
    add_arg("--chronobt", conf.do_chronobt, atof, "ChronoBT. SAT must be DISABLED or this will fail");
    add_arg("--prob", conf.do_probabilistic_hashing, atoi, "Use probabilistic hashing. When set to 0, we are not running in probabilistic mode, but in deterministic mode, i.e. delta is 0 in Ganak mode (not in case we switch to ApproxMC mode via --appmct)");

    // Arjun options
    add_arg("--arjun", do_arjun, atoi, "Use arjun");
    add_arg("--arjunverb", arjun_verb, atoi, "Arjun verb");
    add_arg("--arjungates", arjun_gates, atoi, "Use arjun's gate detection");
    add_arg("--arjunextend", etof_conf.do_extend_indep, atoi, "Extend indep via Arjun's extend system");
    add_arg("--prebackbone", do_pre_backbone, atoi, "Perform backbone before other things");
    add_arg("--puura", do_puura, atoi, "Run Puura");
    add_arg("--puurabackbone", simp_conf.do_backbone_puura, atoi, "Perform backbone in Puura");
    add_arg("--arjuniter1", simp_conf.iter1, atoi, "Arjun's iter1");
    add_arg("--arjuniter2", simp_conf.iter2, atoi, "Arjun's iter2");
    add_arg("--arjunprobe", do_probe_based, atoi, "Probe based arjun");
    add_arg("--arjunsimplev", arjun_simp_level, atoi, "Arjun simp level");
    add_arg("--arjunbackwmaxc", arjun_backw_maxc, atoi, "Arjun backw max confl");
    add_arg("--arjunoraclefindbins", arjun_oracle_find_bins, atoi, "Arjun's oracle should find bins or not");
    add_arg("--arjunautarkies", arjun_autarkies, atoi, "How much autarky for Arjun to do");
    add_arg("--bce", etof_conf.do_bce, atoi, "Do static BCE");
    add_arg("--bveresolvmaxsz", simp_conf.bve_too_large_resolvent, atoi, "Puura BVE max resolvent size in literals. -1 == no limit");
    add_arg("--bvegrowiter1", simp_conf.bve_grow_iter1, atoi, "Puura BVE growth allowance iter1");
    add_arg("--bvegrowiter2", simp_conf.bve_grow_iter2, atoi, "Puura BVE growth allowance iter2");
    add_arg("--extraoracle", simp_conf.oracle_extra, atoi, "Extra oracle at the end of puura");
    add_arg("--resolvsub", simp_conf.do_subs_with_resolvent_clauses, atoi, "Sets relevant CMS option: subsume other clauses with resolvent clauses");
    add_arg("--arjunoraclegetlearnt", simp_conf.oracle_vivify_get_learnts, atoi, "Arjun's oracle should get learnts");
    add_arg("--arjundebugcnf", debug_arjun_cnf, string, "Write debug arjun CNF into this file");
    add_arg("--arjuncmsmult", arjun_cms_glob_mult, atof,  "Pass this multiplier to CMSat through Arjun");
    add_arg("--arjunsamplcutoff", arjun_further_min_cutoff, atoi,  "Only perform further arjun-based minimization in case the minimized indep support is larger or equal to this");
    add_arg("--arjunextendccnr", arjun_extend_ccnr, atoi,  "Filter extend of ccnr gates via CCNR mems, in the millions");
    add_arg("--arjunweakenlim", simp_conf.weaken_limit, atoi,  "Arjun's weaken limitation");

    // TD options
    add_arg("--td", conf.do_td, atoi, "Run TD decompose");
    add_arg("--tdmaxw", conf.td_maxweight, atof, "TD max weight");
    add_arg("--tdminw", conf.td_minweight, atof, "TD min weight");
    add_arg("--tddiv", conf.td_divider, atof, "TD divider");
    add_arg("--tdexpmult", conf.td_exp_mult, atof, "TD exponential multiplier");
    add_arg("--tdcheckagainstind", conf.do_check_td_vs_ind, atoi, "Check TD against indep size");
    add_arg("--tditers", conf.td_iters, atoi, "TD flowcutter iterations (restarts)");
    add_arg("--tdsteps", conf.td_steps, atoll, "TD flowcutter number of steps at most");
    add_arg("--tdlook", conf.td_lookahead, atoi, "-1 means never");
    add_arg("--tdlooktwcut", conf.td_lookahead_tw_cutoff, atoi, "TD lookahead only when TW of current comp is larger than this value");
    add_arg("--tdlookiters", conf.td_lookahead_iters, atoi, "TD lookahead iterations");
    add_arg("--tdcontract", conf.do_td_contract, atoi, "TD contract over opt indep set");
    add_arg("--tdlimit", conf.td_limit, atoi, "If TD is over this, reduce weight to 0.1");
    add_arg("--tdoptindep", conf.do_td_use_opt_indep, atoi, "Use opt indep for TD computation");
    add_arg("--tdmaxdensity", conf.td_max_density, atof, "Max density for TD computation");
    add_arg("--tdmaxedgeratio", conf.td_max_edge_var_ratio, atoi, "Max edge to var ratio for TD computation");
    add_arg("--tduseadj", conf.td_do_use_adj, atoi, "TD should use adjacency matrix for computing TD scores");
    add_arg("--tdreadfile", conf.td_read_file, string, "Read TD scores from this file");
    add_arg("--tdvis", conf.td_visualize_dot_file, string, "Visualize the TD into this file in DOT format");

    // Clause DB options
    add_arg("--rdbclstarget", conf.rdb_cls_target, atoi, "RDB clauses target size (added to this are LBD 3 or lower)");
    add_arg("--rdbeveryn", conf.reduce_db_everyN, atoi, "Reduce the clause DB every N conflicts");
    add_arg("--rdbkeepused", conf.rdb_keep_used, atoi, "RDB keeps clauses that are used");
    add_arg("--consolidateeveryn", conf.consolidate_every_n, atoi, "Consolidate memory after every N learnt clause");
    add_arg("--lbd", conf.base_lbd_cutoff, atoi, "Initial LBD cutoff");
    add_arg("--updatelbdcutoff", conf.do_update_lbd_cutoff, atoi, "Update lbd cutoff");

    // Decision options
    add_arg("--polar", conf.polar_type, atoi, "0=standard_polarity, 1=polar cache, 2=false, 3=true");
    add_arg("--decide", conf.decide, atoi, "ignore or not ignore TD");
    add_arg("--initact", conf.do_init_activity_scores, atoi, "Init activity scores to var freq");
    add_arg("--vsadsadjust", conf.vsads_readjust_every, atoi, "VSADS ajust activity every N");
    add_arg("--actscorediv", conf.act_score_divisor, atof, "Activity score divisor");
    add_arg("--freqscorediv", conf.freq_score_divisor, atof, "Component frequency score divisor");

    // Cache options
    add_arg("--cache", conf.do_use_cache, atoi, "Use (i.e. store and retrieve) cache");
    add_arg("--maxcache", conf.maximum_cache_size_MB, atoll, "Max cache size in MB");
    add_arg("--cachetime", conf.cache_time_update, atoi, "2 = set to mid-point");

    // BuDDy options
    add_arg("--buddy", conf.do_buddy, atoi, "Run BuDDy");
    add_arg("--buddymaxcls", conf.buddy_max_cls, atoi, "Run BuDDy");

    // Vivif options -- inprocessing during Ganak
    add_arg("--vivif", conf.do_vivify, atoi, "Vivify clauses");
    add_arg("--vivifevery", conf.vivif_every, atoi, "Vivify every N conflicts");
    add_arg("--vivifmult", conf.vivif_mult, atof, "How much to multiply timeout for vivif");
    add_arg("--vivifoutern", conf.vivif_outer_every_n, atoi, "How many restarts between outer vivif");
    add_arg("--totusedcutoffvivif", conf.tot_used_cutoff_vivif, atoi, "Total used vivif cutoff");

    // SBVA options
    add_arg("--sbvasteps", etof_conf.num_sbva_steps, atoi, "SBVA steps. 0 = no SBVA");
    add_arg("--sbvaclcut", etof_conf.sbva_cls_cutoff, atoi, "SBVA cls cutoff");
    add_arg("--sbvalitcut", etof_conf.sbva_lits_cutoff, atoi, "SBVA lits cutoff");
    add_arg("--sbvabreak", etof_conf.sbva_tiebreak, atoi, "1 = sbva");

    // SAT solver options
    add_arg("--satsolver", conf.do_use_sat_solver, atoi, "Use SAT solver when all minimal indep set has been set");
    add_arg("--satrst", conf.do_sat_restart, atoi, "Inside SAT solver, perform restarts");
    add_arg("--satrstmult", conf.sat_restart_mult, atoi, "SAT restart multiplier");
    add_arg("--satpolarcache", conf.do_sat_polar_cache, atoi, "Inside SAT solver, use polarity cache");
    add_arg("--satvsids", conf.do_sat_vsids, atoi, "Inside SAT solver, use VSIDS, not VSADS");

    // Opt independent set options
    add_arg("--allindep", etof_conf.all_indep, atoi, "All variables can be made part of the indepedent support. Indep support is given ONLY to help the solver.");
    add_arg("--arjunextendmaxconfl", arjun_extend_max_confl, atoi, "Max number of conflicts per extend operation in Arjun");
    add_arg("--arjunextend", etof_conf.do_extend_indep, atoi, "Max number of conflicts per extend operation in Arjun");
    add_arg("--stripoptindep", strip_opt_indep, atoi, "Strip optional indep support");

    // Analyze candidates options
    add_arg("--analyzecand", conf.analyze_cand_update, atoi, "Update analyze candidates if more than N vars are still undecided from opt indep set");

    // Restart options
    add_arg("--rstfirst", conf.first_restart, atoll, "Run restarts");
    add_arg("--restart", conf.do_restart, atoi, "Run restarts");
    add_arg("--rsttype", conf.restart_type, atoi, "Check count at every step");
    add_arg("--rstcutoff", conf.restart_cutoff_mult, atof, "Multiply cutoff with this");
    add_arg("--rstcheckcnt", conf.do_cube_check_count, atoi, "Check the count of each cube");
    add_arg("--rstreadjust", conf.do_readjust_for_restart, atoi, "Readjust params for restart");
    add_arg("--maxrst", conf.max_num_rst, atoi, "Max number of restarts");
    add_arg("--maxcubesperrst", conf.max_num_cubes_per_restart, atoi,  "Max number of cubes per restart");

    // Multi-threading options
    add_arg("--threads", num_threads, atoi, "Number of threads to use. -1 = all available cores");
    add_arg("--bitsjobs", bits_jobs, atoi, "Number of variables to multi-thread on (8 = 256 jobs)");
    program.add_argument("inputfile").remaining().help("input CNF");
}

void parse_supported_options(int argc, char** argv) {
    add_ganak_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout << "Flexible Weighted Model Counter" << endl << endl
            << "ganak [options] inputfile" << endl;
            cout << program << endl;
            exit(0);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        exit(-1);
    }
    if (conf.do_use_sat_solver && !conf.do_chronobt) {
      cout << "ERROR: When chronobt is disabled, SAT solver cannot be used" << endl;
      exit(-1);
    }
    if (bits_jobs < 0 || bits_jobs > 20) {
      cout << "ERROR: bitsjobs must be between 0 and 20, inclusive" << endl;
      exit(-1);
    }
    if (num_threads < -1) {
      cout << "ERROR: number of threads must not be less than -1" << endl;
      exit(-1);
    }
    if (num_threads > 1024) {
      cout << "ERROR: number of threads must not be more than 1024" << endl;
      exit(-1);
    }
    if (num_threads == 0) {
      cout << "ERROR: number of threads must not be 0" << endl;
      exit(-1);
    }
}

template<class T> void parse_file(const std::string& filename, T* reader) {
  #ifndef USE_ZLIB
  FILE * in;
  if (filename == "-") in = stdin;
  else in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, nullptr, 0, fg);
  #else
  gzFile in;
  if (filename == "-") in = gzdopen(fileno(stdin), "rb");
  else in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, nullptr, 0, fg);
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

void print_vars(vector<uint32_t> vars) {
  std::sort(vars.begin(), vars.end());
  for(const auto& v: vars) cout << v+1 << " ";
}

void setup_ganak(const ArjunNS::SimplifiedCNF& cnf, Ganak& counter) {
  cnf.check_sanity();
  counter.new_vars(cnf.nVars());

  set<uint32_t> tmp;
  for(auto const& s: cnf.sampl_vars) tmp.insert(s+1);
  counter.set_indep_support(tmp);
  if (cnf.get_opt_sampl_vars_set()) {
    tmp.clear();
    for(auto const& s: cnf.opt_sampl_vars) tmp.insert(s+1);
  }
  counter.set_optional_indep_support(tmp);
  if (conf.verb) counter.print_indep_distrib();

  if (cnf.weighted) {
    for(const auto& t: cnf.weights) {
      counter.set_lit_weight(Lit(t.first+1, true), t.second.pos);
      counter.set_lit_weight(Lit(t.first+1, false), t.second.neg);
    }
  }

  for(const auto& cl: cnf.clauses) counter.add_irred_cl(cms_to_ganak_cl(cl));
  for(const auto& cl: cnf.red_clauses) counter.add_red_cl(cms_to_ganak_cl(cl));
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
  arjun.set_oracle_find_bins(arjun_oracle_find_bins);
  arjun.set_cms_glob_mult(arjun_cms_glob_mult);
  arjun.set_autarkies(arjun_autarkies);
  if (do_pre_backbone) arjun.standalone_backbone(cnf);
  arjun.standalone_minimize_indep(cnf, etof_conf.all_indep);
  arjun.set_extend_ccnr(arjun_extend_ccnr);
  if (cnf.get_sampl_vars().size() >= arjun_further_min_cutoff && do_puura) {
    arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);
  } else cnf.renumber_sampling_vars_for_ganak();
  verb_print(1, "Arjun T: " << (cpu_time()-my_time));
}

string print_mpq_as_scientific(const mpq_class& number) {
    mpf_class mpf_value(number);
    std::ostringstream oss;
    oss << std::scientific << setprecision(8) << mpf_value;
    return oss.str();
}

void print_log(const mpfr_t& cnt, string extra = "") {
    mpfr_t log10_val;
    mpfr_init2(log10_val, 256);
    mpfr_set(log10_val, cnt, MPFR_RNDN);
    if (mpfr_sgn(log10_val) < 0) {
      cout << "c s neglog10-estimate" << extra << " ";
      mpfr_neg(log10_val, log10_val, MPFR_RNDN);
    } else {
      cout << "c s log10-estimate" << extra << " ";
    }
    mpfr_log10(log10_val, log10_val, MPFR_RNDN);

    char* tmp = nullptr;
    mpfr_asprintf(&tmp, "%.8Re", log10_val);
    cout << tmp << endl;
    mpfr_free_str(tmp);
    mpfr_clear(log10_val);
}

void print_log(const mpz_class& cnt, string extra = "") {
    mpz_class abs_cnt = cnt;
    if (abs_cnt < 0) {
      cout << "c s neglog10-estimate" << extra << " ";
      abs_cnt *= -1;
    } else {
      cout << "c s log10-estimate" << extra << " ";
    }
    mpfr_t log10_val;
    mpfr_init2(log10_val, 256);
    mpfr_set_z(log10_val, abs_cnt.get_mpz_t(), MPFR_RNDN);
    mpfr_log10(log10_val, log10_val, MPFR_RNDN);

    char* tmp = nullptr;
    mpfr_asprintf(&tmp, "%.8Re", log10_val);
    cout << tmp << endl;
    mpfr_free_str(tmp);
    mpfr_clear(log10_val);
}

// compute collision probability, i.e. 2^(log2(lookups) + log2(elems) - 64)
void compute_collision_prob(mpfr_t& result, const uint64_t lookups, uint64_t elems) {
    mpfr_t lookups2;
    mpfr_init_set_ui(lookups2, lookups, MPFR_RNDN);
    mpfr_log2(lookups2, lookups2, MPFR_RNDN);

    mpfr_t elems2;
    mpfr_init_set_ui(elems2, elems, MPFR_RNDN);
    mpfr_log2(elems2, elems2, MPFR_RNDN);

    mpfr_t e;
    mpfr_init_set_si(e, -64, MPFR_RNDN);
    mpfr_add(e, lookups2, e, MPFR_RNDN);
    mpfr_add(e, elems2, e, MPFR_RNDN);
    // e = log2(lookups) + log2(elems) - 64

    // Compute 2^e
    mpfr_init(result);
    mpfr_exp2(result, e, MPFR_RNDN);

    // Clear temporary variables
    mpfr_clear(lookups2);
    mpfr_clear(elems2);
    mpfr_clear(e);
}

void run_weighted_counter(Ganak& counter, const ArjunNS::SimplifiedCNF& cnf, const double start_time) {
    FF cnt = cnf.multiplier_weight->dup();
    if (!cnf.multiplier_weight->is_zero()) *cnt *= *counter.count(bits_jobs, num_threads);
    cout << "c o Total time [Arjun+GANAK]: " << setprecision(2)
        << std::fixed << (cpu_time() - start_time) << endl;

    string out = "c o type ";
    if (cnf.get_projected()) out+="p";
    if (cnf.weighted) out += "wmc";
    else out += "mc";

    if (!cnt->is_zero()) cout << "s SATISFIABLE" << endl;
    else cout << "s UNSATISFIABLE" << endl;
    if (mode == 0 || mode == 1 || mode == 2 || mode == 6 || mode == 7) {
      std::stringstream ss;
      ss << std::scientific << setprecision(40);
      const CMSat::Field* ptr = cnt.get();
      assert(ptr != nullptr);
      if (mode == 0) {
        // Integer numbers
        if (cnf.get_projected()) cout << "c s type pmc" << endl;
        else cout << "c s type mc" << endl;
        const ArjunNS::FMpz* od = dynamic_cast<const ArjunNS::FMpz*>(ptr);
        print_log(od->val);
        ss << *od;
        if (counter.get_is_approximate()) {
          cout << "c s approx arb int "  << ss.str() << endl;
        } else {
          cout << "c s exact arb int "  << ss.str() << endl;
        }
      } else if (mode == 1) {
        // Rational numbers
        if (cnf.get_projected()) cout << "c s type pwmc" << endl;
        else cout << "c s type wmc" << endl;
        const ArjunNS::FMpq* od = dynamic_cast<const ArjunNS::FMpq*>(ptr);
        mpfr_t r;
        mpfr_init2(r, 256);
        mpfr_set_q(r, od->val.get_mpq_t(), MPFR_RNDN);
        print_log(r);
        mpfr_clear(r);

        cout << "c o exact quadruple float "  << print_mpq_as_scientific(od->val) << endl;
        cout << "c s exact arb frac " << *cnt << endl;
      } else if (mode == 2) {
        // Complex rational numbers
        cout << "c s type amc-complex" << endl;
        const FComplex* od = dynamic_cast<const FComplex*>(ptr);
        mpfr_t r, i;
        mpfr_init2(r, 256);
        mpfr_set_q(r, od->real.get_mpq_t(), MPFR_RNDN);
        mpfr_init2(i, 256);
        mpfr_set_q(i, od->imag.get_mpq_t(), MPFR_RNDN);
        print_log(r, "-real");
        print_log(i, "-imag");
        mpfr_clear(r);
        mpfr_clear(i);

        cout << "c o exact quadruple float " << print_mpq_as_scientific(od->real) << " + "
          << print_mpq_as_scientific(od->imag) << "i" << endl;
        cout << "c s exact arb frac " << *cnt << endl;
      } else if (mode == 6) {
        // Complex MPF numbers
        cout << "c s type amc-complex" << endl;
        const MPFComplex* od = dynamic_cast<const MPFComplex*>(ptr);
        print_log(od->real, "-real");
        print_log(od->imag, "-imag");
        mpfr_printf("c s exact quadruple float %.8Re + %.8Rei\n", od->real, od->imag);
      } else if (mode == 7) {
        // MPFR numbers
        if (cnf.get_projected()) cout << "c s type pwmc" << endl;
        else cout << "c s type wmc" << endl;
        const ArjunNS::FMpfr* od = dynamic_cast<const ArjunNS::FMpfr*>(ptr);
        print_log(od->val);
        mpfr_printf("c s exact quadruple float %.8Re\n", od->val);
      }
    }
    if (counter.get_is_approximate()) {
      cout << "c s pac guarantees epsilon: " << conf.appmc_epsilon << " delta: " << conf.delta << endl;
    } else if (counter.get_num_cache_lookups() == 0 || counter.get_max_cache_elems() == 0) {
      cout << "c s pac guarantees epsilon: 0" << " delta: " << 0 << endl;
    } else if (!conf.do_probabilistic_hashing) {
      cout << "c s pac guarantees epsilon: 0 delta: 0" << endl;
    } else {
      mpfr_t collision_prob;
      compute_collision_prob(collision_prob, counter.get_num_cache_lookups(), counter.get_max_cache_elems());
      cout << "c s pac guarantees epsilon: 0" << " delta: ";
      char* tmp = nullptr;
      mpfr_asprintf(&tmp, "%.8Re", collision_prob);
      cout << tmp << endl;
      mpfr_free_str(tmp);
      mpfr_clear(collision_prob);
    }
}

int main(int argc, char *argv[]) {
  mpf_set_default_prec(256);
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
      if (i+1 < argc) command_line += " ";
  }
  parse_supported_options(argc, argv);
  if (conf.verb) {
    cout << print_version();
    cout << "c o called with: " << command_line << endl;
  }

  switch (mode) {
    case 0:
        fg = std::make_unique<ArjunNS::FGenMpz>();
        break;
    case 1:
        fg = std::make_unique<ArjunNS::FGenMpq>();
        break;
    case 7:
        fg = std::make_unique<ArjunNS::FGenMpfr>();
        break;
    case 2:
        fg = std::make_unique<FGenComplex>();
        break;
    case 6:
        fg = std::make_unique<FGenMPFComplex>();
        break;
    case 3:
        if (poly_nvars == -1) {
          cout << "c o [arjun] ERROR: Must provide number of polynomial vars for mode 3 via --npolyvars" << endl;
          exit(-1);
        }
        fg = std::make_unique<FGenPoly>(poly_nvars);
        break;
    case 4:
        fg = std::make_unique<FGenParity>();
        break;
    case 5:
        if (prime_field == -1) {
          cout << "c o [arjun] ERROR: Must provide prime field for mode 5 via --prime" << endl;
          exit(-1);
        }
        fg = std::make_unique<FGenPrime>(prime_field);
        break;
    default:
        cout << "c o [arjun] ERROR: Unknown mode" << endl;
        exit(-1);
  }
  ArjunNS::SimplifiedCNF cnf(fg);

  // Parse the CNF
  if (!program.is_used("inputfile")) parse_file("-",  &cnf);
  else {
    auto files = program.get<std::vector<std::string>>("inputfile");
    if (files.empty()) {
      cout << "ERROR: you provided --inputfile but no file. Strange. Exiting. " << endl;
      exit(-1);
    } else if (files.size() == 1) {
      const string& fname = files[0];
      parse_file(fname, &cnf);
    } else {
        cout << "[appmc] ERROR: you must only give one CNF as input (or none, and then we read from STDIN)" << endl;
        exit(-1);
    }
  }

  if (cnf.get_weighted() && conf.do_buddy) {
    cout << "ERROR: Cannot run BuDDy with weighted CNF" << endl;
    exit(-1);
  }
  cnf.clean_idiotic_mccomp_weights();
  cnf.check_sanity();
  verb_print(1, "CNF projection set size: " << cnf.get_sampl_vars().size());

  // Run Arjun
  if (!do_arjun) cnf.renumber_sampling_vars_for_ganak();
  else run_arjun(cnf);
  cnf.remove_equiv_weights();
  if (strip_opt_indep) cnf.strip_opt_sampling_vars();
  if (conf.verb >= 2) {
    cout << "c o sampl_vars: "; print_vars(cnf.sampl_vars); cout << endl;
    if (cnf.get_opt_sampl_vars_set()) {
      cout << "c o opt sampl_vars: "; print_vars(cnf.opt_sampl_vars); cout << endl;
    }
  }

  /* // Run BreakID */
  /* vector<map<Lit, Lit>> generators; */
  /* if (cnf.get_sampl_vars().size() >= arjun_further_min_cutoff && conf.do_restart && do_breakid && cnf.clauses.size() > 1) */
  /*   generators = run_breakid(cnf); */

  if (!debug_arjun_cnf.empty()) cnf.write_simpcnf(debug_arjun_cnf, true);

  // Run Ganak
  Ganak counter(conf, fg);
  setup_ganak(cnf, counter);
  run_weighted_counter(counter, cnf, start_time);
  return 0;
}
