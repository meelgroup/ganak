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
#include <charconv>
/* #include <breakid/breakid.hpp> */
#include <arjun/arjun.h>
#include "src/argparse.hpp"
#include "mpoly.hpp"
#include "mparity.hpp"
#include "mcomplex.hpp"
#include "mcomplex-mpfr.hpp"
#include "fmpfi.hpp"
#include "fmpqi.hpp"
#include <approxmc/approxmc.h>
#include "file_read_helper.h"

static constexpr uint32_t max_digit_precision = 1e6;

using std::set;
using namespace GanakInt;
using std::setprecision;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

static int fc_int(const std::string& s) {
    int val = 0;
    std::from_chars(s.data(), s.data() + s.size(), val);
    return val;
}
static double fc_double(const std::string& s) {
    size_t pos;
    double val = std::stod(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("trailing characters in double: " + s);
    return val;
}
static const std::string& fc_string(const std::string& s) { return s; }

using std::string;
using std::vector;
argparse::ArgumentParser program = argparse::ArgumentParser("ganak",
        GANAK::get_version_sha1(),
        argparse::default_arguments::help);

template<typename T, typename F>
void add_arg(const char* name, T& var, F fun, const char* hhelp) {
    using r = std::decay_t<std::invoke_result_t<F, const std::string&>>;
    static_assert(std::is_floating_point_v<r> == std::is_floating_point_v<T>,
        "Floating-point mismatch: use fc_double for floating-point vars, fc_int for integral vars");
    static_assert(std::is_integral_v<r> == std::is_integral_v<T>,
        "Integral/string mismatch: use fc_int for integral vars, fc_string for string vars");
    program.add_argument(name)
        .action([&var, fun](const auto& a) { var = fun(a); })
        .default_value(var)
        .help(hhelp);
}
template<typename T, typename F>
void add_arg2(const char* name1, const char* name2, T& var, F fun, const char* hhelp) {
    using r = std::decay_t<std::invoke_result_t<F, const std::string&>>;
    static_assert(std::is_floating_point_v<r> == std::is_floating_point_v<T>,
        "Floating-point mismatch: use fc_double for floating-point vars, fc_int for integral vars");
    static_assert(std::is_integral_v<r> == std::is_integral_v<T>,
        "Integral/string mismatch: use fc_int for integral vars, fc_string for string vars");
    program.add_argument(name1, name2)
        .action([&var, fun](const auto& a) { var = fun(a); })
        .default_value(var)
        .help(hhelp);
}
template<typename T>
void myflag(const char* name, T& var, const char* hhelp) {
    static_assert(std::is_same_v<T, int>, "myflag var must be int");
    program.add_argument(name)
        .action([&var](const auto&) { var = 1; })
        .default_value(var)
        .flag()
        .help(hhelp);
}

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
int poly_nvars = -1;
int prime_field = -1;
int strip_opt_indep = 0;
FG fg = nullptr;
int num_threads = 1;
int bits_jobs = 10;

// mode
int mode = 0;
int mpfr_precision = 64;

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

    add_arg2("-v", "--verb", conf.verb, fc_int, "Verbosity");
    add_arg2("-s", "--seed", conf.seed, fc_int, "Seed");
    program.add_argument("-v", "--version") \
        .action([&](const auto&) {cout << print_version(); exit(EXIT_SUCCESS);}) \
        .flag()
        .help("Print version and exit");
    add_arg("--mode", mode , fc_int, R"delimiter(0=integer counting,
1=weighted counting over the rationals,
2=complex rational numbers,
3=multivariate polynomials over the rational field,
4=parity counting,
5=counting over a prime field (see --prime),
6=mpfr floating point complex numbers (see --mpfrprec),
7=mpfr floating point real numbers (see --mpfrprec),
8=mpfi intervals (see --mpfrprec)
9=mpqi rational/interval adaptive (see --mpfrprec)
)delimiter");
    add_arg("--prime", prime_field, fc_int, "Prime for prime field counting");
    add_arg("--npolyvars", poly_nvars, fc_int, "Number of variables in the polynomial field");
    add_arg("--delta", conf.delta, fc_double, "Delta");
    /* add_arg("--breakid", do_breakid, fc_int, "Enable BreakID"); */
    add_arg("--appmct", conf.appmc_timeout, fc_double, "after K seconds");
    add_arg("--epsilon", conf.appmc_epsilon, fc_double, "AppMC epsilon");
    add_arg("--chronobt", conf.do_chronobt, fc_int, "ChronoBT. SAT must be DISABLED or this will fail");
    add_arg("--prob", conf.do_probabilistic_hashing, fc_int, "Use probabilistic hashing. When set to 0, we are not running in probabilistic mode, but in deterministic mode, i.e. delta is 0 in Ganak mode (not in case we switch to ApproxMC mode via --appmct)");

    // Arjun options
    add_arg("--arjun", do_arjun, fc_int, "Use arjun");
    add_arg("--arjunverb", arjun_verb, fc_int, "Arjun verb");
    add_arg("--arjungates", arjun_gates, fc_int, "Use arjun's gate detection");
    add_arg("--arjunextend", etof_conf.do_extend_indep, fc_int, "Extend indep via Arjun's extend system");
    add_arg("--prebackbone", do_pre_backbone, fc_int, "Perform backbone before other things");
    add_arg("--puura", do_puura, fc_int, "Run Puura");
    add_arg("--puurabackbone", simp_conf.do_backbone_puura, fc_int, "Perform backbone in Puura");
    add_arg("--puuraautarky", etof_conf.do_autarky, fc_int, "Do autarky in Puura");
    add_arg("--arjuniter1", simp_conf.iter1, fc_int, "Arjun's iter1");
    add_arg("--arjuniter2", simp_conf.iter2, fc_int, "Arjun's iter2");
    add_arg("--arjunprobe", do_probe_based, fc_int, "Probe based arjun");
    add_arg("--arjunsimplev", arjun_simp_level, fc_int, "Arjun simp level");
    add_arg("--arjunbackwmaxc", arjun_backw_maxc, fc_int, "Arjun backw max confl");
    add_arg("--arjunoraclefindbins", arjun_oracle_find_bins, fc_int, "Arjun's oracle should find bins or not");
    add_arg("--arjunoraclemult", simp_conf.oracle_mult, fc_double, "Multiplier for Arjun's oracle timeout when it is called from Puura");
    add_arg("--bce", etof_conf.do_bce, fc_int, "Do static BCE");
    add_arg("--bveresolvmaxsz", simp_conf.bve_too_large_resolvent, fc_int, "Puura BVE max resolvent size in literals. -1 == no limit");
    add_arg("--bvegrowiter1", simp_conf.bve_grow_iter1, fc_int, "Puura BVE growth allowance iter1");
    add_arg("--bvegrowiter2", simp_conf.bve_grow_iter2, fc_int, "Puura BVE growth allowance iter2");
    add_arg("--extraoracle", simp_conf.oracle_extra, fc_int, "Extra oracle at the end of puura");
    add_arg("--resolvsub", simp_conf.do_subs_with_resolvent_clauses, fc_int, "Sets relevant CMS option: subsume other clauses with resolvent clauses");
    add_arg("--arjunoraclegetlearnt", simp_conf.oracle_vivify_get_learnts, fc_int, "Arjun's oracle should get learnts");
    add_arg("--arjundebugcnf", debug_arjun_cnf, fc_string, "Write debug arjun CNF into this file");
    add_arg("--arjuncmsmult", arjun_cms_glob_mult, fc_double,  "Pass this multiplier to CMSat through Arjun");
    add_arg("--arjunsamplcutoff", arjun_further_min_cutoff, fc_int,  "Only perform further arjun-based minimization in case the minimized indep support is larger or equal to this");
    add_arg("--arjunextendccnr", arjun_extend_ccnr, fc_int,  "Filter extend of ccnr gates via CCNR mems, in the millions");
    add_arg("--arjunweakenlim", simp_conf.weaken_limit, fc_int,  "Arjun's weaken limitation");

    // TD options
    add_arg("--td", conf.do_td, fc_int, "Run TD decompose");
    add_arg("--tdmaxw", conf.td_maxweight, fc_double, "TD max weight");
    add_arg("--tdminw", conf.td_minweight, fc_double, "TD min weight");
    add_arg("--tddiv", conf.td_divider, fc_double, "TD divider");
    add_arg("--tdexpmult", conf.td_exp_mult, fc_double, "TD exponential multiplier");
    add_arg("--tdcheckagainstind", conf.do_check_td_vs_ind, fc_int, "Check TD against indep size");
    add_arg("--tditers", conf.td_iters, fc_int, "TD flowcutter iterations (restarts)");
    add_arg("--tdsteps", conf.td_steps, fc_int, "TD flowcutter number of steps at most");
    add_arg("--tdlook", conf.td_lookahead, fc_int, "-1 means never");
    add_arg("--tdlooktwcut", conf.td_lookahead_tw_cutoff, fc_int, "TD lookahead only when TW of current comp is larger than this value");
    add_arg("--tdlookiters", conf.td_lookahead_iters, fc_int, "TD lookahead iterations");
    add_arg("--tdcontract", conf.do_td_contract, fc_int, "TD contract over opt indep set");
    add_arg("--tdlimit", conf.td_limit, fc_int, "If TD is over this, reduce weight to 0.1");
    add_arg("--tdoptindep", conf.do_td_use_opt_indep, fc_int, "Use opt indep for TD computation");
    add_arg("--tdmaxdensity", conf.td_max_density, fc_double, "Max density for TD computation");
    add_arg("--tdmaxedgeratio", conf.td_max_edge_var_ratio, fc_int, "Max edge to var ratio for TD computation");
    add_arg("--tduseadj", conf.td_do_use_adj, fc_int, "TD should use adjacency matrix for computing TD scores");
    add_arg("--tdreadfile", conf.td_read_file, fc_string, "Read TD scores from this file");
    add_arg("--tdvis", conf.td_visualize_dot_file, fc_string, "Visualize the TD into this file in DOT format");

    // Clause DB options
    add_arg("--rdbclstarget", conf.rdb_cls_target, fc_int, "RDB clauses target size (added to this are LBD 3 or lower)");
    add_arg("--rdbeveryn", conf.reduce_db_everyN, fc_int, "Reduce the clause DB every N conflicts");
    add_arg("--rdbkeepused", conf.rdb_keep_used, fc_int, "RDB keeps clauses that are used");
    add_arg("--consolidateeveryn", conf.consolidate_every_n, fc_int, "Consolidate memory after every N learnt clause");
    add_arg("--lbd", conf.base_lbd_cutoff, fc_int, "Initial LBD cutoff");
    add_arg("--updatelbdcutoff", conf.do_update_lbd_cutoff, fc_int, "Update lbd cutoff");

    // Decision options
    add_arg("--polar", conf.polar_type, fc_int, "0=standard_polarity, 1=polar cache, 2=false, 3=true");
    add_arg("--decide", conf.decide, fc_int, "ignore or not ignore TD");
    add_arg("--initact", conf.do_init_activity_scores, fc_int, "Init activity scores to var freq");
    add_arg("--vsadsadjust", conf.vsads_readjust_every, fc_int, "VSADS ajust activity every N");
    add_arg("--actscorediv", conf.act_score_divisor, fc_double, "Activity score divisor");
    add_arg("--freqscorediv", conf.freq_score_divisor, fc_double, "Component frequency score divisor");

    // Cache options
    add_arg("--cache", conf.do_use_cache, fc_int, "Use (i.e. store and retrieve) cache");
    add_arg("--maxcache", conf.maximum_cache_size_MB, fc_int, "Max cache size in MB");
    add_arg("--cachetime", conf.cache_time_update, fc_int, "2 = set to mid-point");
    add_arg("--lru", conf.lru_eviction, fc_int, "Cache eviction: 1 = LRU (evict oldest), 0 = reverse-LRU/default (evict newest)");

    // BuDDy options
    add_arg("--buddy", conf.do_buddy, fc_int, "Run BuDDy");
    add_arg("--buddymaxcls", conf.buddy_max_cls, fc_int, "Run BuDDy");

    // Vivif options -- inprocessing during Ganak
    add_arg("--vivif", conf.do_vivify, fc_int, "Vivify clauses");
    add_arg("--vivifevery", conf.vivif_every, fc_int, "Vivify every N conflicts");
    add_arg("--vivifmult", conf.vivif_mult, fc_double, "How much to multiply timeout for vivif");
    add_arg("--vivifoutern", conf.vivif_outer_every_n, fc_int, "How many restarts between outer vivif");
    add_arg("--totusedcutoffvivif", conf.tot_used_cutoff_vivif, fc_int, "Total used vivif cutoff");

    // SBVA options
    add_arg("--sbvasteps", etof_conf.num_sbva_steps, fc_int, "SBVA steps. 0 = no SBVA");
    add_arg("--sbvaclcut", etof_conf.sbva_cls_cutoff, fc_int, "SBVA cls cutoff");
    add_arg("--sbvalitcut", etof_conf.sbva_lits_cutoff, fc_int, "SBVA lits cutoff");
    add_arg("--sbvabreak", etof_conf.sbva_tiebreak, fc_int, "1 = sbva");

    // SAT solver options
    add_arg("--satsolver", conf.do_use_sat_solver, fc_int, "Use SAT solver when all minimal indep set has been set");
    add_arg("--satrst", conf.do_sat_restart, fc_int, "Inside SAT solver, perform restarts");
    add_arg("--satrstmult", conf.sat_restart_mult, fc_int, "SAT restart multiplier");
    add_arg("--satpolarcache", conf.do_sat_polar_cache, fc_int, "Inside SAT solver, use polarity cache");
    add_arg("--satvsids", conf.do_sat_vsids, fc_int, "Inside SAT solver, use VSIDS, not VSADS");

    // Opt independent set options
    add_arg("--allindep", etof_conf.all_indep, fc_int, "All variables can be made part of the indepedent support. Indep support is given ONLY to help the solver.");
    add_arg("--arjunextendmaxconfl", arjun_extend_max_confl, fc_int, "Max number of conflicts per extend operation in Arjun");
    add_arg("--arjunextend", etof_conf.do_extend_indep, fc_int, "Max number of conflicts per extend operation in Arjun");
    add_arg("--stripoptindep", strip_opt_indep, fc_int, "Strip optional indep support");

    // Analyze candidates options
    add_arg("--analyzecand", conf.analyze_cand_update, fc_int, "Update analyze candidates if more than N vars are still undecided from opt indep set");

    // Restart options
    add_arg("--rstfirst", conf.first_restart, fc_int, "Run restarts");
    add_arg("--restart", conf.do_restart, fc_int, "Run restarts");
    add_arg("--rsttype", conf.restart_type, fc_int, "Check count at every step");
    add_arg("--rstcutoff", conf.restart_cutoff_mult, fc_double, "Multiply cutoff with this");
    add_arg("--rstcheckcnt", conf.do_cube_check_count, fc_int, "Check the count of each cube. 1 = use ganak itself without restart. 2 = use CMS one-by-one counting");
    add_arg("--rstreadjust", conf.do_readjust_for_restart, fc_int, "Readjust params for restart");
    add_arg("--maxrst", conf.max_num_rst, fc_int, "Max number of restarts");
    add_arg("--maxcubesperrst", conf.max_num_cubes_per_restart, fc_int,  "Max number of cubes per restart");

    // Multi-threading options
    add_arg("--threads", num_threads, fc_int, "Number of threads to use. -1 = all available cores");
    add_arg("--bitsjobs", bits_jobs, fc_int, "Number of variables to multi-thread on (8 = 256 jobs)");
    program.add_argument("inputfile").remaining().help("input CNF");

    // Minor options
    add_arg("--mpfrprec", mpfr_precision, fc_int, "MPFR precision in bits");
}

void parse_supported_options(int argc, char** argv) {
    add_ganak_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout << "Flexible Weighted Model Counter" << endl << endl
            << "ganak [options] inputfile" << endl;
            cout << program << endl;
            exit(EXIT_SUCCESS);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    if (conf.do_use_sat_solver && !conf.do_chronobt) {
      cout << "ERROR: When chronobt is disabled, SAT solver cannot be used" << endl;
      exit(EXIT_FAILURE);
    }
    if (bits_jobs < 0 || bits_jobs > 20) {
      cout << "ERROR: bitsjobs must be between 0 and 20, inclusive" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads < -1) {
      cout << "ERROR: number of threads must not be less than -1" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads > 1024) {
      cout << "ERROR: number of threads must not be more than 1024" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads == 0) {
      cout << "ERROR: number of threads must not be 0" << endl;
      exit(EXIT_FAILURE);
    }
    if (mpfr_precision < 2) {
      cout << "ERROR: mpfr precision must be at least 2 bits" << endl;
      exit(EXIT_FAILURE);
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
  cnf.check_cnf_sampl_sanity();
  counter.new_vars(cnf.nVars());

  set<uint32_t> tmp;
  for(auto const& s: cnf.get_sampl_vars()) tmp.insert(s+1);
  counter.set_indep_support(tmp);
  if (cnf.get_opt_sampl_vars_set()) {
    tmp.clear();
    for(auto const& s: cnf.get_opt_sampl_vars()) tmp.insert(s+1);
  }
  counter.set_optional_indep_support(tmp);
  if (conf.verb) counter.print_indep_distrib();

  if (cnf.get_weighted()) {
    for(const auto& t: cnf.get_weights()) {
      counter.set_lit_weight(Lit(t.first+1, true), t.second.pos);
      counter.set_lit_weight(Lit(t.first+1, false), t.second.neg);
    }
  }

  for(const auto& cl: cnf.get_clauses()) counter.add_irred_cl(cms_to_ganak_cl(cl));
  for(const auto& cl: cnf.get_red_clauses()) counter.add_red_cl(cms_to_ganak_cl(cl));
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

double digit_precision_mpfi(mpfi_srcptr v) {
    mpfr_prec_t prec = mpfi_get_prec(v);
    mpfr_t left;
    mpfr_init2(left, prec);
    mpfr_t right;
    mpfr_init2(right, prec);

    mpfi_get_left(left, v);
    mpfi_get_right(right, v);
    if (mpfr_sgn(left) != mpfr_sgn(right)) {
        mpfr_clear(left);
        mpfr_clear(right);
        return 0.0;
    }

    mpfr_t diam;
    mpfr_init2(diam, prec);
    mpfi_diam_rel(diam, v);
    if (mpfr_sgn(diam) == 0) {
        mpfr_clear(diam);
        mpfr_clear(left);
        mpfr_clear(right);
        return max_digit_precision;
    }

    mpfr_log10(diam, diam, MPFR_RNDN);
    double result = -mpfr_get_d(diam, MPFR_RNDN);
    if (result < 0)
        result = 0.0;

    if (result > max_digit_precision)
        result = max_digit_precision;

    mpfr_clear(diam);
    mpfr_clear(left);
    mpfr_clear(right);
    return result;
}

void print_log(const mpfi_t& val, string extra = "") {
    mpfr_t left, right;
    mpfr_init2(left, mpfr_precision);
    mpfr_init2(right, mpfr_precision);
    mpfi_get_left(left, val);
    mpfi_get_right(right, val);
    print_log(left, extra + " left bound");
    print_log(right, extra + " right bound");
    mpfr_clear(left);
    mpfr_clear(right);
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
// result must be initialized with mpfr_init2 before calling
void compute_collision_prob(mpfr_t result, const uint64_t lookups, uint64_t elems) {
    mpfr_t lookups2;
    mpfr_init2(lookups2, 256);
    mpfr_set_ui(lookups2, lookups, MPFR_RNDN);
    mpfr_log2(lookups2, lookups2, MPFR_RNDN);

    mpfr_t elems2;
    mpfr_init2(elems2, 256);
    mpfr_set_ui(elems2, elems, MPFR_RNDN);
    mpfr_log2(elems2, elems2, MPFR_RNDN);

    mpfr_t e;
    mpfr_init2(e, 256);
    mpfr_set_si(e, -64, MPFR_RNDN);
    mpfr_add(e, lookups2, e, MPFR_RNDN);
    mpfr_add(e, elems2, e, MPFR_RNDN);
    // e = log2(lookups) + log2(elems) - 64

    // Compute 2^e
    mpfr_exp2(result, e, MPFR_RNDN);

    // Clear temporary variables
    mpfr_clear(lookups2);
    mpfr_clear(elems2);
    mpfr_clear(e);
}

void run_weighted_counter(Ganak& counter, const ArjunNS::SimplifiedCNF& cnf, const double start_time) {
    FF cnt = cnf.get_multiplier_weight()->dup();
    if (!cnf.get_multiplier_weight()->is_zero()) *cnt *= *counter.count(bits_jobs, num_threads);
    cout << "c o Total time [Arjun+GANAK]: " << setprecision(2)
        << std::fixed << (cpu_time() - start_time) << endl;

    string out = "c o type ";
    if (cnf.get_projected()) out+="p";
    if (cnf.get_weighted()) out += "wmc";
    else out += "mc";

    if (!cnt->is_zero()) cout << "s SATISFIABLE" << endl;
    else cout << "s UNSATISFIABLE" << endl;
    if (mode == 0 || mode == 1 || mode == 2 || mode == 6 || mode == 7 || mode == 8 || mode == 9) {
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
      } else if (mode == 8) {
        // MPFR intervals
        if (cnf.get_projected()) cout << "c s type pwmc" << endl;
        else cout << "c s type wmc" << endl;
        const FMpfi* od = dynamic_cast<const FMpfi*>(ptr);
        assert(od != nullptr);
        print_log(od->val);
        cout << "c s exact quadruple float interval " << *od << endl;
        cout << "c s digit precision of interval: " << digit_precision_mpfi(od->val) << endl;
      } else if (mode == 9) {
        // mpqi rational/interval adaptive
        if (cnf.get_projected()) cout << "c s type pwmc" << endl;
        else cout << "c s type wmc" << endl;
        const FMpqi* od = dynamic_cast<const FMpqi*>(ptr);
        cout << "c s exact arb frac " << *od << endl;
        cout << "c s digit precision: " << digit_precision_mpqi(const_cast<mpqi_ptr>(&od->val)) << endl;
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
      mpfr_init2(collision_prob, 256);
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
        fg = std::make_unique<ArjunNS::FGenMpfr>(mpfr_precision);
        break;
    case 8:
        fg = std::make_unique<FGenMpfi>(mpfr_precision);
        break;
    case 9:
        fg = std::make_unique<FGenMpqi>(mpfr_precision);
        break;
    case 2:
        fg = std::make_unique<FGenComplex>();
        break;
    case 6:
        fg = std::make_unique<FGenMPFComplex>(mpfr_precision);
        break;
    case 3:
        if (poly_nvars == -1) {
          cout << "c o [arjun] ERROR: Must provide number of polynomial vars for mode 3 via --npolyvars" << endl;
          exit(EXIT_FAILURE);
        }
        fg = std::make_unique<FGenPoly>(poly_nvars);
        break;
    case 4:
        fg = std::make_unique<FGenParity>();
        break;
    case 5:
        if (prime_field == -1) {
          cout << "c o [arjun] ERROR: Must provide prime field for mode 5 via --prime" << endl;
          exit(EXIT_FAILURE);
        }
        fg = std::make_unique<FGenPrime>(prime_field);
        break;
    default:
        cout << "c o [arjun] ERROR: Unknown mode" << endl;
        exit(EXIT_FAILURE);
  }
  ArjunNS::SimplifiedCNF cnf(fg);

  // Parse the CNF
  if (!program.is_used("inputfile")) read_in_a_file("-",  &cnf, etof_conf.all_indep, fg);
  else {
    auto files = program.get<std::vector<std::string>>("inputfile");
    if (files.empty()) {
      cout << "ERROR: you provided --inputfile but no file. Strange. Exiting. " << endl;
      exit(EXIT_FAILURE);
    } else if (files.size() == 1) {
      const string& fname = files[0];
      read_in_a_file(fname, &cnf, etof_conf.all_indep, fg);
    } else {
        cout << "ERROR: you must only give one CNF as input (or none, and then we read from STDIN)" << endl;
        cout << "       You provided " << files.size() << " files: ";
        for (const auto& f: files) cout << f << " ";
        cout << endl;
        exit(EXIT_FAILURE);
    }
  }

  if (cnf.get_weighted() && conf.do_buddy) {
    cout << "ERROR: Cannot run BuDDy with weighted CNF" << endl;
    exit(EXIT_FAILURE);
  }
  cnf.clean_idiotic_mccomp_weights();
  cnf.check_cnf_sampl_sanity();
  cnf.check_cnf_vars();
  verb_print(1, "CNF projection set size: " << cnf.get_sampl_vars().size());

  // Run Arjun
  if (!do_arjun) cnf.renumber_sampling_vars_for_ganak();
  else run_arjun(cnf);
  cnf.remove_equiv_weights();
  if (strip_opt_indep) cnf.strip_opt_sampling_vars();
  if (conf.verb >= 2) {
    cout << "c o sampl_vars: "; print_vars(cnf.get_sampl_vars()); cout << endl;
    if (cnf.get_opt_sampl_vars_set()) {
      cout << "c o opt sampl_vars: "; print_vars(cnf.get_opt_sampl_vars()); cout << endl;
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
