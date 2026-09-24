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
#include <iomanip>
#include <gmpxx.h>
#include <mpfr.h>
#include <charconv>
#include <type_traits>
/* #include <breakid/breakid.hpp> */
#include <arjun/arjun.h>
#include "src/argparse.hpp"
#include "mpoly.hpp"
#include "mparity.hpp"
#include <approxmc/approxmc.h>
#include <treedecomp/treedecomp_version.hpp>
#include "file_read_helper.h"

static constexpr uint32_t max_digit_precision = 1e6;

using namespace GanakInt;
using std::setprecision;

#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif

template<class T> static T parse_opt(const std::string& s) {
    if constexpr (std::is_same_v<T, std::string>) {
        return s;
    } else if constexpr (std::is_same_v<T, bool>) {
        return parse_opt<int>(s) != 0;
    } else if constexpr (std::is_floating_point_v<T>) {
        size_t pos = 0;
        double val;
        try { val = std::stod(s, &pos); }
        catch (const std::exception&) { throw std::invalid_argument("not a number: " + s); }
        if (pos != s.size()) throw std::invalid_argument("trailing characters in number: " + s);
        return val;
    } else if constexpr (std::is_integral_v<T>) {
        T val{};
        auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), val);
        if (ec == std::errc::result_out_of_range) throw std::invalid_argument("integer out of range: " + s);
        if (ec != std::errc{}) throw std::invalid_argument("not an integer: " + s);
        if (ptr != s.data() + s.size()) throw std::invalid_argument("trailing characters in integer: " + s);
        return val;
    } else {
        static_assert(sizeof(T) == 0, "parse_opt: unsupported option type");
    }
}

using std::string;
using std::vector;
argparse::ArgumentParser program = argparse::ArgumentParser("ganak",
        GANAK::get_version_sha1(),
        argparse::default_arguments::help);

template<typename T>
void add_arg(const char* name, T& var, const char* hhelp) {
    program.add_argument(name)
        .action([&var](const std::string& a) { var = parse_opt<T>(a); })
        .default_value(var)
        .help(hhelp);
}
template<typename T>
void add_arg2(const char* name1, const char* name2, T& var, const char* hhelp) {
    program.add_argument(name1, name2)
        .action([&var](const std::string& a) { var = parse_opt<T>(a); })
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
int arjun_extend_max_confl = 30000;
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
bool disconnected_allowed = false;
uint32_t arjun_further_min_cutoff = 10;
int arjun_extend_ccnr = 0;
int poly_nvars = -1;
int prime_field = -1;
int strip_opt_indep = 0;
FG fg = nullptr;

// threads
int num_threads = 1;
int bits_jobs = 10;
int debug_threads = 0;

// mode
int mode = 0;
int mpfr_precision = 128;

string print_version()
{
    std::stringstream ss;
    ss << "c o Ganak SHA1: " << GANAK::get_version_sha1() << endl;
    ss << "c o Arjun SHA1: " << ArjunNS::Arjun::get_version_sha1() << endl;
    ss << "c o SBVA SHA1: " << ArjunNS::Arjun::get_sbva_version_sha1() << endl;
    ss << "c o CMS SHA1: " << CMSat::SATSolver::get_version_sha1() << endl;
    ss << "c o CaDiCaL SHA1: " << CMSat::SATSolver::get_cadical_version_sha1() << endl;
    ss << "c o CadiBack SHA1: " << CMSat::SATSolver::get_cadiback_version_sha1() << endl;
    ss << "c o ApproxMC SHA1: " << ApproxMC::AppMC::get_version_sha1() << endl;
    ss << "c o TreeDecomp SHA1: " << TWD::get_version_sha1() << endl;
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

    add_arg2("-v", "--verb", conf.verb, "Verbosity");
    add_arg2("-s", "--seed", conf.seed, "Seed");
    program.add_argument("-v", "--version")
        .action([&](const auto&) {cout << print_version(); exit(EXIT_SUCCESS);})
        .flag()
        .help("Print version and exit");
    add_arg("--mode", mode , R"delimiter(0=integer counting,
1=weighted counting over the rationals,
2=complex rational numbers,
3=multivariate polynomials over the rational field,
4=parity counting,
5=counting over a prime field (see --prime),
6=mpfr floating point complex numbers (see --mpfrprec),
7=mpfr floating point real numbers (see --mpfrprec),
13=multivariate Laurent polynomials over the rational field (see --npolyvars)
)delimiter");
    add_arg("--prime", prime_field, "Prime for prime field counting");
    add_arg("--npolyvars", poly_nvars, "Number of variables in the polynomial field");
    add_arg("--delta", conf.delta, "Delta");
    /* add_arg("--breakid", do_breakid, "Enable BreakID"); */
    add_arg("--appmct", conf.appmc_timeout, "after K seconds");
    add_arg("--epsilon", conf.appmc_epsilon, "AppMC epsilon");
    add_arg("--chronobt", conf.do_chronobt, "ChronoBT. SAT must be DISABLED or this will fail");
    add_arg("--prob", conf.do_probabilistic_hashing, "Use probabilistic hashing. When set to 0, we are not running in probabilistic mode, but in deterministic mode, i.e. delta is 0 in Ganak mode (not in case we switch to ApproxMC mode via --appmct)");
    program.add_argument("--fast")
        .action([&](const auto&) {
          arjun_cms_glob_mult = 0.1;
          simp_conf.oracle_mult = 0.1;
          arjun_backw_maxc = 50;
        })
        .flag()
        .help("Optimize for quick/easy instances (<5 mins)");

    // d-DNNF compilation
    add_arg("--compile", conf.compile_fname, "Compile the search trace into a (Decision-)d-DNNF circuit and write it to this file (d4 .nnf format). Forces a clean single-threaded search (no restarts, exact cache, no BuDDy/vivify, no Arjun/Puura). SAT oracle stays on (witnesses synthesized vars on projected inputs).");

    // Arjun options
    add_arg("--arjun", do_arjun, "Use arjun");
    add_arg("--arjunverb", arjun_verb, "Arjun verb");
    add_arg("--arjungates", arjun_gates, "Use arjun's gate detection");
    add_arg("--arjunextend", etof_conf.do_extend_indep, "Extend indep via Arjun's extend system");
    add_arg("--prebackbone", do_pre_backbone, "Perform backbone before other things");
    add_arg("--puura", do_puura, "Run Puura");
    add_arg("--puurabackbone", simp_conf.do_backbone_puura, "Perform backbone in Puura");
    add_arg("--puurabackbonemaxconfl", simp_conf.backbone_max_confl, "Max conflicts for backbone in Puura (-1 = unlimited)");
    add_arg("--puuraautarky", etof_conf.do_autarky, "Do autarky in Puura");
    add_arg("--arjuniter1", simp_conf.iter1, "Arjun's iter1");
    add_arg("--arjuniter2", simp_conf.iter2, "Arjun's iter2");
    add_arg("--arjunprobe", do_probe_based, "Probe based arjun");
    add_arg("--arjunsimplev", arjun_simp_level, "Arjun simp level");
    add_arg("--arjunbackwmaxc", arjun_backw_maxc, "Arjun backw max confl");
    add_arg("--arjunoraclefindbins", arjun_oracle_find_bins, "Arjun's oracle should find bins or not");
    add_arg("--arjunoraclemult", simp_conf.oracle_mult, "Multiplier for Arjun's oracle timeout when it is called from Puura");
    add_arg("--puuraoraclevivif", simp_conf.oracle_vivify, "Run Puura's main oracle vivification pass");
    add_arg("--puuraoraclesparsify", simp_conf.oracle_sparsify, "Run Puura's main oracle sparsification pass");
    add_arg("--puurabve", simp_conf.do_bve, "Run BVE in Puura");
    add_arg("--bveresolvmaxsz", simp_conf.bve_too_large_resolvent, "Puura BVE max resolvent size in literals. -1 == no limit");
    add_arg("--bveresolvmaxsz2", simp_conf.bve_too_large_resolvent2, "Like --bveresolvmaxsz, for the 2nd elim pass");
    add_arg("--xorgatemaxsize", simp_conf.xor_gate_find_maxsize, "Max clause size for XOR-gate finding");
    add_arg("--bvegrowiter1", simp_conf.bve_grow_iter1, "Puura BVE growth allowance iter1");
    add_arg("--iter2grow", simp_conf.bve_grow_iter2, "Puura BVE growth allowance iter2");
    add_arg("--iter2growlarge", simp_conf.bve_grow_iter2_large, "If >= 0: used instead of --iter2grow when more than --iter2growlargevars vars are left before iter2");
    add_arg("--iter2growlargevars", simp_conf.bve_grow_iter2_large_vars, "Vars-left threshold for --iter2growlarge");
    add_arg("--bveocclim", simp_conf.bve_occ_cutoff, "BVE: refuse a var whose more frequent polarity occurs more than this often (CaDiCaL's elimocclim). 0 = no such limit");
    add_arg("--bveclsmaxsz", simp_conf.bve_cls_max_size, "BVE: refuse a var that occurs in a clause longer than this. 0 = no limit");
    add_arg("--distillremlevel", simp_conf.distill_rem_level, "Clause removal during Puura's distillation. 0 = never, 1 = only on a real conflict, 2 = also when a literal is positively implied. Levels below 2 keep gate clauses that BVE needs to recover definitions");
    add_arg("--extraoracle", simp_conf.oracle_extra, "Extra oracle at the end of puura");
    add_arg("--resolvsub", simp_conf.do_subs_with_resolvent_clauses, "Sets relevant CMS option: subsume other clauses with resolvent clauses");
    add_arg("--arjunoraclegetlearnt", simp_conf.oracle_vivify_get_learnts, "Arjun's oracle should get learnts");
    add_arg("--arjundebugcnf", debug_arjun_cnf, "Write debug arjun CNF into this file");
    add_arg("--arjuncmsmult", arjun_cms_glob_mult,  "Pass this multiplier to CMSat through Arjun");
    add_arg("--arjunsamplcutoff", arjun_further_min_cutoff,  "Only perform further arjun-based minimization in case the minimized indep support is larger or equal to this");
    add_arg("--arjunextendccnr", arjun_extend_ccnr,  "Filter extend of ccnr gates via CCNR mems, in the millions");
    add_arg("--arjunweakenlim", simp_conf.weaken_limit,  "Arjun's weaken limitation");

    // TD options
    add_arg("--td", conf.do_td, "Run TD decompose");
    add_arg("--tdmaxw", conf.td_maxweight, "TD max weight");
    add_arg("--tdminw", conf.td_minweight, "TD min weight");
    add_arg("--tddiv", conf.td_divider, "TD divider");
    add_arg("--tdexpmult", conf.td_exp_mult, "TD exponential multiplier");
    add_arg("--tditers", conf.td_iters, "TD flowcutter iterations (restarts)");
    add_arg("--tdsteps", conf.td_steps, "TD flowcutter number of steps at most");
    add_arg("--tdbandpct", conf.td_band_pct, "TD: a candidate up to this % wider than the narrowest TD seen can still win, by splitting better");
    add_arg("--tdflatpct", conf.td_flat_pct, "TD: if the TD's width is at least this % of the graph's nodes, the graph is too dense for the TD to say anything, and it does not guide the branching. 0 = off");
    add_arg("--tdsepwpct", conf.td_sep_weight_pct, "TD: within one level of the TD, prefer vars that do more separating work (small adhesion in front of a large subtree). In % of one TD level. 0 = off");
    add_arg("--tdsplitwpct", conf.td_split_weight_pct, "TD: scale the TD branching weight by how well the TD splits the graph, by up to this %. 0 = off");
    add_arg("--tddensepct", conf.td_dense_pct, "TD: the split only decides when the width is over this % of the graph's nodes, below it the width alone does. 100 = never");
    add_arg("--tdlook", conf.td_lookahead, "-1 means never");
    add_arg("--tdlooktwcut", conf.td_lookahead_tw_cutoff, "TD lookahead only when TW of current comp is larger than this value");
    add_arg("--tdlookiters", conf.td_lookahead_iters, "TD lookahead iterations");
    add_arg("--tdcontract", conf.do_td_contract, "TD contract over opt indep set");
    add_arg("--tdlimit", conf.td_limit, "If TD is over this, reduce weight to 0.1");
    add_arg("--tdoptindep", conf.do_td_use_opt_indep, "Use opt indep for TD computation");
    add_arg("--tdmaxdensity", conf.td_max_density, "Max density for TD computation");
    add_arg("--tdmaxedgeratio", conf.td_max_edge_var_ratio, "Max edge to var ratio for TD computation");
    add_arg("--tduseadj", conf.td_do_use_adj, "TD should use adjacency matrix for computing TD scores");
    add_arg("--tdreadfile", conf.td_read_file, "Read TD scores from this file");
    add_arg("--tdvis", conf.td_visualize_dot_file, "Visualize the TD into this file in DOT format");
    add_arg("--tddumpcnf", conf.td_dump_cnf_file, "Dump the CNF used to build the primal graph for TD computation, one DIMACS file per component: FILE.0, FILE.1, ...");

    // Clause DB options
    add_arg("--rdbclstarget", conf.rdb_cls_target, "RDB clauses target size (added to this are LBD 3 or lower)");
    add_arg("--rdbeveryn", conf.reduce_db_everyN, "Reduce the clause DB every N conflicts");
    add_arg("--rdbkeepused", conf.rdb_keep_used, "RDB keeps clauses that are used");
    add_arg("--consolidateeveryn", conf.consolidate_every_n, "Consolidate memory after every N learnt clause");
    add_arg("--lbd", conf.base_lbd_cutoff, "Initial LBD cutoff");
    add_arg("--updatelbdcutoff", conf.do_update_lbd_cutoff, "Update lbd cutoff");

    // Decision options
    add_arg("--polar", conf.polar_type, "0=standard_polarity, 1=polar cache, 2=false, 3=true");
    add_arg("--decide", conf.decide, "ignore or not ignore TD");
    add_arg("--initact", conf.do_init_activity_scores, "Init activity scores to var freq");
    add_arg("--vsadsadjust", conf.vsads_readjust_every, "VSADS ajust activity every N");
    add_arg("--actscorediv", conf.act_score_divisor, "Activity score divisor");
    add_arg("--freqscorediv", conf.freq_score_divisor, "Component frequency score divisor");

    // Cache options
    add_arg("--cache", conf.do_use_cache, "Use (i.e. store and retrieve) cache");
    add_arg("--maxcache", conf.maximum_cache_size_MB, "Max cache size in MB");
    add_arg("--cachetime", conf.cache_time_update, "2 = set to mid-point");
    add_arg("--lru", conf.lru_eviction, "Cache eviction: 1 = LRU (evict oldest), 0 = reverse-LRU/default (evict newest)");

    // BuDDy options
    add_arg("--buddy", conf.do_buddy, "Run BuDDy");
    add_arg("--buddymaxcls", conf.buddy_max_cls, "Run BuDDy");

    // Shrinking options
    add_arg("--shrink", conf.do_shrink, "Block-wise secondary UIP shrinking (CaDiCaL-style)");
    add_arg("--bumpreason", conf.do_bump_reason, "Bump reason clause literals during conflict analysis (CaDiCaL-style)");

    // Vivif options -- inprocessing during Ganak
    add_arg("--vivif", conf.do_vivify, "Vivify clauses");
    add_arg("--vivifevery", conf.vivif_every, "Vivify every N conflicts");
    add_arg("--vivifmult", conf.vivif_mult, "How much to multiply timeout for vivif");
    add_arg("--vivifoutern", conf.vivif_outer_every_n, "How many restarts between outer vivif");
    add_arg("--totusedcutoffvivif", conf.tot_used_cutoff_vivif, "Total used vivif cutoff");

    // SBVA options
    add_arg("--sbvasteps", etof_conf.num_sbva_steps, "SBVA steps. 0 = no SBVA");
    add_arg("--sbvaclcut", etof_conf.sbva_cls_cutoff, "SBVA cls cutoff");
    add_arg("--sbvalitcut", etof_conf.sbva_lits_cutoff, "SBVA lits cutoff");
    add_arg("--sbvabreak", etof_conf.sbva_tiebreak, "1 = sbva");
    add_arg("--sbvamaxnewvars", etof_conf.sbva_max_new_vars, "Max number of new variables SBVA may add. 0 = no limit");

    // SAT solver options
    add_arg("--satsolver", conf.do_use_sat_solver, "Use SAT solver when all minimal indep set has been set");
    add_arg("--satrst", conf.do_sat_restart, "Inside SAT solver, perform restarts");
    add_arg("--satrstmult", conf.sat_restart_mult, "SAT restart multiplier");
    add_arg("--satpolarcache", conf.do_sat_polar_cache, "Inside SAT solver, use polarity cache");
    add_arg("--satvsids", conf.do_sat_vsids, "Inside SAT solver, use VSIDS, not VSADS");

    // Opt independent set options
    add_arg("--allindep", etof_conf.all_indep, "All variables can be made part of the indepedent support. Indep support is given ONLY to help the solver.");
    add_arg("--arjunextendmaxconfl", arjun_extend_max_confl, "Max number of conflicts per extend operation in Arjun");
    add_arg("--arjunextend", etof_conf.do_extend_indep, "Max number of conflicts per extend operation in Arjun");
    add_arg("--stripoptindep", strip_opt_indep, "Strip optional indep support");

    // Analyze candidates options
    add_arg("--analyzecand", conf.analyze_cand_update, "Update analyze candidates if more than N vars are still undecided from opt indep set");

    // Restart options
    add_arg("--rstfirst", conf.first_restart, "Run restarts");
    add_arg("--restart", conf.do_restart, "Run restarts");
    add_arg("--rsttype", conf.restart_type, "Check count at every step");
    add_arg("--rstcheckcnt", conf.do_cube_check_count, "Check the count of each cube. 1 = use ganak itself without restart. 2 = use CMS one-by-one counting");
    add_arg("--rstreadjust", conf.do_readjust_for_restart, "Readjust params for restart");
    add_arg("--maxrst", conf.max_num_rst, "Max number of restarts");
    add_arg("--maxcubesperrst", conf.max_num_cubes_per_restart,  "Max number of cubes per restart");
    add_arg("--extendcubes", conf.do_extend_cubes,  "Extend cubes");
    add_arg("--cuberesolve", conf.do_cube_resolve, "Resolution-merge cubes that differ in exactly one blocking literal");
    add_arg("--cubeflp", conf.do_cube_flp, "Failed-literal probing to remove forced blocking literals from cubes");
    add_arg("--smallcubedisable", conf.do_small_cube_disable, "Disable cubes beyond max-num-cubes-per-restart (sorted by LBD)");
    add_arg("--tdwrstdecay", conf.td_weight_restart_decay, "Multiply td_weight by this after each restart (1.0=no decay, 0.5=halve)");

    // Multi-threading options
    add_arg("--threads", num_threads, "Number of threads to use. -1 = all available cores");
    add_arg("--bitsjobs", bits_jobs, "Number of variables to multi-thread on (8 = 256 jobs)");
    add_arg("--debugthreads", debug_threads, "Debug threads -- use thread system, even though only one thread is allowed");
    program.add_argument("inputfile").remaining().help("input CNF");

    // Minor options
    add_arg("--mpfrprec", mpfr_precision, "MPFR precision in bits");
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
        std::string msg = err.what();
        if (msg == "Duplicate argument") {
            std::map<std::string, int> seen;
            for (int i = 1; i < argc; i++)
                if (argv[i][0] == '-') seen[argv[i]]++;
            for (auto& [k, v] : seen)
                if (v > 1) msg += ": " + k;
        }
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!conf.compile_fname.empty()) {
      // d-DNNF needs a single clean DPLL tree: no restarts, exact cache, no
      // opaque leaves (BuDDy), no clause rewriting (vivify), original var
      // numbering (Arjun/Puura off), one thread. SAT oracle stays on; it
      // witnesses synthesized vars on projected inputs and never fires otherwise.
      conf.do_restart = 0;
      conf.do_probabilistic_hashing = 0;
      conf.do_buddy = 0;
      conf.do_vivify = 0;
      do_arjun = 0;
      do_puura = 0;
      num_threads = 1;
      cout << "c o [compile] d-DNNF compilation mode -> " << conf.compile_fname << endl;
    }
    if (conf.do_use_sat_solver && !conf.do_chronobt) {
      cerr << "ERROR: When chronobt is disabled, SAT solver cannot be used" << endl;
      exit(EXIT_FAILURE);
    }
    if (bits_jobs < 0 || bits_jobs > 20) {
      cerr << "ERROR: bitsjobs must be between 0 and 20, inclusive" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads < -1) {
      cerr << "ERROR: number of threads must not be less than -1" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads > 1024) {
      cerr << "ERROR: number of threads must not be more than 1024" << endl;
      exit(EXIT_FAILURE);
    }
    if (num_threads == 0) {
      cerr << "ERROR: number of threads must not be 0" << endl;
      exit(EXIT_FAILURE);
    }
    if (debug_threads > 0 && num_threads > 1) {
      cerr << "ERROR: threads cannot be debugged when num_threads is more than 1" << endl;
      exit(EXIT_FAILURE);
    }
    if (mpfr_precision > 256) {
      cerr << "ERROR: mpfr precision must not be more than 256 bits" << endl;
      exit(EXIT_FAILURE);
    }
    if (mpfr_precision < 2) {
      cerr << "ERROR: mpfr precision must be at least 2 bits" << endl;
      exit(EXIT_FAILURE);
    }
}

void print_vars(vector<uint32_t> vars) {
  std::sort(vars.begin(), vars.end());
  for(const auto& v: vars) cout << v+1 << " ";
}

void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  double const my_time = cpu_time();
  uint64_t lits_before = 0;
  for(const auto& cl: cnf.get_clauses()) lits_before += cl.size();
  ArjunNS::Arjun arjun;
  ArjunNS::Arjun::InterpConf iconf;
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
  arjun.standalone_minimize_indep(cnf, iconf, etof_conf.all_indep);
  arjun.set_extend_ccnr(arjun_extend_ccnr);
  if (cnf.get_sampl_vars().size() >= arjun_further_min_cutoff && do_puura) {
    arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf, iconf);
  } else {
    disconnected_allowed = true;
    verb_print(1, "WARNING. Not performing puura.  "
        << "Number of sampling variables: " << cnf.get_sampl_vars().size()
        << " vs " << arjun_further_min_cutoff
        << " and  --puura is: " << do_puura
        << " components may be disconnected, which will interfere with proper TD weight calculation");
    cnf.renumber_sampling_vars_for_ganak();
  }
  // Preprocessing that hands the counter a much bigger formula is a bug
  uint64_t lits_after = 0;
  for(const auto& cl: cnf.get_clauses()) lits_after += cl.size();
  assert(lits_after <= 10 * lits_before + 100000);
  (void)lits_before; (void)lits_after;
  verb_print(1, "Arjun T: " << (cpu_time()-my_time));
}

const char* mpfr_prec_name(const int prec) {
    if (prec <= 16) return "half float";
    if (prec <= 32) return "single float";
    if (prec <= 64) return "double float";
    if (prec <= 128) return "quadruple float";
    return "octuple float";
}

string print_mpq_as_scientific(const mpq_class& number) {
    mpf_class const mpf_value(number);
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

// compute collision probability, i.e. 2^(log2(lookups) + log2(elems) - 128)
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
    mpfr_set_si(e, -128, MPFR_RNDN);
    mpfr_add(e, lookups2, e, MPFR_RNDN);
    mpfr_add(e, elems2, e, MPFR_RNDN);
    // e = log2(lookups) + log2(elems) - 128

    // Compute 2^e
    mpfr_exp2(result, e, MPFR_RNDN);

    // Clear temporary variables
    mpfr_clear(lookups2);
    mpfr_clear(elems2);
    mpfr_clear(e);
}

void run_weighted_counter(Ganak& counter, const ArjunNS::SimplifiedCNF& cnf, const double start_time) {
    FF cnt = cnf.get_multiplier_weight()->dup();
    if (!cnt->is_zero()) *cnt *= *counter.count(bits_jobs, num_threads, debug_threads);
    cout << "c o Total time [Arjun+GANAK]: " << setprecision(2)
        << std::fixed << (cpu_time() - start_time) << endl;

    string out = "c o type ";
    if (cnf.get_projected()) out+="p";
    if (cnf.get_weighted()) out += "wmc";
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
        const ArjunNS::FComplex* od = dynamic_cast<const ArjunNS::FComplex*>(ptr);
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
        const ArjunNS::MPFComplex* od = dynamic_cast<const ArjunNS::MPFComplex*>(ptr);
        print_log(od->real, "-real");
        print_log(od->imag, "-imag");
        mpfr_printf("c s exact %s %.8Re + %.8Rei\n", mpfr_prec_name(mpfr_precision), od->real, od->imag);
      } else if (mode == 7) {
        // MPFR numbers
        if (cnf.get_projected()) cout << "c s type pwmc" << endl;
        else cout << "c s type wmc" << endl;
        const ArjunNS::FMpfr* od = dynamic_cast<const ArjunNS::FMpfr*>(ptr);
        print_log(od->val);
        mpfr_printf("c s exact %s %.8Re\n", mpfr_prec_name(mpfr_precision), od->val);
      }
    } else if (mode == 3) {
      cout << "c s exact poly " << *cnt << endl;
    } else if (mode == 13) {
      cout << "c s exact laurent " << *cnt << endl;
    } else if (mode == 4) {
      cout << "c s exact parity " << *cnt << endl;
    } else if (mode == 5) {
      cout << "c s exact modprime " << *cnt << endl;
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
    case 2:
        fg = std::make_unique<ArjunNS::FGenComplex>();
        break;
    case 6:
        fg = std::make_unique<ArjunNS::FGenMPFComplex>(mpfr_precision);
        break;
    case 3:
        if (poly_nvars == -1) {
          cout << "c o [arjun] ERROR: Must provide number of polynomial vars for mode 3 via --npolyvars" << endl;
          exit(EXIT_FAILURE);
        }
        fg = std::make_unique<FGenPoly>(poly_nvars);
        break;
    case 13:
        if (poly_nvars == -1) {
          cout << "c o [arjun] ERROR: Must provide number of polynomial vars for mode 13 via --npolyvars" << endl;
          exit(EXIT_FAILURE);
        }
        fg = std::make_unique<LaurentPolyGen>(poly_nvars);
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
      cerr << "ERROR: you provided --inputfile but no file. Strange. Exiting. " << endl;
      exit(EXIT_FAILURE);
    } else if (files.size() == 1) {
      const string& fname = files[0];
      read_in_a_file(fname, &cnf, etof_conf.all_indep, fg);
    } else {
        cerr << "ERROR: you must only give one CNF as input (or none, and then we read from STDIN)" << endl;
        cout << "       You provided " << files.size() << " files: ";
        for (const auto& f: files) cout << f << " ";
        cout << endl;
        exit(EXIT_FAILURE);
    }
  }

  if (cnf.get_weighted() && conf.do_buddy) {
    cerr << "ERROR: Cannot run BuDDy with weighted CNF" << endl;
    exit(EXIT_FAILURE);
  }
  cnf.clean_idiotic_mccomp_weights();
  cnf.check_cnf_sampl_sanity();
  cnf.check_cnf_vars();
  verb_print(1, "CNF projection set size: " << cnf.get_sampl_vars().size());

  // Run Arjun
  if (!do_arjun) {
    cnf.renumber_sampling_vars_for_ganak();
    disconnected_allowed = true;
  } else run_arjun(cnf);
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
  conf.disconnected_allowed = disconnected_allowed;
  Ganak counter(conf, fg);
  setup_ganak(cnf, counter);
  run_weighted_counter(counter, cnf, start_time);
  return 0;
}
