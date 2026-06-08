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

// ddnnf-verify: structural-count + invariant checker for d4 .nnf files from
// `ganak --compile`. C++ port of tests/ddnnf_verify.py (CLI surface only). The
// Python version stays for tiny-N programmatic use from ddnnf_fuzz.py; this
// binary handles the big-circuit cases (multi-million-node
// random_nvN benches) that exhaust Python's stack and dict overhead.
//
// Usage:  ddnnf-verify <in.nnf>                       (default: --check-decomposable
//                                                      and --strict are both ON)
//         ddnnf-verify <in.nnf> --expect-count N
//         ddnnf-verify <in.nnf> --cnf <c.cnf>         (small CNFs only: brute-forces
//                                                      the model set 2^nvars)
//         ddnnf-verify <in.nnf> --no-check-decomposable   (opt out of decomp check)
//         ddnnf-verify <in.nnf> --no-strict           (opt out of strict-form check:
//                                                      root id 1, ids 1..N contiguous,
//                                                      no dead nodes)
//         ddnnf-verify <in.nnf> -q | --quiet          (suppress final "OK" message)
//
// `--check-decomposable` and `--strict` are accepted as no-ops for back-compat
// (they used to be opt-in). To get the raw structural count alone, opt out of
// both checks: `ddnnf-verify in.nnf --no-check-decomposable --no-strict`.

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "ddnnf_io.hpp"

namespace {

using ddnnf_io::Arc;

// Compact bitset over [0, nvars]. Fixed word count is set once per pass so the
// inner loops are word-parallel ORs/ANDs (same trick as ddnnf-cleanup).
struct VarSet {
  std::vector<uint64_t> w;
  void resize(size_t words) { w.assign(words, 0); }
  void set(int v) { w[v >> 6] |= 1ULL << (v & 63); }
  void or_with(const VarSet& o) {
    const size_t n = w.size();
    for (size_t i = 0; i < n; i++) w[i] |= o.w[i];
  }
  // Returns any var in (this & o), else 0.
  int first_intersect_var(const VarSet& o) const {
    const size_t n = w.size();
    for (size_t i = 0; i < n; i++) {
      uint64_t b = w[i] & o.w[i];
      if (b) return (int)(i * 64) + __builtin_ctzll(b);
    }
    return 0;
  }
};

void fail(const std::string& msg) {
  std::cout << "ddnnf-verify: FAIL: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

void die(const std::string& msg) {
  std::cerr << "ddnnf-verify: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

// Iterative post-order DFS body. Visits children before parents; calls
// `on_emit(n)` exactly once per reachable node, with all of n's children
// already emitted. `done` marks emitted nodes (sized to at least max_id+1
// by caller, zero-initialized). Avoids the 8 MB stack blowup at 14M nodes.
template <class OnEmit>
void post_order_visit(int root,
                      const std::vector<std::vector<Arc>>& arcs,
                      std::vector<char>& done,
                      OnEmit on_emit) {
  std::vector<std::pair<int, int>> stack;
  stack.reserve(64);
  if (done[root]) return;
  stack.emplace_back(root, 0);
  while (!stack.empty()) {
    auto& top = stack.back();
    int n = top.first;
    int idx = top.second;
    if (done[n]) { stack.pop_back(); continue; }
    int sz = (n < (int)arcs.size()) ? (int)arcs[n].size() : 0;
    if (idx >= sz) {
      on_emit(n);
      done[n] = 1;
      stack.pop_back();
      continue;
    }
    int c = arcs[n][idx].child;
    top.second = idx + 1;
    if (!done[c]) stack.emplace_back(c, 0);
  }
}

// Structural model count with GMP big-int memo (OR=sum, AND=product, t=1, f=0;
// arc literals ignored). Free vars are explicit OR(v,-v), so no smoothing.
mpz_class structural_count(int root,
                           const std::vector<char>& type,
                           const std::vector<std::vector<Arc>>& arcs) {
  std::vector<mpz_class> memo(type.size());
  std::vector<char> done(type.size(), 0);
  post_order_visit(root, arcs, done, [&](int n) {
    char t = type[n];
    if (t == 't')      memo[n] = 1;
    else if (t == 'f') memo[n] = 0;
    else if (t == 'o') {
      mpz_class r = 0;
      for (const auto& a : arcs[n]) r += memo[a.child];
      memo[n] = std::move(r);
    } else { // 'a'
      mpz_class r = 1;
      for (const auto& a : arcs[n]) r *= memo[a.child];
      memo[n] = std::move(r);
    }
  });
  return memo[root];
}

// Reachable-from-root marker, indexed by node id. Iterative DFS.
std::vector<char> reachable_mask(int root,
                                 const std::vector<std::vector<Arc>>& arcs,
                                 size_t sz) {
  std::vector<char> seen(sz, 0);
  std::vector<int> stack = {root};
  while (!stack.empty()) {
    int n = stack.back(); stack.pop_back();
    if (seen[n]) continue;
    seen[n] = 1;
    if (n < (int)arcs.size())
      for (const auto& a : arcs[n]) stack.push_back(a.child);
  }
  return seen;
}

// Strict d-DNNF check: every reachable AND node's children must have pairwise
// var-disjoint subtree variable sets (where a child's vars = arc lits ∪ that
// child's subtree vars). Returns (true, "") or (false, msg).
std::pair<bool, std::string> check_decomposable(
    int root,
    const std::vector<char>& type,
    const std::vector<std::vector<Arc>>& arcs) {
  int max_var = 0;
  for (const auto& kids : arcs)
    for (const auto& a : kids)
      for (int l : a.lits) { int v = l > 0 ? l : -l; if (v > max_var) max_var = v; }
  const size_t words = (size_t)(max_var / 64) + 1;

  std::vector<VarSet> sv(type.size());
  std::vector<char> done(type.size(), 0);
  post_order_visit(root, arcs, done, [&](int n) {
    sv[n].resize(words);
    for (const auto& a : arcs[n]) {
      for (int l : a.lits) sv[n].set(l > 0 ? l : -l);
      sv[n].or_with(sv[a.child]);
    }
  });

  for (int n = 0; n < (int)type.size(); n++) {
    if (!done[n] || type[n] != 'a') continue;
    VarSet seen; seen.resize(words);
    for (const auto& a : arcs[n]) {
      // child_vars = arc lits ∪ subtree(a.child)
      VarSet cv = sv[a.child];
      for (int l : a.lits) cv.set(l > 0 ? l : -l);
      int v = seen.first_intersect_var(cv);
      if (v != 0) {
        std::ostringstream m;
        m << "AND " << n << " children share var " << v;
        return {false, m.str()};
      }
      seen.or_with(cv);
    }
  }
  return {true, ""};
}

// Strict-form check: matches what ddnnf-cleanup produces -- root id 1, ids
// 1..N contiguous (all declared), no dead (unreachable) nodes.
std::pair<bool, std::string> check_strict(
    int root,
    const std::vector<char>& type,
    const std::vector<std::vector<Arc>>& arcs) {
  if (root != 1) return {false, "root id " + std::to_string(root) + " != 1"};
  // declared = nonzero type slots
  std::vector<int> declared;
  declared.reserve(type.size());
  for (int i = 0; i < (int)type.size(); i++) if (type[i] != 0) declared.push_back(i);
  // contiguous 1..N
  for (size_t k = 0; k < declared.size(); k++) {
    if (declared[k] != (int)k + 1) {
      std::ostringstream m;
      m << "node ids not contiguous 1..N (first gap before id " << declared[k] << ")";
      return {false, m.str()};
    }
  }
  // dead nodes: declared - reachable
  auto seen = reachable_mask(root, arcs, type.size());
  for (int id : declared) if (!seen[id]) {
    std::ostringstream m;
    m << "1 or more dead node(s) present; first: " << id;
    return {false, m.str()};
  }
  // nested AND: a literal-free arc from an AND to another AND is associatively
  // redundant -- the inner AND's conjuncts belong directly in the parent. The
  // streaming compiler emits these (e.g. wrap_lits() witness wrappers conjoined
  // by mk_and()); ddnnf-cleanup's flatten_ands pass removes them. A surviving one
  // means a non-canonical circuit (or a regression in that pass).
  for (int n = 0; n < (int)type.size(); n++) {
    if (!seen[n] || type[n] != 'a') continue;
    for (const auto& a : arcs[n]) {
      if (a.lits.empty() && a.child >= 0 && a.child < (int)type.size()
          && type[a.child] == 'a') {
        std::ostringstream m;
        m << "nested AND: AND " << n << " has a literal-free arc to AND "
          << a.child << " (should be flattened)";
        return {false, m.str()};
      }
    }
  }
  // unary AND: an AND with a single arc conjoins nothing -- it is its child plus
  // the arc's literals, so it belongs folded into its parent edge. wrap_lits()
  // emits these for projected / SAT-witnessed components; ddnnf-cleanup's
  // elide_unary_ands pass removes every non-root one. (The root may stay unary:
  // it carries the top-level literals and has no parent to fold into.)
  for (int n = 0; n < (int)type.size(); n++) {
    if (!seen[n] || type[n] != 'a' || n == root) continue;
    if (arcs[n].size() == 1) {
      std::ostringstream m;
      m << "unary AND: AND " << n << " has a single arc (should be elided into "
           "its parent edge)";
      return {false, m.str()};
    }
  }
  return {true, ""};
}

// ---------------- CNF brute-force model-set comparison ------------------

struct Cnf {
  int nvars = 0;
  std::vector<std::vector<int>> clauses;
};

Cnf read_cnf(const std::string& path) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open cnf: " + path);
  Cnf c;
  std::string line;
  while (std::getline(f, line)) {
    size_t i = 0;
    while (i < line.size() && std::isspace((unsigned char)line[i])) i++;
    if (i >= line.size() || line[i] == 'c') continue;
    if (line[i] == 'p') {
      std::istringstream ss(line);
      std::string p, cnf_kw; int nv, nc;
      ss >> p >> cnf_kw >> nv >> nc;
      c.nvars = nv;
      continue;
    }
    std::istringstream ss(line);
    std::vector<int> cl;
    int x;
    while (ss >> x) { if (x == 0) break; cl.push_back(x); }
    if (!cl.empty()) c.clauses.push_back(std::move(cl));
  }
  return c;
}

// Brute-force model set of a CNF -- only for tiny nvars (we cap at 30 below).
// Returns the set of satisfying assignments as bitmasks (bit v-1 = var v).
std::set<uint64_t> cnf_brute_models(const Cnf& c) {
  std::set<uint64_t> out;
  const uint64_t total = (uint64_t)1 << c.nvars;
  for (uint64_t bits = 0; bits < total; bits++) {
    bool sat = true;
    for (const auto& cl : c.clauses) {
      bool clsat = false;
      for (int l : cl) {
        int v = l > 0 ? l : -l;
        bool val = ((bits >> (v - 1)) & 1) != 0;
        if ((l > 0) == val) { clsat = true; break; }
      }
      if (!clsat) { sat = false; break; }
    }
    if (sat) out.insert(bits);
  }
  return out;
}

// Partial model: (set_mask, val_mask). Bit v-1 of set_mask = "var v is
// constrained"; bit v-1 of val_mask = "constrained to true". A model in
// the circuit's sense.
using PM = std::pair<uint64_t, uint64_t>;

// Extend partial `m` with `lits`; return false on conflict.
bool apply_arc_lits(PM& m, const std::vector<int>& lits) {
  for (int l : lits) {
    int v = l > 0 ? l : -l;
    uint64_t vbit = (uint64_t)1 << (v - 1);
    bool val = l > 0;
    if (m.first & vbit) {
      if ((((m.second >> (v - 1)) & 1) != 0) != val) return false;
    } else {
      m.first |= vbit;
      if (val) m.second |= vbit;
    }
  }
  return true;
}

// Circuit model set: walk the DAG considering arc lits, produce partial models
// per node, then expand free vars. For tiny CNFs only (exponential in
// unconstrained-var count per partial). Mirrors models() in ddnnf_verify.py.
std::set<uint64_t> circuit_models(int root,
                                  const std::vector<char>& type,
                                  const std::vector<std::vector<Arc>>& arcs,
                                  int nvars) {
  std::vector<std::vector<PM>> memo(type.size());
  std::vector<char> done(type.size(), 0);
  post_order_visit(root, arcs, done, [&](int n) {
    char t = type[n];
    auto& res = memo[n];
    if (t == 't') { res.push_back({0, 0}); return; }
    if (t == 'f') return;
    if (t == 'o') {
      for (const auto& a : arcs[n]) {
        for (const auto& cm : memo[a.child]) {
          PM mm = cm;
          if (apply_arc_lits(mm, a.lits)) res.push_back(mm);
        }
      }
      return;
    }
    // AND: cartesian product of children, dropping conflicts.
    res.push_back({0, 0});
    for (const auto& a : arcs[n]) {
      std::vector<PM> child_apps;
      child_apps.reserve(memo[a.child].size());
      for (const auto& cm : memo[a.child]) {
        PM mm = cm;
        if (apply_arc_lits(mm, a.lits)) child_apps.push_back(mm);
      }
      std::vector<PM> next;
      next.reserve(res.size() * child_apps.size());
      for (const auto& base : res) {
        for (const auto& cm : child_apps) {
          uint64_t conflict = base.first & cm.first;
          if ((base.second & conflict) != (cm.second & conflict)) continue;
          next.push_back({base.first | cm.first, base.second | cm.second});
        }
      }
      res = std::move(next);
      if (res.empty()) return;
    }
  });

  std::set<uint64_t> full;
  const uint64_t all_mask = (nvars >= 64) ? ~(uint64_t)0 : (((uint64_t)1 << nvars) - 1);
  for (const auto& m : memo[root]) {
    uint64_t free = (~m.first) & all_mask;
    int nfree = __builtin_popcountll(free);
    std::vector<int> pos;
    pos.reserve(nfree);
    for (int v = 0; v < nvars; v++) if ((free >> v) & 1) pos.push_back(v);
    const uint64_t total = (uint64_t)1 << nfree;
    for (uint64_t k = 0; k < total; k++) {
      uint64_t v = m.second;
      for (int i = 0; i < nfree; i++) if ((k >> i) & 1) v |= ((uint64_t)1 << pos[i]);
      full.insert(v);
    }
  }
  return full;
}

} // namespace

int main(int argc, char** argv) {
  // ---- arg parse: <nnf> + optional flags ----
  std::string nnf_path;
  std::string cnf_path;
  mpz_class expect_count;
  bool have_expect = false;
  bool check_decomp = true;   // default ON; opt out with --no-check-decomposable
  bool strict = true;         // default ON; opt out with --no-strict
  bool quiet = false;

  auto usage = [&]() {
    std::cerr <<
      "Usage: " << argv[0] << " <in.nnf>\n"
      "         [--expect-count N] [--cnf <file>]\n"
      "         [--no-check-decomposable] [--no-strict] [-q|--quiet]\n"
      "\n"
      "--check-decomposable and --strict are ON by default. Opt out with\n"
      "--no-check-decomposable / --no-strict; opt out of BOTH to print just\n"
      "the structural count.\n";
    std::exit(EXIT_FAILURE);
  };

  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") usage();
    else if (a == "-q" || a == "--quiet") quiet = true;
    else if (a == "--check-decomposable") check_decomp = true;   // no-op, kept for back-compat
    else if (a == "--strict") strict = true;                     // no-op, kept for back-compat
    else if (a == "--no-check-decomposable") check_decomp = false;
    else if (a == "--no-strict") strict = false;
    else if (a == "--expect-count") {
      if (i + 1 >= argc) usage();
      expect_count.set_str(argv[++i], 10);
      have_expect = true;
    }
    else if (a == "--cnf") {
      if (i + 1 >= argc) usage();
      cnf_path = argv[++i];
    }
    else if (a.size() > 0 && a[0] == '-') {
      std::cerr << "ddnnf-verify: unknown option: " << a << "\n";
      usage();
    }
    else if (nnf_path.empty()) nnf_path = a;
    else { std::cerr << "ddnnf-verify: unexpected positional: " << a << "\n"; usage(); }
  }
  if (nnf_path.empty()) usage();

  // ---- parse the .nnf ----
  std::ifstream in(nnf_path);
  if (!in.good()) die("cannot open: " + nnf_path);
  std::vector<char> type;
  std::vector<std::vector<Arc>> arcs;
  std::vector<int> order;
  std::string err;
  int root = ddnnf_io::parse_nnf(in, type, arcs, order, err);
  if (root == -1) die(err);

  // ---- structural count (always) ----
  mpz_class c = structural_count(root, type, arcs);

  const bool any_check = have_expect || !cnf_path.empty() || check_decomp || strict;
  if (!any_check) {
    // Bare invocation: print the count, exit 0 (matches Python's behavior).
    std::cout << c.get_str() << "\n";
    return EXIT_SUCCESS;
  }

  if (have_expect && c != expect_count) {
    fail("count " + c.get_str() + " != expected " + expect_count.get_str());
  }

  if (!cnf_path.empty()) {
    Cnf cnf = read_cnf(cnf_path);
    if (cnf.nvars > 30) die("--cnf brute-force only supports nvars <= 30 (got "
                            + std::to_string(cnf.nvars) + ")");
    auto true_models = cnf_brute_models(cnf);
    if (c != (long)true_models.size()) {
      fail("count " + c.get_str() + " != #models(CNF) " + std::to_string(true_models.size()));
    }
    auto circ = circuit_models(root, type, arcs, cnf.nvars);
    if (circ != true_models) {
      size_t extra = 0, missing = 0;
      for (auto x : circ) if (!true_models.count(x)) extra++;
      for (auto x : true_models) if (!circ.count(x)) missing++;
      fail("model set mismatch: +" + std::to_string(extra)
           + " extra, -" + std::to_string(missing) + " missing");
    }
  }

  if (check_decomp) {
    auto [ok, msg] = check_decomposable(root, type, arcs);
    if (!ok) fail("not decomposable: " + msg);
  }

  if (strict) {
    auto [ok, msg] = check_strict(root, type, arcs);
    if (!ok) fail(msg);
  }

  if (!quiet) {
    size_t nodes = 0;
    for (char t : type) if (t != 0) nodes++;
    std::cout << "ddnnf-verify: OK (count=" << c.get_str()
              << ", nodes=" << nodes << ")\n";
  }
  return EXIT_SUCCESS;
}
