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

// ddnnf-verify: structural count + invariant checker for BDNNF circuits
// (spec-bdnnf.md), raw from `ganak --compile` or cleaned by ddnnf-cleanup.
//
// Usage:  ddnnf-verify <in.bdnnf>             streamed structural count only
//         ddnnf-verify <in.bdnnf> [--expect-count N] [--cnf <c.cnf>]
//                      [--no-check-decomposable] [--no-strict] [-q]
// Any checking option loads the circuit and runs every check that is not opted
// out of, plus a check of each property the file claims. --cnf brute-forces the
// model set (small CNFs only). --strict: the circuit is CLEAN (root is the last
// node, no dead nodes) and canonical (no nested or non-root unary ANDs).

#include <algorithm>
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
using ddnnf_io::Circuit;

void fail(const std::string& msg) {
  std::cout << "ddnnf-verify: FAIL: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

void die(const std::string& msg) {
  std::cerr << "ddnnf-verify: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

// One forward scan, no graph built: structural count, plus the per-node property
// checks asked for (memory: one var bitset per node for decomposable/smooth).
struct Streamed {
  mpz_class count;
  uint64_t nnodes = 0;
  bdnnf::KV footer_kv;
  std::string err;
};
Streamed streamed_check(const std::string& path, bool det, bool decomp, bool smooth) {
  bdnnf::Reader r;
  r.open(path);
  Streamed out;
  const size_t n = r.have_footer ? (size_t)r.nnodes + 2 : 2;
  std::vector<mpz_class> val(n);
  val[0] = 0;
  val[1] = 1;
  const bool sets = decomp || smooth;
  const size_t W = (size_t)r.nvars / 64 + 1;
  std::vector<uint64_t> vs(sets ? n * W : 0, 0);
  std::vector<uint64_t> seen(W), av(W), first(W);
  auto row = [&](uint64_t id) { return vs.data() + id * W; };
  auto bit = [](uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
  auto set = [](uint64_t* s, uint32_t v) { s[v >> 6] |= 1ULL << (v & 63); };
  auto fail_once = [&](const std::string& m) { if (out.err.empty()) out.err = m; };
  r.scan([&](uint64_t id, const bdnnf::Node& nd) {
    mpz_class x;
    if (nd.kind == bdnnf::K_AND) {
      x = 1;
      for (uint64_t k : nd.kids) x *= val[k];
      x <<= nd.free.size();
      if (sets) {
        uint64_t* me = row(id);
        auto add = [&](uint32_t v) {
          if (decomp && bit(me, v)) fail_once("AND " + std::to_string(id) + " conjuncts share var " + std::to_string(v));
          set(me, v);
        };
        for (uint64_t k : nd.kids) {
          const uint64_t* kv = row(k);
          for (size_t w = 0; w < W; w++) {
            if (decomp && (me[w] & kv[w]))
              fail_once("AND " + std::to_string(id) + " children share var " +
                        std::to_string(w * 64 + __builtin_ctzll(me[w] & kv[w])));
            me[w] |= kv[w];
          }
        }
        for (int l : nd.units) add((uint32_t)std::abs(l));
        for (uint32_t v : nd.free) add(v);
      }
    } else {
      x = 0;
      for (const auto& a : nd.arcs) x += val[a.child];
      if (det)
        for (size_t i = 0; i < nd.arcs.size(); i++)
          for (size_t j = i + 1; j < nd.arcs.size(); j++) {
            bool clash = false;
            for (int l : nd.arcs[i].lits)
              for (int m : nd.arcs[j].lits) clash |= (l == -m);
            if (!clash)
              fail_once("OR " + std::to_string(id) + " arcs " + std::to_string(i) + " and " +
                        std::to_string(j) + " share no complementary literal (not DETERMINISTIC)");
          }
      if (sets) {
        uint64_t* me = row(id);
        bool have_first = false;
        for (const auto& a : nd.arcs) {
          const uint64_t* cv = row(a.child);
          std::copy(cv, cv + W, av.begin());
          for (int l : a.lits) {
            const uint32_t v = (uint32_t)std::abs(l);
            if (decomp && bit(av.data(), v))
              fail_once("OR " + std::to_string(id) + " arc literal on var " + std::to_string(v) +
                        " also occurs below its child " + std::to_string(a.child));
            set(av.data(), v);
          }
          for (size_t w = 0; w < W; w++) me[w] |= av[w];
          if (smooth && a.child != 0) {
            if (!have_first) { first = av; have_first = true; }
            else if (first != av) fail_once("OR " + std::to_string(id) + " arcs mention different vars (not SMOOTH)");
          }
        }
      }
    }
    val[id] = std::move(x);
  });
  if (smooth && r.root != 0)
    for (uint32_t v = 1; v <= r.nvars; v++)
      if (!bit(row(r.root), v)) { fail_once("root does not mention var " + std::to_string(v) + " (not SMOOTH)"); break; }
  out.count = val[r.root];
  out.nnodes = r.nnodes;
  out.footer_kv = r.footer_kv;
  return out;
}

mpz_class structural_count(const Circuit& c) {
  std::vector<mpz_class> memo(c.type.size());
  std::vector<char> done(c.type.size(), 0);
  ddnnf_io::post_order(c, done, [&](int n) {
    const char t = c.type[n];
    if (t == 't') memo[n] = 1;
    else if (t == 'f') memo[n] = 0;
    else if (t == 'o') {
      mpz_class r = 0;
      for (const auto& a : c.arcs[n]) r += memo[a.child];
      memo[n] = std::move(r);
    } else {
      mpz_class r = 1;
      for (const auto& a : c.arcs[n]) r *= memo[a.child];
      memo[n] = std::move(r);
    }
  });
  return memo[c.root];
}

std::pair<bool, std::string> check_strict(const Circuit& c) {
  if (!(c.props & bdnnf::P_CLEAN)) return {false, "file is not flagged CLEAN"};
  if (c.nfile > 0 && c.root != c.nfile + 1)
    return {false, "root " + std::to_string(c.root) + " is not the last node"};
  std::vector<char> seen(c.type.size(), 0);
  ddnnf_io::post_order(c, seen, [](int) {});
  for (int id = 2; id < c.nfile + 2; id++)
    if (!seen[id]) return {false, "1 or more dead node(s) present; first: " + std::to_string(id)};
  for (int n = 0; n < (int)c.type.size(); n++) {
    if (!seen[n] || c.type[n] != 'a') continue;
    for (const auto& a : c.arcs[n])
      if (a.lits.empty() && c.type[a.child] == 'a') {
        std::ostringstream m;
        m << "nested AND: AND " << n << " has a literal-free arc to AND " << a.child
          << " (should be flattened)";
        return {false, m.str()};
      }
    if (n != c.root && c.arcs[n].size() == 1) {
      std::ostringstream m;
      m << "unary AND: AND " << n << " has a single conjunct (should be elided into its parent edge)";
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
std::set<uint64_t> circuit_models(const ddnnf_io::Circuit& cc, int nvars) {
  const auto& type = cc.type;
  const auto& arcs = cc.arcs;
  const int root = cc.root;
  std::vector<std::vector<PM>> memo(type.size());
  std::vector<char> done(type.size(), 0);
  ddnnf_io::post_order(cc, done, [&](int n) {
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
  std::string path;
  std::string cnf_path;
  mpz_class expect_count;
  bool have_expect = false;
  bool check_decomp = true;
  bool strict = true;
  bool quiet = false;
  bool any_opt = false;

  auto usage = [&]() {
    std::cerr <<
      "Usage: " << argv[0] << " <in.bdnnf>\n"
      "         [--expect-count N] [--cnf <file>]\n"
      "         [--no-check-decomposable] [--no-strict] [-q|--quiet]\n"
      "\n"
      "With no option, prints the structural count (streamed, no graph built).\n";
    std::exit(EXIT_FAILURE);
  };

  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") usage();
    else if (a == "-q" || a == "--quiet") { quiet = true; any_opt = true; }
    else if (a == "--check-decomposable") { check_decomp = true; any_opt = true; }
    else if (a == "--strict") { strict = true; any_opt = true; }
    else if (a == "--no-check-decomposable") { check_decomp = false; any_opt = true; }
    else if (a == "--no-strict") { strict = false; any_opt = true; }
    else if (a == "--expect-count") {
      if (i + 1 >= argc) usage();
      expect_count.set_str(argv[++i], 10);
      have_expect = true;
      any_opt = true;
    }
    else if (a == "--cnf") {
      if (i + 1 >= argc) usage();
      cnf_path = argv[++i];
      any_opt = true;
    }
    else if (!a.empty() && a[0] == '-' && a != "-") {
      std::cerr << "ddnnf-verify: unknown option: " << a << "\n";
      usage();
    }
    else if (path.empty()) path = a;
    else { std::cerr << "ddnnf-verify: unexpected positional: " << a << "\n"; usage(); }
  }
  if (path.empty()) usage();

  if (!any_opt) {
    try {
      std::cout << streamed_check(path, false, false, false).count.get_str() << "\n";
    } catch (const bdnnf::Error& e) {
      die(path + ": " + e.what());
    }
    return EXIT_SUCCESS;
  }

  uint64_t props = 0;
  try {
    bdnnf::Reader pr;
    pr.open(path);
    props = pr.props;
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  // The graph is only needed for the model set and the strict (canonical) form.
  if (cnf_path.empty() && !strict) {
    Streamed st;
    try {
      st = streamed_check(path, props & bdnnf::P_DETERMINISTIC,
                          check_decomp || (props & bdnnf::P_DECOMPOSABLE), props & bdnnf::P_SMOOTH);
    } catch (const bdnnf::Error& e) { die(e.what()); }
    if (!st.err.empty()) fail(st.err);
    const mpz_class& sc = st.count;
    const uint64_t nn = st.nnodes;
    const std::string fc = bdnnf::Reader().get(st.footer_kv, "count");
    if (!fc.empty() && mpz_class(fc, 10) != sc) fail("count " + sc.get_str() + " != footer count " + fc);
    if (have_expect && sc != expect_count)
      fail("count " + sc.get_str() + " != expected " + expect_count.get_str());
    if (!quiet) std::cout << "ddnnf-verify: OK (count=" << sc.get_str() << ", nodes=" << nn << ", streamed)\n";
    return EXIT_SUCCESS;
  }

  Circuit c;
  std::string err;
  if (!ddnnf_io::load(path, c, err)) die(err);

  const mpz_class cnt = structural_count(c);
  {
    mpz_class sc;
    try { sc = streamed_check(path, false, false, false).count; } catch (const bdnnf::Error& e) { die(e.what()); }
    if (sc != cnt) fail("streamed count " + sc.get_str() + " != graph count " + cnt.get_str());
  }
  const std::string fcount = bdnnf::Reader().get(c.footer_kv, "count");
  if (!fcount.empty() && mpz_class(fcount, 10) != cnt)
    fail("count " + cnt.get_str() + " != footer count " + fcount);

  if (have_expect && cnt != expect_count)
    fail("count " + cnt.get_str() + " != expected " + expect_count.get_str());

  if (!cnf_path.empty()) {
    Cnf cnf = read_cnf(cnf_path);
    if (cnf.nvars > 30) die("--cnf brute-force only supports nvars <= 30 (got "
                            + std::to_string(cnf.nvars) + ")");
    auto true_models = cnf_brute_models(cnf);
    if (cnt != (long)true_models.size())
      fail("count " + cnt.get_str() + " != #models(CNF) " + std::to_string(true_models.size()));
    auto circ = circuit_models(c, cnf.nvars);
    if (circ != true_models) {
      size_t extra = 0, missing = 0;
      for (auto x : circ) if (!true_models.count(x)) extra++;
      for (auto x : true_models) if (!circ.count(x)) missing++;
      fail("model set mismatch: +" + std::to_string(extra)
           + " extra, -" + std::to_string(missing) + " missing");
    }
  }

  if (check_decomp || (c.props & bdnnf::P_DECOMPOSABLE)) {
    auto [ok, msg] = ddnnf_io::check_decomposable(c);
    if (!ok) fail("not decomposable: " + msg);
  }
  if (c.props & bdnnf::P_SMOOTH) {
    auto [ok, msg] = ddnnf_io::check_smooth(c);
    if (!ok) fail("flagged SMOOTH but " + msg);
  }
  if (c.props & bdnnf::P_DETERMINISTIC) {
    auto [ok, msg] = ddnnf_io::check_deterministic(c);
    if (!ok) fail("flagged DETERMINISTIC but " + msg);
  }
  if (strict) {
    auto [ok, msg] = check_strict(c);
    if (!ok) fail(msg);
  }

  if (!quiet)
    std::cout << "ddnnf-verify: OK (count=" << cnt.get_str() << ", nodes=" << c.nfile << ")\n";
  return EXIT_SUCCESS;
}
