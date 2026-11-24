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

// bdnnf verify: checks a BDNNF circuit (spec-bdnnf.md) in streamed passes, with
// one var bitset per node: the format, the structural count, and every
// property the file claims (DECOMPOSABLE, DETERMINISTIC, SMOOTH, CLEAN); with
// --cnf, that every model of the file satisfies the CNF.
//
// Usage:  bdnnf verify <in.bdnnf>               streamed structural count only
//         bdnnf verify <in.bdnnf> [--expect-count N] [--cnf F] [-q]

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gmpxx.h>

#include "bdnnf.hpp"

namespace {

void fail(const std::string& msg) {
  std::cout << "bdnnf verify: FAIL: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

void die(const std::string& msg) {
  std::cerr << "bdnnf verify: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

// One forward scan, no graph built: structural count, plus the per-node property
// checks asked for (memory: one var bitset per node for decomposable/smooth).
struct Streamed {
  mpz_class count;
  uint64_t nnodes = 0;
  uint64_t root = 0;
  uint64_t props = 0;
  bdnnf::KV footer_kv;
  std::string err;
};
Streamed streamed_check(const std::string& path, bool det, bool decomp, bool smooth) {
  bdnnf::Reader r;
  r.open(path);
  Streamed out;
  const size_t n = (size_t)r.nnodes + 2;
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
  // Projected: arcs must also differ on the projection, or the projected count
  // adds a projection twice. A var outside it that is a function of the
  // projection (header key projection-defined) separates them too.
  std::vector<char> sep(r.nvars + 1, 1);
  const bool projected = r.has(r.header_kv, "projection");
  if (projected) {
    sep.assign(r.nvars + 1, 0);
    for (const char* key : {"projection", "projection-defined"})
      for (uint32_t v : bdnnf::parse_var_list(r.get(r.header_kv, key)))
        if (v >= 1 && v <= r.nvars) sep[v] = 1;
  }
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
            bool clash = false, sep_clash = false;
            for (int l : nd.arcs[i].lits)
              for (int m : nd.arcs[j].lits)
                if (l == -m) {
                  clash = true;
                  sep_clash |= sep[(uint32_t)std::abs(l)] != 0;
                }
            if (!clash)
              fail_once("OR " + std::to_string(id) + " arcs " + std::to_string(i) + " and " +
                        std::to_string(j) + " share no complementary literal (not DETERMINISTIC)");
            else if (!sep_clash)
              fail_once("OR " + std::to_string(id) + " arcs " + std::to_string(i) + " and " +
                        std::to_string(j) + " differ only outside the projection: the projected count adds "
                        "their common projections twice");
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
  // A defined var is given by its definition: the circuit must not mention it.
  if (!r.defs.empty()) {
    std::vector<uint64_t> all(W, 0);
    for (size_t id = 2; id < n && sets; id++)
      for (size_t q = 0; q < W; q++) all[q] |= row(id)[q];
    for (const auto& [v, rf] : r.defs.defs)
      if (sets && bit(all.data(), v)) fail_once("defined var " + std::to_string(v) + " is mentioned by the circuit");
    std::vector<signed char> assign(r.nvars + 1, 0);
    try { bdnnf::eval_defs(r.defs, r.nvars, assign); } catch (const bdnnf::Error& e) { fail_once(e.what()); }
  }
  if (smooth && r.root != 0)
    for (uint32_t v = 1; v <= r.nvars; v++)
      if (!bit(row(r.root), v)) { fail_once("root does not mention var " + std::to_string(v) + " (not SMOOTH)"); break; }
  out.count = val[r.root];
  out.nnodes = r.nnodes;
  out.root = r.root;
  out.props = r.props;
  out.footer_kv = r.footer_kv;
  return out;
}

// Unweighted projected count, right for a file not flagged SMOOTH too
// (unlike the structural count).
mpz_class model_count(const std::string& path) {
  bdnnf::Reader r;
  r.open(path);
  const uint32_t nv = r.nvars;
  std::vector<char> in_proj(nv + 1, 1);
  if (r.has(r.header_kv, "projection")) {
    in_proj.assign(nv + 1, 0);
    for (uint32_t v : bdnnf::parse_var_list(r.get(r.header_kv, "projection")))
      if (v >= 1 && v <= nv) in_proj[v] = 1;
  }
  for (const auto& [v, rf] : r.defs.defs) in_proj[v] = 0;
  in_proj[0] = 0;
  const size_t n = (size_t)r.nnodes + 2;
  const size_t W = (size_t)nv / 64 + 1;
  std::vector<mpz_class> val(n);
  std::vector<uint64_t> vs(n * W, 0);
  std::vector<char> dead(n, 0);
  val[1] = 1;
  dead[0] = 1;
  auto row = [&](uint64_t id) { return vs.data() + id * W; };
  auto set = [](uint64_t* s, uint32_t v) { s[v >> 6] |= 1ULL << (v & 63); };
  auto conflicting = [](const std::vector<int>& ls) {
    for (size_t i = 0; i < ls.size(); i++)
      for (size_t j = i + 1; j < ls.size(); j++) if (ls[i] == -ls[j]) return true;
    return false;
  };
  auto pv = [&](int l) { return in_proj[(uint32_t)std::abs(l)] != 0; };
  std::vector<uint64_t> av(W);
  r.scan([&](uint64_t id, const bdnnf::Node& nd) {
    uint64_t* me = row(id);
    if (nd.kind == bdnnf::K_AND) {
      bool d = conflicting(nd.units);
      for (uint64_t k : nd.kids) d |= dead[k] != 0;
      if (d) { dead[id] = 1; return; }
      mpz_class x = 1;
      for (uint64_t k : nd.kids) {
        x *= val[k];
        for (size_t q = 0; q < W; q++) me[q] |= row(k)[q];
      }
      for (int l : nd.units) if (pv(l)) set(me, (uint32_t)std::abs(l));
      for (uint32_t v : nd.free) if (in_proj[v]) { x <<= 1; set(me, v); }
      val[id] = std::move(x);
      return;
    }
    bool any = false;
    for (const auto& a : nd.arcs) {
      if (dead[a.child] || conflicting(a.lits)) continue;
      any = true;
      for (size_t q = 0; q < W; q++) me[q] |= row(a.child)[q];
      for (int l : a.lits) if (pv(l)) set(me, (uint32_t)std::abs(l));
    }
    if (!any) { dead[id] = 1; return; }
    mpz_class x = 0;
    for (const auto& a : nd.arcs) {
      if (dead[a.child] || conflicting(a.lits)) continue;
      std::copy(row(a.child), row(a.child) + W, av.begin());
      for (int l : a.lits) if (pv(l)) set(av.data(), (uint32_t)std::abs(l));
      size_t gap = 0;
      for (size_t q = 0; q < W; q++) gap += __builtin_popcountll(me[q] & ~av[q]);
      x += mpz_class(val[a.child]) << gap;
    }
    val[id] = std::move(x);
  });
  if (dead[r.root]) return 0;
  mpz_class res = val[r.root];
  size_t gap = 0;
  for (uint32_t v = 1; v <= nv; v++)
    if (in_proj[v] && !((row(r.root)[v >> 6] >> (v & 63)) & 1)) gap++;
  return res << gap;
}

// --cnf: every model of the file satisfies every clause, i.e. the circuit and
// the clause's negation have no common model. On a decomposable circuit that is
// a bottom-up SAT pass, 64 clauses per pass, over all vars: the projection does
// not matter, a var on no path takes any value. A clause with a defined var is
// skipped: check the circuit against the CNF it was compiled from (ganak
// --arjundebugcnf), which has none, and the definitions with ganak
// --compilecheckthorough.
struct CnfCheck {
  uint64_t clauses = 0, skipped = 0;
  std::string err;
};

CnfCheck check_cnf(bdnnf::Reader& r, const std::string& cnf_path) {
  const uint32_t nv = r.nvars;
  const size_t n = (size_t)r.nnodes + 2;
  std::vector<char> defined((size_t)nv + 1, 0);
  for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;

  CnfCheck out;
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  std::vector<uint64_t> kill(2 * (size_t)nv + 2, 0), sm(n, 0);
  std::vector<std::vector<int>> qs;
  std::vector<std::vector<int>> cls;  // of the queries in qs

  auto run_chunk = [&]() {
    if (qs.empty() || !out.err.empty()) return;
    const size_t k = qs.size();
    const uint64_t all = k == 64 ? ~0ULL : (1ULL << k) - 1;
    for (size_t j = 0; j < k; j++)
      for (int l : qs[j]) kill[lidx(-l)] |= 1ULL << j;
    auto killed = [&](const std::vector<int32_t>& ls) {
      uint64_t m = 0;
      for (int l : ls) m |= kill[lidx(l)];
      return m;
    };
    sm[0] = 0;
    sm[1] = all;
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      uint64_t m;
      if (nd.kind == bdnnf::K_AND) {
        m = all & ~killed(nd.units);
        for (uint64_t c : nd.kids) m &= sm[c];
      } else {
        m = 0;
        for (const auto& a : nd.arcs) m |= sm[a.child] & ~killed(a.lits);
      }
      sm[id] = m;
    });
    const uint64_t sat = sm[r.root];
    for (size_t j = 0; j < k && out.err.empty(); j++) {
      if (!((sat >> j) & 1)) continue;
      out.err = "clause (";
      for (int l : cls[j]) out.err += std::to_string(l) + " ";
      out.err += "0) of " + cnf_path + " is not implied: the circuit has a model with";
      for (int l : qs[j]) out.err += " " + std::to_string(l);
    }
    for (size_t j = 0; j < k; j++)
      for (int l : qs[j]) kill[lidx(-l)] = 0;
    qs.clear();
    cls.clear();
  };

  std::ifstream f(cnf_path);
  if (!f.good()) throw bdnnf::Error("cannot open " + cnf_path);
  std::string line;
  std::vector<int> cl;
  auto add_clause = [&]() {
    out.clauses++;
    bool taut = false, has_def = false;
    for (int l : cl) {
      has_def |= defined[(uint32_t)std::abs(l)] != 0;
      taut |= std::find(cl.begin(), cl.end(), -l) != cl.end();
    }
    if (taut || has_def) {
      out.skipped += has_def ? 1 : 0;
      cl.clear();
      return;
    }
    std::vector<int> q;
    for (int l : cl) if (std::find(q.begin(), q.end(), -l) == q.end()) q.push_back(-l);
    qs.push_back(std::move(q));
    cls.push_back(cl);
    cl.clear();
    if (qs.size() == 64) run_chunk();
  };
  while (std::getline(f, line) && out.err.empty()) {
    if (line.empty() || line[0] == 'c' || line[0] == 'p') continue;
    std::istringstream ss(line);
    long x;
    while (ss >> x) {
      if (x == 0) {
        if (!cl.empty()) add_clause();
        continue;
      }
      if ((uint64_t)std::labs(x) > nv)
        throw bdnnf::Error("CNF var " + std::to_string(std::labs(x)) + " is beyond the circuit's " +
                           std::to_string(nv) + " vars");
      cl.push_back((int)x);
    }
  }
  if (!cl.empty()) add_clause();
  run_chunk();
  return out;
}

// CLEAN: root last, no dead nodes (backward scan, one bit per node).
std::string check_clean(const std::string& path, const Streamed& st) {
  if (st.nnodes > 0 && st.root != st.nnodes + 1)
    return "root " + std::to_string(st.root) + " is not the last node";
  std::vector<char> reach(st.nnodes + 2, 0);
  reach[st.root] = 1;
  bdnnf::Reader r;
  r.open(path);
  std::string err;
  r.scan_backward([&](uint64_t id, const bdnnf::Node& nd) {
    if (!reach[id]) {
      if (err.empty()) err = "1 or more dead node(s) present; first: " + std::to_string(id);
      return;
    }
    if (nd.kind == bdnnf::K_AND) for (uint64_t k : nd.kids) reach[k] = 1;
    else for (const auto& a : nd.arcs) reach[a.child] = 1;
  });
  return err;
}

}  // namespace

int bdnnf_verify(int argc, char** argv) {
  std::string path;
  mpz_class expect_count;
  bool have_expect = false;
  bool quiet = false;
  bool any_opt = false;
  std::vector<std::string> cnfs;

  auto usage = [&]() {
    std::cerr <<
      "Usage: bdnnf verify <in.bdnnf>\n"
      "         [--expect-count N] [--cnf F] [-q|--quiet]\n"
      "\n"
      "With no option, prints the structural count. --cnf F (repeatable): every\n"
      "model of the file must satisfy every clause of F without a defined var\n"
      "(checked exactly, not by sampling).\n";
    std::exit(EXIT_FAILURE);
  };

  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") usage();
    else if (a == "-q" || a == "--quiet") { quiet = true; any_opt = true; }
    else if (a == "--expect-count") {
      if (i + 1 >= argc) usage();
      const std::string n = argv[++i];
      // set_str alone skips whitespace and needs a separate result check
      if (n.empty() || !std::all_of(n.begin(), n.end(), [](unsigned char ch) { return ch >= '0' && ch <= '9'; })
          || expect_count.set_str(n, 10) != 0) {
        std::cerr << "bdnnf verify: --expect-count needs a non-negative integer, got: " << n << "\n";
        usage();
      }
      have_expect = true;
      any_opt = true;
    }
    else if (a == "--cnf") {
      if (i + 1 >= argc) usage();
      cnfs.emplace_back(argv[++i]);
      any_opt = true;
    }
    else if (!a.empty() && a[0] == '-' && a != "-") {
      std::cerr << "bdnnf verify: unknown option: " << a << "\n";
      usage();
    }
    else if (path.empty()) path = a;
    else { std::cerr << "bdnnf verify: unexpected positional: " << a << "\n"; usage(); }
  }
  if (path.empty()) usage();

  try {
    if (!any_opt) {
      std::cout << streamed_check(path, false, false, false).count.get_str() << "\n";
      return EXIT_SUCCESS;
    }
    uint64_t props;
    {
      bdnnf::Reader pr;
      pr.open(path);
      props = pr.props;
    }
    const Streamed st = streamed_check(path, props & bdnnf::P_DETERMINISTIC,
        true, props & bdnnf::P_SMOOTH);
    if (!st.err.empty()) fail(st.err);
    const std::string fc = bdnnf::Reader().get(st.footer_kv, "count");
    if (!fc.empty() && mpz_class(fc, 10) != st.count)
      fail("count " + st.count.get_str() + " != footer count " + fc);
    if (have_expect && st.count != expect_count)
      fail("count " + st.count.get_str() + " != expected " + expect_count.get_str());
    // The structural count is the model count only for a d-DNNF that is smooth too
    const uint64_t all3 = bdnnf::P_DECOMPOSABLE | bdnnf::P_DETERMINISTIC | bdnnf::P_SMOOTH;
    if (have_expect && (props & all3) != all3) {
      const mpz_class mc = model_count(path);
      if (mc != expect_count)
        fail("model count " + mc.get_str() + " != expected " + expect_count.get_str() +
             " (structural count " + st.count.get_str() + ")");
    }
    if (props & bdnnf::P_CLEAN) {
      const std::string e = check_clean(path, st);
      if (!e.empty()) fail(e);
    }
    if (!cnfs.empty()) {
      bdnnf::Reader r;
      r.open(path);
      // small circuits stay in memory for the passes; a scan decodes the definitions
      if (std::filesystem::file_size(path) <= (64ULL << 20)) r.cache_all();
      else r.scan([](uint64_t, const bdnnf::Node&) {});
      for (const auto& cnf : cnfs) {
        const CnfCheck cc = check_cnf(r, cnf);
        if (!cc.err.empty()) fail(cc.err);
        if (!quiet)
          std::cout << "bdnnf verify: every model satisfies the " << cc.clauses - cc.skipped << " clauses of "
                    << cnf << (cc.skipped ? " without a defined var (" + std::to_string(cc.skipped) +
                               " have one)" : std::string()) << "\n";
      }
    }
    if (!quiet)
      std::cout << "bdnnf verify: OK (count=" << st.count.get_str() << ", nodes=" << st.nnodes << ")\n";
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
