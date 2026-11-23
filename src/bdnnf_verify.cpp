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
// property the file claims (DECOMPOSABLE, DETERMINISTIC, SMOOTH, CLEAN).
//
// Usage:  bdnnf verify <in.bdnnf>               streamed structural count only
//         bdnnf verify <in.bdnnf> [--expect-count N] [-q]

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
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

  auto usage = [&]() {
    std::cerr <<
      "Usage: bdnnf verify <in.bdnnf>\n"
      "         [--expect-count N] [-q|--quiet]\n"
      "\n"
      "With no option, prints the structural count.\n";
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
    if (!quiet)
      std::cout << "bdnnf verify: OK (count=" << st.count.get_str() << ", nodes=" << st.nnodes << ")\n";
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
