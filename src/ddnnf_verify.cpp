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

// ddnnf-verify: checks a BDNNF circuit (spec-bdnnf.md) in streamed passes, with
// one var bitset per node: the format, the structural count, and every
// property the file claims (DECOMPOSABLE, DETERMINISTIC, SMOOTH, CLEAN).
//
// Usage:  ddnnf-verify <in.bdnnf>               streamed structural count only
//         ddnnf-verify <in.bdnnf> [--expect-count N] [--no-check-decomposable]
//                      [--no-strict] [-q]
// --strict (default with any option): CLEAN, no dead nodes, root last, and
// canonical (no AND under an AND, no non-root single-conjunct AND).

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
  uint64_t root = 0;
  uint64_t props = 0;
  bdnnf::KV footer_kv;
  std::string err;
  // Per node: is an AND; AND with an AND child (nested); AND with one conjunct.
  std::vector<char> is_and, nested, unary;
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
  out.is_and.assign(n, 0);
  out.nested.assign(n, 0);
  out.unary.assign(n, 0);
  r.scan([&](uint64_t id, const bdnnf::Node& nd) {
    mpz_class x;
    if (nd.kind == bdnnf::K_AND) {
      out.is_and[id] = 1;
      for (uint64_t k : nd.kids) out.nested[id] |= out.is_and[k];
      out.unary[id] = nd.kids.size() + (nd.units.empty() ? 0 : 1) + nd.free.size() == 1;
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
  out.root = r.root;
  out.props = r.props;
  out.footer_kv = r.footer_kv;
  return out;
}

// Backward scan from the root: reachable nodes (one bit each), then the strict
// form on them.
std::string check_strict(const std::string& path, const Streamed& st) {
  if (!(st.props & bdnnf::P_CLEAN)) return "file is not flagged CLEAN";
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
    if (st.nested[id] && err.empty())
      err = "nested AND: AND " + std::to_string(id) + " has an AND child (should be flattened)";
    if (st.unary[id] && id != st.root && err.empty())
      err = "unary AND: AND " + std::to_string(id) + " has a single conjunct (should be elided into its parent edge)";
    if (nd.kind == bdnnf::K_AND) for (uint64_t k : nd.kids) reach[k] = 1;
    else for (const auto& a : nd.arcs) reach[a.child] = 1;
  });
  return err;
}

}  // namespace

int main(int argc, char** argv) {
  std::string path;
  mpz_class expect_count;
  bool have_expect = false;
  bool check_decomp = true;
  bool strict = true;
  bool quiet = false;
  bool any_opt = false;

  auto usage = [&]() {
    std::cerr <<
      "Usage: " << argv[0] << " <in.bdnnf>\n"
      "         [--expect-count N] [--no-check-decomposable] [--no-strict] [-q|--quiet]\n"
      "\n"
      "With no option, prints the structural count.\n";
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
    else if (!a.empty() && a[0] == '-' && a != "-") {
      std::cerr << "ddnnf-verify: unknown option: " << a << "\n";
      usage();
    }
    else if (path.empty()) path = a;
    else { std::cerr << "ddnnf-verify: unexpected positional: " << a << "\n"; usage(); }
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
        check_decomp || (props & bdnnf::P_DECOMPOSABLE), props & bdnnf::P_SMOOTH);
    if (!st.err.empty()) fail(st.err);
    const std::string fc = bdnnf::Reader().get(st.footer_kv, "count");
    if (!fc.empty() && mpz_class(fc, 10) != st.count)
      fail("count " + st.count.get_str() + " != footer count " + fc);
    if (have_expect && st.count != expect_count)
      fail("count " + st.count.get_str() + " != expected " + expect_count.get_str());
    if (strict || (props & bdnnf::P_CLEAN)) {
      const std::string e = check_strict(path, st);
      if (!e.empty() && (strict || e.rfind("file is not", 0) != 0)) fail(e);
    }
    if (!quiet)
      std::cout << "ddnnf-verify: OK (count=" << st.count.get_str() << ", nodes=" << st.nnodes << ")\n";
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
