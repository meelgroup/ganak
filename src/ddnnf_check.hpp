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

#pragma once

// On-the-fly checker for --compile (--compilecheck 1, always on with
// SLOW_DEBUG): sees each node exactly as it is written, independently of how the
// compiler built it, and aborts at the first one that is wrong, while the search
// that built it is still on the stack. Memory: one var bitset and one count per
// node.

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <gmpxx.h>

namespace GanakInt {

class DDNNFChecker {
public:
  DDNNFChecker(uint32_t _nvars, bool _expect_smooth)
      : nvars(_nvars), W(_nvars / 64 + 1), expect_smooth(_expect_smooth) {
    vs.assign(2 * W, 0);
    cnt.resize(2);
    cnt[0] = 0;
    cnt[1] = 1;
  }

  template <class Kids, class Units, class Free>
  void on_and(uint64_t id, const Kids& kids, const Units& units, const Free& free) {
    uint64_t* me = new_node(id);
    mpz_class c = 1;
    for (auto k : kids) {
      check_ref(id, (uint64_t)k);
      const uint64_t* kv = row((uint64_t)k);
      for (size_t q = 0; q < W; q++) {
        if (me[q] & kv[q]) fail(id, "AND children share var " + std::to_string(low(me[q] & kv[q], q)));
        me[q] |= kv[q];
      }
      c *= cnt[(uint64_t)k];
    }
    for (auto l : units) add_var(id, me, (uint32_t)std::abs((int)l), "unit");
    for (auto v : free) add_var(id, me, (uint32_t)v, "free var");
    c <<= (unsigned long)std::size(free);
    cnt[id] = std::move(c);
  }

  template <class Arcs>
  void on_or(uint64_t id, const Arcs& arcs) {
    uint64_t* me = new_node(id);
    mpz_class c = 0;
    bool have_first = false;
    for (size_t i = 0; i < arcs.size(); i++) {
      const uint64_t ch = (uint64_t)arcs[i].child;
      check_ref(id, ch);
      av.assign(row(ch), row(ch) + W);
      for (int l : arcs[i].lits) {
        const uint32_t v = (uint32_t)std::abs(l);
        if (v == 0 || v > nvars) fail(id, "literal out of range");
        if (bit(av.data(), v)) fail(id, "arc " + std::to_string(i) + " literal " + std::to_string(l) +
                                         " is on a var also below its child or on the arc twice");
        av[v >> 6] |= 1ULL << (v & 63);
      }
      for (size_t j = 0; j < i; j++) {
        bool clash = false;
        for (int l : arcs[i].lits)
          for (int m : arcs[j].lits) clash |= (l == -m);
        if (!clash) fail(id, "arcs " + std::to_string(j) + " and " + std::to_string(i) +
                             " share no complementary literal (not deterministic)");
      }
      if (ch != 0) {
        if (!have_first) { first = av; have_first = true; }
        else if (expect_smooth && av != first)
          fail(id, "arcs cover different vars (a literal missing from an arc?)");
      }
      for (size_t q = 0; q < W; q++) me[q] |= av[q];
      c += cnt[ch];
    }
    cnt[id] = std::move(c);
  }

  // The node built for a decision level must count what ganak counted for it.
  void expect_count(uint64_t id, const std::string& ganak_count) const {
    if (id < cnt.size() && cnt[id].get_str() != ganak_count)
      fail(id, "counts " + cnt[id].get_str() + " but ganak counted " + ganak_count + " for its level");
  }

private:
  uint32_t nvars;
  size_t W;
  bool expect_smooth;
  std::vector<uint64_t> vs;
  std::vector<mpz_class> cnt;
  std::vector<uint64_t> av, first;

  uint64_t* row(uint64_t id) { return vs.data() + id * W; }
  static bool bit(const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; }
  static uint32_t low(uint64_t w, size_t q) { return (uint32_t)(q * 64 + __builtin_ctzll(w)); }

  uint64_t* new_node(uint64_t id) {
    if (id != cnt.size()) fail(id, "node ids are not consecutive");
    vs.resize((id + 1) * W, 0);
    cnt.emplace_back();
    return row(id);
  }
  void check_ref(uint64_t id, uint64_t ch) const {
    if (ch >= id) fail(id, "child " + std::to_string(ch) + " is not an earlier node");
  }
  void add_var(uint64_t id, uint64_t* me, uint32_t v, const char* what) {
    if (v == 0 || v > nvars) fail(id, std::string(what) + " out of range");
    if (bit(me, v)) fail(id, std::string("AND ") + what + " " + std::to_string(v) + " shares its var with another conjunct");
    me[v >> 6] |= 1ULL << (v & 63);
  }
  [[noreturn]] static void fail(uint64_t id, const std::string& msg) {
    std::cerr << "*** d-DNNF CHECK FAILED at node " << id << ": " << msg << std::endl;
    std::abort();
  }
};

}  // namespace GanakInt
