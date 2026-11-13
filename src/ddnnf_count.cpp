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

// ddnnf-count: exact (weighted, projected) model count of a BDNNF circuit in one
// streamed pass: one var bitset and one rational per node. Smooths on the fly,
// so it needs a decomposable, deterministic circuit but not a SMOOTH one.
// Weights ('c p weight', Arjun's semantics) and projection ('c p show') come
// from the CNF; literals of vars outside the projection weigh 1.
//
// --sample N also draws N models (proportionally to the unweighted smoothed
// counts, so uniformly over the projection when the circuit is decomposable and
// deterministic; other vars not fixed on the path get random values) and checks
// each against the CNF's clauses: a circuit too weak fails this, one too strong
// fails the count. Random access through the frame index: no graph in memory.
//
// --assume L1,L2,... conditions on literals of projection vars (weight 0 on
// their negations): the count must equal ganak's on the CNF plus those units.
//
// Usage:  ddnnf-count <in.bdnnf> [--cnf <file.cnf>] [--assume L1,L2,...]
//                     [--sample N] [--seed S]

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gmpxx.h>

#include "bdnnf.hpp"

namespace {

[[noreturn]] void die(const std::string& msg) {
  std::cerr << "ddnnf-count: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

// "a", "a/b", "-0.25", "1e-3"-free decimals.
mpq_class parse_num(const std::string& s) {
  mpq_class r;
  const auto dot = s.find('.');
  if (dot == std::string::npos) {
    if (r.set_str(s, 10) != 0) die("bad weight: " + s);
  } else {
    std::string digits = s.substr(0, dot) + s.substr(dot + 1);
    if (digits.empty() || digits == "-" || digits == "+") digits += "0";
    mpz_class num(digits, 10), den = 1;
    for (size_t i = dot + 1; i < s.size(); i++) den *= 10;
    r = mpq_class(num, den);
  }
  r.canonicalize();
  return r;
}

struct Weights {
  bool weighted = false;
  std::map<int, mpq_class> w;
  bool projected = false;
  std::vector<char> in_proj;  // by var
};

Weights read_cnf(const std::string& path, uint32_t nvars) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open " + path);
  Weights W;
  std::vector<char> proj(nvars + 1, 0);
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream ss(line);
    std::string a, b, c;
    ss >> a >> b;
    if (a != "c") continue;
    if (b == "ind") {
      W.projected = true;
      long v;
      while (ss >> v && v != 0) if (v > 0 && (uint64_t)v <= nvars) proj[v] = 1;
      continue;
    }
    if (b != "p") continue;
    ss >> c;
    if (c == "weight") {
      long l;
      std::string ws;
      ss >> l >> ws;
      const mpq_class x = parse_num(ws);
      W.weighted = true;
      if (!W.w.count((int)l) && !W.w.count((int)-l)) W.w[(int)-l] = 1 - x;
      W.w[(int)l] = x;
    } else if (c == "show") {
      W.projected = true;
      long v;
      while (ss >> v && v != 0) if (v > 0 && (uint64_t)v <= nvars) proj[v] = 1;
    }
  }
  if (!W.projected) proj.assign(nvars + 1, 1);
  W.in_proj = std::move(proj);
  return W;
}

std::vector<std::vector<int>> read_clauses(const std::string& path) {
  std::ifstream f(path);
  std::vector<std::vector<int>> cls;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == 'c' || line[0] == 'p') continue;
    std::istringstream ss(line);
    std::vector<int> cl;
    long x;
    while (ss >> x && x != 0) cl.push_back((int)x);
    if (!cl.empty()) cls.push_back(std::move(cl));
  }
  return cls;
}

// Random access to nodes through the frame index, one decoded frame cached.
struct NodeAccess {
  bdnnf::Reader& r;
  size_t cur = SIZE_MAX;
  std::vector<uint8_t> raw;
  std::vector<size_t> starts;
  explicit NodeAccess(bdnnf::Reader& rd) : r(rd) {}
  void get(uint64_t id, bdnnf::Node& n) {
    const size_t fi = r.frame_of(id);
    if (fi != cur) {
      r.load_frame(fi, raw);
      starts.clear();
      bdnnf::Cursor c(raw.data(), raw.size());
      for (uint32_t i = 0; i < r.frames[fi].nnodes; i++) {
        starts.push_back((size_t)(c.p - raw.data()));
        bdnnf::decode_node(c, r.frames[fi].first_id + i, r.nvars, n);
      }
      cur = fi;
    }
    const size_t k = (size_t)(id - r.frames[fi].first_id);
    bdnnf::Cursor c(raw.data() + starts[k], raw.size() - starts[k]);
    bdnnf::decode_node(c, id, r.nvars, n);
  }
};

}  // namespace

int main(int argc, char** argv) {
  std::string path, cnf;
  uint64_t samples = 0, seed = 1;
  std::vector<int> assume;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--cnf" && i + 1 < argc) cnf = argv[++i];
    else if (a == "--sample" && i + 1 < argc) samples = std::stoull(argv[++i]);
    else if (a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
    else if (a == "--assume" && i + 1 < argc) {
      std::string t = argv[++i];
      for (char& ch : t) if (ch == ',') ch = ' ';
      std::istringstream ss(t);
      long l;
      while (ss >> l) if (l != 0) assume.push_back((int)l);
    }
    else if (path.empty()) path = a;
    else die("usage: ddnnf-count <in.bdnnf> [--cnf <file.cnf>]");
  }
  if (path.empty()) die("usage: ddnnf-count <in.bdnnf> [--cnf <file.cnf>]");
  try {
    bdnnf::Reader r;
    r.open(path);
    const uint32_t nv = r.nvars;
    Weights W;
    if (!cnf.empty()) W = read_cnf(cnf, nv);
    else W.in_proj.assign(nv + 1, 1);
    for (const auto& [v, rf] : r.defs.defs)
      if (W.in_proj[v]) die("var " + std::to_string(v) + " is in the projection but defined: "
                            "its count is not a d-DNNF query");
    for (int l : assume) {
      if ((uint32_t)std::abs(l) > nv || !W.in_proj[(uint32_t)std::abs(l)])
        die("--assume " + std::to_string(l) + ": not a projection var");
      if (!W.weighted) W.w.clear();
      W.weighted = true;
      if (!W.w.count(l)) W.w[l] = 1;
      W.w[-l] = 0;
    }
    const mpq_class one = 1;
    auto w = [&](int l) -> mpq_class {
      const uint32_t v = (uint32_t)std::abs(l);
      if (!W.in_proj[v] || !W.weighted) return one;
      auto it = W.w.find(l);
      return it == W.w.end() ? one : it->second;
    };
    auto both = [&](uint32_t v) -> mpq_class { return w((int)v) + w(-(int)v); };

    const size_t n = (size_t)r.nnodes + 2;
    const size_t WD = (size_t)nv / 64 + 1;
    std::vector<mpq_class> val(n);
    std::vector<mpz_class> ucnt(samples ? n : 0);
    if (samples) { ucnt[0] = 0; ucnt[1] = 1; }
    auto pop_diff = [&](const uint64_t* a, const uint64_t* b) {
      unsigned long k = 0;
      for (size_t q = 0; q < WD; q++) k += (unsigned long)__builtin_popcountll(a[q] & ~b[q]);
      return k;
    };
    std::vector<uint64_t> vs(n * WD, 0);
    std::vector<char> dead(n, 0);
    val[0] = 0;
    val[1] = 1;
    dead[0] = 1;
    auto row = [&](uint64_t id) { return vs.data() + id * WD; };
    auto set = [](uint64_t* s, uint32_t v) { s[v >> 6] |= 1ULL << (v & 63); };
    auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
    std::vector<uint64_t> av;
    std::vector<std::vector<uint64_t>> arc_sets;
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      uint64_t* me = row(id);
      if (nd.kind == bdnnf::K_AND) {
        bool d = false;
        for (uint64_t k : nd.kids) d |= dead[k] != 0;
        for (size_t i = 0; i < nd.units.size() && !d; i++)
          for (size_t j = i + 1; j < nd.units.size(); j++) d |= nd.units[i] == -nd.units[j];
        if (d) { dead[id] = 1; val[id] = 0; if (samples) ucnt[id] = 0; return; }
        mpq_class x = 1;
        for (uint64_t k : nd.kids) {
          x *= val[k];
          const uint64_t* kv = row(k);
          for (size_t q = 0; q < WD; q++) me[q] |= kv[q];
        }
        for (int l : nd.units) {
          x *= w(l);
          if (W.in_proj[(uint32_t)std::abs(l)]) set(me, (uint32_t)std::abs(l));
        }
        for (uint32_t v : nd.free)
          if (W.in_proj[v]) { x *= both(v); set(me, v); }
        val[id] = std::move(x);
        if (samples) {
          mpz_class u = 1;
          for (uint64_t k : nd.kids) u *= ucnt[k];
          unsigned long nf = 0;
          for (uint32_t v : nd.free) nf += W.in_proj[v] ? 1 : 0;
          u <<= nf;
          ucnt[id] = std::move(u);
        }
        return;
      }
      arc_sets.clear();
      std::vector<size_t> live;
      for (size_t i = 0; i < nd.arcs.size(); i++) {
        const auto& a = nd.arcs[i];
        bool d = dead[a.child] != 0;
        for (size_t p = 0; p < a.lits.size() && !d; p++)
          for (size_t q = p + 1; q < a.lits.size(); q++) d |= a.lits[p] == -a.lits[q];
        if (d) continue;
        const uint64_t* cv = row(a.child);
        av.assign(cv, cv + WD);
        for (int l : a.lits)
          if (W.in_proj[(uint32_t)std::abs(l)]) set(av.data(), (uint32_t)std::abs(l));
        for (size_t q = 0; q < WD; q++) me[q] |= av[q];
        arc_sets.push_back(av);
        live.push_back(i);
      }
      if (live.empty()) {
        dead[id] = 1;
        std::fill(me, me + WD, 0);
        val[id] = 0;
        if (samples) ucnt[id] = 0;
        return;
      }
      if (samples) {
        mpz_class u = 0;
        for (size_t k = 0; k < live.size(); k++) {
          mpz_class y = ucnt[nd.arcs[live[k]].child];
          y <<= pop_diff(me, arc_sets[k].data());
          u += y;
        }
        ucnt[id] = std::move(u);
      }
      mpq_class x = 0;
      for (size_t k = 0; k < live.size(); k++) {
        const auto& a = nd.arcs[live[k]];
        mpq_class y = val[a.child];
        for (int l : a.lits) y *= w(l);
        for (size_t q = 0; q < WD; q++)
          for (uint64_t m = me[q] & ~arc_sets[k][q]; m; m &= m - 1)
            y *= both((uint32_t)(q * 64 + __builtin_ctzll(m)));
        x += y;
      }
      val[id] = std::move(x);
    });
    mpq_class res = val[r.root];
    if (!dead[r.root])
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] && !has(row(r.root), v)) res *= both(v);
    res.canonicalize();
    std::cout << res.get_str() << std::endl;
    if (samples && !dead[r.root]) {
      if (cnf.empty()) die("--sample needs --cnf");
      const auto cls = read_clauses(cnf);
      NodeAccess acc(r);
      bdnnf::Node nd;
      gmp_randclass rnd(gmp_randinit_default);
      rnd.seed(seed);
      std::vector<signed char> val2(nv + 1);
      std::vector<uint64_t> stk;
      std::vector<mpz_class> wts;
      for (uint64_t smp = 0; smp < samples; smp++) {
        std::fill(val2.begin(), val2.end(), 0);
        auto setl = [&](int l) {
          const signed char want = l > 0 ? 1 : -1;
          if (val2[std::abs(l)] == -want) die("FAIL: sampled path has conflicting lits");
          val2[std::abs(l)] = want;
        };
        stk.assign(1, r.root);
        std::vector<uint64_t> path;
        while (!stk.empty()) {
          const uint64_t x = stk.back();
          stk.pop_back();
          path.push_back(x);
          if (x <= 1) { if (x == 0) die("FAIL: sampled into FALSE"); continue; }
          acc.get(x, nd);
          if (nd.kind == bdnnf::K_AND) {
            for (int l : nd.units) setl(l);
            for (uint64_t k : nd.kids) stk.push_back(k);
            continue;
          }
          wts.clear();
          mpz_class tot = 0;
          for (const auto& a : nd.arcs) {
            bool d = dead[a.child] != 0;
            for (size_t p = 0; p < a.lits.size() && !d; p++)
              for (size_t q = p + 1; q < a.lits.size(); q++) d |= a.lits[p] == -a.lits[q];
            mpz_class y = 0;
            if (!d) {
              av.assign(row(a.child), row(a.child) + WD);
              for (int l : a.lits)
                if (W.in_proj[(uint32_t)std::abs(l)]) set(av.data(), (uint32_t)std::abs(l));
              y = ucnt[a.child];
              y <<= pop_diff(row(x), av.data());
            }
            tot += y;
            wts.push_back(std::move(y));
          }
          mpz_class pick = rnd.get_z_range(tot);
          size_t i = 0;
          while (pick >= wts[i]) { pick -= wts[i]; i++; }
          for (int l : nd.arcs[i].lits) setl(l);
          stk.push_back(nd.arcs[i].child);
        }
        std::vector<char> rnd_set(nv + 1, 0);
        std::vector<char> defined(nv + 1, 0);
        for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
        for (uint32_t v = 1; v <= nv; v++)
          if (!val2[v] && !defined[v]) { val2[v] = (rnd.get_z_bits(1) == 1) ? 1 : -1; rnd_set[v] = 1; }
        bdnnf::eval_defs(r.defs, nv, val2);
        for (const auto& cl : cls) {
          bool sat = false;
          for (int l : cl) sat |= (val2[std::abs(l)] == (l > 0 ? 1 : -1));
          if (sat) continue;
          std::ostringstream m;
          m << "FAIL: sample " << smp << " violates clause";
          for (int l : cl) m << " " << l << (rnd_set[std::abs(l)] ? "(not on path)" : "");
          m << "; nodes visited:";
          for (uint64_t x : path) m << " " << x;
          die(m.str());
        }
      }
      std::cerr << "ddnnf-count: " << samples << " sampled models satisfy the CNF" << std::endl;
    }
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
