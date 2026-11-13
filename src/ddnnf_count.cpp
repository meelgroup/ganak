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
// counts, so uniformly when the circuit is decomposable and deterministic) and
// checks each against the CNF's clauses: a circuit too weak fails this, one too
// strong fails the count. Needs the graph in memory.
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
#include "ddnnf_io.hpp"

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

// Returns an error message, or "".
std::string sample_check(const std::string& path, const std::string& cnf, uint64_t n, uint64_t seed) {
  ddnnf_io::Circuit c;
  std::string err;
  if (!ddnnf_io::load(path, c, err)) return err;
  const uint32_t nv = c.nvars;
  const size_t WD = nv / 64 + 1;
  const size_t N = c.type.size();
  std::vector<std::vector<uint64_t>> vs(N, std::vector<uint64_t>(WD, 0));
  std::vector<mpz_class> cnt(N);
  std::vector<char> done(N, 0);
  auto pop = [&](const std::vector<uint64_t>& a) { int k = 0; for (auto x : a) k += __builtin_popcountll(x); return k; };
  std::vector<std::vector<uint64_t>> arcvs;
  ddnnf_io::post_order(c, done, [&](int x) {
    if (c.type[x] == 't') { cnt[x] = 1; return; }
    if (c.type[x] == 'f') { cnt[x] = 0; return; }
    auto& me = vs[x];
    if (c.type[x] == 'a') {
      cnt[x] = 1;
      for (const auto& a : c.arcs[x]) {
        cnt[x] *= cnt[a.child];
        for (size_t q = 0; q < WD; q++) me[q] |= vs[a.child][q];
        for (int l : a.lits) me[std::abs(l) >> 6] |= 1ULL << (std::abs(l) & 63);
      }
    } else {
      cnt[x] = 0;
      arcvs.assign(c.arcs[x].size(), std::vector<uint64_t>(WD));
      for (size_t i = 0; i < c.arcs[x].size(); i++) {
        const auto& a = c.arcs[x][i];
        arcvs[i] = vs[a.child];
        for (int l : a.lits) arcvs[i][std::abs(l) >> 6] |= 1ULL << (std::abs(l) & 63);
        for (size_t q = 0; q < WD; q++) me[q] |= arcvs[i][q];
      }
      for (size_t i = 0; i < c.arcs[x].size(); i++) {
        mpz_class y = cnt[c.arcs[x][i].child];
        y <<= (unsigned long)(pop(me) - pop(arcvs[i]));
        cnt[x] += y;
      }
    }
  });
  if (cnt[c.root] == 0) return "";
  const auto cls = read_clauses(cnf);
  gmp_randclass rnd(gmp_randinit_default);
  rnd.seed(seed);
  std::vector<signed char> val(nv + 1);
  for (uint64_t s = 0; s < n; s++) {
    std::fill(val.begin(), val.end(), 0);
    std::vector<int> st{c.root};
    while (!st.empty()) {
      const int x = st.back();
      st.pop_back();
      auto setl = [&](int l) {
        const signed char want = l > 0 ? 1 : -1;
        if (val[std::abs(l)] == -want) return false;
        val[std::abs(l)] = want;
        return true;
      };
      if (c.type[x] == 'a') {
        for (const auto& a : c.arcs[x]) {
          for (int l : a.lits) if (!setl(l)) return "sampled model has conflicting lits";
          st.push_back(a.child);
        }
      } else if (c.type[x] == 'o') {
        std::vector<mpz_class> wts;
        mpz_class tot = 0;
        for (const auto& a : c.arcs[x]) {
          std::vector<uint64_t> av = vs[a.child];
          for (int l : a.lits) av[std::abs(l) >> 6] |= 1ULL << (std::abs(l) & 63);
          mpz_class y = cnt[a.child];
          y <<= (unsigned long)(pop(vs[x]) - pop(av));
          wts.push_back(y);
          tot += y;
        }
        mpz_class pick = rnd.get_z_range(tot);
        size_t i = 0;
        while (pick >= wts[i]) { pick -= wts[i]; i++; }
        for (int l : c.arcs[x][i].lits) if (!setl(l)) return "sampled model has conflicting lits";
        st.push_back(c.arcs[x][i].child);
      } else if (c.type[x] == 'f') {
        return "sampled into FALSE";
      }
    }
    std::vector<char> rnd_set(nv + 1, 0);
    for (uint32_t v = 1; v <= nv; v++)
      if (!val[v]) { val[v] = (rnd.get_z_bits(1) == 1) ? 1 : -1; rnd_set[v] = 1; }
    for (const auto& cl : cls) {
      bool sat = false;
      for (int l : cl) sat |= (val[std::abs(l)] == (l > 0 ? 1 : -1));
      if (!sat) {
        std::ostringstream m;
        m << "sample " << s << " violates clause";
        for (int l : cl) m << " " << l << (rnd_set[std::abs(l)] ? "(not on path)" : "");
        m << "; model:";
        for (uint32_t v = 1; v <= nv; v++) m << " " << (val[v] > 0 ? (int)v : -(int)v);
        return m.str();
      }
    }
  }
  return "";
}

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
        if (d) { dead[id] = 1; val[id] = 0; return; }
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
        return;
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
    if (samples) {
      if (cnf.empty()) die("--sample needs --cnf");
      const std::string e = sample_check(path, cnf, samples, seed);
      if (!e.empty()) die("FAIL: " + e);
      std::cerr << "ddnnf-count: " << samples << " sampled models satisfy the CNF" << std::endl;
    }
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
