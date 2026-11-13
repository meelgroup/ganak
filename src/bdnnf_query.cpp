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

// bdnnf count / bdnnf sample: the exact (weighted, projected) model count of a
// BDNNF circuit in one streamed pass, one var bitset and one rational per node.
// Smooths on the fly, so it needs a decomposable, deterministic circuit but not
// a SMOOTH one. Weights ('c p weight', Arjun's semantics) come from --cnf; the
// projection is the file's "projection" header key (a --cnf must agree).
// --assume L1,L2,... conditions on literals of projection vars.
//
// sample draws N models, each with probability proportional to its weight
// (uniformly when unweighted) over the projection, top-down through the frame
// index; non-projection vars come from the path, the definitions, or at
// random. Models are printed as DIMACS lines over the projection (--full: all
// vars), to stdout or -o. With --cnf, each is checked against the clauses.

#include <algorithm>
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
  std::cerr << "bdnnf: " << msg << std::endl;
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

// a * (product of fs except fs[i]) for every i, without dividing.
void excl_products(const std::vector<mpq_class>& fs, const mpq_class& a, std::vector<mpq_class>& out) {
  const size_t k = fs.size();
  out.assign(k, a);
  mpq_class acc = 1;
  for (size_t i = 0; i < k; i++) { out[i] *= acc; acc *= fs[i]; }
  acc = 1;
  for (size_t i = k; i-- > 0;) { out[i] *= acc; acc *= fs[i]; }
}

// Weighted count of the models containing each projection literal (index
// 2v for v, 2v+1 for -v): the derivatives of the count, top-down in one
// backward scan. A literal's count is w(l) times the derivative by w(l).
template <class WFn>
std::vector<mpq_class> lit_counts(bdnnf::Reader& r, const std::vector<mpq_class>& val,
                                  const std::vector<char>& dead, const std::vector<uint64_t>& vs,
                                  size_t WD, const std::vector<char>& in_proj, WFn&& w) {
  const uint32_t nv = r.nvars;
  std::vector<mpq_class> cnt(2 * (size_t)nv + 2), dv(val.size());
  auto row = [&](uint64_t id) { return vs.data() + id * WD; };
  auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  auto both = [&](uint32_t v) -> mpq_class { return w((int)v) + w(-(int)v); };
  std::vector<mpq_class> fs, ex;
  std::vector<uint32_t> gap;
  // Each var of `gap` is a factor (w(v) + w(-v)) at fs[off + i].
  auto add_gap = [&](size_t off) {
    for (size_t i = 0; i < gap.size(); i++) {
      cnt[lidx((int)gap[i])] += ex[off + i] * w((int)gap[i]);
      cnt[lidx(-(int)gap[i])] += ex[off + i] * w(-(int)gap[i]);
    }
  };
  if (dead[r.root]) return cnt;
  fs.assign(1, val[r.root]);
  gap.clear();
  for (uint32_t v = 1; v <= nv; v++)
    if (in_proj[v] && !has(row(r.root), v)) { gap.push_back(v); fs.push_back(both(v)); }
  excl_products(fs, 1, ex);
  dv[r.root] = ex[0];
  add_gap(1);
  std::vector<uint64_t> av(WD);
  r.scan_backward([&](uint64_t id, const bdnnf::Node& nd) {
    if (dead[id] || dv[id] == 0) return;
    const mpq_class& g = dv[id];
    if (nd.kind == bdnnf::K_AND) {
      fs.clear();
      for (uint64_t k : nd.kids) fs.push_back(val[k]);
      gap.clear();
      for (uint32_t v : nd.free) if (in_proj[v]) { gap.push_back(v); fs.push_back(both(v)); }
      mpq_class u = 1;
      for (int l : nd.units) u *= w(l);
      excl_products(fs, g * u, ex);
      for (size_t i = 0; i < nd.kids.size(); i++) dv[nd.kids[i]] += ex[i];
      add_gap(nd.kids.size());
      for (int l : nd.units) if (in_proj[(uint32_t)std::abs(l)]) cnt[lidx(l)] += g * val[id];
      return;
    }
    const uint64_t* me = row(id);
    for (const auto& a : nd.arcs) {
      if (dead[a.child]) continue;
      bool conflict = false;
      for (size_t p = 0; p < a.lits.size(); p++)
        for (size_t q = p + 1; q < a.lits.size(); q++) conflict |= a.lits[p] == -a.lits[q];
      if (conflict) continue;
      const uint64_t* cv = row(a.child);
      av.assign(cv, cv + WD);
      mpq_class u = 1;
      for (int l : a.lits) {
        u *= w(l);
        const uint32_t v = (uint32_t)std::abs(l);
        if (in_proj[v]) av[v >> 6] |= 1ULL << (v & 63);
      }
      fs.assign(1, val[a.child]);
      gap.clear();
      for (size_t q = 0; q < WD; q++)
        for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1) {
          gap.push_back((uint32_t)(q * 64 + __builtin_ctzll(m)));
          fs.push_back(both(gap.back()));
        }
      excl_products(fs, g * u, ex);
      dv[a.child] += ex[0];
      add_gap(1);
      const mpq_class y = ex[0] * val[a.child];
      for (int l : a.lits) if (in_proj[(uint32_t)std::abs(l)]) cnt[lidx(l)] += y;
    }
  });
  return cnt;
}

std::vector<std::vector<int>> read_queries(const std::string& path, const std::vector<char>& in_proj) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open " + path);
  std::vector<std::vector<int>> qs;
  std::string line;
  while (std::getline(f, line)) {
    if (!line.empty() && line[0] == 'c') continue;
    std::istringstream ss(line);
    std::vector<int> q;
    long l;
    while (ss >> l && l != 0) {
      if ((size_t)std::abs(l) >= in_proj.size() || !in_proj[(size_t)std::abs(l)])
        die("query literal " + std::to_string(l) + " is not of a projection var");
      q.push_back((int)l);
    }
    qs.push_back(std::move(q));
  }
  return qs;
}

// Answers the queries (assumption sets) in chunks, one forward scan per chunk:
// --sat with one bit per query and node, 64 per scan; counts with one
// rational per query and node, as many per scan as fit kMaxVals.
template <class WFn>
void batch(bdnnf::Reader& r, const std::vector<char>& dead, const std::vector<uint64_t>& vs,
           size_t WD, const std::vector<char>& in_proj, WFn&& w,
           const std::vector<std::vector<int>>& qs, bool sat, const std::vector<int>& assume,
           std::vector<std::string>& out) {
  constexpr size_t kMaxVals = 1 << 22;
  const uint32_t nv = r.nvars;
  const size_t n = dead.size();
  auto row = [&](uint64_t id) { return vs.data() + id * WD; };
  auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
  auto conflicting = [](const std::vector<int>& ls) {
    for (size_t p = 0; p < ls.size(); p++)
      for (size_t q = p + 1; q < ls.size(); q++) if (ls[p] == -ls[q]) return true;
    return false;
  };
  const size_t K = sat ? 64 : std::max<size_t>(1, std::min<size_t>(64, kMaxVals / n));
  // kill[lit]: the queries of the chunk that assume -lit
  std::vector<uint64_t> kill(2 * (size_t)nv + 2);
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  std::vector<uint64_t> sm(sat ? n : 0);
  std::vector<mpq_class> val(sat ? 0 : n * K);
  std::vector<uint32_t> gap;
  std::vector<uint64_t> av(WD);
  for (size_t q0 = 0; q0 < qs.size(); q0 += K) {
    const size_t k = std::min(K, qs.size() - q0);
    const uint64_t all = k == 64 ? ~0ULL : (1ULL << k) - 1;
    std::fill(kill.begin(), kill.end(), 0);
    for (size_t j = 0; j < k; j++)
      for (int l : qs[q0 + j]) kill[lidx(-l)] |= 1ULL << j;
    if (sat) for (int l : assume) kill[lidx(-l)] |= all;
    auto wj = [&](int l, size_t j) -> mpq_class {
      if (!in_proj[(uint32_t)std::abs(l)]) return 1;
      return (kill[lidx(l)] >> j) & 1 ? mpq_class(0) : w(l);
    };
    auto bothj = [&](uint32_t v, size_t j) -> mpq_class { return wj((int)v, j) + wj(-(int)v, j); };
    auto killed = [&](const std::vector<int>& ls) {
      uint64_t m = 0;
      for (int l : ls) if (in_proj[(uint32_t)std::abs(l)]) m |= kill[lidx(l)];
      return m;
    };
    if (sat) { sm[0] = 0; sm[1] = all; }
    else for (size_t j = 0; j < K; j++) { val[j] = 0; val[K + j] = 1; }
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      if (dead[id]) {
        if (sat) sm[id] = 0;
        else for (size_t j = 0; j < k; j++) val[id * K + j] = 0;
        return;
      }
      if (nd.kind == bdnnf::K_AND) {
        if (sat) {
          uint64_t m = all & ~killed(nd.units);
          for (uint64_t c : nd.kids) m &= sm[c];
          sm[id] = m;
          return;
        }
        for (size_t j = 0; j < k; j++) {
          mpq_class x = 1;
          for (uint64_t c : nd.kids) x *= val[c * K + j];
          for (int l : nd.units) x *= wj(l, j);
          for (uint32_t v : nd.free) if (in_proj[v]) x *= bothj(v, j);
          val[id * K + j] = std::move(x);
        }
        return;
      }
      if (sat) {
        uint64_t m = 0;
        for (const auto& a : nd.arcs)
          if (!dead[a.child] && !conflicting(a.lits)) m |= sm[a.child] & ~killed(a.lits);
        sm[id] = m;
        return;
      }
      const uint64_t* me = row(id);
      for (size_t j = 0; j < k; j++) val[id * K + j] = 0;
      for (const auto& a : nd.arcs) {
        if (dead[a.child] || conflicting(a.lits)) continue;
        const uint64_t* cv = row(a.child);
        av.assign(cv, cv + WD);
        for (int l : a.lits) {
          const uint32_t v = (uint32_t)std::abs(l);
          if (in_proj[v]) av[v >> 6] |= 1ULL << (v & 63);
        }
        gap.clear();
        for (size_t q = 0; q < WD; q++)
          for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1) gap.push_back((uint32_t)(q * 64 + __builtin_ctzll(m)));
        for (size_t j = 0; j < k; j++) {
          mpq_class y = val[a.child * K + j];
          for (int l : a.lits) y *= wj(l, j);
          for (uint32_t v : gap) y *= bothj(v, j);
          val[id * K + j] += y;
        }
      }
    });
    for (size_t j = 0; j < k; j++) {
      if (sat) {
        const bool ok = !dead[r.root] && ((sm[r.root] >> j) & 1);
        bool conflict = false;
        for (int l : qs[q0 + j]) conflict |= (kill[lidx(l)] >> j) & 1;
        for (int l : assume) conflict |= (kill[lidx(l)] >> j) & 1;
        out.push_back(ok && !conflict ? "SAT" : "UNSAT");
        continue;
      }
      mpq_class x = dead[r.root] ? mpq_class(0) : val[r.root * K + j];
      if (!dead[r.root])
        for (uint32_t v = 1; v <= nv; v++)
          if (in_proj[v] && !has(row(r.root), v)) x *= bothj(v, j);
      x.canonicalize();
      out.push_back(x.get_str());
    }
  }
}

}  // namespace

int bdnnf_query(int argc, char** argv) {
  const std::string cmd = argv[0];
  const bool enumerating = cmd == "enum";
  const bool atomic = cmd == "atomic";
  const bool sampling = cmd == "sample" || enumerating || atomic;
  const bool lits = cmd == "lits" || cmd == "core";
  const bool core = cmd == "core";
  const bool is_batch = cmd == "batch";
  const std::string usage = atomic
    ? "usage: bdnnf atomic <in.bdnnf> [--cross] [--cnf <file.cnf>] [--assume L1,L2,...]"
    : enumerating
    ? "usage: bdnnf enum <in.bdnnf> [-n N] [--start I] [--cnf <file.cnf>] [--assume L1,L2,...] [--full] [-o out]"
    : sampling
    ? "usage: bdnnf sample <in.bdnnf> [-n N] [--seed S] [--cnf <file.cnf>] [--assume L1,L2,...] [--full] [-o out]"
    : is_batch
    ? "usage: bdnnf batch <in.bdnnf> <queries> [--sat] [--cnf <file.cnf>]   (a query: a line of literals)"
    : "usage: bdnnf " + cmd + " <in.bdnnf> [--cnf <file.cnf>] [--assume L1,L2,...]";
  std::string path, cnf, out_path;
  uint64_t num = enumerating ? UINT64_MAX : atomic ? 64 : 1, seed = 1;
  mpz_class start = 0;
  bool full = false, sat = false, cross = false;
  std::string queries;
  std::vector<int> assume;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--cnf" && i + 1 < argc) cnf = argv[++i];
    else if (sampling && (a == "-n" || a == "--num") && i + 1 < argc) num = std::stoull(argv[++i]);
    else if (sampling && !enumerating && a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
    else if (enumerating && a == "--start" && i + 1 < argc) start = mpz_class(argv[++i], 10);
    else if (sampling && a == "--full") full = true;
    else if (sampling && a == "-o" && i + 1 < argc) out_path = argv[++i];
    else if (a == "--assume" && i + 1 < argc) {
      std::string t = argv[++i];
      for (char& ch : t) if (ch == ',') ch = ' ';
      std::istringstream ss(t);
      long l;
      while (ss >> l) if (l != 0) assume.push_back((int)l);
    }
    else if (is_batch && a == "--sat") sat = true;
    else if (atomic && a == "--cross") cross = true;
    else if (path.empty() && (a == "-" || a.empty() || a[0] != '-')) path = a;
    else if (is_batch && queries.empty() && (a.empty() || a[0] != '-')) queries = a;
    else die(usage);
  }
  if (path.empty() || (is_batch && queries.empty())) die(usage);
  try {
    bdnnf::Reader r;
    r.open(path);
    const uint32_t nv = r.nvars;
    std::vector<char> file_proj(nv + 1, 1);
    const std::string hp = r.get(r.header_kv, "projection");
    if (!hp.empty()) {
      file_proj.assign(nv + 1, 0);
      for (uint32_t v : bdnnf::parse_var_list(hp)) {
        if (v == 0 || v > nv) die("projection var out of range: " + std::to_string(v));
        file_proj[v] = 1;
      }
    }
    Weights W;
    if (!cnf.empty()) {
      W = read_cnf(cnf, nv);
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] != file_proj[v])
          die("the CNF's projection differs from the circuit's (var " + std::to_string(v) + ")");
    } else {
      W.projected = !hp.empty();
      W.in_proj = file_proj;
    }
    for (const auto& [v, rf] : r.defs.defs)
      if (W.in_proj[v]) die("var " + std::to_string(v) + " is in the projection but defined: "
                            "its count is not a d-DNNF query");
    if (core || enumerating || atomic) {
      W.weighted = false;
      W.w.clear();
    }
    for (int l : assume) {
      if ((uint32_t)std::abs(l) > nv || !W.in_proj[(uint32_t)std::abs(l)])
        die("--assume " + std::to_string(l) + ": not a projection var");
      if (!W.weighted) W.w.clear();
      W.weighted = true;
      if (!W.w.count(l)) W.w[l] = 1;
      W.w[-l] = 0;
    }
    if (sampling && !enumerating)
      for (const auto& [l, x] : W.w)
        if (x < 0 && W.in_proj[(uint32_t)std::abs(l)])
          die("negative weight on literal " + std::to_string(l) + ": cannot sample");
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
    auto conflicting = [](const std::vector<int>& ls) {
      for (size_t p = 0; p < ls.size(); p++)
        for (size_t q = p + 1; q < ls.size(); q++) if (ls[p] == -ls[q]) return true;
      return false;
    };
    std::vector<uint64_t> av;
    // weight of an OR arc, smoothed over the vars of its node `me`
    auto arc_val = [&](const uint64_t* me, const bdnnf::Arc& a) -> mpq_class {
      if (dead[a.child] || conflicting(a.lits)) return 0;
      const uint64_t* cv = row(a.child);
      av.assign(cv, cv + WD);
      mpq_class y = val[a.child];
      for (int l : a.lits) {
        y *= w(l);
        if (W.in_proj[(uint32_t)std::abs(l)]) set(av.data(), (uint32_t)std::abs(l));
      }
      for (size_t q = 0; q < WD; q++)
        for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1)
          y *= both((uint32_t)(q * 64 + __builtin_ctzll(m)));
      return y;
    };
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      uint64_t* me = row(id);
      if (nd.kind == bdnnf::K_AND) {
        bool d = conflicting(nd.units);
        for (uint64_t k : nd.kids) d |= dead[k] != 0;
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
      bool any = false;
      for (const auto& a : nd.arcs) {
        if (dead[a.child] || conflicting(a.lits)) continue;
        any = true;
        const uint64_t* cv = row(a.child);
        for (size_t q = 0; q < WD; q++) me[q] |= cv[q];
        for (int l : a.lits)
          if (W.in_proj[(uint32_t)std::abs(l)]) set(me, (uint32_t)std::abs(l));
      }
      if (!any) { dead[id] = 1; val[id] = 0; return; }
      mpq_class x = 0;
      for (const auto& a : nd.arcs) x += arc_val(me, a);
      val[id] = std::move(x);
    });
    mpq_class res = val[r.root];
    if (!dead[r.root])
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] && !has(row(r.root), v)) res *= both(v);
    res.canonicalize();
    if (is_batch) {
      std::vector<std::string> ans;
      batch(r, dead, vs, WD, W.in_proj, w, read_queries(queries, W.in_proj), sat, assume, ans);
      std::string o;
      for (const auto& a : ans) { o += a; o += '\n'; }
      std::cout << o << std::flush;
      return EXIT_SUCCESS;
    }
    if (lits) {
      const auto cnt = lit_counts(r, val, dead, vs, WD, W.in_proj, w);
      std::string line;
      if (core) {
        if (res == 0) die("no model");
        for (uint32_t v = 1; v <= nv; v++) {
          if (!W.in_proj[v]) continue;
          if (cnt[2 * v] == res) line += std::to_string(v) + " ";
          else if (cnt[2 * v + 1] == res) line += "-" + std::to_string(v) + " ";
        }
        std::cout << line << "0" << std::endl;
        return EXIT_SUCCESS;
      }
      for (uint32_t v = 1; v <= nv; v++) {
        if (!W.in_proj[v]) continue;
        mpq_class a = cnt[2 * v], b = cnt[2 * v + 1];
        a.canonicalize();
        b.canonicalize();
        line = std::to_string(v) + " " + a.get_str() + "\n-" + std::to_string(v) + " " + b.get_str() + "\n";
        std::cout << line;
      }
      std::cout << std::flush;
      return EXIT_SUCCESS;
    }
    if (!sampling) {
      std::cout << res.get_str() << std::endl;
      return EXIT_SUCCESS;
    }
    if (res == 0 && !enumerating) die(atomic ? "no model" : "no model of non-zero weight");
    if (enumerating) {
      const mpz_class left = start < res.get_num() ? mpz_class(res.get_num() - start) : mpz_class(0);
      if (left < num) num = left.get_ui();
    }

    const auto cls = cnf.empty() ? std::vector<std::vector<int>>{} : read_clauses(cnf);
    NodeAccess acc(r);
    bdnnf::Node nd;
    gmp_randclass rnd(gmp_randinit_default);
    rnd.seed(seed);
    std::vector<mpq_class> ys;
    std::vector<mpz_class> zs;
    // index i with probability ys[i] / sum(ys)
    auto pick = [&]() -> size_t {
      mpz_class den = 1, tot = 0;
      for (const auto& y : ys) mpz_lcm(den.get_mpz_t(), den.get_mpz_t(), y.get_den_mpz_t());
      zs.clear();
      for (const auto& y : ys) {
        zs.push_back(y.get_num() * (den / y.get_den()));
        tot += zs.back();
      }
      if (tot <= 0) die("sampled into a node of weight 0");
      mpz_class p = rnd.get_z_range(tot);
      size_t i = 0;
      while (p >= zs[i]) { p -= zs[i]; i++; }
      return i;
    };
    std::vector<char> defined(nv + 1, 0);
    for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
    std::vector<signed char> val2(nv + 1);
    std::vector<char> rnd_set(nv + 1);
    std::vector<uint64_t> stk, path_nodes;
    std::vector<mpz_class> istk;
    std::vector<uint64_t> cav(WD);
    std::vector<uint64_t> sig(atomic ? nv + 1 : 0);  // atomic: var values in the samples
    // Unranking: the models of a node are numbered in mixed radix over its
    // factors (children, then free vars, each var a digit of radix w(v)+w(-v)).
    auto digit = [](mpz_class& i, const mpz_class& radix) {
      mpz_class d;
      mpz_fdiv_qr(i.get_mpz_t(), d.get_mpz_t(), i.get_mpz_t(), radix.get_mpz_t());
      return d;
    };
    std::ofstream of;
    if (!out_path.empty() && out_path != "-") {
      of.open(out_path);
      if (!of) die("cannot open " + out_path);
    }
    std::ostream& out = of.is_open() ? of : std::cout;
    std::string line;
    for (uint64_t smp = 0; smp < num; smp++) {
      std::fill(val2.begin(), val2.end(), 0);
      std::fill(rnd_set.begin(), rnd_set.end(), 0);
      auto setl = [&](int l) {
        const signed char want = l > 0 ? 1 : -1;
        if (val2[std::abs(l)] == -want) die("FAIL: sampled path has conflicting lits");
        val2[std::abs(l)] = want;
      };
      // a free var of the unranked model: radix 2, or 1 when one literal has weight 0
      auto set_free = [&](uint32_t v, mpz_class& i) {
        const mpz_class rad = both(v).get_num();
        setl(rad == 1 ? (w((int)v) == 0 ? -(int)v : (int)v) : (digit(i, rad) == 1 ? (int)v : -(int)v));
      };
      stk.assign(1, r.root);
      istk.clear();
      if (enumerating) {
        mpz_class i = start + smp;
        for (uint32_t v = 1; v <= nv; v++)
          if (W.in_proj[v] && !has(row(r.root), v)) set_free(v, i);
        istk.push_back(i);
      }
      path_nodes.clear();
      while (!stk.empty()) {
        const uint64_t x = stk.back();
        stk.pop_back();
        mpz_class i;
        if (enumerating) { i = istk.back(); istk.pop_back(); }
        path_nodes.push_back(x);
        if (x <= 1) { if (x == 0) die("FAIL: sampled into FALSE"); continue; }
        acc.get(x, nd);
        if (nd.kind == bdnnf::K_AND) {
          for (int l : nd.units) setl(l);
          for (uint64_t k : nd.kids) {
            stk.push_back(k);
            if (enumerating) istk.push_back(digit(i, val[k].get_num()));
          }
          if (enumerating) for (uint32_t v : nd.free) if (W.in_proj[v]) set_free(v, i);
          continue;
        }
        ys.clear();
        for (const auto& a : nd.arcs) ys.push_back(arc_val(row(x), a));
        size_t ai = 0;
        if (!enumerating) ai = pick();
        else
          while (i >= ys[ai].get_num()) { i -= ys[ai].get_num(); ai++; }
        const auto& a = nd.arcs[ai];
        for (int l : a.lits) setl(l);
        stk.push_back(a.child);
        if (enumerating) {
          istk.push_back(digit(i, val[a.child].get_num()));
          const uint64_t* cv = row(a.child);
          cav.assign(cv, cv + WD);
          for (int l : a.lits) {
            const uint32_t v = (uint32_t)std::abs(l);
            if (W.in_proj[v]) set(cav.data(), v);
          }
          const uint64_t* me = row(x);
          for (size_t q = 0; q < WD; q++)
            for (uint64_t m = me[q] & ~cav[q]; m; m &= m - 1) set_free((uint32_t)(q * 64 + __builtin_ctzll(m)), i);
        }
      }
      // Projection vars off the path are free: each is drawn by its weight.
      for (uint32_t v = 1; v <= nv; v++) {
        if (val2[v] || defined[v]) continue;
        rnd_set[v] = 1;
        if (enumerating) {
          if (W.in_proj[v]) die("FAIL: model " + mpz_class(start + smp).get_str() + " leaves var " + std::to_string(v) + " open");
          val2[v] = -1;
        } else if (W.in_proj[v]) {
          ys.assign({w(-(int)v), w((int)v)});
          val2[v] = pick() == 1 ? 1 : -1;
        } else val2[v] = rnd.get_z_bits(1) == 1 ? 1 : -1;
      }
      bdnnf::eval_defs(r.defs, nv, val2);
      for (const auto& cl : cls) {
        bool ok = false;
        for (int l : cl) ok |= (val2[std::abs(l)] == (l > 0 ? 1 : -1));
        if (ok) continue;
        std::ostringstream m;
        m << "FAIL: sample " << smp << " violates clause";
        for (int l : cl) m << " " << l << (rnd_set[std::abs(l)] ? "(not on path)" : "");
        m << "; nodes visited:";
        for (uint64_t x : path_nodes) m << " " << x;
        die(m.str());
      }
      if (atomic) {
        for (uint32_t v = 1; v <= nv; v++) if (val2[v] > 0) sig[v] |= 1ULL << smp;
        continue;
      }
      line.clear();
      for (uint32_t v = 1; v <= nv; v++) {
        if (!full && !W.in_proj[v]) continue;
        if (val2[v] < 0) line += '-';
        line += std::to_string(v);
        line += ' ';
      }
      line += "0\n";
      out << line;
    }
    if (atomic) {
      const auto cnt = lit_counts(r, val, dead, vs, WD, W.in_proj, w);
      // Candidates: same literal count, same values in every sample. With
      // --cross each var is taken by its literal true in the first sample.
      std::map<std::pair<std::string, uint64_t>, std::vector<int>> groups;
      for (uint32_t v = 1; v <= nv; v++) {
        if (!W.in_proj[v]) continue;
        const int l = cross && !(sig[v] & 1) ? -(int)v : (int)v;
        const uint64_t sg = l > 0 ? sig[v] : ~sig[v];
        groups[{mpq_class(cnt[2 * v + (l < 0 ? 1 : 0)]).get_str(), sg}].push_back(l);
      }
      std::vector<std::vector<int>> cand, sets;
      for (auto& [k, g] : groups) if (g.size() >= 2) cand.push_back(std::move(g));
      // Confirm each member y against the first, r: {r, -y} and {-r, y} UNSAT.
      while (!cand.empty()) {
        std::vector<std::vector<int>> qs, next;
        for (const auto& g : cand)
          for (size_t i = 1; i < g.size(); i++) {
            qs.push_back({g[0], -g[i]});
            qs.push_back({-g[0], g[i]});
          }
        std::vector<std::string> ans;
        batch(r, dead, vs, WD, W.in_proj, w, qs, true, assume, ans);
        size_t qi = 0;
        for (const auto& g : cand) {
          std::vector<int> eq{g[0]}, rest;
          for (size_t i = 1; i < g.size(); i++, qi += 2)
            (ans[qi] == "UNSAT" && ans[qi + 1] == "UNSAT" ? eq : rest).push_back(g[i]);
          if (eq.size() >= 2) sets.push_back(std::move(eq));
          if (rest.size() >= 2) next.push_back(std::move(rest));
        }
        cand = std::move(next);
      }
      for (auto& st : sets) {
        if (st[0] < 0) for (int& l : st) l = -l;
        std::sort(st.begin(), st.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
      }
      std::sort(sets.begin(), sets.end(), [](const auto& a, const auto& b) { return std::abs(a[0]) < std::abs(b[0]); });
      for (const auto& st : sets) {
        line.clear();
        for (int l : st) line += std::to_string(l) + " ";
        out << line << "0\n";
      }
    }
    out.flush();
    if (!out) die("write error");
    if (!cls.empty() && !atomic)
      std::cerr << "bdnnf " << cmd << ": " << num << " models satisfy the CNF" << std::endl;
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}
