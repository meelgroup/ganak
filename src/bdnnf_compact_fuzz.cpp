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

// Fuzzer of bdnnf compact: random d-DNNFs with dead nodes and (nested)
// duplicates, compacted, then checked: every reachable node keeps its value
// under random literal weights and its var scope, the output has no dead or
// identical nodes, compacting again changes no byte, and `bdnnf verify`
// accepts the input and the output.
// Usage: bdnnf-compact-fuzz [--num N] [--seed S] [--bdnnf PATH]   (no --num: forever)

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

#include "bdnnf_compact.hpp"

namespace {

constexpr uint64_t MOD = (1ULL << 61) - 1;
__extension__ typedef unsigned __int128 u128;
uint64_t mulm(uint64_t a, uint64_t b) {
  const u128 x = (u128)a * b;
  uint64_t r = (uint64_t)(x & MOD) + (uint64_t)(x >> 61);
  return r >= MOD ? r - MOD : r;
}
uint64_t addm(uint64_t a, uint64_t b) { return (a + b) % MOD; }

struct Gen {
  uint32_t nvars = 0;
  std::vector<bdnnf::Node> nodes;  // node id i + 2
  std::vector<uint64_t> scope;     // var mask per node id (0, 1 included)
  std::vector<uint64_t> cls;       // first node built with the same structure
  std::map<uint64_t, std::vector<uint64_t>> members;
  uint64_t root = 0, props = 0;
  bdnnf::KV hkv, fkv;
  bdnnf::Defs defs;
};

uint64_t vbit(int l) { return 1ULL << (std::abs(l) - 1); }

Gen generate(std::mt19937_64& rng) {
  auto rnd = [&](uint64_t n) { return n ? rng() % n : 0; };
  Gen g;
  g.nvars = 1 + (uint32_t)rnd(60);
  const uint32_t ndef = rnd(4) == 0 ? (uint32_t)rnd(std::min<uint32_t>(4, g.nvars)) : 0;
  const uint32_t cv = g.nvars - ndef;  // circuit vars 1..cv
  g.scope = {0, 0};
  g.cls = {0, 1};
  const uint64_t budget = 1 + rnd(rnd(4) == 0 ? 600 : 80);
  std::map<uint64_t, std::vector<uint64_t>> by_scope;
  auto vars_of = [](uint64_t m) {
    std::vector<int> vs;
    for (int v = 1; m; v++, m >>= 1) if (m & 1) vs.push_back(v);
    return vs;
  };
  auto sortl = [](std::vector<int32_t>& ls) {
    std::sort(ls.begin(), ls.end(), [](int a, int b) { return bdnnf::lit_code(a) < bdnnf::lit_code(b); });
  };
  auto emit = [&](const bdnnf::Node& nd, uint64_t sc, uint64_t c) {
    const uint64_t id = g.nodes.size() + 2;
    g.nodes.push_back(nd);
    g.scope.push_back(sc);
    g.cls.push_back(c == 0 ? id : c);
    g.members[g.cls.back()].push_back(id);
    by_scope[sc].push_back(id);
    return id;
  };
  // a node over (a subset of, if not smooth) the vars of S
  std::function<uint64_t(uint64_t, int)> build = [&](uint64_t S, int depth) -> uint64_t {
    auto it = by_scope.find(S);
    if (it != by_scope.end() && rnd(3) == 0) {
      const uint64_t src = it->second[rnd(it->second.size())];
      if (rnd(2)) return src;
      // a copy, its children replaced by copies of theirs, order shuffled
      bdnnf::Node nd = g.nodes[src - 2];
      auto swap_kid = [&](uint64_t& k) {
        if (k < 2 || rnd(2)) return;
        const auto& m = g.members[g.cls[k]];
        k = m[rnd(m.size())];
      };
      for (auto& k : nd.kids) swap_kid(k);
      for (auto& a : nd.arcs) swap_kid(a.child);
      std::shuffle(nd.kids.begin(), nd.kids.end(), rng);
      std::shuffle(nd.arcs.begin(), nd.arcs.end(), rng);
      return emit(nd, g.scope[src], g.cls[src]);
    }
    auto vs = vars_of(S);
    std::shuffle(vs.begin(), vs.end(), rng);
    bdnnf::Node nd;
    uint64_t sc = 0;
    if (S == 0 || depth > 12 || g.nodes.size() >= budget || rnd(6) == 0) {
      if (rnd(10) == 0) return rnd(2);
      nd.kind = bdnnf::K_AND;
      for (int v : vs) {
        if (rnd(8) == 0) continue;
        if (rnd(4) == 0) nd.free.push_back((uint32_t)v);
        else nd.units.push_back(rnd(2) ? v : -v);
        sc |= vbit(v);
      }
    } else if (rnd(2)) {
      nd.kind = bdnnf::K_AND;
      const uint64_t parts = 1 + rnd(4);
      std::vector<uint64_t> P(parts, 0);
      for (int v : vs) {
        const uint64_t r = rnd(parts + 1);
        if (r == parts) {
          if (rnd(3) == 0) nd.free.push_back((uint32_t)v);
          else nd.units.push_back(rnd(2) ? v : -v);
          sc |= vbit(v);
        } else P[r] |= vbit(v);
      }
      for (uint64_t part : P) {
        if (part == 0 && rnd(2)) continue;
        const uint64_t k = build(part, depth + 1);
        nd.kids.push_back(k);
        sc |= g.scope[k];
      }
    } else {
      nd.kind = bdnnf::K_OR;
      const size_t nd_vars = std::min<size_t>(vs.size(), 1 + rnd(2));
      const std::vector<int> D(vs.begin(), vs.begin() + (long)nd_vars);
      const std::vector<int> rest(vs.begin() + (long)nd_vars, vs.end());
      uint64_t Dm = 0;
      for (int v : D) Dm |= vbit(v);
      std::vector<uint32_t> asg(1U << D.size());
      for (uint32_t i = 0; i < asg.size(); i++) asg[i] = i;
      std::shuffle(asg.begin(), asg.end(), rng);
      asg.resize(1 + rnd(asg.size()));
      const bool smooth_arcs = rnd(2);
      for (uint32_t a : asg) {
        bdnnf::Arc arc;
        for (size_t j = 0; j < D.size(); j++) arc.lits.push_back((a >> j & 1) ? D[j] : -D[j]);
        uint64_t cs = S & ~Dm;
        for (int v : rest) {
          if (rnd(6) == 0) {
            arc.lits.push_back(rnd(2) ? v : -v);
            cs &= ~vbit(v);
          } else if (!smooth_arcs && rnd(6) == 0) cs &= ~vbit(v);
        }
        sortl(arc.lits);
        arc.child = build(cs, depth + 1);
        uint64_t s = g.scope[arc.child];
        for (int l : arc.lits) s |= vbit(l);
        sc |= s;
        nd.arcs.push_back(arc);
      }
    }
    sortl(nd.units);
    std::sort(nd.free.begin(), nd.free.end());
    return emit(nd, sc, 0);
  };
  const uint64_t all = cv == 64 ? ~0ULL : (1ULL << cv) - 1;
  auto random_scope = [&]() {
    uint64_t m = 0;
    for (uint32_t v = 1; v <= cv; v++) if (rnd(2)) m |= vbit((int)v);
    return m;
  };
  for (uint64_t i = 0, nd0 = rnd(3); i < nd0; i++) build(random_scope(), 6);
  g.root = build(rnd(2) ? all : random_scope(), 0);
  for (uint64_t i = 0, nd1 = rnd(3) == 0 ? rnd(3) : 0; i < nd1; i++) build(random_scope(), 6);
  if (rnd(10) == 0) g.root = rnd(2);
  const uint64_t last = g.nodes.size() + 1;

  // the claims: DECOMPOSABLE and DETERMINISTIC by construction; SMOOTH and
  // CLEAN when they happen to hold (bdnnf verify checks SMOOTH on dead nodes too)
  g.props = bdnnf::P_DECOMPOSABLE | bdnnf::P_DETERMINISTIC;
  std::vector<char> reach(g.nodes.size() + 2, 0);
  if (g.root >= 2) reach[g.root] = 1;
  bool smooth = ndef == 0 && g.scope[g.root] == (g.nvars == 64 ? ~0ULL : (1ULL << g.nvars) - 1);
  for (uint64_t id = last; id >= 2; id--) {
    const auto& nd = g.nodes[id - 2];
    if (reach[id]) {
      for (uint64_t k : nd.kids) reach[k] = 1;
      for (const auto& a : nd.arcs) reach[a.child] = 1;
    }
    uint64_t first = 0;
    bool have = false;
    for (const auto& a : nd.arcs) {
      uint64_t s = g.scope[a.child];
      for (int l : a.lits) s |= vbit(l);
      if (have && s != first) smooth = false;
      first = s;
      have = true;
    }
  }
  if (smooth && rnd(2)) g.props |= bdnnf::P_SMOOTH;
  bool clean = g.root == last;
  for (uint64_t id = 2; id <= last; id++) clean &= reach[id] != 0;
  if (clean && rnd(2)) g.props |= bdnnf::P_CLEAN;

  if (rnd(2)) g.hkv.emplace_back("producer", "bdnnf-compact-fuzz");
  if (rnd(2)) g.hkv.emplace_back("projection", "");
  if (rnd(2)) g.fkv.emplace_back("note", std::to_string(rng()));
  if (ndef) {
    auto aref = [&]() -> uint64_t {
      const uint64_t k = rnd(4) == 0 || g.defs.nodes.empty()
          ? (rnd(8) == 0 ? 0 : 1 + rnd(cv))
          : g.nvars + 1 + rnd(g.defs.nodes.size());
      return (k << 1) | rnd(2);
    };
    for (uint64_t i = 0, nn = rnd(6); i < nn; i++) g.defs.nodes.emplace_back(aref(), aref());
    for (uint32_t v = cv + 1; v <= g.nvars; v++) g.defs.defs.emplace_back(v, aref());
  }
  return g;
}

void write(const Gen& g, const std::string& path, std::mt19937_64& rng) {
  bdnnf::Writer w;
  w.zstd_level = (int)(rng() % 2);
  w.frame_target = 16 + rng() % 2048;
  w.open(path, g.nvars, g.hkv);
  for (const auto& nd : g.nodes) {
    if (nd.kind == bdnnf::K_AND) w.add_and(nd.kids, nd.units, nd.free);
    else {
      w.begin_or();
      for (const auto& a : nd.arcs) w.or_arc(a.lits, a.child);
      w.end_or();
    }
  }
  w.finish(g.root, g.props, g.fkv, &g.defs);
}

std::vector<bdnnf::Node> load(bdnnf::Reader& r) {
  std::vector<bdnnf::Node> v;
  r.scan([&](uint64_t, const bdnnf::Node& n) { v.push_back(n); });
  return v;
}

// per node: the arithmetic value under weights wt (literal code -> weight)
std::vector<uint64_t> values(const std::vector<bdnnf::Node>& ns, const std::vector<uint64_t>& wt) {
  std::vector<uint64_t> val{0, 1};
  for (const auto& nd : ns) {
    uint64_t x;
    if (nd.kind == bdnnf::K_AND) {
      x = 1;
      for (uint64_t k : nd.kids) x = mulm(x, val[k]);
      for (int l : nd.units) x = mulm(x, wt[bdnnf::lit_code(l)]);
      for (uint32_t v : nd.free) x = mulm(x, addm(wt[2 * v], wt[2 * v + 1]));
    } else {
      x = 0;
      for (const auto& a : nd.arcs) {
        uint64_t y = val[a.child];
        for (int l : a.lits) y = mulm(y, wt[bdnnf::lit_code(l)]);
        x = addm(x, y);
      }
    }
    val.push_back(x);
  }
  return val;
}

std::vector<uint64_t> scopes(const std::vector<bdnnf::Node>& ns) {
  std::vector<uint64_t> s{0, 0};
  for (const auto& nd : ns) {
    uint64_t m = 0;
    for (uint64_t k : nd.kids) m |= s[k];
    for (int l : nd.units) m |= vbit(l);
    for (uint32_t v : nd.free) m |= vbit((int)v);
    for (const auto& a : nd.arcs) {
      m |= s[a.child];
      for (int l : a.lits) m |= vbit(l);
    }
    s.push_back(m);
  }
  return s;
}

using Canon = std::tuple<uint8_t, std::vector<uint64_t>, std::vector<int32_t>, std::vector<uint32_t>,
                         std::vector<std::pair<std::vector<int32_t>, uint64_t>>>;
Canon canon(const bdnnf::Node& nd) {
  Canon c{nd.kind, nd.kids, nd.units, nd.free, {}};
  std::sort(std::get<1>(c).begin(), std::get<1>(c).end());
  for (const auto& a : nd.arcs) std::get<4>(c).emplace_back(a.lits, a.child);
  std::sort(std::get<4>(c).begin(), std::get<4>(c).end());
  return c;
}

std::string slurp(const std::string& p) {
  std::ifstream f(p, std::ios::binary);
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

}  // namespace

int main(int argc, char** argv) {
  uint64_t num = UINT64_MAX, seed = std::random_device{}();
  std::string exe = argv[0];
  std::string bdnnf_bin = exe.substr(0, exe.find_last_of('/') + 1) + "bdnnf";
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--num" && i + 1 < argc) num = std::stoull(argv[++i]);
    else if (a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
    else if (a == "--bdnnf" && i + 1 < argc) bdnnf_bin = argv[++i];
    else {
      std::cerr << "Usage: bdnnf-compact-fuzz [--num N] [--seed S] [--bdnnf PATH]\n";
      return EXIT_FAILURE;
    }
  }
  const char* td = getenv("TMPDIR");
  const std::string dir = td && *td ? td : "/tmp";
  const std::string tag = dir + "/compact_fuzz_" + std::to_string(getpid());
  const std::string fin = tag + "_in.bdnnf", fout = tag + "_out.bdnnf", fout2 = tag + "_out2.bdnnf";
  std::cout << "bdnnf-compact-fuzz: seed " << seed << std::endl;
  uint64_t in_nodes = 0, reached = 0, out_nodes = 0;
  for (uint64_t t = 0; t < num; t++) {
    const uint64_t s = seed + t;
    std::mt19937_64 rng(s);
    auto fail = [&](const std::string& msg) {
      std::cout << "FAIL (case seed " << s << "): " << msg << "\n  input kept in " << fin
                << "\n  rerun: " << argv[0] << " --num 1 --seed " << s << std::endl;
      exit(EXIT_FAILURE);
    };
    const Gen g = generate(rng);
    write(g, fin, rng);
    if (std::system((bdnnf_bin + " verify -q " + fin + " > /dev/null 2>&1").c_str()) != 0)
      fail("generator bug: bdnnf verify rejects the input");

    std::vector<uint64_t> map;
    bdnnf::CompactStats st;
    try {
      bdnnf::Reader r;
      r.open(fin);
      bdnnf::Writer w;
      w.open(fout, r.nvars, r.header_kv);
      st = bdnnf::compact(r, w, &map);
    } catch (const bdnnf::Error& e) {
      fail(std::string("compact threw: ") + e.what());
    }
    if (std::system((bdnnf_bin + " verify -q " + fout + " > /dev/null 2>&1").c_str()) != 0)
      fail("bdnnf verify rejects the output " + fout);

    bdnnf::Reader ri, ro;
    ri.open(fin);
    ro.open(fout);
    const auto ni = load(ri), no = load(ro);
    if (ro.nvars != ri.nvars || ro.header_kv != ri.header_kv || ro.footer_kv != ri.footer_kv)
      fail("header or footer changed");
    if (ro.props != (ri.props | bdnnf::P_CLEAN)) fail("props not the input's plus CLEAN");
    if (ro.defs.nodes != ri.defs.nodes || ro.defs.defs != ri.defs.defs) fail("definitions changed");
    if (no.size() != st.out_nodes) fail("written node count differs from the stats");

    std::vector<char> reach(ni.size() + 2, 0);
    if (ri.root >= 2) reach[ri.root] = 1;
    for (uint64_t id = ni.size() + 1; id >= 2; id--) {
      if (!reach[id]) continue;
      for (uint64_t k : ni[id - 2].kids) reach[k] = 1;
      for (const auto& a : ni[id - 2].arcs) reach[a.child] = 1;
    }
    if (map[ri.root] != ro.root) fail("root not mapped to the output root");
    const auto si = scopes(ni), so = scopes(no);
    for (int round = 0; round < 2; round++) {
      std::vector<uint64_t> wt(2 * (size_t)ri.nvars + 2);
      for (auto& x : wt) x = rng() % MOD;
      const auto vi = values(ni, wt), vo = values(no, wt);
      for (uint64_t id = 0; id < ni.size() + 2; id++) {
        if (id >= 2 && !reach[id]) continue;
        const uint64_t m = id < 2 ? id : map[id];
        if (m >= no.size() + 2) fail("node " + std::to_string(id) + " mapped out of range");
        if (vi[id] != vo[m]) fail("node " + std::to_string(id) + " -> " + std::to_string(m) + " changes value");
        if (si[id] != so[m]) fail("node " + std::to_string(id) + " -> " + std::to_string(m) + " changes scope");
      }
    }
    std::set<Canon> seen;
    for (size_t i = 0; i < no.size(); i++)
      if (!seen.insert(canon(no[i])).second) fail("output node " + std::to_string(i + 2) + " is a duplicate");

    try {
      bdnnf::Reader r;
      r.open(fout);
      bdnnf::Writer w;
      w.open(fout2, r.nvars, r.header_kv);
      bdnnf::compact(r, w);
    } catch (const bdnnf::Error& e) {
      fail(std::string("second compact threw: ") + e.what());
    }
    if (slurp(fout) != slurp(fout2)) fail("compacting the output again changes it");

    in_nodes += st.in_nodes;
    reached += st.reachable;
    out_nodes += st.out_nodes;
    if ((t + 1) % 1000 == 0)
      std::cout << (t + 1) << " cases OK, nodes: " << in_nodes << " in, " << reached << " reachable, "
                << out_nodes << " out" << std::endl;
  }
  std::cout << "bdnnf-compact-fuzz: " << num << " cases OK, nodes: " << in_nodes << " in, " << reached
            << " reachable, " << out_nodes << " out" << std::endl;
  for (const auto& p : {fin, fout, fout2}) unlink(p.c_str());
  return EXIT_SUCCESS;
}
