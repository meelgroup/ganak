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

// Compaction of a BDNNF circuit: drops the nodes the root does not reach and
// merges structurally identical ones (the same kind, units, free vars, and,
// after renumbering, the same children and arcs in any order). Both keep every
// property, so the output claims the input's flags plus CLEAN.

#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "bdnnf.hpp"

namespace bdnnf {

struct CompactStats {
  uint64_t in_nodes = 0, reachable = 0, out_nodes = 0;
};

// old_to_new (if given) gets, for every reachable node, its id in the output.
inline CompactStats compact(Reader& r, Writer& w, std::vector<uint64_t>* old_to_new = nullptr) {
  CompactStats st;
  st.in_nodes = r.nnodes;
  std::vector<char> reach(r.nnodes + 2, 0);
  if (r.root >= 2) reach[r.root] = 1;
  r.scan_backward([&](uint64_t id, const Node& n) {
    if (!reach[id]) return;
    st.reachable++;
    if (n.kind == K_AND) for (uint64_t c : n.kids) reach[c] = 1;
    else for (const auto& a : n.arcs) reach[a.child] = 1;
  });
  std::vector<uint64_t> nid(r.nnodes + 2, 0);
  nid.at(1) = 1;
  std::unordered_map<std::string, uint64_t> seen;
  std::string key;
  std::vector<std::string> parts;
  std::vector<uint8_t> b;
  auto put = [&](uint64_t v) {
    b.clear();
    put_varint(b, v);
    key.append(b.begin(), b.end());
  };
  auto lits = [&](const std::vector<int32_t>& ls) {
    put(ls.size());
    for (int32_t l : ls) put(lit_code(l));
  };
  Node out;
  r.scan([&](uint64_t id, const Node& n) {
    if (!reach[id]) return;
    key.assign(1, (char)n.kind);
    if (n.kind == K_AND) {
      out.kids.clear();
      for (uint64_t c : n.kids) out.kids.push_back(nid[c]);
      std::vector<uint64_t> ks = out.kids;
      std::sort(ks.begin(), ks.end());
      put(ks.size());
      for (uint64_t c : ks) put(c);
      lits(n.units);
      put(n.free.size());
      for (uint32_t v : n.free) put(v);
    } else {
      out.arcs = n.arcs;
      parts.clear();
      for (auto& a : out.arcs) {
        a.child = nid[a.child];
        std::string k0;
        std::swap(key, k0);
        lits(a.lits);
        put(a.child);
        std::swap(key, k0);
        parts.push_back(std::move(k0));
      }
      std::sort(parts.begin(), parts.end());
      put(parts.size());
      for (const auto& p : parts) key += p;
    }
    auto [it, fresh] = seen.try_emplace(key, 0);
    if (fresh) {
      if (n.kind == K_AND) it->second = w.add_and(out.kids, n.units, n.free);
      else {
        w.begin_or();
        for (const auto& a : out.arcs) w.or_arc(a.lits, a.child);
        it->second = w.end_or();
      }
    }
    nid[id] = it->second;
  });
  if (old_to_new) *old_to_new = nid;
  st.out_nodes = w.num_nodes();
  w.finish(nid[r.root], r.props | P_CLEAN, r.footer_kv, &r.defs);
  return st;
}

}  // namespace bdnnf
