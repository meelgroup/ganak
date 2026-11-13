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

// In-memory circuit for bdnnf to-d4 and dot, loaded from BDNNF
// (spec-bdnnf.md). The model is d4-shaped: node id -> 'f'|'t'|'a'|'o', arcs
// carry literals, ids 0/1 are FALSE/TRUE. An AND's units become one arc to
// TRUE, and each free var v becomes one shared synthetic node
// OR(v -> TRUE, -v -> TRUE) appended after the file's nodes.

#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bdnnf.hpp"

namespace ddnnf_io {

struct Arc {
  int child;
  std::vector<int> lits;
};

struct Circuit {
  std::vector<char> type;
  std::vector<std::vector<Arc>> arcs;
  int root = 0;
  uint32_t nvars = 0;
  uint64_t props = 0;
  int nfile = 0;  // ids [2, nfile + 2) are nodes stored in the file
  bdnnf::KV header_kv;
  bdnnf::KV footer_kv;
  bdnnf::Defs defs;  // passed through unchanged

};

inline bool load(const std::string& path, Circuit& c, std::string& err) {
  try {
    bdnnf::Reader r;
    r.open(path);
    if (r.nnodes + 2 > (uint64_t)INT32_MAX / 2) throw bdnnf::Error("circuit too large for this tool");
    c.nvars = r.nvars;
    c.nfile = (int)r.nnodes;
    const size_t n = (size_t)r.nnodes + 2;
    c.type.assign(n, 0);
    c.arcs.assign(n, {});
    c.type[0] = 'f';
    c.type[1] = 't';
    std::vector<int> free_or(r.nvars + 1, -1);
    auto free_node = [&](uint32_t v) {
      if (free_or[v] < 0) {
        free_or[v] = (int)c.type.size();
        c.type.push_back('o');
        c.arcs.push_back({Arc{1, {(int)v}}, Arc{1, {-(int)v}}});
      }
      return free_or[v];
    };
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      std::vector<Arc> out;
      if (nd.kind == bdnnf::K_AND) {
        c.type[id] = 'a';
        for (uint64_t k : nd.kids) out.push_back(Arc{(int)k, {}});
        if (!nd.units.empty()) out.push_back(Arc{1, {nd.units.begin(), nd.units.end()}});
        for (uint32_t v : nd.free) out.push_back(Arc{free_node(v), {}});
      } else {
        c.type[id] = 'o';
        for (const auto& a : nd.arcs) out.push_back(Arc{(int)a.child, {a.lits.begin(), a.lits.end()}});
      }
      c.arcs[id] = std::move(out);
    });
    c.root = (int)r.root;
    c.props = r.props;
    c.header_kv = r.header_kv;
    c.footer_kv = r.footer_kv;
    c.defs = r.defs;
  } catch (const bdnnf::Error& e) {
    err = path + ": " + e.what();
    return false;
  }
  return true;
}

// d4 .nnf for legacy consumers: reachable nodes only, root = 1, ids
// contiguous in BFS order, declarations before arcs.
inline void write_d4(std::ostream& out, const Circuit& c) {
  std::vector<int> newid(c.type.size(), 0);
  std::vector<int> order;
  std::queue<int> q;
  newid[c.root] = 1;
  order.push_back(c.root);
  q.push(c.root);
  while (!q.empty()) {
    const int n = q.front();
    q.pop();
    for (const auto& a : c.arcs[n]) {
      if (newid[a.child] == 0) {
        newid[a.child] = (int)order.size() + 1;
        order.push_back(a.child);
        q.push(a.child);
      }
    }
  }
  for (int old : order) out << c.type[old] << " " << newid[old] << " 0\n";
  for (int old : order)
    for (const auto& a : c.arcs[old]) {
      out << newid[old] << " " << newid[a.child];
      for (int l : a.lits) out << " " << l;
      out << " 0\n";
    }
}

// Iterative post-order (children first) over the part reachable from `root`;
// fn(n) is called once per node. Recursion would blow the stack on big circuits.
template <class Fn>
void post_order(const Circuit& c, std::vector<char>& done, Fn&& fn) {
  std::vector<std::pair<int, int>> stack;
  if (done[c.root]) return;
  stack.emplace_back(c.root, 0);
  while (!stack.empty()) {
    auto& top = stack.back();
    const int n = top.first;
    if (done[n]) { stack.pop_back(); continue; }
    const int idx = top.second;
    if (idx >= (int)c.arcs[n].size()) {
      fn(n);
      done[n] = 1;
      stack.pop_back();
      continue;
    }
    const int ch = c.arcs[n][idx].child;
    top.second = idx + 1;
    if (!done[ch]) stack.emplace_back(ch, 0);
  }
}

}  // namespace ddnnf_io
