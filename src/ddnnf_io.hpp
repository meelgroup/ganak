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

// Header-only reader for the d4 .nnf format from `ganak --compile`. Used by the
// standalone ddnnf-cleanup and ddnnf2dot tools (no ganak/arjun/GMP link deps).

#include <cctype>
#include <cstddef>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

namespace ddnnf_io {

struct Arc {
  int child;
  std::vector<int> lits; // DIMACS literals true on this arc
};

// Parse a d4 .nnf stream into `type` (id -> 'f'|'t'|'a'|'o', 0 for unused
// ids), `arcs` (outgoing per id), and `order` (declaration order). Both `type`
// and `arcs` are sized to max_id+1 and indexed by node id. Ganak emits dense
// 0..N-1 ids (and ddnnf-cleanup renumbers root=1 contiguous), so a vector is
// the right shape: avoids the per-node hash-map overhead that dominates memory
// at multi-million-node scale. Returns the root (first-declared node), or -1
// with `err` set on malformed input or no nodes.
inline int parse_nnf(std::istream& in,
                     std::vector<char>& type,
                     std::vector<std::vector<Arc>>& arcs,
                     std::vector<int>& order,
                     std::string& err) {
  // Single pass with lazy resize. Ganak emits dense 0..N-1 ids and
  // ddnnf-cleanup renumbers root=1 contiguous, so the resize cost is
  // amortized O(1) per node in practice (we never see a gap or out-of-order
  // id and so push_back-equivalent semantics dominate). Critical at
  // multi-million-node scale: collecting intermediate buffers first would
  // double peak memory.
  type.clear();
  arcs.clear();
  order.clear();
  auto ensure_size = [&](int id) {
    if ((int)type.size() <= id) {
      type.resize((std::size_t)id + 1, 0);
      arcs.resize((std::size_t)id + 1);
    }
  };
  int root = -1;
  std::string line;
  while (std::getline(in, line)) {
    size_t i = 0;
    while (i < line.size() && std::isspace((unsigned char)line[i])) i++;
    if (i >= line.size()) continue;                 // blank
    if (line[i] == 'c') continue;                   // comment

    std::istringstream ss(line);
    if (line[i] == 'f' || line[i] == 't' || line[i] == 'a' || line[i] == 'o') {
      char t;
      int id;
      ss >> t >> id;
      if (!ss) { err = "malformed node declaration: " + line; return -1; }
      ensure_size(id);
      type[id] = t;
      order.push_back(id);
      if (root == -1) root = id;                    // first declared node is the root
    } else {
      std::vector<int> nums;
      int x;
      while (ss >> x) nums.push_back(x);
      if (!nums.empty() && nums.back() == 0) nums.pop_back();
      if (nums.size() < 2) { err = "malformed arc line: " + line; return -1; }
      int parent = nums[0];
      Arc a;
      a.child = nums[1];
      a.lits.assign(nums.begin() + 2, nums.end());
      ensure_size(parent > a.child ? parent : a.child);
      arcs[parent].push_back(std::move(a));
    }
  }
  if (root == -1) { err = "input has no nodes"; return -1; }
  return root;
}

} // namespace ddnnf_io
