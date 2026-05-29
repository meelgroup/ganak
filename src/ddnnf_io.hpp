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

// Shared, header-only reader for the d4 .nnf circuit format emitted by
// `ganak --compile`. Used by the standalone ddnnf-cleanup and ddnnf2dot tools;
// keeps them self-contained (no ganak/arjun/GMP link deps).

#include <cctype>
#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace ddnnf_io {

struct Arc {
  int child;
  std::vector<int> lits; // DIMACS literals true on this arc
};

// Parse a d4 .nnf stream. Fills `type` (id -> 'f'|'t'|'a'|'o'), `arcs` (id ->
// outgoing arcs) and `order` (node ids in declaration order). Returns the root
// id (the first-declared node, the d4 convention), or -1 with `err` set on
// malformed input or a file that declares no nodes.
inline int parse_nnf(std::istream& in,
                     std::unordered_map<int, char>& type,
                     std::unordered_map<int, std::vector<Arc>>& arcs,
                     std::vector<int>& order,
                     std::string& err) {
  int root = -1;
  std::string line;
  while (std::getline(in, line)) {
    // First non-space token decides the line kind.
    size_t i = 0;
    while (i < line.size() && std::isspace((unsigned char)line[i])) i++;
    if (i >= line.size()) continue;                 // blank
    if (line[i] == 'c') continue;                   // comment

    std::istringstream ss(line);
    if (line[i] == 'f' || line[i] == 't' || line[i] == 'a' || line[i] == 'o') {
      // Node declaration:  <type> <id> 0
      char t;
      int id;
      ss >> t >> id;
      if (!ss) { err = "malformed node declaration: " + line; return -1; }
      type[id] = t;
      arcs.emplace(id, std::vector<Arc>{});
      order.push_back(id);
      if (root == -1) root = id;                    // first declared node is the root
    } else {
      // Arc line:  <parent> <child> [lits...] 0
      std::vector<int> nums;
      int x;
      while (ss >> x) nums.push_back(x);
      if (!nums.empty() && nums.back() == 0) nums.pop_back();  // drop terminator
      if (nums.size() < 2) { err = "malformed arc line: " + line; return -1; }
      Arc a;
      a.child = nums[1];
      a.lits.assign(nums.begin() + 2, nums.end());
      arcs[nums[0]].push_back(std::move(a));
    }
  }
  if (root == -1) err = "input has no nodes";
  return root;
}

} // namespace ddnnf_io
