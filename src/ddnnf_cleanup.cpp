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

// ddnnf-cleanup: post-process a d4 .nnf circuit emitted by `ganak --compile`.
//
// The streaming compiler writes a valid-but-loose file: the root is declared
// first (the d4 convention), but the file keeps the node ids it used internally
// and may contain UNREACHABLE ("dead") nodes (orphaned when an AND
// short-circuited to FALSE). Those are harmless to a root-rooted count, but a
// strict d4 reader (and the non-relaxed verifier) wants a clean circuit.
//
// This tool reads such a file, keeps only the nodes reachable from the root,
// renumbers them root=1 and contiguous (in breadth-first order from the root --
// the same layout the in-memory compiler used to produce), and writes the
// classic two-section d4 file (all node declarations, then all arc lines). The
// structural model count is unchanged.
//
// Usage:  ddnnf-cleanup <in.nnf> [out.nnf]      ("-" or omitted out => stdout)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "ddnnf_io.hpp"

namespace {

using ddnnf_io::Arc;

[[noreturn]] void die(const std::string& msg) {
  std::cerr << "ddnnf-cleanup: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <in.nnf> [out.nnf]   (\"-\"/omitted out => stdout)"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string in_path = argv[1];
  const std::string out_path = (argc == 3) ? argv[2] : "-";

  std::ifstream in(in_path);
  if (!in.good()) die("cannot open input file: " + in_path);

  std::unordered_map<int, char> type;             // node id -> 'f'|'t'|'a'|'o'
  std::unordered_map<int, std::vector<Arc>> arcs; // node id -> outgoing arcs
  std::vector<int> decl_order;                    // declaration order (unused; we renumber by BFS)
  std::string err;
  const int root = ddnnf_io::parse_nnf(in, type, arcs, decl_order, err);
  if (root == -1) die(err);

  // Breadth-first renumber from the root: root => 1, then children in arc order.
  // Unvisited nodes (dead) never receive an id and are dropped.
  std::unordered_map<int, int> newid;
  std::vector<int> order;                           // new-id order -> old id
  std::queue<int> q;
  newid[root] = 1;
  order.push_back(root);
  q.push(root);
  while (!q.empty()) {
    int nid = q.front();
    q.pop();
    auto it = arcs.find(nid);
    if (it == arcs.end()) continue;
    for (const auto& a : it->second) {
      if (newid.find(a.child) == newid.end()) {
        newid[a.child] = (int)order.size() + 1;
        order.push_back(a.child);
        q.push(a.child);
      }
    }
  }

  std::ofstream fout;
  std::ostream* out = &std::cout;
  if (out_path != "-") {
    fout.open(out_path);
    if (!fout.good()) die("cannot open output file: " + out_path);
    out = &fout;
  }

  // Section 1: declarations in new-id order.
  for (int oldid : order) *out << type[oldid] << " " << newid[oldid] << " 0\n";
  // Section 2: arcs in new-id order, with ids remapped.
  for (int oldid : order) {
    for (const auto& a : arcs[oldid]) {
      *out << newid[oldid] << " " << newid[a.child];
      for (int l : a.lits) *out << " " << l;
      *out << " 0\n";
    }
  }
  out->flush();

  std::cerr << "ddnnf-cleanup: declared=" << type.size()
            << " reachable=" << order.size()
            << " dropped=" << (type.size() - order.size()) << std::endl;
  return EXIT_SUCCESS;
}
