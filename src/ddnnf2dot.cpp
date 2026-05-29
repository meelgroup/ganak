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

// ddnnf2dot: render a d4 .nnf circuit (raw from `ganak --compile`, or cleaned by
// `ddnnf-cleanup`) as a Graphviz DOT graph for visualization.
//
//   OR  node ('o') -> blue ellipse  labelled  "∨" + id
//   AND node ('a') -> yellow box     labelled  "∧" + id
//   TRUE ('t')      -> "⊤" ;  FALSE ('f') -> "⊥"
//   arc            -> edge labelled with its decided+propagated literals.
//
// Nodes are emitted in declaration order (root first, the d4 convention), arcs in
// file order. This is a pure pretty-printer: it does not require the circuit to be
// clean, and it changes nothing structurally.
//
// Usage:  ddnnf2dot <in.nnf> [out.dot]      ("-" or omitted out => stdout)

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Arc {
  int child;
  std::vector<int> lits;
};

[[noreturn]] void die(const std::string& msg) {
  std::cerr << "ddnnf2dot: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <in.nnf> [out.dot]   (\"-\"/omitted out => stdout)"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string in_path = argv[1];
  const std::string out_path = (argc == 3) ? argv[2] : "-";

  std::ifstream in(in_path);
  if (!in.good()) die("cannot open input file: " + in_path);

  std::unordered_map<int, char> type;             // node id -> 'f'|'t'|'a'|'o'
  std::unordered_map<int, std::vector<Arc>> arcs; // node id -> outgoing arcs
  std::vector<int> order;                         // node ids in declaration order
  int root = -1;

  std::string line;
  while (std::getline(in, line)) {
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
      if (!ss) die("malformed node declaration: " + line);
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
      if (nums.size() < 2) die("malformed arc line: " + line);
      Arc a;
      a.child = nums[1];
      a.lits.assign(nums.begin() + 2, nums.end());
      arcs[nums[0]].push_back(std::move(a));
    }
  }
  if (root == -1) die("input has no nodes");

  std::ofstream fout;
  std::ostream* outp = &std::cout;
  if (out_path != "-") {
    fout.open(out_path);
    if (!fout.good()) die("cannot open output file: " + out_path);
    outp = &fout;
  }
  std::ostream& out = *outp;

  out << "digraph ddnnf {\n";
  out << "  rankdir=TB;\n";
  out << "  node [fontname=\"monospace\", fontsize=11];\n";
  out << "  edge [fontname=\"monospace\", fontsize=10];\n";

  // Node declarations (in file order).
  for (int id : order) {
    out << "  n" << id << " [";
    switch (type[id]) {
      case 'o':
        out << "shape=ellipse, style=filled, fillcolor=\"#cfe8ff\", label=\"∨\\n" << id << "\"";
        break;
      case 'a':
        out << "shape=box, style=filled, fillcolor=\"#fff3c4\", label=\"∧\\n" << id << "\"";
        break;
      case 't':
        out << "shape=plaintext, label=\"⊤\"";
        break;
      default: // 'f'
        out << "shape=plaintext, label=\"⊥\"";
        break;
    }
    out << "];\n";
  }

  // Arcs, labelled with their decided+propagated literals.
  for (int id : order) {
    for (const auto& a : arcs[id]) {
      out << "  n" << id << " -> n" << a.child;
      if (!a.lits.empty()) {
        out << " [label=\"";
        for (size_t k = 0; k < a.lits.size(); k++) out << (k ? " " : "") << a.lits[k];
        out << "\"]";
      }
      out << ";\n";
    }
  }

  out << "  // root = n" << root << "\n";
  out << "}\n";
  out.flush();

  return EXIT_SUCCESS;
}
