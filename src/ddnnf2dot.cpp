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
//   OR  node ('o')  -> blue ellipse "∨"
//   AND node ('a')  -> yellow box   "∧"
//   FALSE ('f')      -> "⊥"
//   arc            -> edge. An arc carrying literals is labelled in bold blue as
//                     "1 [2,-3]" -- the decided literal, then the propagated ones
//                     in brackets (just "1" when nothing propagated). An AND node's
//                     decomposition arcs carry no literals, so we instead label
//                     each with the set of variables in that child's component --
//                     this shows how the AND cut the problem into disjoint parts.
//
// The single TRUE ('t') sink is NOT drawn as one shared ⊤ node. Instead, every
// literal labelling an arc into ⊤ becomes its own (shared) green leaf ("3", "¬1",
// ...). A terminal arc with several literals fans out through a fresh AND node to
// each literal's leaf. This declutters the otherwise-hub-like ⊤. (A ⊤ node is only
// kept for the degenerate cases: a literal-free terminal arc, or a tautology whose
// whole circuit is ⊤.) Node labels carry no ids -- variables appear only on edges
// and leaves.
//
// Nodes are emitted in declaration order (root first, the d4 convention), arcs in
// file order. This is a pure pretty-printer: it does not require the circuit to be
// clean, and it changes nothing structurally.
//
// Usage:  ddnnf2dot <in.nnf> [out.dot]      ("-" or omitted out => stdout)

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
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

  auto is_true = [&](int id) {
    auto it = type.find(id);
    return it != type.end() && it->second == 't';
  };

  // Variables appearing in a node's reachable subtree (on arc literals). Used to
  // describe how an AND node splits its component across children. Memoized DFS
  // over the DAG; std::map keeps cached references stable across insertions.
  std::map<int, std::vector<int>> sub_vars;
  std::function<const std::vector<int>&(int)> vars_of = [&](int nid) -> const std::vector<int>& {
    auto it = sub_vars.find(nid);
    if (it != sub_vars.end()) return it->second;
    std::set<int> acc;
    for (const auto& a : arcs[nid]) {
      for (int l : a.lits) acc.insert(std::abs(l));
      for (int v : vars_of(a.child)) acc.insert(v);
    }
    auto& vec = sub_vars[nid];
    vec.assign(acc.begin(), acc.end());
    return vec;
  };

  // Per-literal TRUE sinks. Instead of one shared ⊤ with many incoming arcs, every
  // literal that labels an arc into ⊤ gets its own (shared) leaf node. An arc with
  // several literals (a decision + propagated lits) fans out through a fresh AND
  // node to each literal's sink. We build the edges first (discovering the needed
  // sinks/fan-out nodes), then declare everything.
  std::map<int, std::string> sink_name;   // literal -> dot node name
  std::vector<int> sink_order;            // literals in first-seen order
  auto sink_for = [&](int l) -> const std::string& {
    auto it = sink_name.find(l);
    if (it != sink_name.end()) return it->second;
    std::string nm = std::string("lit") + (l > 0 ? "p" : "n") + std::to_string(std::abs(l));
    sink_order.push_back(l);
    return sink_name.emplace(l, std::move(nm)).first->second;
  };
  std::vector<std::string> fanout_nodes;  // fresh AND nodes for multi-literal arcs
  bool need_true = is_true(root);         // a lone ⊤ root still needs one node
  std::ostringstream edges;

  // The bold-blue edge label. lits[0] is the DECIDED literal; the rest were unit-
  // propagated. Rendered as "1 [2,-3,4]" (decided, then propagated in brackets), or
  // just "1" when nothing propagated.
  auto blue_lits = [](std::ostream& os, const std::vector<int>& lits) {
    os << " [label=<<FONT COLOR=\"#1565c0\"><B>" << lits[0];
    if (lits.size() > 1) {
      os << " [";
      for (size_t k = 1; k < lits.size(); k++) os << (k > 1 ? "," : "") << lits[k];
      os << "]";
    }
    os << "</B></FONT>>]";
  };

  for (int id : order) {
    if (is_true(id)) continue;            // ⊤ nodes are replaced by per-literal sinks
    const bool parent_and = (type[id] == 'a');
    for (const auto& a : arcs[id]) {
      if (is_true(a.child)) {
        // Terminal arc -> per-literal sink(s). Edges that land DIRECTLY on a green
        // leaf carry no label (the leaf already names the literal). A multi-literal
        // arc goes via a fresh AND node, and that n->AND edge is not a leaf edge, so
        // it keeps its bold-blue literals (and the AND->leaf edges stay bare).
        if (a.lits.empty()) {
          need_true = true;
          edges << "  n" << id << " -> T;\n";
        } else if (a.lits.size() == 1) {
          edges << "  n" << id << " -> " << sink_for(a.lits[0]) << ";\n";
        } else {
          std::string an = "and" + std::to_string(fanout_nodes.size());
          fanout_nodes.push_back(an);
          edges << "  n" << id << " -> " << an;
          blue_lits(edges, a.lits);
          edges << ";\n";
          for (int l : a.lits) edges << "  " << an << " -> " << sink_for(l) << ";\n";
        }
        continue;
      }
      // Non-terminal arc.
      edges << "  n" << id << " -> n" << a.child;
      if (!a.lits.empty()) {
        // Decision / implied literals, bold blue.
        blue_lits(edges, a.lits);
      } else if (parent_and) {
        // AND decomposition: show this child component's variables.
        const auto& vs = vars_of(a.child);
        if (!vs.empty()) {
          edges << " [label=<<FONT COLOR=\"#777777\">{";
          for (size_t k = 0; k < vs.size(); k++) edges << (k ? "," : "") << vs[k];
          edges << "}</FONT>>, fontsize=9]";
        }
      }
      edges << ";\n";
    }
  }

  out << "digraph ddnnf {\n";
  out << "  rankdir=TB;\n";
  out << "  node [fontname=\"monospace\", fontsize=11];\n";
  out << "  edge [fontname=\"monospace\", fontsize=10];\n";

  // Operator / FALSE node declarations (in file order; ⊤ nodes dropped). Labels
  // are just the operator symbol -- variables only ever appear on edges/sinks.
  for (int id : order) {
    if (is_true(id)) continue;
    out << "  n" << id << " [";
    if (type[id] == 'o')      out << "shape=ellipse, style=filled, fillcolor=\"#cfe8ff\", label=\"∨\"";
    else if (type[id] == 'a') out << "shape=box, style=filled, fillcolor=\"#fff3c4\", label=\"∧\"";
    else                      out << "shape=plaintext, label=\"⊥\"";  // 'f'
    out << "];\n";
  }
  // Fan-out AND nodes (for multi-literal terminal arcs).
  for (const auto& an : fanout_nodes)
    out << "  " << an << " [shape=box, style=filled, fillcolor=\"#fff3c4\", label=\"∧\"];\n";
  // Per-literal TRUE sinks (green leaves), e.g. "3" or "¬1".
  for (int l : sink_order)
    out << "  " << sink_name[l] << " [shape=box, style=\"rounded,filled\", fillcolor=\"#c8e6c9\", "
        << "label=\"" << (l > 0 ? "" : "¬") << std::abs(l) << "\"];\n";
  if (need_true) out << "  T [shape=plaintext, label=\"⊤\"];\n";

  out << edges.str();
  out << "  // root = " << (is_true(root) ? "T" : "n" + std::to_string(root)) << "\n";
  out << "}\n";
  out.flush();

  return EXIT_SUCCESS;
}
