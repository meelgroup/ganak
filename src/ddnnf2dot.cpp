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

// ddnnf2dot: render a d4 .nnf circuit (raw from `ganak --compile` or cleaned by
// `ddnnf-cleanup`) as a Graphviz DOT graph.
//
//   OR ('o') -> blue ellipse "∨";  AND ('a') -> yellow box "∧";  FALSE ('f') -> "⊥"
//   No operator arc carries literals. An OR branch that fixes literal(s) and
//     continues to a child operator routes through a fresh AND node: OR -> AND,
//     AND -> child (gray "{vars}"), AND -> one green "{v}" literal-leaf per fixed
//     literal -- so there are no literal-carrying OR->operator edges. AND arcs
//     likewise pull each fixed literal out as its own green literal-leaf child
//     (gray "{v}" edge), so every operator stays a visible conjunction/
//     disjunction of its drawn operands.
//
// The shared TRUE ('t') sink is not drawn; instead each literal on an arc into ⊤
// becomes its own green leaf ("3", "¬1"). A single-literal terminal OR arc lands
// directly on that leaf; a multi-literal one fans out through a fresh AND with a
// "{v}" leaf per literal. A ⊤ node is kept only for the degenerate cases (a
// literal-free terminal arc, or a whole-circuit tautology). Operator/FALSE nodes
// carry the .nnf node id as a small gray "#<id>" under the symbol; variables
// still appear only on edges/leaves.
//
// Pure pretty-printer (no structural change); does not require a clean circuit.
//
// Usage:  ddnnf2dot <in.nnf> [out.dot]      ("-" or omitted out => stdout)

#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "ddnnf_io.hpp"

namespace {

using ddnnf_io::Arc;

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

  std::vector<char> type;                          // node id -> 'f'|'t'|'a'|'o' (0 = no node)
  std::vector<std::vector<Arc>> arcs;              // node id -> outgoing arcs
  std::vector<int> order;                          // node ids in declaration order
  std::string err;
  const int root = ddnnf_io::parse_nnf(in, type, arcs, order, err);
  if (root == -1) die(err);

  std::ofstream fout;
  std::ostream* outp = &std::cout;
  if (out_path != "-") {
    fout.open(out_path);
    if (!fout.good()) die("cannot open output file: " + out_path);
    outp = &fout;
  }
  std::ostream& out = *outp;

  auto is_true = [&](int id) {
    return id >= 0 && id < (int)type.size() && type[id] == 't';
  };

  // Variables in a node's reachable subtree (on arc lits); labels AND children.
  // Memoized DFS; std::map keeps cached references stable across insertions.
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

  // Per-literal TRUE sinks: each literal on an arc into ⊤ gets its own shared leaf;
  // a multi-literal arc fans out through a fresh AND node. Build edges first
  // (discovering the needed sinks/fan-out nodes), then declare everything.
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

  // Bold-blue decided/implied literals: "1[2,-3,4]", or just "1". HTML fragment.
  auto blue_frag = [](std::ostream& os, const std::vector<int>& lits) {
    os << "<FONT COLOR=\"#1565c0\" POINT-SIZE=\"9\"><B>" << lits[0];
    if (lits.size() > 1) {
      os << "[";
      for (size_t k = 1; k < lits.size(); k++) os << (k > 1 ? "," : "") << lits[k];
      os << "]";
    }
    os << "</B></FONT>";
  };
  // Gray variable-set "{2,3,4}" of an arc: its lit vars plus the child subtree's
  // vars (none for a ⊤ child). HTML fragment. Empty set -> nothing written.
  auto gray_frag = [&](std::ostream& os, const std::vector<int>& lits, int child) {
    std::set<int> vs;
    for (int l : lits) vs.insert(l > 0 ? l : -l);
    if (!is_true(child)) for (int v : vars_of(child)) vs.insert(v);
    if (vs.empty()) return false;
    os << "<FONT COLOR=\"#777777\">{";
    bool first = true;
    for (int v : vs) { os << (first ? "" : ",") << v; first = false; }
    os << "}</FONT>";
    return true;
  };

  // Every AND out-edge is a decomposition arc, so it ALWAYS carries the child
  // component's variable set "{...}" (gray). AND arcs never show blue literals:
  // any literal fixed on an AND arc is split off as its own literal-leaf child
  // (see lit_leaf), so the structural edge only needs the child's var-set.
  // Non-AND (OR) edges keep the plain blue decided-literal label, and terminal
  // arcs into ⊤ keep landing on per-literal leaves.
  auto and_label = [&](std::ostream& os, const std::vector<int>& lits, int child) {
    std::ostringstream body;
    const bool has_blue = !lits.empty();
    if (has_blue) blue_frag(body, lits);
    std::ostringstream g;
    const bool has_gray = gray_frag(g, lits, child);
    if (has_blue && has_gray) body << "<BR/>";
    if (has_gray) body << g.str();
    const std::string b = body.str();
    if (b.empty()) return;
    os << " [label=<" << b << ">";
    if (!has_blue) os << ", fontsize=9";   // gray-only labels match decomp style
    os << "]";
  };
  // A literal fixed on an AND arc, rendered as its own green literal-leaf child
  // with a gray "{v}" decomposition-style edge (so every AND out-edge is a
  // labelled conjunct, and the AND is a visible conjunction of its operands).
  auto lit_leaf_from = [&](std::ostream& os, const std::string& from, int l) {
    os << "  " << from << " -> " << sink_for(l)
       << " [label=<<FONT COLOR=\"#777777\">{" << std::abs(l)
       << "}</FONT>>, fontsize=9];\n";
  };
  auto lit_leaf = [&](std::ostream& os, int from_id, int l) {
    lit_leaf_from(os, "n" + std::to_string(from_id), l);
  };

  for (int id : order) {
    if (is_true(id)) continue;            // ⊤ replaced by per-literal sinks
    const bool parent_and = (type[id] == 'a');
    for (const auto& a : arcs[id]) {
      if (is_true(a.child)) {
        // Terminal arc -> per-literal sink(s). For an AND parent each fixed
        // literal becomes its own gray-"{v}" leaf child (no fan-out node, no
        // blue). For an OR parent the leaf names the literal so a single-literal
        // edge stays bare; a multi-literal edge fans out through a fresh AND.
        if (a.lits.empty()) {
          need_true = true;
          edges << "  n" << id << " -> T;\n";
        } else if (parent_and) {
          for (int l : a.lits) lit_leaf(edges, id, l);
        } else if (a.lits.size() == 1) {
          edges << "  n" << id << " -> " << sink_for(a.lits[0]) << ";\n";
        } else {
          // OR branch fixing several literals (all terminal): fan out through a
          // fresh AND with one gray-"{v}" literal-leaf per literal. The decided
          // literal (lits[0]) labels the OR->AND edge in bold blue so the OR's
          // branch variable stays visible.
          std::string an = "and" + std::to_string(fanout_nodes.size());
          fanout_nodes.push_back(an);
          const std::vector<int> dec{a.lits[0]};
          edges << "  n" << id << " -> " << an << " [label=<";
          blue_frag(edges, dec);
          edges << ">];\n";
          for (int l : a.lits) lit_leaf_from(edges, an, l);
        }
        continue;
      }
      // Non-terminal arc.
      if (parent_and) {
        // AND decomposition: the structural edge shows just the child
        // component's variables; any decided literals fixed on this arc are
        // pulled out as their own literal-leaf children below, so the AND
        // visibly conjoins them. (A single-child arc carrying the literal only
        // in its label would look like the AND "ANDs nothing".)
        edges << "  n" << id << " -> n" << a.child;
        static const std::vector<int> none;
        and_label(edges, none, a.child);
        edges << ";\n";
        for (int l : a.lits) lit_leaf(edges, id, l);
        continue;
      }
      // OR parent.
      if (a.lits.empty()) {
        // Bare structural edge to the child operator.
        edges << "  n" << id << " -> n" << a.child << ";\n";
      } else {
        // OR branch that fixes literal(s): the OR->child edge becomes
        // OR -> AND -> child, with the *decided* literal (lits[0]) kept as a
        // bold-blue label on the OR->AND edge so the OR's branch variable stays
        // visible. The AND then conjoins the child (gray "{vars}") with one
        // green "{v}" literal-leaf per fixed literal -- no literal-carrying
        // OR->operator edge remains.
        std::string an = "and" + std::to_string(fanout_nodes.size());
        fanout_nodes.push_back(an);
        const std::vector<int> dec{a.lits[0]};
        edges << "  n" << id << " -> " << an << " [label=<";
        blue_frag(edges, dec);
        edges << ">];\n";
        edges << "  " << an << " -> n" << a.child;
        static const std::vector<int> none;
        and_label(edges, none, a.child);
        edges << ";\n";
        for (int l : a.lits) lit_leaf_from(edges, an, l);
      }
    }
  }

  out << "digraph ddnnf {\n";
  out << "  rankdir=TB;\n";
  out << "  node [fontname=\"monospace\", fontsize=11];\n";
  out << "  edge [fontname=\"monospace\", fontsize=10];\n";

  // Operator / FALSE declarations (file order; ⊤ dropped). Labels are just the
  // operator symbol -- variables appear only on edges/sinks.
  for (int id : order) {
    if (is_true(id)) continue;
    // Two-line label: the operator symbol on top, the .nnf node id below it
    // (a smaller gray "#<id>") so nodes can be referenced back to the circuit.
    const std::string idtag = "<BR/><FONT COLOR=\"#777777\" POINT-SIZE=\"8\">#"
                              + std::to_string(id) + "</FONT>";
    out << "  n" << id << " [";
    if (type[id] == 'o')      out << "shape=ellipse, style=filled, fillcolor=\"#cfe8ff\", label=<∨" << idtag << ">";
    else if (type[id] == 'a') out << "shape=box, style=filled, fillcolor=\"#fff3c4\", label=<∧" << idtag << ">";
    else                      out << "shape=plaintext, label=<⊥" << idtag << ">";  // 'f'
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
