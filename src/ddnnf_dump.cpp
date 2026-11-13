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

// ddnnf-dump: print a BDNNF file as text (spec-bdnnf.md, section 7), streaming.
//
// Usage:  ddnnf-dump <in.bdnnf | -> [out.txt]      ("-" or omitted out => stdout)

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "bdnnf.hpp"

namespace {

std::string ref(uint64_t r) {
  if (r == 0) return "F";
  if (r == 1) return "T";
  return std::to_string(r);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <in.bdnnf | -> [out.txt]" << std::endl;
    return EXIT_FAILURE;
  }
  std::ios::sync_with_stdio(false);
  std::ofstream fout;
  if (argc == 3 && std::string(argv[2]) != "-") {
    fout.open(argv[2]);
    if (!fout.good()) {
      std::cerr << "ddnnf-dump: cannot open " << argv[2] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::ostream& out = fout.is_open() ? fout : std::cout;
  try {
    bdnnf::Reader r;
    r.open(argv[1]);
    out << "bdnnf " << bdnnf::VERSION << " nvars " << r.nvars << "\n";
    for (const auto& [k, v] : r.header_kv) out << "meta " << k << " " << v << "\n";
    r.scan([&](uint64_t id, const bdnnf::Node& n) {
      out << id;
      if (n.kind == bdnnf::K_AND) {
        out << " and";
        for (uint64_t k : n.kids) out << " " << ref(k);
        if (!n.units.empty()) {
          out << " units";
          for (int l : n.units) out << " " << l;
        }
        if (!n.free.empty()) {
          out << " free";
          for (uint32_t v : n.free) out << " " << v;
        }
      } else {
        out << " or";
        for (size_t i = 0; i < n.arcs.size(); i++) {
          out << (i ? " ;" : "");
          for (int l : n.arcs[i].lits) out << " " << l;
          out << " -> " << ref(n.arcs[i].child);
        }
      }
      out << "\n";
    });
    for (size_t j = 0; j < r.defs.nodes.size(); j++)
      out << "defnode " << r.nvars + 1 + j << " " << r.defs.nodes[j].first << " " << r.defs.nodes[j].second << "\n";
    for (const auto& [v, rf] : r.defs.defs) out << "def " << v << " " << rf << "\n";
    out << "root " << ref(r.root) << "\n";
    out << "props";
    const char* names[] = {"DECOMPOSABLE", "DETERMINISTIC", "SMOOTH", "CLEAN"};
    for (int b = 0; b < 4; b++) if (r.props & (1ULL << b)) out << " " << names[b];
    out << "\n";
    for (const auto& [k, v] : r.footer_kv) out << "meta " << k << " " << v << "\n";
  } catch (const bdnnf::Error& e) {
    out.flush();
    std::cerr << "ddnnf-dump: " << argv[1] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  out.flush();
  return EXIT_SUCCESS;
}
