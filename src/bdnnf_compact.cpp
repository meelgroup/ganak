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

// Usage:  bdnnf compact <in.bdnnf> <out.bdnnf>

#include <cstdlib>
#include <iostream>
#include <string>

#include "bdnnf_compact.hpp"

int bdnnf_compact(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: bdnnf compact <in.bdnnf> <out.bdnnf>\n"
                 "Drops the nodes the root does not reach and merges identical ones.\n";
    return EXIT_FAILURE;
  }
  if (std::string(argv[1]) == argv[2]) {
    std::cerr << "bdnnf compact: input and output must differ" << std::endl;
    return EXIT_FAILURE;
  }
  try {
    bdnnf::Reader r;
    r.open(argv[1]);
    bdnnf::Writer w;
    w.open(argv[2], r.nvars, r.header_kv);
    const auto st = bdnnf::compact(r, w);
    std::cerr << "bdnnf compact: " << st.in_nodes << " nodes, " << st.reachable << " reachable, "
              << st.out_nodes << " written" << std::endl;
  } catch (const bdnnf::Error& e) {
    std::cerr << "bdnnf compact: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
