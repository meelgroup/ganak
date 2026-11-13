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

// ddnnf-to-d4: convert a BDNNF file to the d4 .nnf text format, for legacy
// consumers (ddnnife and co.). Root = 1, ids contiguous, reachable nodes only;
// free variables become OR(v, -v) nodes. Run ddnnf-cleanup first if the
// consumer needs strict decomposability.
//
// Usage:  ddnnf-to-d4 <in.bdnnf> [out.nnf]      ("-" or omitted out => stdout)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "ddnnf_io.hpp"

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <in.bdnnf> [out.nnf]" << std::endl;
    return EXIT_FAILURE;
  }
  ddnnf_io::Circuit c;
  std::string err;
  if (!ddnnf_io::load(argv[1], c, err)) {
    std::cerr << "ddnnf-to-d4: " << err << std::endl;
    return EXIT_FAILURE;
  }
  if (!c.defs.empty()) {
    std::cerr << "ddnnf-to-d4: " << argv[1] << " defines vars by AIGs (a definitions frame), "
                 "which the d4 format cannot hold" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string out_path = argc == 3 ? argv[2] : "-";
  std::ofstream fout;
  if (out_path != "-") {
    fout.open(out_path);
    if (!fout.good()) {
      std::cerr << "ddnnf-to-d4: cannot open " << out_path << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::ostream& out = out_path == "-" ? std::cout : fout;
  ddnnf_io::write_d4(out, c);
  out.flush();
  if (!out.good()) {
    std::cerr << "ddnnf-to-d4: write error" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
