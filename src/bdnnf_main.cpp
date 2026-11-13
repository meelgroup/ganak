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

// bdnnf: queries and conversions of BDNNF circuits (spec-bdnnf.md).

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

int bdnnf_verify(int argc, char** argv);
int bdnnf_query(int argc, char** argv);
int bdnnf_dump(int argc, char** argv);
int bdnnf_to_d4(int argc, char** argv);
int bdnnf_dot(int argc, char** argv);

namespace {
struct Cmd {
  const char* name;
  int (*fn)(int, char**);
  const char* help;
};
const Cmd cmds[] = {
  {"verify", bdnnf_verify, "check the format and every claimed property"},
  {"count", bdnnf_query, "weighted/projected model count, under --assume"},
  {"sample", bdnnf_query, "weighted/uniform random models, under --assume"},
  {"lits", bdnnf_query, "weighted/projected count of the models with each literal"},
  {"core", bdnnf_query, "literals true in every model (core; negated: dead)"},
  {"dump", bdnnf_dump, "print as text"},
  {"to-d4", bdnnf_to_d4, "convert to d4 .nnf"},
  {"dot", bdnnf_dot, "render as Graphviz DOT"},
};
}  // namespace

int main(int argc, char** argv) {
  // stderr: stdout carries the answers
  std::cerr << "c o bdnnf: queries on BDNNF d-DNNF circuits, from Ganak\n"
               "c o Its queries follow ddnnife by Sundermann, Raab, Hess, Thum and Schaefer,\n"
               "c o   \"Reusing d-DNNFs for Efficient Feature-Model Counting\", ACM TOSEM 2024,\n"
               "c o   https://github.com/practical-kc/ddnnife\n";
  if (argc >= 2)
    for (const auto& c : cmds)
      if (std::strcmp(argv[1], c.name) == 0) return c.fn(argc - 1, argv + 1);
  std::cerr << "Usage: bdnnf <command> [args]   (bdnnf <command> -h for its args)\n\n";
  for (const auto& c : cmds) std::cerr << "  " << c.name << std::string(9 - std::strlen(c.name), ' ')
                                       << c.help << "\n";
  return argc >= 2 && (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0)
             ? EXIT_SUCCESS : EXIT_FAILURE;
}
