/******************************************
Copyright (C) 2026 Authors of GANAK, see AUTHORS file

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

// The ganak-model language (spec-ganak-model.md): parser, CNF encoder, writers.

#pragma once

#include <gmpxx.h>

#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace gm {

struct Loc {
  int line = 0, col = 0;
};

struct Diag {
  Loc loc;
  std::string msg;
};

enum class Op { Const, Ref, Not, And, Or, Xor, Iff, Implies, Ite, AtLeast, AtMost, Exactly, Flip };

struct Expr {
  Op op = Op::Const;
  Loc loc;
  bool val = false;
  std::string name;
  int node = -1;  // Ref: the node; Flip: its anonymous event
  int k = 0;
  mpq_class p;
  std::vector<int> kids;
};

// Event, ChoiceValue and Flip nodes are the random variables.
enum class Kind { Event, ChoiceValue, Flip, Def };

struct Node {
  std::string name, label, text;
  Kind kind = Kind::Event;
  Loc loc;
  mpq_class p;
  int expr = -1;
  std::vector<int> rules;
  int self = -1;
  int choice = -1;
  bool anon = false;
  int var = 0;
};

struct Rule {
  int body = -1, head = -1, flip = -1;
  mpq_class p;
  Loc loc;
};

struct Choice {
  std::string name, label;
  Loc loc;
  std::vector<int> values;
};

enum class SKind { Constraint, Query, Require };

struct Stmt {
  SKind kind = SKind::Query;
  Loc loc;
  std::string text;
  int expr = -1;
  std::vector<int> given;
  std::string op;
  mpq_class bound;
  int lit = 0;
  std::vector<int> given_lits;
};

struct Model {
  std::vector<Node> nodes;
  std::map<std::string, int> by_name;
  std::vector<Expr> exprs;
  std::vector<Rule> rules;
  std::vector<Choice> choices;
  std::map<std::string, int> choice_by_name;
  std::vector<Stmt> stmts;
  std::vector<Diag> errors;
  bool ok() const { return errors.empty(); }
  bool is_random(int n) const { return nodes[n].kind != Kind::Def; }
};

struct Cnf {
  int nvars = 0;
  std::vector<std::vector<int>> clauses;
  std::map<int, std::pair<mpq_class, mpq_class>> w;  // var -> (w(v), w(-v))
  std::vector<int> show;
  std::vector<int> random;
};

Model parse(const std::string& src);
Cnf encode(Model& m);

bool parse_prob(const std::string& s, mpq_class& out);
std::string fmt_q(const mpq_class& q);
std::string approx(const mpq_class& q, int digits = 6);
std::string json_str(const std::string& s);

void write_cnf(std::ostream& o, const Model& m, const Cnf& c);
void write_json(std::ostream& o, const Model& m, const Cnf* c);
void write_dot(std::ostream& o, const Model& m);

// Values of every node, given the values of the random nodes (indexed by node).
std::vector<char> eval(const Model& m, const std::vector<char>& random_val);
bool eval_expr(const Model& m, int e, const std::vector<char>& node_val);

}  // namespace gm
