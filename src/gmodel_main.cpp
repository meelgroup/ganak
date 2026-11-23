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

// ganak-model: probabilistic models (spec-ganak-model.md) to weighted projected
// CNF; --solve answers the queries with Ganak.

#include <cryptominisat5/solvertypesmini.h>
#include <arjun/arjun.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ganak.hpp"
#include "gmodel.hpp"

using std::cerr;
using std::cout;
using std::string;
using std::vector;

namespace {

const char* usage =
  "usage: ganak-model [options] <model.gm | ->\n"
  "  (no option)      write the weighted projected CNF to stdout\n"
  "  -o FILE          write the CNF to FILE\n"
  "  --json FILE      write the model, its statements and errors as JSON (- = stdout)\n"
  "  --dot FILE       write the model graph as DOT (- = stdout)\n"
  "  --solve          answer every query and requirement with Ganak\n"
  "  --check          only report errors\n"
  "  --arjun 0|1      --solve: preprocess with Arjun (default 1)\n"
  "Exit code: 1 on errors, 3 if a requirement fails (--solve).\n";

int do_arjun = 1;

CMSat::Lit mklit(int l) { return CMSat::Lit((uint32_t)std::abs(l) - 1, l < 0); }

mpq_class count(const gm::Cnf& c, const vector<int>& assume) {
  std::unique_ptr<CMSat::FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
  ArjunNS::SimplifiedCNF cnf(fg);
  cnf.new_vars(c.nvars);
  vector<CMSat::Lit> cl;
  for (const auto& x : c.clauses) {
    cl.clear();
    for (const int l : x) cl.push_back(mklit(l));
    cnf.add_clause(cl);
  }
  for (const int a : assume) cnf.add_clause({mklit(a)});
  cnf.set_weighted(true);
  cnf.set_projected(true);
  for (const auto& [v, w] : c.w) {
    cnf.set_lit_weight(mklit(v), ArjunNS::FMpq(w.first));
    cnf.set_lit_weight(mklit(-v), ArjunNS::FMpq(w.second));
  }
  vector<uint32_t> samp;
  for (const int v : c.random) samp.push_back((uint32_t)v - 1);
  cnf.set_sampl_vars(samp);

  GanakInt::CounterConfiguration conf;
  conf.verb = 0;
  ArjunNS::Arjun arjun;
  arjun.set_verb(0);
  ArjunNS::Arjun::InterpConf iconf;
  if (do_arjun) {
    arjun.standalone_minimize_indep(cnf, iconf, false);
    if (cnf.get_sampl_vars().size() >= 10) {
      ArjunNS::SimpConf simp_conf;
      ArjunNS::Arjun::ElimToFileConf etof_conf;
      arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf, iconf);
    } else {
      conf.disconnected_allowed = true;
      cnf.renumber_sampling_vars_for_ganak();
    }
  } else {
    conf.disconnected_allowed = true;
    cnf.renumber_sampling_vars_for_ganak();
  }
  cnf.remove_equiv_weights();
  std::unique_ptr<CMSat::Field> cnt = cnf.get_multiplier_weight()->dup();
  if (!cnt->is_zero()) {
    Ganak counter(conf, fg);
    setup_ganak(cnf, counter);
    *cnt *= *counter.count();
  }
  const auto* q = dynamic_cast<const ArjunNS::FMpq*>(cnt.get());
  if (!q) {
    cerr << "ganak-model: internal error: the count is not a rational\n";
    exit(EXIT_FAILURE);
  }
  return q->val;
}

string pretty(const gm::Stmt& s) {
  string t = s.text;
  const string g = " given ";
  const auto p = t.find(g);
  if (p != string::npos) t = t.substr(0, p) + " | " + t.substr(p + g.size());
  return "P(" + t + ")";
}

int solve(const gm::Model& m, const gm::Cnf& c) {
  std::map<vector<int>, mpq_class> memo;
  auto cnt = [&](vector<int> a) {
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    auto it = memo.find(a);
    if (it != memo.end()) return it->second;
    const mpq_class r = count(c, a);
    memo[a] = r;
    return r;
  };
  const mpq_class z = cnt({});
  bool has_constraints = false;
  for (const auto& s : m.stmts) has_constraints |= s.kind == gm::SKind::Constraint;
  if (has_constraints) cout << "P(constraints) = " << gm::fmt_q(z) << " ~ " << gm::approx(z) << "\n";
  if (z == 0) {
    cout << "The constraints can never hold: every probability is undefined.\n";
    return EXIT_FAILURE;
  }
  int ret = EXIT_SUCCESS;
  for (const auto& s : m.stmts) {
    if (s.kind == gm::SKind::Constraint) continue;
    const mpq_class den = cnt(s.given_lits);
    if (den == 0) {
      cout << pretty(s) << " is undefined: the condition has probability 0\n";
      if (s.kind == gm::SKind::Require) ret = 3;
      continue;
    }
    vector<int> a = s.given_lits;
    a.push_back(s.lit);
    mpq_class p = cnt(a) / den;
    p.canonicalize();
    if (s.kind == gm::SKind::Query) {
      cout << pretty(s) << " = " << gm::fmt_q(p) << " ~ " << gm::approx(p) << "\n";
      continue;
    }
    const bool ok = s.op == "<" ? p < s.bound : s.op == "<=" ? p <= s.bound : s.op == ">" ? p > s.bound : p >= s.bound;
    cout << "require " << pretty(s) << " " << s.op << " " << gm::approx(s.bound, 3) << ": " << (ok ? "PASS" : "FAIL")
         << " (" << gm::fmt_q(p) << " ~ " << gm::approx(p) << ")\n";
    if (!ok) ret = 3;
  }
  return ret;
}

bool write_to(const string& path, const std::function<void(std::ostream&)>& f) {
  if (path == "-") {
    f(cout);
    return true;
  }
  std::ofstream o(path);
  if (!o) {
    cerr << "ganak-model: cannot write " << path << "\n";
    return false;
  }
  f(o);
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  string in, out, json, dot;
  bool do_solve = false, check = false;
  for (int i = 1; i < argc; i++) {
    const string a = argv[i];
    auto val = [&]() -> string {
      if (i + 1 >= argc) {
        cerr << "ganak-model: " << a << " needs a value\n" << usage;
        exit(EXIT_FAILURE);
      }
      return argv[++i];
    };
    if (a == "-h" || a == "--help") {
      cout << usage;
      return EXIT_SUCCESS;
    } else if (a == "-o") out = val();
    else if (a == "--json") json = val();
    else if (a == "--dot") dot = val();
    else if (a == "--solve") do_solve = true;
    else if (a == "--check") check = true;
    else if (a == "--arjun") do_arjun = std::stoi(val());
    else if (a.size() > 1 && a[0] == '-') {
      cerr << "ganak-model: unknown option " << a << "\n" << usage;
      return EXIT_FAILURE;
    } else if (in.empty()) in = a;
    else {
      cerr << "ganak-model: more than one input\n" << usage;
      return EXIT_FAILURE;
    }
  }
  if (in.empty()) {
    cerr << usage;
    return EXIT_FAILURE;
  }
  std::stringstream buf;
  if (in == "-") buf << std::cin.rdbuf();
  else {
    std::ifstream f(in);
    if (!f) {
      cerr << "ganak-model: cannot open " << in << "\n";
      return EXIT_FAILURE;
    }
    buf << f.rdbuf();
  }
  const string src = buf.str();
  const string fname = in == "-" ? "<stdin>" : in;

  gm::Model m = gm::parse(src);
  gm::Cnf c;
  const bool parsed = m.ok();
  if (parsed) c = gm::encode(m);
  if (!json.empty() && !write_to(json, [&](std::ostream& o) { gm::write_json(o, m, parsed ? &c : nullptr); }))
    return EXIT_FAILURE;
  for (const auto& e : m.errors) cerr << fname << ":" << e.loc.line << ":" << e.loc.col << ": error: " << e.msg << "\n";
  if (!m.ok()) return EXIT_FAILURE;
  if (check) return EXIT_SUCCESS;
  if (!dot.empty() && !write_to(dot, [&](std::ostream& o) { gm::write_dot(o, m); })) return EXIT_FAILURE;
  if (!out.empty() || (!do_solve && json != "-" && dot != "-"))
    if (!write_to(out.empty() ? "-" : out, [&](std::ostream& o) { gm::write_cnf(o, m, c); })) return EXIT_FAILURE;
  return do_solve ? solve(m, c) : EXIT_SUCCESS;
}
