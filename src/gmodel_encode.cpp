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

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <set>
#include <sstream>

#include "gmodel.hpp"

namespace gm {

bool parse_prob(const std::string& str, mpq_class& out) {
  std::string s = str;
  bool pct = false;
  if (!s.empty() && s.back() == '%') {
    pct = true;
    s.pop_back();
  }
  if (s.empty()) return false;
  auto digits = [](const std::string& x) {
    return !x.empty() && std::all_of(x.begin(), x.end(), [](char c) { return std::isdigit((unsigned char)c); });
  };
  const auto slash = s.find('/');
  if (slash != std::string::npos) {
    const std::string a = s.substr(0, slash), b = s.substr(slash + 1);
    if (!digits(a) || !digits(b)) return false;
    mpz_class den(b, 10);
    if (den == 0) return false;
    out = mpq_class(mpz_class(a, 10), den);
  } else {
    std::string mant = s;
    long ex = 0;
    const auto e = s.find_first_of("eE");
    if (e != std::string::npos) {
      mant = s.substr(0, e);
      std::string es = s.substr(e + 1);
      bool neg = false;
      if (!es.empty() && (es[0] == '+' || es[0] == '-')) {
        neg = es[0] == '-';
        es = es.substr(1);
      }
      if (!digits(es) || es.size() > 6) return false;
      ex = std::stol(es) * (neg ? -1 : 1);
    }
    const auto dot = mant.find('.');
    std::string ip = mant, fp;
    if (dot != std::string::npos) {
      ip = mant.substr(0, dot);
      fp = mant.substr(dot + 1);
    }
    if (ip.empty() && fp.empty()) return false;
    if ((!ip.empty() && !digits(ip)) || (!fp.empty() && !digits(fp))) return false;
    mpz_class num(ip + fp, 10);
    ex -= (long)fp.size();
    mpz_class p10;
    mpz_ui_pow_ui(p10.get_mpz_t(), 10, (unsigned long)std::labs(ex));
    out = ex >= 0 ? mpq_class(num * p10) : mpq_class(num, p10);
  }
  if (pct) out /= 100;
  out.canonicalize();
  return true;
}

std::string fmt_q(const mpq_class& q) { return q.get_str(); }

std::string approx(const mpq_class& q, int digits) {
  if (q == 0) return "0";
  mpf_class f(q, 512);
  mp_exp_t ex;
  std::string s = f.get_str(ex, 10, (size_t)digits);
  std::string sign;
  if (!s.empty() && s[0] == '-') {
    sign = "-";
    s = s.substr(1);
  }
  while ((int)s.size() < digits) s += '0';
  const long e = (long)ex - 1;
  std::string out = sign + s.substr(0, 1);
  if (digits > 1) out += "." + s.substr(1);
  out += "e";
  out += e < 0 ? "-" : "+";
  const std::string es = std::to_string(std::labs(e));
  out += (es.size() < 2 ? "0" : "") + es;
  return out;
}

std::string json_str(const std::string& s) {
  std::string o = "\"";
  for (const char c : s) {
    if (c == '"') o += "\\\"";
    else if (c == '\\') o += "\\\\";
    else if (c == '\n') o += "\\n";
    else if (c == '\t') o += "\\t";
    else if ((unsigned char)c < 0x20) {
      char buf[8];
      snprintf(buf, sizeof buf, "\\u%04x", c);
      o += buf;
    } else o += c;
  }
  return o + "\"";
}

namespace {

constexpr int TL = INT_MAX;

struct Enc {
  Model& m;
  Cnf& c;
  std::map<std::vector<int>, int> hash;
  std::vector<int> nlit;
  std::vector<char> st;
  std::vector<int> stack;
  std::set<int> owned;

  Enc(Model& model, Cnf& cnf) : m(model), c(cnf), nlit(model.nodes.size(), 0), st(model.nodes.size(), 0) {}

  int fresh() { return ++c.nvars; }
  void cl(std::vector<int> v) { c.clauses.push_back(std::move(v)); }

  int mk_and(std::vector<int> ls) {
    std::vector<int> v;
    for (const int l : ls) {
      if (l == -TL) return -TL;
      if (l != TL) v.push_back(l);
    }
    std::sort(v.begin(), v.end(), [](int a, int b) { return std::abs(a) != std::abs(b) ? std::abs(a) < std::abs(b) : a < b; });
    v.erase(std::unique(v.begin(), v.end()), v.end());
    for (size_t k = 1; k < v.size(); k++)
      if (v[k] == -v[k - 1]) return -TL;
    if (v.empty()) return TL;
    if (v.size() == 1) return v[0];
    std::vector<int> key{0};
    key.insert(key.end(), v.begin(), v.end());
    auto it = hash.find(key);
    if (it != hash.end()) return it->second;
    const int x = fresh();
    std::vector<int> big{x};
    for (const int l : v) {
      cl({-x, l});
      big.push_back(-l);
    }
    cl(big);
    hash[key] = x;
    return x;
  }
  int mk_or(std::vector<int> ls) {
    for (int& l : ls) l = -l;
    return -mk_and(std::move(ls));
  }
  int mk_xor(int a, int b) {
    if (std::abs(a) == TL) return a == TL ? -b : b;
    if (std::abs(b) == TL) return b == TL ? -a : a;
    bool neg = false;
    if (a < 0) {
      a = -a;
      neg = !neg;
    }
    if (b < 0) {
      b = -b;
      neg = !neg;
    }
    if (a == b) return neg ? TL : -TL;
    if (a > b) std::swap(a, b);
    const std::vector<int> key{1, a, b};
    int x;
    auto it = hash.find(key);
    if (it != hash.end()) x = it->second;
    else {
      x = fresh();
      cl({-x, a, b});
      cl({-x, -a, -b});
      cl({x, -a, b});
      cl({x, a, -b});
      hash[key] = x;
    }
    return neg ? -x : x;
  }
  int atleast(const std::vector<int>& xs, int k) {
    if (k <= 0) return TL;
    if (k > (int)xs.size()) return -TL;
    std::vector<int> s(k + 1, -TL);
    s[0] = TL;
    for (const int x : xs)
      for (int j = k; j >= 1; j--) s[j] = mk_or({s[j], mk_and({x, s[j - 1]})});
    return s[k];
  }

  int lit(int e) {
    const Expr& x = m.exprs[e];
    switch (x.op) {
      case Op::Const: return x.val ? TL : -TL;
      case Op::Ref:
      case Op::Flip: return node(x.node);
      case Op::Not: return -lit(x.kids[0]);
      case Op::And: return mk_and({lit(x.kids[0]), lit(x.kids[1])});
      case Op::Or: return mk_or({lit(x.kids[0]), lit(x.kids[1])});
      case Op::Xor: return mk_xor(lit(x.kids[0]), lit(x.kids[1]));
      case Op::Iff: return -mk_xor(lit(x.kids[0]), lit(x.kids[1]));
      case Op::Implies: return mk_or({-lit(x.kids[0]), lit(x.kids[1])});
      case Op::Ite: {
        const int cc = lit(x.kids[0]), a = lit(x.kids[1]), b = lit(x.kids[2]);
        return mk_or({mk_and({cc, a}), mk_and({-cc, b})});
      }
      case Op::AtLeast:
      case Op::AtMost:
      case Op::Exactly: {
        std::vector<int> xs;
        for (const int k : x.kids) xs.push_back(lit(k));
        if (x.op == Op::AtLeast) return atleast(xs, x.k);
        if (x.op == Op::AtMost) return -atleast(xs, x.k + 1);
        return mk_and({atleast(xs, x.k), -atleast(xs, x.k + 1)});
      }
    }
    return -TL;
  }

  int named(int l) {
    if (l > 0 && l != TL && !owned.count(l) && l > (int)c.random.size()) {
      owned.insert(l);
      return l;
    }
    const int v = fresh();
    if (l == TL) cl({v});
    else if (l == -TL) cl({-v});
    else {
      cl({-v, l});
      cl({v, -l});
    }
    return v;
  }

  int node(int n) {
    if (m.is_random(n)) return nlit[n];
    if (st[n] == 2) return nlit[n];
    if (st[n] == 1) {
      std::string path;
      auto it = std::find(stack.begin(), stack.end(), n);
      for (; it != stack.end(); ++it) path += m.nodes[*it].name + " -> ";
      path += m.nodes[n].name;
      m.errors.push_back({m.nodes[n].loc, "cycle: " + path});
      return -TL;
    }
    st[n] = 1;
    stack.push_back(n);
    const int l = m.nodes[n].expr >= 0 ? lit(m.nodes[n].expr) : -TL;
    stack.pop_back();
    if (st[n] == 2) return nlit[n];
    const int v = named(l);
    nlit[n] = v;
    m.nodes[n].var = v;
    st[n] = 2;
    return v;
  }

  int stmt_lit(int e) {
    const int l = lit(e);
    if (std::abs(l) != TL) return l;
    const int v = fresh();
    cl({l == TL ? v : -v});
    return v;
  }

  void run() {
    for (size_t n = 0; n < m.nodes.size(); n++) {
      if (!m.is_random((int)n)) continue;
      const int v = fresh();
      nlit[n] = v;
      m.nodes[n].var = v;
      c.random.push_back(v);
      const Node& x = m.nodes[n];
      c.w[v] = x.kind == Kind::ChoiceValue ? std::make_pair(x.p, mpq_class(1)) : std::make_pair(x.p, mpq_class(1 - x.p));
    }
    for (size_t n = 0; n < m.nodes.size(); n++)
      if (!m.is_random((int)n)) node((int)n);
    for (const Choice& ch : m.choices) {
      std::vector<int> any;
      for (const int v : ch.values) any.push_back(m.nodes[v].var);
      cl(any);
      for (size_t a = 0; a < any.size(); a++)
        for (size_t b = a + 1; b < any.size(); b++) cl({-any[a], -any[b]});
    }
    for (Stmt& s : m.stmts) {
      if (s.kind == SKind::Constraint) {
        const int l = lit(s.expr);
        s.lit = l;
        if (l == TL) continue;
        if (l == -TL) {
          const int v = fresh();
          cl({v});
          cl({-v});
          s.lit = v;
        } else cl({l});
        continue;
      }
      s.lit = stmt_lit(s.expr);
      for (const int g : s.given) s.given_lits.push_back(stmt_lit(g));
    }
    std::set<int> show(c.random.begin(), c.random.end());
    for (const Node& x : m.nodes)
      if (x.var) show.insert(x.var);
    for (const Stmt& s : m.stmts) {
      if (s.kind == SKind::Constraint) continue;
      show.insert(std::abs(s.lit));
      for (const int g : s.given_lits) show.insert(std::abs(g));
    }
    c.show.assign(show.begin(), show.end());
  }
};

void deps(const Model& m, int e, int pol, std::map<int, int>& out) {
  const Expr& x = m.exprs[e];
  switch (x.op) {
    case Op::Const: return;
    case Op::Ref:
    case Op::Flip:
      if (x.node >= 0) out[x.node] |= pol;
      return;
    case Op::Not: deps(m, x.kids[0], pol == 1 ? 2 : pol == 2 ? 1 : pol, out); return;
    case Op::Implies:
      deps(m, x.kids[0], pol == 1 ? 2 : pol == 2 ? 1 : pol, out);
      deps(m, x.kids[1], pol, out);
      return;
    case Op::Xor:
    case Op::Iff:
      for (const int k : x.kids) deps(m, k, 3, out);
      return;
    case Op::Ite:
      deps(m, x.kids[0], 3, out);
      deps(m, x.kids[1], pol, out);
      deps(m, x.kids[2], pol, out);
      return;
    default:
      for (const int k : x.kids) deps(m, k, pol, out);
  }
}

const char* kind_name(Kind k) {
  switch (k) {
    case Kind::Event: return "event";
    case Kind::ChoiceValue: return "value";
    case Kind::Flip: return "flip";
    case Kind::Def: return "def";
  }
  return "?";
}

}  // namespace

Cnf encode(Model& m) {
  Cnf c;
  Enc e(m, c);
  e.run();
  return c;
}

void write_cnf(std::ostream& o, const Model& m, const Cnf& c) {
  o << "c t pwmc\n";
  o << "p cnf " << c.nvars << " " << c.clauses.size() << "\n";
  for (const auto& [v, w] : c.w) {
    o << "c p weight " << v << " " << fmt_q(w.first) << " 0\n";
    o << "c p weight " << -v << " " << fmt_q(w.second) << " 0\n";
  }
  o << "c p show";
  for (const int v : c.show) o << " " << v;
  o << " 0\n";
  for (const Node& x : m.nodes)
    if (x.var) o << "c " << x.var << " " << x.name << "\n";
  int qi = 0;
  for (const Stmt& s : m.stmts) {
    if (s.kind == SKind::Constraint) continue;
    qi++;
    o << "c gm " << (s.kind == SKind::Query ? "query " : "require ") << qi << " " << s.lit;
    if (!s.given_lits.empty()) {
      o << " given";
      for (const int g : s.given_lits) o << " " << g;
    }
    if (s.kind == SKind::Require) o << " " << s.op << " " << fmt_q(s.bound);
    o << "\nc gm text " << qi << " " << s.text << "\n";
  }
  for (const auto& cl : c.clauses) {
    for (const int l : cl) o << l << " ";
    o << "0\n";
  }
}

void write_json(std::ostream& o, const Model& m, const Cnf* c) {
  o << "{\"errors\":[";
  for (size_t i = 0; i < m.errors.size(); i++)
    o << (i ? "," : "") << "{\"line\":" << m.errors[i].loc.line << ",\"col\":" << m.errors[i].loc.col
      << ",\"msg\":" << json_str(m.errors[i].msg) << "}";
  o << "],\"nvars\":" << (c ? c->nvars : 0) << ",\"nodes\":[";
  for (size_t n = 0; n < m.nodes.size(); n++) {
    const Node& x = m.nodes[n];
    o << (n ? "," : "") << "\n{\"name\":" << json_str(x.name) << ",\"label\":" << json_str(x.label)
      << ",\"kind\":\"" << kind_name(x.kind) << "\",\"var\":" << x.var << ",\"line\":" << x.loc.line
      << ",\"anon\":" << (x.anon ? "true" : "false");
    if (m.is_random((int)n)) o << ",\"p\":" << json_str(fmt_q(x.p));
    if (x.choice >= 0) o << ",\"choice\":" << json_str(m.choices[x.choice].name);
    if (!x.text.empty()) o << ",\"text\":" << json_str(x.text);
    o << ",\"deps\":[";
    if (x.expr >= 0) {
      std::map<int, int> d;
      deps(m, x.expr, 1, d);
      bool first = true;
      for (const auto& [dn, pol] : d) {
        o << (first ? "" : ",") << "{\"name\":" << json_str(m.nodes[dn].name) << ",\"pol\":\""
          << (pol == 1 ? "+" : pol == 2 ? "-" : "*") << "\"}";
        first = false;
      }
    }
    o << "]}";
  }
  o << "],\"choices\":[";
  for (size_t i = 0; i < m.choices.size(); i++) {
    const Choice& ch = m.choices[i];
    o << (i ? "," : "") << "{\"name\":" << json_str(ch.name) << ",\"label\":" << json_str(ch.label) << ",\"values\":[";
    for (size_t k = 0; k < ch.values.size(); k++) o << (k ? "," : "") << json_str(m.nodes[ch.values[k]].name);
    o << "]}";
  }
  o << "],\"stmts\":[";
  for (size_t i = 0; i < m.stmts.size(); i++) {
    const Stmt& s = m.stmts[i];
    o << (i ? "," : "") << "\n{\"kind\":\""
      << (s.kind == SKind::Constraint ? "constraint" : s.kind == SKind::Query ? "query" : "require")
      << "\",\"text\":" << json_str(s.text) << ",\"line\":" << s.loc.line << ",\"lit\":" << s.lit << ",\"given\":[";
    for (size_t k = 0; k < s.given_lits.size(); k++) o << (k ? "," : "") << s.given_lits[k];
    o << "]";
    if (s.kind == SKind::Require) o << ",\"op\":" << json_str(s.op) << ",\"bound\":" << json_str(fmt_q(s.bound));
    std::map<int, int> d;
    deps(m, s.expr, 1, d);
    o << ",\"deps\":[";
    bool first = true;
    for (const auto& [dn, pol] : d) {
      o << (first ? "" : ",") << "{\"name\":" << json_str(m.nodes[dn].name) << ",\"pol\":\""
        << (pol == 1 ? "+" : pol == 2 ? "-" : "*") << "\"}";
      first = false;
    }
    o << "]}";
  }
  o << "]}\n";
}

void write_dot(std::ostream& o, const Model& m) {
  o << "digraph model {\n  rankdir=TB;\n  node [fontname=\"Helvetica\", fontsize=11];\n";
  for (size_t n = 0; n < m.nodes.size(); n++) {
    const Node& x = m.nodes[n];
    std::string label = x.label.empty() ? x.name : x.label;
    if (m.is_random((int)n)) label += "\\n" + approx(x.p, 3);
    std::string esc;
    for (const char ch : label) esc += ch == '"' ? std::string("\\\"") : std::string(1, ch);
    o << "  n" << n << " [label=\"" << esc << "\", shape="
      << (m.is_random((int)n) ? (x.anon ? "circle" : "ellipse") : "box") << "];\n";
  }
  for (size_t n = 0; n < m.nodes.size(); n++) {
    if (m.nodes[n].expr < 0) continue;
    std::map<int, int> d;
    deps(m, m.nodes[n].expr, 1, d);
    for (const auto& [dn, pol] : d)
      o << "  n" << n << " -> n" << dn << (pol == 2 ? " [style=dashed, arrowhead=odot]" : pol == 3 ? " [style=dotted]" : "")
        << ";\n";
  }
  o << "}\n";
}

bool eval_expr(const Model& m, int e, const std::vector<char>& v) {
  const Expr& x = m.exprs[e];
  auto k = [&](size_t i) { return eval_expr(m, x.kids[i], v); };
  switch (x.op) {
    case Op::Const: return x.val;
    case Op::Ref:
    case Op::Flip: return v[x.node];
    case Op::Not: return !k(0);
    case Op::And: return k(0) && k(1);
    case Op::Or: return k(0) || k(1);
    case Op::Xor: return k(0) != k(1);
    case Op::Iff: return k(0) == k(1);
    case Op::Implies: return !k(0) || k(1);
    case Op::Ite: return k(0) ? k(1) : k(2);
    case Op::AtLeast:
    case Op::AtMost:
    case Op::Exactly: {
      int cnt = 0;
      for (size_t i = 0; i < x.kids.size(); i++) cnt += k(i);
      return x.op == Op::AtLeast ? cnt >= x.k : x.op == Op::AtMost ? cnt <= x.k : cnt == x.k;
    }
  }
  return false;
}

namespace {
void eval_node(const Model& m, int n, std::vector<char>& v, std::vector<char>& done) {
  if (done[n]) return;
  done[n] = 1;
  if (m.nodes[n].expr < 0) return;
  std::vector<int> todo{m.nodes[n].expr};
  while (!todo.empty()) {
    const int e = todo.back();
    todo.pop_back();
    const Expr& x = m.exprs[e];
    if ((x.op == Op::Ref || x.op == Op::Flip) && x.node >= 0) eval_node(m, x.node, v, done);
    for (const int k : x.kids) todo.push_back(k);
  }
  v[n] = eval_expr(m, m.nodes[n].expr, v);
}
}  // namespace

std::vector<char> eval(const Model& m, const std::vector<char>& random_val) {
  std::vector<char> v(m.nodes.size(), 0), done(m.nodes.size(), 0);
  for (size_t n = 0; n < m.nodes.size(); n++)
    if (m.is_random((int)n)) {
      v[n] = random_val[n];
      done[n] = 1;
    }
  for (size_t n = 0; n < m.nodes.size(); n++) eval_node(m, (int)n, v, done);
  return v;
}

}  // namespace gm
