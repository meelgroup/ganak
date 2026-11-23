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
#include <set>

#include "gmodel.hpp"

namespace gm {

namespace {

enum class T { Name, Str, Num, Sym, NL, End };

struct Tok {
  T t = T::End;
  std::string s;
  Loc loc;
  size_t beg = 0, end = 0;
};

struct ParseError {
  Loc loc;
  std::string msg;
};

const std::set<std::string> keywords = {
  "event", "choice", "rule", "constraint", "query", "require", "given", "if", "then", "else",
  "true", "false", "not", "and", "or", "xor", "flip", "atleast", "atmost", "exactly", "of"};

bool id_start(char c) { return std::isalpha((unsigned char)c) || c == '_'; }
bool id_char(char c) { return std::isalnum((unsigned char)c) || c == '_'; }

std::vector<Tok> lex(const std::string& s, std::vector<Diag>& errs) {
  std::vector<Tok> out;
  size_t i = 0;
  int line = 1;
  size_t line_start = 0;
  auto loc = [&](size_t at) { return Loc{line, (int)(at - line_start) + 1}; };
  auto push = [&](T t, size_t b, size_t e, std::string text) {
    out.push_back(Tok{t, std::move(text), loc(b), b, e});
  };
  while (i < s.size()) {
    const char c = s[i];
    if (c == '\n') {
      push(T::NL, i, i + 1, "newline");
      i++;
      line++;
      line_start = i;
    } else if (c == '#') {
      while (i < s.size() && s[i] != '\n') i++;
    } else if (std::isspace((unsigned char)c)) {
      i++;
    } else if (id_start(c)) {
      const size_t b = i;
      while (i < s.size() && id_char(s[i])) i++;
      while (i + 1 < s.size() && s[i] == '.' && id_start(s[i + 1])) {
        i++;
        while (i < s.size() && id_char(s[i])) i++;
      }
      push(T::Name, b, i, s.substr(b, i - b));
    } else if (std::isdigit((unsigned char)c) || (c == '.' && i + 1 < s.size() && std::isdigit((unsigned char)s[i + 1]))) {
      const size_t b = i;
      while (i < s.size() && std::isdigit((unsigned char)s[i])) i++;
      if (i < s.size() && s[i] == '.') {
        i++;
        while (i < s.size() && std::isdigit((unsigned char)s[i])) i++;
      }
      if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        size_t j = i + 1;
        if (j < s.size() && (s[j] == '+' || s[j] == '-')) j++;
        if (j < s.size() && std::isdigit((unsigned char)s[j])) {
          i = j;
          while (i < s.size() && std::isdigit((unsigned char)s[i])) i++;
        }
      }
      if (i + 1 < s.size() && s[i] == '/' && std::isdigit((unsigned char)s[i + 1])) {
        i++;
        while (i < s.size() && std::isdigit((unsigned char)s[i])) i++;
      }
      if (i < s.size() && s[i] == '%') i++;
      push(T::Num, b, i, s.substr(b, i - b));
    } else if (c == '"') {
      const size_t b = i++;
      std::string v;
      bool closed = false;
      while (i < s.size() && s[i] != '\n') {
        if (s[i] == '\\' && i + 1 < s.size() && (s[i + 1] == '"' || s[i + 1] == '\\')) {
          v += s[i + 1];
          i += 2;
        } else if (s[i] == '"') {
          i++;
          closed = true;
          break;
        } else {
          v += s[i++];
        }
      }
      if (!closed) errs.push_back({loc(b), "unterminated label"});
      push(T::Str, b, i, v);
    } else {
      static const char* syms[] = {"<->", "->", "<=", ">=", "(", ")", "{", "}", ",", ";", "=", "!", "~",
                                   "&", "|", "^", "<", ">"};
      bool found = false;
      for (const char* sym : syms) {
        const size_t n = std::char_traits<char>::length(sym);
        if (s.compare(i, n, sym) == 0) {
          push(T::Sym, i, i + n, sym);
          i += n;
          found = true;
          break;
        }
      }
      if (!found) {
        errs.push_back({loc(i), std::string("unexpected character '") + c + "'"});
        i++;
      }
    }
  }
  out.push_back(Tok{T::End, "end of input", loc(i), i, i});
  return out;
}

size_t edit_distance(const std::string& a, const std::string& b) {
  std::vector<size_t> d(b.size() + 1);
  for (size_t j = 0; j <= b.size(); j++) d[j] = j;
  for (size_t i = 1; i <= a.size(); i++) {
    size_t prev = d[0];
    d[0] = i;
    for (size_t j = 1; j <= b.size(); j++) {
      const size_t cur = d[j];
      d[j] = std::min({d[j] + 1, d[j - 1] + 1, prev + (a[i - 1] == b[j - 1] ? 0 : 1)});
      prev = cur;
    }
  }
  return d[b.size()];
}

struct PendingRule {
  std::string head;
  Loc head_loc;
};

struct Parser {
  const std::string& src;
  std::vector<Tok> toks;
  size_t i = 0;
  int depth = 0;
  bool cont = true;
  size_t last_end = 0;
  Model& m;
  std::vector<PendingRule> pending;

  Parser(const std::string& s, Model& model) : src(s), m(model) { toks = lex(s, m.errors); }

  static bool is_cont(const Tok& t) {
    if (t.t == T::Sym) return t.s != ")" && t.s != "}" && t.s != ";";
    if (t.t == T::Name)
      return t.s == "given" || t.s == "then" || t.s == "else" || t.s == "if" || t.s == "of" || t.s == "not" ||
             t.s == "and" || t.s == "or" || t.s == "xor";
    return false;
  }

  void skip() {
    while (toks[i].t == T::NL && (depth > 0 || cont)) i++;
  }
  const Tok& peek() {
    skip();
    return toks[i];
  }
  Tok take() {
    skip();
    Tok t = toks[i];
    if (t.t == T::End) return t;
    i++;
    cont = is_cont(t);
    if (t.t == T::Sym) {
      if (t.s == "(" || t.s == "{") depth++;
      if ((t.s == ")" || t.s == "}") && depth > 0) depth--;
    }
    last_end = t.end;
    return t;
  }
  [[noreturn]] void fail(const Tok& t, const std::string& msg) { throw ParseError{t.loc, msg}; }
  static std::string what(const Tok& t) {
    if (t.t == T::NL || t.t == T::End) return t.s;
    if (t.t == T::Str) return "label \"" + t.s + "\"";
    return "'" + t.s + "'";
  }
  bool is_sym(const Tok& t, const char* s) const { return t.t == T::Sym && t.s == s; }
  bool is_kw(const Tok& t, const char* s) const { return t.t == T::Name && t.s == s; }
  void expect_sym(const char* s) {
    const Tok& t = peek();
    if (!is_sym(t, s)) fail(t, std::string("expected '") + s + "', found " + what(t));
    take();
  }

  int add_expr(Expr e) {
    m.exprs.push_back(std::move(e));
    return (int)m.exprs.size() - 1;
  }
  int bin(Op op, int a, int b, Loc loc) {
    Expr e;
    e.op = op;
    e.loc = loc;
    e.kids = {a, b};
    return add_expr(std::move(e));
  }

  mpq_class prob(const Tok& t) {
    mpq_class q;
    if (t.t != T::Num || !parse_prob(t.s, q)) fail(t, "expected a probability, found " + what(t));
    if (q < 0 || q > 1) fail(t, "probability " + t.s + " is not between 0 and 1");
    return q;
  }

  std::string decl_name() {
    const Tok t = peek();
    if (t.t != T::Name || keywords.count(t.s)) fail(t, "expected a name, found " + what(t));
    if (t.s.find('.') != std::string::npos) fail(t, "a declared name cannot contain '.': " + t.s);
    take();
    return t.s;
  }

  std::string opt_label() {
    if (peek().t == T::Str) return take().s;
    return "";
  }

  int declare(const std::string& name, Loc loc, Kind kind) {
    auto it = m.by_name.find(name);
    if (it != m.by_name.end())
      throw ParseError{loc, "'" + name + "' is already declared at line " + std::to_string(m.nodes[it->second].loc.line)};
    if (m.choice_by_name.count(name))
      throw ParseError{loc, "'" + name + "' is already declared as a choice at line " +
                                std::to_string(m.choices[m.choice_by_name[name]].loc.line)};
    Node n;
    n.name = name;
    n.kind = kind;
    n.loc = loc;
    m.nodes.push_back(n);
    m.by_name[name] = (int)m.nodes.size() - 1;
    return (int)m.nodes.size() - 1;
  }

  // ---- expressions ----

  int expr() {
    const Tok& t = peek();
    if (is_kw(t, "if")) {
      const Loc loc = take().loc;
      const int c = expr();
      if (!is_kw(peek(), "then")) fail(peek(), "expected 'then', found " + what(peek()));
      take();
      const int a = expr();
      if (!is_kw(peek(), "else")) fail(peek(), "expected 'else', found " + what(peek()));
      take();
      const int b = expr();
      Expr e;
      e.op = Op::Ite;
      e.loc = loc;
      e.kids = {c, a, b};
      return add_expr(std::move(e));
    }
    int l = imp();
    while (is_sym(peek(), "<->")) {
      const Loc loc = take().loc;
      l = bin(Op::Iff, l, imp(), loc);
    }
    return l;
  }
  int imp() {
    const int l = disj();
    if (is_sym(peek(), "->")) {
      const Loc loc = take().loc;
      return bin(Op::Implies, l, imp(), loc);
    }
    return l;
  }
  int disj() {
    int l = xr();
    while (is_sym(peek(), "|") || is_kw(peek(), "or")) {
      const Loc loc = take().loc;
      l = bin(Op::Or, l, xr(), loc);
    }
    return l;
  }
  int xr() {
    int l = conj();
    while (is_sym(peek(), "^") || is_kw(peek(), "xor")) {
      const Loc loc = take().loc;
      l = bin(Op::Xor, l, conj(), loc);
    }
    return l;
  }
  int conj() {
    int l = neg();
    while (is_sym(peek(), "&") || is_kw(peek(), "and")) {
      const Loc loc = take().loc;
      l = bin(Op::And, l, neg(), loc);
    }
    return l;
  }
  int neg() {
    const Tok& t = peek();
    if (is_sym(t, "!") || is_sym(t, "~") || is_kw(t, "not")) {
      const Loc loc = take().loc;
      Expr e;
      e.op = Op::Not;
      e.loc = loc;
      e.kids = {neg()};
      return add_expr(std::move(e));
    }
    return atom();
  }
  int atom() {
    const Tok t = peek();
    Expr e;
    e.loc = t.loc;
    if (is_sym(t, "(")) {
      take();
      const int x = expr();
      expect_sym(")");
      return x;
    }
    if (is_kw(t, "if")) return expr();
    if (t.t != T::Name) fail(t, "expected an expression, found " + what(t));
    if (t.s == "true" || t.s == "false") {
      take();
      e.op = Op::Const;
      e.val = t.s == "true";
      return add_expr(std::move(e));
    }
    if (t.s == "flip") {
      take();
      expect_sym("(");
      const Tok pt = take();
      e.op = Op::Flip;
      e.p = prob(pt);
      expect_sym(")");
      const std::string name = "flip@" + std::to_string(t.loc.line) + ":" + std::to_string(t.loc.col);
      const int n = declare(name, t.loc, Kind::Flip);
      m.nodes[n].p = e.p;
      m.nodes[n].anon = true;
      m.nodes[n].label = "flip(" + pt.s + ")";
      e.node = n;
      return add_expr(std::move(e));
    }
    if (t.s == "atleast" || t.s == "atmost" || t.s == "exactly") {
      take();
      e.op = t.s == "atleast" ? Op::AtLeast : t.s == "atmost" ? Op::AtMost : Op::Exactly;
      const Tok kt = take();
      if (kt.t != T::Num || !std::all_of(kt.s.begin(), kt.s.end(), [](char c) { return std::isdigit((unsigned char)c); }) ||
          kt.s.size() > 6)
        fail(kt, "expected a count after '" + t.s + "', found " + what(kt));
      e.k = std::stoi(kt.s);
      if (is_kw(peek(), "of")) take();
      expect_sym("(");
      e.kids.push_back(expr());
      while (is_sym(peek(), ",")) {
        take();
        e.kids.push_back(expr());
      }
      expect_sym(")");
      return add_expr(std::move(e));
    }
    if (keywords.count(t.s)) fail(t, "expected an expression, found '" + t.s + "'");
    take();
    e.op = Op::Ref;
    e.name = t.s;
    return add_expr(std::move(e));
  }

  // ---- statements ----

  std::string text_from(size_t beg) {
    std::string out;
    bool sp = false;
    for (size_t k = beg; k < last_end && k < src.size(); k++) {
      const char c = src[k];
      if (c == '#') {
        while (k < last_end && src[k] != '\n') k++;
        sp = true;
        continue;
      }
      if (std::isspace((unsigned char)c)) sp = true;
      else {
        if (sp && !out.empty()) out += ' ';
        sp = false;
        out += c;
      }
    }
    return out;
  }

  void givens(Stmt& st) {
    if (!is_kw(peek(), "given")) return;
    take();
    st.given.push_back(expr());
    while (is_sym(peek(), ",")) {
      take();
      st.given.push_back(expr());
    }
  }

  void end_stmt() {
    skip();
    const Tok& t = toks[i];
    if (t.t == T::End) return;
    if (t.t == T::NL || is_sym(t, ";")) {
      i++;
      return;
    }
    fail(t, "unexpected " + what(t) + " (expected the end of the statement)");
  }

  void stmt() {
    const Tok t = peek();
    if (t.t == T::End) return;
    if (is_sym(t, ";")) {
      take();
      return;
    }
    if (t.t != T::Name) fail(t, "expected a statement, found " + what(t));
    if (t.s == "event") {
      take();
      const Tok nt = peek();
      const std::string name = decl_name();
      const std::string label = opt_label();
      const mpq_class p = prob(take());
      const int n = declare(name, nt.loc, Kind::Event);
      m.nodes[n].label = label;
      m.nodes[n].p = p;
    } else if (t.s == "choice") {
      take();
      const Tok nt = peek();
      Choice ch;
      ch.name = decl_name();
      ch.label = opt_label();
      ch.loc = nt.loc;
      if (m.by_name.count(ch.name) || m.choice_by_name.count(ch.name))
        fail(nt, "'" + ch.name + "' is already declared");
      const int ci = (int)m.choices.size();
      m.choice_by_name[ch.name] = ci;
      m.choices.push_back(ch);
      expect_sym("{");
      while (!is_sym(peek(), "}")) {
        const Tok vt = peek();
        const std::string v = decl_name();
        const mpq_class p = prob(take());
        const int n = declare(ch.name + "." + v, vt.loc, Kind::ChoiceValue);
        m.nodes[n].p = p;
        m.nodes[n].choice = ci;
        m.choices[ci].values.push_back(n);
        if (!is_sym(peek(), ",")) break;
        take();
      }
      expect_sym("}");
      if (m.choices[ci].values.empty()) fail(nt, "choice '" + ch.name + "' has no values");
    } else if (t.s == "rule") {
      take();
      const size_t beg = peek().beg;
      const int e = expr();
      const Expr& x = m.exprs[e];
      if (x.op != Op::Implies || m.exprs[x.kids[1]].op != Op::Ref)
        fail(t, "a rule is 'rule BODY -> NAME [probability]'");
      Rule r;
      r.body = x.kids[0];
      r.loc = t.loc;
      r.p = 1;
      const Expr& head = m.exprs[x.kids[1]];
      if (head.name.find('.') != std::string::npos) throw ParseError{head.loc, "a rule cannot define " + head.name};
      if (peek().t == T::Num) r.p = prob(take());
      (void)beg;
      m.rules.push_back(r);
      pending.push_back({head.name, head.loc});
    } else if (t.s == "constraint" || t.s == "query" || t.s == "require") {
      take();
      Stmt st;
      st.kind = t.s == "constraint" ? SKind::Constraint : t.s == "query" ? SKind::Query : SKind::Require;
      st.loc = t.loc;
      const size_t beg = peek().beg;
      st.expr = expr();
      if (st.kind != SKind::Constraint) givens(st);
      st.text = text_from(beg);
      if (st.kind == SKind::Require) {
        const Tok ot = take();
        if (!(ot.t == T::Sym && (ot.s == "<" || ot.s == "<=" || ot.s == ">" || ot.s == ">=")))
          fail(ot, "expected <, <=, > or >= in the requirement, found " + what(ot));
        st.op = ot.s;
        st.bound = prob(take());
      }
      m.stmts.push_back(st);
    } else {
      if (keywords.count(t.s)) fail(t, "expected a statement, found '" + t.s + "'");
      const std::string name = decl_name();
      const std::string label = opt_label();
      const Tok eq = peek();
      if (!is_sym(eq, "="))
        fail(eq, "expected '=' after '" + name + "' (or start the line with event, choice, rule, constraint, query, require)");
      take();
      const size_t beg = peek().beg;
      const int e = expr();
      const int n = declare(name, t.loc, Kind::Def);
      m.nodes[n].label = label;
      m.nodes[n].expr = e;
      m.nodes[n].text = text_from(beg);
    }
    end_stmt();
  }

  void run() {
    while (true) {
      depth = 0;
      cont = true;
      if (peek().t == T::End) break;
      try {
        stmt();
      } catch (const ParseError& e) {
        m.errors.push_back({e.loc, e.msg});
        if (m.errors.size() >= 50) return;
        while (toks[i].t != T::NL && toks[i].t != T::End && !(toks[i].t == T::Sym && toks[i].s == ";")) i++;
        if (toks[i].t != T::End) i++;
      }
    }
  }
};

int ref_expr(Model& m, int node, Loc loc) {
  Expr e;
  e.op = Op::Ref;
  e.node = node;
  e.name = m.nodes[node].name;
  e.loc = loc;
  m.exprs.push_back(e);
  return (int)m.exprs.size() - 1;
}

void resolve(Model& m, const std::vector<PendingRule>& pending) {
  auto err = [&](Loc l, const std::string& s) { m.errors.push_back({l, s}); };
  std::set<int> heads;
  for (size_t r = 0; r < m.rules.size(); r++) {
    const std::string& h = pending[r].head;
    if (m.choice_by_name.count(h)) {
      err(pending[r].head_loc, "a rule cannot define the choice '" + h + "'");
      continue;
    }
    auto it = m.by_name.find(h);
    int n;
    if (it == m.by_name.end()) {
      Node x;
      x.name = h;
      x.kind = Kind::Def;
      x.loc = m.rules[r].loc;
      m.nodes.push_back(x);
      n = (int)m.nodes.size() - 1;
      m.by_name[h] = n;
      heads.insert(n);
    } else {
      n = it->second;
      Node& x = m.nodes[n];
      if (!heads.count(n)) {
        if (x.kind == Kind::Def) {
          err(m.rules[r].loc, "'" + h + "' is defined with '=' at line " + std::to_string(x.loc.line) +
                                  ", so it cannot also have rules");
          continue;
        }
        if (x.kind != Kind::Event || x.anon) {
          err(m.rules[r].loc, "a rule cannot define '" + h + "'");
          continue;
        }
        Node self;
        self.name = h + ".self";
        self.label = (x.label.empty() ? h : x.label) + ", by itself";
        self.kind = Kind::Event;
        self.p = x.p;
        self.loc = x.loc;
        self.anon = true;
        m.nodes.push_back(self);
        const int si = (int)m.nodes.size() - 1;
        m.by_name[self.name] = si;
        m.nodes[n].kind = Kind::Def;
        m.nodes[n].self = si;
        heads.insert(n);
      }
    }
    m.rules[r].head = n;
    m.nodes[n].rules.push_back((int)r);
    if (m.rules[r].p != 1) {
      const Loc l = m.rules[r].loc;
      Node f;
      f.name = "rule@" + std::to_string(l.line) + ":" + std::to_string(l.col);
      f.label = "rule " + h + " (" + approx(m.rules[r].p, 3) + ")";
      f.kind = Kind::Flip;
      f.p = m.rules[r].p;
      f.loc = l;
      f.anon = true;
      m.nodes.push_back(f);
      m.rules[r].flip = (int)m.nodes.size() - 1;
      m.by_name[f.name] = m.rules[r].flip;
    }
  }
  for (const int n : heads) {
    std::vector<int> terms;
    std::string text;
    for (const int r : m.nodes[n].rules) {
      const Rule& ru = m.rules[r];
      int t = ru.body;
      if (ru.flip >= 0) {
        Expr e;
        e.op = Op::And;
        e.loc = ru.loc;
        e.kids = {ru.body, ref_expr(m, ru.flip, ru.loc)};
        m.exprs.push_back(e);
        t = (int)m.exprs.size() - 1;
      }
      terms.push_back(t);
    }
    if (m.nodes[n].self >= 0) terms.push_back(ref_expr(m, m.nodes[n].self, m.nodes[n].loc));
    int acc = terms[0];
    for (size_t k = 1; k < terms.size(); k++) {
      Expr e;
      e.op = Op::Or;
      e.loc = m.nodes[n].loc;
      e.kids = {acc, terms[k]};
      m.exprs.push_back(e);
      acc = (int)m.exprs.size() - 1;
    }
    m.nodes[n].expr = acc;
    m.nodes[n].text = std::to_string(m.nodes[n].rules.size()) + (m.nodes[n].rules.size() == 1 ? " rule" : " rules") +
                      (m.nodes[n].self >= 0 ? " + itself" : "");
  }
  for (Expr& e : m.exprs) {
    if (e.op != Op::Ref || e.node >= 0) continue;
    auto it = m.by_name.find(e.name);
    if (it != m.by_name.end()) {
      e.node = it->second;
      continue;
    }
    if (m.choice_by_name.count(e.name)) {
      const Choice& c = m.choices[m.choice_by_name[e.name]];
      err(e.loc, "'" + e.name + "' is a choice: use one of its values" +
                     (c.values.empty() ? std::string() : ", e.g. " + m.nodes[c.values[0]].name));
      continue;
    }
    std::string best;
    size_t bd = 3;
    for (const auto& [name, idx] : m.by_name) {
      if (m.nodes[idx].anon) continue;
      const size_t d = edit_distance(name, e.name);
      if (d < bd) {
        bd = d;
        best = name;
      }
    }
    err(e.loc, "unknown name '" + e.name + "'" + (best.empty() ? "" : " (did you mean '" + best + "'?)"));
  }
  for (const Choice& c : m.choices) {
    if (c.values.empty()) continue;  // its declaration failed, already reported
    mpq_class s = 0;
    for (const int v : c.values) s += m.nodes[v].p;
    if (s != 1) err(c.loc, "the probabilities of choice '" + c.name + "' sum to " + fmt_q(s) + ", not 1");
  }
}

}  // namespace

Model parse(const std::string& src) {
  Model m;
  Parser p(src, m);
  p.run();
  resolve(m, p.pending);
  std::stable_sort(m.errors.begin(), m.errors.end(), [](const Diag& a, const Diag& b) {
    return a.loc.line != b.loc.line ? a.loc.line < b.loc.line : a.loc.col < b.loc.col;
  });
  return m;
}

}  // namespace gm
