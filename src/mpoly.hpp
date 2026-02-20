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

#pragma once
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/mpoly_types.h>
#include <cryptominisat5/solvertypesmini.h>
#include <memory>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cctype>

struct AutoPoly {
  ~AutoPoly() { fmpq_mpoly_ctx_clear(ctx); }
  AutoPoly(int nvars) {
    fmpq_mpoly_ctx_init(ctx, nvars, ORD_DEGREVLEX);
  }
  fmpq_mpoly_ctx_t ctx;
};

class FPoly final : public CMSat::Field {
public:
  fmpq_mpoly_t val;
  std::shared_ptr<AutoPoly> ctx;

  ~FPoly() final {
    fmpq_mpoly_clear(val, ctx->ctx);
  }

  FPoly(const fmpq_mpoly_t& v,
      const std::shared_ptr<AutoPoly>& c): ctx(c) {
    memcpy(val, v, sizeof(fmpq_mpoly_t));
  }

  FPoly(const FPoly& other): ctx(other.ctx) {
    ctx = other.ctx;
    fmpq_mpoly_init(val, ctx->ctx);
    fmpq_mpoly_set(val, other.val, ctx->ctx);
  }
  const auto& get_val() const { return val; }

  Field& operator=(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_set(val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator+=(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_add(val, val, od.val, ctx->ctx);
    return *this;
  }

  std::unique_ptr<Field> add(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_t v;
    fmpq_mpoly_init(v, ctx->ctx);
    fmpq_mpoly_add(v, val, od.val, ctx->ctx);
    return std::make_unique<FPoly>(v, ctx);
  }

  Field& operator-=(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_sub(val, val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator*=(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_mul(val, val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator/=(const Field& other) final {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_div(val, val, od.val, ctx->ctx);
    return *this;
  }

  bool operator==(const Field& other) const final {
      const auto& od = dynamic_cast<const FPoly&>(other);
      return fmpq_mpoly_equal(val, od.val, ctx->ctx);
  }

  std::ostream& display(std::ostream& os) const final {
    auto varnames = nvars_names();
    char* str = fmpq_mpoly_get_str_pretty(val, (const char**)varnames, ctx->ctx);
    os << str;
    free(str);
    del_nvars_names(varnames);
    return os;
  }

  std::unique_ptr<Field> dup() const final {
    fmpq_mpoly_t v;
    fmpq_mpoly_init(v, ctx->ctx);
    fmpq_mpoly_set(v, val, ctx->ctx);
    return std::make_unique<FPoly>(v, ctx);
  }

  bool is_zero() const final {
    return fmpq_mpoly_is_zero(val, ctx->ctx);
  }

  bool is_one() const final {
    return fmpq_mpoly_is_one(val, ctx->ctx);
  }

  std::string rem_trail_space(const std::string& str) {
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    if (end == std::string::npos) return "";
    return str.substr(0, end + 1);
  }

  void del_nvars_names(char** names) const {
    uint32_t nvars = fmpq_mpoly_ctx_nvars(ctx->ctx);
    for(uint32_t i = 0; i < nvars; i++) delete[] names[i];
    delete[] names;
  }
  char** nvars_names() const {
    uint32_t nvars = fmpq_mpoly_ctx_nvars(ctx->ctx);
    char** vars = new char*[nvars];
    for(uint32_t i = 0; i < nvars; i++) {
      std::string x = "x";
      x += std::to_string(i);
      vars[i] = new char[x.size()+1];
      memcpy(vars[i], x.c_str(), x.size()+1);
    }
    return vars;
  }

  bool parse(const std::string& str, const uint32_t line_no) final {
    auto str2 = rem_trail_space(str);
    if (str2.empty() || str2.back() != '0') {
      std::cerr << "Error parsing polynomial on line " << line_no
        << " -- it should have a 0 at the end?"
        << " -- poly: " << str
        << " -- stripped poly: " << str2 << std::endl;
      exit(EXIT_FAILURE);
    }
    str2.pop_back();
    auto varnames = nvars_names();
    auto ret = fmpq_mpoly_set_str_pretty(val, str2.c_str(), (const char**)varnames, ctx->ctx);
    del_nvars_names(varnames);
    if (ret == -1) {
      std::cerr << "Error parsing polynomial on line " << line_no
        << " -- poly: " << str << std::endl;
      exit(EXIT_FAILURE);
    }
    return true;
  }

  void set_zero() final { fmpq_mpoly_zero(val, ctx->ctx); }
  void set_one() final { fmpq_mpoly_one(val, ctx->ctx); }

  uint64_t bytes_used() const final {
    return sizeof(FPoly) + fmpq_mpoly_length(val, ctx->ctx) * 100;
  }
};

class FGenPoly final : public CMSat::FieldGen {
public:
  std::shared_ptr<AutoPoly> ctx;

  FGenPoly(int nvars) : ctx(std::make_shared<AutoPoly>(nvars)) {}

  std::unique_ptr<CMSat::Field> zero() const final {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx->ctx);
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<CMSat::Field> one() const final {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx->ctx);
    fmpq_mpoly_one(val, ctx->ctx);
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<FieldGen> dup() const final {
      return std::make_unique<FGenPoly>(*this);
  }

  bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
      const auto& ad = dynamic_cast<const FPoly&>(a);
      const auto& bd = dynamic_cast<const FPoly&>(b);
      return fmpq_mpoly_cmp(ad.val, bd.val, ctx->ctx) == -1;
  }

  bool weighted() const final { return true; }
};

// LaurentPoly: multivariate Laurent polynomial over Q.
//
// Representation: (poly, shift) where the mathematical value is
//   z^shift * poly(z),
// with poly having only non-negative exponents (handled by FLINT) and
// shift in Z^m recording the per-variable minimal exponent so that all
// negative exponents are absorbed into shift.
//
// Operations:
//   operator*= : multiply polynomials, add shift vectors.
//   operator+=/-= : align shifts to their componentwise minimum, then
//                   multiply each poly by the corresponding correction
//                   monomial before adding/subtracting.
//
// Parse/display format uses variable names z0, z1, z2, ...
// Exponents may be negative; parentheses around them are optional.
// Example: "1/2*z0^2 - 3*z1^-1 + z0^(-2) 0"
// (The trailing '0' follows the same convention as FPoly.)
class LaurentPoly final : public CMSat::Field {
public:
  fmpq_mpoly_t poly;         // non-negative-exponent part
  std::vector<slong> shift;  // shift[i] = minimal exponent of z_i
  std::shared_ptr<AutoPoly> ctx;

  explicit LaurentPoly(const std::shared_ptr<AutoPoly>& c)
      : shift((size_t)fmpq_mpoly_ctx_nvars(c->ctx), 0), ctx(c) {
    fmpq_mpoly_init(poly, ctx->ctx);
  }

  LaurentPoly(const LaurentPoly& o) : shift(o.shift), ctx(o.ctx) {
    fmpq_mpoly_init(poly, ctx->ctx);
    fmpq_mpoly_set(poly, o.poly, ctx->ctx);
  }

  ~LaurentPoly() final { fmpq_mpoly_clear(poly, ctx->ctx); }

  int nvars() const { return (int)shift.size(); }

  // Multiply the internal polynomial (in-place) by z^d, where d[i] >= 0.
  void mul_by_nonneg_monomial(const std::vector<slong>& d) {
    int nv = nvars();
    std::vector<ulong> ue(nv);
    for (int i = 0; i < nv; i++) ue[i] = (ulong)d[i];
    fmpq_mpoly_t mono;
    fmpq_mpoly_init(mono, ctx->ctx);
    fmpq_t one;
    fmpq_init(one);
    fmpq_one(one);
    fmpq_mpoly_set_coeff_fmpq_ui(mono, one, ue.data(), ctx->ctx);
    fmpq_clear(one);
    fmpq_mpoly_mul(poly, poly, mono, ctx->ctx);
    fmpq_mpoly_clear(mono, ctx->ctx);
  }

  static bool all_zero(const std::vector<slong>& v) {
    for (auto x : v) if (x != 0) return false;
    return true;
  }

  Field& operator=(const Field& other) final {
    const auto& od = dynamic_cast<const LaurentPoly&>(other);
    if (this == &od) return *this;
    shift = od.shift;
    fmpq_mpoly_set(poly, od.poly, ctx->ctx);
    return *this;
  }

  Field& operator+=(const Field& other) final {
    const auto& od = dynamic_cast<const LaurentPoly&>(other);
    int nv = nvars();
    std::vector<slong> ns(nv), d1(nv), d2(nv);
    for (int i = 0; i < nv; i++) {
      ns[i] = std::min(shift[i], od.shift[i]);
      d1[i] = shift[i] - ns[i];    // >= 0
      d2[i] = od.shift[i] - ns[i]; // >= 0
    }
    if (!all_zero(d1)) mul_by_nonneg_monomial(d1);
    if (!all_zero(d2)) {
      LaurentPoly tmp(od);
      tmp.mul_by_nonneg_monomial(d2);
      fmpq_mpoly_add(poly, poly, tmp.poly, ctx->ctx);
    } else {
      fmpq_mpoly_add(poly, poly, od.poly, ctx->ctx);
    }
    shift = ns;
    return *this;
  }

  std::unique_ptr<Field> add(const Field& other) final {
    auto res = std::make_unique<LaurentPoly>(*this);
    *res += other;
    return res;
  }

  Field& operator-=(const Field& other) final {
    const auto& od = dynamic_cast<const LaurentPoly&>(other);
    int nv = nvars();
    std::vector<slong> ns(nv), d1(nv), d2(nv);
    for (int i = 0; i < nv; i++) {
      ns[i] = std::min(shift[i], od.shift[i]);
      d1[i] = shift[i] - ns[i];
      d2[i] = od.shift[i] - ns[i];
    }
    if (!all_zero(d1)) mul_by_nonneg_monomial(d1);
    if (!all_zero(d2)) {
      LaurentPoly tmp(od);
      tmp.mul_by_nonneg_monomial(d2);
      fmpq_mpoly_sub(poly, poly, tmp.poly, ctx->ctx);
    } else {
      fmpq_mpoly_sub(poly, poly, od.poly, ctx->ctx);
    }
    shift = ns;
    return *this;
  }

  // multiply: (p1,s1)*(p2,s2) = (p1*p2, s1+s2)
  Field& operator*=(const Field& other) final {
    const auto& od = dynamic_cast<const LaurentPoly&>(other);
    fmpq_mpoly_mul(poly, poly, od.poly, ctx->ctx);
    for (int i = 0; i < nvars(); i++) shift[i] += od.shift[i];
    return *this;
  }

  Field& operator/=(const Field& other) final {
    const auto& od = dynamic_cast<const LaurentPoly&>(other);
    fmpq_mpoly_div(poly, poly, od.poly, ctx->ctx);
    for (int i = 0; i < nvars(); i++) shift[i] -= od.shift[i];
    return *this;
  }

  // Equality via subtraction: handles non-normalised representations.
  bool operator==(const Field& other) const final {
    LaurentPoly diff(*this);
    diff -= other;
    return diff.is_zero();
  }

  bool is_zero() const final { return fmpq_mpoly_is_zero(poly, ctx->ctx); }

  bool is_one() const final {
    for (auto s : shift) if (s != 0) return false;
    return fmpq_mpoly_is_one(poly, ctx->ctx);
  }

  void set_zero() final {
    fmpq_mpoly_zero(poly, ctx->ctx);
    std::fill(shift.begin(), shift.end(), 0);
  }

  void set_one() final {
    fmpq_mpoly_one(poly, ctx->ctx);
    std::fill(shift.begin(), shift.end(), 0);
  }

  std::unique_ptr<Field> dup() const final {
    return std::make_unique<LaurentPoly>(*this);
  }

  uint64_t bytes_used() const final {
    return sizeof(LaurentPoly)
        + fmpq_mpoly_length(poly, ctx->ctx) * 100
        + shift.size() * sizeof(slong);
  }

  // Display as: c1*z0^e0*z1^e1 + c2*z0^f0 - ...
  // Coefficients ±1 are suppressed unless the term is a constant.
  // Exponents of 1 are suppressed.  Negative exponents are printed as-is.
  std::ostream& display(std::ostream& os) const final {
    int nv = nvars();
    slong nterms = fmpq_mpoly_length(poly, ctx->ctx);
    if (nterms == 0) { os << "0"; return os; }

    fmpq_t coeff, minus_one;
    fmpq_init(coeff);
    fmpq_init(minus_one);
    fmpq_set_si(minus_one, -1, 1);
    std::vector<ulong> uexp(nv);

    for (slong t = 0; t < nterms; t++) {
      fmpq_mpoly_get_term_coeff_fmpq(coeff, poly, t, ctx->ctx);
      fmpq_mpoly_get_term_exp_ui(uexp.data(), poly, t, ctx->ctx);

      // Actual Laurent exponents (poly exponent + shift)
      bool has_vars = false;
      std::vector<slong> aexp(nv);
      for (int i = 0; i < nv; i++) {
        aexp[i] = (slong)uexp[i] + shift[i];
        if (aexp[i] != 0) has_vars = true;
      }

      bool is_neg = (fmpq_sgn(coeff) < 0);
      bool abs_is_one = fmpq_is_one(coeff) || fmpq_equal(coeff, minus_one);

      if (t > 0) os << (is_neg ? " - " : " + ");
      else if (is_neg) os << "-";

      // Print |coeff| unless it's 1 and there are variable factors
      if (!abs_is_one || !has_vars) {
        fmpq_t ac; fmpq_init(ac);
        fmpq_abs(ac, coeff);
        char* cs = fmpq_get_str(nullptr, 10, ac);
        os << cs;
        free(cs);
        fmpq_clear(ac);
        if (has_vars) os << "*";
      }

      bool first_v = true;
      for (int i = 0; i < nv; i++) {
        if (aexp[i] == 0) continue;
        if (!first_v) os << "*";
        first_v = false;
        os << "z" << i;
        if (aexp[i] != 1) os << "^" << aexp[i];
      }
    }

    fmpq_clear(coeff);
    fmpq_clear(minus_one);
    return os;
  }

  // Parse a Laurent polynomial string.
  //
  // Grammar (whitespace freely permitted between tokens):
  //   expr     ::= [sign] term { sign term }*
  //   term     ::= [rational '*'] var_product
  //              | rational                   (constant)
  //   rational ::= integer ['/' integer]
  //   var_product ::= var_power { '*' var_power }*
  //   var_power   ::= 'z' index ['^' ['('] [sign] integer [')'] ]
  //   sign        ::= '+' | '-'
  //   index       ::= digit+
  //   integer     ::= digit+
  //
  // The string must end with a trailing '0' (stripped before parsing),
  // following the same convention as FPoly.
  bool parse(const std::string& str, const uint32_t line_no) final {
    std::string s = rem_trail_space(str);
    if (s.empty() || s.back() != '0') {
      std::cerr << "Error parsing Laurent polynomial on line " << line_no
          << ": expected trailing '0', got: '" << str << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
    s.pop_back();
    s = rem_trail_space(s);

    fmpq_mpoly_zero(poly, ctx->ctx);
    std::fill(shift.begin(), shift.end(), 0);
    if (s.empty()) return true;

    int nv = nvars();

    struct RawTerm {
      slong num = 1, den = 1;
      std::vector<slong> exp;
    };
    std::vector<RawTerm> terms;

    // Helper: skip spaces/tabs
    auto skip_ws = [&](size_t& p) {
      while (p < s.size() && (s[p] == ' ' || s[p] == '\t')) ++p;
    };

    // Helper: parse a run of decimal digits into a non-negative slong
    auto parse_digits = [&](size_t& p, slong& val) -> bool {
      if (p >= s.size() || !isdigit((unsigned char)s[p])) return false;
      val = 0;
      while (p < s.size() && isdigit((unsigned char)s[p]))
        val = val * 10 + (s[p++] - '0');
      return true;
    };

    size_t pos = 0;
    while (pos < s.size()) {
      skip_ws(pos);
      if (pos >= s.size()) break;

      RawTerm term;
      term.exp.resize(nv, 0);

      // Leading sign
      int sign = 1;
      if (s[pos] == '+')      { ++pos; }
      else if (s[pos] == '-') { sign = -1; ++pos; }
      skip_ws(pos);

      // Try to parse a numeric coefficient (integer or p/q)
      bool have_coeff = false;
      if (pos < s.size() && isdigit((unsigned char)s[pos])) {
        size_t rewind = pos;
        slong n; parse_digits(pos, n);

        if (pos < s.size() && s[pos] == '/') {
          // rational p/q
          ++pos;
          slong d; parse_digits(pos, d);
          term.num = sign * n; term.den = d;
          have_coeff = true; sign = 1;
        } else if (pos >= s.size() || s[pos]=='+' || s[pos]=='-'
                   || s[pos]==' '  || s[pos]=='\t') {
          // standalone integer constant – push as a term with zero exponents
          term.num = sign * n; term.den = 1;
          terms.push_back(term);
          continue;
        } else if (s[pos] == '*') {
          // coefficient followed by '*'
          ++pos;
          term.num = sign * n; term.den = 1;
          have_coeff = true; sign = 1;
        } else {
          // digits were not a coefficient (unexpected char) – rewind
          pos = rewind;
        }
      }

      if (!have_coeff) { term.num = sign; term.den = 1; }

      skip_ws(pos);

      // Parse one or more variable powers: z<n> [^[sign]int] [* ...]
      bool parsed_var = false;
      while (pos < s.size() && s[pos] == 'z') {
        ++pos;
        slong vidx;
        if (!parse_digits(pos, vidx)) {
          std::cerr << "Error: expected variable index after 'z' on line "
              << line_no << std::endl;
          exit(EXIT_FAILURE);
        }
        if (vidx >= nv) {
          std::cerr << "Error: variable z" << vidx << " out of range ("
              << nv << " vars) on line " << line_no << std::endl;
          exit(EXIT_FAILURE);
        }

        slong exp_val = 1;
        if (pos < s.size() && s[pos] == '^') {
          ++pos;
          bool has_paren = (pos < s.size() && s[pos] == '(');
          if (has_paren) ++pos;
          int esign = 1;
          if      (pos < s.size() && s[pos] == '-') { esign = -1; ++pos; }
          else if (pos < s.size() && s[pos] == '+') {              ++pos; }
          slong e;
          if (!parse_digits(pos, e)) {
            std::cerr << "Error: expected exponent after '^' on line "
                << line_no << std::endl;
            exit(EXIT_FAILURE);
          }
          if (has_paren && pos < s.size() && s[pos] == ')') ++pos;
          exp_val = esign * e;
        }

        term.exp[vidx] += exp_val;
        parsed_var = true;

        skip_ws(pos);
        // Consume '*' only if the next non-space token is another variable
        if (pos < s.size() && s[pos] == '*') {
          size_t peek = pos + 1;
          while (peek < s.size() && (s[peek]==' ' || s[peek]=='\t')) ++peek;
          if (peek < s.size() && s[peek] == 'z') {
            ++pos; skip_ws(pos);
          } else {
            break;
          }
        }
      }

      if (!have_coeff && !parsed_var) {
        std::cerr << "Error: empty term in Laurent polynomial on line "
            << line_no << std::endl;
        exit(EXIT_FAILURE);
      }

      terms.push_back(term);
    }

    if (terms.empty()) return true;

    // Compute per-variable shift = componentwise minimum of all exponents
    std::vector<slong> new_shift(nv, 0);
    for (auto& t : terms)
      for (int i = 0; i < nv; i++)
        new_shift[i] = std::min(new_shift[i], t.exp[i]);

    // Build the FLINT polynomial with non-negative adjusted exponents.
    // Use get+add+set to combine any like terms.
    std::vector<ulong> ue(nv);
    fmpq_t cur_c, add_c;
    fmpq_init(cur_c); fmpq_init(add_c);
    for (auto& t : terms) {
      for (int i = 0; i < nv; i++) ue[i] = (ulong)(t.exp[i] - new_shift[i]);
      fmpq_mpoly_get_coeff_fmpq_ui(cur_c, poly, ue.data(), ctx->ctx);
      fmpq_set_si(add_c, t.num, (ulong)t.den);
      fmpq_add(cur_c, cur_c, add_c);
      fmpq_mpoly_set_coeff_fmpq_ui(poly, cur_c, ue.data(), ctx->ctx);
    }
    fmpq_clear(cur_c); fmpq_clear(add_c);

    shift = new_shift;
    return true;
  }

private:
  std::string rem_trail_space(const std::string& str) {
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    if (end == std::string::npos) return "";
    return str.substr(0, end + 1);
  }
};

class LaurentPolyGen final : public CMSat::FieldGen {
public:
  std::shared_ptr<AutoPoly> ctx;

  explicit LaurentPolyGen(int nvars)
      : ctx(std::make_shared<AutoPoly>(nvars)) {}

  std::unique_ptr<CMSat::Field> zero() const final {
    return std::make_unique<LaurentPoly>(ctx);
  }

  std::unique_ptr<CMSat::Field> one() const final {
    auto p = std::make_unique<LaurentPoly>(ctx);
    p->set_one();
    return p;
  }

  std::unique_ptr<FieldGen> dup() const final {
    return std::make_unique<LaurentPolyGen>(*this);
  }

  bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
    const auto& ad = dynamic_cast<const LaurentPoly&>(a);
    const auto& bd = dynamic_cast<const LaurentPoly&>(b);
    return fmpq_mpoly_cmp(ad.poly, bd.poly, ctx->ctx) == -1;
  }

  bool weighted() const final { return true; }
};
