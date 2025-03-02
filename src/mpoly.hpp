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

struct AutoPoly {
  ~AutoPoly() { fmpq_mpoly_ctx_clear(ctx); }
  AutoPoly(int nvars) {
    fmpq_mpoly_ctx_init(ctx, nvars, ORD_DEGREVLEX);
  }
  fmpq_mpoly_ctx_t ctx;
};

class FPoly : public CMSat::Field {
public:
  fmpq_mpoly_t val;
  std::shared_ptr<AutoPoly> ctx;

  ~FPoly() override {
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

  Field& operator=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_set(val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator+=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_add(val, val, od.val, ctx->ctx);
    return *this;
  }

  std::unique_ptr<Field> add(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_t v;
    fmpq_mpoly_init(v, ctx->ctx);
    fmpq_mpoly_add(v, val, od.val, ctx->ctx);
    return std::make_unique<FPoly>(v, ctx);
  }

  Field& operator-=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_sub(val, val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator*=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_mul(val, val, od.val, ctx->ctx);
    return *this;
  }

  Field& operator/=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_div(val, val, od.val, ctx->ctx);
    return *this;
  }

  bool operator==(const Field& other) const override {
      const auto& od = dynamic_cast<const FPoly&>(other);
      return fmpq_mpoly_equal(val, od.val, ctx->ctx);
  }

  std::ostream& display(std::ostream& os) const override {
    auto varnames = nvars_names();
    char* str = fmpq_mpoly_get_str_pretty(val, (const char**)varnames, ctx->ctx);
    os << str;
    free(str);
    del_nvars_names(varnames);
    return os;
  }

  std::unique_ptr<Field> dup() const override {
    fmpq_mpoly_t v;
    fmpq_mpoly_init(v, ctx->ctx);
    fmpq_mpoly_set(v, val, ctx->ctx);
    return std::make_unique<FPoly>(v, ctx);
  }

  bool is_zero() const override {
    return fmpq_mpoly_is_zero(val, ctx->ctx);
  }

  bool is_one() const override {
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

  bool parse(const std::string& str, const uint32_t line_no) override {
    auto str2 = rem_trail_space(str);
    if (str2.empty() || str2.back() != '0') {
      std::cerr << "Error parsing polynomial on line " << line_no
        << " -- it should have a 0 at the end?"
        << " -- poly: " << str
        << " -- stripped poly: " << str2 << std::endl;
      exit(-1);
    }
    str2.pop_back();
    auto varnames = nvars_names();
    auto ret = fmpq_mpoly_set_str_pretty(val, str2.c_str(), (const char**)varnames, ctx->ctx);
    del_nvars_names(varnames);
    if (ret == -1) {
      std::cerr << "Error parsing polynomial on line " << line_no
        << " -- poly: " << str << std::endl;
      exit(-1);
    }
    return true;
  }

  void set_zero() override {
    fmpq_mpoly_zero(val, ctx->ctx);
  }
  void set_one() override {
    fmpq_mpoly_one(val, ctx->ctx);
  }

  uint64_t bytes_used() const override {
    return 0;
  }
};

class FGenPoly : public CMSat::FieldGen {
public:
  std::shared_ptr<AutoPoly> ctx;

  FGenPoly(int nvars) : ctx(std::make_shared<AutoPoly>(nvars)) {}

  std::unique_ptr<CMSat::Field> zero() const override {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx->ctx);
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<CMSat::Field> one() const override {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx->ctx);
    fmpq_mpoly_one(val, ctx->ctx);
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<FieldGen> dup() const override {
      return std::make_unique<FGenPoly>(*this);
  }

  bool weighted() const override { return true; }
};
