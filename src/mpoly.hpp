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
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_mpoly.h>
#include <flint/mpoly_types.h>
#include <cryptominisat5/solvertypesmini.h>
#include <memory>
#include <cstdlib>

inline void gen_poly() {
    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX);
}

class FPoly : public CMSat::Field {
public:
  fmpq_mpoly_t val;
  std::shared_ptr<fmpq_mpoly_ctx_t> ctx;

  ~FPoly() override {
    fmpq_mpoly_clear(val, ctx.get());
  }

  FPoly(const fmpq_mpoly_t& v,
      const std::shared_ptr<fmpq_mpoly_ctx_t>& c): ctx(c) {
    fmpq_mpoly_init(val, ctx.get());
    fmpq_mpoly_set(val, v, ctx.get());
  }

  FPoly(const FPoly& other): ctx(other.ctx) {
    ctx = other.ctx;
    fmpq_mpoly_init(val, ctx.get());
    fmpq_mpoly_set(val, other.val, ctx.get());
  }
  const auto& get_val() const { return val; }

  Field& operator=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_set(val, od.val, ctx.get());
    return *this;
  }

  Field& operator+=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_add(val, val, od.val, ctx.get());
    return *this;
  }

  Field& operator-=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_sub(val, val, od.val, ctx.get());
    return *this;
  }

  Field& operator*=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_mul(val, val, od.val, ctx.get());
    return *this;
  }

  Field& operator/=(const Field& other) override {
    const auto& od = dynamic_cast<const FPoly&>(other);
    fmpq_mpoly_div(val, val, od.val, ctx.get());
    return *this;
  }

  bool operator==(const Field& other) const override {
      const auto& od = dynamic_cast<const FPoly&>(other);
      return fmpq_mpoly_equal(val, od.val, ctx.get());
  }

  std::ostream& display(std::ostream& os) const override {
  char* str;
  fmpq_mpoly_get_str_pretty(val, (const char**)&str, ctx.get());
  os << *str;
  free(str);
  return os;
  }

  std::unique_ptr<Field> dup() const override {
    return std::make_unique<FPoly>(val, ctx);
  }

  bool is_zero() const override {
    return fmpq_mpoly_is_zero(val, ctx.get());
  }

  bool is_one() const override {
    return fmpq_mpoly_is_one(val, ctx.get());
  }

  bool parse(const std::string& str, const uint32_t line_no) override {
    const char* vars[] = {"x", "y"};
    fmpq_mpoly_set_str_pretty(val, str.c_str(), vars, ctx.get());
    return true;
  }

  void set_zero() override {
    fmpq_mpoly_zero(val, ctx.get());
  }
  void set_one() override {
    fmpq_mpoly_one(val, ctx.get());
  }

  uint64_t bytes_used() const override {
    return 0;
  }
};


class FGenPoly : public CMSat::FieldGen {
public:
  std::shared_ptr<fmpq_mpoly_ctx_t> ctx;

  FGenPoly(const FGenPoly& other) = default;
  ~FGenPoly() override = default;

  std::unique_ptr<CMSat::Field> zero() const override {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx.get());
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<CMSat::Field> one() const override {
    fmpq_mpoly_t val;
    fmpq_mpoly_init(val, ctx.get());
    fmpq_mpoly_one(val, ctx.get());
    return std::make_unique<FPoly>(val, ctx);
  }

  std::unique_ptr<FieldGen> dup() const override {
      return std::make_unique<FGenPoly>(*this);
  }

  bool weighted() const override { return true; }
};

