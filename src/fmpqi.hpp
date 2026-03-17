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

#pragma once

#include "fmpfi.hpp"
#include "mpqi.h"
#include <arjun/arjun.h>
#include <cryptominisat5/solvertypesmini.h>
#include <gmpxx.h>
#include <memory>

class FMpqi final : public CMSat::Field {
public:
    mpqi_t val;

    ~FMpqi() final { mpqi_clear(&val); }
    FMpqi() = delete;

    explicit FMpqi(mpfr_prec_t prec) {
        mpqi_init2(&val, prec);
    }

    explicit FMpqi(long _val, mpfr_prec_t prec) {
        mpqi_init2(&val, prec);
        mpqi_set_d(&val, (double)_val);
    }

    explicit FMpqi(double _val, mpfr_prec_t prec) {
        mpqi_init2(&val, prec);
        mpqi_set_d(&val, _val);
    }

    explicit FMpqi(const mpqi_t& _val) {
        mpqi_init2(&val, _val.prec);
        mpqi_set(&val, &_val);
    }

    explicit FMpqi(const FMpqi& other) {
        mpqi_init2(&val, other.val.prec);
        mpqi_set(&val, &other.val);
    }

    Field& operator=(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        mpqi_set(&val, &od.val);
        return *this;
    }

    Field& operator+=(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        mpqi_add(&val, &val, &od.val);
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        auto ret = std::make_unique<FMpqi>(val);
        mpqi_add(&ret->val, &ret->val, &od.val);
        return ret;
    }

    Field& operator-=(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        mpqi_sub(&val, &val, &od.val);
        return *this;
    }

    Field& operator*=(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        mpqi_mul(&val, &val, &od.val);
        return *this;
    }

    Field& operator/=(const Field& other) final {
        const auto& od = static_cast<const FMpqi&>(other);
        if (mpqi_has_zero(&od.val))
            throw std::runtime_error("Division by mpqi containing zero");
        mpqi_div(&val, &val, &od.val);
        return *this;
    }

    bool operator==(const Field& other) const final {
        const auto& od = static_cast<const FMpqi&>(other);
        if (val.qsize > 0 && od.val.qsize > 0)
            return mpq_equal(val.qval, od.val.qval) != 0;
        if (val.qsize == 0 && od.val.qsize == 0)
            return mpfi_cmp(val.mval, od.val.mval) == 0;
        // Mixed: promote the rational to a finite-precision point interval and
        // compare.  Can return false for mathematically equal values if the
        // interval side was computed with rounding error (i.e. the interval is
        // not a perfect point at the rational value).  This is intentional:
        // once a value has been approximated as an interval it is no longer
        // exactly represented, so strict equality with a rational is impossible.
        const mpqi_t* rat = val.qsize > 0 ? &val : &od.val;
        const mpqi_t* itv = val.qsize > 0 ? &od.val : &val;
        mpfi_t tmp;
        mpfi_init2(tmp, rat->prec);
        mpfi_set_q(tmp, rat->qval);
        bool result = mpfi_cmp(tmp, itv->mval) == 0;
        mpfi_clear(tmp);
        return result;
    }

    std::ostream& display(std::ostream& os) const final {
        if (val.qsize > 0) {
            char* str = mpq_get_str(nullptr, 10, val.qval);
            os << str;
            free(str);
        } else {
            mpfr_prec_t prec = mpfi_get_prec(val.mval);
            mpfr_t left, right;
            mpfr_init2(left, prec);
            mpfr_init2(right, prec);
            mpfi_get_left(left, val.mval);
            mpfi_get_right(right, val.mval);
            char* l_str = nullptr;
            char* r_str = nullptr;
            mpfr_asprintf(&l_str, "%.8Re", left);
            mpfr_asprintf(&r_str, "%.8Re", right);
            os << "[ " << l_str << " " << r_str << " ]";
            mpfr_free_str(l_str);
            mpfr_free_str(r_str);
            mpfr_clear(left);
            mpfr_clear(right);
        }
        return os;
    }

    std::unique_ptr<Field> dup() const final {
        return std::make_unique<FMpqi>(val);
    }

    bool is_zero() const final {
        return mpqi_is_zero(&val);
    }

    bool is_one() const final {
        if (val.qsize > 0)
            return mpq_cmp_ui(val.qval, 1, 1) == 0;
        else
            return mpfi_cmp_si(val.mval, 1) == 0;
    }

    bool parse(const std::string& str, const uint32_t line_no) final {
        uint32_t at = 0;
        ArjunNS::FMpq val_pre;
        if (!val_pre.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        mpqi_set_q(&val, val_pre.get_val().get_mpq_t());
        return true;
    }

    void set_zero() final { mpqi_set_d(&val, 0.0); }
    void set_one() final { mpqi_set_d(&val, 1.0); }

    uint64_t bytes_used() const final {
        if (val.qsize > 0)
            return sizeof(FMpqi) + val.qsize;
        else
            return sizeof(FMpqi) + mpfi_memory_usage(val.mval);
    }
};

class FGenMpqi final : public CMSat::FieldGen {
public:
    mpfr_prec_t prec;
    ~FGenMpqi() final = default;
    FGenMpqi(const FGenMpqi& other) : prec(other.prec) {}
    FGenMpqi& operator=(const FGenMpqi& other) {
        if (this != &other) prec = other.prec;
        return *this;
    }
    explicit FGenMpqi(mpfr_prec_t _prec) : prec(_prec) {}

    std::unique_ptr<CMSat::Field> zero() const final {
        return std::make_unique<FMpqi>((long)0, prec);
    }

    std::unique_ptr<CMSat::Field> one() const final {
        return std::make_unique<FMpqi>((long)1, prec);
    }

    std::unique_ptr<FieldGen> dup() const final {
        return std::make_unique<FGenMpqi>(prec);
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
        const auto& ad = static_cast<const FMpqi&>(a);
        const auto& bd = static_cast<const FMpqi&>(b);
        if (ad.val.qsize > 0 && bd.val.qsize > 0)
            return mpq_cmp(ad.val.qval, bd.val.qval) > 0;
        mpfr_t a_left, b_left;
        mpfr_init2(a_left, prec);
        mpfr_init2(b_left, prec);
        mpqi_left(a_left, &ad.val);
        mpqi_left(b_left, &bd.val);
        bool result = mpfr_cmp(a_left, b_left) > 0;
        mpfr_clear(a_left);
        mpfr_clear(b_left);
        return result;
    }

    bool weighted() const final { return true; }
    bool exact() const final { return false; }
};
