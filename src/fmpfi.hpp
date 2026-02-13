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

#include <arjun/arjun.h>
#include <cryptominisat5/solvertypesmini.h>
#include <gmpxx.h>
#include <memory>
#include <mpfi.h>

inline unsigned int mpfi_memory_usage(const mpfi_t& x) {
    mpfr_prec_t precision = mpfi_get_prec(x);
    size_t limb_size = mp_bits_per_limb;
    size_t num_limbs = (precision + limb_size - 1) / limb_size;
    // Two MPFR endpoints (left and right)
    unsigned int base_size = 2 * sizeof(__mpfr_struct);
    unsigned int mantissa_size = 2 * num_limbs * sizeof(mp_limb_t);
    return base_size + mantissa_size;
}

class FMpfi final : public CMSat::Field {
public:
    uint16_t prec;
    mpfi_t val;

    ~FMpfi() final { mpfi_clear(val); }
    FMpfi() = delete;

    explicit FMpfi(uint16_t _prec) : prec(_prec) {
        mpfi_init2(val, prec);
        mpfi_set_si(val, 0);
    }

    explicit FMpfi(const int _val, uint16_t _prec) : prec(_prec) {
        mpfi_init2(val, prec);
        mpfi_set_si(val, _val);
    }

    explicit FMpfi(const double _val, uint16_t _prec) : prec(_prec) {
        mpfi_init2(val, prec);
        mpfi_set_d(val, _val);
    }

    explicit FMpfi(const mpfi_t& _val, uint16_t _prec) : prec(_prec) {
        mpfi_init2(val, prec);
        mpfi_set(val, _val);
    }

    explicit FMpfi(const FMpfi& other, uint16_t _prec) : prec(_prec) {
        mpfi_init2(val, prec);
        mpfi_set(val, other.val);
    }

    const mpfi_t& get_val() const { return val; }

    Field& operator=(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        mpfi_set(val, od.val);
        return *this;
    }

    Field& operator+=(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        mpfi_add(val, val, od.val);
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        mpfi_t res;
        mpfi_init2(res, prec);
        mpfi_add(res, val, od.val);
        auto ret = std::make_unique<FMpfi>(res, prec);
        mpfi_clear(res);
        return ret;
    }

    Field& operator-=(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        mpfi_sub(val, val, od.val);
        return *this;
    }

    Field& operator*=(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        mpfi_mul(val, val, od.val);
        return *this;
    }

    Field& operator/=(const Field& other) final {
        const auto& od = static_cast<const FMpfi&>(other);
        if (mpfi_has_zero(od.val))
            throw std::runtime_error("Division by interval containing zero");
        mpfi_div(val, val, od.val);
        return *this;
    }

    bool operator==(const Field& other) const final {
        const auto& od = static_cast<const FMpfi&>(other);
        return mpfi_cmp(val, od.val) == 0;
    }

    std::ostream& display(std::ostream& os) const final {
        mpfr_t left, right;
        mpfr_init2(left, prec);
        mpfr_init2(right, prec);
        mpfi_get_left(left, val);
        mpfi_get_right(right, val);
        char* l_str = nullptr;
        char* r_str = nullptr;
        mpfr_asprintf(&l_str, "%.8Re", left);
        mpfr_asprintf(&r_str, "%.8Re", right);
        os << "[" << l_str << ", " << r_str << "]";
        mpfr_free_str(l_str);
        mpfr_free_str(r_str);
        mpfr_clear(left);
        mpfr_clear(right);
        return os;
    }

    std::unique_ptr<Field> dup() const final {
        return std::make_unique<FMpfi>(val, prec);
    }

    bool is_zero() const final {
        return mpfi_is_zero(val) != 0;
    }

    bool is_one() const final {
        return mpfi_cmp_si(val, 1) == 0;
    }

    bool parse(const std::string& str, const uint32_t line_no) final {
        uint32_t at = 0;
        ArjunNS::FMpq val_pre;
        if (!val_pre.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        mpfi_set_q(val, val_pre.get_val().get_mpq_t());
        return true;
    }

    void set_zero() final { mpfi_set_si(val, 0); }
    void set_one() final { mpfi_set_si(val, 1); }

    uint64_t bytes_used() const final {
        return sizeof(FMpfi) + mpfi_memory_usage(val);
    }
};

class FGenMpfi final : public CMSat::FieldGen {
public:
    uint16_t prec;
    ~FGenMpfi() final = default;
    FGenMpfi(const FGenMpfi& other) : prec(other.prec) {}
    FGenMpfi& operator=(const FGenMpfi& other) {
        if (this != &other) prec = other.prec;
        return *this;
    }
    explicit FGenMpfi(uint16_t _prec) : prec(_prec) {}

    std::unique_ptr<CMSat::Field> zero() const final {
        return std::make_unique<FMpfi>(0, prec);
    }

    std::unique_ptr<CMSat::Field> one() const final {
        return std::make_unique<FMpfi>(1, prec);
    }

    std::unique_ptr<FieldGen> dup() const final {
        return std::make_unique<FGenMpfi>(prec);
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
        const auto& ad = static_cast<const FMpfi&>(a);
        const auto& bd = static_cast<const FMpfi&>(b);
        // Compare using left endpoints of the intervals
        mpfr_t a_left, b_left;
        mpfr_init2(a_left, prec);
        mpfr_init2(b_left, prec);
        mpfi_get_left(a_left, ad.val);
        mpfi_get_left(b_left, bd.val);
        bool result = mpfr_cmp(a_left, b_left) > 0;
        mpfr_clear(a_left);
        mpfr_clear(b_left);
        return result;
    }

    bool weighted() const final { return true; }
};
