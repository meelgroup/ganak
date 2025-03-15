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

#include "arjun/arjun.h"
#include <cryptominisat5/solvertypesmini.h>
#include <gmpxx.h>
#include <memory>
#include <cstdlib>
#include <mpfr.h>

class MPFComplex : public CMSat::Field {
public:
    mpfr_t real;
    mpfr_t imag;
    MPFComplex() {
      mpfr_init2(real, 256);
      mpfr_init2(imag, 256);
      mpfr_set_si(real, 0, MPFR_RNDN);
      mpfr_set_si(imag, 0, MPFR_RNDN);
    }
    MPFComplex(int r, int i) {
      mpfr_init2(real, 256);
      mpfr_init2(imag, 256);
      mpfr_set_si(real, r, MPFR_RNDN);
      mpfr_set_si(imag, i, MPFR_RNDN);
    }
    MPFComplex(const mpfr_t& _real, const mpfr_t& _imag) {
      mpfr_init2(real, 256);
      mpfr_init2(imag, 256);
      mpfr_set(real, _real, MPFR_RNDN);
      mpfr_set(imag, _imag, MPFR_RNDN);
    }
    MPFComplex(const MPFComplex& other) : MPFComplex(other.real, other.imag) {}
    ~MPFComplex() override {
      mpfr_clear(real);
      mpfr_clear(imag);
    }

    Field& operator=(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        mpfr_set(real, od.real, MPFR_RNDN);
        mpfr_set(imag, od.imag, MPFR_RNDN);
        return *this;
    }

    Field& operator+=(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        mpfr_add(real, real, od.real, MPFR_RNDN);
        mpfr_add(imag, imag, od.imag, MPFR_RNDN);
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        mpfr_t r;
        mpfr_t i;
        mpfr_init2(r, 256);
        mpfr_init2(i, 256);
        mpfr_add(r, real, od.real, MPFR_RNDN);
        mpfr_add(i, imag, od.imag, MPFR_RNDN);
        auto ret = std::make_unique<MPFComplex>(r, i);
        mpfr_clear(r);
        mpfr_clear(i);
        return ret;
    }

    Field& operator-=(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        mpfr_sub(real, real, od.real, MPFR_RNDN);
        mpfr_sub(imag, imag, od.imag, MPFR_RNDN);
        return *this;
    }

    Field& operator*=(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        mpfr_t r;
        mpfr_init2(r, 256);

        mpfr_t tmp;
        mpfr_init2(tmp, 256);
        mpfr_t tmp2;
        mpfr_init2(tmp2, 256);

        mpfr_mul(tmp, real, od.real, MPFR_RNDN);
        mpfr_mul(tmp2, imag, od.imag, MPFR_RNDN);
        mpfr_sub(r, tmp, tmp2, MPFR_RNDN);

        mpfr_mul(tmp, real, od.imag, MPFR_RNDN);
        mpfr_mul(tmp2, imag, od.real, MPFR_RNDN);
        mpfr_add(imag, tmp, tmp2, MPFR_RNDN);

        mpfr_set(real, r, MPFR_RNDN);
        mpfr_clear(r);
        mpfr_clear(tmp);
        mpfr_clear(tmp2);

        /* mpq_class r = real; */
        /* mpq_class i = imag; */
        /* real = r*od.real-i*od.imag; */
        /* imag = r*od.imag+i*od.real; */
        return *this;
    }

    Field& operator/=(const Field& other) override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        if (od.is_zero()) throw std::runtime_error("Division by zero");
        mpfr_t div;
        mpfr_init2(div, 256);
        mpfr_t r;
        mpfr_init2(r, 256);

        mpfr_t tmp;
        mpfr_init2(tmp, 256);
        mpfr_t tmp2;
        mpfr_init2(tmp2, 256);

        // div = od.imag*od.imag+od.real*od.real;
        mpfr_mul(div, od.imag, od.imag, MPFR_RNDN);
        mpfr_mul(tmp, od.real, od.real, MPFR_RNDN);
        mpfr_add(div, div, tmp, MPFR_RNDN);

        mpfr_mul(tmp, real, od.real, MPFR_RNDN);
        mpfr_mul(tmp2, imag, od.imag, MPFR_RNDN);
        mpfr_add(r, tmp, tmp2, MPFR_RNDN);
        mpfr_div(r, r, div, MPFR_RNDN);

        mpfr_mul(tmp, imag, od.real, MPFR_RNDN);
        mpfr_mul(tmp2, real, od.imag, MPFR_RNDN);
        mpfr_sub(imag, tmp, tmp2, MPFR_RNDN);
        mpfr_div(imag, imag, div, MPFR_RNDN);
        mpfr_set(real, r, MPFR_RNDN);

        mpfr_clear(div);
        mpfr_clear(r);
        mpfr_clear(tmp);
        mpfr_clear(tmp2);

        /* mpq_class div = od.imag*od.imag+od.real*od.real; */
        /* mpq_class r = real; */
        /* mpq_class i = imag; */
        /* real = r*od.real+i*od.imag; */
        /* real /= div; */
        /* imag = i*od.real-r*od.imag; */
        /* imag /= div; */
        return *this;
    }

    bool operator==(const Field& other) const override {
        const auto& od = dynamic_cast<const MPFComplex&>(other);
        return mpfr_equal_p(real, od.real) && mpfr_equal_p(imag, od.imag);
    }

    std::ostream& display(std::ostream& os) const override {
      char* tmp = nullptr;
      mpfr_asprintf(&tmp, "%.6Rf + %.6Rfi", real, imag);
      os << tmp;
      mpfr_free_str(tmp);
      return os;
    }

    std::unique_ptr<Field> dup() const override {
        return std::make_unique<MPFComplex>(real, imag);
    }

    bool is_zero() const override {
        return mpfr_zero_p(real) && mpfr_zero_p(imag);
    }

    bool is_one() const override {
        return mpfr_cmp_si(real, 1) == 0 && mpfr_zero_p(imag);
    }

    void set_zero() override {
      mpfr_set_si(real, 0, MPFR_RNDN);
      mpfr_set_si(imag, 0, MPFR_RNDN);
    }

    void set_one() override {
      mpfr_set_si(real, 1, MPFR_RNDN);
      mpfr_set_si(imag, 0, MPFR_RNDN);
    }

    inline uint64_t helper(const mpz_class& v) const {
      return v.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    }

    unsigned int mpfr_memory_usage(const mpfr_t& x) const {
        // Get the precision of the mpfr_t variable (in bits)
        mpfr_prec_t precision = mpfr_get_prec(x);

        // Determine the size of a limb (in bits)
        size_t limb_size = sizeof(mp_limb_t) * 8;

        // Calculate the number of limbs required for the significand
        size_t num_limbs = (precision + limb_size - 1) / limb_size;

        // Calculate the memory used by the significand (in bytes)
        size_t significand_memory = num_limbs * sizeof(mp_limb_t);

        // Estimate memory used by the exponent and other fields
        // This is an approximation and may vary depending on the MPFR implementation
        size_t overhead_memory = sizeof(mpfr_exp_t) + sizeof(int) + sizeof(mpfr_prec_t);

        // Total memory usage
        size_t total_memory = significand_memory + overhead_memory;

        return (unsigned int)total_memory;
    }

    uint64_t bytes_used() const override {
      return sizeof(MPFComplex) + mpfr_memory_usage(real) + mpfr_memory_usage(imag);
    }

    bool parse(const std::string& str, const uint32_t line_no) override {
        uint32_t at = 0;
        ArjunNS::FMpq _real;
        if (!_real.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        ArjunNS::FMpq _imag;
        if (!_imag.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        mpfr_set_q(real, _real.get_val().get_mpq_t(), MPFR_RNDN);
        mpfr_set_q(imag, _imag.get_val().get_mpq_t(), MPFR_RNDN);
        return true;
    }
};

class FGenMPFComplex : public CMSat::FieldGen {
public:
    ~FGenMPFComplex() override = default;
    std::unique_ptr<CMSat::Field> zero() const override {
        return std::make_unique<MPFComplex>();
    }

    std::unique_ptr<CMSat::Field> one() const override {
        return std::make_unique<MPFComplex>(1, 0);
    }

    std::unique_ptr<FieldGen> dup() const override {
        return std::make_unique<FGenMPFComplex>();
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const override {
      const auto& ad = dynamic_cast<const MPFComplex&>(a);
      const auto& bd = dynamic_cast<const MPFComplex&>(b);
      return mpfr_greaterequal_p(ad.real, bd.real) && mpfr_greaterequal_p(ad.imag, bd.imag);
    }

    bool weighted() const override { return true; }
};
