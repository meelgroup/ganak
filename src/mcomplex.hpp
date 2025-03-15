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

#include <cryptominisat5/solvertypesmini.h>
#include <gmpxx.h>
#include <memory>
#include <cstdlib>
#include <arjun/arjun.h>

class FComplex : public CMSat::Field {
public:
    mpq_class real;
    mpq_class imag;
    FComplex() : real(0), imag(0) {}
    FComplex(const mpq_class& _real, const mpq_class& _imag) : real(_real), imag(_imag) {}
    FComplex(const FComplex& other) : real(other.real), imag(other.imag) {}

    Field& operator=(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        real = od.real;
        imag = od.imag;
        return *this;
    }

    Field& operator+=(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        real += od.real;
        imag += od.imag;
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        return std::make_unique<FComplex>(real+od.real, imag+od.imag);
    }

    Field& operator-=(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        real -= od.real;
        imag -= od.imag;
        return *this;
    }

    Field& operator*=(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        mpq_class r = real;
        mpq_class i = imag;
        real = r*od.real-i*od.imag;
        imag = r*od.imag+i*od.real;
        return *this;
    }

    Field& operator/=(const Field& other) override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        if (od.is_zero()) throw std::runtime_error("Division by zero");
        mpq_class div = od.imag*od.imag+od.real*od.real;
        mpq_class r = real;
        mpq_class i = imag;
        real = r*od.real+i*od.imag;
        real /= div;
        imag = i*od.real-r*od.imag;
        imag /= div;
        return *this;
    }

    bool operator==(const Field& other) const override {
        const auto& od = dynamic_cast<const FComplex&>(other);
        return real == od.real && imag == od.imag;
    }

    std::ostream& display(std::ostream& os) const override {
        os << real << " + " << imag << "i";
        return os;
    }

    std::unique_ptr<Field> dup() const override {
        return std::make_unique<FComplex>(real, imag);
    }

    bool is_zero() const override {
        return real == 0 && imag == 0;
    }

    bool is_one() const override {
        return real == 1 && imag == 0;
    }

    void set_zero() override {
        real = 0;
        imag = 0;
    }

    void set_one() override {
        real = 1;
        imag = 0;
    }

    inline uint64_t helper(const mpz_class& v) const {
      return v.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    }

    uint64_t bytes_used() const override {
      return sizeof(FComplex) +
          helper(imag.get_num()) + helper(imag.get_den()) +
          helper(real.get_num()) + helper(real.get_den());
    }

    bool parse(const std::string& str, const uint32_t line_no) override {
        uint32_t at = 0;
        ArjunNS::FMpq _real;
        if (!_real.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        ArjunNS::FMpq _imag;
        if (!_imag.parse_mpq(str, at, line_no)) return false;
        skip_whitespace(str, at);
        real = _real.get_val();
        imag = _imag.get_val();
        return true;
    }
};

class FGenComplex : public CMSat::FieldGen {
public:
    ~FGenComplex() override = default;
    std::unique_ptr<CMSat::Field> zero() const override {
        return std::make_unique<FComplex>();
    }

    std::unique_ptr<CMSat::Field> one() const override {
        return std::make_unique<FComplex>(1, 0);
    }

    std::unique_ptr<FieldGen> dup() const override {
        return std::make_unique<FGenComplex>();
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const override {
      const auto& ad = dynamic_cast<const FComplex&>(a);
      const auto& bd = dynamic_cast<const FComplex&>(b);
      return ad.real > bd.real || (ad.real == bd.real && ad.imag > bd.imag);
    }

    bool weighted() const override { return true; }
};
