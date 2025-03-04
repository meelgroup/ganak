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
#include <cryptominisat5/solvertypesmini.h>
#include <memory>
#include <cstdlib>


class FParity : public CMSat::Field {
public:
    bool val;
    FParity() : val(0) {}
    FParity(const bool _val) : val(_val) {}
    FParity(const FParity& other) : val(other.val) {}

    Field& operator=(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        val = od.val;
        return *this;
    }

    Field& operator+=(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        val += od.val;
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        return std::make_unique<FParity>(val + od.val);
    }

    Field& operator-=(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        val -= od.val;
        return *this;
    }

    Field& operator*=(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        val &= od.val;
        return *this;
    }

    Field& operator/=(const Field& other) override {
        const auto& od = dynamic_cast<const FParity&>(other);
        if (od.val == 0) throw std::runtime_error("Division by zero");
        val /= od.val;
        return *this;
    }

    bool operator==(const Field& other) const override {
        const auto& od = dynamic_cast<const FParity&>(other);
        return val == od.val;
    }

    std::ostream& display(std::ostream& os) const override {
        os << val;
        return os;
    }

    std::unique_ptr<Field> dup() const override {
        return std::make_unique<FParity>(val);
    }

    bool is_zero() const override {
        return val == 0;
    }

    bool is_one() const override {
        return val == 1;
    }

    bool parse(const std::string& str, const uint32_t line_no) override {
        uint32_t at = 0;
        mpz_class head;
        if (!parse_int(head, str, at, line_no)) return false;
        val = head.get_ui();
        return check_end_of_weight(str, at, line_no);
    }

    void set_zero() override { val = 0; }
    void set_one() override { val = 1; }

    inline uint64_t helper(const mpz_class& v) const {
      return v.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    }

    uint64_t bytes_used() const override {
      return 1;
    }
};

class FGenParity : public CMSat::FieldGen {
public:
    ~FGenParity() override = default;
    std::unique_ptr<CMSat::Field> zero() const override {
        return std::make_unique<FParity>(0);
    }

    std::unique_ptr<CMSat::Field> one() const override {
        return std::make_unique<FParity>(1);
    }

    std::unique_ptr<FieldGen> dup() const override {
        return std::make_unique<FGenParity>();
    }

    bool weighted() const override { return true; }
};
