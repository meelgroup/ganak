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
#include <arjun/arjun.h>
#include <memory>

// Integers mod a prime p < 2^31. val is always in [0, p), so a product fits
// in 64 bits and == compares elements, not representations.
class FPrime final : public CMSat::Field {
public:
    int64_t val;
    int64_t field;
    FPrime(const int64_t _val, const int64_t _field) : val(norm(_val, _field)), field(_field) {}
    FPrime(const FPrime& other) : val(other.val), field(other.field) {}

    static int64_t norm(int64_t v, const int64_t p) {
        v %= p;
        return v < 0 ? v + p : v;
    }

    // a^-1 mod p, by the extended Euclidean algorithm
    static int64_t inverse(const int64_t a, const int64_t p) {
        int64_t t = 0, new_t = 1, r = p, new_r = a;
        while (new_r != 0) {
            const int64_t q = r / new_r;
            t -= q*new_t; std::swap(t, new_t);
            r -= q*new_r; std::swap(r, new_r);
        }
        if (r != 1) throw std::runtime_error("Not invertible mod the prime");
        return norm(t, p);
    }

    Field& operator=(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        val = od.val;
        return *this;
    }

    Field& operator+=(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        val += od.val;
        if (val >= field) val -= field;
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        return std::make_unique<FPrime>(val + od.val, field);
    }

    Field& operator-=(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        val -= od.val;
        if (val < 0) val += field;
        return *this;
    }

    Field& operator*=(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        val = (val * od.val) % field;
        return *this;
    }

    Field& operator/=(const Field& other) final {
        const auto& od = static_cast<const FPrime&>(other);
        if (od.val == 0) throw std::runtime_error("Division by zero");
        val = (val * inverse(od.val, field)) % field;
        return *this;
    }

    bool operator==(const Field& other) const final {
        const auto& od = static_cast<const FPrime&>(other);
        return val == od.val;
    }

    std::ostream& display(std::ostream& os) const final {
        os << val << " mod " << field;
        return os;
    }

    std::unique_ptr<Field> dup() const final {
        return std::make_unique<FPrime>(val, field);
    }

    bool is_zero() const final {
        return val == 0;
    }

    bool is_one() const final {
        return val == 1;
    }

    // Integers, decimals and fractions, as mode 1 reads them: a/b is a*b^-1 mod p
    bool parse(const std::string& str, const uint32_t line_no) final {
        uint32_t at = 0;
        ArjunNS::FMpq q;
        if (!q.parse_mpq(str, at, line_no)) return false;
        // floor remainders: non-negative for negative weights too
        const auto p = (unsigned long)field;
        const int64_t num = (int64_t)mpz_fdiv_ui(q.get_val().get_num_mpz_t(), p);
        const int64_t den = (int64_t)mpz_fdiv_ui(q.get_val().get_den_mpz_t(), p);
        if (den == 0) {
            std::cerr << "PARSE ERROR! The weight's denominator is divisible by the prime "
                << field << " at line " << line_no << std::endl;
            return false;
        }
        val = (num * inverse(den, field)) % field;
        return check_end_of_weight(str, at, line_no);
    }

    void set_zero() final { val = 0; }
    void set_one() final { val = 1; }

    uint64_t bytes_used() const final {
      return ArjunNS::heap_bytes(sizeof(FPrime));
    }
};

class FGenPrime final : public CMSat::FieldGen {
public:
    int64_t field;
    FGenPrime(int64_t _field) : field(_field) {}
    ~FGenPrime() final = default;
    std::unique_ptr<CMSat::Field> zero() const final {
        return std::make_unique<FPrime>(0, field);
    }

    std::unique_ptr<CMSat::Field> one() const final {
        return std::make_unique<FPrime>(1, field);
    }

    std::unique_ptr<FieldGen> dup() const final {
        return std::make_unique<FGenPrime>(field);
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
      const auto& ad = static_cast<const FPrime&>(a);
      const auto& bd = static_cast<const FPrime&>(b);
      return ad.val > bd.val;
    }

    bool weighted() const final { return true; }
    bool exact() const final { return true; }
};

class FParity final : public CMSat::Field {
public:
    bool val;
    FParity() : val(false) {}
    FParity(const bool _val) : val(_val) {}
    FParity(const FParity& other) : val(other.val) {}

    Field& operator=(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        val = od.val;
        return *this;
    }

    Field& operator+=(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        val ^= od.val;
        return *this;
    }

    std::unique_ptr<Field> add(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        return std::make_unique<FParity>(val ^ od.val);
    }

    Field& operator-=(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        val ^= od.val;
        return *this;
    }

    Field& operator*=(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        val &= od.val;
        return *this;
    }

    Field& operator/=(const Field& other) final {
        const auto& od = static_cast<const FParity&>(other);
        if (od.val == 0) throw std::runtime_error("Division by zero");
        return *this;
    }

    bool operator==(const Field& other) const final {
        const auto& od = static_cast<const FParity&>(other);
        return val == od.val;
    }

    std::ostream& display(std::ostream& os) const final {
        os << val;
        return os;
    }

    std::unique_ptr<Field> dup() const final {
        return std::make_unique<FParity>(val);
    }

    bool is_zero() const final {
        return val == false;
    }

    bool is_one() const final {
        return val == true;
    }

    bool parse(const std::string& str, const uint32_t line_no) final {
        uint32_t at = 0;
        mpz_class head;
        if (!parse_int(head, str, at, line_no)) return false;
        val = (head.get_ui() % 2) == 1;
        return check_end_of_weight(str, at, line_no);
    }

    void set_zero() final { val = false; }
    void set_one() final { val = true; }

    uint64_t bytes_used() const final {
      return ArjunNS::heap_bytes(sizeof(FParity));
    }
};

class FGenParity final : public CMSat::FieldGen {
public:
    ~FGenParity() final = default;
    std::unique_ptr<CMSat::Field> zero() const final {
        return std::make_unique<FParity>(false);
    }

    std::unique_ptr<CMSat::Field> one() const final {
        return std::make_unique<FParity>(true);
    }

    std::unique_ptr<FieldGen> dup() const final {
        return std::make_unique<FGenParity>();
    }

    bool larger_than(const CMSat::Field& a, const CMSat::Field& b) const final {
        const auto& ad = static_cast<const FParity&>(a);
        const auto& bd = static_cast<const FParity&>(b);
        return ad.val > bd.val;
    }

    bool weighted() const final { return true; }
    bool exact() const final { return true; }
};
