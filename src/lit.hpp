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

#include <vector>
#include <set>
#include <cstdint>
#include <ostream>

namespace GanakInt {

class Lit {
public:
  Lit(uint32_t val) = delete;
  Lit(int val) = delete;
  explicit constexpr Lit() : value_(0) { }
  explicit constexpr Lit(uint32_t var, bool sign) : value_((var << 1) + (uint32_t) sign) {}

  uint32_t var() const { return (value_ >> 1U); }
  constexpr bool operator<(const Lit other) const { return value_ < other.value_; }
  constexpr bool operator>(const Lit other) const { return value_ > other.value_; }
  constexpr int to_visual_int() const { return ((int) value_ >> 1) * ((sign()) ? 1 : -1); }
  constexpr void inc() {++value_;}
  constexpr static Lit toLit(uint32_t data) {
    Lit l;
    l.copyRaw(data);
    return l;
  }
  constexpr void copyRaw(uint32_t v) { value_ = v; }

  // True if it's NON-NEGATED and False if it's NEGATED
  constexpr bool sign() const { return (bool) (value_ & 0x01); }
  constexpr bool operator!=(const Lit &rL2) const { return value_ != rL2.value_; }
  constexpr bool operator==(const Lit &rL2) const { return value_ == rL2.value_; }

  constexpr Lit neg() const {
    Lit l;
    l.value_ = value_ ^ 1;
    return l;
  }
  int val() const { return (sign() ? (int)var() : -1*(int)var()); }
  constexpr uint32_t raw() const { return value_;}

private:
  uint32_t value_;
};

inline std::ostream& operator<<(std::ostream& os, const Lit lit)
{
  if (lit.raw() == 0) os << "UNDEF";
  else os << lit.to_visual_int();
  return os;
}

inline std::ostream& operator<<(std::ostream& co, const std::vector<Lit>& lits)
{
  for (uint32_t i = 0; i < lits.size(); i++) {
    co << lits[i];
    if (i != lits.size()-1) co << " ";
  }
  return co;
}

inline std::ostream& operator<<(std::ostream& co, const std::set<Lit>& lits)
{
  size_t i = 0;
  for (const auto& l: lits) {
    co << l;
    if (i != lits.size()-1) co << " ";
    i++;
  }
  return co;
}

}
