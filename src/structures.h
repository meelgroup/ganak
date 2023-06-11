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
#include <cassert>
#include <iostream>
#include <iomanip>
#include "primitive_types.h"
#include "common.h"

using std::vector;

constexpr int32_t INVALID_DL = -1;

typedef uint8_t TriValue;
#define   F_TRI  0
#define   T_TRI  1
#define   X_TRI  2


class Lit {
public:
  Lit(uint32_t val) = delete;
  Lit(int val) = delete;
  explicit constexpr Lit() : value_(0) { }
  explicit constexpr Lit(VariableIndex var, bool sign) : value_((var << 1) + (uint32_t) sign)
  {}

  uint32_t var() const {
    return (value_ >> 1U);
  }

  constexpr bool operator<(const Lit other) const {
    return value_ < other.value_;
  }

  constexpr bool operator>(const Lit other) const {
    return value_ > other.value_;
  }

  // Does NOT do what you'd expect
  int toInteger() const {
    return ((int) value_ >> 1) * ((sign()) ? 1 : -1);
  }

  uint32_t toPosInt() const {
    return value_;
  }

  void inc(){++value_;}

  static Lit toLit(uint32_t data) {
    Lit l;
    l.copyRaw(data);
    return l;
  }

  void copyRaw(uint32_t v) { value_ = v; }

  // True if it's NON-NEGATED and False if it's NEGATED
  bool sign() const {
    return (bool) (value_ & 0x01);
  }

  bool operator!=(const Lit &rL2) const {
    return value_ != rL2.value_;
  }

  bool operator==(const Lit &rL2) const {
    return value_ == rL2.value_;
  }

  Lit neg() const {
    Lit l;
    l.value_ = value_ ^ 1;
    return l;
  }

  int val() const {
    return (sign() ? var() : -1*var());
  }
  uint32_t raw() const { return value_;}

private:
  uint32_t value_;
};

inline std::ostream& operator<<(std::ostream& os, const Lit lit)
{
  if (lit.raw() == 0) os << "UNDEF";
  else os << lit.toInteger();
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

static const Lit NOT_A_LIT(0, false);
#define SENTINEL_LIT NOT_A_LIT

struct ClOffsBlckL {
  ClOffsBlckL() {}
  ClOffsBlckL(ClauseOfs _ofs, Lit l) : ofs(_ofs), blckLit(l) {}
  ClauseOfs ofs;
  Lit blckLit;
};

class LitWatchList {
public:
  vector<Lit> binary_links_;
  vector<ClOffsBlckL> watch_list_;
  uint32_t last_irred_bin = 0;
  double activity = 0.0;

  void removeWatchLinkTo(ClauseOfs offs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (it->ofs == offs) {
        *it = watch_list_.back();
        watch_list_.pop_back();
        return;
      }
  }

  void replaceWatchLinkTo(ClauseOfs off, ClauseOfs replace_ofs) {
    bool found = false;
    for (auto& w: watch_list_) {
      if (w.ofs == off) { w.ofs = replace_ofs; found = true; break; }
    }
    assert(found && "Should have found watch!!!");
  }

  void addWatchLinkTo(ClauseIndex offs, Lit blockedLit) {
    watch_list_.push_back(ClOffsBlckL(offs, blockedLit));
  }

  void addBinLinkTo(Lit lit, bool red) {
    binary_links_.push_back(lit);
    if (!red) last_irred_bin = binary_links_.size();
  }

  void resetWatchList(){
    watch_list_.clear();
  }

  bool hasBinaryLinkTo(Lit lit) const {
    for (const auto& l : binary_links_) {
      if (l == lit) return true;
    }
    return false;
  }
};

enum class AnteType {
  clause, lit, fake, decision
};

class Antecedent {
  uint32_t val_;
  AnteType type = AnteType::decision;

public:
  Antecedent() {
    type = AnteType::decision;
  }

  explicit Antecedent(const ClauseOfs cl_ofs) {
     val_ = cl_ofs;
     type = AnteType::clause;
   }
  explicit Antecedent(const Lit idLit) {
    val_ = idLit.raw();
    type = AnteType::lit;
  }

  bool isNull() const {return type == AnteType::decision;}
  bool isAnt() const {return !isNull();}
  bool isFake() const {return type == AnteType::fake;}
  bool isAClause() const { return type == AnteType::clause; }
  bool isALit() const { return type == AnteType::lit; }

  ClauseOfs asCl() const {
    SLOW_DEBUG_DO(assert(isAClause()));
    return val_;
  }

  Lit asLit() const {
    SLOW_DEBUG_DO(assert(isALit()));
    Lit idLit;
    idLit.copyRaw(val_);
    return idLit;
  }

  bool operator==(const Antecedent& other) const {
    return val_ == other.val_;
  }
  bool operator!=(const Antecedent& other) const {
    return val_ != other.val_;
  }

  static Antecedent fakeAnte() {
    Antecedent ante;
    ante.type = AnteType::fake;
    return ante;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Antecedent& val)
{
  std::stringstream s;
  if (val.isNull()) {
    s << std::setw(5) << "DEC  " << std::setw(10) << "";
  } else if (val.isFake()) {
    s << std::setw(5) <<"fake " << std::setw(10) << "";
  } else if (val.isAClause()) {
    s << "CL:  " << std::setw(10) << val.asCl();
  } else if (val.isALit()) {
    s << "Lit: " << std::setw(10) << val.asLit();
  } else {assert(false);}
  os << s.str();
  return os;
}

struct Variable {
  Antecedent ante;
  int32_t decision_level = INVALID_DL;
  bool last_polarity = false;
  bool set_once = false; //it has once been set to some value
  uint32_t sublevel;
};

class Clause {
public:
  Clause(bool _red, uint32_t _sz):
    sz(_sz), red(_red)  {}

  void increaseScore() {
    used = 1;
    total_used++;
  }
  uint32_t sz;
  uint32_t total_used = 0;
  uint8_t lbd = 0;
  uint8_t used:1 = 1;
  uint8_t red:1 = 0;
  uint8_t freed:1 = 0;
  uint8_t reloced:1 = 0;

  void update_lbd(uint32_t _lbd) {
    if (_lbd > 250) return;
    if (_lbd < lbd) lbd = _lbd;
  }
  Lit* getData() const {
    return (Lit*)((char*)this + sizeof(Clause));
  }
  Lit* begin() {
    return getData();
  }
  const Lit* begin() const {
    return getData();
  }
  const Lit* end() const {
    return begin()+sz;
  }
  Lit* end() {
    return begin()+sz;
  }
  Lit& operator[](uint32_t at) {
    return *(getData()+at);
  }
  const Lit& operator[](uint32_t at) const {
    return *(getData()+at);
  }
};
