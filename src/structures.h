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

  int toInt() const {
    return ((int) value_ >> 1) * ((sign()) ? 1 : -1);
  }

  uint32_t toPosInt() const {
    return value_;
  }

  void inc(){++value_;}

  void copyRaw(uint32_t v) {
    value_ = v;
  }

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
  else os << lit.toInt();
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
  vector<ClauseOfs> occ;
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
#ifdef SLOW_DEBUG
    found = false;
    for (auto& o: occ) {
      if (o == off) { o = replace_ofs; found = true; break; }
    }
    assert(found && "Should have found occ!!!");
#endif
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

class Antecedent {
  uint32_t val_; // stuffed  value. LSB indicates long clause/binary cl
  bool fake = false; // for PROBE and BPROP

public:
  Antecedent() {
    val_ = 1;
    fake = false;
  }

  explicit Antecedent(const ClauseOfs cl_ofs) {
     val_ = (cl_ofs << 1) | 1;
     fake = false;
   }
  explicit Antecedent(const Lit idLit, bool _fake = false) :
    fake (_fake) {
    if (_fake) return;
    val_ = (idLit.raw() << 1);
  }

  bool isDecision() const {return isAClause() && asCl() == NOT_A_CLAUSE;}
  bool isFake() const {return fake;}
  bool isAClause() const { return !fake && (val_ & 0x01); }
  ClauseOfs asCl() const {
    SLOW_DEBUG_DO(assert(isAClause()));
    return val_ >> 1;
  }

  Lit asLit() const {
    SLOW_DEBUG_DO(assert(!fake && !isAClause()));
    Lit idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }

  // Has an antecedent?
  bool isAnt() const {
    if (fake) return true;
    //Note that literals and clause offsets both start
    // at a higher-than 0 index, so if it's been set to be an antecdent, it'll be
    // different than 1
    return val_ != 1;
  }

  bool operator==(const Antecedent& other) const {
    return val_ == other.val_;
  }
  bool operator!=(const Antecedent& other) const {
    return val_ != other.val_;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Antecedent& val)
{
  std::stringstream s;
  if (val.isAClause() && val.asCl() == NOT_A_CLAUSE) {
    s << std::setw(5) << "DEC  " << std::setw(10) << "";
  } else if (!val.isAnt()) {
    s << std::setw(5) << "???  " << std::setw(10) << "";
  } else if (val.isFake()) {
    s << std::setw(5) <<"fake " << std::setw(10) << "";
  } else if (val.isAClause()) {
    s << "CL:  " << std::setw(10) << val.asCl();
  } else {
    s << "Lit: " << std::setw(10) << val.asLit();
  }
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

class ClHeader {
public:
  ClHeader(uint8_t _lbd, bool _red): lbd(_lbd), red(_red)  {}

  void increaseScore() {
    used = 1;
    total_used++;
  }
  uint32_t total_used = 0;
  uint8_t lbd;
  uint8_t used:1 = 1;
  uint8_t marked_deleted:1 = 0;
  uint8_t red:1 = 0;

  void update_lbd(uint32_t _lbd) {
    if (_lbd > 250) return;
    if (_lbd < lbd) lbd = _lbd;
  }

  constexpr static uint32_t overheadInLits() {
    return sizeof(ClHeader)/sizeof(Lit) + (bool)(sizeof(ClHeader)%sizeof(Lit));
  }
};
