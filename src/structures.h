/*
 * structures.h
 *
 *  Created on: Jun 25, 2012
 *      Author: Marc Thurley
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <iostream>
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

  constexpr Lit() : value_(0) { }
  constexpr Lit(VariableIndex var, bool sign) : value_((var << 1) + (uint32_t) sign)
  {}

  VariableIndex var() const {
    return (value_ >> 1);
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

  void inc(){++value_;}

  void copyRaw(uint32_t v) {
    value_ = v;
  }

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
    if (lit.raw() == 0) {
        os << "UNDEF";
    } else {
        os << lit.toInt();
    }
    return os;
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
  vector<Lit> binary_links_ = vector<Lit>(1,SENTINEL_LIT);
  vector<ClOffsBlckL> watch_list_;
  uint32_t last_irred_bin = 0;
  double activity = 0.0;

  void removeWatchLinkTo(ClauseOfs clause_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (it->ofs == clause_ofs) {
        *it = watch_list_.back();
        watch_list_.pop_back();
        return;
      }
  }

  void replaceWatchLinkTo(ClauseOfs clause_ofs, ClauseOfs replace_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (it->ofs == clause_ofs) {
        it->ofs = replace_ofs;
        return;
      }
  }

  void addWatchLinkTo(ClauseIndex clause_ofs, Lit blockedLit) {
    watch_list_.push_back(ClOffsBlckL(clause_ofs, blockedLit));
  }

  void addBinLinkTo(Lit lit, bool irred) {
    binary_links_.back() = lit;
    binary_links_.push_back(SENTINEL_LIT);
    if (irred) last_irred_bin = binary_links_.size();
  }

  void resetWatchList(){
    watch_list_.clear();
  }

  bool hasBinaryLinkTo(Lit lit) {
    for (auto l : binary_links_) {
      if (l == lit)
        return true;
    }
    return false;
  }

  bool hasBinaryLinks() {
    return !binary_links_.empty();
  }
};

class Antecedent {
  uint32_t val_; // stuffed  value. LSB indicates long clause/binary cl

public:
  Antecedent() {
    val_ = 1;
  }

  Antecedent(const ClauseOfs cl_ofs) {
     val_ = (cl_ofs << 1) | 1;
   }
  Antecedent(const Lit idLit) {
    val_ = (idLit.raw() << 1);
  }

  bool isAClause() const {
    return val_ & 0x01;
  }

  ClauseOfs asCl() const {
    SLOW_DEBUG_DO(assert(isAClause()));
    return val_ >> 1;
  }

  Lit asLit() const {
    SLOW_DEBUG_DO(assert(!isAClause()));
    Lit idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }

  // Has an antecedent?
  bool isAnt() const {
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


struct Variable {
  Antecedent ante;
  int32_t decision_level = INVALID_DL;
  bool last_polarity = false;
  bool bprop = false;
  bool set_once = false; //it has once been set to some value
};

// for now Clause Header is just a dummy
// we keep it for possible later changes
class ClauseHeader {
  uint32_t creation_time_; // number of conflicts seen at creation time
  uint32_t score_;
  uint32_t length_;
public:

  void increaseScore() {
    score_++;
  }
  void decayScore() {
      score_ >>= 1;
  }
  uint32_t score() {
      return score_;
  }

  uint32_t creation_time() {
      return creation_time_;
  }
  uint32_t length(){ return length_;}
  void set_length(uint32_t length) {length_ = length;}

  void set_creation_time(uint32_t time) {
    creation_time_ = time;
  }
  static uint32_t overheadInLits() {return sizeof(ClauseHeader)/sizeof(Lit);}
};

#endif /* STRUCTURES_H_ */
