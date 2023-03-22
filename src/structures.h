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

using std::vector;

constexpr int32_t INVALID_DL = -1;

typedef uint8_t TriValue;
#define   F_TRI  0
#define   T_TRI  1
#define   X_TRI  2

class Lit {
public:

  Lit() {
    value_ = 0;
  }
  Lit(int lit) {
    value_ = (abs(lit) << 1) + (unsigned) (lit > 0);
  }

  Lit(VariableIndex var, bool sign) {
    value_ = (var << 1) + (unsigned) sign;
  }

  VariableIndex var() const {
    return (value_ >> 1);
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

  const Lit neg() const {
    return Lit(var(), !sign());
  }

  int val() const {
    return (sign() ? var() : -1*var());
  }
  unsigned raw() const { return value_;}

private:
  unsigned value_;
};

inline std::ostream& operator<<(std::ostream& os, const Lit lit)
{
    if (lit == 0) {
        os << "UNDEF";
    } else {
        os << lit.toInt();
    }
    return os;
}

static const Lit NOT_A_LIT(0, false);
#define SENTINEL_LIT NOT_A_LIT

class LitWatchList {
public:
  vector<Lit> binary_links_ = vector<Lit>(1,SENTINEL_LIT);
  vector<ClauseOfs> watch_list_ = vector<ClauseOfs>(1,SENTINEL_CL);
  float activity_score_ = 0.0f;

  void increaseActivity(unsigned u = 1){
    activity_score_+= u;
  }

  void removeWatchLinkTo(ClauseOfs clause_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
          if (*it == clause_ofs) {
            *it = watch_list_.back();
            watch_list_.pop_back();
            return;
          }
  }

  void replaceWatchLinkTo(ClauseOfs clause_ofs, ClauseOfs replace_ofs) {
        for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
          if (*it == clause_ofs) {
            *it = replace_ofs;
            return;
          }
  }

  void addWatchLinkTo(ClauseIndex clause_ofs) {
    watch_list_.push_back(clause_ofs);
  }

  void addBinLinkTo(Lit lit) {
    binary_links_.back() = lit;
    binary_links_.push_back(SENTINEL_LIT);
  }

  void resetWatchList(){
        watch_list_.clear();
        watch_list_.push_back(SENTINEL_CL);
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
  uint32_t val_;

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
      return val_ >> 1;
    }

  Lit asLit() {
    Lit idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }
  // A NON-Antecedent will only be A NOT_A_CLAUSE Clause Id
  bool isAnt() {
    return val_ != 1; //i.e. NOT a NOT_A_CLAUSE;
  }

  bool operator==(const Antecedent& other) const {
    return val_ == other.val_;
  }
};


struct Variable {
  Antecedent ante;
  int32_t decision_level = INVALID_DL;
  bool polarity = false;
  bool set = false;
};

// for now Clause Header is just a dummy
// we keep it for possible later changes
class ClauseHeader {
  unsigned creation_time_; // number of conflicts seen at creation time
  unsigned score_;
  unsigned length_;
public:

  void increaseScore() {
    score_++;
  }
  void decayScore() {
      score_ >>= 1;
  }
  unsigned score() {
      return score_;
  }

  unsigned creation_time() {
      return creation_time_;
  }
  unsigned length(){ return length_;}
  void set_length(unsigned length) {length_ = length;}

  void set_creation_time(unsigned time) {
    creation_time_ = time;
  }
  static unsigned overheadInLits() {return sizeof(ClauseHeader)/sizeof(Lit);}
};

#endif /* STRUCTURES_H_ */
