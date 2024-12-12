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

#include <climits>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <cstdint>
#include "lit.hpp"
#include "common.hpp"
#include "Vec.hpp"

namespace GanakInt {

using ClauseIndex = uint32_t;
using ClauseOfs = uint32_t;
using CacheEntryID = uint32_t;
using TriValue = uint8_t;

constexpr uint32_t sentinel = 0;
constexpr int32_t INVALID_DL = -1;
constexpr uint8_t F_TRI = 0;
constexpr uint8_t T_TRI = 1;
constexpr uint8_t X_TRI = 2;
constexpr Lit NOT_A_LIT(0, false);
constexpr auto SENTINEL_LIT = NOT_A_LIT;

struct ClOffsBlckL {
  ClOffsBlckL() = default;
  ClOffsBlckL(ClauseOfs _ofs, Lit l) : ofs(_ofs), blckLit(l) {}
  ClauseOfs ofs;
  Lit blckLit;

  bool operator==(const ClOffsBlckL& other) const {
    return ofs == other.ofs && blckLit == other.blckLit;
  }

  bool operator!=(const ClOffsBlckL& other) const {
    return !(this->operator==(other));
  }
};

struct BinCl {
  BinCl(uint32_t val) = delete;
  BinCl(int val) = delete;
  explicit BinCl(Lit _lit, bool _red) {
    v = _lit.raw() << 1 | (uint32_t)_red;
  }
  Lit lit() const { return Lit::toLit(v >> 1); }
  bool red() const { return v&(1U); }
  bool irred() const { return !red(); }
  uint32_t v;

};

class LitWatchList {
public:
  vec<BinCl> binaries;
  vec<ClOffsBlckL> watch_list_;
  float activity = 0.0;

  void del_c(ClauseOfs offs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++) {
      if (it->ofs == offs) {
        *it = watch_list_.back();
        watch_list_.pop_back();
        return;
      }
    }
    assert(false && "should have found it!");
  }

  inline void add_cl(ClauseIndex offs, Lit blocked_lit) {
    watch_list_.push_back(ClOffsBlckL(offs, blocked_lit));
  }

  void add_bin(Lit lit, bool red) {
    binaries.push_back(BinCl(lit, red));
  }
};

enum class AnteType { clause, lit, decision };

class Antecedent {
  uint32_t val = 0;
  AnteType type = AnteType::decision;

public:
  Antecedent() {
    type = AnteType::decision;
  }

  explicit Antecedent(const ClauseOfs cl_ofs) {
     val = cl_ofs;
     type = AnteType::clause;
   }
  explicit Antecedent(const Lit idLit) {
    val = idLit.raw();
    type = AnteType::lit;
  }

  bool isAClause() const { return type == AnteType::clause; }
  bool isALit() const { return type == AnteType::lit; }
  bool isNull() const {return type == AnteType::decision;}
  bool isAnt() const {return !isNull();}

  ClauseOfs as_cl() const {
    SLOW_DEBUG_DO(assert(isAClause()));
    return val;
  }

  Lit as_lit() const {
    SLOW_DEBUG_DO(assert(isALit()));
    Lit idLit;
    idLit.copyRaw(val);
    return idLit;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Antecedent& val)
{
  std::stringstream s;
  if (val.isNull()) {
    s << std::setw(5) << "DEC  " << std::setw(10) << "";
  } else if (val.isAClause()) {
    s << "CL offs:  " << std::setw(10) << val.as_cl();
  } else if (val.isALit()) {
    s << "Lit: " << std::setw(10) << val.as_lit();
  } else {assert(false);}
  os << s.str();
  return os;
}

template<typename T>
struct Cube {
  Cube() = default;
  Cube(const std::vector<Lit>& _cnf, const T& _cnt, bool _symm = false) : cnf(_cnf), cnt(_cnt), symm(_symm) {}
  std::vector<Lit> cnf;
  T cnt;
  bool enabled = true;
  bool symm = false;
  uint32_t lbd = 100;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Cube<T>& c) {
  os << "CNF: " << c.cnf << " cnt: " << std::setprecision(40) << c.cnt << " enabled: " << (int)c.enabled << " symm: " << (int)c.symm << " lbd: " << c.lbd;
  return os;
}

struct VarData {
  Antecedent ante;
  int32_t decision_level = INVALID_DL;
  uint32_t sublevel;
  bool last_polarity = false;
};

class Clause {
public:
  Clause(bool _red, uint32_t _sz):
    sz(_sz), red(_red)  {}

  void set_used() {
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
  uint8_t vivifed:1 = 0;
  auto size() const { return sz; }
  void resize(const uint32_t sz2) {sz = sz2;}
  void update_lbd(uint32_t _lbd) {
    if (_lbd > 250) return;
    if (_lbd < lbd) lbd = _lbd;
  }
  Lit* data() const {
    return (Lit*)((char*)this + sizeof(Clause));
  }
  Lit* begin() { return data(); }
  Lit* end() { return begin()+sz; }
  Lit& operator[](uint32_t at) { return *(data()+at); }
  const Lit* begin() const { return data(); }
  const Lit* end() const { return begin()+sz; }
  const Lit& operator[](uint32_t at) const { return *(data()+at); }
};

inline std::ostream& operator<<(std::ostream& os, const Clause& cl) {
  for(const auto& l: cl) os << l << " ";
  os << "0"
    << " (red: " << (int)cl.red << " lbd: " << (int)cl.lbd
    << " used: " << (int)cl.used << " total_used: " << (int)cl.total_used
    << " vivifed: " << (int)cl.vivifed
    << " freed: " << (int)cl.freed;
  return os;
}

}
