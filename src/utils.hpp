/******************************************
Copyright (C) 2021 Tuukka Korhonen and Matti Jarvisalo

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
#include <string>
#include <iostream>
#include <cassert>
#include <limits>
#include <algorithm>
#include <chrono>
#include <random>

#include "bitset.hpp"

namespace sspp {

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::pair;

inline bool IsInt(const string& s, int64_t lb=std::numeric_limits<int64_t>::min(), int64_t ub=std::numeric_limits<int64_t>::max()) {
  try {
    int x = std::stoll(s);
    return lb <= x && x <= ub;
  } catch (...) {
    return false;
  }
}

inline bool IsDouble(const string& s, double lb=std::numeric_limits<double>::min(), double ub=std::numeric_limits<double>::max()) {
	try {
		double x = std::stod(s);
		return lb <= x && x <= ub;
	} catch (...) {
		return false;
	}
}

template<typename T>
void Shuffle(vector<T>& vec, std::mt19937& gen) {
	std::shuffle(vec.begin(), vec.end(), gen);
}

inline bool RandBool(std::mt19937& gen) {
	return std::uniform_int_distribution<int>(0,1)(gen);
}

template<typename T>
inline T RandInt(T a, T b, std::mt19937& gen) {
	return std::uniform_int_distribution<T>(a,b)(gen);
}

template<typename T>
T Power2(int p) {
	assert(p >= 0);
	if (p == 0) return 1;
	if (p%2 == 0) {
		T x = Power2<T>(p/2);
		return x*x;
	} else {
		return Power2<T>(p-1)*2;
	}
}

template<typename T>
bool IsSorted(const vector<T>& vec) {
	for (size_t i = 1; i < vec.size(); i++) {
		if (vec[i] < vec[i-1]) return false;
	}
	return true;
}

template<typename T>
void Append(std::vector<T>& a, const std::vector<T>& b) {
  a.reserve(a.size() + b.size());
  for (const T& x : b) {
    a.push_back(x);
  }
}

template<typename T>
std::vector<T> PermInverse(const std::vector<T>& perm) {
  std::vector<T> ret(perm.size());
  for (int i = 0; i < (int)perm.size(); i++) {
    ret[perm[i]] = i;
  }
  return ret;
}

inline Bitset ToBitset(const std::vector<int>& a, int n) {
  assert(n>=0);
  Bitset bs(n);
  for (int x : a) {
    assert(x>=0&&x<n);
    bs.SetTrue(x);
  }
  return bs;
}

template<typename T>
void SortAndDedup(vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

template<typename T>
void SwapDel(vector<T>& vec, size_t i) {
	assert(i < vec.size());
	std::swap(vec[i], vec.back());
	vec.pop_back();
}

template<typename T>
void ShiftDel(vector<T>& vec, size_t i) {
	assert(i < vec.size());
	for (; i+1 < vec.size(); i++) {
		vec[i] = std::move(vec[i+1]);
	}
	vec.pop_back();
}

template<typename T>
int Ind(const std::vector<T>& a, const T& x) {
	int ind = std::lower_bound(a.begin(), a.end(), x) - a.begin();
	assert(a[ind] == x);
	return ind;
}

template<typename T>
bool BS(const std::vector<T>& a, const T x) {
  return std::binary_search(a.begin(), a.end(), x);
}

class Timer {
 private:
  bool timing;
  std::chrono::duration<double> elapsedTime;
  std::chrono::time_point<std::chrono::steady_clock> startTime;
 public:
  Timer();
  void start();
  void stop();
  void clear();
  double get() const;
};

inline Timer::Timer() {
  timing = false;
  elapsedTime = std::chrono::duration<double>(std::chrono::duration_values<double>::zero());
}

inline void Timer::start() {
  if (timing) return;
  timing = true;
  startTime = std::chrono::steady_clock::now();
}

inline void Timer::stop() {
  if (!timing) return;
  timing = false;
  std::chrono::time_point<std::chrono::steady_clock> endTime = std::chrono::steady_clock::now();
  elapsedTime += (endTime - startTime);
}

inline double Timer::get() const {
  if (timing) {
  	auto tela = elapsedTime;
  	tela += (std::chrono::steady_clock::now() - startTime);
    return tela.count();
  }
  else {
    return elapsedTime.count();
  }
}

} // namespace sspp
