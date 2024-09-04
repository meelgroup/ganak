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

#include "utils.hpp"
#include "bitset.hpp"

namespace sspp {

template<typename T>
class StaticSet {
public:
  StaticSet();
  explicit StaticSet(const std::vector<T>& values);
  explicit StaticSet(const std::vector<std::pair<T, T> >& values);
  void Init(const std::vector<T>& values);
  void Init(const std::vector<std::pair<T, T> >& values);
  int Rank(T value) const;
  T Kth(int k) const;
  int Size() const;
  std::vector<T> Values() const;
private:
  std::vector<T> values_;
};

using Edge = pair<int, int>;

class Graph {
public:
  explicit Graph(int n);
  explicit Graph(std::vector<Edge> edges);
  void AddEdge(int v, int u);
  void AddEdge(Edge e);

  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  std::vector<Edge> Edges() const;
  std::vector<int> Vertices() const;

  Bitset Neighbors(const Bitset& vs) const;
  const std::vector<int>& Neighbors(int v) const;

  Graph(const Graph& rhs) = default;
  Graph& operator=(const Graph& rhs) = default;

private:
  int n_, m_;
  StaticSet<int> vertex_map_;
  std::vector<std::vector<int> > adj_list_;
  std::vector<Bitset> adj_mat2_;
};

class TreeDecomposition {
 public:
   explicit TreeDecomposition(int bs_, int n_);
   const vector<int>& Neighbors(int b) const;
   int nverts() const;
   int nbags() const;
   void AddEdge(int a, int b);
   void SetBag(int v, vector<int> bag);
   int Width() const;
   bool Verify(const Graph& graph) const;
   bool InBag(int b, int v) const;
   int Centroid() const;
  vector<int> GetOrd() const;
 private:
   int bs, n, width;
   Graph tree;
   vector<vector<int>> bags;
  void OdDes(int b, int p, int d, vector<int>& ret) const;
   int CenDfs(int x, int p, int& cen) const;
   bool dfs(int x, int v, int p, vector<int>& u) const;
};

// Implementation
template<typename T>
StaticSet<T>::StaticSet(const std::vector<T>& values) {
  Init(values);
}

template<typename T>
StaticSet<T>::StaticSet() : StaticSet(std::vector<T>()) { }

template<typename T>
StaticSet<T>::StaticSet(const std::vector<std::pair<T, T> >& values) {
  Init(values);
}

template<typename T>
void StaticSet<T>::Init(const std::vector<T>& values) {
  values_ = values;
  SortAndDedup(values_);
}

template<typename T>
void StaticSet<T>::Init(const std::vector<std::pair<T, T> >& values) {
  values_.clear();
  for (const std::pair<T, T>& value : values) {
    values_.push_back(value.first);
    values_.push_back(value.second);
  }
  SortAndDedup(values_);
}

template<typename T>
int StaticSet<T>::Rank(T value) const {
  return std::lower_bound(values_.begin(), values_.end(), value) - values_.begin();
}

template<typename T>
T StaticSet<T>::Kth(int k) const {
  return values_[k];
}

template<typename T>
int StaticSet<T>::Size() const {
  return values_.size();
}

template<typename T>
std::vector<T> StaticSet<T>::Values() const {
  return values_;
}
} // namespace sspp
