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

typedef pair<int, int> Edge;

class Graph {
public:
  explicit Graph(int n);
  explicit Graph(std::vector<Edge> edges);
  void AddEdge(int v, int u);
  void AddEdge(Edge e);
  void AddEdges(const std::vector<Edge>& edges);

  void RemoveEdge(int v, int u);
  void RemoveEdgesBetween(int v, const std::vector<int>& vs);

  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  std::vector<Edge> Edges() const;
  std::vector<int> Vertices() const;

  int Degeneracy() const;
  Bitset Neighbors(const Bitset& vs) const;
  const std::vector<int>& Neighbors(int v) const;
  std::vector<std::vector<int> > Comps(const std::vector<int>& separator) const;
  std::vector<std::vector<int> > NComps(const std::vector<int>& separator) const;
  std::vector<Bitset> NComps(const Bitset& bs) const;
  std::vector<int> FindCompAndMark(int v, std::vector<char>& block) const;

  std::vector<Edge> EdgesIn(const std::vector<int>& vs) const;
  std::vector<Edge> FillEdges(const std::vector<int>& clq) const;
  std::vector<Edge> FillEdges(const Graph& other) const;
  std::vector<Edge> FillEdges(Bitset bs) const;
  void FillBS(Bitset bs);
  int FillSize(Bitset bs) const;

  int MapBack(int v) const;
  std::vector<int> MapBack(std::vector<int> vs) const;
  Edge MapBack(Edge e) const;
  std::vector<Edge> MapBack(std::vector<Edge> es) const;
  std::pair<int, int> MapBack(int v, int u) const;
  int MapInto(int v) const;
  std::vector<int> MapInto(std::vector<int> vs) const;

  void InheritMap(const Graph& parent);

  std::vector<std::vector<int>> CompNeighs(const std::vector<int>& block) const;
  std::vector<Bitset> CompNeighsBit(const Bitset& block) const;

  Bitset AnotherComp(int x, const Bitset& minsep) const;

  bool IsMinsep(const std::vector<int>& separator) const;
  bool IsMinsep(const Bitset& separator) const;
  bool HasNFullComps(const Bitset& separator, int n) const;

  void Dfs2(int v, Bitset& sep, Bitset& vis, std::vector<int>& f) const;
  std::vector<Bitset> BitComps(Bitset vis) const;
  void Dfs22(int v, Bitset& sep, Bitset& vis, std::vector<int>& f, const Bitset& good) const;
  void Dfs2Bit(Bitset& vis, Bitset& ne) const;

  std::vector<Bitset> adj_mat2_;

  std::vector<int> Distances(const std::vector<int>& start) const;
  std::vector<std::vector<int>> DistanceMatrix() const;

  int MaximalIS(const Bitset& vs) const;

  Graph(const Graph& rhs) = default;
  Graph& operator=(const Graph& rhs) = default;

private:
  int n_, m_;
  StaticSet<int> vertex_map_;
  std::vector<std::vector<int> > adj_list_;
  void Dfs(int v, std::vector<char>& blocked, std::vector<int>& component) const;
};

class TreeDecomposition {
 public:
   explicit TreeDecomposition(int bs_, int n_);
   const vector<vector<int>>& Bags() const;
   const vector<int>& Neighbors(int b) const;
   int nverts() const;
   int nbags() const;
   void AddEdge(int a, int b);
   void SetBag(int v, vector<int> bag);
   int Width() const;
   bool Verify(const Graph& graph) const;
   bool InBag(int b, int v) const;
   Graph Chordal() const;
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
