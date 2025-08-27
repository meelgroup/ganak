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
  explicit StaticSet(const vector<T>& values);
  explicit StaticSet(const vector<std::pair<T, T> >& values);
  void Init(const vector<T>& values);
  void Init(const vector<std::pair<T, T> >& values);
  int Rank(T value) const;
  T Kth(int k) const;
  int Size() const;
  const vector<T>& Values() const;
private:
  vector<T> values_;
};

using Edge = pair<int, int>;

class Graph {
public:
  explicit Graph(int n);
  explicit Graph(vector<Edge> edges);
  void addEdge(int v, int u);
  void addEdge(Edge e);

  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  vector<Edge> edges() const;

  Bitset Neighbors(const Bitset& vs) const;
  const vector<int>& Neighbors(int v) const;

  Graph(const Graph& rhs) = default;
  Graph& operator=(const Graph& rhs) = default;

private:
  int n_, m_;
  StaticSet<int> vertex_map_;
  vector<vector<int> > adj_list_;
  vector<Bitset> adj_mat2_;
};

class TreeDecomposition {
public:
  explicit TreeDecomposition(int nBags, int nVars);
  const vector<int>& neighbor_bags(int b) const;
  int numVars() const;
  int numBags() const;
  void addEdge(int a, int b);
  void setBag(int v, const vector<int>& bag);
  int Width() const;
  bool InBag(int b, int v) const;
  int getCentroid() const;
  void visualizeTree(const std::string& fname) const;
  vector<int> getOrd(int& centroid) const;
private:
  int nBags; // number of bags in the tree decomposition
  int nVars; // number of vertices in the original graph
  int width; // width of the tree decomposition
  Graph tree; // the tree of bags
  vector<vector<int>> bags;
  void OdDes(int b, int p, int d, vector<int>& ret, vector<int>& bagDepths) const;
  int CenDfs(int x, int p, int& cen) const;
  bool dfs(int x, int v, int p, vector<int>& u) const;
  bool bagsConnected(int start) const {
    if (nBags == 0) return true;
    vector<int> visited(nBags, 0);

    auto dfs = [&](const int b, auto&& dfs_ref) -> void {
      assert(b >= 0 && b < nBags);
      visited[b] = 1;
      for (int b2 : neighbor_bags(b)) {
        if (!visited[b2]) dfs_ref(b2, dfs_ref);
      }
    };
    dfs(start, dfs);

    for (int i = 0; i < nBags; i++) if (!visited[i]) return false;
    return true;
  }
};

template<typename T>
StaticSet<T>::StaticSet(const vector<T>& values) {
  Init(values);
}

template<typename T>
StaticSet<T>::StaticSet() : StaticSet(vector<T>()) { }

template<typename T>
StaticSet<T>::StaticSet(const vector<std::pair<T, T> >& values) {
  Init(values);
}

template<typename T>
void StaticSet<T>::Init(const vector<T>& values) {
  values_ = values;
  SortAndDedup(values_);
}

template<typename T>
void StaticSet<T>::Init(const vector<std::pair<T, T> >& values) {
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
const vector<T>& StaticSet<T>::Values() const {
  return values_;
}
} // namespace sspp
