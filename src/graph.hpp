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
  explicit Graph(int vars, const vector<vector<Lit>>& clauses);
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

  bool IsSimp(int v) const;

  int Degeneracy() const;
  
  bool IsConnected() const;
  bool IsConnectedOrIsolated() const;
  
  Bitset Neighbors(const Bitset& vs) const;
  const std::vector<int>& Neighbors(int v) const;
  std::vector<std::vector<int> > Components(const std::vector<int>& separator) const;
  std::vector<std::vector<int> > NComponents(const std::vector<int>& separator) const;
  std::vector<Bitset> NComponents(const Bitset& bs) const;
  std::vector<int> FindComponentAndMark(int v, std::vector<char>& block) const;
  
  bool IsClique(const std::vector<int>& clique) const;
  bool IsAlmostClique(const std::vector<int>& clq) const;
  bool IsClique(Bitset bs) const;

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
  bool HasNFullComponents(const Bitset& separator, int n) const;
  
  void Dfs2(int v, Bitset& sep, Bitset& vis, std::vector<int>& f) const;
  bool IsFull(int v, Bitset sep, Bitset vis) const;
  bool IsFull2(int v, Bitset sep, Bitset& vis) const;
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