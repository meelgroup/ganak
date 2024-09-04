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

#include "graph.hpp"
#include "utils.hpp"

namespace sspp {

Graph::Graph(int n)
  : n_(n), m_(0), adj_list_(n) {
  adj_mat2_.resize(n_);
  std::vector<int> identity(n);
  for (int i = 0; i < n; i++) {
    identity[i] = i;
    adj_mat2_[i] = Bitset(n);
    adj_mat2_[i].SetTrue(i);
  }
  vertex_map_.Init(identity);
}

Graph::Graph(std::vector<Edge> edges) : vertex_map_(edges) {
  n_ = vertex_map_.Size();
  m_ = 0;
  adj_list_.resize(n_);
  adj_mat2_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_mat2_[i] = Bitset(n_);
    adj_mat2_[i].SetTrue(i);
  }
  for (auto edge : edges) {
    AddEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
}

int Graph::n() const {
  return n_;
}

int Graph::m() const {
  return m_;
}

bool Graph::HasEdge(int v, int u) const {
  return adj_mat2_[v].Get(u);
}

bool Graph::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

std::vector<Edge> Graph::Edges() const {
  std::vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

std::vector<int> Graph::Vertices() const {
  std::vector<int> ret(n_);
  for (int i=0;i<n_;i++){
    ret[i] = i;
  }
  return ret;
}

const std::vector<int>& Graph::Neighbors(int v) const {
  return adj_list_[v];
}

Bitset Graph::Neighbors(const Bitset& vs) const {
  Bitset nbs(n_);
  for (int v : vs) {
    nbs |= adj_mat2_[v];
  }
  nbs.TurnOff(vs);
  return nbs;
}

void Graph::AddEdge(int v, int u) {
  if (HasEdge(v, u)) return;
  assert(v != u);
  m_++;
  adj_list_[v].push_back(u);
  adj_list_[u].push_back(v);
  adj_mat2_[v].SetTrue(u);
  adj_mat2_[u].SetTrue(v);
}

void Graph::AddEdge(Edge e) {
  AddEdge(e.first, e.second);
}

void Graph::AddEdges(const std::vector<Edge>& edges) {
  for (auto& edge : edges) AddEdge(edge);
}

std::vector<Edge> Graph::EdgesIn(const std::vector<int>& vs) const {
  std::vector<char> is(n_);
  for (int v : vs) {
    is[v] = true;
  }
  std::vector<Edge> edges;
  for (int v : vs) {
    if (adj_list_[v].size() <= vs.size()) { // Two cases for optimization
      for (int nv : adj_list_[v]) {
        if (is[nv] && nv > v) edges.push_back({v, nv});
      }
    }
    else {
      for (int nv : vs) {
        if (adj_mat2_[v].Get(nv) && nv > v) edges.push_back({v, nv});
      }
    }
  }
  return edges;
}

int Graph::FillSize(Bitset bs) const {
  int chunks = bs.chunks_;
  int ans = 0;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        ans += __builtin_popcountll(bs.data_[j] & (~adj_mat2_[v].data_[j]));
      }
    }
  }
  return ans;
}

void Graph::FillBS(Bitset bs) {
  int chunks = bs.chunks_;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        uint64_t td = bs.data_[j] & (~adj_mat2_[v].data_[j]);
        while (td) {
          int u = j*BITS + __builtin_ctzll(td);
          td &= ~-td;
          AddEdge(v, u);
        }
      }
    }
  }
}

int Graph::MapBack(int v) const {
  return vertex_map_.Kth(v);
}
std::vector<int> Graph::MapBack(std::vector<int> vs) const {
  for (int& v : vs) {
    v = MapBack(v);
  }
  return vs;
}
std::pair<int, int> Graph::MapBack(int v, int u) const {
  return {MapBack(v), MapBack(u)};
}
int Graph::MapInto(int v) const {
  return vertex_map_.Rank(v);
}
std::vector<int> Graph::MapInto(std::vector<int> vs) const {
  for (int& v : vs) {
    v = MapInto(v);
  }
  return vs;
}
Edge Graph::MapBack(Edge e) const {
  return {MapBack(e.first), MapBack(e.second)};
}
std::vector<Edge> Graph::MapBack(std::vector<Edge> es) const {
  for (Edge& e : es) {
    e = MapBack(e);
  }
  return es;
}

int Graph::MaximalIS(const Bitset& vs) const {
  Bitset is(n_);
  int ans = 0;
  for (int v : vs) {
    if (!is.Intersects(adj_mat2_[v])) {
      is.SetTrue(v);
      ans++;
    }
  }
  return ans;
}

TreeDecomposition::TreeDecomposition(int bs_, int n_)
 : bs(bs_), n(n_), width(-1), tree(bs+1), bags(bs+1) {}

void TreeDecomposition::AddEdge(int a, int b) {
  tree.AddEdge(a, b);
}

void TreeDecomposition::SetBag(int v, vector<int> bag) {
  assert(v >= 1 && v <= bs);
  assert(bags[v].empty());
  bags[v] = bag;
  SortAndDedup(bags[v]);
  width = std::max(width, (int)bags[v].size()-1);
  for (int u : bags[v]) {
    assert(0 <= u && u < n);
  }
}

int TreeDecomposition::Width() const {
  return width;
}

bool TreeDecomposition::InBag(int b, int v) const {
  assert(1 <= b && b <= bs && 0 <= v && v < n);
  return BS(bags[b], v);
}

const vector<vector<int>>& TreeDecomposition::Bags() const {
  return bags;
}

Graph TreeDecomposition::Chordal() const {
  Graph ret(n);
  for (const auto& bag : bags) {
    for (int i = 0; i < (int)bag.size(); i++) {
      for (int j = i+1; j < (int)bag.size(); j++) {
        assert(bag[i] >= 0 && bag[i] < n && bag[j] >= 0 && bag[j] < n);
        ret.AddEdge(bag[i], bag[j]);
      }
    }
  }
  return ret;
}

int TreeDecomposition::nbags() const {
  return bs;
}

int TreeDecomposition::nverts() const {
  return n;
}

const vector<int>& TreeDecomposition::Neighbors(int b) const {
  assert(b >= 1 && b <= bs);
  return tree.Neighbors(b);
}

int TreeDecomposition::CenDfs(int b, int p, int& cen) const {
  assert(b >= 1 && b <= bs);
  assert(p >= 0 && p <= bs);
  assert(cen == 0);
  int intro = 0;
  for (int nb : Neighbors(b)) {
    if (nb == p) continue;
    int cintro = CenDfs(nb, b, cen);
    intro += cintro;
    if (cintro >= n/2) {
      assert(cen);
      return intro;
    }
  }
  for (int v : bags[b]) {
    if (p == 0 || !InBag(p, v)) {
      intro++;
    }
  }
  if (intro >= n/2) {
    cen = b;
  }
  return intro;
}

int TreeDecomposition::Centroid() const {
  int cen = 0;
  CenDfs(1, 0, cen);
  assert(cen >= 1 && cen <= bs);
  return cen;
}

void TreeDecomposition::OdDes(int b, int p, int d, vector<int>& ret) const {
  assert(b >= 1 && b <= bs);
  assert(p >= 0 && p <= bs);
  assert(d >= 1);
  bool new_vs = false;
  for (int v : bags[b]) {
    if (ret[v] == 0) {
      new_vs = true;
    } else {
      assert(ret[v] <= d);
      assert(binary_search(bags[p].begin(), bags[p].end(), v));
    }
  }
  if (new_vs) {
    d++;
    for (int v : bags[b]) {
      if (ret[v] == 0) {
        ret[v] = d;
      }
    }
  }
  for (int nb : Neighbors(b)) {
    if (nb == p) continue;
    OdDes(nb, b, d, ret);
  }
}

vector<int> TreeDecomposition::GetOrd() const {
  int centroid = Centroid();
  assert(centroid >= 1 && centroid <= bs);
  vector<int> ret(n);
  OdDes(centroid, 0, 1, ret);
  for (int i = 0; i < n; i++) assert(ret[i] > 0);
  return ret;
}
} // namespace sspp
