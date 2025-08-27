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
#include "common.hpp"
#include "utils.hpp"

namespace sspp {

Graph::Graph(int n)
  : n_(n), m_(0), adj_list_(n) {
  adj_mat2_.resize(n_);
  vector<int> identity(n);
  for (int i = 0; i < n; i++) {
    identity[i] = i;
    adj_mat2_[i] = Bitset(n);
    adj_mat2_[i].SetTrue(i);
  }
  vertex_map_.Init(identity);
}

Graph::Graph(vector<Edge> edges) : vertex_map_(edges) {
  n_ = vertex_map_.Size();
  m_ = 0;
  adj_list_.resize(n_);
  adj_mat2_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_mat2_[i] = Bitset(n_);
    adj_mat2_[i].SetTrue(i);
  }
  for (auto edge : edges) {
    addEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
}

int Graph::n() const { return n_; }
int Graph::m() const { return m_; }

bool Graph::HasEdge(int v, int u) const {
  return adj_mat2_[v].Get(u);
}

bool Graph::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

vector<Edge> Graph::edges() const {
  vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

const vector<int>& Graph::Neighbors(int v) const {
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

void Graph::addEdge(int v, int u) {
  if (HasEdge(v, u)) return;
  assert(v != u);
  m_++;
  adj_list_[v].push_back(u);
  adj_list_[u].push_back(v);
  adj_mat2_[v].SetTrue(u);
  adj_mat2_[u].SetTrue(v);
}

void Graph::addEdge(Edge e) {
  addEdge(e.first, e.second);
}

TreeDecomposition::TreeDecomposition(int nBags_, int nVars_)
 : nBags(nBags_), nVars(nVars_), width(-1), tree(nBags), bags(nBags) {}

void TreeDecomposition::addEdge(int a, int b) {
  tree.addEdge(a, b);
}

void TreeDecomposition::setBag(int b, const vector<int>& bag) {
  assert(b >= 0 && b < nBags);
  assert(bags[b].empty());
  bags[b] = bag;
  SortAndDedup(bags[b]);
  width = std::max(width, (int)bags[b].size());
  for (int v : bags[b]) assert(0 <= v && v < nVars);
}

int TreeDecomposition::Width() const {
  return width;
}

bool TreeDecomposition::InBag(int b, int v) const {
  assert(0 <= b && b < nBags);
  assert(0 <= v && v < nVars);
  return binary_search(bags[b], v);
}

int TreeDecomposition::numBags() const { return nBags; }
int TreeDecomposition::numVars() const { return nVars; }

const vector<int>& TreeDecomposition::neighbor_bags(int b) const {
  assert(b >= 0 && b < nBags);
  return tree.Neighbors(b);
}

int TreeDecomposition::CenDfs(int bag, int parent, int& cen) const {
  assert(bag >= 0 && bag < nBags);
  assert(cen == -1);
  int intro = 0;
  for (int nb : neighbor_bags(bag)) {
    if (nb == parent) continue;
    int cintro = CenDfs(nb, bag, cen);
    intro += cintro;
    if (cintro >= nVars/2) {
      assert(cen >= 0);
      return intro;
    }
  }
  for (int v : bags[bag]) {
    if (parent == -1 || !InBag(parent, v)) intro++;
  }
  if (intro >= nVars/2) cen = bag;
  return intro;
}

int TreeDecomposition::getCentroid() const {
  int cen = -1;
  CenDfs(0, -1, cen);
  assert(cen >= 0 && cen < nBags);
  return cen;
}

/**
    bag: Current bag in the tree decomposition.
    pparent: Parent bag.
    depth: Current depth
    ret: Output vector storing the order of each variable.
*/
void TreeDecomposition::OdDes(int bag, int parent, int depth, vector<int>& ret, vector<int>& bagDepths) const {
  VERBOSE_DEBUG_DO(cout << "c o OdDes: bag=" << bag << ", parent=" << parent << ", depth=" << depth << endl);
  assert(bag >= 0 && bag < nBags);
  assert(depth >= 1);
  if (bagDepths[bag] != -1 && bagDepths[bag] <= depth) return;

  depth++;
  bagDepths[bag] = depth;
  for (int v : bags[bag]) if (ret[v] > depth || ret[v] == -1) ret[v] = depth;
  for (int b : neighbor_bags(bag)) {
    if (b == parent) continue;
    OdDes(b, bag, depth, ret, bagDepths);
  }
}

// Gets the order of vertices in the tree decomposition
//    Assigns an incremental order to vertices in the graph based on their
//    appearance in the tree decomposition. The order starts from the centroid
//    and propagates outward, ensuring: Vertices in parent bags are processed
//    before their children. Newly discovered vertices get a higher (later)
//    order.
vector<int> TreeDecomposition::getOrd() const {
  int centroid = getCentroid();
  SLOW_DEBUG_DO(assert(bagsConnected(centroid)));
  VERBOSE_DEBUG_DO(cout << "c o Tree decomposition centroid: " << centroid << endl);
  assert(centroid >= 0 && centroid < nBags);
  vector<int> ret(nVars, -1);
  vector<int> bagDepths(nBags, -1);
  OdDes(centroid, -1, 1, ret, bagDepths);

  for (int i = 0; i < nVars; i++) if (ret[i] == 0) {
    cout << "c o Vertex " << i << " not assigned an order!" << endl;
  }
  for (int i = 0; i < nVars; i++) assert(ret[i] > 0);
  return ret;
}

}
