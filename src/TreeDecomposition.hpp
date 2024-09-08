/******************************************
Copyright (C) 2023 Kenji Hashimoto

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
#include <algorithm>

#include "bitset.hpp"
#include "structures.hpp"
namespace TWD {

class Graph {
public:
  Graph() : nodes(0), edges(0) { }
  Graph(int vars);

  void init(int n);
  void clear();

  int numNodes() { return this->nodes; }
  int numEdges() { return this->edges; }

  const vector<vector<int>>& get_adj_list() const;
  void addEdge(int v1, int v2);
  void contract(int v, int max_edges);
  bool hasEdge(int v1, int v2) { return adj_mat[v1].Get(v2); }
  const std::vector<int> Neighbors(int v) const { return adj_list[v]; }

protected:
  int nodes;
  int edges;
  std::vector<std::vector<int>> adj_list;
  std::vector<sspp::Bitset> adj_mat;
};

class TreeDecomposition : public Graph {
public:
  TreeDecomposition();

  void initBags() { bags.clear(); bags.resize(nodes); }
  std::vector<std::vector<int>>& Bags() { return bags; }
  bool inBag(int v, int x) const { return std::binary_search(bags[v].begin(), bags[v].end(), x); }

  void setWidth(int width) { tw = width; }
  int width() const { return tw; }

  void setNumGraphNodes(int n) { gnodes = n; }

  int centroid(int npvars, int verb = 0);
  std::vector<int> distanceFromCentroid(int npvars);
  double start_time;

private:
  int findCentroid(int v, int parent, int& centroid) const;
  void computeDistance(int v, int parent, int depth, std::vector<int>& distance);

  std::vector<std::vector<int>> bags;
  int tw;
  int gnodes;
  int cent;
};
}
