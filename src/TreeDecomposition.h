/*
 * TreeDecomposition.h
 *
 *  Created on: 2022/03/29
 *      Author: k-hasimt
 */

#ifndef PREPROCESSOR_TREEDECOMPOSITION_H_
#define PREPROCESSOR_TREEDECOMPOSITION_H_

#include <vector>
#include <algorithm>

#include "bitset.hpp"
#include "structures.h"

class Graph {
public:
  Graph() : nodes(0), edges(0) { }
  Graph(int vars);

  void init(int n);
  void clear();

  int numNodes() { return this->nodes; }
  int numEdges() { return this->edges; }

  void addEdge(int v1, int v2);
  bool hasEdge(int v1, int v2) { return adj_mat[v1].Get(v2); }
  const std::vector<int> Neighbors(int v) const { return adj_list[v]; }

  bool isSimplical(int v) { return isClique(adj_list[v]); }
  bool isClique(const std::vector<int>& adj);

  // Debug
  void toDimacs(std::ostream& out, bool withHeader=true);

protected:
  int nodes;
  int edges;
  std::vector<std::vector<int>> adj_list;
  std::vector<sspp::Bitset> adj_mat;
};

class TreeDecomposition : public Graph {
public:
  TreeDecomposition() : tw(0), gnodes(0), cent(-1) { }

  void initBags() { bags.clear(); bags.resize(nodes); }
  void setBag(int v, std::vector<int> bag) { bags[v] = bag; }
  std::vector<std::vector<int>>& Bags() { return bags; }
  bool inBag(int v, int x) const { return std::binary_search(bags[v].begin(), bags[v].end(), x); }

  void setWidth(int width) { tw = width; }
  int width() const { return tw; }

  void setNumGraphNodes(int n) { gnodes = n; }

  int centroid(int npvars);
  std::vector<int> distanceFromCentroid(int npvars);

  // Debug
  void toDimacs(std::ostream& out, int cent, int npvars);

private:
  int findCentroid(int v, int parent, int& centroid) const;
  void computeDistance(int v, int parent, int depth, std::vector<int>& distance);

  std::vector<std::vector<int>> bags;
  int tw;
  int gnodes;
  int cent;
};

#endif /* PREPROCESSOR_TREEDECOMPOSITION_H_ */
