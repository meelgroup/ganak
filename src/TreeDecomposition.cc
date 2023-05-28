/*
 * TreeDecomposition.cc
 *
 *  Created on: 2022/03/29
 *      Author: k-hasimt
 */

#include "TreeDecomposition.h"

#include <iostream>
#include <assert.h>
#include "structures.h"
using namespace std;

// Graph
Graph::Graph(int vars)
{
	clear();
	init(vars);
}

void Graph::init(int n)
{
	nodes = n;

	adj_list.clear();
	adj_list.resize(nodes, {});

	adj_mat.clear();
	adj_mat.resize(n);
	for(int i=0; i<n; i++) {
		adj_mat[i] = sspp::Bitset(n);
		adj_mat[i].SetFalse(i);
	}
}
void Graph::clear()
{
	nodes = 0;
	edges = 0;
	adj_list.clear();
	adj_mat.clear();
}
void Graph::addEdge(int v1, int v2)
{
	if(adj_mat[v1].Get(v2))
		return;

	adj_list[v1].push_back(v2);
	adj_list[v2].push_back(v1);
	adj_mat[v1].SetTrue(v2);
	adj_mat[v2].SetTrue(v1);
	edges++;
}

const vector<vector<int>>& Graph::get_adj_list() const {
	return adj_list;
}

bool Graph::isClique(const vector<int>& adj)
{
	for(size_t i=0; i<adj.size(); i++)
		for(size_t j=i+1; j<adj.size(); j++)
			if(!hasEdge(adj[i], adj[j])) return false;

	return true;
}

void Graph::toDimacs(ostream& out, bool withHeader)
{
	if(withHeader)
		out << "p tw "<< nodes << " " << edges << endl;

	for(int i=0; i<nodes; i++) {
		for(auto j : adj_list[i]) {
			if(i>=j) continue;
			out << (i+1) << " " << (j+1) << endl;
		}
	}
}

// TreeDecomposition
int TreeDecomposition::centroid(int npvars) {
	int centroid = -1;
	findCentroid(0, -1, centroid);

	int size = bags[centroid].size();
	int npvars_in_bags = 0;
	for(int i = 0; i < size; i++)
		if(bags[centroid][i] < npvars) npvars_in_bags++;

	cent = centroid;
	printf("c o centroid bag size %d  #npvars %d\n", size, npvars_in_bags);
	return centroid;
}

int TreeDecomposition::findCentroid(int v, int parent, int &centroid) const
{
	int intros = 0;

	for (auto ch : Neighbors(v)) {
		if (ch == parent) continue;
		intros += findCentroid(ch, v, centroid);
		if (centroid != -1)
			return intros;
	}

	for (auto x : bags[v])
		if (parent == -1 || !inBag(parent, x)) {
			intros++;
		}

	if (intros >= gnodes / 2)
		centroid = v;

	return intros;
}
vector<int> TreeDecomposition::distanceFromCentroid(int npvars)
{
	int cen = cent == -1 ? centroid(npvars) : cent;
	vector<int> distance(gnodes, -1);

	computeDistance(cen, -1, 0, distance);
	return distance;
}
void TreeDecomposition::computeDistance(int v, int parent, int depth, vector<int>& distance)
{
	bool BagHasNewVar = false;

	for (auto x : bags[v]) {
		if (parent == -1 || distance[x] == -1) {
			distance[x] = depth;
			BagHasNewVar = true;
		}
	}

	int d = BagHasNewVar ? depth+1 : depth;

	for (auto ch : Neighbors(v)) {
		if (ch == parent) continue;
		computeDistance(ch, v, d, distance);
	}

}

void TreeDecomposition::toDimacs(ostream& out, int cent2, int npvars)
{
	out << "s td " << bags.size() << " " << (tw+1) << " " << gnodes << endl;
	out << "c pvars " << npvars << endl;
	if(cent2 >= 0) out << "c centroid " << cent2 << endl;

	int count = 1;
	for(const auto& bag : bags) {
		out << "b " << count ;
		for(auto v : bag) {
			out << " " << (v+1);
		}
		out << endl;
		count++;
	}

	Graph::toDimacs(out, false);
}
