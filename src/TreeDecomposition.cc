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

// TreeDecomposition
int TreeDecomposition::centroid(int npvars, int verb) {
	int centroid = -1;
	findCentroid(0, -1, centroid);

	int size = bags[centroid].size();
	int npvars_in_bags = 0;
	for(int i = 0; i < size; i++)
		if(bags[centroid][i] < npvars) npvars_in_bags++;

	cent = centroid;
	if (verb > 0) printf("c o centroid bag size %d  #npvars %d\n", size, npvars_in_bags);
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
