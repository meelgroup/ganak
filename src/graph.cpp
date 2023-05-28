#include "graph.hpp"

#include "utils.hpp"

#include <queue>

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

Graph::Graph(int vars, const vector<vector<Lit>>& clauses) : Graph((int)vars+1) {
	for (const auto& clause : clauses) {
		for (int i = 0; i < (int)clause.size(); i++) {
			for (int j = i+1; j < (int)clause.size(); j++) {
				Var v1 = VarOf(clause[i]);
				Var v2 = VarOf(clause[j]);
				assert(v1 >= 1 && v1 <= vars && v2 >= 1 && v2 <= vars);
				if (v1 != v2) {
					AddEdge(v1, v2);
				}
			}
		}
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

bool Graph::IsConnected() const {
  auto cs = Components({});
  return (cs.size() == 1) && ((int)cs[0].size() == n_);
}

bool Graph::IsConnectedOrIsolated() const {
  auto cs = Components({});
  int f = 0;
  for (const auto& c : cs) {
    if ((int)c.size() > 1) f++;
  }
  return f <= 1;
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

void Graph::RemoveEdge(int v, int u) {
  assert(HasEdge(v, u) && HasEdge(u, v));
  m_--;
  adj_mat2_[v].SetFalse(u);
  adj_mat2_[u].SetFalse(v);
  int fo = 0;
  for (int i = 0; i < (int)adj_list_[v].size(); i++) {
    if (adj_list_[v][i] == u) {
      std::swap(adj_list_[v][i], adj_list_[v].back());
      adj_list_[v].pop_back();
      fo++;
      break;
    }
  }
  for (int i = 0; i < (int)adj_list_[u].size(); i++) {
    if (adj_list_[u][i] == v) {
      std::swap(adj_list_[u][i], adj_list_[u].back());
      adj_list_[u].pop_back();
      fo++;
      break;
    }
  }
  assert(fo == 2);
}

int Graph::Degeneracy() const {
  std::vector<std::vector<int>> q(n_);
  std::vector<int> dg(n_);
  int vs = 0;
  for (int i=0;i<n_;i++) {
    dg[i] = adj_list_[i].size();
    if (dg[i] > 0) {
      q[dg[i]].push_back(i);
      vs++;
    }
  }
  int mt = 0;
  int t = 0;
  for (int it=0;it<vs;it++) {
    int x = -1;
    while (x == -1) {
      if (q[t].empty()) {
        t++;
      } else {
        x = q[t].back();
        q[t].pop_back();
        if (dg[x] == -1) {
          x = -1;
        } else {
          assert(dg[x] == t);
        }
      }
    }
    mt = std::max(mt, t);
    assert(x>=0&&x<n_&&dg[x]==t&&t>=0);
    for (int nx : adj_list_[x]) {
      if (dg[nx] >= 0) {
        assert(dg[nx] >= 1);
        dg[nx]--;
        dg[x]--;
        q[dg[nx]].push_back(nx);
      }
    }
    assert(dg[x] == 0);
    dg[x] = -1;
    t = std::max(0, t-1);
  }
  return mt;
}

void Graph::Dfs(int v, std::vector<char>& block, std::vector<int>& component) const {
  block[v] = true;
  component.push_back(v);
  for (int nv : adj_list_[v]) {
    if (!block[nv]) {
      Dfs(nv, block, component);
    }
  }
}

std::vector<int> Graph::FindComponentAndMark(int v, std::vector<char>& block) const {
  std::vector<int> component;
  Dfs(v, block, component);
  return component;
}

std::vector<std::vector<int> > Graph::Components(const std::vector<int>& separator) const {
  std::vector<char> blocked(n_);
  for (int v : separator) {
    blocked[v] = true;
  }
  std::vector<std::vector<int> > components;
  for (int i = 0; i < n_; i++) {
    if (!blocked[i]) {
      components.push_back(FindComponentAndMark(i, blocked));
    }
  }
  return components;
}

std::vector<std::vector<int> > Graph::NComponents(const std::vector<int>& separator) const {
  std::vector<char> blocked(n_);
  for (int v : separator) {
    blocked[v] = true;
  }
  std::vector<std::vector<int> > components;
  for (int v : separator) {
    for (int nv : adj_list_[v]) {
      if (!blocked[nv]) {
        components.push_back(FindComponentAndMark(nv, blocked));
      }
    }
  }
  return components;
}

bool Graph::IsClique(const std::vector<int>& clique) const {
  for (int i = 0; i < (int)clique.size(); i++) {
    for (int ii = i + 1; ii < (int)clique.size(); ii++) {
      if (!HasEdge(clique[i], clique[ii])) return false;
    }
  }
  return true;
}

bool Graph::IsFull(int v, Bitset sep, Bitset vis) const {
  static std::vector<int> q;
  q.resize(n_);
  q[0] = v;
  vis.SetFalse(v);
  int i = 0;
  int s = 1;
  int chunks = vis.Chunks();
  while (i < s) {
    int x = q[i++];
    for (int j = 0; j < chunks; j++) {
      uint64_t go = vis.data_[j] & adj_mat2_[x].data_[j];
      while (go) {
        vis.data_[j] &= (~(go&-go));
        q[s++] = __builtin_ctzll(go) + j*BITS;
        go &= ~-go;
      }
    }
    bool has = false;
    for (int j = 0; j < chunks; j++) {
      sep.data_[j] &= (~adj_mat2_[x].data_[j]);
      if (sep.data_[j]) has = true;
    }
    if (!has) return true;
  }
  return false;
}

bool Graph::IsFull2(int v, Bitset sep, Bitset& vis) const {
  static std::vector<int> q;
  q.resize(n_);
  q[0] = v;
  vis.SetFalse(v);
  int i = 0;
  int s = 1;
  int chunks = vis.Chunks();
  while (i < s) {
    int x = q[i++];
    for (int j = 0; j < chunks; j++) {
      uint64_t go = vis.data_[j] & adj_mat2_[x].data_[j];
      while (go) {
        vis.data_[j] &= (~(go&-go));
        q[s++] = __builtin_ctzll(go) + j*BITS;
        go &= ~-go;
      }
    }
    bool has = false;
    for (int j = 0; j < chunks; j++) {
      sep.data_[j] &= (~adj_mat2_[x].data_[j]);
      if (sep.data_[j]) has = true;
    }
    if (!has) return true;
  }
  return false;
}

void Graph::Dfs22(int v, Bitset& sep, Bitset& vis, std::vector<int>& f, const Bitset& good) const {
  vis.SetTrue(v);
  int chunks = vis.Chunks();
  Bitset ne(n_);
  ne.SetTrue(v);
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        gv &= ~-gv;
      }
      uint64_t gs = sep.data_[j] & ne.data_[j];
      while (gs) {
        sep.data_[j] &= (~(gs&-gs));
        if (good.data_[j] & (gs&-gs)) {
          f.push_back(__builtin_ctzll(gs) + j*BITS);
        }
        gs &= ~-gs;
      }
    }
  }
}

void Graph::Dfs2(int v, Bitset& sep, Bitset& vis, std::vector<int>& f) const {
  vis.SetTrue(v);
  int chunks = vis.Chunks();
  Bitset ne(n_);
  ne.SetTrue(v);
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        gv &= ~-gv;
      }
      uint64_t gs = sep.data_[j] & ne.data_[j];
      while (gs) {
        sep.data_[j] &= (~(gs&-gs));
        f.push_back(__builtin_ctzll(gs) + j*BITS);
        gs &= ~-gs;
      }
    }
  }
}

std::vector<Bitset> Graph::NComponents(const Bitset& separator) const {
  Bitset vis(n_);
  vis.FillUpTo(n_);
  vis.TurnOff(separator);
  Bitset nbs = Neighbors(separator);
  std::vector<Bitset> ret;
  for (Bitset comp : BitComps(vis)) {
    if (comp.Intersects(nbs)) {
      ret.push_back(comp);
    }
  }
  return ret;
}

std::vector<Bitset> Graph::BitComps(Bitset vis) const {
  Bitset ne(n_);
  int chunks = vis.Chunks();
  std::vector<Bitset> ret;
  bool fo = false;
  while (1) {
    if (!fo) {
      for (int j = 0; j < chunks; j++) {
        if (vis.data_[j]) {
          int x = __builtin_ctzll(vis.data_[j]) + j*BITS;
          ne.SetTrue(x);
          ret.push_back(Bitset(n_));
          fo = true;
          break;
        }
      }
      if (!fo) return ret;
    }
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        ret.back().SetTrue(x);
        gv &= ~-gv;
      }
    }
  }
}

std::vector<std::vector<int>> Graph::CompNeighs(const std::vector<int>& block) const {
  Bitset sep(n_);
  Bitset vis(n_);
  for (int i = 0; i < n_; i++) {
    if (adj_list_[i].size() > 0) {
      vis.SetTrue(i);
    }
  }
  for (int x : block) {
    sep.SetTrue(x);
    vis.SetFalse(x);
  }
  Bitset sb = sep;
  std::vector<std::vector<int>> ret;
  std::vector<int> f;
  f.reserve(block.size());
  for (int i = 0; i < n_; i++) {
    if (vis.Get(i)) {
      Dfs2(i, sep, vis, f);
      if (!f.empty()) {
        ret.push_back(f);
        f.clear();
      }
      sep = sb;
    }
  }
  return ret;
}

bool Graph::HasNFullComponents(const Bitset& separator, int n) const {
  Bitset vis(n_);
  vis.FillTrue();
  vis.TurnOff(separator);
  std::vector<int> f;
  f.reserve(separator.Popcount());
  int cnt = 0;
  Bitset ne(n_);
  for (int i = 0; i < n_; i++) {
    if (vis.Get(i) && adj_list_[i].size() > 0) {
      ne.SetTrue(i);
      Dfs2Bit(vis, ne);
      if (ne.Subsumes(separator)) {
        cnt++;
      }
      if (cnt >= n) return true;
      ne.Clear();
    }
  }
  return false;
}

bool Graph::IsMinsep(const Bitset& separator) const {
  return HasNFullComponents(separator, 2);
}

bool Graph::IsMinsep(const std::vector<int>& separator) const {
  return HasNFullComponents(ToBitset(separator, n_), 2);
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

Bitset Graph::AnotherComp(int x, const Bitset& minsep) const {
  // could be optimized
  int fulls = 0;
  bool ffx = false;
  bool fnx = false;
  Bitset ret;
  Bitset vis(n_);
  vis.FillUpTo(n_);
  vis.TurnOff(minsep);
  for (auto comp : BitComps(vis)) {
    if (Neighbors(comp).Popcount() == minsep.Popcount()) {
      fulls++;
      if (!comp.Get(x)) {
        fnx = true;
        ret = comp;
      } else {
        ffx = true;
      }
    }
  }
  assert(fulls == 2);
  assert(ffx && fnx);
  return ret;
}

std::vector<Edge> Graph::FillEdges(const std::vector<int>& clq) const {
  std::vector<Edge> ret;
  for (int i=0;i<(int)clq.size();i++){
    for (int ii=i+1;ii<(int)clq.size();ii++){
      if (!HasEdge(clq[i], clq[ii])) {
        ret.push_back(std::minmax(clq[i], clq[ii]));
      }
    }
  }
  return ret;
}

std::vector<Edge> Graph::FillEdges(const Graph& other) const {
  std::vector<Edge> ret;
  for (auto e : other.Edges()) {
    if (!HasEdge(e)) {
      ret.push_back(e);
    }
  }
  return ret;
}

std::vector<Edge> Graph::FillEdges(Bitset bs) const {
  int chunks = bs.chunks_;
  std::vector<Edge> ret;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        uint64_t td = bs.data_[j] & (~adj_mat2_[v].data_[j]);
        while (td) {
          int u = j*BITS + __builtin_ctzll(td);
          td &= ~-td;
          ret.push_back({v, u});
        }
      }
    }
  }
  return ret;
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

bool Graph::IsClique(Bitset bs) const {
  int chunks = bs.chunks_;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        if (bs.data_[j] & (~adj_mat2_[v].data_[j])) {
          return false;
        }
      }
    }
  }
  return true;
}

void Graph::RemoveEdgesBetween(int v, const std::vector<int>& vs) {
  for (int u : vs) {
    assert(u != v);
    RemoveEdge(u, v);
  }
}

bool Graph::IsAlmostClique(const std::vector<int>& clq) const {
  std::vector<int> rm;
  for (int i=0;i<(int)clq.size();i++){
    for (int ii=i+1;ii<(int)clq.size();ii++){
      if (!HasEdge(clq[i], clq[ii])) {
        if (rm.size() == 0) {
          rm = {clq[i], clq[ii]};
        } else if (rm.size() == 1) {
          if (clq[i] != rm[0] && clq[ii] != rm[0]) {
            return false;
          }
        } else if (rm.size() == 2) {
          if (clq[i] == rm[0]) {
            assert(clq[ii] != rm[1]);
            rm = {rm[0]};
          } else if(clq[i] == rm[1]) {
            assert(clq[ii] != rm[0]);
            rm = {rm[1]};
          } else if(clq[ii] == rm[0]) {
            assert(clq[i] != rm[1]);
            rm = {rm[0]};
          } else if(clq[ii] == rm[1]) {
            assert(clq[i] != rm[0]);
            rm = {rm[1]};
          } else {
            return false;
          }
        } else {
          assert(0);
        }
      }
    }
  }
  return true;
}

std::vector<int> Graph::Distances(const std::vector<int>& start) const {
  assert(start.size() > 0);
  std::vector<int> d(n_);
  for (int i=0;i<n_;i++){
    d[i] = n_+1;
  }
  std::queue<int> q;
  for (int v : start) {
    d[v] = 0;
    q.push(v);
  }
  while (!q.empty()) {
    int x = q.front();
    q.pop();
    for (int nx : adj_list_[x]){
      if (d[nx] == n_+1) {
        d[nx] = d[x]+1;
        q.push(nx);
      }
    }
  }
  return d;
}

std::vector<std::vector<int>> Graph::DistanceMatrix() const {
  std::vector<std::vector<int>> ret;
  for (int i=0;i<n_;i++){
    ret.push_back(Distances({i}));
    assert((int)ret.back().size() == n_);
  }
  return ret;
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
void Graph::InheritMap(const Graph& parent) {
  vertex_map_ = StaticSet<int>(parent.MapBack(vertex_map_.Values()));
}

std::vector<Bitset> Graph::CompNeighsBit(const Bitset& block) const {
  Bitset vis = ~block;
  std::vector<Bitset> ret;
  Bitset ne(n_);
  Bitset sep(n_);
  for (int i=0;i<n_;i++){
    if (vis.Get(i)) {
      ne.CopyFrom(adj_mat2_[i]);
      Dfs2Bit(vis, ne);
      sep.SetAnd(block, ne);
      if (sep.Popcount() > 0) {
        ret.push_back(sep);
      }
    }
  }
  return ret;
}


void Graph::Dfs2Bit(Bitset& vis, Bitset& ne) const {
  int chunks = vis.Chunks();
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        gv &= ~-gv;
        ne |= adj_mat2_[x];
      }
    }
  }
}

bool Graph::IsSimp(int v) const {
	return IsClique(adj_mat2_[v]);
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
	width = max(width, (int)bags[v].size()-1);
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

bool TreeDecomposition::dfs(int x, int v, int p, vector<int>& u) const {
	assert(InBag(x, v));
	assert(u[x] != v);
	u[x] = v;
	bool ok = true;
	for (int nx : tree.Neighbors(x)) {
		if (InBag(nx, v) && nx != p) {
			if (u[nx] == v) {
				return false;
			}
			ok = dfs(nx, v, x, u);
			if (!ok) {
				return false;
			}
		}
	}
	return true;
}

bool TreeDecomposition::Verify(const Graph& graph) const {
	assert(n == graph.n());
	vector<vector<char>> aps(n);
	for (int i = 0; i < n; i++) {
		aps[i].resize(n);
	}
	for (auto bag : bags) {
		for (int v : bag) {
			for (int u : bag) {
				aps[v][u] = 1;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		if (aps[i][i] == 0) return false;
	}
	for (auto e : graph.Edges()) {
		if (aps[e.F][e.S] == 0) return false;
	}
	vector<int> u(bs+1);
	for (int i = 1; i <= bs; i++) {
		u[i] = -1;
	}
	for (int i = 0; i < n; i++) {
		bool f = false;
		for (int j = 1; j <= bs; j++) {
			if (InBag(j, i)) {
				if (!f) {
					bool ok = dfs(j, i, 0, u);
					if (!ok) return false;
					f = true;
				}
				if (u[j] != i) {
					return false;
				}
			}
		}
	}
	return true;
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
  for (int i = 0; i < n; i++) {
    assert(ret[i] > 0);
  }
  return ret;
}
} // namespace sspp