/*
 Copyright (c) 2016, Ben Strasser
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice, this list
 of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice, this
 list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 Kenji Hashimoto has to say:
 *   This file is a degraded copy of pace.cpp in flow-cutter-pace17.
 *   I made a small change to the source code for our use (e.g., timeout).
 *
 *  Editted on: 2022/03/29
 *      Author: k-hasimt
 */

//#include "flow-cutter-pace17/src/id_func.hpp"
//#include "flow-cutter-pace17/src/list_graph.hpp"
/* #include "flow-cutter-pace17/src/multi_arc.hpp" */
//#include "flow-cutter-pace17/src/sort_arc.hpp"
//#include "flow-cutter-pace17/src/chain.hpp"
/* #include "flow-cutter-pace17/src/union_find.hpp" */
//#include "flow-cutter-pace17/src/node_flow_cutter.hpp"
#include "flow-cutter-pace17/src/separator.hpp"
#include "flow-cutter-pace17/src/id_multi_func.hpp"
#include "flow-cutter-pace17/src/filter.hpp"
#include "flow-cutter-pace17/src/preorder.hpp"
#include "flow-cutter-pace17/src/contraction_graph.hpp"
//#include "flow-cutter-pace17/src/tree_decomposition.hpp"
#include "flow-cutter-pace17/src/greedy_order.hpp"
#include "flow-cutter-pace17/src/min_max.hpp"
#include "flow-cutter-pace17/src/heap.hpp"

#include "IFlowCutter.hpp"
#include "TreeDecomposition.hpp"
#include "time_mem.hpp"

// #include <limits>
// #include <signal.h>
// #include <stdlib.h>
// #include <string.h>
// #include <string>
// #include <sstream>
#include <queue>

#include <sys/time.h>
// #include <unistd.h>

#include <iostream>
using namespace std;
using namespace TWD;


IFlowCutter::IFlowCutter(int n, int m, int _verb) :
    nodes(n), best_bag_size(numeric_limits<int>::max()), head(2*m, n), tail(2*m, n),
    verb(_verb) {
      start_time = cpu_time();
    }

void IFlowCutter::importGraph(const Graph& g)
{
  int next_arc = 0;
  for(int i=0; i<nodes; i++)
    for(auto j : g.Neighbors(i)) {
      if(i >= j) continue;
      assert(next_arc < head.preimage_count());
      head[next_arc] = j;
      tail[next_arc] = i;
      next_arc++;
      head[next_arc] = i;
      tail[next_arc] = j;
      next_arc++;
    }
}

void ignore_return_value(long long){}

// This hack is actually standard compilant
template <class T, class S, class C>
S& access_internal_vector(std::priority_queue<T, S, C>& q) {
  struct Hacked : private priority_queue<T, S, C> {
    static S& access(priority_queue<T, S, C>& q) {
      return q.*&Hacked::c;
    }
  };
  return Hacked::access(q);
}

void IFlowCutter::print_comment(std::string msg){
  if (verb == 0) return;
  cout << "c o [td] " << msg << endl;
}

template<class Tail, class Head>
void check_multilevel_partition_invariants(
    [[maybe_unused]] const Tail&tail, [[maybe_unused]] const Head&head,
    [[maybe_unused]] const std::vector<Cell>&multilevel_partition){
  #ifdef SLOW_DEBUG
  const int node_count = tail.image_count();
  const int arc_count = tail.preimage_count();

  auto is_child_of = [&](int c, int p){
    for(;;){
      if(c == p)
        return true;
      if(c == -1)
        return false;
      c = multilevel_partition[c].parent_cell;
    }
  };

  auto are_ordered = [&](int a, int b){
    return is_child_of(a, b) || is_child_of(b, a);
  };

  ArrayIDFunc<int> cell_of_node(node_count);
  cell_of_node.fill(-1);

  for(int i=0; i<(int)multilevel_partition.size(); ++i){
    for(auto&y:multilevel_partition[i].separator_node_list){
      assert(cell_of_node(y) == -1);
      cell_of_node[y] = i;
    }
  }

  for(auto x:cell_of_node)
    assert(x != -1);

  for(int xy = 0; xy < arc_count; ++xy){
    int x = cell_of_node(tail(xy)), y = cell_of_node(head(xy));
    assert(are_ordered(x, y));
  }
  #endif
}

template <
    class T,
    class Container = std::vector<T>,
    class Compare = std::less<typename Container::value_type>>
class my_priority_queue : public std::priority_queue<T, Container, Compare> {
public:
  T top_and_pop() {
    std::pop_heap(c.begin(), c.end(), comp);
    T value = std::move(c.back());
    c.pop_back();
    return value;
  }

protected:
  using std::priority_queue<T, Container, Compare>::c;
  using std::priority_queue<T, Container, Compare>::comp;
};

template<class Tail, class Head, class ComputeSeparator, class OnNewMP>
void compute_multilevel_partition(const Tail&tail, const Head&head, const ComputeSeparator&compute_separator, int smallest_known_treewidth, const OnNewMP&on_new_multilevel_partition){

  const int node_count = tail.image_count();
  const int arc_count = tail.preimage_count();

  std::vector<Cell>closed_cells;
  my_priority_queue<Cell>open_cells;

  {
    Cell top_level_cell;
    top_level_cell.separator_node_list.resize(node_count);
    for(int i=0; i<node_count; ++i)
      top_level_cell.separator_node_list[i] = i;
    //top_level_cell.boundary_node_list = {};
    top_level_cell.parent_cell = -1;

    open_cells.push(std::move(top_level_cell));
  }

  int max_closed_bag_size = 0;
  int max_open_bag_size = node_count;

  auto check_if_better = [&]{
    int current_tree_width = std::max(max_closed_bag_size, max_open_bag_size);

    if(current_tree_width < smallest_known_treewidth){
      smallest_known_treewidth = current_tree_width;

      std::vector<Cell>cells = closed_cells;
      for(auto&q:access_internal_vector(open_cells))
        cells.push_back(q);
      check_multilevel_partition_invariants(tail, head, cells);
      on_new_multilevel_partition(cells, open_cells.empty() || max_closed_bag_size>=max_open_bag_size);
    }
  };

  check_if_better();


  ArrayIDFunc<int>node_to_sub_node(node_count);
  node_to_sub_node.fill(-1);

  auto inv_tail = invert_sorted_id_id_func(tail);

  BitIDFunc in_child_cell(node_count);
  in_child_cell.fill(false);

  while(!open_cells.empty()){

    #ifdef SLOW_DEBUG

    int real_max_closed_bag_size = 0;
    for(auto&x:closed_cells)
      max_to(real_max_closed_bag_size, (int)x.bag_size());
    assert(max_closed_bag_size == real_max_closed_bag_size);

    int real_max_open_bag_size = 0;
    for(auto&x:access_internal_vector(open_cells))
      max_to(real_max_open_bag_size, (int)x.bag_size());
    assert(max_open_bag_size == real_max_open_bag_size);

    #endif

    Cell current_cell = open_cells.top_and_pop();
    bool must_recompute_max_open_bag_size = ((int)current_cell.bag_size() == max_open_bag_size);

    int closed_cell_id = closed_cells.size();

    if((int)current_cell.bag_size() > max_closed_bag_size){

      auto interior_node_list = std::move(current_cell.separator_node_list);
      int interior_node_count = interior_node_list.size();

      ArrayIDFunc<int>sub_node_to_node(interior_node_count);

      int next_sub_id = 0;
      for(int x:interior_node_list){
        node_to_sub_node[x] = next_sub_id;
        sub_node_to_node[next_sub_id] = x;
        ++next_sub_id;
      }

      auto is_node_interior = id_func(
        node_count,
        [&](int x)->bool{
          return node_to_sub_node(x) != -1;
        }
      );

      auto is_arc_interior = id_func(
        arc_count,
        [&](int xy)->bool{
          return is_node_interior(tail(xy)) && is_node_interior(head(xy));
        }
      );

      int interior_arc_count = count_true(is_arc_interior);
      auto sub_tail = keep_if(is_arc_interior, interior_arc_count, tail);
      auto sub_head = keep_if(is_arc_interior, interior_arc_count, head);

      for(auto&x:sub_tail)
        x = node_to_sub_node(x);
      sub_tail.set_image_count(interior_node_count);

      for(auto&x:sub_head)
        x = node_to_sub_node(x);
      sub_head.set_image_count(interior_node_count);

      auto sub_separator = compute_separator(sub_tail, sub_head);

      BitIDFunc is_in_sub_separator(interior_node_count);
      is_in_sub_separator.fill(false);
      for(auto x:sub_separator)
        is_in_sub_separator.set(x, true);

      UnionFind uf(interior_node_count);

      for(int xy=0; xy<interior_arc_count; ++xy){
        int x = sub_tail(xy);
        int y = sub_head(xy);
        if(!is_in_sub_separator(x) && !is_in_sub_separator(y))
          uf.unite(x, y);
      }

      std::vector<std::vector<int>>nodes_of_representative(interior_node_count);
      for(int x=0; x<interior_node_count; ++x)
        if(!is_in_sub_separator(x))
          nodes_of_representative[uf(x)].push_back(x);

      auto&separator = sub_separator;
      for(auto&x:separator)
        x = sub_node_to_node(x);

      for(int x=0; x<interior_node_count; ++x){
        if(!nodes_of_representative[x].empty()){
          Cell new_cell;

          auto&new_cell_interior_node_list = nodes_of_representative[x];
          for(auto&x2:new_cell_interior_node_list) x2 = sub_node_to_node(x2);

          new_cell.parent_cell = closed_cell_id;

          new_cell.separator_node_list = std::move(new_cell_interior_node_list);

          new_cell.boundary_node_list = current_cell.boundary_node_list;
          new_cell.boundary_node_list.insert(new_cell.boundary_node_list.end(), separator.begin(), separator.end());

          {
            for(auto x2:new_cell.separator_node_list) in_child_cell.set(x2, true);
            new_cell.boundary_node_list.erase(
              std::remove_if(
                new_cell.boundary_node_list.begin(),
                new_cell.boundary_node_list.end(),
                [&](int x2)->bool{
                  for(auto xy:inv_tail(x2))
                    if(in_child_cell(head(xy)))
                      return false;
                  return true;
                }
              ),
              new_cell.boundary_node_list.end()
            );
            for(auto x2:new_cell.separator_node_list) in_child_cell.set(x2, false);
          }

          new_cell.separator_node_list.shrink_to_fit();
          new_cell.boundary_node_list.shrink_to_fit();

          if((int)new_cell.bag_size() > max_open_bag_size)
            max_open_bag_size = new_cell.bag_size();

          open_cells.push(std::move(new_cell));
        }
      }

      current_cell.separator_node_list = std::move(separator);
      current_cell.separator_node_list.shrink_to_fit();

      for(int x2:interior_node_list) node_to_sub_node[x2] = -1;
    }

    if((int)current_cell.bag_size() > max_closed_bag_size)
      max_closed_bag_size = current_cell.bag_size();

    if(must_recompute_max_open_bag_size){
      max_open_bag_size = 0;
      for(auto&x2:access_internal_vector(open_cells))
        if((int)x2.bag_size() > max_open_bag_size) max_open_bag_size = x2.bag_size();
    }

    closed_cells.push_back(std::move(current_cell));

    check_if_better();

    if(max_closed_bag_size >= smallest_known_treewidth){
      return;
    }

    if(max_closed_bag_size >= max_open_bag_size){
      return;
    }
  }
}

int IFlowCutter::compute_max_bag_size_of_order(const ArrayIDIDFunc&order){
  auto inv_order = inverse_permutation(order);
  int current_tail = -1;
  int current_tail_up_deg = 0;
  int max_up_deg = 0;
  compute_chordal_supergraph(
    chain(tail, inv_order), chain(head, inv_order),
    [&](int x, int /*y*/){
      if(current_tail != x){
        current_tail = x;
        max_to(max_up_deg, current_tail_up_deg);
        current_tail_up_deg = 0;
      }
      ++current_tail_up_deg;
    }
  );
  return max_up_deg+1;
}

void IFlowCutter::test_new_order(const ArrayIDIDFunc&order, TreeDecomposition& td){
  int x = compute_max_bag_size_of_order(order);
  {
    if(x < best_bag_size){
      best_bag_size = x;
      td = output_tree_decompostion_of_order(tail, head, order);
    }
  }
}

TreeDecomposition IFlowCutter::output_tree_decompostion_of_order(
    ArrayIDIDFunc tail_loc, ArrayIDIDFunc head_loc, const ArrayIDIDFunc&order){
  TreeDecomposition better_td;

  const int node_count = tail_loc.image_count();

  auto inv_order = inverse_permutation(order);
  tail_loc = chain(tail_loc, inv_order);
  head_loc = chain(head_loc, inv_order);

  vector<vector<int>>nodes_in_bag;
  ArrayIDFunc<vector<int>>bags_of_node(node_count);

  auto is_left_subset_of_right = [](const std::vector<int>&l, const std::vector<int>&r){
    auto i = l.begin(), j = r.begin();

    for(;;){
      if(i == l.end())
        return true;
      if(j == r.end())
        return false;

      if(*i < *j)
        return false;
      if(*i == *j)
        ++i;
      ++j;
    }
  };

  auto compute_intersection_size = [](const std::vector<int>&l, const std::vector<int>&r){
    auto i = l.begin(), j = r.begin();
    int n = 0;
    for(;;){
      if(i == l.end() || j == r.end())
        return n;
      if(*i < *j)
        ++i;
      else if(*i > *j)
        ++j;
      else{
        ++i;
        ++j;
        ++n;
      }
    }
  };

  auto on_new_potential_maximal_clique = [&](int lowest_node_in_clique, std::vector<int>clique){
    for(auto b:bags_of_node(lowest_node_in_clique))
      if(is_left_subset_of_right(clique, nodes_in_bag[b]))
        return;
    int bag_id = nodes_in_bag.size();
    for(auto x:clique)
      bags_of_node[x].push_back(bag_id);
    nodes_in_bag.push_back(std::move(clique));
  };


  {
    BitIDFunc is_root(node_count);
    is_root.fill(true);
    std::vector<int>upper_neighborhood_of_z;
    int z = -1;
    compute_chordal_supergraph(tail_loc, head_loc, [&](int x, int y){
        is_root.set(x, false);
        if(z != -1 && z != x){
          upper_neighborhood_of_z.push_back(z);
          sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
          on_new_potential_maximal_clique(z, std::move(upper_neighborhood_of_z));
          upper_neighborhood_of_z.clear();
        }
        z = x;
        upper_neighborhood_of_z.push_back(y);
      }
    );
    if(z != -1){
      upper_neighborhood_of_z.push_back(z);
      sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
      on_new_potential_maximal_clique(z, std::move(upper_neighborhood_of_z));
    }

    for(int x=0; x<node_count; ++x){
      if(is_root(x)){
        on_new_potential_maximal_clique(x, {x});
      }
    }
  }

  int bag_count = nodes_in_bag.size();

  int maximum_bag_size = 0;
  for(auto&b:nodes_in_bag)
    if((int)b.size() > maximum_bag_size)
      maximum_bag_size = b.size();

  better_td.init(bag_count);
  better_td.setWidth(maximum_bag_size-1);
  better_td.setNumGraphNodes(node_count);
  if (verb > 0) {
  cout << "c o [td] #bags " << bag_count << ", tw " << maximum_bag_size
      << ", elapsed " << (cpu_time()-start_time)  << " s"<< endl; // << endl;// << ", #vars " << node_count << endl;
  }
  better_td.initBags();

  for(int i=0; i<bag_count; ++i) {
    vector<vector<int>>& bags = better_td.Bags();
    for(auto x: nodes_in_bag[i])
      bags[i].push_back(order(x));
  }

  {
    std::vector<int>tail2, head2, weight;

    for(int b=0; b<bag_count; ++b){
      vector<int>neighbor_bags;
      for(auto x:nodes_in_bag[b]){
        vector<int>tmp;
        std::set_union(
          bags_of_node[x].begin(), bags_of_node[x].end(),
          neighbor_bags.begin(), neighbor_bags.end(),
          std::back_inserter(tmp)
        );
        neighbor_bags.swap(tmp);
      }
      for(auto p:neighbor_bags){
        if(p != b){
          tail2.push_back(b);
          head2.push_back(p);
          weight.push_back(compute_intersection_size(nodes_in_bag[b], nodes_in_bag[p]));
        }
      }
    }

    int arc_count = tail2.size();

    auto out_arc = invert_id_id_func(
      id_id_func(
        arc_count, bag_count,
        [&](unsigned a){return tail2[a];}
      )
    );

    BitIDFunc in_tree(bag_count);
    in_tree.fill(false);
    max_id_heap<int>q(arc_count);

    for(int b=0; b<bag_count; ++b){
      if(!in_tree(b)){
        if(b != 0)
          better_td.addEdge(0, b);
        in_tree.set(b, true);
        for(int a:out_arc(b))
          q.push(a, weight[a]);
        while(!q.empty()){
          int xy = q.pop();

          int x = tail2[xy];
          int y = head2[xy];

          assert(in_tree(x));

          if(!in_tree(y)){
            better_td.addEdge(x, y);
            in_tree.set(y, true);
            for(int yz:out_arc(y)){
              assert(!q.contains(yz));
              q.push(yz, weight[yz]);
            }
          }
        }
      }
    }
  }

  return better_td;
}

TreeDecomposition IFlowCutter::output_tree_decompostion_of_multilevel_partition(
    const ArrayIDIDFunc&/*tail*/, const ArrayIDIDFunc&/*head*/, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list){
  TreeDecomposition better_td;

  int bag_count = cell_list.size();

  better_td.init(bag_count);
  better_td.setWidth(get_treewidth_of_multilevel_partition(cell_list)-1);
  better_td.setNumGraphNodes(get_node_count_of_multilevel_partition(cell_list));
  if (verb > 0)
    cout << "c o [td] #bags " << bag_count
      << " tw " << get_treewidth_of_multilevel_partition(cell_list)-1
      << " elapsed " << cpu_time()-start_time  << " s"
      << endl; // << ", #vars " << get_node_count_of_multilevel_partition(cell_list) << endl;
  better_td.initBags();

  for(int i=0; i<bag_count; ++i) {
    vector<vector<int>>& bags = better_td.Bags();
    for(auto&x:cell_list[i].separator_node_list)
      bags[i].push_back(to_input_node_id(x));
    for(auto&x:cell_list[i].boundary_node_list)
      bags[i].push_back(to_input_node_id(x));
  }

  for(int i=0; i<bag_count; ++i){
    if(cell_list[i].parent_cell != -1)
      better_td.addEdge(i, cell_list[i].parent_cell);
  }

  return better_td;
}

TreeDecomposition IFlowCutter::constructTD(const int64_t conf_steps, const int conf_iters)
{
  TreeDecomposition td;
  ArrayIDIDFunc preorder, inv_preorder;
  double t = cpu_time();

  /* int random_seed = 0; */
  try{
    {
      preorder = compute_preorder(compute_successor_function(tail, head));
      for(int i=0; i<tail.image_count(); ++i)
        preorder[i] = i;
      inv_preorder = inverse_permutation(preorder);
      tail = chain(std::move(tail), inv_preorder);
      head = chain(std::move(head), inv_preorder);
    }

    {
      auto p = sort_arcs_first_by_tail_second_by_head(tail, head);
      tail = chain(p, std::move(tail));
      head = chain(p, std::move(head));
    }

    const int node_count = tail.image_count();

    auto on_new_multilevel_partition = [&](const std::vector<Cell>&multilevel_partition, bool /*must_print*/){
      int tw = get_treewidth_of_multilevel_partition(multilevel_partition);
      {
        /* update best tree decomposition*/
        td = output_tree_decompostion_of_multilevel_partition(tail, head, preorder, multilevel_partition);
        best_bag_size = tw;
      }
    };


    {
      try{
        std::minstd_rand rand_gen;
        rand_gen.seed(0);

        if(node_count > 500000)
        {
          print_comment("start F1 with 0.1 min balance and edge_first");
          flow_cutter::Config config;
          config.cutter_count = 1;
          config.random_seed = rand_gen();
          config.min_small_side_size = 0.1;
          config.max_cut_size = 500;
          config.separator_selection = flow_cutter::Config::SeparatorSelection::edge_first;
          compute_multilevel_partition(tail, head, flow_cutter::ComputeSeparator(config), best_bag_size, on_new_multilevel_partition);
        }

        int64_t steps = conf_steps;
        int64_t next_step_print = steps-1e4;
        if(node_count < 50000){
          print_comment("min degree heuristic");
          test_new_order(chain(compute_greedy_min_degree_order(tail, head), inv_preorder), td);
        }

        if(node_count < 1000){
          print_comment("min shortcut heuristic");
          test_new_order(chain(compute_greedy_min_shortcut_order(tail, head), inv_preorder), td);
        }

        {
          print_comment("run with 0.0/0.1/0.2 min balance and node_min_expansion in endless loop with varying seed");
          flow_cutter::Config config;
          config.cutter_count = 1;
          config.random_seed = rand_gen();
          config.max_cut_size = 10000;
          config.separator_selection = flow_cutter::Config::SeparatorSelection::node_min_expansion;

          for(int i=2; i < conf_iters && steps > 0;++i){
            /* cout << "nodes: " << nodes << " preimage count: " << head.preimage_count_ */
            /*   << " mul: " << ((int64_t)nodes * std::sqrt((int64_t)head.preimage_count_))/50 */
            /*   << endl; */
            steps -= (std::sqrt((int64_t)nodes) * std::sqrt((int64_t)head.preimage_count_))/50;
            config.random_seed = rand_gen();
            if(i % 10 == 0) ++config.cutter_count;
            if(i % 100 == 0) config.cutter_count = 1;

            switch(i % 3){
            case 2: config.min_small_side_size = 0.2; break;
            case 1: config.min_small_side_size = 0.1; break;
            case 0: config.min_small_side_size = 0.0; break;
            }

            compute_multilevel_partition(tail, head, flow_cutter::ComputeSeparator(config), best_bag_size, on_new_multilevel_partition);

            if (i % 100 == 99 || steps < next_step_print) {
              if (verb) {
                cout << "c o [td] iter " << i << " width: " << td.width()
                  << " stepsK remain: " << steps/1000 << " T: " << (cpu_time()-t) << endl;
              }
              next_step_print -= 1e5;
            }
          }
        }
      }catch(...){
      }

    }
  }catch(...){
  }

  return td;
}

