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

#pragma once

#include "flow-cutter-pace17/src/array_id_func.hpp"
#include "flow-cutter-pace17/src/cell.hpp"
#include "TreeDecomposition.hpp"

class IFlowCutter {
public:
  IFlowCutter(int n, int m, int verb = 0);

  void importGraph(const Graph& g);
  TreeDecomposition constructTD();
  auto num_nodes() const { return nodes; }

private:
  void print_comment(std::string msg);
  int compute_max_bag_size_of_order(const ArrayIDIDFunc&order);
  void test_new_order(const ArrayIDIDFunc&order, TreeDecomposition&td);

  TreeDecomposition output_tree_decompostion_of_order(ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order);
  TreeDecomposition output_tree_decompostion_of_multilevel_partition(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list);

  int nodes;
  int best_bag_size;

  ArrayIDIDFunc head, tail;
  int verb = 0;
  double start_time;
};
