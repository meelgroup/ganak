/*
 * FlowCutter.h
 *  This is an interface to flow_cutter_pace17 only for our use.
 *
 *  Created on: 2022/03/29
 *      Author: k-hasimt
 */

#ifndef PREPROCESSOR_IFLOWCUTTER_H_
#define PREPROCESSOR_IFLOWCUTTER_H_

#include "flow-cutter-pace17/src/array_id_func.h"
#include "flow-cutter-pace17/src/cell.h"
#include "TreeDecomposition.h"

class IFlowCutter {
public:
  IFlowCutter(int n, int m, double timeout=120);

  void importGraph(const Graph& g);
  TreeDecomposition constructTD();

  void setTimeout(double _timeout) { this->timeout = _timeout; }

private:
  int compute_max_bag_size_of_order(const ArrayIDIDFunc&order);
  void test_new_order(const ArrayIDIDFunc&order, TreeDecomposition&td);

  TreeDecomposition output_tree_decompostion_of_order(ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order);
  TreeDecomposition output_tree_decompostion_of_multilevel_partition(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list);

  int nodes;
  int best_bag_size;

  ArrayIDIDFunc head, tail;

  double timeout;
};



#endif /* PREPROCESSOR_IFLOWCUTTER_H_ */
