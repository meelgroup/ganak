#ifndef TREE_DECOMPOSITION_H
#define TREE_DECOMPOSITION_H

#include "array_id_func.h"
#include "cell.h"
#include <string>
#include <ostream>

void print_tree_decompostion_of_order(std::ostream&out, ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order);
void print_tree_decompostion_of_multilevel_partition(std::ostream&out, const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list);

#endif
