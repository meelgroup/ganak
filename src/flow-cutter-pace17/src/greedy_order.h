#ifndef GREEDY_ORDER_H
#define GREEDY_ORDER_H

#include "array_id_func.h"

ArrayIDIDFunc compute_greedy_min_degree_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);
ArrayIDIDFunc compute_greedy_min_shortcut_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);

#endif
