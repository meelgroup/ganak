#ifndef LIST_GRAPH_H
#define LIST_GRAPH_H

#include "array_id_func.h"

#include <tuple>
#include <string>

struct ListGraph{
	ListGraph()=default;
	ListGraph(int node_count, int arc_count)
		:head(arc_count, node_count), tail(arc_count, node_count){}

	int node_count()const{ return head.image_count(); }
	int arc_count()const{ return head.preimage_count(); }

	ArrayIDIDFunc head, tail;
};

ListGraph uncached_load_pace_graph(const std::string&file_name);


#endif
