#include "list_graph.h"
#include "io_helper.h"
#include "multi_arc.h"
#include "id_multi_func.h"

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>

static
ListGraph load_pace_graph_impl(std::istream&in){
	ListGraph graph;
	std::string line;
	int line_num = 0;
	int next_arc = 0;

	bool was_header_read = false;
	while(std::getline(in, line)){
		++line_num;
		if(line.empty() || line[0] == 'c')
			continue;

		std::istringstream lin(line);
		if(!was_header_read){
			was_header_read = true;
			std::string p, sp;
			int node_count;
			int arc_count;
			if(!(lin >> p >> sp >> node_count >> arc_count))
				throw std::runtime_error("Can not parse header in pace file.");
			if(p != "p" || sp != "tw" || node_count < 0 || arc_count < 0)
				throw std::runtime_error("Invalid header in pace file.");
			graph = ListGraph(node_count, 2*arc_count);
		}else{
			int h, t;
			if(!(lin >> t >> h))
				throw std::runtime_error("Can not parse line num "+std::to_string(line_num)+" \""+line+"\" in pace file.");
			--h;
			--t;
			if(h < 0 || h >= graph.node_count() || t < 0 || t >= graph.node_count())
				throw std::runtime_error("Invalid arc in line num "+std::to_string(line_num)+" \""+line+"\" in pace file.");
			if(next_arc < graph.arc_count()){
				graph.head[next_arc] = h;
				graph.tail[next_arc] = t;
			}
			++next_arc;
			if(next_arc < graph.arc_count()){
				graph.head[next_arc] = t;
				graph.tail[next_arc] = h;
			}
			++next_arc;
		}
	}

	if(next_arc != graph.arc_count())
		throw std::runtime_error("The arc count in the header ("+std::to_string(graph.arc_count())+") does not correspond with the actual number of arcs ("+std::to_string(next_arc)+").");

	return graph; // NVRO
}

ListGraph uncached_load_pace_graph(const std::string&file_name){
	return load_uncached_text_file(file_name, load_pace_graph_impl);
}

