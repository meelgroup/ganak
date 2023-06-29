
#ifndef FLOW_CUTTER_CONFIG_H
#define FLOW_CUTTER_CONFIG_H
#include <string>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <cassert>

namespace flow_cutter{
	struct Config{
		int cutter_count;
		int random_seed;
		int max_cut_size;
		float min_small_side_size;

		enum class SkipNonMaximumSides{
			skip,
			no_skip
		};
		SkipNonMaximumSides skip_non_maximum_sides;

		enum class SeparatorSelection{
			node_min_expansion,
			edge_min_expansion,
			node_first,
			edge_first
		};
		SeparatorSelection separator_selection;

		enum class GraphSearchAlgorithm{
			pseudo_depth_first_search,
			breadth_first_search,
			depth_first_search
		};
		GraphSearchAlgorithm graph_search_algorithm;

		enum class AvoidAugmentingPath{
			avoid_and_pick_best,
			do_not_avoid,
			avoid_and_pick_oldest,
			avoid_and_pick_random
		};
		AvoidAugmentingPath avoid_augmenting_path;

		enum class PierceRating{
			max_target_minus_source_hop_dist,
			min_source_hop_dist,
			max_target_hop_dist,
			random,
			oldest
		};
		PierceRating pierce_rating;

		Config():
			cutter_count(3),
			random_seed(5489),
			max_cut_size(1000),
			min_small_side_size(0.2),
			skip_non_maximum_sides(SkipNonMaximumSides::skip),
			separator_selection(SeparatorSelection::node_min_expansion),
			graph_search_algorithm(GraphSearchAlgorithm::pseudo_depth_first_search),
			avoid_augmenting_path(AvoidAugmentingPath::avoid_and_pick_best),
			pierce_rating(PierceRating::max_target_minus_source_hop_dist){}

		void set(const std::string&var, const std::string&val){
			int val_id = -1;
			if(var == "SkipNonMaximumSides" || var == "skip_non_maximum_sides"){
				if(val == "skip" || val_id == static_cast<int>(SkipNonMaximumSides::skip))
					skip_non_maximum_sides = SkipNonMaximumSides::skip;
				else if(val == "no_skip" || val_id == static_cast<int>(SkipNonMaximumSides::no_skip))
					skip_non_maximum_sides = SkipNonMaximumSides::no_skip;
				else throw std::runtime_error("Unknown config value "+val+" for variable SkipNonMaximumSides; valid are skip, no_skip");
			}else if(var == "SeparatorSelection" || var == "separator_selection"){
				if(val == "node_min_expansion" || val_id == static_cast<int>(SeparatorSelection::node_min_expansion))
					separator_selection = SeparatorSelection::node_min_expansion;
				else if(val == "edge_min_expansion" || val_id == static_cast<int>(SeparatorSelection::edge_min_expansion))
					separator_selection = SeparatorSelection::edge_min_expansion;
				else if(val == "node_first" || val_id == static_cast<int>(SeparatorSelection::node_first))
					separator_selection = SeparatorSelection::node_first;
				else if(val == "edge_first" || val_id == static_cast<int>(SeparatorSelection::edge_first))
					separator_selection = SeparatorSelection::edge_first;
				else throw std::runtime_error("Unknown config value "+val+" for variable SeparatorSelection; valid are node_min_expansion, edge_min_expansion, node_first, edge_first");
			}else if(var == "GraphSearchAlgorithm" || var == "graph_search_algorithm"){
				if(val == "pseudo_depth_first_search" || val_id == static_cast<int>(GraphSearchAlgorithm::pseudo_depth_first_search))
					graph_search_algorithm = GraphSearchAlgorithm::pseudo_depth_first_search;
				else if(val == "breadth_first_search" || val_id == static_cast<int>(GraphSearchAlgorithm::breadth_first_search))
					graph_search_algorithm = GraphSearchAlgorithm::breadth_first_search;
				else if(val == "depth_first_search" || val_id == static_cast<int>(GraphSearchAlgorithm::depth_first_search))
					graph_search_algorithm = GraphSearchAlgorithm::depth_first_search;
				else throw std::runtime_error("Unknown config value "+val+" for variable GraphSearchAlgorithm; valid are pseudo_depth_first_search, breadth_first_search, depth_first_search");
			}else if(var == "AvoidAugmentingPath" || var == "avoid_augmenting_path"){
				if(val == "avoid_and_pick_best" || val_id == static_cast<int>(AvoidAugmentingPath::avoid_and_pick_best))
					avoid_augmenting_path = AvoidAugmentingPath::avoid_and_pick_best;
				else if(val == "do_not_avoid" || val_id == static_cast<int>(AvoidAugmentingPath::do_not_avoid))
					avoid_augmenting_path = AvoidAugmentingPath::do_not_avoid;
				else if(val == "avoid_and_pick_oldest" || val_id == static_cast<int>(AvoidAugmentingPath::avoid_and_pick_oldest))
					avoid_augmenting_path = AvoidAugmentingPath::avoid_and_pick_oldest;
				else if(val == "avoid_and_pick_random" || val_id == static_cast<int>(AvoidAugmentingPath::avoid_and_pick_random))
					avoid_augmenting_path = AvoidAugmentingPath::avoid_and_pick_random;
				else throw std::runtime_error("Unknown config value "+val+" for variable AvoidAugmentingPath; valid are avoid_and_pick_best, do_not_avoid, avoid_and_pick_oldest, avoid_and_pick_random");
			}else if(var == "PierceRating" || var == "pierce_rating"){
				if(val == "max_target_minus_source_hop_dist" || val_id == static_cast<int>(PierceRating::max_target_minus_source_hop_dist))
					pierce_rating = PierceRating::max_target_minus_source_hop_dist;
				else if(val == "min_source_hop_dist" || val_id == static_cast<int>(PierceRating::min_source_hop_dist))
					pierce_rating = PierceRating::min_source_hop_dist;
				else if(val == "max_target_hop_dist" || val_id == static_cast<int>(PierceRating::max_target_hop_dist))
					pierce_rating = PierceRating::max_target_hop_dist;
				else if(val == "random" || val_id == static_cast<int>(PierceRating::random))
					pierce_rating = PierceRating::random;
				else if(val == "oldest" || val_id == static_cast<int>(PierceRating::oldest))
					pierce_rating = PierceRating::oldest;
				else throw std::runtime_error("Unknown config value "+val+" for variable PierceRating; valid are max_target_minus_source_hop_dist, min_source_hop_dist, max_target_hop_dist, random, oldest");
			}else if(var == "cutter_count"){
				int x = std::stoi(val);
				if(!(x>0))
					throw std::runtime_error("Value for \"cutter_count\" must fullfill \"x>0\"");
				cutter_count = x;
			}else if(var == "random_seed"){
				random_seed = std::stoi(val);
			}else if(var == "max_cut_size"){
				int x = std::stoi(val);
				if(!(x>=1))
					throw std::runtime_error("Value for \"max_cut_size\" must fullfill \"x>=1\"");
				max_cut_size = x;
			}else if(var == "min_small_side_size"){
				double x = std::stod(val);
				if(!(0.5>=x&&x>=0.0))
					throw std::runtime_error("Value for \"min_small_side_size\" must fullfill \"0.5>=x&&x>=0.0\"");
				min_small_side_size = x;
			}else throw std::runtime_error("Unknown config variable "+var+"; valid are SkipNonMaximumSides, SeparatorSelection, GraphSearchAlgorithm, AvoidAugmentingPath, PierceRating, cutter_count, random_seed, max_cut_size, min_small_side_size");
		}
		std::string get(const std::string&var)const{
			if(var == "SkipNonMaximumSides" || var == "skip_non_maximum_sides"){
				if(skip_non_maximum_sides == SkipNonMaximumSides::skip) return "skip";
				else if(skip_non_maximum_sides == SkipNonMaximumSides::no_skip) return "no_skip";
				else {assert(false); return "";}
			}else if(var == "SeparatorSelection" || var == "separator_selection"){
				if(separator_selection == SeparatorSelection::node_min_expansion) return "node_min_expansion";
				else if(separator_selection == SeparatorSelection::edge_min_expansion) return "edge_min_expansion";
				else if(separator_selection == SeparatorSelection::node_first) return "node_first";
				else if(separator_selection == SeparatorSelection::edge_first) return "edge_first";
				else {assert(false); return "";}
			}else if(var == "GraphSearchAlgorithm" || var == "graph_search_algorithm"){
				if(graph_search_algorithm == GraphSearchAlgorithm::pseudo_depth_first_search) return "pseudo_depth_first_search";
				else if(graph_search_algorithm == GraphSearchAlgorithm::breadth_first_search) return "breadth_first_search";
				else if(graph_search_algorithm == GraphSearchAlgorithm::depth_first_search) return "depth_first_search";
				else {assert(false); return "";}
			}else if(var == "AvoidAugmentingPath" || var == "avoid_augmenting_path"){
				if(avoid_augmenting_path == AvoidAugmentingPath::avoid_and_pick_best) return "avoid_and_pick_best";
				else if(avoid_augmenting_path == AvoidAugmentingPath::do_not_avoid) return "do_not_avoid";
				else if(avoid_augmenting_path == AvoidAugmentingPath::avoid_and_pick_oldest) return "avoid_and_pick_oldest";
				else if(avoid_augmenting_path == AvoidAugmentingPath::avoid_and_pick_random) return "avoid_and_pick_random";
				else {assert(false); return "";}
			}else if(var == "PierceRating" || var == "pierce_rating"){
				if(pierce_rating == PierceRating::max_target_minus_source_hop_dist) return "max_target_minus_source_hop_dist";
				else if(pierce_rating == PierceRating::min_source_hop_dist) return "min_source_hop_dist";
				else if(pierce_rating == PierceRating::max_target_hop_dist) return "max_target_hop_dist";
				else if(pierce_rating == PierceRating::random) return "random";
				else if(pierce_rating == PierceRating::oldest) return "oldest";
				else {assert(false); return "";}
			}else if(var == "cutter_count"){
				return std::to_string(cutter_count);
			}else if(var == "random_seed"){
				return std::to_string(random_seed);
			}else if(var == "max_cut_size"){
				return std::to_string(max_cut_size);
			}else if(var == "min_small_side_size"){
				return std::to_string(min_small_side_size);
			}else throw std::runtime_error("Unknown config variable "+var+"; valid are SkipNonMaximumSides,SeparatorSelection,GraphSearchAlgorithm,AvoidAugmentingPath,PierceRating, cutter_count, random_seed, max_cut_size, min_small_side_size");
		}
		std::string get_config()const{
			std::ostringstream out;
			out
				<< std::setw(30) << "SkipNonMaximumSides" << " : " << get("SkipNonMaximumSides") << '\n'
				<< std::setw(30) << "SeparatorSelection" << " : " << get("SeparatorSelection") << '\n'
				<< std::setw(30) << "GraphSearchAlgorithm" << " : " << get("GraphSearchAlgorithm") << '\n'
				<< std::setw(30) << "AvoidAugmentingPath" << " : " << get("AvoidAugmentingPath") << '\n'
				<< std::setw(30) << "PierceRating" << " : " << get("PierceRating") << '\n'
				<< std::setw(30) << "cutter_count" << " : " << get("cutter_count") << '\n'
				<< std::setw(30) << "random_seed" << " : " << get("random_seed") << '\n'
				<< std::setw(30) << "max_cut_size" << " : " << get("max_cut_size") << '\n'
				<< std::setw(30) << "min_small_side_size" << " : " << get("min_small_side_size") << '\n';
			return out.str();
		}

	};
}
#endif
