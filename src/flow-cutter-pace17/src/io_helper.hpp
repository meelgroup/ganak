#ifndef IO_HELPER_H
#define IO_HELPER_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

template<class SaveFunc, class ...Args>
void save_text_file(const std::string&file_name, const SaveFunc&save, Args&&...args){
	if(file_name == "-"){
		save(std::cout, std::forward<Args>(args)...);
		std::cout << std::flush;
	} else if(file_name == "-null"){
	} else {
		std::ofstream out(file_name);
		if(!out)
			throw std::runtime_error("Could not open "+file_name+" for text writing");
		save(out, std::forward<Args>(args)...);
	}
}

template<class LoadFunc>
auto load_uncached_text_file(const std::string&file_name, const LoadFunc&load)->decltype(load(std::cin)){
	if(file_name == "-"){
		return load(std::cin);
	} else {
		std::ifstream in(file_name);
		if(!in)
			throw std::runtime_error("Could not load "+file_name+" for text reading");
		return load(in);
	}
}

#endif
