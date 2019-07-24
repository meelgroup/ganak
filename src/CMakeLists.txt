add_subdirectory (component_types)
add_subdirectory (clhash)

SET(SOURCES alt_component_analyzer.cpp 
	component_analyzer.cpp
	component_cache.cpp
	component_management.cpp
	instance.cpp
	new_component_analyzer.cpp
	solver.cpp
	statistics.cpp)

SET(LIB	clhash
		component_types)

ADD_LIBRARY(libganak ${SOURCES} ${LIB})

#SET_TARGET_PROPERTIES(libganak PROPERTIES LINK_SEARCH_START_STATIC 1)
target_link_libraries (libganak ${LIB}) 
#SET_TARGET_PROPERTIES(libganak PROPERTIES LINK_SEARCH_START_STATIC 0)
