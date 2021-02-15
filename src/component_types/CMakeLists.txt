

SET(component_types_SOURCES
        component_archetype.cpp
		ClHashComponent.h packed_component.cpp packed_component.h sparsegraph_component.h)

ADD_LIBRARY(component_types ${component_types_SOURCES})