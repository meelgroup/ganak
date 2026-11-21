CXX ?= c++
CXXFLAGS ?= -O3
all: resume enumerate

resume: resume.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++20 $< $(LDFLAGS) -lgmpxx -lgmp -o $@

enumerate: enumerate.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++20 $< $(LDFLAGS) -lcadical -o $@

.PHONY: all
