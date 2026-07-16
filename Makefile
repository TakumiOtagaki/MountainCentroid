CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic

.PHONY: exact

exact: bin/exact_mountain_centroid

bin/exact_mountain_centroid: cpp/exact_mountain_centroid.cpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@
