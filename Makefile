CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic

.PHONY: constrained

constrained: bin/sequence_constrained_mountain_centroid

bin/sequence_constrained_mountain_centroid: cpp/sequence_constrained_mountain_centroid.cpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@
