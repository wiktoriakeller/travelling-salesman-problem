#pragma once
#include<vector>

struct Chromosome {
	std::vector<int> tour;
	float fitness;
	float path;

	Chromosome(int capacity) : tour(capacity) {}
};