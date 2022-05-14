#pragma once
#include <tuple>
#include <limits>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

class TSPGreedy {
public:
	static float FindShortestPath(int startCity, float** matrix, int size);
	static std::tuple<int, float> FindNearestNeighbour(int city, std::map<int, bool>& visitedCities, float** marix, int size);
};