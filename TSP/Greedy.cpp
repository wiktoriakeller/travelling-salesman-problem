#include "Greedy.h"

std::tuple<int, float> TSPGreedy::FindNearestNeighbour(int city, std::map<int, bool>& visitedCities,
	float** matrix, int size) {
	float minPath = std::numeric_limits<float>::max();
	int nearestNeighbour;

	for (int i = 0; i < size; ++i) {
		if (matrix[city][i] < minPath && !visitedCities[i] && matrix[city][i] != 0) {
			minPath = matrix[city][i];
			nearestNeighbour = i;
		}
	}
	visitedCities[nearestNeighbour] = true;
	return std::make_tuple(nearestNeighbour, minPath);
}

float TSPGreedy::FindShortestPath(int startCity, float** matrix, int size) {
	std::vector<int> path;
	std::map<int, bool> visitedCities;
	int numOfVisitedCities = 0;
	float shortestPath = 0;

	for (int i = 0; i < size; ++i) {
		visitedCities[i] = false;
	}

	int currentCity = startCity;
	visitedCities[startCity] = true;
	numOfVisitedCities++;
	path.emplace_back(startCity);

	while (numOfVisitedCities != size) {
		auto neighbourData = FindNearestNeighbour(currentCity, visitedCities, matrix, size);
		numOfVisitedCities++;
		currentCity = std::get<0>(neighbourData);
		shortestPath += std::get<1>(neighbourData);
		path.emplace_back(currentCity);
	}

	//path.emplace_back(startCity);

	shortestPath += matrix[currentCity][startCity];
	//for (int i = 0; i < path.size(); ++i) {
	//	std::cout << path[i] + 1 << " -> ";
	//}
	//std::cout << "\n" << shortestPath << "\n";

	return shortestPath;
}