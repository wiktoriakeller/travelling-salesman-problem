#pragma once
#include <cmath>
#include <array>
#include "InstanceGenerator.h"
#include "GeneticAlgorithm.h"
#include "Greedy.h"
#include "ParallelGeneticAlgorithm.h"

class Manager {
public:
	void TSPSolver();
	void ParallelTSPSolver();
	void TestParameters();
	void RunTests();
	void RunTestsParallel();

private:
	float** PointsToMatrix(std::vector<Point>& points, int& size);
	std::vector<Point> LoadPointsFromFile(std::string filePath);
	std::vector<Point> PrepareFile(std::string filePath);
	void DeleteMatrix(float** matrix, int size);
	std::vector<std::string> Split(std::string word, char mark);
	void SavePointsToFile(std::string filePath, std::vector<int> chromosome, std::vector<Point> points);
};
