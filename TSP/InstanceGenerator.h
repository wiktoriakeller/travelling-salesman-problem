#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Point.h"

class InstanceGenerator {
public:
	static void GenerateInstance(int numberOfCities, int canvasSize, std::string filePath);

private:
	static bool Overlaps(int x, int y, const std::vector<Point>& points);
	static void SaveToFile(const std::vector<Point>& points, std::string filePath);
};