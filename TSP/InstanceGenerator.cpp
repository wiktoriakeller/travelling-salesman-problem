#include "InstanceGenerator.h"

void InstanceGenerator::GenerateInstance(int numberOfCities, int canvasSize, std::string filePath) {
	std::vector<Point> points;
	srand(time(NULL));

	for (int i = 0; i < numberOfCities; i++) {
		int x = rand() % canvasSize;
		int y = rand() % canvasSize;

		if (Overlaps(x, y, points)) {
			i--;
		}
		else {
			points.emplace_back(Point(x, y));
		}
	}

	SaveToFile(points, filePath);
}

bool InstanceGenerator::Overlaps(int x, int y, const std::vector<Point>& points) {
	for (int i = 0; i < points.size(); i++) {
		if (points[i].GetX() == x && points[i].GetY() == y) {
			return true;
		}
	}

	return false;
}

void InstanceGenerator::SaveToFile(const std::vector<Point>& points, std::string filePath) {
	std::ofstream outputFile(filePath);

	outputFile << points.size() << "\n";

	for (int i = 0; i < points.size(); i++) {
		outputFile << i + 1 << " " << points[i].GetX() << " " << points[i].GetY() << "\n";
	}

	outputFile.close();
}
