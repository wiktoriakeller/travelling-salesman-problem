#include "Manager.h"
#include <time.h>

void Manager::TSPSolver() {
	int size;
	std::vector<Point> points = LoadPointsFromFile("Instances\\210.txt");
	float** matrix = PointsToMatrix(points, size);
	GeneticAlgorithm genetic(10000, size, 4, 10, 89, 13, matrix);
	genetic.Run(1500);
	DeleteMatrix(matrix, size);
}

void Manager::ParallelTSPSolver() {
	//omp_set_nested(1);
	//omp_set_num_threads(8);

	int size;
	std::vector<Point> points = LoadPointsFromFile("Instances\\210.txt");
	float** matrix = PointsToMatrix(points, size);
	ParallelGeneticAlgorithm parallel(1000, size, 4, 18, 89, 13, matrix);
	parallel.Run(1500);
	DeleteMatrix(matrix, size);
}

float** Manager::PointsToMatrix(std::vector<Point>& points, int& size) {
	size = points.size();

	float** cities = new float* [size];
	for (int i = 0; i < size; ++i) {
		cities[i] = new float[size];
		for (int j = 0; j < size; ++j) {
			cities[i][j] = 0;
		}
	}

	for (int i = 0; i < size; ++i) {
		for (int j = i + 1; j < size; ++j) {
			float distance = sqrt(pow(points[i].GetX() - points[j].GetX(), 2) + pow(points[i].GetY() - points[j].GetY(), 2));
			cities[i][j] = distance;
			cities[j][i] = distance;
		}
	}

	return cities;
}

std::vector<Point> Manager::LoadPointsFromFile(std::string filePath) {
	std::vector<Point> points;
	std::ifstream pointsFile(filePath);

	if (pointsFile.is_open()) {
		int size;
		int number;
		int x, y;

		pointsFile >> size;

		while (pointsFile >> number >> x >> y) {
			points.emplace_back(Point(x, y));
		}
	}

	pointsFile.close();
	return points;
}

std::vector<Point> Manager::PrepareFile(std::string filePath) {
	std::vector<Point> points;
	std::ifstream pointsFile(filePath);

	if (pointsFile.is_open()) {
		std::string line;
		std::vector<std::string> splitted;
		std::getline(pointsFile, line);

		int size;
		int number;
		int x, y;

		while (line != "EOF") {
			splitted = Split(line, ' ');

			if (splitted[0].compare("NAME:") != 0 && splitted[0].compare("TYPE:") != 0 && splitted[0].compare("COMMENT:") != 0
				&& splitted[0].compare("EDGE_WEIGHT_TYPE:") != 0 && splitted[0].compare("NODE_COORD_SECTION") != 0
				&& splitted[0].compare("DIMENSION:") != 0) {

				int numbers[3] = { 0, 0, 0 };
				int j = 0;
				for (int i = 0; i < splitted.size(); i++) {
					if (splitted[i] != "" && splitted[i] != " ") {
						numbers[j++] = std::stoi(splitted[i]);
					}
				}

				points.emplace_back(Point(numbers[1], numbers[2]));
			}
			else if (splitted[0].compare("DIMENSION:") == 0) {
				number = std::stoi(splitted[1]);
			}

			std::getline(pointsFile, line);
		}
	}

	pointsFile.close();
	return points;
}

void Manager::DeleteMatrix(float** matrix, int size) {
	for (int i = 0; i < size; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

std::vector<std::string> Manager::Split(std::string word, char mark) {
	std::string tmp = "";
	std::vector<std::string> splitted;

	for (int i = 0; i < word.size(); i++) {
		if (word[i] == mark) {
			splitted.push_back(tmp);
			tmp.clear();
		}
		else {
			tmp += word[i];
		}
	}
	splitted.push_back(tmp);

	return splitted;
}

void Manager::SavePointsToFile(std::string filePath, std::vector<int> chromosome, std::vector<Point> points) {
	std::ofstream pointsFile(filePath);
	for (auto& i : chromosome) {
		pointsFile << points[i].GetX() << " " << points[i].GetY() << " " << i << '\n';
	}
	pointsFile.close();
}

void Manager::RunTestsParallel() {
	std::vector<std::string> files = { "30.txt", "60.txt", "90.txt", "120.txt", "150.txt",
		"180.txt", "210.txt", "240.txt", "270.txt", "300.txt",
		"330.txt", "360.txt", "390.txt", "420.txt", "450.txt", "berlin52.txt", "kroD100.txt", "kroB150.txt" };

	std::vector<float> bestFound = {
		4188.85, 5756.4, 7191.84, 8067.61, 9168.92, 10012.1, 11853.9, 11532.7, 12276.4, 13201.4, 13571.7, 14268, 14728.2, 15310.8, 15751,
		7544.37, 21399.8, 26127.4
	};

	std::string path = "Instances\\";
	clock_t start, stop;
	std::ofstream result("Results\\parallel-8-watkow.txt");

	for (int i = 0; i < files.size(); i++) {
		int size;
		std::vector<Point> points = LoadPointsFromFile(path + files[i]);
		float** matrix = PointsToMatrix(points, size);

		start = clock();
		ParallelGeneticAlgorithm genetic(10000, size, 4, 20, 89, 13, matrix);
		genetic.Run(1500, bestFound[i]);
		stop = clock();
		double elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC;

		result << files[i] << "\n";
		result << "Best path " << "Best Fitness " << "Time " << "\n";
		result << genetic.bestPath << " " << genetic.bestFitness << " " << elapsed << "\n";

		std::cout << files[i] << "\n";
		std::cout << "Best path " << "Best Fitness " << "Time " << "\n";
		std::cout << genetic.bestPath << " " << genetic.bestFitness << " " << elapsed << "\n";
		DeleteMatrix(matrix, size);
	}
}

void Manager::RunTests() {
	//std::vector<std::string> files = { "30.txt", "60.txt", "90.txt", "120.txt", "150.txt",
	//	"180.txt", "210.txt", "240.txt", "270.txt", "300.txt",
	//	"330.txt", "360.txt", "390.txt", "420.txt", "450.txt" };

	std::vector<std::string> files = { "berlin52.txt", "cities20.txt", 
	"cities40.txt", "kroA200.txt", "kroB150.txt", "kroD100.txt", "tsp250.txt",
		"tsp500.txt", "tsp1000.txt", "pr299.txt"};

	std::string path = "Instances\\";
	clock_t start, stop;
	std::ofstream result("Results\\testSynch4.txt");

	for (int i = 0; i < files.size(); i++) {
		int size;
		std::vector<Point> points = LoadPointsFromFile(path + files[i]);
		float** matrix = PointsToMatrix(points, size);
		
		start = clock();
		GeneticAlgorithm genetic(10000, size, 4, 10, 89, 13, matrix);
		genetic.RunTest(1500);
		stop = clock();
		double elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC;
		
		result << files[i] << "\n";
		result << "Best path " << "Best Fitness " << "Time " << "\n";
		result << genetic.bestPath << " " << genetic.bestFitness << " " << elapsed << "\n";
		
		std::cout << files[i] << "\n";
		std::cout << "Best path " << "Best Fitness " << "Time " << "\n";
		std::cout << genetic.bestPath << " " << genetic.bestFitness << " " << elapsed << "\n";
	}
}

void Manager::TestParameters() {
	int startPopulation = 10000;
	int mutationRate = 4;
	int startNumberOfPairs = 40;
	int chanceToUseCloseCity = 75;
	int twoOptIterations = 10;
	std::ofstream result("Results\\data1.txt");
	float score;
	float bestScore;
	float TheBestScore;
	std::array<std::string, 2> filePaths = { "Instances\\tsp500.txt", "Instances\\tsp250.txt" };
	for (int n = 0; n < filePaths.size(); n++) {
		int size;
		TheBestScore = 1000000000;
		result << filePaths[n] << "\n";
		std::cout << filePaths[n] << "\n";
		std::vector<Point> points = LoadPointsFromFile(filePaths[n]);
		float** matrix = PointsToMatrix(points, size);

		for (int i = startPopulation; i <= 25000; i += 5000) {
			for (int j = startNumberOfPairs; j <= 100; j += 10) {
				for (int g = chanceToUseCloseCity; g <= 100; g += 5) {
					for (int h = twoOptIterations; h <= 40; h += 5) {
						for (int a = 0; a < 2; a++) {
							GeneticAlgorithm genetic(i, size, mutationRate, j, g, h, matrix);
							score = genetic.RunFixedTime(120);
							if (a == 0) {
								bestScore = score;
							}
							else {
								if (score < bestScore) {
									bestScore = score;
								}
							}
						}
						result << bestScore << " MutationRate: " << mutationRate << " twoOptIterations: " << h << " chanceToUseCloseCity: " << g << " NumberOfPairs: " << j << " Population: " << i << "\n";
						std::cout << bestScore << '\n';
						if (bestScore < TheBestScore) {
							TheBestScore = bestScore;
						}
					}
				}
			}
		}
		result << TheBestScore << "\n\n\n";
		DeleteMatrix(matrix, size);
	}
}