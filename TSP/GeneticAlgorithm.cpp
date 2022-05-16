#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm(int populationSize, int numberOfCities, float mutationRate, int numberOfParentPairs,
	int chanceToUseCloseCity, int twoOptIterations, float** cities) : population(populationSize, numberOfCities) {
	srand(time(NULL));

	this->populationSize = populationSize;
	this->numberOfCities = numberOfCities;
	this->mutationRate = mutationRate;
	this->cities = cities;
	this->numberOfParentPairs = numberOfParentPairs;
	this->chanceToUseCloseCity = chanceToUseCloseCity;
	this->twoOptIterations = twoOptIterations;

	bestFitness = std::numeric_limits<float>::min();
	bestPath = std::numeric_limits<float>::max();
	bestWrittenFitness = std::numeric_limits<float>::min();
	InitializePopulation();
}

void GeneticAlgorithm::SetPopulationSize(int newSize) {
	populationSize = newSize;
}

int GeneticAlgorithm::GetPopulationSize() const {
	return populationSize;
}

void GeneticAlgorithm::SetMutationRate(float newRate) {
	mutationRate = newRate;
}

float GeneticAlgorithm::GetMutationRate() const {
	return mutationRate;
}

void GeneticAlgorithm::SetNumberOfCities(int newNumberOfCities) {
	numberOfCities = newNumberOfCities;
}

int GeneticAlgorithm::GetNumberOfCities() const {
	return numberOfCities;
}

int** GeneticAlgorithm::SetCitiesMatrix(int** matrix) {
	return matrix;
}

void GeneticAlgorithm::ClearPopulation() {
	if (!population.empty()) {
		bestChromosome.clear();
		bestFitness = std::numeric_limits<float>::min();
		bestPath = std::numeric_limits<float>::max();
		bestWrittenFitness = std::numeric_limits<float>::min();
		population.clear();
	}
}

static std::mutex i_mutex;
static std::mutex best_mutex;

void GeneticAlgorithm::GenerateTour(int* i) {
	int iLocal;
	int startCity;
	int nextCity;
	std::vector<bool> visitedCities(numberOfCities, false);
	while (true) {
		{
			std::lock_guard<std::mutex> lock(i_mutex);
			iLocal = *i;
			(*i)++;
		}
		if (iLocal >= populationSize) {
			break;
		}

		startCity = rand() % numberOfCities;
		nextCity = startCity;
		visitedCities[startCity] = true;

		for (int j = 0; j < numberOfCities - 1; ++j) {
			population[iLocal].tour[j] = nextCity;
			if ((rand() % 100) < chanceToUseCloseCity) {
				nextCity = FindNearestNeighbour(population[iLocal].tour[j], visitedCities);
			}
			else {
				std::vector<int> onlyUnvisitedCities;
				onlyUnvisitedCities.reserve(numberOfCities);
				for (int g = 0; g < visitedCities.size(); ++g) {
					if (!visitedCities[g]) {
						onlyUnvisitedCities.emplace_back(g);
					}
				}

				int randomCity = rand() % onlyUnvisitedCities.size();
				nextCity = onlyUnvisitedCities[randomCity];
				visitedCities[nextCity] = true;
			}
		}
		population[iLocal].tour[numberOfCities - 1] = nextCity;
		CalculateFitness(iLocal);


		if (population[iLocal].fitness > bestFitness) {
			std::lock_guard<std::mutex> lock(best_mutex);
			bestFitness = population[iLocal].fitness;
			bestPath = population[iLocal].path;
			bestChromosome = population[iLocal].tour;
		}


		for (int j = 0; j < visitedCities.size(); ++j) {
			visitedCities[j] = false;
		}
	}
}

void GeneticAlgorithm::InitializePopulation() {
	int i = 0;
	std::thread th1(&GeneticAlgorithm::GenerateTour, this, &i);
	std::thread th2(&GeneticAlgorithm::GenerateTour, this, &i);
	std::thread th3(&GeneticAlgorithm::GenerateTour, this, &i);
	std::thread th4(&GeneticAlgorithm::GenerateTour, this, &i);
	th1.join();
	th2.join();
	th3.join();
	th4.join();
}

int GeneticAlgorithm::FindNearestNeighbour(int city, std::vector<bool>& visitedCities) {
	float minPath = std::numeric_limits<float>::max();
	int nearestNeighbour;

	for (int i = 0; i < numberOfCities; ++i) {
		if (cities[city][i] < minPath && !visitedCities[i] && i != city) {
			minPath = cities[city][i];
			nearestNeighbour = i;
		}
	}

	visitedCities[nearestNeighbour] = true;
	return nearestNeighbour;
}

float GeneticAlgorithm::CalculateDistance(std::vector<int>& chromosome) {
	float path = 0;
	for (int i = 0; i < numberOfCities - 1; ++i) {
		path += cities[chromosome[i]][chromosome[i + 1]];
	}

	path += cities[chromosome[numberOfCities - 1]][chromosome[0]];
	return path;
}

void GeneticAlgorithm::CalculateFitness(int index) {
	float path;
	path = CalculateDistance(population[index].tour);
	population[index].path = path;
	population[index].fitness = 1 / path;
}

void GeneticAlgorithm::MakeNextGeneration() {
	int bestParentIndex = 0;
	int bestSecondParentIndex = 0;
	int worstParentIndex = 0;
	int worstSecondParentIndex = 0;

	int bestTour;

	if (numberOfParentPairs == 1) {
		//choosing best, second best, worst and second worst parent
		bestTour = 0;
		for (int i = 1; i < populationSize; i++) {
			if (population[i].fitness > population[bestSecondParentIndex].fitness) {
				if (population[i].fitness > population[bestParentIndex].fitness) {
					bestParentIndex = i;
				}
				else {
					bestSecondParentIndex = i;
				}
			}

			if (population[i].fitness < population[worstSecondParentIndex].fitness) {
				if (population[i].fitness < population[worstParentIndex].fitness) {
					worstParentIndex = i;
				}
				else {
					worstSecondParentIndex = i;
				}
			}
		}

		bestTour = bestParentIndex;

		population[worstParentIndex].tour = Crossover(population[bestParentIndex].tour, population[bestSecondParentIndex].tour);
		Mutate(population[worstParentIndex].tour);
		TwoOpt(population[worstParentIndex].tour);
		CalculateFitness(worstParentIndex);

		if (population[worstParentIndex].fitness > bestFitness) {
			bestFitness = population[worstParentIndex].fitness;
			bestPath = population[worstParentIndex].path;
			bestChromosome = population[worstParentIndex].tour;
		}

		if (population[worstParentIndex].fitness > population[bestTour].fitness) {
			bestTour = worstParentIndex;
		}

		population[worstSecondParentIndex].tour = Crossover(population[bestSecondParentIndex].tour, population[bestParentIndex].tour);
		Mutate(population[worstSecondParentIndex].tour);
		TwoOpt(population[worstSecondParentIndex].tour);
		CalculateFitness(worstSecondParentIndex);

		if (population[worstSecondParentIndex].fitness > bestFitness) {
			bestFitness = population[worstSecondParentIndex].fitness;
			bestPath = population[worstSecondParentIndex].path;
			bestChromosome = population[worstSecondParentIndex].tour;
		}

		if (population[worstSecondParentIndex].fitness > population[bestTour].fitness) {
			bestTour = worstSecondParentIndex;
		}
	}
	else {
		std::sort(population.begin(), population.end(), CompareFitness);
		bestTour = 0;

		for (int i = 0; i < numberOfParentPairs * 2; i += 2) {

			bestParentIndex = i;
			bestSecondParentIndex = i + 1;
			worstParentIndex = populationSize - i - 1;
			worstSecondParentIndex = populationSize - i - 2;

			population[worstParentIndex].tour = Crossover(population[bestParentIndex].tour, population[bestSecondParentIndex].tour);
			Mutate(population[worstParentIndex].tour);
			TwoOpt(population[worstParentIndex].tour);
			CalculateFitness(worstParentIndex);

			if (population[worstParentIndex].fitness > bestFitness) {
				bestFitness = population[worstParentIndex].fitness;
				bestPath = population[worstParentIndex].path;
				bestChromosome = population[worstParentIndex].tour;
			}

			if (population[worstParentIndex].fitness > population[bestTour].fitness) {
				bestTour = worstParentIndex;
			}

			population[worstSecondParentIndex].tour = Crossover(population[bestSecondParentIndex].tour, population[bestParentIndex].tour);
			Mutate(population[worstSecondParentIndex].tour);
			TwoOpt(population[worstSecondParentIndex].tour);
			CalculateFitness(worstSecondParentIndex);

			if (population[worstSecondParentIndex].fitness > bestFitness) {
				bestFitness = population[worstSecondParentIndex].fitness;
				bestPath = population[worstSecondParentIndex].path;
				bestChromosome = population[worstSecondParentIndex].tour;
			}

			if (population[worstSecondParentIndex].fitness > population[bestTour].fitness) {
				bestTour = worstSecondParentIndex;
			}
		}
	}

	TwoOpt(population[bestTour].tour);

	if (population[bestTour].fitness > bestFitness) {
		bestFitness = population[bestTour].fitness;
		bestPath = population[bestTour].path;
		bestChromosome = population[bestTour].tour;
	}
}

bool GeneticAlgorithm::CompareFitness(Chromosome& a, Chromosome& b) {
	return a.fitness > b.fitness;
}

void GeneticAlgorithm::TwoOpt(std::vector<int>& chromosome) {
	for (int n = 0; n < twoOptIterations; ++n) {
		float minchange = 0;
		for (int i = 0; i < chromosome.size() - 2; ++i) {
			for (int j = i + 2; j < chromosome.size() - 1; ++j) {
				float change = GetDistance(chromosome[i], chromosome[j])
					+ GetDistance(chromosome[i + 1], chromosome[j + 1])
					- GetDistance(chromosome[i], chromosome[i + 1])
					- GetDistance(chromosome[j], chromosome[j + 1]);

				if (change < 0) {
					std::swap(chromosome[j], chromosome[i + 1]);
				}

				if (change < minchange) {
					minchange = change;
				}
			}
		}

		if (minchange >= 0) {
			break;
		}
	}
}

float GeneticAlgorithm::GetDistance(int a, int b) {
	return cities[a][b];
}

void GeneticAlgorithm::Mutate(std::vector<int>& chromosome) {
	int randomCityIndex = rand() % chromosome.size();
	int randomShift = rand() % (chromosome.size() - 1) + 1;
	int direction = rand() % 2;
	int i = randomCityIndex + 1;

	if (direction == 0) {
		i = randomCityIndex - 1;
	}

	while (randomShift > 0) {
		if (i >= chromosome.size()) {
			i -= chromosome.size();
		}

		if (i < 0) {
			i = chromosome.size() - 1;
		}

		std::swap(chromosome[randomCityIndex], chromosome[i]);

		randomCityIndex = i;
		randomShift--;

		if (direction == 1) {
			i++;
		}
		else {
			i--;
		}
	}
}

std::vector<int> GeneticAlgorithm::Crossover(std::vector<int>& parent1, std::vector<int>& parent2) {
	//PMX crossover
	std::vector<int> offspring(parent1.size(), -1);
	int crossoverPoint2 = rand() % parent1.size();
	int crossoverPoint1 = rand() % (parent1.size() - crossoverPoint2);

	if (crossoverPoint2 < crossoverPoint1) {
		std::swap(crossoverPoint1, crossoverPoint2);
	}

	for (int i = crossoverPoint1; i <= crossoverPoint2; ++i) {
		offspring[i] = parent1[i];
	}

	for (int i = crossoverPoint1; i <= crossoverPoint2; ++i) {
		if (FindCityIndex(offspring, parent2[i], crossoverPoint1, crossoverPoint2) == -1) {
			int j = parent1[i];
			int index = FindCityIndex(parent2, j, 0, parent2.size() - 1);
			while (index >= crossoverPoint1 && index <= crossoverPoint2) {
				j = parent1[index];
				index = FindCityIndex(parent2, j, 0, parent2.size() - 1);
			}
			offspring[index] = parent2[i];
		}
	}

	for (int i = 0; i < offspring.size(); ++i) {
		if (offspring[i] == -1) {
			offspring[i] = parent2[i];
		}
	}

	return offspring;
}

int GeneticAlgorithm::FindCityIndex(std::vector<int>& chromosome, int city, int start, int end) {
	for (int i = start; i <= end; ++i) {
		if (chromosome[i] == city) {
			return i;
		}
	}
	return -1;
}

void GeneticAlgorithm::Run(int numberOfIterations) {
	for (int i = 0; i < numberOfIterations; ++i) {
		PrintBestChromosome();
		MakeNextGeneration();
	}
}

void GeneticAlgorithm::RunTest(int numberOfIterations) {
	std::vector<float> allBestPaths(numberOfIterations, -1);
	int cutoff = 150;

	for (int i = 0; i < numberOfIterations; ++i) {
		MakeNextGeneration();
		allBestPaths[i] = bestPath;

		if (i >= cutoff && abs(allBestPaths[i - cutoff] - bestPath) <= 60)
			break;
	}
}

float GeneticAlgorithm::RunFixedTime(double seconds) {
	time_t start, end;
	double elapsed;
	start = time(NULL);
	bool run = true;

	while (run) {
		end = time(NULL);
		elapsed = difftime(end, start);
		if (elapsed >= seconds) {
			run = false;
		}
		else {
			MakeNextGeneration();
		}
	}
	return bestPath;
}

void GeneticAlgorithm::PrintChromosome(std::vector<int>& chromosome) {
	for (int j = 0; j < chromosome.size(); ++j) {
		std::cout << chromosome[j] + 1 << ", ";
	}
	std::cout << chromosome[0] + 1 << "\n";
}

void GeneticAlgorithm::PrintBestChromosome() {
	if (bestWrittenFitness != bestFitness) {
		std::cout << "Best fitness: " << bestFitness << "\nBest Path: " << bestPath << "\nBest Chromosome: ";
		PrintChromosome(bestChromosome);
		std::cout << "\n";
		bestWrittenFitness = bestFitness;
	}
}

void GeneticAlgorithm::PrintPopulation() {
	for (int i = 0; i < populationSize; i++) {
		PrintChromosome(population[i].tour);
		std::cout << "\n";
	}
}
