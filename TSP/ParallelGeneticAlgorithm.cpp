#include "ParallelGeneticAlgorithm.h"

ParallelGeneticAlgorithm::ParallelGeneticAlgorithm(int populationSize, int numberOfCities, float mutationRate, int numberOfParentPairs,
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

void ParallelGeneticAlgorithm::SetPopulationSize(int newSize) {
	populationSize = newSize;
}

int ParallelGeneticAlgorithm::GetPopulationSize() const {
	return populationSize;
}

void ParallelGeneticAlgorithm::SetMutationRate(float newRate) {
	mutationRate = newRate;
}

float ParallelGeneticAlgorithm::GetMutationRate() const {
	return mutationRate;
}

void ParallelGeneticAlgorithm::SetNumberOfCities(int newNumberOfCities) {
	numberOfCities = newNumberOfCities;
}

int ParallelGeneticAlgorithm::GetNumberOfCities() const {
	return numberOfCities;
}

int** ParallelGeneticAlgorithm::SetCitiesMatrix(int** matrix) {
	return matrix;
}

float ParallelGeneticAlgorithm::GetDistance(int a, int b) {
	return cities[a][b];
}

void ParallelGeneticAlgorithm::Run(int numberOfIterations, float bestFound) {
	//omp_set_nested(1);
	omp_set_num_threads(8);

	for (int i = 0; i < numberOfIterations; ++i) {
		int bestTour = 0;
		std::sort(population.begin(), population.end(), CompareFitness);

		#pragma omp parallel 
		{
			int bestParentIndex;
			int bestSecondParentIndex;
			int worstParentIndex;
			int worstSecondParentIndex;

			float localBestFitness = bestFitness;
			int localBestTour = -1;

			#pragma omp for
			for (int j = 0; j < numberOfParentPairs * 2; j += 2) {
				bestParentIndex = j;
				bestSecondParentIndex = j + 1;
				worstParentIndex = populationSize - j - 1;
				worstSecondParentIndex = populationSize - j - 2;

				population[worstParentIndex].tour = 
					Crossover(population[bestParentIndex].tour, 
						population[bestSecondParentIndex].tour);

				Mutate(population[worstParentIndex].tour);
				TwoOpt(population[worstParentIndex].tour);
				
				float path1 = 0.0f;
				int chrom1 = worstParentIndex;
				#pragma omp parallel 
				{
					float localPath = 0.0f;
					#pragma omp for
					for (int g = 0; g < numberOfCities - 1; g++) {
						localPath += cities[(population[chrom1].tour)[g]][(population[chrom1].tour)[g + 1]];
					}
					localPath += cities[(population[chrom1].tour)[numberOfCities - 1]][(population[chrom1].tour)[0]];

					#pragma omp atomic
					path1 += localPath;
				}

				population[chrom1].path = path1;
				population[chrom1].fitness = 1 / path1;

				population[worstSecondParentIndex].tour = 
					Crossover(population[bestSecondParentIndex].tour, 
						population[bestParentIndex].tour);

				Mutate(population[worstSecondParentIndex].tour);
				TwoOpt(population[worstSecondParentIndex].tour);
				
				float path2 = 0.0f;
				int chrom2 = worstSecondParentIndex;
				#pragma omp parallel 
				{
					float localPath = 0.0f;
					#pragma omp for
					for (int g = 0; g < numberOfCities - 1; g++) {
						localPath += cities[(population[chrom2].tour)[g]][(population[chrom2].tour)[g + 1]];
					}
					localPath += cities[(population[chrom2].tour)[numberOfCities - 1]][(population[chrom2].tour)[0]];

					#pragma omp atomic
					path2 += localPath;
				}

				population[chrom2].path = path2;
				population[chrom2].fitness = 1 / path2;

				if (population[worstParentIndex].fitness > localBestFitness) {
					localBestTour = worstParentIndex;
					localBestFitness = population[worstParentIndex].fitness;
				}

				if (population[worstSecondParentIndex].fitness > localBestFitness) {
					localBestTour = worstSecondParentIndex;
					localBestFitness = population[worstSecondParentIndex].fitness;
				}
			}

			#pragma omp critical
			{
				if (localBestTour != -1 && population[localBestTour].fitness > bestFitness) {
					bestFitness = population[localBestTour].fitness;
					bestPath = population[localBestTour].path;
					bestChromosome = population[localBestTour].tour;
					bestTour = localBestTour;
				}
			}
		}

		TwoOpt(population[bestTour].tour);
		if (population[bestTour].fitness > bestFitness) {
			bestFitness = population[bestTour].fitness;
			bestPath = population[bestTour].path;
			bestChromosome = population[bestTour].tour;
		}

		if (bestFound != -1.0f && bestFound >= bestPath)
			break;

		//PrintBestChromosome();
	}
}

void ParallelGeneticAlgorithm::ClearPopulation() {
	if (!population.empty()) {
		bestChromosome.clear();
		bestFitness = std::numeric_limits<float>::min();
		bestPath = std::numeric_limits<float>::max();
		bestWrittenFitness = std::numeric_limits<float>::min();
		population.clear();
	}
}

void ParallelGeneticAlgorithm::InitializePopulation() {
	#pragma omp parallel
	{
		int startCity;
		int nextCity;
		#pragma omp for
		for (int i = 0; i < populationSize; ++i) {
			std::vector<bool> visitedCities(numberOfCities, false);
			startCity = rand() % numberOfCities;
			nextCity = startCity;
			visitedCities[startCity] = true;

			for (int j = 0; j < numberOfCities - 1; ++j) {
				population[i].tour[j] = nextCity;
				if ((rand() % 100) < chanceToUseCloseCity) {
					nextCity = FindNearestNeighbour(population[i].tour[j], visitedCities);
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

			population[i].tour[numberOfCities - 1] = nextCity;
			CalculateFitness(i);

			#pragma omp critical 
			{
				if (population[i].fitness > bestFitness) {
					bestFitness = population[i].fitness;
					bestPath = population[i].path;
					bestChromosome = population[i].tour;
				}
			}
		}
	}
}

bool ParallelGeneticAlgorithm::CompareFitness(Chromosome& a, Chromosome& b) {
	return a.fitness > b.fitness;
}

int ParallelGeneticAlgorithm::FindNearestNeighbour(int city, std::vector<bool>& visitedCities) {
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

float ParallelGeneticAlgorithm::CalculateDistance(std::vector<int>& chromosome) {
	float path = 0;
	for (int i = 0; i < numberOfCities - 1; ++i) {
		path += cities[chromosome[i]][chromosome[i + 1]];
	}

	path += cities[chromosome[numberOfCities - 1]][chromosome[0]];
	return path;
}

void ParallelGeneticAlgorithm::CalculateFitness(int index) {
	float path;
	path = CalculateDistance(population[index].tour);
	population[index].path = path;
	population[index].fitness = 1 / path;
}

void ParallelGeneticAlgorithm::TwoOpt(std::vector<int>& chromosome) {
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

void ParallelGeneticAlgorithm::Mutate(std::vector<int>& chromosome) {
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

std::vector<int> ParallelGeneticAlgorithm::Crossover(std::vector<int>& parent1, std::vector<int>& parent2) {
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

int ParallelGeneticAlgorithm::FindCityIndex(std::vector<int>& chromosome, int city, int start, int end) {
	for (int i = start; i <= end; ++i) {
		if (chromosome[i] == city) {
			return i;
		}
	}
	return -1;
}

void ParallelGeneticAlgorithm::PrintChromosome(std::vector<int>& chromosome) {
	for (int j = 0; j < chromosome.size(); ++j) {
		std::cout << chromosome[j] + 1 << " ";
	}
	std::cout << chromosome[0] + 1 << "\n";
}

void ParallelGeneticAlgorithm::PrintBestChromosome() {
	if (bestWrittenFitness != bestFitness) {
		std::cout << "Best fitness: " << bestFitness << "\nBest Path: " << bestPath << "\nBest Chromosome: ";
		PrintChromosome(bestChromosome);
		std::cout << "\n";
		bestWrittenFitness = bestFitness;
	}
}

void ParallelGeneticAlgorithm::PrintPopulation() {
	for (int i = 0; i < populationSize; i++) {
		PrintChromosome(population[i].tour);
		std::cout << "\n";
	}
}

bool ParallelGeneticAlgorithm::IsPathValid(std::vector<int>& path) {
	std::vector<int> visitedCities(path.size(), 0);
	for (int j = 0; j < path.size(); j++) {
		visitedCities[path[j]]++;
	}

	bool valid = true;
	for (int j = 0; j < path.size(); j++) {
		if (visitedCities[j] != 1) {
			std::cout << "Invalid city!\n";
			valid = false;
			break;
		}
	}

	if (!valid)
		std::cout << "Invalid chromosome!\n";

	return valid;
}