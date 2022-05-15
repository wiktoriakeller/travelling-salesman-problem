#pragma once
#include <vector>
#include <algorithm>   
#include <iostream>
#include <ctime>
#include <random>
#include <limits>
#include <map>
#include <functional>
#include <mutex>
#include <thread>
#include <omp.h>
#include "Greedy.h"
#include "Chromosome.h"

class ParallelGeneticAlgorithm {
public:
	ParallelGeneticAlgorithm(int populationSize, int numberOfCities, float mutationRate, int numberOfParentPairs,
		int chanceToUseCloseCity, int twoOptIterations, float** cities);
	void SetPopulationSize(int newSize);
	int GetPopulationSize() const;
	void SetMutationRate(float newRate);
	float GetMutationRate() const;
	void SetNumberOfCities(int newNumberOfCities);
	int GetNumberOfCities() const;
	int** SetCitiesMatrix(int** matrix);
	void InitializePopulation();
	void PrintPopulation();
	float Run(int numberOfIterations);
	static bool CompareFitness(Chromosome& a, Chromosome& b);
	std::vector<int> GetBestChromosome();

private:
	int numberOfParentPairs;
	int populationSize;
	int numberOfCities;
	int mutationRate;
	int chanceToUseCloseCity;
	int twoOptIterations;
	float** cities;
	float bestFitness;
	float bestPath;
	float bestWrittenFitness;
	std::vector<Chromosome> population;
	std::vector<int> bestChromosome;

	void ClearPopulation();
	void CalculateFitness(int index);
	float CalculateDistance(std::vector<int>& chromosome);
	int FindCityIndex(std::vector<int>& chromosome, int city, int start, int end);
	std::vector<int> Crossover(std::vector<int>& parent1, std::vector<int>& parent2);
	int FindNearestNeighbour(int city, std::vector<bool>& visitedCities);
	void Mutate(std::vector<int>& chromosome);
	void TwoOpt(std::vector<int>& chromosome);
	float GetDistance(int a, int b);
	void PrintBestChromosome();
	void PrintChromosome(std::vector<int>& chromosome);
	bool IsPathValid(std::vector<int>& chromosome);
};