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
#include "Greedy.h"
#include "Chromosome.h"

class GeneticAlgorithm {
public:
	float bestFitness;
	float bestPath;
	std::vector<int> bestChromosome;

	GeneticAlgorithm(int populationSize, int numberOfCities, float mutationRate, int numberOfParentPairs,
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
	void Run(int numberOfIterations);
	void RunTest(int numberOfIterations);
	float RunFixedTime(double seconds);
	static bool CompareFitness(Chromosome& a, Chromosome& b);

private:
	int numberOfParentPairs;
	int populationSize;
	int numberOfCities;
	int mutationRate;
	int chanceToUseCloseCity;
	int twoOptIterations;
	float** cities;
	float bestWrittenFitness;
	std::vector<Chromosome> population;

	void ClearPopulation();
	void CalculateFitness(int index);
	float CalculateDistance(std::vector<int>& chromosome);
	int FindCityIndex(std::vector<int>& chromosome, int city, int start, int end);
	std::vector<int> Crossover(std::vector<int>& parent1, std::vector<int>& parent2);
	int FindNearestNeighbour(int city, std::vector<bool>& visitedCities);
	void Mutate(std::vector<int>& chromosome);
	void TwoOpt(std::vector<int>& chromosome);
	void MakeNextGeneration();
	float GetDistance(int a, int b);
	void PrintBestChromosome();
	void PrintChromosome(std::vector<int>& chromosome);
	void GenerateTour(int* i);
};