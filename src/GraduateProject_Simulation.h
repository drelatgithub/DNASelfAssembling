#pragma once

#include <random>

constexpr double k_B = 1.381E-23; // SI units
constexpr double N_A = 6.023E23; // SI units

constexpr double cal2J = 4.184; // thermochemical calorie

extern std::uniform_real_distribution<> judge;
extern double T;
extern int *cluster_series;
extern int cluster_size;
int simulationPrepare();
int simulationProcess();

int showStats(int step, int totalSteps, int step_stat);
int maxCorrectSize();
