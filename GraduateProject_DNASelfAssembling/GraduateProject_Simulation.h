#pragma once

#define k_B 1.381E-23 // SI units
#define N_A 6.023E23 // SI units

#define cal2J 4.184 // thermochemical calorie

extern uniform_real_distribution<> judge;
extern double T;
int simulationPrepare();
int simulationProcess();

int showStats(int step, int totalSteps, int step_stat);
int maxCorrectSize();
