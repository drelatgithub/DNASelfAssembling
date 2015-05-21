#include"stdafx.h"

bool couldPlaceMolecule(int x, int y, int z){
	int i, j, k;
	ppos npx;
	for (i = -1; i <= 1; i++){
		for (j = -1; j <= 1; j++){
			for (k = -1; k <= 1; k++){
				npx = ppos(x + i, y + j, z + k);
				if (i*i + j*j + k*k != 3 && stage[npx.x][npx.y][npx.z] >= 0){
					return false;
				}
			}
		}
	}
	return true;
}
int simulationPrepare(){
	/*
		The following algorithm used has LOW efficiency, and needs to be improved.

		We put the DNA molecules in the lattice randomly in the following manner:
		1: generate a random number as the position in the lattice;
		2: if the position is not taken, put the molecule in it, otherwise repeat 1;
		3: repeat 1 - 2 for all molecules.

		The ratio of useful to useless repetitions could be extremely LOW when the lattice becomes crowded.
		Assume the lattice with n different positions in which m is already taken,
		then the expectation value for the repetitions required to put a new molecule is n/(n-m).
		The overall time complexity to put k molecules is O(n ln(n/(n-k)))
		where the best case is O(k), and the worst case is O(\infty).

		As for this case, every DNA molecule would in average consume 4 of the lattice positions,
		thus making the stage easily get crowded.

		====== This information should be removed after revising. ======
	*/
	uniform_int_distribution<long> transDis_x(0, _Nx - 1);
	uniform_int_distribution<long> transDis_y(0, _Ny - 1);
	uniform_int_distribution<long> transDis_z(0, _Nz - 1);
	uniform_int_distribution<> orntDis(0, 23);
	int x, y, z;
	for (int m = 0; m < N; m++){
		do{
			x = transDis_x(gen);
			y = transDis_y(gen);
			z = transDis_z(gen);
		} while (!couldPlaceMolecule(x, y, z));
		stage[x][y][z] = m;
		mol[m].put(x, y, z, orntDis(gen));
	}

	return 0;
}

const double delta_H_NN_init = 0.2E3;
const double delta_S_NN_init = -5.7;
double T; // Kelvin

double energy_local(int s){
	/*
	This energy includes the repulsive energy (100K) and the hybridization free energy for the matched patches,
	allowing for single internal mismatches.

	The free energy is determined using the nearest-neighbor parametrization (pH 7),
	taking into account terminal A-T penalties,
	                    internal mismatches,
						dangling ends, and
						temperature and salt concentration dependence,
	but do not consider loops or bulges.

	References:
	1. A. Reinhardt and D. Frenkel, Numerical Evidence for Nucleated Self-Assembly of DNA Brick Structures, Phys. Rev. Lett. 112, 238103 (2014)
	2. J. SantaLucia, Jr. and D. Hicks, The Thermodynamics of DNA Structural Motifs, Annu. Rev. Biophys. Biomol. Struct. 33, 415 (2004)
	*/
	int i, j, k, nearests;
	ppos npx;

	double E_repulsive; // Joules
	nearests = 0;
	for (i = -1; i <= 1; i += 2){
		for (j = -1; j <= 1; j += 2){
			for (k = -1; k <= 1; k += 2){
				npx = mol[s].px + ppos(i, j, k);
				if (stage[npx.x][npx.y][npx.z] >= 0)nearests++;
			}
		}
	}
	E_repulsive = nearests*k_B*100;

	/*
	For a given number of basepairs in a certain duplex, and a given salt concentration,
	the free energy corrections by sodium dependence is constant, thus ignored when calculating energy differences.
	*/
	double E_patches_cal = 0;
	for (i = -1; i <= 1; i += 2){
		for (j = -1; j <= 1; j += 2){
			for (k = -1; k <= 1; k += 2){
				npx = mol[s].px + ppos(i, j, k);
				int s1 = stage[npx.x][npx.y][npx.z];
				if (s1 >= 0){
					int ornt0 = 4 * (i + 1) / 2 + 2 * (j + 1) / 2 + (k + 1) / 2;
					int ornt1 = 7 - ornt0;
					int n0ps = -1, n1ps = -1; // patch serial for the relative bonding. -1 if not exist.
					n0ps = mol[s].findPatchSerial(ornt0);
					n1ps = mol[s].findPatchSerial(ornt1);
					if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps]){
						int mismatches = 0, internalMismatches = 0;
						int m;
						for (m = 0; m < 8; m++){
							if (ntPair(mol[s].patch[n0ps][m]) != mol[s1].patch[n1ps][7 - m]){
								mismatches++;
								if (m > 0 && m < 7)internalMismatches++;
							}
						}
						switch (mismatches){
						case 0:
							E_patches_cal += delta_H_NN_init - T*delta_S_NN_init;
							for (m = 0; m < 7; m++){

							}
							break;
						case 1:
							if (internalMismatches == 1){

							}
						}
					}
				}
			}
		}
	}

	cout << E_repulsive << endl;
	

	return 0;
}
int moveStep(int s){
	double E0 = energy_local(s);
	return 0;
}
int simulationProcess(){
	int totalSteps = 10000;
	int step = 0;
	T = 300;

	for (int i = 0; i < N; i++){
		cout << i << '\t';
		moveStep(i);
	}

	return 0;
}