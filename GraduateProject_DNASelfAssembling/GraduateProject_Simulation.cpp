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
/*
	The algorithm used in this function has LOW efficiency, and needs to be improved.
	
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
int simulationPrepare(){
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

#define k_B 1.381E-23
#define N_A 6.023E23


int moveStep(int s){

	return 0;
}
int simulationProcess(){
	int totalSteps = 10000;
	int step = 0;

	return 0;
}