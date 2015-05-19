#include"stdafx.h"

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

	This information should be removed after revising.
*/
int simulationPrepare(){
	static uniform_int_distribution<long> transDis_x(0, _Nx - 1);
	static uniform_int_distribution<long> transDis_y(0, _Ny - 1);
	static uniform_int_distribution<long> transDis_z(0, _Nz - 1);
	static uniform_int_distribution<> orntDis(0, 23);
	int x, y, z;
	for (int m = 0; m < N; m++){
		do{
			x = transDis_x(gen);
			y = transDis_y(gen);
			z = transDis_z(gen);
		} while (stage[x][y][z] >= 0);
		stage[x][y][z] = m;
		mol[m].put(x, y, z, orntDis(gen));
	}

	return 0;
}