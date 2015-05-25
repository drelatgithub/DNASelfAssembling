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
bool couldPlaceMolecule(ppos px){
	return couldPlaceMolecule(px.x, px.y, px.z);
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

double T; // Kelvin

double energy_local_patch(int s, int n0ps, int s1, int n1ps){
	double E_patch_kcal = 0;
	int mismatches = 0, internalMismatches = 0;
	int m;
	int a, b, c;
	for (m = 0; m < 8; m++){
		if (ntSerialPair[mol[s].patch[n0ps][m]] != mol[s1].patch[n1ps][7 - m]){
			mismatches++;
			if (m > 0 && m < 7)internalMismatches++;
		}
	}
	int s53 = isPatchConsequent[n0ps] ? s : s1, ps53 = isPatchConsequent[n0ps] ? n0ps : n1ps;
	int s35 = isPatchConsequent[n0ps] ? s1 : s, ps35 = isPatchConsequent[n0ps] ? n1ps : n0ps;
	switch (mismatches){
	case 0:
		for (m = 0; m < 7; m++){
			a = mol[s53].patch[ps53][m];
			b = mol[s53].patch[ps53][m + 1];
			c = mol[s35].patch[ps35][6 - m];
			E_patch_kcal += NN_dH_kcal[a][b][c] - T*NN_dS_cal[a][b][c] / 1000;
		}
		break;
	case 1:
		if (internalMismatches == 1){
			for (m = 0; m < 7; m++){
				a = mol[s53].patch[ps53][m];
				b = mol[s53].patch[ps53][m + 1];
				c = mol[s35].patch[ps35][6 - m];
				E_patch_kcal += NN_dH_kcal[a][b][c] - T*NN_dS_cal[a][b][c] / 1000;
				if (ntSerialPair[b] != c){
					a = mol[s35].patch[ps35][5 - m];
					b = mol[s35].patch[ps35][6 - m];
					c = mol[s53].patch[ps53][m + 1];
					E_patch_kcal += NN_dH_kcal[a][b][c] - T*NN_dS_cal[a][b][c] / 1000;
					m++;
				}
			}
		}
	}
	return E_patch_kcal;
}
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
	static ppos npx;

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
	1. For a given number of basepairs in a certain duplex, and a given salt concentration,
	   the free energy corrections by sodium dependence is constant, thus ignored when calculating energy differences.
	2. The initiation energy is not applicable because the DNA brick would not allow for a duplex to start together.
	3. The terminal A-T penalty is not applicable because the DNA brick would not allow for a duplex to end together.
	*/
	double E_patches_kcal = 0; // 1000 calories per mol
	for (i = -1; i <= 1; i += 2){
		for (j = -1; j <= 1; j += 2){
			for (k = -1; k <= 1; k += 2){
				npx = mol[s].px + ppos(i, j, k);
				int s1 = stage[npx.x][npx.y][npx.z];
				if (s1 >= 0){
					int ornt0 = 4 * (i + 1) / 2 + 2 * (j + 1) / 2 + (k + 1) / 2;
					int ornt1 = 7 - ornt0;
					int n0ps = -1, n1ps = -1; // patch serial for the relative bonding. -1 if not exist.
					n0ps = findPatchSerial[mol[s].ornt][ornt0];
					n1ps = findPatchSerial[mol[s1].ornt][ornt1];
					if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps]){
						E_patches_kcal += energy_local_patch(s, n0ps, s1, n1ps);
					}
				}
			}
		}
	}

	double E_total = E_repulsive + E_patches_kcal * 1000 * cal2J / N_A;
	return E_total;
}
int moveStep(int s){
	double E0 = energy_local(s);
	static uniform_int_distribution<> translateOrRotate(0, 1);
	static uniform_int_distribution<> anotherOrnt(0, 22); // only 23 possible orientations in order to exclude the original orientation
	static uniform_int_distribution<> anotherCoor(-1, 1);
	static uniform_real_distribution<> judge(0, 1); // Metropolis Criterion
	static ppos opx, npx;
	if (translateOrRotate(gen)){
		int oOrnt = mol[s].ornt;
		int nOrnt = anotherOrnt(gen);
		nOrnt = (nOrnt >= oOrnt) ? (nOrnt + 1) : nOrnt;
		mol[s].ornt = nOrnt;
		double E1 = energy_local(s);
		if (E1 > E0 && judge(gen) > exp(-(E1 - E0) / k_B / T)){
			mol[s].ornt = oOrnt;
		}
	}
	else{
		opx = mol[s].px;
		npx = opx + ppos(anotherCoor(gen), anotherCoor(gen), anotherCoor(gen));
		if (stage[npx.x][npx.y][npx.z] == -1){
			mol[s].px = npx;
			stage[opx.x][opx.y][opx.z] = -1;
			stage[npx.x][npx.y][npx.z] = s;
			double E1 = energy_local(s);
			if (E1 > E0 && judge(gen) > exp(-(E1 - E0) / k_B / T)){
				mol[s].px = opx;
				stage[opx.x][opx.y][opx.z] = s;
				stage[npx.x][npx.y][npx.z] = -1;
			}
		}
	}
	return 0;
}
int simulationProcess(){
	long totalSteps = 10000000;
	long step;
	long step_stat = 5000;
	T = 310;

	for (step = 1; step <= totalSteps; step++){
		if (step % step_stat == 0){
			showStats(step, totalSteps, step_stat);
		}
		for (int i = 0; i < N; i++){
			moveStep(i);
		}
	}

	return 0;
}

int timeDisplay(double time_sec){
	int time_min = (int)(time_sec / 60);
	time_sec -= time_min * 60;
	int time_hour = time_min / 60;
	time_min -= time_hour * 60;
	int time_day = time_hour / 24;
	time_hour -= time_day * 24;
	if (time_day != 0){
		cout << time_day << "d";
	}
	if (time_hour != 0){
		cout << ' ' << time_hour << "h";
	}
	if (time_min != 0){
		cout << ' ' << time_min << "m";
	}
	cout << ' ' << (int)time_sec << "s" << endl;
	return 0;
}
int showStats(int step, int totalSteps, int step_stat){
	t_end = clock();
	double time_used = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	double time_per_step = time_used / step;
	double time_remaining = time_per_step*(totalSteps - step);
	int maxSize = maxCorrectSize();
	static int historyMax = 0;
	if (historyMax < maxSize)historyMax = maxSize;
	cout << endl << step << '/' << totalSteps << " steps" << endl;
	cout << "max size: " << maxSize << endl;
	cout << "time used: "; timeDisplay(time_used);
	cout << "time per " << step_stat << " steps: " << time_per_step * step_stat << endl;
	cout << "estimated remaining time: "; timeDisplay(time_remaining);
	cout << endl;
	ofstream coorOut("F:\coordinates.txt");
	coorOut << T << endl << step << endl << maxSize << endl << historyMax << endl;
	for (int i = 0; i < N; i++){
		coorOut << mol[i].px.x << ' ' << mol[i].px.y << ' ' << mol[i].px.z << endl;
	}
	return 0;
}

int anotherMoleculeCombined(int *mark, int previousSerial){
	int ok = 1;
	mark[previousSerial] = 1;
	for (int i = 0; i < 4; i++){
		if (mol[previousSerial].correctbond[i] >= 0 && !mark[mol[previousSerial].correctbond[i]]){
			ppos dpx = mol[mol[previousSerial].correctbond[i]].px - mol[previousSerial].px;
			if ((dpx.x == 1 || dpx.x == _Nx - 1) && (dpx.y == 1 || dpx.y == _Ny - 1) && (dpx.z == 1 || dpx.z == _Nz - 1)){
				ok += anotherMoleculeCombined(mark, mol[previousSerial].correctbond[i]);
			}
		}
	}
	return ok;
}
int maxCorrectSize(){
	int *mark = new int[N];
	int max = 0, temp;
	int p;
	for (p = 0; p < N; p++){
		mark[p] = 0;
	}
	p = 0;
	while (p < N){
		if (!mark[p]){
			temp = anotherMoleculeCombined(mark, p);
			cout << temp << '\t';
			if (max < temp)max = temp;
		}
		p++;
	}
	cout << endl;

	delete[]mark;
	return max;
}