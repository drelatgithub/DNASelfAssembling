#include"stdafx.h"

bool couldPlaceMolecule(int x, int y, int z){
	int i, j, k;
	static ppos npx;
	for (i = -1; i <= 1; i++){
		for (j = -1; j <= 1; j++){
			for (k = -1; k <= 1; k++){
				npx.set(x + i, y + j, z + k);
				npx.adjust();
				if (i*i + j*j + k*k != 3 && stage[npx.x][npx.y][npx.z] >= 0){
					return false;
				}
			}
		}
	}
	return true;
}
bool couldPlaceMolecule(const ppos &px){
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
	int m;
	for (m = 0; m < N; m++){
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
uniform_real_distribution<> judge(0, 1);

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
	for (i = -1; i <= 1; i++){
		for (j = -1; j <= 1; j++){
			for (k = -1; k <= 1; k++){
				if (i || j || k){
					npx = mol[s].px.plus(i, j, k);
					npx.adjust();
					if (stage[npx.x][npx.y][npx.z] >= 0)nearests++;
				}
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
				npx = mol[s].px.plus(i, j, k);
				npx.adjust();
				int s1 = stage[npx.x][npx.y][npx.z];
				if (s1 >= 0){
					int ornt0 = 4 * (i + 1) / 2 + 2 * (j + 1) / 2 + (k + 1) / 2;
					int ornt1 = 7 - ornt0;
					int n0ps = -1, n1ps = -1; // patch serial for the relative bonding. -1 if not exist.
					n0ps = findPatchSerial[mol[s].ornt][ornt0];
					n1ps = findPatchSerial[mol[s1].ornt][ornt1];
					if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps] && couldPatchInteract_ornt(mol[s].ornt,n0ps,mol[s1].ornt,n1ps)){
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
	static uniform_int_distribution<> anotherOrntAxis(0, 2);
	static uniform_int_distribution<> anotherOrntTurns(1, 3);
	static uniform_int_distribution<> anotherCoor(-1, 1);
	static ppos opx, npx;
	if (translateOrRotate(gen)){
		int oOrnt = mol[s].ornt;
		int nOrntAxis = anotherOrntAxis(gen);
		int nOrntTurns = anotherOrntTurns(gen);
		int nOrnt = orntRot[oOrnt][nOrntAxis][nOrntTurns];
		mol[s].ornt = nOrnt;
		double E1 = energy_local(s);
		if (E1 > E0 && judge(gen) > exp(-(E1 - E0) / k_B / T)){
			mol[s].ornt = oOrnt;
		}
	}
	else{
		opx = mol[s].px;
		npx = opx.plus(anotherCoor(gen), anotherCoor(gen), anotherCoor(gen));
		npx.adjust();
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
ppos posAfterMove(int dpxx, int dpxy, int dpxz, int cs, int axis, int turns, int dx, int dy, int dz){
	int ndpxx = dpxx, ndpxy = dpxy, ndpxz = dpxz;
	if (turns){
		switch (axis){
		case 0:
			switch (turns){
			case 1:ndpxy = -dpxz; ndpxz = dpxy; break;
			case 2:ndpxy = -dpxy; ndpxz = -dpxz; break;
			case 3:ndpxy = dpxz; ndpxz = -dpxy;
			}break;
		case 1:
			switch (turns){
			case 1:ndpxz = -dpxx; ndpxx = dpxz; break;
			case 2:ndpxz = -dpxz; ndpxx = -dpxx; break;
			case 3:ndpxz = dpxx; ndpxx = -dpxz;
			}break;
		case 2:
			switch (turns){
			case 1:ndpxx = -dpxy; ndpxy = dpxx; break;
			case 2:ndpxx = -dpxx; ndpxy = -dpxy; break;
			case 3:ndpxx = dpxy; ndpxy = -dpxx;
			}
		}
	}
	return mol[cs].px.plus(ndpxx + dx, ndpxy + dy, ndpxz + dz);
}
double interactEnergy(int s0, int s1){ // never use stage[][][] in this function
	static ppos dpx;
	double E = 100 * k_B;
	dpx = mol[s1].px - mol[s0].px;
	dpx.adjust();
	if (dpx.x > 1)dpx.x -= _Nx;
	if (dpx.y > 1)dpx.y -= _Ny;
	if (dpx.z > 1)dpx.z -= _Nz;
	if (dpx.x*dpx.x + dpx.y*dpx.y + dpx.z*dpx.z > 3){
		E = 0.0;
		return E;
	}
	if (dpx.x && dpx.y && dpx.z){
		int ornt0 = 4 * (dpx.x + 1) / 2 + 2 * (dpx.y + 1) / 2 + (dpx.z + 1) / 2;
		int ornt1 = 7 - ornt0;
		int n0ps = -1, n1ps = -1; // patch serial for the relative bonding. -1 if not exist.
		n0ps = findPatchSerial[mol[s0].ornt][ornt0];
		n1ps = findPatchSerial[mol[s1].ornt][ornt1];
		if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps] && couldPatchInteract_ornt(mol[s0].ornt, n0ps, mol[s1].ornt, n1ps)){
			E += energy_local_patch(s0, n0ps, s1, n1ps) * 1000 * cal2J / N_A;
		}
	}
	return E;
}
bool recruit(int s, int dpxx, int dpxy, int dpxz, int cs, int axis, int turns, int dx, int dy, int dz){
	/* By introducing seemingly redundant parameter dpxx, dpxy and dpxz
	   (which is simply the position difference between molecule s and cs), 
	   we could be sure that the position difference would always have a change between -1 and 1,
	   thus avoiding the unphysical rotation near the surface due to periodic boundary conditions.

	   note: do not use static variables in this function unless required, because the function is recursive.
	*/
	mol[s].status = 2;
	int i, j, k;
	int s1;
	int oOrnt, nOrnt;
	double E0, Ef, Er;
	double pf, pr;
	ppos px_cnt = mol[s].px,
		 px_fwd = posAfterMove(dpxx, dpxy, dpxz, cs, axis, turns, dx, dy, dz),
		 px_rev = posAfterMove(dpxx, dpxy, dpxz, cs, axis, turns ? (4 - turns) : 0, -dx, -dy, -dz),
		 npx;
	px_fwd.adjust();
	px_rev.adjust();
	mol[s].px_cnt_backup = px_cnt;
	mol[s].px_fwd_dest = px_fwd;
	mol[s].ornt_cnt_backup = oOrnt = mol[s].ornt;
	mol[s].ornt_fwd_dest = nOrnt = orntRot[mol[s].ornt][axis][turns];
	bool accept = true;
	for (i = -1; i <= 1; i++){
		for (j = -1; j <= 1; j++){
			for (k = -1; k <= 1; k++){
				npx = mol[s].px.plus(i, j, k);
				npx.adjust();
				if ((i || j || k) && stage[npx.x][npx.y][npx.z] >= 0){ // propose a link
					s1 = stage[npx.x][npx.y][npx.z];
					if (px_fwd == npx){ // forward move cause overlap
						switch (mol[s1].status){
						case 0: // recruit and accept (temp)
							if (mol[s].links >= 8){
								cout << "Warning: too many links for particle " << s << endl;
								break;
							}
							mol[s].link_with[mol[s].links] = s1;
							if (px_rev == npx)mol[s].link_with_pRatio[mol[s].links] = 1; // reverse also overlap
							else{ // reverse not overlap
								E0 = interactEnergy(s, s1);
								mol[s].put(px_rev.x, px_rev.y, px_rev.z, orntRot[oOrnt][axis][turns ? (4 - turns) : 0]);
								Er = interactEnergy(s, s1);
								mol[s].put(px_cnt.x, px_cnt.y, px_cnt.z, oOrnt);
								pr = (Er > E0) ? (1 - exp(-(Er - E0) / k_B / T)) : 0;
								mol[s].link_with_pRatio[mol[s].links] = pr / 1.0;
							}
							mol[s].links++;
							accept = accept && recruit(s1, dpxx + i, dpxy + j, dpxz + k, cs, axis, turns, dx, dy, dz);
							break;
						case 1: // unrecruit and deny
							accept = false;
							break;
						case 2: ;// do nothing because already recruited
						}
					}
					else{ // forward move do not cause overlap
						switch (mol[s1].status){
						case 0: // try to recruit
							E0 = interactEnergy(s, s1);
							mol[s].put(px_fwd.x, px_fwd.y, px_fwd.z, nOrnt);
							Ef = interactEnergy(s, s1);
							mol[s].put(px_cnt.x, px_cnt.y, px_cnt.z, oOrnt);
							pf = (Ef > E0) ? (1 - exp(-(Ef - E0) / k_B / T)) : 0;
							if (judge(gen) < pf){ // link formed (where pf must be over 0.0)
								if (mol[s].links >= 8){
									cout << "Warning: too many links for particle " << s << endl;
									break;
								}
								mol[s].link_with[mol[s].links] = s1;
								if (px_rev == npx)mol[s].link_with_pRatio[mol[s].links] = 1; // reverse overlap
								else{ // reverse not overlap
									mol[s].put(px_rev.x, px_rev.y, px_rev.z, orntRot[oOrnt][axis][turns ? (4 - turns) : 0]);
									Er = interactEnergy(s, s1);
									mol[s].put(px_cnt.x, px_cnt.y, px_cnt.z, oOrnt);
									pr = (Er > E0) ? (1 - exp(-(Er - E0) / k_B / T)) : 0;
									mol[s].link_with_pRatio[mol[s].links] = pr / pf;
								}
								mol[s].links++;
								accept = accept && recruit(s1, dpxx + i, dpxy + j, dpxz + k, cs, axis, turns, dx, dy, dz);
							}
							break;
						case 1: case 2: ;// do nothing
						}
					}
					//calculate energy
				}
			}
		}
	}
	// do not prevent far away collisions
	return accept;
}
double groupInteractEnergy(){
	int m, i, j, k, s;
	static ppos npx;
	double E = 0;
	for (m = 0; m < N; m++){
		if (mol[m].status == 2){
			for (i = -1; i <= 1; i++){
				for (j = -1; j <= 1; j++){
					for (k = -1; k <= 1; k++){
						if (i || j || k){
							npx = mol[m].px.plus(i, j, k);
							npx.adjust();
							s = stage[npx.x][npx.y][npx.z];
							if (s >= 0){
								if (mol[s].status != 2){
									E += interactEnergy(m, s);
								}
							}
						}
					}
				}
			}
		}
	}
	return E;
}
int groupMove(bool forward){
	int m;
	for (m = 0; m < N; m++){
		if (mol[m].status == 2){
			stage[mol[m].px.x][mol[m].px.y][mol[m].px.z] = -1;
		}
	}
	for (m = 0; m < N; m++){
		if (mol[m].status == 2){
			if (forward){
				mol[m].px = mol[m].px_fwd_dest;
				mol[m].ornt = mol[m].ornt_fwd_dest;
			}
			else{
				mol[m].px = mol[m].px_cnt_backup;
				mol[m].ornt = mol[m].ornt_cnt_backup;
			}
			stage[mol[m].px.x][mol[m].px.y][mol[m].px.z] = m;
		}
	}
	return 0;
}
int moveStep_Group(int s){
	int m, i;
	static uniform_int_distribution<> translateOrRotate(0, 1);
	static uniform_int_distribution<> anotherOrntAxis(0, 2);
	static uniform_int_distribution<> anotherOrntTurns(1, 3);
	static uniform_int_distribution<> anotherCoor(-1, 1);
	int nOrntAxis = 0, nOrntTurns = 0;
	static ppos dpx;
	int s1;
	double E0, E1, p_acc;
	bool b_acc = true;
	dpx.set(0, 0, 0);
	if (translateOrRotate(gen)){
		nOrntAxis = anotherOrntAxis(gen);
		nOrntTurns = anotherOrntTurns(gen);
	}
	else{
		dpx.set(anotherCoor(gen), anotherCoor(gen), anotherCoor(gen)); // without adjustments
	}
	if (nOrntTurns || dpx.x || dpx.y || dpx.z){ // at least one move
		if (recruit(s, 0, 0, 0, s, nOrntAxis, nOrntTurns, dpx.x, dpx.y, dpx.z)){
			for (m = 0; m < N; m++){
				if (mol[m].status == 2){
					s1 = stage[mol[m].px_fwd_dest.x][mol[m].px_fwd_dest.y][mol[m].px_fwd_dest.z];
					if (s1 >= 0 && mol[s1].status != 2){
						b_acc = false; // prevent overlapping
					}
				}
			}
			if (b_acc){
				E0 = groupInteractEnergy();
				groupMove(true);
				E1 = groupInteractEnergy();
				p_acc = (E1 > E0) ? exp(-(E1 - E0) / k_B / T) : 1;
				for (m = 0; m < N; m++){
					if (mol[m].status == 2){
						for (i = 0; i < mol[m].links; i++){
							p_acc *= mol[m].link_with_pRatio[i];
						}
					}
				}
				if (judge(gen) < p_acc); // accept the move
				else{ // return to original position
					groupMove(false);
				}
			}
		}
		for (m = 0; m < N; m++){
			if (mol[m].status == 2){
				mol[m].status = 1; // mark as used
			}
		}
	} // (else) no move proposed, and leave the molecule marked 0 (unused), to facilitate future moves if recruited
	return 0;
}
int simulationProcess(){
	long totalSteps = 10000000;
	long step;
	long step_stat = 400;
	T = 315;

	for (step = 1; step <= totalSteps; step++){
		for (int i = 0; i < N; i++){
			mol[i].status = 0; // marked as unused
			mol[i].links = 0; // unlinked
		}
		if (step % step_stat == 0){
			showStats(step, totalSteps, step_stat);
		}
		for (int i = 0; i < N; i++){
			//moveStep(i);
			if (mol[i].status == 0){
				moveStep_Group(i);
			}
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
	ofstream coorOut("F:\\coordinates.txt");
	coorOut << T << endl << step << endl << maxSize << endl << historyMax << endl;
	for (int i = 0; i < N; i++){
		coorOut << mol[i].px.x << ' ' << mol[i].px.y << ' ' << mol[i].px.z << endl;
	}
	return 0;
}

int anotherMoleculeCombined(int *mark, int previousSerial){
	int ok = 1;
	mark[previousSerial] = 1;
	ppos dpx;
	for (int i = 0; i < 4; i++){
		if (mol[previousSerial].correctbond[i] >= 0 && !mark[mol[previousSerial].correctbond[i]]){
			dpx = mol[mol[previousSerial].correctbond[i]].px - mol[previousSerial].px;
			dpx.adjust();
			if ((dpx.x == 1 || dpx.x == _Nx - 1) && (dpx.y == 1 || dpx.y == _Ny - 1) && (dpx.z == 1 || dpx.z == _Nz - 1)){
				ok += anotherMoleculeCombined(mark, mol[previousSerial].correctbond[i]);
			}
		}
	}
	return ok;
}
int maxCorrectSize(){
	int *mark = new int[N];
	int *clustersize = new int[N];
	int clusternum = 0;
	int sum_clustersize = 0;
	double avg_clustersize, stdevp_clustersize = 0.0;
	int max = 0, temp;
	int p;
	for (p = 0; p < N; p++){
		mark[p] = 0;
	}
	p = 0;
	while (p < N){
		if (!mark[p]){
			temp = anotherMoleculeCombined(mark, p);
			clustersize[clusternum++] = temp;
			if (max < temp)max = temp;
		}
		p++;
	}
	for (p = 0; p < clusternum; p++){
		sum_clustersize += clustersize[p];
	}
	avg_clustersize = (double)sum_clustersize / clusternum;
	for (p = 0; p < clusternum; p++){
		stdevp_clustersize += (clustersize[p] - avg_clustersize)*(clustersize[p] - avg_clustersize);
	}
	stdevp_clustersize = sqrt(stdevp_clustersize / clusternum);
	cout << "clusters: " << clusternum << endl;
	cout << "avg size: " << avg_clustersize << endl;
	cout << "std error: " << stdevp_clustersize << endl;

	delete[]mark;
	delete[]clustersize;
	return max;
}