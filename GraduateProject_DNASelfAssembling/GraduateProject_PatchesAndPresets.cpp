#include"stdafx.h"

atom_in_unitcell atom_in_unitcell::set(int nx, int ny, int nz, int nsd, int nornt0, int nornt1){
	px = ppos(nx, ny, nz);
	sd = nsd;
	int mainOrnt = nornt0;
	int subOrnt = 0, subOrnt_as_pow = 7 - nornt0^nornt1;
	for (; subOrnt_as_pow > 1; subOrnt_as_pow >>= 1)subOrnt++;
	ornt = mainOrnt + 8 * subOrnt;
	return *this;
}
atom_in_unitcell unitcell[8];
int unitcell_init(){
	unitcell[0].set(0, 0, 0, 3, 7, 4);
	unitcell[1].set(0, 2, 2, 2, 7, 4);
	unitcell[2].set(2, 0, 2, 2, 4, 7);
	unitcell[3].set(2, 2, 0, 2, 4, 7);
	unitcell[4].set(1, 1, 1, 1, 5, 6);
	unitcell[5].set(1, 3, 3, 1, 5, 6);
	unitcell[6].set(3, 1, 3, 1, 6, 5);
	unitcell[7].set(3, 3, 1, 1, 6, 5);
	return 0;
}

DNAmol::DNAmol(){
	px.x = px.y = px.z = 0;
	ornt = 0;
	int i, j;
	for (i = 0; i < 4; i++){
		correctbond[i] = -1;
		for (j = 0; j < 8; j++)patch[i][j] = 3; // 'T' as unpaired nucleotides.
	}
}
DNAmol DNAmol::put(int nx, int ny, int nz, int nornt){
	px = ppos(nx, ny, nz);
	ornt = nornt;
	return *this;
}
DNAmol DNAmol::put(const atom_in_unitcell &a, int x_init, int y_init, int z_init){
	px = a.px + ppos(x_init, y_init, z_init);
	ornt = a.ornt;
	return *this;
}
int DNAmol::displayPatch()const{
	int i, j;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 8; j++){
			cout.put(nt[patch[i][j]]);
		}
	}
	cout << endl;
	return 0;
}
int DNAmol::displayPatch(ofstream &out)const{
	int i, j;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 8; j++){
			out.put(nt[patch[i][j]]);
		}
	}
	out << endl;
	return 0;
}
int stage[_Nx][_Ny][_Nz]; // To store DNAmol serial. -1 if not occupied. Using periodic boundary conditions.
DNAmol *mol;
int N; // Total molecules.

char nt[4] = { 'A', 'C', 'G', 'T' };
int ntSerial(char which_nt){
	switch (which_nt){
	case'A':return 0;
	case'C':return 1;
	case'G':return 2;
	case'T':return 3;
	default:cout << "Wrong nucleotide!" << endl; return -1;
	}
}
int ntSerialPair[4] = { 3, 2, 1, 0 };
int couldPatchInteract[4][4] = {
	{ 0, 1, 1, 0 },
	{ 1, 0, 0, 1 },
	{ 1, 0, 0, 1 },
	{ 0, 1, 1, 0 }
};
int isPatchConsequent[4] = { 1, 0, 0, 1 };

int stageClear(){
	int i, j, k;
	for (i = 0; i < _Nx; i++){
		for (j = 0; j < _Ny; j++){
			for (k = 0; k < _Nz; k++){
				stage[i][j][k] = -1;
			}
		}
	}
	return 0;
}
int molPreset(){
	int _N;
	N = 0;

	// 2x2x2 unitcells
	int i, j, k, m;
	_N = 2 * 2 * 2 * 8;
	mol = new DNAmol[_N];
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			for (k = 0; k < 2; k++){
				for (m = 0; m < 8; m++){
					mol[N].put(unitcell[m], i * 4 + 1, j * 4 + 1, k * 4 + 1);
					stage[mol[N].px.x][mol[N].px.y][mol[N].px.z] = N;
					N++;
				}
			}
		}
	}
	if (N != _N){
		cout << "The number of molecules is incorrect." << endl;
	}

	return 0;
}

int correctbonding(int n0serial, int n01x, int n01y, int n01z){
	ppos npx = mol[n0serial].px + ppos(n01x, n01y, n01z);
	int n1serial = stage[npx.x][npx.y][npx.z];
	int ornt0 = 4 * (n01x + 1) / 2 + 2 * (n01y + 1) / 2 + (n01z + 1) / 2;
	int ornt1 = 7 - ornt0;
	int n0ps = -1, n1ps = -1; // patch serial for the right bonding. -1 if not exist.
	int temp_ntSerial;
	n0ps = findPatchSerial[mol[n0serial].ornt][ornt0];
	n1ps = findPatchSerial[mol[n1serial].ornt][ornt1];
	if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps]){
		mol[n0serial].correctbond[n0ps] = n1serial;
		mol[n1serial].correctbond[n1ps] = n0serial;
		// patch assign
		static uniform_int_distribution<> patchDis(0, 3);
		for (int i = 0; i < 8; i++){
			temp_ntSerial = patchDis(gen);
			mol[n0serial].patch[n0ps][i] = temp_ntSerial;
			mol[n1serial].patch[n1ps][7 - i] = ntSerialPair[temp_ntSerial];
		}
		return 0;
	}
	return 1;
}
int nucleotidePreset(){
	int m, i, j, k;
	ppos npx;
	for (m = 0; m < N; m++){
		for (i = -1; i <= 1; i+=2){
			for (j = -1; j <= 1; j+=2){
				for (k = -1; k <= 1; k+=2){
					npx = mol[m].px + ppos(i, j, k);
					if (stage[npx.x][npx.y][npx.z]>m){
						correctbonding(m, i, j, k);
					}
				}
			}
		}
	}
	ofstream out("F:\\DNAPatches.txt");
	for (m = 0; m < N; m++){
		mol[m].displayPatch(out);
	}
	out.close();
	return 0;
}
int patchPreset(){
	stageClear();
	molPreset();
	nucleotidePreset();
	stageClear();
	return 0;
}