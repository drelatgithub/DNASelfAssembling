#include "GraduateProject_PatchesAndPresets.h"

#include <iostream>
#include <fstream>

#include "globals.h"
#include "GraduateProject_PositionAndOrientation.h"

atom_in_unitcell& atom_in_unitcell::set(int nx, int ny, int nz, int nsd, int nornt0, int nornt1){
	px.set(nx, ny, nz);
	sd = nsd;
	ornt = bpornt2ornt[nornt0][nornt1];
	return *this;
}
atom_in_unitcell unitcell[8];
int unitcell_init(){
	unitcell[0].set(0, 0, 0, 3, 7, 1);
	unitcell[1].set(0, 2, 2, 1, 7, 1);
	unitcell[2].set(2, 0, 2, 1, 4, 2);
	unitcell[3].set(2, 2, 0, 1, 4, 2);
	unitcell[4].set(1, 1, 1, 0, 5, 0);
	unitcell[5].set(1, 3, 3, 0, 5, 0);
	unitcell[6].set(3, 1, 3, 0, 6, 3);
	unitcell[7].set(3, 3, 1, 0, 6, 3);
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
DNAmol& DNAmol::put(int nx, int ny, int nz, int nornt){
	px.set(nx, ny, nz);
	px.adjust();
	ornt = nornt;
	return *this;
}
DNAmol& DNAmol::put(const atom_in_unitcell &a, int x_init, int y_init, int z_init){
	px = a.px.plus(x_init, y_init, z_init);
	px.adjust();
	ornt = a.ornt;
	return *this;
}
int DNAmol::displayPatch()const{
	int i, j;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 8; j++){
			std::cout.put(nt[patch[i][j]]);
		}
	}
	std::cout << std::endl;
	return 0;
}
int DNAmol::displayPatch(std::ofstream &out)const{
	int i, j;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 8; j++){
			out.put(nt[patch[i][j]]);
		}
	}
	out << std::endl;
	return 0;
}
int stage[_Nx][_Ny][_Nz]; // To store DNAmol serial. -1 if not occupied. Using periodic boundary conditions.
DNAmol *mol;
int N; // Total molecules.

char nt[4] = { 'A', 'C', 'G', 'T' };
short ntSerial(char which_nt){
	switch (which_nt){
	case'A':return 0;
	case'C':return 1;
	case'G':return 2;
	case'T':return 3;
    default:std::cout << "Wrong nucleotide!" << std::endl; return -1;
	}
}
short ntSerialPair[4] = { 3, 2, 1, 0 };
int couldPatchInteract[4][4] = {
	{ 0, 1, 1, 0 },
	{ 1, 0, 0, 1 },
	{ 1, 0, 0, 1 },
	{ 0, 1, 1, 0 }
};
short couldPatchInteract_paraJudge[4] = { 1, 3, 0, 2 };
bool couldPatchInteract_ornt(int ornt0, int n0ps, int ornt1, int n1ps){
	return(ornt2bpornt[ornt0][couldPatchInteract_paraJudge[n0ps]] + ornt2bpornt[ornt1][couldPatchInteract_paraJudge[n1ps]] == 7);
}
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

	// u_x * u_y * u_z unitcells
	int u_x = 3, u_y = 3, u_z = 3;
	int i, j, k, m;
	_N = u_x * u_y * u_z * 8;
	mol = new DNAmol[_N];
	for (i = 0; i < u_x; i++){
		for (j = 0; j < u_y; j++){
			for (k = 0; k < u_z; k++){
				for (m = 0; m < 8; m++){
					mol[N].put(unitcell[m], i * 4 + 1, j * 4 + 1, k * 4 + 1);
					stage[mol[N].px.x][mol[N].px.y][mol[N].px.z] = N;
					N++;
				}
			}
		}
	}
	if (N != _N){
		std::cout << "The number of molecules is incorrect." << std::endl;
	}

	return 0;
}

int correctbonding(int n0serial, int n01x, int n01y, int n01z){
	ppos npx = mol[n0serial].px.plus(n01x, n01y, n01z);
	npx.adjust();
	int n1serial = stage[npx.x][npx.y][npx.z];
	int ornt0 = 4 * (n01x + 1) / 2 + 2 * (n01y + 1) / 2 + (n01z + 1) / 2;
	int ornt1 = 7 - ornt0;
	int n0ps = -1, n1ps = -1; // patch serial for the right bonding. -1 if not exist.
	int temp_ntSerial;
	n0ps = findPatchSerial[mol[n0serial].ornt][ornt0];
	n1ps = findPatchSerial[mol[n1serial].ornt][ornt1];
	if (n0ps >= 0 && n1ps >= 0 && couldPatchInteract[n0ps][n1ps] && couldPatchInteract_ornt(mol[n0serial].ornt,n0ps,mol[n1serial].ornt,n1ps)){
		mol[n0serial].correctbond[n0ps] = n1serial;
		mol[n1serial].correctbond[n1ps] = n0serial;
		// patch assign
		static std::uniform_int_distribution<> patchDis(0, 3);
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
	static ppos npx;
	for (m = 0; m < N; m++){
		for (i = -1; i <= 1; i+=2){
			for (j = -1; j <= 1; j+=2){
				for (k = -1; k <= 1; k+=2){
					npx = mol[m].px.plus(i, j, k);
					npx.adjust();
					if (stage[npx.x][npx.y][npx.z]>m){
						correctbonding(m, i, j, k);
					}
				}
			}
		}
	}
	std::ofstream out("F:\\DNAPatches.txt");
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