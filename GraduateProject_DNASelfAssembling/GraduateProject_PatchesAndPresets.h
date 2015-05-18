#pragma once

struct atom_in_unitcell{
	ppos px;
	int sd; // Shared Dimension
	int ornt; // 0, 1, ..., 23
	atom_in_unitcell set(int nx, int ny, int nz, int nsd, int nornt0, int nornt1);
};
extern atom_in_unitcell unitcell[];
int unitcell_init();

struct DNAmol{
	ppos px;
	int ornt; // 0, 1, 2, ..., 23
	char patch[4][8]; // 32-nt patch sequences are divided into 4 parts, each ranging from the center to the ends.
	int correctbond[4]; // -1 by default and if not bonded.
	DNAmol();
	DNAmol put(int nx, int ny, int nz, int nornt);
	DNAmol put(const atom_in_unitcell &a, int x_init, int y_init, int z_init);
	int displayPatch()const;
};
extern DNAmol *mol;
extern int stage[_Nx][_Ny][_Nz]; // To store DNAmol serial. -1 if not occupied. Using periodic boundary conditions.
extern int N; // Total molecules.

int patchPreset();