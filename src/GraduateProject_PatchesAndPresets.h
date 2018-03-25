#pragma once

struct atom_in_unitcell{
	ppos px;
	int sd; // Shared Dimension
	int ornt; // 0, 1, ..., 23
	atom_in_unitcell& set(int nx, int ny, int nz, int nsd, int nornt0, int nornt1);
};
extern atom_in_unitcell unitcell[];
int unitcell_init();

struct DNAmol{
	ppos px;
	int ornt; // 0, 1, 2, ..., 23
	short patch[4][8]; // 32-nt patch sequences are divided into 4 parts, each ranging from the center to the ends.
	short correctbond[4]; // -1 by default and if not bonded.
	short status; // 0 - not mentioned; 1 - used; 2 - in current cluster
	int links;
	int link_with[8];
	double link_with_pRatio[8];
	ppos px_cnt_backup;
	ppos px_fwd_dest;
	int ornt_cnt_backup, ornt_fwd_dest;
	double cluster_eff_rad_diff2;
	int outputTrace; // marked as different values for output use
	DNAmol();
	DNAmol& put(int nx, int ny, int nz, int nornt);
	DNAmol& put(const atom_in_unitcell &a, int x_init, int y_init, int z_init);
	int displayPatch()const;
	int displayPatch(ofstream &out)const;
};
extern DNAmol *mol;
extern int stage[_Nx][_Ny][_Nz]; // To store DNAmol serial. -1 if not occupied. Using periodic boundary conditions.
extern int N; // Total molecules.

extern char nt[4]; // A C G T
short ntSerial(char which_nt);
extern short ntSerialPair[4];
extern int couldPatchInteract[4][4];
extern short couldPatchInteract_paraJudge[4];
bool couldPatchInteract_ornt(int ornt0, int n0ps, int ornt1, int n1ps);
extern int isPatchConsequent[4]; // If the nucleotide sequence in a patch is 5'-3', it is consequent.

int patchPreset();