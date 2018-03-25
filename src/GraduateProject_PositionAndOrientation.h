#pragma once

#include <iostream>

// periodic stage size
constexpr int _Nx = 30;
constexpr int _Ny = 30;
constexpr int _Nz = 30;

struct ppos{
	int x, y, z;
	ppos& adjust();
	ppos(){}
	ppos(int nx, int ny, int nz){ x = nx; y = ny; z = nz; }
	ppos& set(int nx, int ny, int nz);
	ppos operator+(const ppos &a)const;
	ppos operator-()const;
	ppos operator-(const ppos &a)const;
	bool operator==(const ppos &a)const;
	ppos plus(int ax, int ay, int az)const;
	ppos crossproduct(int ax, int ay, int az)const;
	friend std::ostream& operator<<(std::ostream &os, const ppos &px);
};
std::ostream& operator<<(std::ostream &os, const ppos &px);
int mod2(const ppos &a);


//ornt2bpornt[whichOrnt][bpSerial] gives an integer between 0 - 7
extern short ornt2bpornt[][4];
int ornt2bpornt_init();
//bpornt2ornt[serial 0 ornt][serial 1 ornt] gives the orientation, or -1 if not valid
extern short bpornt2ornt[8][8];
int bpornt2ornt_init();
//findPatchSerial[whichOrnt][bpOrnt] gives the patch serial 0 - 3, or -1 if not exist
extern short findPatchSerial[][8];
int findPatchSerial_init();
//bporntRot[bpOrnt][axis x, y, z][times of 90-degree rotation] gives the new base-pair orientation
extern short bporntRot[8][3][4];
int bporntRot_init();
//orntRot[whichOrnt][axis x, y, z][times of 90-degree rotation] gives the new molecule orientation
extern short orntRot[24][3][4];
int orntRot_init();