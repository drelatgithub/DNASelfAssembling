#pragma once

// periodic stage size
#define _Nx 12
#define _Ny 12
#define _Nz 12

struct ppos{
	int x, y, z;
	ppos adjust();
	ppos(){adjust(); }
	ppos(int nx, int ny, int nz){ x = nx; y = ny; z = nz; adjust(); }
	ppos operator+(const ppos &a)const;
	ppos operator-()const;
	ppos operator-(const ppos &a)const;
	friend ostream& operator<<(ostream &os, const ppos &px);
};
ostream& operator<<(ostream &os, const ppos &px);


//ornt2bpornt[whichOrnt][bpSerial] gives an integer between 0 - 7
extern int ornt2bpornt[][4];
int ornt2bpornt_init();
