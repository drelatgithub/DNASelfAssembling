#include"stdafx.h"

ppos ppos::operator+(const ppos &a)const{
	return ppos(x + a.x, y + a.y, z + a.z);
}
ppos ppos::operator-()const{
	return ppos(-x, -y, -z);
}
ppos ppos::operator-(const ppos &a)const{
	return ppos(x - a.x, y - a.y, z - a.z);
}
ppos ppos::adjust(){
	x = (x < 0) ? (x%_Nx + _Nx) % _Nx : x%_Nx;
	y = (y < 0) ? (y%_Ny + _Ny) % _Ny : y%_Ny;
	z = (z < 0) ? (z%_Nz + _Nz) % _Nz : z%_Nz;
	return *this;
}
ostream& operator<<(ostream &os, const ppos &px){
	os << '(' << px.x << ',' << px.y << ',' << px.z << ')';
	return os;
}


int ornt2bpornt[24][4];
int ornt2bpornt_init(){
	int mainOrnt, subOrnt;
	int i, j, a[3], b[3];
	for (i = 0; i < 24; i++){
		mainOrnt = i % 8; // 0, 1, ..., 7
		subOrnt = i / 8; // 0, 1, 2
		ornt2bpornt[i][0] = mainOrnt;
		ornt2bpornt[i][1] = mainOrnt ^ (7 - (1 << subOrnt));
		/*
		The 8-nt domains are defined in the following manner:

		                   Domain_2(1)   Domain_1(0)
		combined_head(5')                             tail(3')
		                   Domain_3(2)   Domain_4(3)

		Base-pair orientations 2 and 3 are determined so that rotation along 0->1->2->0 would give the angular speed towards 3.
		Reference: Y. Ke, L. L. Ong, W. M. Shih, and P. Yin, Three-Dimensional Structures Self-Assembled from DNA Bricks, Science 338, 1177 (2012)
		*/
		ornt2bpornt[i][2] = ornt2bpornt[i][3] = 0;
		for (j = 0; j < 3; j++){
			a[j] = (ornt2bpornt[i][0] >> (2 - j)) % 2;
			b[j] = (ornt2bpornt[i][1] >> (2 - j)) % 2;
		}
		for (j = 0; j < 3; j++){
			ornt2bpornt[i][2] <<= 1;
			ornt2bpornt[i][3] <<= 1;
			ornt2bpornt[i][2] += b[(j + 1) % 3] * a[(j + 2) % 3] - a[(j + 1) % 3] * b[(j + 2) % 3] + (-a[j % 3] - b[j % 3] + a[(j + 1) % 3] - b[(j + 1) % 3] - a[(j + 2) % 3] + b[(j + 2) % 3]) / 2 + 1;
			ornt2bpornt[i][3] += a[(j + 1) % 3] * b[(j + 2) % 3] - b[(j + 1) % 3] * a[(j + 2) % 3] + (-a[j % 3] - b[j % 3] - a[(j + 1) % 3] + b[(j + 1) % 3] + a[(j + 2) % 3] - b[(j + 2) % 3]) / 2 + 1;
		}
	}
	return 0;
}
int findPatchSerial[24][8];
int findPatchSerial_init(){
	int i, j;
	for (i = 0; i < 24; i++){
		for (j = 0; j < 8; j++)findPatchSerial[i][j] = -1;
	}
	for (i = 0; i < 24; i++){
		for (j = 0; j < 4; j++)findPatchSerial[i][ornt2bpornt[i][j]] = j;
	}
	return 0;
}