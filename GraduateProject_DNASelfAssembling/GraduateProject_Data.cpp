#include"stdafx.h"

const double NN_dH_kcal[4][4][4] = {
	{ // AX/TY
		{ 1.2, 2.3, -0.6, -7.6 },
		{ 5.3, 0.0, -8.4, 0.7 },
		{ -0.7, -7.8, -3.1, 1.0 },
		{ -7.2, -1.2, -2.5, -2.7 }
	},
	{ // CX/GY
		{ -0.9, 1.9, -0.7, -8.5 },
		{ 0.6, -1.5, -8.0, -0.8 },
		{ -4.0, -10.6, -4.9, -4.1 },
		{ -7.8, -1.5, -2.8, -5.0 }
	},
	{ // GX/CY
		{ -2.9, 5.2, -0.6, -8.2 },
		{ -0.7, 3.6, -9.8, 2.3 },
		{ 0.5, -8.0, -6.0, 3.3 },
		{ -8.4, 5.2, -4.4, -2.2 }
	},
	{ // TX/AY
		{ 4.7, 3.4, 0.7, -7.2 },
		{ 7.6, 6.1, -8.2, 1.2 },
		{ 3.0, -8.5, 1.6, -0.1 },
		{ -7.6, 1.0, -1.3, 0.2 }
	}
};
const double NN_dS_cal[4][4][4] = {
	{ // AX/TY
		{ 1.7, 4.6, -2.3, -21.3 },
		{ 14.6, -4.4, -22.4, 0.2 },
		{ -2.3, -21.0, -9.5, 0.9 },
		{ -20.4, -6.2, -8.3, -10.8 }
	},
	{ // CX/GY
		{ -4.2, 3.7, -2.3, -22.7 },
		{ -0.6, -7.2, -19.9, -4.5 },
		{ -13.2, -27.2, -15.3, -11.7 },
		{ -21.0, -6.1, -8.0, -15.8 }
	},
	{ // GX/CY
		{ -9.8, 14.2, -1.0, -22.2 },
		{ -3.8, 8.9, -24.4, 5.4 },
		{ 3.2, -19.9, -15.8, 10.4 },
		{ -22.4, 13.5, -12.3, -8.4 }
	},
	{ // TX/AY
		{ 12.9, 8.0, 0.7, -21.3 },
		{ 20.2, 16.4, -22.2, 0.7 },
		{ 7.4, -22.7, 3.6, -1.7 },
		{ -21.3, 0.7, -5.3, -1.5 }
	}
};
int NN_dG_check(){
	double check[4][4][4];
	int i, j, k;
	double T37 = 273.15 + 37;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){
			for (k = 0; k < 4; k++){
				check[i][j][k] = NN_dH_kcal[i][j][k] - T37*NN_dS_cal[i][j][k] / 1000;
				cout << "  \t" << check[i][j][k];
			}
			cout << endl;
		}
		cout << endl;
	}
	return 0;
}
const double NN_dH_init_kcal = 0.2;
const double NN_dS_init_cal = -5.7;
const double NN_dH_ATPenalty_kcal = 2.2;
const double NN_dS_ATPenalty_cal = 6.9;
/*
	References:

	1. J. SantaLucia, Jr. and D. Hicks, The Thermodynamics of DNA Structural Motifs, Annu. Rev. Biophys. Biomol. Struct. 33, 415 (2004)
	2. Allawi HT, SantaLucia J Jr. Thermodynamics and NMR of internal G-T mismatches in DNA. Biochemistry 36:10581-94 (1997)
	3. Allawi HT, SantaLucia J Jr. Nearestneighbor thermodynamics of internal A-C mismatches in DNA: sequence dependence and pH effects. Biochemistry 37:9435-44 (1998)
	4. Allawi HT, SantaLucia J Jr. Nearestneighbor thermodynamics parameters for internal G-A mismatches in DNA. Biochemistry 37:2170-79 (1998)
	5. Allawi HT, SantaLucia J Jr. Thermodynamics of internal C-T mismatches in DNA. Nucleic Acids Res. 26:2694-701 (1998)
	6. Peyret N, Seneviratne PA, Allawi HT, SantaLucia J Jr. Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A-A, C-C, G-G, and T-T mismatches. Biochemistry 38:3468-77 (1999)
*/

extern const double NN_dH_DE_kcal[2][4][4] = {
	{
		{ 0.2, 0.6, -1.1, -6.9 },
		{ -6.3, -4.4, -5.1, -4.0 },
		{ -3.7, -4.0, -3.9, -4.9 },
		{ -2.9, -4.1, -4.2, -0.2 }
	},
	{
		{ -0.5, 4.7, -4.1, -3.8 },
		{ -5.9, -2.6, -3.2, -5.2 },
		{ -2.1, -0.2, -3.9, -4.4 },
		{ -0.7, 4.4, -1.6, 2.9 }
	}
};
extern const double NN_dS_DE_cal[2][4][4] = {
	{
		{ 2.3, 3.3, -1.6, -20.0 },
		{ -17.1, -12.6, -14.0, -10.9 },
		{ -10.0, -11.9, -10.9, -13.8 },
		{ -7.6, -13.0, -15.0, -0.5 }
	},
	{
		{ -1.1, 14.2, -13.1, -12.6 },
		{ -16.5, -7.4, -10.4, -15.0 },
		{ -3.9, -0.1, -11.2, -13.1 },
		{ -0.8, 14.9, -3.6, 10.4 }
	}
};
int NN_dG_DE_check(){
	double check[2][4][4];
	int i, j, k;
	double T37 = 273.15 + 37;
	for (i = 0; i < 2; i++){
		for (j = 0; j < 4; j++){
			for (k = 0; k < 4; k++){
				check[i][j][k] = NN_dH_DE_kcal[i][j][k] - T37*NN_dS_DE_cal[i][j][k] / 1000;
				cout << "  \t" << check[i][j][k];
			}
			cout << endl;
		}
		cout << endl;
	}
	return 0;
}
/*
	References:

	1. J. SantaLucia, Jr. and D. Hicks, The Thermodynamics of DNA Structural Motifs, Annu. Rev. Biophys. Biomol. Struct. 33, 415 (2004)
	2. Bommarito S, Peyret N, SantaLucia J Jr. Thermodynamic parameters for DNA sequences with dangling ends. Nucleic Acids Res. 28:1929-34 (2000)
*/
