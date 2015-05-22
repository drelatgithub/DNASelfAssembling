// GraduateProject_DNASelfAssembling.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

mt19937 gen(123456);

int programInitiation(){
	ornt2bpornt_init();
	unitcell_init();
	return 0;
}
int _tmain(int argc, _TCHAR* argv[])
{
	programInitiation();
	patchPreset();
	simulationPrepare();
	simulationProcess();

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
	}

	return 0;
}

