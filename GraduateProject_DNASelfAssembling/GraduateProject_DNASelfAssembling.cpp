// GraduateProject_DNASelfAssembling.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

random_device rd;
mt19937 gen(rd());
clock_t t_start;
clock_t t_end;

int programInitiation(){
	ornt2bpornt_init();
	bpornt2ornt_init();
	findPatchSerial_init();
	bporntRot_init();
	orntRot_init();

	unitcell_init();
	patchPreset();

	return 0;
}
int _tmain(int argc, _TCHAR* argv[]){
	programInitiation();
	t_start = clock();
	simulationPrepare();
	simulationProcess();


	return 0;
}
