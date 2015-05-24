// GraduateProject_DNASelfAssembling.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

mt19937 gen(123456);
clock_t t_start;
clock_t t_end;

int programInitiation(){
	ornt2bpornt_init();
	unitcell_init();
	return 0;
}
int _tmain(int argc, _TCHAR* argv[]){
	t_start = clock();
	programInitiation();
	patchPreset();
	simulationPrepare();
	simulationProcess();

	return 0;
}
