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

	return 0;
}

