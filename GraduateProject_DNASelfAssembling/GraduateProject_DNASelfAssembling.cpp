// GraduateProject_DNASelfAssembling.cpp : �������̨Ӧ�ó������ڵ㡣
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

	for (int i = 0; i < N; i++){
		cout << i << '\t';
		mol[i].displayPatch();
	}

	return 0;
}

