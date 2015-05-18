// GraduateProject_DNASelfAssembling.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


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
		cout << i;
		for (int j = 0; j < 4; j++){
			cout << '\t' << mol[i].correctbond[j];
		}
		cout << endl;
	}
	return 0;
}

