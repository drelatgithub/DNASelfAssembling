#pragma once

extern const double NN_dH_kcal[4][4][4];
extern const double NN_dS_cal[4][4][4];
/*
	NN_dH_kcal[a][b][c], NN_dS_cal[a][b][c]

	5'  ---  a  ---  b  ---  3'
	         |       |
	3'  ---  a' ---  c  ---  5'

	where a' is the paired nucleotide of a.
	If c = b' then the index [a][b][b'] and [b'][a'][a] must refer to the same value.
*/
//int NN_dG_check();
extern const double NN_dH_init_kcal;
extern const double NN_dS_init_cal;
extern const double NN_dH_ATPenalty_kcal;
extern const double NN_dS_ATPenalty_cal;

extern const double NN_dH_DE_kcal[2][4][4];
extern const double NN_dS_DE_cal[2][4][4];
/*
	NN_dH_DE_kcal[0,1][a][b], NN_dS_DE_cal[0,1][a][b]

	0 mode: 5' dangling ends

		5'  ---  b  ---  a  ---  3'
		                 |
						 a' ---  5'

	1 mode: 3' dangling ends

		                 a' ---  3'
						 |
		3'  ---  b  ---  a  ---  5'
*/
//int NN_dG_DE_check();
