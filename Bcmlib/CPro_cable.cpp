#include "stdafx.h"

#include "cpro_lpr.h"
#include "cpro_cable.h"

#ifndef ___ABRIDGE_PROFILE_MODE___
#include "ccable3d.h"
#include "CCells.h"
#endif
#include "unit_mes.h"

#ifdef ___WINDOWS_LOG_MESSAGE___
#define  Message(Msg)    theMainFrame->Message(Msg)
#else
#define  Message(Msg)    printf(Msg);  printf("\n")
#endif

#ifdef ___PRO_CABLE2D_cpp___
/////////////////////////////////////
//...число типов проводов и их имена;
int    GetCABLESampleCount(void) { return(NUM_CABLE_SAMPLES); }
char * GetCABLESampleName (int N_sm)
{
  static char * S[] = { "ќптический кабель",			//_CABLE1
								"ќптический грозотрос",		//_CABLE2
								"ќптический провод",			//_CABLE3
  };
  return S[N_sm];
}

/*================================================*/
/*        “јЅЋ»÷џ ќѕ“ќ-¬ќЋќ ќЌЌџ’  јЅ≈Ћ≈…         */
/*================================================*/
//////////////////////////////////////////
//...таблицы параметров опто-волоконных кабелей;
Table * cable_OKLG(int id_static_char)
{
	static double Diam [] = { 0., 13.8, 13.8, 17.0, 18.6, 13.8, 15.7, 14.6, 15.4, 15.0,	15.4,	
											15.6,	16.1, 14.0, 13.4, 12.0, 
											11.7, 13.0,  15.9, 15.4, 13.25, 9.1, 11.1,  
											 9.2, 11.45, 12.1, 13.3, 13.9, 13.9, 14.6, 15.0, 15.0, 15.95};
	static double F_eff[] = { 0., ((int)(M_PI*Diam[1]*Diam[1]*2.5))*.1,  ((int)(M_PI*Diam[2]*Diam[2]*2.5))*.1,  ((int)(M_PI*Diam[3]*Diam[3]*2.5))*.1, ((int)(M_PI*Diam[4]*Diam[4]*2.5))*.1, ((int)(M_PI*Diam[5]*Diam[5]*2.5))*.1, 
											((int)(M_PI*Diam[6]*Diam[6]*2.5))*.1,  ((int)(M_PI*Diam[7]*Diam[7]*2.5))*.1,  ((int)(M_PI*Diam[8]*Diam[8]*2.5))*.1, ((int)(M_PI*Diam[9]*Diam[9]*2.5))*.1, ((int)(M_PI*Diam[10]*Diam[10]*2.5))*.1,
											((int)(M_PI*Diam[11]*Diam[11]*2.5))*.1,((int)(M_PI*Diam[12]*Diam[12]*2.5))*.1, 110.8,                               100.0,                                
											75.4 , 52.9,  58.0,  93.1, 142.4,  107.2,  50.,    74.5, 
											44.43, 76.02, 80.5, 100.3, 109.35, 109.35, 121.07, 126.47, 126.47, 145.18};
	static double E_crp[] = { 0., 17100., 16300., 14130., 18920., 16300.,  16220., 19410., 17450., 22050., 20170., 
											19650., 18110., 75200., 108000., 148000., 
											145400., 177700., 141000., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double E_ini[] = { 0., 17100., 16300., 14130., 18920., 16300.,  16220., 19410., 17450., 22050., 20170., 
											19650., 18110., 85600., 108000., 151200., 
											145400., 177700., 141000., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double E_fin[] = { 0., 17100., 17620., 14130., 18920., 17620.,  16220., 20990., 18860., 22050., 20170., 
											19650., 18110., 92000, 108000., 162000., 
											145400., 177700., 141000., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double M_cable[] = { 0., 153., 155., 231., 265., 155., 206., 180., 200., 184., 204., 
											  204., 216., 473., 496., 528.6, 
											  373., 495., 628., 755., 484., 418., 623., 
											  308., 351., 376., 394., 497., 549., 534., 537., 480., 619.};
	static double F_break[] = { 0., 46600., 43600., 56300., 89500., 43600., 55000., 56600., 58000., 65900.,	65500., 
											  65400., 64400., 67300., 75400., 91200., 
											  54230., 72500., 77500., 96826., 65000., 61200., 78300., 
											  (int)(5272.*M_G+.5), (int)(5370.*M_G+.5), (int)(5705.*M_G+.5), (int)(5489.*M_G+.5), (int)(7530.*M_G+.5), 
											  (int)(8730.*M_G+.5), (int)(7954.*M_G+.5), (int)(7716.*M_G+.5), (int)(6433.*M_G+.5), (int)(8916.*M_G+.5)};
	static double A_temp [] = { 0., 1.62e-6, 1.44e-6, 1.44e-6, 0.95e-6, 1.44e-6, 1.57e-6, 0.92e-6, 1.33e-6, 0.39e-6, 0.84e-6, 
											  0.86e-6, 1.10e-6, 1.72e-5, 1.57e-5, 1.3e-5, 
											  1.50e-5,  1.43e-5,  1.71e-5,  1.45e-5,  1.6e-5,   1.2e-5,   1.2e-5,
											  12.95e-6, 15.97e-6, 15.91e-6, 17.79e-6, 16.17e-6, 15.22e-6, 16.50e-6, 16.98e-6, 18.31e-6, 16.95e-6};
	static char * S0[] = {"ќ Ћ∆-01-6-16-10/125-0.36/0.22-3.5/18-19.5",  "ќ Ћ∆-01-6-20-10/125-0,36/0,22-3,5/18-20",		"ќ Ћ∆-01-6-38-10/125-0,36/0,22-3,5/18-25",    "ќ Ћ∆-01-6-16-10/125-0,36/0,22-3,5/18-40,0",  "ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-20",  
								 "ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-25,0","ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-26,0","ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-26,0","ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-30,0","ќ Ћ∆-“-01-6-36-10/125-0,36/0,22-3,5/18-30,0",
								 "ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-30,0","ќ Ћ∆-“-01-6-70-10/125-0,36/0,22-3,5/18-30,0","OPGW-L1-24K(16G652B+8G655)-74AY/37ACS", 	    "OPGW-C1-24K(16G652B+8G655)-51AL3/49A20SA",   "OPGW-C1-24K(16G652B+8G655)-75A20SA", 
								 "OPGW 17A37z", "OPGW 30C50z", "OPGW 59L63z", "ј—-70/72", "ј∆—-70/39", "—-50", "—-70", 
								 "ќ √“ц-1-24-(G652)-9,2/52",  "ќ √“ц-1-24-(G652)-11,5/53", "ќ √“ц-1-24-(G652)-12,1/56", "ќ √“ц-1-24-(G652)-13,3/54", "ќ √“ц-1-24-(G652)-13,9/74", 
								 "ќ √“ц-1-24-(G652)-13,9/85", "ќ √“ц-1-24-(G652)-14,6/78", "ќ √“ц-1-24-(G652)-15/75",   "ќ √“ц-1-24-(G652)-15/63",   "ќ √“ц-1-24-(G652)-16/87",
	};
	Table * table = GetCableDataBaseTable(32, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = Diam;
		table->table[1].parm_values = F_eff;
		table->table[2].parm_values = E_crp;
		table->table[3].parm_values = E_ini;
		table->table[4].parm_values = E_fin;
		table->table[5].parm_values = M_cable;
		table->table[6].parm_values = F_break;
		table->table[7].parm_values = A_temp;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = Diam[k];
		((double *)(table->table[1].parm_values))[k] = F_eff[k];
		((double *)(table->table[2].parm_values))[k] = E_crp[k];
		((double *)(table->table[3].parm_values))[k] = E_ini[k];
		((double *)(table->table[4].parm_values))[k] = E_fin[k];
		((double *)(table->table[5].parm_values))[k] = M_cable[k];
		((double *)(table->table[6].parm_values))[k] = F_break[k];
		((double *)(table->table[7].parm_values))[k] = A_temp[k];
	}
	return table;
}

Table * cable_OPGW(int id_static_char) //...пришлось добавить в список cable_OKLG;
{
	static double Diam [] = { 0., 11.7, 13.0,  15.9, 15.4, 13.25, 9.1, 11.1,  
											 9.2, 11.45, 12.1, 13.3, 13.9, 13.9, 14.6,   15.0,   15.0,   15.95};
	static double F_eff[] = { 0., 52.9,  58.0,  93.1, 142.4, 107.2,  50.,    74.5, 
											44.43, 76.02, 80.5, 100.3, 109.35, 109.35, 121.07, 126.47, 126.47, 145.18};
	static double E_crp[] = { 0., 145400., 177700., 141000., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double E_ini[] = { 0., 145400., 177700., 141000., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double E_fin[] = { 0., 148220., 177700., 143730., 133273., 114544., 200000., 200000., 
											162590., 106360., 107130.,  88050., 104030., 116470., 100370., 95390., 83880., 95710.};
	static double M_cable[] = { 0., 373., 495., 628., 755., 484., 418., 623., 
											  308., 351., 376., 394., 497., 549., 534., 537., 480., 619.};
	static double F_break[] = { 0., 54230., 72500., 77500., 96826., 65000., 61200., 78300., 
											  (int)(5272.*M_G+.5), (int)(5370.*M_G+.5), (int)(5705.*M_G+.5), (int)(5489.*M_G+.5), (int)(7530.*M_G+.5), 
											  (int)(8730.*M_G+.5), (int)(7954.*M_G+.5), (int)(7716.*M_G+.5), (int)(6433.*M_G+.5), (int)(8916.*M_G+.5)};
	static double A_temp [] = { 0., 1.50e-5,  1.43e-5,  1.71e-5,  1.45e-5,  1.6e-5,   1.2e-5,   1.2e-5,
											  12.95e-6, 15.97e-6, 15.91e-6, 17.79e-6, 16.17e-6, 15.22e-6, 16.50e-6, 16.98e-6, 18.31e-6, 16.95e-6};
	static char * S0[] = { "OPGW 17A37z", "OPGW 30C50z", "OPGW 59L63z", "ј—-70/72", "ј∆—-70/39", "—-50", "—-70", 
								  "ќ √“ц-1-24-(G652)-9,2/52",  "ќ √“ц-1-24-(G652)-11,5/53", "ќ √“ц-1-24-(G652)-12,1/56", "ќ √“ц-1-24-(G652)-13,3/54", "ќ √“ц-1-24-(G652)-13,9/74", 
								  "ќ √“ц-1-24-(G652)-13,9/85", "ќ √“ц-1-24-(G652)-14,6/78", "ќ √“ц-1-24-(G652)-15/75",   "ќ √“ц-1-24-(G652)-15/63",   "ќ √“ц-1-24-(G652)-16/87",
	};
	Table * table = GetCableDataBaseTable(17, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = Diam;
		table->table[1].parm_values = F_eff;
		table->table[2].parm_values = E_crp;
		table->table[3].parm_values = E_ini;
		table->table[4].parm_values = E_fin;
		table->table[5].parm_values = M_cable;
		table->table[6].parm_values = F_break;
		table->table[7].parm_values = A_temp;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = Diam[k];
		((double *)(table->table[1].parm_values))[k] = F_eff[k];
		((double *)(table->table[2].parm_values))[k] = E_crp[k];
		((double *)(table->table[3].parm_values))[k] = E_ini[k];
		((double *)(table->table[4].parm_values))[k] = E_fin[k];
		((double *)(table->table[5].parm_values))[k] = M_cable[k];
		((double *)(table->table[6].parm_values))[k] = F_break[k];
		((double *)(table->table[7].parm_values))[k] = A_temp[k];
	}
	return table;
}

Table * cable_OKGT(int id_static_char) //...пришлось добавить в список cable_OPGW и cable_OKLG;
{
	static double Diam [] = { 0.,  9.2,   11.45,  12.1,    13.3,  13.9,   13.9,   14.6,   15.0,   15.0,   15.95};
	static double F_eff[] = { 0., 44.43,  76.02,  80.5,   100.3, 109.35, 109.35, 121.07, 126.47, 126.47, 145.18};
	static double E_crp[] = { 0., 162590., 106360., 107130., 88050., 104030., 116470., 100370.,  95390.,  83880.,  95710.};
	static double E_ini[] = { 0., 162590., 106360., 107130., 88050., 104030., 116470., 100370.,  95390.,  83880.,  95710.};
	static double E_fin[] = { 0., 162590., 106360., 107130., 88050., 104030., 116470., 100370.,  95390.,  83880.,  95710.};
	static double M_cable[] = { 0., 308., 351., 376., 394., 497., 549., 534., 537., 480., 619.};
	static double F_break[] = { 0., (int)(5272.*M_G+.5), (int)(5370.*M_G+.5), (int)(5705.*M_G+.5), (int)(5489.*M_G+.5), (int)(7530.*M_G+.5), 
											  (int)(8730.*M_G+.5), (int)(7954.*M_G+.5), (int)(7716.*M_G+.5), (int)(6433.*M_G+.5), (int)(8916.*M_G+.5)};
	static double A_temp [] = { 0., 12.95e-6, 15.97e-6, 15.91e-6, 17.79e-6, 16.17e-6, 15.22e-6, 16.50e-6, 16.98e-6, 18.31e-6, 16.95e-6};
	static char * S0[] = { "ќ √“ц-1-24-(G652)-9,2/52",  "ќ √“ц-1-24-(G652)-11,5/53", "ќ √“ц-1-24-(G652)-12,1/56", "ќ √“ц-1-24-(G652)-13,3/54", "ќ √“ц-1-24-(G652)-13,9/74", 
								  "ќ √“ц-1-24-(G652)-13,9/85", "ќ √“ц-1-24-(G652)-14,6/78", "ќ √“ц-1-24-(G652)-15/75",   "ќ √“ц-1-24-(G652)-15/63",   "ќ √“ц-1-24-(G652)-16/87",
	};
	Table * table = GetCableDataBaseTable(10, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = Diam;
		table->table[1].parm_values = F_eff;
		table->table[2].parm_values = E_crp;
		table->table[3].parm_values = E_ini;
		table->table[4].parm_values = E_fin;
		table->table[5].parm_values = M_cable;
		table->table[6].parm_values = F_break;
		table->table[7].parm_values = A_temp;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = Diam[k];
		((double *)(table->table[1].parm_values))[k] = F_eff[k];
		((double *)(table->table[2].parm_values))[k] = E_crp[k];
		((double *)(table->table[3].parm_values))[k] = E_ini[k];
		((double *)(table->table[4].parm_values))[k] = E_fin[k];
		((double *)(table->table[5].parm_values))[k] = M_cable[k];
		((double *)(table->table[6].parm_values))[k] = F_break[k];
		((double *)(table->table[7].parm_values))[k] = A_temp[k];
	}
	return table;
}

///////////////////////////////////////////////////////////////
// ...база данных стандартных кабелей (только основной список);
int cable_DB(void ** cable_list, char * cable, double * par)
{
//////////////////////////////////////////////////////////////////////
//...ищем марку кабел€ в списке базы данных и заполн€ем его параметры;
	char * sss, * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	memset(par,  0, NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
	setlocale(LC_COLLATE, "rus"); 
	for (int k = 0; cable_list && cable_list[k]; k++) { 
		Table * table = (Table *)GetTable(cable_list[k]);
		for (int i = 1; table && i <= table->N_group; i++) 
			if ((sss = GetProfileName(cable_list[k], i)) != NULL && ! _strnicmp(cable, sss, strlen(sss))) { //...сравниваем без учета регистра;
				for (int j = 0; j < NUM_OF_CABLE_BASE_ELEMENTS; j++) 
				par[j] = GetTableParam(cable_list[k], i, j);
				setlocale(LC_COLLATE, sloc); free(sloc);
				return(1);
			}
	}
	setlocale(LC_COLLATE, sloc); free(sloc);
	return(0);
}

////////////////////////////////////////////////////////////////////////
// ...база данных стандартных кабелей (включа€ дополнительные проверки);
void CableBase(void ** cable_list, char * cable, double * par)
{
	static double cable_par[][NUM_OF_CABLE_BASE_ELEMENTS] = {
		{13.8, ((int)(M_PI*13.8*13.8*2.5))*.1, 17100., 17100., 17100., 153., 46600., 1.62e-6}, //[0] //...ќ Ћ∆-01-6-16-10/125-0.36/0.22-3.5/18-19.5 (ќ Ћ∆-19, в —аратове, T_max = 0.42*T_break)
		{13.8, ((int)(M_PI*13.8*13.8*2.5))*.1, 16300., 16300., 17620., 155., 43600., 1.44e-6}, //[1] //...ќ Ћ∆-01-6-20-10/125-0,36/0,22-3,5/18-20 (ќ Ћ∆-20, как ќ Ћ∆“-20)
		{17.0, ((int)(M_PI*17.0*17.0*2.5))*.1, 14130., 14130., 14130., 231., 56300., 1.44e-6}, //[2] //...ќ Ћ∆-01-6-38-10/125-0,36/0,22-3,5/18-25 (ќ Ћ∆-25)
		{13.8, ((int)(M_PI*13.8*13.8*2.5))*.1, 16300., 16300., 17620., 155., 43600., 1.44e-6}, //[3] //...ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-20 (ќ Ћ∆“-20)
		{15.7, ((int)(M_PI*15.7*15.7*2.5))*.1, 16220., 16220., 16220., 206., 55000., 1.57e-6}, //[4] //...ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-25,0 (ќ Ћ∆“-25)
		{14.6, ((int)(M_PI*14.6*14.6*2.5))*.1, 19410., 19410., 20990., 180., 56600., 0.92e-6}, //[5] //...ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-26,0 (ќ Ћ∆“-26)
		{15.4, ((int)(M_PI*15.4*15.4*2.5))*.1, 17450., 17450., 18860., 200., 58000., 1.33e-6}, //[6] //...ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-26,0 (ќ Ћ∆“-26)
		{11.7, 52.9,   145400., 145400., 148220., 373., 54230., 1.50e-5}, //[7] //...OPGW 17A37z
		{13.0, 58.0,   177700., 177700., 177700., 495., 72500., 1.43e-5}, //[8] //...OPGW 30C50z (OPGW A17/S47 F8 ???)
		{15.9, 93.1,   141000., 141000., 143730., 628., 77500., 1.71e-5}, //[9] //...OPGW 59L63z
		{15.4, 142.4,  133273., 133273., 133273., 755., 96826., 1.45e-5}, //[10] //...ј—-70/72
		{13.25,107.2,  114544., 114544., 114544., 484., 65000., 1.6e-5}, //[11] //...ј∆—-70/39
		{9.1,  50.,    200000., 200000., 200000., 418., 61200., 1.2e-5}, //[12] //...—-50
		{11.1, 74.5,   200000., 200000., 200000., 623., 78300., 1.2e-5}, //[13] //...—-70
		{9.2,  44.43,  162590., 162590., 162590., 308., (int)(5272.*M_G+.5),12.95e-6}, //[14] //...ќ √“ц-1-24-(G652)-9,2/52 
		{11.45,76.02,  106360., 106360., 106360., 351., (int)(5370.*M_G+.5),15.97e-6}, //[15] //...ќ √“ц-1-24-(G652)-11,5/53
		{12.1, 80.5,   107130., 107130., 107130., 376., (int)(5705.*M_G+.5),15.91e-6}, //[16] //...ќ √“ц-1-24-(G652)-12,1/56
		{13.3, 100.3,   88050.,  88050.,  88050., 394., (int)(5489.*M_G+.5),17.79e-6}, //[17] //...ќ √“ц-1-24-(G652)-13,3/54
		{13.9, 109.35, 104030., 104030., 104030., 497., (int)(7530.*M_G+.5),16.17e-6}, //[18] //...ќ √“ц-1-24-(G652)-13,9/74
		{13.9, 109.35, 116470., 116470., 116470., 549., (int)(8730.*M_G+.5),15.22e-6}, //[19] //...ќ √“ц-1-24-(G652)-13,9/85
		{14.6, 121.07, 100370., 100370., 100370., 534., (int)(7954.*M_G+.5),16.50e-6}, //[20] //...ќ √“ц-1-24-(G652)-14,6/78
		{15.0, 126.47,  95390.,  95390.,  95390., 537., (int)(7716.*M_G+.5),16.98e-6}, //[21] //...ќ √“ц-1-24-(G652)-15/75  
		{15.0, 126.47,  83880.,  83880.,  83880., 480., (int)(6433.*M_G+.5),18.31e-6}, //[22] //...ќ √“ц-1-24-(G652)-15/63  
		{15.95,145.18,  95710.,  95710.,  95710., 619., (int)(8916.*M_G+.5),16.95e-6}, //[23] //...ќ √“ц-1-24-(G652)-16/87
	};

/////////////////////////////////////////////////////
//...логика выбора соответствующих параметров кабел€;
	char * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	setlocale(LC_COLLATE, "rus"); 
	if (cable/* && ! cable_DB(cable_list, cable, par)*/) { //...сравниваем без учета регистра;
		if (! _stricmp (cable, "ќ Ћ∆-01-6-16-10/125-0.36/0.22-3.5/18-19.5"))
			memcpy(par, cable_par[0], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-01-6-20-10/125-0,36/0,22-3,5/18-20"))
			memcpy(par, cable_par[1], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-01-6-38-10/125-0,36/0,22-3,5/18-25"))
			memcpy(par, cable_par[2], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-20"))
			memcpy(par, cable_par[3], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-25,0"))
			memcpy(par, cable_par[4], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-“-01-6-20-10/125-0,36/0,22-3,5/18-26,0"))
			memcpy(par, cable_par[5], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ Ћ∆-“-01-6-40-10/125-0,36/0,22-3,5/18-26,0"))
			memcpy(par, cable_par[6], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "OPGW 17A37z"))
			memcpy(par, cable_par[7], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "OPGW 30C50z"))
			memcpy(par, cable_par[8], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "OPGW 59L63z"))
			memcpy(par, cable_par[9], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ј—-70/72"))
			memcpy(par, cable_par[10], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ј∆—-70/39"))
			memcpy(par, cable_par[11], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "—-50"))
			memcpy(par, cable_par[12], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "—-70"))
			memcpy(par, cable_par[13], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-9,2/52"))
			memcpy(par, cable_par[14], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-11,5/53"))
			memcpy(par, cable_par[15], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-12,1/56"))
			memcpy(par, cable_par[16], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-13,3/54"))
			memcpy(par, cable_par[17], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-13,9/74"))
			memcpy(par, cable_par[18], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-13,9/85"))
			memcpy(par, cable_par[19], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-14,6/78"))
			memcpy(par, cable_par[20], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-15/75"))
			memcpy(par, cable_par[21], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-15/63"))
			memcpy(par, cable_par[22], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! _stricmp (cable, "ќ √“ц-1-24-(G652)-16/87"))
			memcpy(par, cable_par[23], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
		else
		if (! cable_DB(cable_list, cable, par))
			memcpy(par, cable_par[0], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));
	}
	else
	if (! cable)
		memcpy(par, cable_par[0], NUM_OF_CABLE_BASE_ELEMENTS*sizeof(double));

	setlocale(LC_COLLATE, sloc); free(sloc);
}

////////////////////////////////////////////////////////////////////
//...функци€ образовани€ контекста дл€ базы опто-волоконных кабелей;
void * CreateCABLEContext(int N_sm)
{
	if (N_sm < 0 || N_sm >= NUM_CABLE_SAMPLES) N_sm = 0;

	Context * cont = (Context *)new_struct(sizeof(Context));
	if (! cont) return(NULL);

	cont->N           = N_sm+SHIFT_CABLE_SAMPLES;
	cont->static_char = ADDITIONAL_STATE;
	cont->sample_name = GetCABLESampleName(N_sm);
	cont->units       = UNIT_WIRE;

	switch (cont->N-SHIFT_CABLE_SAMPLES) {
      case _CABLE1: cont->table = cable_OKLG(cont->static_char); break;
      case _CABLE2: cont->table = cable_OPGW(cont->static_char); break;
      case _CABLE3: {
			static char * S0[] = {""};
			cont->table = GetCableDataBaseTable(1, S0, cont->static_char); 
		}	break;
  }
  set_default    (cont);
  set_table_units(cont->table, UNIT_WIRE);
  
  SetTableIndex(cont, 1);
  SetUnits     (cont, UNIT_WIRE); return (void *)cont;
}

/*
/*===============================================================================================*/
/*                                 »Ќ»÷»јЋ»«ј÷»я ќЅ–ј«÷ј                                         */
/*===============================================================================================*/
///////////////////////////////////////////////////////////////////
//...пpедваpительна€ инициализаци€ сечений дл€ базы данных кабелей;
#ifndef ___ABRIDGE_PROFILE_MODE___
int cable2D_init(void * context, CGrid * block_nd)
{
	Context * cont = (Context *)context;

	if (! is_Cable2D(cont)  || ! cont->table) return(0);

	cont->left = cont->right = cont->bottom = cont->top = cont->back = cont->front = 0.;

	DeleteSample(cont->sm);

	switch (cont->N-SHIFT_CABLE_SAMPLES) {
		case _CABLE1: 
		case _CABLE2: { 
		} break;
	}
	return(1);
// return(cont->sm != NULL);
}
#endif

#endif
#undef  Message
