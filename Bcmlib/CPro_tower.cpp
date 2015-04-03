#include "stdafx.h"

#include "cpro_lpr.h"
#include "cpro_tower.h"

#ifndef ___ABRIDGE_PROFILE_MODE___
#include "ctower3d.h"
#include "ccebasic.h"
#endif
#include "unit_mes.h"

#ifdef ___WINDOWS_LOG_MESSAGE___
#define  Message(Msg)    theMainFrame->Message(Msg)
#else
#define  Message(Msg)    printf(Msg);  printf("\n")
#endif

#ifdef ___PRO_TOWER3D_cpp___
/////////////////////////////////////
//...число типов проводов и их имена;
int    GetTOWERSampleCount(void) { return(NUM_TOWER_SAMPLES); }
char * GetTOWERSampleName (int N_sm)
{
  static char * S[] = { "Дополнительный список опор",     //_TOWER1
	                     "Железобетонные анкерно-угловые", //_TOWER2
								"Железобетонные промежуточные",   //_TOWER3
                        "Металлические анкерно-угловые",  //_TOWER4
                        "Металлические промежуточные",    //_TOWER5
                        "Деревянные",                     //_TOWER6
  };
  return S[N_sm];
}

/*===============================================*/
/*                 ТАБЛИЦЫ ОПОР                  */
/*===============================================*/
///////////////////////////////
//...доплнительный список опор;
Table * tower_DOP(int id_static_char)
{
	static double h1  [] = { 0., 10.70, 10.50, 10.50, 17.50, 17.50, 17.50, 19.025, 14.10, 14.10, 17.30, 
										  10.50, 10.50, 13.70, 13.70, 13.70, 17.70, 20.20,  14.50, 10.50, 10.50, 
										  25.70, 25.70, 25.70, 15.00, 15.00, 25.50, 25.50,  25.50, 22.50, 23.61,
										  23.61, 23.61, 25.50, 25.50, 25.50, 17.30, 17.10,  17.10, 23.00, 23.50,
										  23.50, 20.20, 20.20, 20.20, 20.20, 23.20, 23.20,  23.20, 23.20, 23.20,
										  13.80, 25.70, 23.00, 19.30, 21.00, 17.00, 16.00,  19.00, 10.50, 15.50,
										  15.50,  6.00, 20.00, 22.00, 13.50, 13.50, 22.50,  22.50, 25.50, 13.70,
										  23.20, 13.50, 14.10, 13.70, 13.70, 14.10, 13.70,  17.30, 17.70, 13.20,
										  7.00};
	static double h2  [] = { 0., 6.50, 8.10, 4.00, 3.50, 3.50, 4.00, 3.50, 3.60, 3.60, 3.60, 
										  3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 6.00, 4.00, 4.00, 4.00, 
										  MAX_HIT, MAX_HIT, MAX_HIT, 5.38, 5.38, MAX_HIT, MAX_HIT, MAX_HIT, 6.00, 6.00,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, 6.00, 6.00, 6.00, 6.00, 7.00, 7.00, MAX_HIT, MAX_HIT, MAX_HIT,
										  7.00, MAX_HIT, MAX_HIT, MAX_HIT, 5.40, 6.00, 5.50, 6.00, 4.00, 3.00,
										  3.00, 0.00, 6.30, 6.50, 3.00, 3.00, 6.00, 6.00, 6.50, 7.00,
										  7.00, 3.00, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, MAX_HIT,
										  0.00};
	static double h3  [] = { 0., 7.00, MAX_HIT, 4.00, 3.50, 3.50, 4.00, 3.50, MAX_HIT, MAX_HIT, MAX_HIT, 
										  3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 6.20, MAX_HIT, MAX_HIT, 4.00, 
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 6.50, 6.50,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, 7.00, 7.00, 7.00, 7.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 6.00, MAX_HIT, 3.00,
										  3.00, MAX_HIT, MAX_HIT, 6.50, 3.00, 3.00, 6.50, 6.50, MAX_HIT, MAX_HIT,
										  MAX_HIT, 3.00, MAX_HIT, 3.60, 3.60, MAX_HIT, 3.60, MAX_HIT, 3.60, MAX_HIT,
										  MAX_HIT};
	static double h4  [] = { 0., 9.30, 0.00, 3.00, 7.45, 6.50, 8.45, 5.625, 4.40, 4.40, 4.40, 
										  4.40, 4.40, 4.40, 4.40, 4.40, 4.40, 8.20,  2.00, 5.40, 5.40, 
										  MAX_HIT, 4.45, 4.45, 2.74, 2.14, MAX_HIT,  4.45, 4.45, 5.80, 5.80,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 4.45, 4.45, MAX_HIT, 4.0, 4.50,
										  4.50, 7.00, 7.00, 7.00, 7.00, 0.00, 0.00, 6.80, 6.80, 6.80,
										  0.00, MAX_HIT, 4.50, 3.80, 3.30, 3.00, 2.50, 4.50, 6.20, 7.40,
										  7.40, 0.00, 0.00, 6.00, 2.15, 2.15, 7.20, 7.20, 4.00, 4.00,
										  4.00, 1.75, 4.40, 4.40, 4.40, 4.40, 4.40, 4.40, 4.40, 6.80,
										  0.00};
	static double a1  [] = { 0., 6.00, 6.60, 3.50, 3.50, 3.10, 5.00, 3.00, 4.50, 4.50, 4.50, 
										  2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 5.40, 2.00, 3.10, 3.10, 
										  7.60, 7.60, 7.60, 3.00, 3.00, 9.05, 9.05, 9.05, 4.40, 7.60,
										  7.60, 7.60, 9.05, 9.05, 9.05, 7.60, 9.05, 9.05, 3.20, 3.20,
										  3.20, 5.70, 5.70, 5.70, 5.70, 6.20, 6.20, 8.30, 8.30, 8.30,
										  6.20, 7.60, 5.00, 3.20, 2.50, 2.10, 2.80, 2.10, 5.00, 2.75,
										  2.75, 2.50, 4.20, 4.00, 2.00, 2.00, 5.40, 5.40, 4.20, 6.20,
										  6.20, 2.00, 4.50, 2.60, 2.60, 4.50, 2.60, 4.50, 2.60, 8.30,
										  10.00};
	static double a2  [] = { 0., 8.90, 4.00, 5.00, 3.90, 3.45, 3.35, 3.35, MAX_HIT, MAX_HIT, MAX_HIT, 
										  4.50, 4.50, 4.50, 4.50, 4.50, 4.50, 7.40, MAX_HIT, MAX_HIT, 4.60, 
										  4.00, 4.00, 4.00, MAX_HIT, MAX_HIT, 4.40, 4.40, 4.40, 6.40, 4.00,
										  4.00, 4.00, 4.40, 4.40, 4.40, 4.00, 4.40, 4.40, 4.20, 4.20,
										  4.20, 6.40, 6.40, 6.40, 6.40, 3.00, 3.00, 4.00, 4.00, 4.00,
										  3.00, 4.00, 6.50, 4.20, MAX_HIT, MAX_HIT, MAX_HIT, 4.20, MAX_HIT, 3.50,
										  3.50, 0.00, 0.00, 6.00, 3.50, 3.50, 7.80, 7.80, MAX_HIT, MAX_HIT,
										  MAX_HIT, 2.75, MAX_HIT, 4.50, 4.50, MAX_HIT, 4.50, MAX_HIT, 4.50, 4.00,
										  10.00};
	static double a3  [] = { 0., 5.10, MAX_HIT, 3.50, 4.30, 3.80, 5.00, 3.70, MAX_HIT, MAX_HIT, MAX_HIT, 
										  2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 4.80, MAX_HIT, MAX_HIT, 3.10, 
										  MAX_HIT, 4.60, 4.60, MAX_HIT, MAX_HIT, MAX_HIT, 5.80, 5.80, 3.80, 3.80,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 5.80, MAX_HIT, MAX_HIT, 5.50,
										  5.50, 5.10, 5.10, 5.10, 5.10, MAX_HIT, MAX_HIT, 6.00, 6.00, 6.00,
										  MAX_HIT, MAX_HIT, 7.20, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.10, MAX_HIT, 4.25,
										  4.25, MAX_HIT, MAX_HIT, 4.00, 2.00, 2.00, 4.80, 4.80, MAX_HIT, MAX_HIT,
										  MAX_HIT, 2.00, MAX_HIT, 2.50, 2.50, MAX_HIT, 2.50, MAX_HIT, 2.50, 6.00,
										  MAX_HIT};
	static double b1  [] = { 0., 6.00, 6.60, 3.50, 3.50, 3.10, 5.00, 3.00, 4.50, 4.50, 4.50, 
										  2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 5.40, 3.50, 4.60, 3.10, 
										  7.60, 7.60, 7.60, 5.00, 4.50, 9.05, 9.05, 9.05, 4.40, 7.60,
										  7.60, 7.60, 9.05, 9.05, 9.05, 7.60, 9.05, 9.05, 3.20, 3.20,
										  3.20, 5.70, 5.70, 5.70, 5.70, 6.20, 6.20, 8.30, 8.30, 8.30,
										  6.20, 7.60, 5.00, 3.20, 4.50, 4.20, 4.80, 2.10, 5.00, 2.75,
										  2.75, 2.50, 4.20, 4.00, 2.00, 2.00, 5.40, 5.40, 5.80, 6.20,
										  6.20, 2.00, 4.50, 2.60, 2.60, 4.50, 2.60, 4.50, 2.60, 8.30,
										 -5.00};
	static double b2  [] = { 0., 8.90, 5.50, 5.00, 3.90, 3.45, 3.35, 3.35, 2.50, 2.50, 2.50, 
										  4.50, 4.50, 4.50, 4.50, 4.50, 4.50, 7.40, 2.00, 3.10, 4.60, 
										  4.00, 4.00, 4.00, 3.00, 3.00, 4.40, 4.40, 4.40, 6.40, 4.00,
										  4.00, 4.00, 4.40, 4.40, 4.40, 4.00, 4.40, 4.40, 4.20, 4.20,
										  4.20, 6.40, 6.40, 6.40, 6.40, 5.50, 5.50, 4.00, 4.00, 4.00,
										  5.50, 4.00, 6.50, 4.20, 2.50, 2.10, 2.80, 4.20, 3.50, 3.50,
										  3.50, MAX_HIT, MAX_HIT, 6.00, 3.50, 3.50, 7.80, 7.80, 3.80, 5.50,
										  5.50, 2.75, 2.50, 4.50, 4.50, 2.50, 4.50, 2.50, 4.50, 4.00,
										 -5.00};
	static double b3  [] = { 0., 5.10, MAX_HIT, 3.50, 4.30, 3.80, 5.00, 3.70, MAX_HIT, MAX_HIT, MAX_HIT, 
										  2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 4.80, MAX_HIT, MAX_HIT, 3.10, 
										  MAX_HIT, 4.60, 4.60, MAX_HIT, MAX_HIT, MAX_HIT, 5.80, 5.80, 3.80, 3.80,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 5.80, MAX_HIT, MAX_HIT, 5.50,
										  5.50, 5.10, 5.10, 5.10, 5.10, MAX_HIT, MAX_HIT, 6.00, 6.00, 6.00,
										  MAX_HIT, MAX_HIT, 7.20, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.10, MAX_HIT, 4.25,
										  4.25, MAX_HIT, MAX_HIT, 4.00, 2.00, 2.00, 4.80, 4.80, MAX_HIT, MAX_HIT,
										  MAX_HIT, 2.00, MAX_HIT, 2.50, 2.50, MAX_HIT, 2.50, MAX_HIT, 2.50, 6.00,
										  MAX_HIT};
	static double c1  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										  6.50, 6.50, 6.50, MAX_HIT, MAX_HIT, 6.90, 6.90, 6.90, MAX_HIT, 6.50,
										  6.50, 6.50, 6.90, 6.90, 6.90, 6.50, 6.90, 6.90, 4.20, 6.75,
										  6.75, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, 6.50, 6.50, 4.20, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, MAX_HIT, 5.306,6.16, 6.16, 5.306,6.16, 6.16, 7.54, MAX_HIT};
	static double c2  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										  6.50, 6.50, 6.50, MAX_HIT, MAX_HIT, 6.90, 6.90, 6.90, MAX_HIT, 6.50,
										  6.50, 6.50, 6.90, 6.90, 6.90, 6.50, 6.90, 6.90, 4.20, 6.75,
										  6.75, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, 6.50, 6.50, 4.20, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT, MAX_HIT, 5.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT,
										  MAX_HIT};
	static int N_sheme[] = { 0, 
							_TOPO_TWO_CIRCLE_TREE,		_TOPO_LEFT_GROZO_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		 _TOPO_ONE_CIRCLE_TREE,		  _TOPO_ONE_CIRCLE_TREE,
							_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		 _TOPO_ONE_CIRCLE_TREE,		  _TOPO_TWO_CIRCLE_TREE,
							_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B, _TOPO_TWO_CIRCLE_TREE,		  _TOPO_ONE_CIRCUIT_DUPLEX_A, 
							_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_A, _TOPO_ONE_CIRCUIT_DUPLEX_B, _TOPO_ONE_CIRCUIT_DUPLEX_B,
							_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_LEFT_GROZO_TREE,		_TOPO_LEFT_GROZO_TREE,		_TOPO_HORIZONTAL_CROSS_TREE,_TOPO_HORIZONTAL_CROSS_TREE,_TOPO_HORIZONTAL_CROSS_TREE,		
							_TOPO_LEFT_GROZO_TREE,		_TOPO_ONE_CIRCUIT_DUPLEX_A,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		 _TOPO_ONE_CIRCLE_TREE,		  _TOPO_TWO_CIRCLE_TREE,
							_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		 _TOPO_ONE_CIRCLE_TREE,		  _TOPO_ONE_CIRCLE_TREE,
							_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,	 	 _TOPO_TWO_CIRCLE_TREE,		  _TOPO_HORIZONTAL_CROSS_TREE,
							_TOPO_ORDINARY_DUPLEX,
	};	
	static char * S0[] = { "У330-2", "У220-3", "УС110-8", "УШЛБ-12", "УТЛБ-12",	"ТЛБ-12", "ПЛБ-12",		"У1М",      "У3М",		"У5М", 
								  "У2М-2",  "У6М-2",  "У2М",     "У4М",     "УМ",			"У6М",    "У39",			"ПБ110-15", "У110-3",	"У110-4",
								  "П21Ма",	"П21",	 "П21М",		"П220а",	  "П220",		"П22Ма",  "П22",			"П22М",		"П26М",		"П21М-1",	
								  "П21М-2",	"П21М-3", "П22М-1",	"П22М-2",  "П22М-3",		"СП21М",  "СП22М",		"СП22Ма",	"ПБ330-7н", "ПБ16",
								  "ПУ30",	"У32",	 "У32М",		"У38",	  "У38М",		"У33",	  "У33М",		"У34",		"У34М",	   "У35",	 
								  "CУ33",	"П21Ма-1","ПБ500-5н","ПВС500-2","ПЕ110-6+НА","ПЕ110-5Н","ПБ220-3и",	"П110-6Н",	"УС110-3",	"УТЛБ-10", 
								  "ПЛБ-10",	"СН220",	 "КТС",		"ПС220-6т","ПБ26",		"ПБ28",	  "ЦП28",		"ЦП28М",		"ЦП23",		"ЦУ37-2", 
								  "ЦУ37",	"ДП110-2", "У1",	   "У2",      "У2к",       "У3",      "У4",        "У5",       "У6",			"У35п",
								  "портал"};

	Table * table = GetTowerDataBaseTable(81, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = h1;
		table->table[1].parm_values = h2;
		table->table[2].parm_values = h3;
		table->table[3].parm_values = h4;
		table->table[4].parm_values = a1;
		table->table[5].parm_values = a2;
		table->table[6].parm_values = a3;
		table->table[7].parm_values = b1;
		table->table[8].parm_values = b2;
		table->table[9].parm_values = b3;
		table->table[10].parm_values = c1;
		table->table[11].parm_values = c2;
		table->table[12].parm_values = N_sheme;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = h1[k];
		((double *)(table->table[1].parm_values))[k] = h2[k];
		((double *)(table->table[2].parm_values))[k] = h3[k];
		((double *)(table->table[3].parm_values))[k] = h4[k];
		((double *)(table->table[4].parm_values))[k] = a1[k];
		((double *)(table->table[5].parm_values))[k] = a2[k];
		((double *)(table->table[6].parm_values))[k] = a3[k];
		((double *)(table->table[7].parm_values))[k] = b1[k];
		((double *)(table->table[8].parm_values))[k] = b2[k]; 
		((double *)(table->table[9].parm_values))[k] = b3[k];
		((double *)(table->table[10].parm_values))[k] = c1[k];
		((double *)(table->table[11].parm_values))[k] = c2[k];
		((int    *)(table->table[12].parm_values))[k] = N_sheme[k];
	}
	return table;
}

//////////////////////////////////////////
//...таблицы геометрических размеров опор;
Table * tower_UB(int id_static_char)
{
	static double h1  [] = { 0., 10.30, 8.00, 13.50, 10.00, 12.50, 16.20, 8.50, MAX_HIT};
	static double h2  [] = { 0., 3.00, MAX_HIT, 3.00, 4.00, 4.00, 4.00, 4.00, MAX_HIT};
	static double h3  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT};
	static double h4  [] = { 0., 2.35, MAX_HIT, 2.35, 2.00, 2.50, 2.70, 6.70, MAX_HIT};
	static double a1  [] = { 0., 1.60, MAX_HIT, 1.60, 1.75, 3.10, 3.10, 3.10, MAX_HIT};
	static double a2  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT};
	static double a3  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT};
	static double b1  [] = { 0., 1.60, 1.00, 1.60, 1.75, 3.10, 3.10, 3.10, MAX_HIT};
	static double b2  [] = { 0., 1.00, MAX_HIT, 1.00, 1.75, 2.30, 2.30, 2.30, 7.00};
	static double b3  [] = { 0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT};
	static double c1  [] = { 0., MAX_HIT, MAX_HIT, 6.20, MAX_HIT, 6.10, 7.50, 5.10, MAX_HIT};
	static double c2  [] = { 0., MAX_HIT, MAX_HIT, 9.80, 6.00, 6.10, 7.50, 5.10, 10.1};
	static int N_sheme[] = { 0, 
							_TOPO_ONE_CIRCLE_TREE, _TOPO_SIMPLE_TRIPLEX, _TOPO_ONE_CIRCLE_TREE, _TOPO_ONE_CIRCLE_TREE, _TOPO_ONE_CIRCLE_TREE, _TOPO_ONE_CIRCLE_TREE, _TOPO_ONE_CIRCLE_TREE, _TOPO_ONE_CIRCUIT_TRIPLEX
	};
	static char * S0[] = { "УБ35-1В", "УБ35-3В(???)", "УСБ35-1В", "УБ35-1", "УБ110-1", "УСБ110-1", "УСБ110-3", "УБ500-1(???)"};
	Table * table = GetTowerDataBaseTable(8, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = h1;
		table->table[1].parm_values = h2;
		table->table[2].parm_values = h3;
		table->table[3].parm_values = h4;
		table->table[4].parm_values = a1;
		table->table[5].parm_values = a2;
		table->table[6].parm_values = a3;
		table->table[7].parm_values = b1;
		table->table[8].parm_values = b2;
		table->table[9].parm_values = b3;
		table->table[10].parm_values = c1;
		table->table[11].parm_values = c2;
		table->table[12].parm_values = N_sheme;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = h1[k];
		((double *)(table->table[1].parm_values))[k] = h2[k];
		((double *)(table->table[2].parm_values))[k] = h3[k];
		((double *)(table->table[3].parm_values))[k] = h4[k];
		((double *)(table->table[4].parm_values))[k] = a1[k];
		((double *)(table->table[5].parm_values))[k] = a2[k];
		((double *)(table->table[6].parm_values))[k] = a3[k];
		((double *)(table->table[7].parm_values))[k] = b1[k];
		((double *)(table->table[8].parm_values))[k] = b2[k]; 
		((double *)(table->table[9].parm_values))[k] = b3[k];
		((double *)(table->table[10].parm_values))[k] = c1[k];
		((double *)(table->table[11].parm_values))[k] = c2[k];
		((int    *)(table->table[12].parm_values))[k] = N_sheme[k];
	}
	return table;

}

Table * tower_PB(int id_static_char)
{
	static double h1  [] = {0., 10.80, 10.40, 10.30, 12.40, 11.00, 10.30, 15.50, 14.50, 12.50, 10.50, 
				 						 14.50, 14.50, 14.50, 14.50, 13.50, 13.50, 11.50, 13.50, 14.50, 15.50, 
										 13.50, 14.50, 16.00, 17.50, 19.50, 22.90, 17.50, 17.50, 12.50, 18.50, 
										 12.50, 23.00};
	static double h2  [] = {0., 2.50, 2.50, 3.00, 3.00, 3.50, 3.00, 3.00, 4.00, 3.00, 4.00, 
				 						 4.00, 2.00/*3.00???*/, 3.00, 4.00, 3.00, 3.00, 4.00, 4.00, 3.00, 3.00, 
										 4.00, 7.00, 5.50, 5.50, 4.50, 4.10, 3.00, 3.00, 4.00, 3.00, 
										 3.00, MAX_HIT};
	static double h3  [] = {0., MAX_HIT, 2.50,    MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 3.00,    4.00, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 3.00,    3.00,    4.00,    4.00,    3.00,    3.00, 
										 2.00,    MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 2.50, MAX_HIT};
	static double h4  [] = {0., 2.35, 0.50, 2.35, 0.50,  0.50, 0.00, 2.80, 2.80, 3.70, 3.70, 
				 						 1.10, 2.00, 2.00, 2.00, 2.70,    3.00,    2.70,    3.00,    4.00, 3.00, 
										 3.00, 2.50, 2.50, 2.50, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 4.00, 3.00, 
										 1.20, MAX_HIT};
	static double a1  [] = {0., 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.75, 1.75, 1.75, 
				 						 1.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 
										 2.50, 2.80, 2.80, 2.80, 4.20, 4.20, 2.50, 4.80, 2.00, 2.00, 
										 3.10, 5.30};
	static double a2  [] = {0., MAX_HIT, 1.70, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.50, 2.50, 
				 						 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 3.50, 3.50, 3.50, 3.50, 3.50, 4.00, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, 5.60};
	static double a3  [] = {0., MAX_HIT, 1.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 1.75, 1.75, 
				 						 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double b1  [] = {0., 1.70, 1.00, 1.70, 1.70, 1.70, 2.00, 2.50, 1.75, 1.75, 1.75, 
				 						 2.50, 3.50, 3.50, 3.50, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 
										 4.00, 4.80, 4.80, 4.80, 4.20, 4.20, 2.50, 2.80, 4.00, 4.00, 
										 3.10, 5.30};
	static double b2  [] = {0., 1.00, 1.70, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 2.50, 2.50, 
				 						 1.75, 2.00, 2.00, 2.00, 3.50, 3.50, 3.50, 3.50, 3.50, 4.00, 
										 2.50, 2.80, 2.80, 2.80, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.80, 2.00, 
										 2.30, 5.60};
	static double b3  [] = {0., MAX_HIT, 1.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 1.75, 1.75, 
				 						 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double c1  [] = {0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
				 						 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 8.40, 8.40, 5.32, 5.96, MAX_HIT, MAX_HIT, 
										 9.00, 7.60};
	static double c2  [] = {0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 5.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
				 						 6.50, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 6.50, MAX_HIT, 
										 9.00, 7.60};
	static int N_sheme[] = { 0, 
							_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCLE_TREE,	 _TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,	 _TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,	_TOPO_TWO_CIRCLE_TREE,
							_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCLE_TREE,	_TOPO_TWO_CIRCLE_TREE,	 _TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,	 _TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,	_TOPO_TWO_CIRCLE_TREE,
							_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCUIT_DUPLEX,_TOPO_ONE_CIRCUIT_DUPLEX,	_TOPO_ONE_CIRCUIT_DUPLEX,_TOPO_ONE_CIRCUIT_DUPLEX,	_TOPO_ONE_CIRCLE_TREE,	_TOPO_ONE_CIRCLE_TREE,
							_TOPO_ONE_CIRCUIT_TRIPLEX, _TOPO_ONE_CIRCUIT_DUPLEX_A
	};
	static char * S0[] = { "ПБ35-1В",  "ПБ35-2В",				  "ПБ35-3В", "ПБ35-5В", "ПБ35-7В", "ПУСБ35-1В", "ПБ35-1",   "ПБ35-3",					"ПБ35-2",    "ПБ35-4",
								  "ПУСБ35",   "ПБ110-1",				  "ПБ110-3", "ПБ110-5", "ПБ110-2", "ПБ110-4",   "ПБ110-6",  "ПБ110-8(III-IV р.г.)", "ПБ110-8",   "ПБ110-10",
								  "ПБ150-1",  "ПБ220-1(III-IV р.г.)", "ПБ220-1", "ПБ220-3", "ПБ330-1", "ПБ330-3",   "ПСБ150-1", "ПСБ220-1",					"ПУСБ110-1", "ПСБ110-1", 
								  "КСБ110-1(???)", "ПБ500-1"};
	Table * table = GetTowerDataBaseTable(32, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = h1;
		table->table[1].parm_values = h2;
		table->table[2].parm_values = h3;
		table->table[3].parm_values = h4;
		table->table[4].parm_values = a1;
		table->table[5].parm_values = a2;
		table->table[6].parm_values = a3;
		table->table[7].parm_values = b1;
		table->table[8].parm_values = b2;
		table->table[9].parm_values = b3;
		table->table[10].parm_values = c1;
		table->table[11].parm_values = c2;
		table->table[12].parm_values = N_sheme;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = h1[k];
		((double *)(table->table[1].parm_values))[k] = h2[k];
		((double *)(table->table[2].parm_values))[k] = h3[k];
		((double *)(table->table[3].parm_values))[k] = h4[k];
		((double *)(table->table[4].parm_values))[k] = a1[k];
		((double *)(table->table[5].parm_values))[k] = a2[k];
		((double *)(table->table[6].parm_values))[k] = a3[k];
		((double *)(table->table[7].parm_values))[k] = b1[k];
		((double *)(table->table[8].parm_values))[k] = b2[k]; 
		((double *)(table->table[9].parm_values))[k] = b3[k];
		((double *)(table->table[10].parm_values))[k] = c1[k];
		((double *)(table->table[11].parm_values))[k] = c2[k];
		((int    *)(table->table[12].parm_values))[k] = N_sheme[k];
	}
	return table;
}

Table * tower_US(int id_static_char)
{
	static double h1  [] = {0., 10.00, 10.50, 10.50, 10.50, 10.50, 10.50, 15.50, 15.50, 10.70, 10.70, 
				 						 10.70, 19.70, 10.50, 10.50, 10.50, 15.50, 15.50, 10.50, 10.50, 0.60, 
										 0.60,  0.60};
	static double h2  [] = {0., 3.00, 3.00, 4.00, 4.00, 6.50, 6.50, 6.50, 6.50, 7.00, 8.60, 
				 						 6.50, 6.50, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 8.10, 10.0, 
										 10.0, 10.00};
	static double h3  [] = {0., MAX_HIT, 3.00, MAX_HIT, 4.00, MAX_HIT, 6.50, MAX_HIT, 6.50, MAX_HIT, MAX_HIT, 
				 						 7.00, 8.60, MAX_HIT, 4.00, MAX_HIT, MAX_HIT, 4.00, 4.00,    MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double h4  [] = {0., 3.95, 3.95, 6.20, 6.20, 8.10, 8.10, 8.10, 8.10, 9.30, 0.00, 
				 						 9.30, 9.50, 5.40, 5.40, 6.20, 6.20, 6.20, 6.20, 0.00, 16.00, 
										 25.0, 31.0};
	static double a1  [] = {0., 2.80, 2.80, 3.50, 3.50, 6.60, 4.60, 6.60, 5.90, 8.00, 8.00, 
				 						 6.00/*8.00-длина траверсы*/, 8.00, 3.10, 3.10, 5.00, 3.50, 3.50, MAX_HIT, 6.60, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double a2  [] = {0., MAX_HIT, 3.50, MAX_HIT, 5.00, MAX_HIT, 6.60, MAX_HIT, 6.60, MAX_HIT, 4.30, 
				 						 8.90, 11.00, MAX_HIT, 4.60, MAX_HIT, MAX_HIT, 5.00, 5.00, 4.00, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double a3  [] = {0., MAX_HIT, 2.80, MAX_HIT, 3.50, MAX_HIT, 4.60, MAX_HIT, 5.90, MAX_HIT, MAX_HIT, 
				 						 5.10/*7.10-длина траверсы*/, 7.10, MAX_HIT, 3.10, MAX_HIT, MAX_HIT, 3.50, 3.50, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double b1  [] = {0., 3.50, 2.80, 5.00, 3.50, 6.60, 4.60, 6.60, 5.90, 8.00, 8.00, 
				 						 6.00/*8.00-длина траверсы*/, 8.00, 4.60, 3.10, 5.00, 5.00, 3.50, 3.50, 6.60, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double b2  [] = {0., 2.80, 3.50, 3.50, 5.00, 4.60/*5.90-длина траверсы*/, 6.60, 5.90, 6.60, 5.10/*7.10-длина траверсы*/, 7.00, 
				 						 8.90, 11.00, 3.10, 4.60, 3.50, 3.50, 5.00, 5.00, 5.50, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double b3  [] = {0., MAX_HIT, 2.80, MAX_HIT, 3.50, MAX_HIT, 4.60, MAX_HIT, 5.90, MAX_HIT, MAX_HIT, 
				 						 5.10/*7.10-длина траверсы*/, 7.10, MAX_HIT, 3.10, MAX_HIT, MAX_HIT, 3.50, 3.50, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double c1  [] = {0., 4.20, 4.20, 4.80, 4.80, 5.20, 5.20, 4.10, 4.10, 6.24, 6.24, 
				 						 6.85, 9.55, 4.10, 4.10, 4.80, 3.50, 3.50, 4.80, 5.20, MAX_HIT, 
										 MAX_HIT, MAX_HIT};
	static double c2  [] = {0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
				 						 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 20.00, 
										 25.00, 20.00};
	static int N_sheme[] = { 0,
							_TOPO_ONE_CIRCLE_TREE, _TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE, _TOPO_TWO_CIRCLE_TREE, _TOPO_ONE_CIRCLE_TREE,		 _TOPO_TWO_CIRCLE_TREE,		 _TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE, _TOPO_LEFT_GROZO_TREE, 
				 		   _TOPO_TWO_CIRCLE_TREE, _TOPO_ONE_CIRCUIT_DUPLEX_B, _TOPO_ONE_CIRCLE_TREE, _TOPO_TWO_CIRCLE_TREE, _TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_TRIPLEX,_TOPO_ONE_CIRCUIT_TRIPLEX, _TOPO_LEFT_GROZO_TREE, _TOPO_ORDINARY_TRIPLEX, 
							_TOPO_ORDINARY_TRIPLEX,    _TOPO_ORDINARY_TRIPLEX
	};
	static char * S0[] = { "У35-1т",  "У35-2т", "У110-1",  "У110-2",  "У220-1",  "У220-2",  "УС220-5", "УС220-6", "У330-1", "У330-3",
								  "У330-2", "УС330-2(???)", "У110-3Н", "У110-4Н", "УС110-3(???)", "УС110-5(???)", "УС110-6(???)", "УС110-7(???)", "У220-3",      "АУ16(???)",	   
								  "АУ25(???)",    "АУ30(???)"};
	Table * table = GetTowerDataBaseTable(22, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = h1;
		table->table[1].parm_values = h2;
		table->table[2].parm_values = h3;
		table->table[3].parm_values = h4;
		table->table[4].parm_values = a1;
		table->table[5].parm_values = a2;
		table->table[6].parm_values = a3;
		table->table[7].parm_values = b1;
		table->table[8].parm_values = b2;
		table->table[9].parm_values = b3;
		table->table[10].parm_values = c1;
		table->table[11].parm_values = c2;
		table->table[12].parm_values = N_sheme;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = h1[k];
		((double *)(table->table[1].parm_values))[k] = h2[k];
		((double *)(table->table[2].parm_values))[k] = h3[k];
		((double *)(table->table[3].parm_values))[k] = h4[k];
		((double *)(table->table[4].parm_values))[k] = a1[k];
		((double *)(table->table[5].parm_values))[k] = a2[k];
		((double *)(table->table[6].parm_values))[k] = a3[k];
		((double *)(table->table[7].parm_values))[k] = b1[k];
		((double *)(table->table[8].parm_values))[k] = b2[k]; 
		((double *)(table->table[9].parm_values))[k] = b3[k];
		((double *)(table->table[10].parm_values))[k] = c1[k];
		((double *)(table->table[11].parm_values))[k] = c2[k];
		((int    *)(table->table[12].parm_values))[k] = N_sheme[k];
	}
	return table;
}

Table * tower_PUS(int id_static_char)
{
	static double h1  [] = {0., 15.00, 14.00, 19.00, 19.00, 19.00, 22.00, 19.00, 19.00, 19.00, 19.00, 
										 19.00, 16.50, 20.50, 17.50, 22.50, 25.50, 22.50, 22.50, 22.50, 25.50, 
										 22.50, 20.60, 17.50, 25.50, 22.50, 5.50, 11.00, 15.00, 15.00, 15.00, 
										 15.00, 17.00, 12.00, 19.00, 19.00, 22.00, 19.00, 19.00, 19.00, 25.50, 
										 25.50, 22.50, 27.00, 27.00, 27.00, 27.00, 27.00, 27.00, 27.00, 27.00,
										 27.00, 32.00, 32.00, 32.70, 32.00, 19.00, 19.00, 19.00, 19.00};
	static double h2  [] = {0., 3.00, 3.00, 4.00, 4.00, 6.00, 4.00, 4.00, 4.00, 6.00, 6.00, 
										 6.00, 6.50, 6.50, 6.50, 6.50, 6.50, 6.50, 6.50, 6.50, 7.50, 
										 6.50, 7.50, 6.50, 7.70, 7.70, MAX_HIT, 3.00, 4.00, 6.00, 4.00, 
										 6.00, 4.00, 4.00, 6.00, 6.00, 6.00, 4.00, 6.00, 6.00, 6.50, 
										 6.50, 6.50, 1.50, 1.50, 1.50, 1.50, 1.50, 15.00, 15.00, 1.80,
										 1.80, 1.00, 2.40, 1.14, 2.00, 3.60, 3.60, 5.40, 5.40};
	static double h3  [] = {0., MAX_HIT, 3.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 4.00, 4.00, 6.00, MAX_HIT, 
										 6.00, MAX_HIT, MAX_HIT, 6.50, MAX_HIT, MAX_HIT, 6.50, MAX_HIT, 6.50, MAX_HIT, 
										 7.50, MAX_HIT, 7.50, MAX_HIT, MAX_HIT, MAX_HIT, 3.00, MAX_HIT, MAX_HIT, 4.00, 
										 6.00, MAX_HIT, 4.00, MAX_HIT, 6.00, MAX_HIT, MAX_HIT, MAX_HIT, 6.00, MAX_HIT, 
										 MAX_HIT, 6.50, 3.50, 3.50, 3.50, 3.50, 3.50, 8.00, 8.00, 3.50,
										 3.50, 6.60, 5.20, MAX_HIT, MAX_HIT, 3.60, 3.60, 5.40, 5.40};
	static double h4  [] = {0., 2.00, 2.00, 2.00, 2.00, 3.00, 4.00, 4.00, 4.00, 4.00, 3.00, 
										 4.00, 4.00, 4.00, 5.50, 3.60, 4.00, 6.00, 3.20, 9.20, 4.70, 
										 7.00, 4.70, 7.00, 5.30, 6.80, 25.50, 2.00, 6.00, 4.00, 4.00, 
										 4.00, 4.00, MAX_HIT, 2.00, 3.00, 3.00, 2.00, 4.50, 4.50, 4.00, 
										 4.00, 5.50, 5.00, 5.00, 5.00, 5.00, 5.00, 6.00, 6.00, 5.30,
										 5.30, 7.60, 7.60, 6.44, 5.60, 3.30, 3.30, 3.30, 3.30};
	static double a1  [] = {0., 2.00, 2.00, 2.00, 2.10, 2.10, 5.20, 2.00, 2.10, 2.10, 2.60, 
										 2.60, 5.90, 3.90, 4.20, 4.00, 6.10, 4.00, 5.50, 4.00, 5.80, 
										 5.60, 5.80, 5.60, 6.00, 6.40, 5.75, 2.00, 2.10, 2.10, 2.10, 
										 2.10, 5.20, 2.10, 2.60, 2.60, 5.20, 2.10, 3.40, 3.40, 5.90, 
										 3.90, 4.20, 3.90, 3.90, 4.20, 4.20, 4.20, 6.20, 6.20, 4.95,
										 5.65, 4.90, 4.90, 19.50, 20.00, 2.50, 2.50, 2.50, 2.50};
	static double a2  [] = {0., MAX_HIT, 3.30, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 4.10, 4.20, 4.20, MAX_HIT, 
										 4.20, MAX_HIT, MAX_HIT, 6.40, MAX_HIT, MAX_HIT, 6.00, MAX_HIT, 6.00, MAX_HIT, 
										 8.80, MAX_HIT, 8.80, MAX_HIT, 9.60, 4.00, 3.30, MAX_HIT, MAX_HIT, 4.20, 
										 4.20, MAX_HIT, 4.20, MAX_HIT, 4.20, MAX_HIT, MAX_HIT, MAX_HIT, 4.60, MAX_HIT, 
										 MAX_HIT, 6.40, 8.10, 8.10, 8.60, 8.60, 8.60, 12.0, 12.0, 9.85, 
										 7.90, 18.50, 18.50, 19.50, 20.00, 3.50, 3.50, 4.50, 4.50};
	static double a3  [] = {0., MAX_HIT, 2.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.00, 2.10, 2.10, MAX_HIT, 
										 2.60, MAX_HIT, MAX_HIT, 3.50, MAX_HIT, MAX_HIT, 4.00, MAX_HIT, 4.00, MAX_HIT, 
										 4.90, MAX_HIT, 4.90, MAX_HIT, 5.60, 8.70, 2.00, MAX_HIT, MAX_HIT, 2.10, 
										 2.10, MAX_HIT, 2.10, MAX_HIT, 2.60, MAX_HIT, MAX_HIT, MAX_HIT, 3.40, MAX_HIT, 
										 MAX_HIT, 3.50/*4.50???*/, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 7.60, 7.60, MAX_HIT, 										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.50, 2.50, 2.50, 2.50};
	static double b1  [] = {0., 3.30, 2.00, 4.10, 4.20, 4.20, 5.20, 2.00, 2.10, 2.10, 4.20, 
										 2.60, 5.90, 6.10, 4.20, 6.00, 6.10, 4.00, 7.50, 5.50, 8.30, 
										 5.60, 8.30, 5.60, 9.60, 6.90, 5.75, 2.00, 4.20, 4.20, 2.10, 
										 2.10, 5.20, 2.10, 4.20, 2.60, 5.20, 4.20, 4.60, 3.40, 5.90, 
										 6.10, 4.20, 3.90, 3.90, 4.20, 4.20, 4.20, 6.20, 6.20, 4.95,
										 5.65, 4.90, 4.90, 19.50, 20.00, 2.50, 2.50, 2.50, 2.50};
	static double b2  [] = {0., 2.00, 3.30, 2.00, 2.10, 2.10, 2.60, 4.10, 4.20, 4.20, 2.60, 
										 4.20, 3.30, 3.50, 6.40, 4.00, 3.50, 6.00, 5.50, 7.50, 4.80, 
										 8.80, 4.80, 8.80, 5.60, 9.60, 4.00, 3.30, 2.10, 2.10, 4.20, 
										 4.20, 2.60, 2.10, 2.60, 4.20, 2.60, 2.10, 3.40, 4.60, 3.30, 
										 3.50, 6.40, 8.10, 8.10, 8.60, 8.60, 8.60, 12.0, 12.0, 9.85,
										 7.90, 18.50, 18.50, 19.50, 20.00, 3.50, 3.50, 4.50, 4.50};
	static double b3  [] = {0., MAX_HIT, 2.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.00, 2.10, 2.10, MAX_HIT, 
										 2.60, MAX_HIT, MAX_HIT, 3.50, MAX_HIT, MAX_HIT, 4.00, MAX_HIT, 5.50, MAX_HIT, 
										 4.90, MAX_HIT, 4.90, MAX_HIT, 5.60, 8.70, 2.00, MAX_HIT, MAX_HIT, 2.10, 
										 2.10, MAX_HIT, 4.20, MAX_HIT, 2.60, MAX_HIT, MAX_HIT, MAX_HIT, 3.40, MAX_HIT, 
										 MAX_HIT, 3.50, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 7.60, 7.60, MAX_HIT, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 2.50, 2.50, 2.50, 2.50};
	static double c1  [] = {0., 1.80, 1.80, 2.50, 2.80, 2.80, 6.00, 2.50, 2.80, 2.80, 2.80, 
										 2.80, 5.00, 4.42, 4.82, 2.10, 7.00, 4.10, 4.55, 4.55, 5.42, 
										 5.75, 4.82, 5.17, 5.33, 5.75, 5.42, 1.50, 2.40, 2.40, 2.40, 
										 2.40, 4.80, 2.10, 2.80, 2.75, 6.00, MAX_HIT, 3.30, 3.30, 7.00, 
										 5.00, 5.40, 8.70, 8.70, 9.20, 9.20, 9.20, 12.00, 12.00, 10.40, 
										 8.00, 13.60, 13.60, 4.77, 4.80,2.908,2.908,2.908,2.908};
	static double c2  [] = {0., MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 12.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, 10.00, MAX_HIT, MAX_HIT, MAX_HIT, 14.00, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, MAX_HIT, 
										 MAX_HIT, 9.60, 2.06, MAX_HIT, 2.75, 12.00, MAX_HIT, MAX_HIT, MAX_HIT, 14.00, 
										 MAX_HIT, MAX_HIT, 8.70, 8.70, 9.20, 9.20, 9.20, 8.00, 8.00, 10.40, 
										 11.90, 13.60, 13.60, MAX_HIT, MAX_HIT,2.408,2.408,2.408,2.408};
	static int N_sheme[] = { 0,
							_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE, 
				 		   _TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE, 
						   _TOPO_TWO_CIRCLE_TREE,		_TOPO_SPLIT_CIRCLE_TREE,	_TOPO_TWO_CIRCLE_TREE,		_TOPO_SPLIT_CIRCLE_TREE,	_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCLE_RUMKA,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE, 
						   _TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_TREE,		_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_SPLIT_CIRCLE_TREE,	_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCLE_RUMKA,		_TOPO_SPLIT_CIRCLE_TREE,	_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_TRIPLEX, _TOPO_ONE_CIRCLE_TREE, 
						   _TOPO_ONE_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCLE_RUMKA,		_TOPO_ONE_CIRCLE_RUMKA,		_TOPO_ONE_CIRCUIT_DUPLEX_B, 
						   _TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ONE_CIRCUIT_DUPLEX_B,_TOPO_ORDINARY_DUPLEX,		_TOPO_ORDINARY_DUPLEX,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE,		_TOPO_TWO_CIRCLE_TREE
	};
	static char * S0[] = { "П35-1",   "П35-2",   "П110-1",  "П110-3",  "П110-5",   "П110-7",   "П110-2",   "П110-4",   "П110-6",   "П150-1",
								  "П150-2",  "ПС220-1", "ПС220-3", "ПС220-2", "ПС220-5",  "ПС220-7",  "ПС220-6",  "ПУС220-1", "ПУС220-2", "П330-3",
								  "П330-2",  "ПС330-3(???)", "ПС330-2", "ПС330-5(???)", "ПС330-6(???)",  "ПС330-7(???)",  "ПС35-2",   "ПС110-3",  "ПС110-5",  "ПС110-4",
								  "ПС110-6", "ПС110-7", "ПС35-4(???)",  "ПС110-9(???)", "ПС110-10", "ПС110-11(???)", "ПС110-13(???)", "ПУС110-1(???)", "ПУС110-2(???)", "П220-1",  
								  "П220-3",  "П220-2",  "ПБ1",     "ПБ2",		 "ПБ3",      "ПБ4",      "ПБ5",      "Р1",       "Р2",		 "ПУБ-2", 
								  "ПУБ-5",   "ПО",      "ПМО",     "ПП",	    "ПУ5",	    "П2",	    "П4",	    "П6",	    "П8",
	};
	Table * table = GetTowerDataBaseTable(59, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = h1;
		table->table[1].parm_values = h2;
		table->table[2].parm_values = h3;
		table->table[3].parm_values = h4;
		table->table[4].parm_values = a1;
		table->table[5].parm_values = a2;
		table->table[6].parm_values = a3;
		table->table[7].parm_values = b1;
		table->table[8].parm_values = b2;
		table->table[9].parm_values = b3;
		table->table[10].parm_values = c1;
		table->table[11].parm_values = c2;
		table->table[12].parm_values = N_sheme;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = h1[k];
		((double *)(table->table[1].parm_values))[k] = h2[k];
		((double *)(table->table[2].parm_values))[k] = h3[k];
		((double *)(table->table[3].parm_values))[k] = h4[k];
		((double *)(table->table[4].parm_values))[k] = a1[k];
		((double *)(table->table[5].parm_values))[k] = a2[k];
		((double *)(table->table[6].parm_values))[k] = a3[k];
		((double *)(table->table[7].parm_values))[k] = b1[k];
		((double *)(table->table[8].parm_values))[k] = b2[k]; 
		((double *)(table->table[9].parm_values))[k] = b3[k];
		((double *)(table->table[10].parm_values))[k] = c1[k];
		((double *)(table->table[11].parm_values))[k] = c2[k];
		((int    *)(table->table[12].parm_values))[k] = N_sheme[k];
	}
	return table;
}

/////////////////////////////////////////
//...функция заполнения параметров опоры;
void set_tower_param(void * context, int i, double * par)
{
	int N_circuit = -1, l; 
	if ((int)par[0] <  _NUM_TOWER_TOPOS && 
		 (int)par[0] != _TOPO_TWO_CIRCLE_TREE) N_circuit = 1; else
	if ((int)par[0] <  _NUM_TOWER_TOPOS && 
		 (int)par[0] == _TOPO_TWO_CIRCLE_TREE) N_circuit = 2;
	if ((int)par[0] == _TOPO_ONE_CIRCLE_TREE ||
		 (int)par[0] == _TOPO_TWO_CIRCLE_TREE) {
		par[1] =  GetTableParam(context, i, 0); //...h1
		par[2] = -GetTableParam(context, i, 4); //...a1
		par[3] =  GetTableParam(context, i, 7); //...b1
		for (l = 0; l < N_circuit; l++) {
			par[4+l*3] = par[4+(l-1)*3] + GetTableParam(context, i, 1+l); //...hl
			if ((par[5+l*3] = GetTableParam(context, i, 5+l)) == MAX_HIT) par[5+l*3] = 0.; else par[5+l*3] = -par[5+l*3];//...al
			if ((par[6+l*3] = GetTableParam(context, i, 8+l)) == MAX_HIT) par[6+l*3] = 0.;//...bl
		}
		par[4+l*3] = par[4+(l-1)*3] + GetTableParam(context, i, 3); //...h_last
		return;
	}
	if ((int)par[0] == _TOPO_LEFT_GROZO_TREE) {
		par[1] =  GetTableParam(context, i, 0); //...h1
		par[2] = -GetTableParam(context, i, 4); //...a1
		par[3] =  GetTableParam(context, i, 7); //...b1
		for (l = 0; l < N_circuit; l++) {
			par[4+l*3] = par[4+(l-1)*3] + GetTableParam(context, i, 1+l); //...hl
			if ((par[6+l*3] = GetTableParam(context, i, 8+l)) == MAX_HIT) par[6+l*3] = 0.;//...bl
		}
		par[4+l*3] = par[7+l*3] = par[4+(l-1)*3] + GetTableParam(context, i, 3); //...h_last
		if ((par[5+l*3] = GetTableParam(context, i, 4+l)) == MAX_HIT) par[5+l*3] = 0.; else par[5+l*3] = -par[5+l*3];//...a_last
		return;
	}
	if ((int)par[0] == _TOPO_SPLIT_CIRCLE_TREE) {
		par[1] = GetTableParam(context, i, 0); //...h1
		par[2] = -3.5; //...a1
		par[3] =  3.5; //...b1
		par[4] =  4.0+par[1]; //...hl
		par[5] = -5.0; //...al
		par[6] =  5.0; //...bl
		par[7] =  4.0+par[4];
		par[8] = -3.5;
		par[9] =  3.5;
		par[10] =  3.0+par[7];
		par[12] =  3.5;
		par[13] =  4.0+par[10];
		par[15] =  5.0;
		par[16] =  4.0+par[13];
		par[18] =  3.5;
		par[19] =  6.2+par[16]; //...h_last
		return;
	}
	if ((int)par[0] == _TOPO_HORIZONTAL_CROSS_TREE) {
		par[1] =  GetTableParam(context, i, 0); //...h1
		par[2] = -GetTableParam(context, i, 4); //...a1
		par[3] =  GetTableParam(context, i, 7); //...b1
		par[4] =  GetTableParam(context, i, 3)+par[1]; //...h4
		par[5] = -GetTableParam(context, i, 5); //...a2
		par[6] =  GetTableParam(context, i, 8); //...b2
		par[7] = -GetTableParam(context, i, 6); //...a3
		par[8] =  GetTableParam(context, i, 9); //...b3
		return;
	}
	if ((int)par[0] == _TOPO_ORDINARY_DUPLEX) {
		par[2] = par[5] = par[10] = -GetTableParam(context, i, 4)*.5; //...a1
		par[3] = par[6] = par[11] =  GetTableParam(context, i, 5)*.5; //...a2
		par[7] = -GetTableParam(context, i, 7)*.5; //...b1
		par[8] =  GetTableParam(context, i, 8)*.5; //...b2
		par[4] = GetTableParam(context, i, 0); //...h1
		par[9] = par[14] = GetTableParam(context, i, 1)+GetTableParam(context, i, 3)+par[4]; //...h2+h4
		return;
	}
	if ((int)par[0] == _TOPO_ONE_CIRCUIT_DUPLEX) {
		par[2] = -(par[3] = GetTableParam(context, i, 10)*.5); //...c1
		par[4] =  GetTableParam(context, i, 0); //...h1
		par[5] = -(par[6] = GetTableParam(context, i, 10)*.5); //...c1
		par[7] = -GetTableParam(context, i, 4); //...a1
		par[8] =  GetTableParam(context, i, 7); //...b1
		if ((par[9] =  GetTableParam(context, i, 3)) == MAX_HIT) par[9] = par[4]; 
		else par[9] += par[4]; //...h4
		par[14] =  (par[9] += GetTableParam(context, i, 1)); //...h2
		par[10] = -(par[11] = GetTableParam(context, i, 10)*.5); //...c1
		return;
	}
	if ((int)par[0] == _TOPO_ONE_CIRCUIT_DUPLEX_A) {
		par[2] = -GetTableParam(context, i, 10); //...c1
		par[3] =  GetTableParam(context, i, 11); //...c2
		par[4] = par[9]  = par[14] = GetTableParam(context, i, 0); //...h1
		par[5] = par[10] = -GetTableParam(context, i, 5); //...a2
		par[6] = par[11] =  GetTableParam(context, i, 8); //...b2
		par[7] = -GetTableParam(context, i, 4); //...a1
		par[8] =  GetTableParam(context, i, 7); //...b1
		return;
	}
	if ((int)par[0] == _TOPO_ONE_CIRCUIT_DUPLEX_B) {
		par[2] = -GetTableParam(context, i, 10); //...c1
		par[3] =  GetTableParam(context, i, 11); //...c2
		par[4] =  GetTableParam(context, i, 0);  //...h1
		par[5] = -GetTableParam(context, i, 5); //...a2
		par[6] =  GetTableParam(context, i, 8); //...b2
		par[7] = -GetTableParam(context, i, 4); //...a1
		par[8] =  GetTableParam(context, i, 7); //...b1
		par[9] = par[14] = GetTableParam(context, i, 3)+ par[4]; //...h4
		if ((par[10] = GetTableParam(context, i, 6)) == MAX_HIT) par[10] = par[5]; else par[10] = -par[10]; //...a3
		if ((par[11] = GetTableParam(context, i, 9)) == MAX_HIT) par[11] = par[6];//...b3
		return;
	}
	if ((int)par[0] == _TOPO_ONE_CIRCLE_RUMKA) {
		par[1] = GetTableParam(context, i, 1); //...h2
		par[2] = GetTableParam(context, i, 4); //...a1
		par[3] = GetTableParam(context, i, 2); //...h3
		par[4] = GetTableParam(context, i, 0); //...h1
		par[5] = -(par[6] = GetTableParam(context, i, 4)); //...а1
		par[7] = -(GetTableParam(context, i, 5)+(par[10] = -GetTableParam(context, i, 4))); //...a2-a1
		par[8] =  (GetTableParam(context, i, 8)-(par[11] =  GetTableParam(context, i, 7))); //...b2-b1
		par[12] =-(GetTableParam(context, i, 6)-GetTableParam(context, i, 4)); //...a3-a1
		par[13] = (GetTableParam(context, i, 9)-GetTableParam(context, i, 7)); //...b3-b1
		par[14] = par[9] = GetTableParam(context, i, 3)+par[4]; //...h4
		return;
	}
	if ((int)par[0] == _TOPO_ORDINARY_TRIPLEX) {
		return;
	}
}

////////////////////////////////
//...дополнительный список опор;
int tower_DB1(void * context, char * tower, double * par)
{
//////////////////////////////////////////////////////////////////////
//...ищем опору в дополнительном списке опор и заполняем ее параметры;
	char * sss;
	memset(par,  0, NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));

	char * str = strstr(tower, "+"), temp;
	if (str && (str[1] < '0' || str[1] > '9'))
	while ((str = strstr(str+1, "+")) != NULL && (str[1] < '0' || str[1] > '9')); //...пытаемся найти признак базы;
	
	if (str) {
		temp = str[0]; 
		str[0] = '\x0';
	}

	Table * table = (Table *)GetTable(context);
	for (int i = 1; table && i <= table->N_group; i++) //...сравниваем с учетом регистра;
	if ((sss = GetProfileName(context, i)) != NULL && tower) {
		if (! ::strcmp(tower, sss)) {
			par[0] = GetTableParam(context, i, 12);
			set_tower_param		 (context, i, par);

			if (str) str[0] = temp;
			return(1);
		}
	}
	if (str) str[0] = temp;
	return(0);
}

//////////////////////////
//...основной список опор;
int tower_DB(void ** tower_list, char * tower, double * par, int k_excl)
{
	int kV = strstr(tower, "110") ? 1 :
				strstr(tower, "220") ? 2 :
				strstr(tower, "330") ? 3 : 0;

	char * str = strstr(tower, "+"), temp;
	if (str && (str[1] < '0' || str[1] > '9'))
	while ((str = strstr(str+1, "+")) != NULL && (str[1] < '0' || str[1] > '9')); //...пытаемся найти признак базы;

	if (str) {
		temp = str[0]; 
		str[0] = '\x0';
	}

//////////////////////////////////////////////////////////////
//...ищем опору в списке базы данных и заполняем ее параметры;
	char * sss;
	memset(par,  0, NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));

	for (int k = 0; tower_list && tower_list[k]; k++) if (k != k_excl) { //...сравниваем без учета регистра;
		Table * table = (Table *)GetTable(tower_list[k]);
		for (int i = 1; table && i <= table->N_group; i++) 
			if ((sss = GetProfileName(tower_list[k], i)) != NULL && tower) {
				if (! _strnicmp(tower, sss, strlen(sss)) &&
//				if (! _stricmp (tower, sss) &&
					(! kV || kV == 1 && strstr(sss, "110") || kV == 2 && strstr(sss, "220") || kV == 3 && strstr(sss, "330"))) {

					par[0] = GetTableParam(tower_list[k], i, 12);
					set_tower_param		 (tower_list[k], i, par);

					if (str) str[0] = temp;
					return(1);
				}
			}
	}
	if (str) str[0] = temp;
	return(0);
}

///////////////////////////////
//...список нестандартных опор;
int tower_DB2(void * context, char * tower, double * par)
{
/////////////////////////////////////////////////////////////////////
//...ищем опору в списке нестандартных опор и заполняем ее параметры;
	char * sss;
	memset(par, 0, NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));

	char * str = strstr(tower, "+"), temp;
	if (str && (str[1] < '0' || str[1] > '9'))
	while ((str = strstr(str+1, "+")) != NULL && (str[1] < '0' || str[1] > '9')); //...пытаемся найти признак базы;

	if (str) {
		temp = str[0]; 
		str[0] = '\x0';
	}

	Table * table = (Table *)GetTable(context);
	for (int i = 1; table && i <= table->N_group; i++) //...сравниваем с учетом регистра;
	if ((sss = (char *)(long long)GetTableParam(context, i, 0)) == NULL && ! tower || sss && tower) {
		if (! sss && ! tower || sss && tower && ! ::strcmp(tower, sss)) {
			for ( int j = 1; j <= NUM_OF_TOWER_BASE_ELEMENTS; j++)
				par[j-1] = GetTableParam(table, i, j);

			if (str) str[0] = temp;
			return(1);
		}
	}
	if (str) str[0] = temp;
	return(0);
}

/////////////////////////////////////////////////////////////////////
// ...база данных стандартных опор (включая дополнительные проверки);
void TowerBase(void ** tower_list, char * tower, double * par, void * context)
{
	static double tower_par[][NUM_OF_TOWER_BASE_ELEMENTS] = {
		{_TOPO_TWO_CIRCLE_TREE, 20.2, -5.7, 5.7, 26.2, -6.4, 6.4, 33.2, -5.1, 5.1,  40.2 },					//[0]//...У32, У36(???)
		{_TOPO_ORDINARY_TRIPLEX, 2.0, -4.5, 8.5, 17.0, 2.0, -4.5, 8.5, 0., -6.5, 6.5, 0., 21.5, 0., 0.},//[1]//...У???(столбы)
		{_TOPO_ONE_CIRCLE_TREE, 12.7, -2.6, 4.5, 16.3, -0.0, 2.5, 20.7, -0.0, 0.0,  0.0 },					//[2]//...УБ110-5(???)
/////////////////////////////
//...специальная серия У и П;
		{_TOPO_ORDINARY_DUPLEX, 0., -5., 5., 7., -5., 5., 2.5, -2.5, 7., -5., 5., 0., 0., 7.},				//[3]//...ПС
		{_TOPO_SIMPLE_TREE, 10.5, 0., 0., 10.5, 0., 0., 10.5, 0., 0., 10.5 },									//[4]//...У
		{_TOPO_SIMPLE_DUPLEX, 0., 2.5, 2.5, 10.5, 2.5, 2.5, -2.5, -2.5, 10.5, 2.5, 2.5, 0., 0., 10.5 },	//[5]//...П
		{_TOPO_TWO_CIRCLE_TREE, 14.5, -10.0, 10.0, 14.5, -0.0, 0.0, 14.5, -0.0, 0.0,  14.5 },				//[6]//...консоль
	};
	double base = 0.;

/////////////////////////
//...параметры подставки;
	char * str, * sss, * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	setlocale(LC_COLLATE, "rus"); //...сравниваем без учета регистра;
	if ((sss = tower) != NULL)
		while ((str = strstr(sss, "+")) != NULL && strlen(str) > 1)
		if (str[1] < '0' || str[1] > '9') sss = str+1; else { //...все заменяем на десятичную точку;
			char * ppp = str;
			while ((ppp = strstr(ppp, ",")) != NULL) { ppp[0] = localeconv()->decimal_point[0]; ppp++;}

			ppp = str;
			while ((ppp = strstr(ppp, ".")) != NULL) { ppp[0] = localeconv()->decimal_point[0]; ppp++;}
			
			base += strtod(str+1, & sss);

			ppp = str;
			while ((ppp = strstr(ppp, localeconv()->decimal_point)) != NULL) { ppp[0] = '.'; ppp++;}
		}

///////////////////////////////////////////////////////////////////////
//...логика выбора соответствующих параметров опоры и способа подвески;
		if (tower) {
		if (! _strnicmp(tower, "У32", 3) || ! _strnicmp(tower, "У36", 3))
			memcpy(par, tower_par[0], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
		else
		if (! _strnicmp(tower, "СУ33", 4))
			memcpy(par, tower_par[1], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
		else
		if (! _strnicmp(tower, "УБ110-5", 7))
			memcpy(par, tower_par[2], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
		else
//////////////////////////////////////////////////////////
//...прибавляем специальные типы опор после всех проверок;
		if (! tower_DB2(context, tower, par) && ! tower_DB1(tower_list[0], tower, par) && ! tower_DB (tower_list, tower, par, 0)) {
			if (strstr(tower, "ПС"))
				memcpy(par, tower_par[3], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
			else
			if (! _strnicmp(tower, "У", 1))
				memcpy(par, tower_par[4], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
			else
			if (! _strnicmp(tower, "П", 1))
				memcpy(par, tower_par[5], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
			else
			if (! _strnicmp(tower, "консоль", 7))
				memcpy(par, tower_par[6], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
			else
				memcpy(par, tower_par[3], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
		}
	}
	else
	if (! tower)
//		memcpy(par, tower_par[0], NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));
		memset(par,            0, NUM_OF_TOWER_BASE_ELEMENTS*sizeof(double));

	setlocale(LC_COLLATE, sloc); free(sloc); 

///////////////////////////////////////////////
//...прибавляем параметры подставки к траверсе;
	switch ((int)par[0]) {
      case _TOPO_SIMPLE_TREE:
			par[1] += base;
		break;
      case _TOPO_ONE_CIRCLE_TREE:
      case _TOPO_LEFT_GROZO_TREE:
			par[1]  += base;
			par[4]  += base;
			par[7]  += base;
		break;
      case _TOPO_TWO_CIRCLE_TREE:
			par[1]  += base;
			par[4]  += base;
			par[7]  += base;
			par[10] += base;
		break;
 		case _TOPO_SPLIT_CIRCLE_TREE: 
			par[1]  += base;
			par[4]  += base;
			par[7]  += base;
			par[10] += base;
			par[13] += base;
			par[16] += base;
			par[19] += base;
		break;
		case _TOPO_HORIZONTAL_CROSS_TREE:
			par[1] += base;
			par[4] += base;
		break;
      case _TOPO_SIMPLE_DUPLEX:
			par[4] += base;
		break;
		case _TOPO_ORDINARY_DUPLEX:
		case _TOPO_ONE_CIRCUIT_DUPLEX:
		case _TOPO_ONE_CIRCUIT_DUPLEX_A:
		case _TOPO_ONE_CIRCUIT_DUPLEX_B:
		case _TOPO_ONE_CIRCLE_RUMKA:
			if (par[1]) 
			par[1] += base;
			par[4] += base;
			par[9] += base;
			if (par[14]) 
			par[14] += base;
		break;
		case _TOPO_SIMPLE_TRIPLEX:
		case _TOPO_ORDINARY_TRIPLEX:
		case _TOPO_ONE_CIRCUIT_TRIPLEX:
			par[4]  += base;
			par[12] += base;
		break;
	}
}

////////////////////////////////////////////////////////
//...функция образования контекста для базы данных опор;
void * CreateTOWERContext(int N_sm)
{
	if (N_sm < 0 || N_sm >= NUM_TOWER_SAMPLES) N_sm = 0;

	Context * cont = (Context *)new_struct(sizeof(Context));
	if (! cont) return(NULL);

	cont->N           = N_sm+SHIFT_TOWER_SAMPLES;
	cont->static_char = ADDITIONAL_STATE;
	cont->sample_name = GetTOWERSampleName(N_sm);
	cont->units       = UNIT_SI;

	switch (cont->N-SHIFT_TOWER_SAMPLES) {
      case _TOWER1: cont->table = tower_DOP(cont->static_char); break;
 		case _TOWER2: cont->table = tower_UB (cont->static_char); break;
      case _TOWER3: cont->table = tower_PB (cont->static_char); break;
      case _TOWER4: cont->table = tower_US (cont->static_char); break;
      case _TOWER5: cont->table = tower_PUS(cont->static_char); break;
  }
  set_default    (cont);
  set_table_units(cont->table, UNIT_SI);
  
  SetTableIndex(cont, 1);
  SetUnits     (cont, UNIT_SI); return (void *)cont;
}

/*
/*===============================================================================================*/
/*                                 ИНИЦИАЛИЗАЦИЯ ОБРАЗЦА                                         */
/*===============================================================================================*/
////////////////////////////////////////////////////////////////////
//...пpедваpительная инициализация сечений для базы данных проводов;
#ifndef ___ABRIDGE_PROFILE_MODE___
int tower3D_init(void * context, CGrid * block_nd)
{
	Context * cont = (Context *)context;
//	int  units_tab = cont->units;

	if (! is_Tower3D(cont)  || ! cont->table) return(0);
//	if (cont->table->index > 0) units_tab = cont->table->table_units[cont->table->index-1];

//	Table * tab = cont->table;

	cont->left = cont->right = cont->bottom = cont->top = cont->back = cont->front = 0.;

	DeleteSample(cont->sm);

	switch (cont->N-SHIFT_TOWER_SAMPLES) {
		case _TOWER1: 
		case _TOWER2:
		case _TOWER3:
		case _TOWER4: { 
		} break;
	}
	return(1);
//	return(cont->sm != NULL);
}
#endif

/*
////////////////////////////////////////
//...pешатель функций фоpмы для 3D опор;
int tower3D_solv(void * context, int & sec, int & hund, int id_solv)
{
    Context  * cont = (Context *)context; sec = hund = 0;
    if ((is_Skin2D(cont) || is_SkGp3D(cont) || is_SkIb3D(cont)) && cont->sm && ! cont->sm->solv) {
        CSample * sm = cont->sm;

///////////////////////////
//...solving of the probem;
        clock_t start = clock();
        int  counting = (cont->N == _NUM_SK2D_SAMPLES-1+SHIFT_SK2D_SAMPLES ? 
                         ANALYTICAL_COUNTING : 
                             BASIC_COUNTING);

        if (sm->counting_kernel(counting) != OK_COUNTING) {
            Message("Error in sample counting...");

            if (sm->solver)
                sm->solver->release_matrix();

            return(0);
        }
        Message("O'K");
        if (id_solv & MASK9) cont->sm->solv = 1;

////////////////////////////////////////////
//...фиксируем время и выходим из программы;
        timing_process(start, & hund, & sec);
    }
    return(1);
}
*/
/*
//////////////////////////////////////////////////////////
//...результаты для 2D образца (задачи пограничного слоя);
char * GetFuncTOWERName(int num)
{ static char *L[]={ "Full displacement",    //_FTOWER1
                     "Classic displacement", //_FTOWER2
                     "Full normal stress",   //_FTOWER3
                     "Classic normal stress",//_FTOWER4
                     "Full shear stress",    //_FTOWER5
                     "Classic shear stress"  //_FTOWER6
                   };
  return L[num];
}

void solver_Tower(void * context, CGrid * nd, double * F, double * par, int id_F)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm && nd && F) {
      int id_block, l;

////////////////////////////////////
//...обpабатываем все входные точки;
      if (nd->X && nd->Y && cont->sm->B)
      for (l = 0; l < nd->N; l++) {
           if (nd->hit) id_block = nd->hit[l];
           else         id_block = -1;
//////////////////////////////////
//...search for appropriate block;
           if (id_block < 0 && Sph3D_struc_in(cont->sm, id_block, nd->X[l], nd->Y[l], 0.) &&
               cont->sm->B[id_block].link && par  &&
               cont->sm->B[id_block].link[0] == 6 && cont->sm->bar && 
             ((CCeDemo *) cont->sm->bar)->in_test_inclusion2D(nd->X[l], nd->Y[l], par) != -cont->sm->B[id_block].link[5])
               id_block = cont->sm->B[id_block].link[6];
//////////////////////
//...data calculation;
           if (id_block >= 0 && id_block < cont->sm->N) {
               double FF[3] = {0., 0., 0.};
               switch (id_F) {
                      case _FK1: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 2);
                            F[l] = FF[0];
                            break;
                      case _FK2: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 0);
                            F[l] = FF[0];
                            break;
                      case _FK3: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 5);
                            F[l] = FF[0];
                            break;
                      case _FK4: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 3);
                            F[l] = FF[0];
                            break;
                      case _FK5: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 5);
                            F[l] = FF[1];
                            break;
                      case _FK6: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 3);
                            F[l] = FF[1];
                            break;
                      case _NUM_SK_FUNCTIONS: 
                            F[l] = id_block;
                      default:
                            F[l] = 0.;
               }
           }
           if (nd->hit) nd->hit[l] = id_block;
      }
  }
}
*/

/////////////////////////////////////////////////////////////////
//...интерфейсные функции для взаимодействия с внешними модулями;
/*double get_full_square(void * context)
{
   CWire2D * sm = (CWire2D *)GetSample(context);
   double S = sm ? sm->get_full_square() : MAX_HIT;
   SetSample(context, sm);
   return S;
}*/
#endif
#undef  Message
