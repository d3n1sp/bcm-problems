#include "stdafx.h"
#include "cpro_lpr.h"

#ifdef ___PRO_LPR_cpp___
///////////////////////////
//...готовые шаблоны таблиц;
Table * GetLineTable(int N_group)
{
	Shablon records[] = {{			 INT_TYPE_RECORD, NULL},//(0) -- SpanNumber
								{			CHAR_TYPE_RECORD, NULL},//(1) -- TowerCode 
								{			CHAR_TYPE_RECORD, NULL},//(2) -- TowerExtNumber
								{COMMENT_CHAR_TYPE_RECORD, NULL},//(3) -- TowerFoto
								{		 DANGLE_TYPE_RECORD, NULL},//(4) -- GPS_N 
								{		 DANGLE_TYPE_RECORD, NULL},//(5) -- GPS_E 
								{		DLENGTH_TYPE_RECORD, NULL},//(6) -- GPS_SpanLength 
								{		 DANGLE_TYPE_RECORD, NULL},//(7) -- RotationAngle 
								{		DLENGTH_TYPE_RECORD, NULL},//(8) -- SpanLength 
								{		DLENGTH_TYPE_RECORD, NULL},//(9) -- TowerJump 
								{		DLENGTH_TYPE_RECORD, NULL},//(10) -- EarthDiverg 
								{		DLENGTH_TYPE_RECORD, NULL},//(11) -- AverLength 
								{			CHAR_TYPE_RECORD, NULL},//(12) -- Mufta 
								{			 INT_TYPE_RECORD, NULL},//(13) -- N_circuit 
								{			 INT_TYPE_RECORD, NULL},//(14) -- N_grozo 
								{			CHAR_TYPE_RECORD, NULL},//(15) -- Shifr_Muft 
								{			CHAR_TYPE_RECORD, NULL},//(16) -- Izol_Marka 
								{			 INT_TYPE_RECORD, NULL},//(17) -- Izol_Type
								{			 INT_TYPE_RECORD, NULL},//(18) -- Izol_Count
								{			 INT_TYPE_RECORD, NULL},//(19) -- Kreplenie
								{			 INT_TYPE_RECORD, NULL},//(20) -- Krep_Shema
								{			 INT_TYPE_RECORD, NULL},//(21) -- ClimateIndex
								{			 INT_TYPE_RECORD, NULL},//(22) -- TerrainIndex 
								{COMMENT_CHAR_TYPE_RECORD, NULL},//(23) -- Comment
	};
	Table * table = get_shablon_table(24, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++)
		((int *)(table->table[NUMCIRCUIT].parm_values))[k] = 1;
	return table;
}

Table * GetTowerTable(int N_group)
{
	Shablon records[] = {{   CHAR_TYPE_RECORD, NULL},//(0) -- TowerCode
								{DLENGTH_TYPE_RECORD, NULL},//(1) -- геометрические размеры опор
								{	  ERR_TYPE_RECORD, NULL},//(2) -- признак продолжения...
	};
	Table * table = get_shablon_table(NUM_OF_TOWER_BASE_ELEMENTS+1, N_group, records, NULL);
	return  table;
}

Table * GetLineHangTable(int N_group)
{
	Shablon records[] = {{	  INT_TYPE_RECORD, NULL},//(0) -- SpanNumber
								{    INT_TYPE_RECORD, NULL},//(1) -- Cable_cross
								{DLENGTH_TYPE_RECORD, NULL},//(2) -- Cable_point
								{DLENGTH_TYPE_RECORD, NULL},//(3) -- Cable_shift
								{    INT_TYPE_RECORD, NULL},//(4) -- Grozo_cross 
								{DLENGTH_TYPE_RECORD, NULL},//(5) -- Insulator0 
								{DWEIGHT_TYPE_RECORD, NULL},//(6) -- Insulator0W 
								{    INT_TYPE_RECORD, NULL},//(7) -- Wire1_cross
								{DLENGTH_TYPE_RECORD, NULL},//(8) -- Insulator1 
								{DWEIGHT_TYPE_RECORD, NULL},//(9) -- Insulator1W
								{    INT_TYPE_RECORD, NULL},//(10) -- Wire2_cross
								{DLENGTH_TYPE_RECORD, NULL},//(11) -- Insulator2 
								{DWEIGHT_TYPE_RECORD, NULL},//(12) -- Insulator2W
								{    INT_TYPE_RECORD, NULL},//(13) -- Wire3_cross
								{DLENGTH_TYPE_RECORD, NULL},//(14) -- Insulator3 
								{DWEIGHT_TYPE_RECORD, NULL},//(15) -- Insulator3W
								{    INT_TYPE_RECORD, NULL},//(16) -- Wire4_cross
								{DLENGTH_TYPE_RECORD, NULL},//(17) -- Insulator4 
								{DWEIGHT_TYPE_RECORD, NULL},//(18) -- Insulator4W
								{    INT_TYPE_RECORD, NULL},//(19) -- Wire5_cross
								{DLENGTH_TYPE_RECORD, NULL},//(20) -- Insulator5 
								{DWEIGHT_TYPE_RECORD, NULL},//(21) -- Insulator5W
								{    INT_TYPE_RECORD, NULL},//(22) -- Wire6_cross
								{DLENGTH_TYPE_RECORD, NULL},//(23) -- Insulator6 
								{DWEIGHT_TYPE_RECORD, NULL},//(24) -- Insulator6W
	};
	Table * table = get_shablon_table(25, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((int *)(table->table[1].parm_values))[k] = ERR_STATE;
	}
	return table;
}

Table * GetLineSectTable(int N_group)
{
	Shablon records[] = {{			 INT_TYPE_RECORD, NULL},//(0) -- SpanNumber
								{		DLENGTH_TYPE_RECORD, NULL},//(1) -- L_to_Obj (<0 -- измерение от правой опоры, т.е. надо прибавить длину пролета)
								{		DLENGTH_TYPE_RECORD, NULL},//(2) -- H_Obj
								{			CHAR_TYPE_RECORD, NULL},//(3) -- Name
								{			CHAR_TYPE_RECORD, NULL},//(4) -- Data
								{		DTEMPER_TYPE_RECORD, NULL},//(5) -- Temperatura
								{		DLENGTH_TYPE_RECORD, NULL},//(6) -- H_L_Podves
								{		DLENGTH_TYPE_RECORD, NULL},//(7) -- H_R_Podves
								{		DLENGTH_TYPE_RECORD, NULL},//(8) -- Proves
								{		DLENGTH_TYPE_RECORD, NULL},//(9) -- H_Provod
								{		DLENGTH_TYPE_RECORD, NULL},//(10) -- W_Obj
								{			 INT_TYPE_RECORD, NULL},//(11) -- Type
								{		 DANGLE_TYPE_RECORD, NULL},//(12) -- UU(PlanRotationAngle)
								{		DLENGTH_TYPE_RECORD, NULL},//(13) -- DH1(Pop_h)
								{		DLENGTH_TYPE_RECORD, NULL},//(14) -- DH2(Pop_dh)
								{		DLENGTH_TYPE_RECORD, NULL},//(15) -- L1(PlanSectionDXTowerA)
								{		DLENGTH_TYPE_RECORD, NULL},//(16) -- L2(PlanInverseSectionDXTowerB)
								{		DLENGTH_TYPE_RECORD, NULL},//(17) -- H1(Pop_Ha)
								{		DLENGTH_TYPE_RECORD, NULL},//(18) -- H2(Pop_Hb)
								{COMMENT_CHAR_TYPE_RECORD, NULL},//(19) -- Comment
	};
	Table * table = get_shablon_table(20, N_group, records, NULL);
	return table;
}

Table * GetLineProfileTable(int N_group, int N_points)
{
	Shablon records[] = {{	  INT_TYPE_RECORD, NULL},//(0) -- SpanNumber
								{	  INT_TYPE_RECORD, NULL},//(1) -- PointsNumber
								{DLENGTH_TYPE_RECORD, NULL},//(2) -- L_point/H_point
								{	  ERR_TYPE_RECORD, NULL},//(3) -- признак продолжения...
	};
	Table * table = get_shablon_table(N_points+2, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++)
		((int *)(table->table[1].parm_values))[k] = N_points;
	return table;
}

Table * GetLineClimateTable(int N_group)
{
	Shablon records[] = {{          INT_TYPE_RECORD, NULL},//(0) -- Num_rush_zone
								{          INT_TYPE_RECORD, NULL},//(1) -- Num_glaze_zone
								{      DLENGTH_TYPE_RECORD, NULL},//(2) -- C_glaze
								{      DLENGTH_TYPE_RECORD, NULL},//(3) -- C_rush_glaze
								{      DTEMPER_TYPE_RECORD, NULL},//(4) -- T_glaze
								{		 DTEMPER_TYPE_RECORD, NULL},//(5) -- T_rush
								{      DTEMPER_TYPE_RECORD, NULL},//(6) -- T_aver
								{      DTEMPER_TYPE_RECORD, NULL},//(7) -- T_abs_max
								{      DTEMPER_TYPE_RECORD, NULL},//(8) -- T_abs_min
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(9) -- Q_rush
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(10) -- Q_rush_glaze
								{    DVELOCITY_TYPE_RECORD, NULL},//(11) -- V_max_wind
								{    DVELOCITY_TYPE_RECORD, NULL},//(12) -- V_glaze
								{DMASS_DENSITY_TYPE_RECORD, NULL},//(13) -- Ro_glaze
								{          INT_TYPE_RECORD, NULL},//(14) -- Condition_index
	};
	Table * table = get_shablon_table(15, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((int    *)(table->table[0].parm_values)) [k] = ERR_STATE;
		((int    *)(table->table[1].parm_values)) [k] = ERR_STATE;
		((double *)(table->table[13].parm_values))[k] = 0.9;
	}
	return table;
}

Table * GetLineModeTable(int N_group)
{
	Shablon records[] = {{			  INT_TYPE_RECORD, NULL},//(0) -- SpanNumber
								{			 CHAR_TYPE_RECORD, NULL},//(1) -- WireMark
								{			  INT_TYPE_RECORD, NULL},//(2) -- CurrentWire
								{			  INT_TYPE_RECORD, NULL},//(3) -- CurrentWireProfile
								{		 DLENGTH_TYPE_RECORD, NULL},//(4) -- Diam
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(5) -- Mode_p1
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(6) -- Mode_p3
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(7) -- Mode_p6
								{  DFORCE_LINE_TYPE_RECORD, NULL},//(8) -- Mode_p7
								{		  DFORCE_TYPE_RECORD, NULL},//(9) -- T_break
								{		  DFORCE_TYPE_RECORD, NULL},//(10) -- T_max
								{		  DFORCE_TYPE_RECORD, NULL},//(11) -- T_opt								
								{		DPERCENT_TYPE_RECORD, NULL},//(12) -- T_max_tension(%)
								{		DPERCENT_TYPE_RECORD, NULL},//(13) -- T_max_P7(%)	 -- не используется;
								{		DPERCENT_TYPE_RECORD, NULL},//(14) -- T_max_tmin(%) -- не используется;

//								{		  DFORCE_TYPE_RECORD, NULL},//(10) -- T_abs_max
//								{		  DOUBLE_TYPE_RECORD, NULL},//(11) -- T_max_P7(%)
//								{		  DOUBLE_TYPE_RECORD, NULL},//(12) -- T_max_tmin(%)
//								{		  DOUBLE_TYPE_RECORD, NULL},//(13) -- T_opt_P1(%)
//								{		  DOUBLE_TYPE_RECORD, NULL},//(14) -- T_max_tension(%)

								{DMASS_DENSITY_TYPE_RECORD, NULL},//(15) -- EF_creep
								{DMASS_DENSITY_TYPE_RECORD, NULL},//(16) -- EF
								{DMASS_DENSITY_TYPE_RECORD, NULL},//(17) -- EF_fin
								{DMASS_DENSITY_TYPE_RECORD, NULL},//(18) -- Alpha
								{      DSQUARE_TYPE_RECORD, NULL},//(19) -- EffectSquare
								{		  DOUBLE_TYPE_RECORD, NULL},//(20) -- ResidualStrain
								{			  INT_TYPE_RECORD, NULL},//(21) -- Type
								{			  INT_TYPE_RECORD, NULL},//(22) -- NP_AMB
								{		  DOUBLE_TYPE_RECORD, NULL},//(23) -- DP_AMB
								{		  DOUBLE_TYPE_RECORD, NULL},//(24) -- RS_AMB
								{		  DOUBLE_TYPE_RECORD, NULL},//(25) -- RD_AMB
								{			  INT_TYPE_RECORD, NULL},//(26) -- NP_Stal
								{		  DOUBLE_TYPE_RECORD, NULL},//(27) -- DP_Stal
								{		  DOUBLE_TYPE_RECORD, NULL},//(28) -- RS_Stal
								{		  DOUBLE_TYPE_RECORD, NULL},//(29) -- RD_Stal
								{ COMMENT_CHAR_TYPE_RECORD, NULL},//(30) -- Comment
	};
	Table * table = get_shablon_table(31, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((int *)(table->table[2].parm_values))[k] = ERR_STATE;
		((int *)(table->table[3].parm_values))[k] = ERR_STATE;
	}
	return table;
}

Table * GetLineTensionTable(int N_group)
{
	Shablon records[] = {{DLENGTH_TYPE_RECORD, NULL},//(0) -- Cable_D 
								{ DFORCE_TYPE_RECORD, NULL},//(1) -- Cable_T 
								{DLENGTH_TYPE_RECORD, NULL},//(2) -- Cable_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(3) -- Cable_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(4) -- Cable_I_W

								{DLENGTH_TYPE_RECORD, NULL},//(5) -- Grozo_D 
								{ DFORCE_TYPE_RECORD, NULL},//(6) -- Grozo_T 
								{DLENGTH_TYPE_RECORD, NULL},//(7) -- Grozo_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(8) -- Grozo_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(9) -- Grozo_I_W

								{DLENGTH_TYPE_RECORD, NULL},//(10) -- Wire1_D 
								{ DFORCE_TYPE_RECORD, NULL},//(11) -- Wire1_T 
								{DLENGTH_TYPE_RECORD, NULL},//(12) -- Wire1_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(13) -- Wire1_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(14) -- Wire1_I_W

								{DLENGTH_TYPE_RECORD, NULL},//(15) -- Wire2_D 
								{ DFORCE_TYPE_RECORD, NULL},//(16) -- Wire2_T
								{DLENGTH_TYPE_RECORD, NULL},//(17) -- Wire2_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(18) -- Wire2_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(19) -- Wire2_I_W
								
								{DLENGTH_TYPE_RECORD, NULL},//(20) -- Wire3_D 
								{ DFORCE_TYPE_RECORD, NULL},//(21) -- Wire3_T
								{DLENGTH_TYPE_RECORD, NULL},//(22) -- Wire3_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(23) -- Wire3_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(24) -- Wire3_I_W
								
								{DLENGTH_TYPE_RECORD, NULL},//(25) -- Wire4_D 
								{ DFORCE_TYPE_RECORD, NULL},//(26) -- Wire4_T
								{DLENGTH_TYPE_RECORD, NULL},//(27) -- Wire4_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(28) -- Wire4_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(29) -- Wire4_I_W
								
								{DLENGTH_TYPE_RECORD, NULL},//(30) -- Wire5_D 
								{ DFORCE_TYPE_RECORD, NULL},//(31) -- Wire5_T
								{DLENGTH_TYPE_RECORD, NULL},//(32) -- Wire5_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(33) -- Wire5_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(34) -- Wire5_I_W
								
								{DLENGTH_TYPE_RECORD, NULL},//(35) -- Wire6_D 
								{ DFORCE_TYPE_RECORD, NULL},//(36) -- Wire6_T 
								{DLENGTH_TYPE_RECORD, NULL},//(37) -- Wire6_I_dL 
								{DLENGTH_TYPE_RECORD, NULL},//(38) -- Wire6_I_L 
								{DWEIGHT_TYPE_RECORD, NULL},//(39) -- Wire6_I_W
	};
	Table * table = get_shablon_table(40, N_group, records, NULL);
	return table;
}

Table * GetLineDampherTable(int N_group)
{
	Shablon records[] = {{   CHAR_TYPE_RECORD, NULL},//(0) -- SchemeCode
								{	  INT_TYPE_RECORD, NULL},//(1) -- ClosingType
								{   CHAR_TYPE_RECORD, NULL},//(2) -- FittingMark1
								{   CHAR_TYPE_RECORD, NULL},//(3) -- FittingMark2
								{ DOUBLE_TYPE_RECORD, NULL},//(4) -- MaxStressPercent
								{DLENGTH_TYPE_RECORD, NULL},//(5) -- MaxSpanLength

								{   CHAR_TYPE_RECORD, NULL},//(6) -- DampherMark1
								{DLENGTH_TYPE_RECORD, NULL},//(7) -- Dampher1

								{   CHAR_TYPE_RECORD, NULL},//(8) -- DampherMark2
								{DLENGTH_TYPE_RECORD, NULL},//(9) -- Dampher2

								{   CHAR_TYPE_RECORD, NULL},//(10) -- DampherMark3
								{DLENGTH_TYPE_RECORD, NULL},//(11) -- Dampher3

								{   CHAR_TYPE_RECORD, NULL},//(12) -- DampherMark4
								{DLENGTH_TYPE_RECORD, NULL},//(13) -- Dampher4
	};
	Table * table = get_shablon_table(14, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((int *)(table->table[NUMMAXSPAN-4].parm_values))[k] = PZS_PZS_STATE;
	}
	return table;
}

Table * GetLineDampherTable_Aeol(int N_group)
{
	Shablon records[] = {{	  INT_TYPE_RECORD, NULL},//(0) -- Fitting
								{	  INT_TYPE_RECORD, NULL},//(1) -- Category
								{ DOUBLE_TYPE_RECORD, NULL},//(2) -- StressPercent
								{	  INT_TYPE_RECORD, NULL},//(3) -- ClosingType
								{	  INT_TYPE_RECORD, NULL},//(4) -- DampherScheme
								{   CHAR_TYPE_RECORD, NULL},//(5) -- DampherMark1
								{	 CHAR_TYPE_RECORD, NULL},//(6) -- DampherMark2
								{DLENGTH_TYPE_RECORD, NULL},//(7) -- LeftDampher1
								{DLENGTH_TYPE_RECORD, NULL},//(8) -- LeftDampher2
								{DLENGTH_TYPE_RECORD, NULL},//(9) -- RightDampher1
								{DLENGTH_TYPE_RECORD, NULL},//(10) -- RightDampher2
								{	  INT_TYPE_RECORD, NULL},//(11) -- LoopDampher
								{DLENGTH_TYPE_RECORD, NULL},//(12) -- LoopLength
								{DLENGTH_TYPE_RECORD, NULL},//(13) -- LoopSag
								{DLENGTH_TYPE_RECORD, NULL},//(14) -- MinLength
								{DLENGTH_TYPE_RECORD, NULL},//(15) -- MaxLength
	};
	Table * table = get_shablon_table(16, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((int *)(table->table[0].parm_values))[k] = ERR_STATE;
		((int *)(table->table[1].parm_values))[k] = 1;
		((int *)(table->table[3].parm_values))[k] = ERR_STATE;
		((int *)(table->table[4].parm_values))[k] = ERR_STATE;
		((int *)(table->table[11].parm_values))[k] = ERR_STATE;
	}
	return table;
}

///////////////////////////////////////////////////
//...шаблон таблиц для задания предельных нагрузок;
Table * GetTenLimitTable(int N_group, char * table_names[])
{
	Shablon records[] = {{	  CHAR_TYPE_RECORD, "Марка"},						//(0) -- Марка
								{  DFORCE_TYPE_RECORD, "T_break"},					//(1) -- T_break
								{  DFORCE_TYPE_RECORD, "Предельное тяжение"},   //(2) -- T_max
								{  DFORCE_TYPE_RECORD, "Оптимальное тяжение"},  //(3) -- T_opt
								{DPERCENT_TYPE_RECORD, "Предельное напряжение"},//(4) -- T_max_abs/T_break*100
	};
	Table * table = get_shablon_table(5, N_group, records, table_names, INITIAL_STATE);
	return table;
}

//////////////////////////////
//...шаблон монтажной таблицы;
Table * GetMontageTable(int N_group, int N, double * temper)
{
	Shablon records[] = {{DFORCE_TYPE_RECORD, NULL},//(0) -- Cable_T for temper[0] и т.д.
								{	 ERR_TYPE_RECORD, NULL},//(1) -- признак продолжения...
	};
	Table * table = get_shablon_table(N, N_group, records, NULL);

///////////////////////////////////////////////////////////
//...вписываем температуру, как динамическое имя параметра;
	for (int i = 0; table && temper && i < table->N; i++) {
		char buff[201], * sss;
		sprintf(buff, "%g", temper[i]);
		sss = (char *)new_struct(((int)strlen(buff)+1)*sizeof(char));
		::strcpy(sss, buff);
		table->table[i].parm_name = sss;
	}
	return table;
}

////////////////////////////////////////////
//...шаблон таблиц для базы данных проводов;
Table * GetWireDataBaseTable(int N_group, char * table_names[], int id_static_char)
{
	Shablon records[] = {{					 CHAR_TYPE_RECORD, "Число проволок"},	   //(0) -- N_wire															
								{		  DLENGTH_CHAR_TYPE_RECORD, "Диаметр проволок"},	//(1) -- diameter
								{					 CHAR_TYPE_RECORD, "Число повивов"},		//(2) -- N_layer
								{					 CHAR_TYPE_RECORD, "Кратн. шага скрутки"},//(3) -- H_twist
								{				 DWEIGHT_TYPE_RECORD, "Масса провода"},		//(4) -- M_wire
								{				 DWEIGHT_TYPE_RECORD, "Масса смазки"},		   //(5) -- M_lube
								{				  DFORCE_TYPE_RECORD, "Разрывное усилие"},   //(6) -- F_break
								{		  DSTRESS_CHAR_TYPE_RECORD, "Модуль упругости"},   //(7) -- E_module
								{					 CHAR_TYPE_RECORD, "Коэф. Пуассона"},	   //(8) -- nju
								{DTEMPER_EXTENT_CHAR_TYPE_RECORD, "Коэф. темп. расш."},  //(9) -- A_temp
	};
	Table * table = get_shablon_table(10, N_group, records, table_names, id_static_char);
	return  table;
}

///////////////////////////////////////////
//...шаблон таблиц для базы данных кабелей;
Table * GetCableDataBaseTable(int N_group, char * table_names[], int id_static_char)
{
	Shablon records[] = {{		  DLENGTH_TYPE_RECORD, "Diam"},	  //(0) -- Diam
								{		  DSQUARE_TYPE_RECORD, "F_effect"},//(1) -- F_eff
								{		  DSTRESS_TYPE_RECORD, "E_creep"}, //(2) -- E_creep
								{		  DSTRESS_TYPE_RECORD, "E_ini"},	  //(3) -- E_ini
								{		  DSTRESS_TYPE_RECORD, "E_fin"},	  //(4) -- E_fin
								{		  DWEIGHT_TYPE_RECORD, "M_cable"}, //(5) -- M_cable
								{			DFORCE_TYPE_RECORD, "F_break"}, //(6) -- F_break
								{DTEMPER_EXTENT_TYPE_RECORD, "A_temp"},  //(7) -- A_temp
	};
	Table * table = get_shablon_table(8, N_group, records, table_names, id_static_char);
	return  table;
}

////////////////////////////////////////
//...шаблон таблиц для базы данных опор;
Table * GetTowerDataBaseTable(int N_group, char * table_names[], int id_static_char)
{
	Shablon records[] = {{DLENGTH_TYPE_RECORD, "h1"},			 //(0) -- h1
								{DLENGTH_TYPE_RECORD, "h2"},			 //(1) -- h2
								{DLENGTH_TYPE_RECORD, "h3"},			 //(2) -- h3
								{DLENGTH_TYPE_RECORD, "h4"},			 //(3) -- h4
								{DLENGTH_TYPE_RECORD, "a1"},			 //(4) -- a1
								{DLENGTH_TYPE_RECORD, "a2"},			 //(5) -- a2
								{DLENGTH_TYPE_RECORD, "a3"},			 //(6) -- a3
								{DLENGTH_TYPE_RECORD, "b1"},			 //(7) -- b1
								{DLENGTH_TYPE_RECORD, "b2"},			 //(8) -- b2
								{DLENGTH_TYPE_RECORD, "b3"},			 //(9) -- b3
								{DLENGTH_TYPE_RECORD, "c1"},			 //(10) -- c1
								{DLENGTH_TYPE_RECORD, "c2"},			 //(11) -- c2
								{    INT_TYPE_RECORD, "Схема опоры"},//(12) -- N_sheme
	};
	Table * table = get_shablon_table(13, N_group, records, table_names, id_static_char);
	return  table;
}

////////////////////////////////////////////////////////////////
//...шаблон таблиц для базы данных схем защиты провода и кабеля;
Table * GetDampherDataBaseTable(int N_group)
{
	Shablon records[] = {{   CHAR_TYPE_RECORD, NULL},//(0) -- FittingMark1
								{   CHAR_TYPE_RECORD, NULL},//(1) -- FittingMark2
								{ DOUBLE_TYPE_RECORD, NULL},//(2) -- StressPercent
								{DLENGTH_TYPE_RECORD, NULL},//(3) -- SpanLength
								{   CHAR_TYPE_RECORD, NULL},//(4) -- DampherMark1
								{   CHAR_TYPE_RECORD, NULL},//(5) -- DampherMark2
								{DLENGTH_TYPE_RECORD, NULL},//(6) -- Dampher1
								{DLENGTH_TYPE_RECORD, NULL},//(7) -- Dampher2
								{DLENGTH_TYPE_RECORD, NULL},//(8) -- Dampher3
								{DLENGTH_TYPE_RECORD, NULL},//(9) -- Dampher4
	};
	Table * table = get_shablon_table(10, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++) {
		((double *)(table->table[2].parm_values))[k] = MAX_HIT;
		((double *)(table->table[3].parm_values))[k] = MAX_HIT;
		((double *)(table->table[6].parm_values))[k] = MAX_HIT;
		((double *)(table->table[7].parm_values))[k] = MAX_HIT;
		((double *)(table->table[8].parm_values))[k] = MAX_HIT;
		((double *)(table->table[9].parm_values))[k] = MAX_HIT;
	}
	return table;
}

/////////////////////////////////////////////////////////////
//...шаблон таблиц для ввода схем защиты кабеля и грозотроса;
Table * GetDampherSchemeTable(int N_group, int id_name)
{
	Shablon records[] = {{   CHAR_TYPE_RECORD, "Номер схемы"},//(0) -- Scheme number
								{DLENGTH_TYPE_RECORD, "Макс. длина"},//(1) -- MaxSpanLength
								{   CHAR_TYPE_RECORD, "Левый  гаситель"},//(2) -- DampherMark1
								{   CHAR_TYPE_RECORD, "Правый гаситель"},//(3) -- DampherMark2
								{   CHAR_TYPE_RECORD, "Доп. гаситель 1"},//(4) -- DampherMark3
								{   CHAR_TYPE_RECORD, "Доп. гаситель 2"},//(5) -- DampherMark4
	};
	if (! id_name) 
		for (int k = 0; k < 6; k++)  records[k].parm_name = NULL;

	Table * table = get_shablon_table(6, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++)
		((double *)(table->table[1].parm_values))[k] = MAX_HIT;

	return table;
}

////////////////////////////////////////////
//...извлекаем эксплуатационный номер опоры;
char * operating_num(Table * table, int i, char * buff, int num_project)
{
	char * num, * sss, * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	setlocale(LC_COLLATE, "rus"); //...сравниваем без учета регистра;
	if (num_project == OK_STATE) {
		sprintf(buff, "%i", ((int **)table->table[0].parm_values)[i+1]);  //...проектный номер;
		num   = buff;
	}
	else
	if (num_project == ADDITIONAL_STATE) {
		sprintf(buff, "%i/%s", ((int  **)table->table[0].parm_values)[i+1],
									  ((char **)table->table[2].parm_values)[i+1]);  //...проектный/эксплуатац. номер;
		num   = buff;
	}
	else
	if ((sss = ((char **)table->table[2].parm_values)[i+1]) == NULL ||
		strlen(sss) == NULL) {
		sprintf(buff, "%i", ((int **)table->table[0].parm_values)[i+1]);  //...проектный номер;
		num   = buff;
	}
	else if (_stricmp("БН", sss) == NULL) { //...в скобочках проектный номер;
		sprintf(buff, "%s(%i)", sss, ((int **)table->table[0].parm_values)[i+1]);
		num   = buff;
	}
	else num = sss; //...эксплуатационный номер;
	setlocale(LC_COLLATE, sloc); free(sloc);
	return num;
}

///////////////////////////////
//...определяем анкерную опору;
int id_anker_num(Table * table, int i, int id_tower)
{
	if (! table) return(1);
	char * sss, * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	int m; i = min(max(0, i), table->N_group-1); 
	setlocale(LC_COLLATE, "rus"); //...сравниваем без учета регистра;
	if (id_tower) m = (sss = ((char **)table->table[1].parm_values)[i+1]) == NULL || ! _strnicmp(sss, "У", 1) || ! _strnicmp(sss, "АУ", 1) || ! _strnicmp(sss, "опора", 5);
	else			  m = (sss = ((char **)table->table[NUMCIRCUIT-1].parm_values)[i+1]) != NULL && ! _strnicmp(sss, "анкер", 5);
	setlocale(LC_COLLATE, sloc); free(sloc);
	return m;
}

/////////////////////////////////////////////////////////
// ...корректируем позиции пустого комментария в таблице;
void comment_correct(Table * table, int id_static_char) 
{
	if (id_static_char < ADDITIONAL_STATE)
	for (int i = 0; table && i < table->N; i++) 
		if (table->table[i].type == COMMENT_CHAR_TYPE_RECORD)
			for (int j = 1; j <= table->N_group; j++) 
				if (! ((char **)table->table[i].parm_values)[j])
						((char **)table->table[i].parm_values)[j] = new_struct(sizeof(char));
//					SetTableParam(table, j, i, (unsigned long)"");
}
void comment_correct(void * context) { 
	comment_correct((Table *)GetTable(context), GetIdStaticChar(context));
}

//////////////////////////////////////////////////////////////////////////////////
// ...корректируем таблицы подвески, сдвигая пролет на анкерный или промежуточный;
void hunging_correct(Table * table, Table * tab_ank, int * anker_num, int id_anker) 
{
	for (int j = 0; table && j < table->N_group; j++) {
		int i = abs((int)GetTableParam(table, j+1, 0))-1;

		if (anker_num && anker_num[0]) {
			if (i < anker_num[1]) i = anker_num[1];
			while (i < anker_num[anker_num[0]+1] && 
					(id_anker && ! id_anker_num(tab_ank, i) || ! id_anker && id_anker_num(tab_ank, i))) i++;

			SetTableParam(table, j+1, 0, i+1.);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// ...выбрасываем повторяющиеся строки (отдельно в климатической таблице);
void table_double_correct(Table * table, int excl_pos, int id_static_char, int id_sequent) 
{
	for (int i = table ? table->N_group : 0; i > 1; i--)
	for (int j = i-1; j > 0; j--) {
		int  first = 0;
		for (int k = 0; ! first && k < table->N; k++)
		if ( excl_pos <= ERR_STATE || k != excl_pos) {
			if (table->table[k].type < CHAR_TYPE_RECORD); else
			if (table->table[k].type <  INT_TYPE_RECORD) {
				char * str1 = (char *)(long long)GetTableParam(table, i, k),
					  * str2 = (char *)(long long)GetTableParam(table, j, k);
				if (str1 && ! str2 || ! str1 && str2 || str1 && str2 && ::strcmp(str1, str2)) first = 1;
			}
			else
			if (table->table[k].type < NUM_OF_TYPE_RECORD) {
				if (fabs(GetTableParam(table, i, k)-GetTableParam(table, j, k)) >= EE_ker) first = 1;
			}
			if (first && ERR_STATE-excl_pos-1 == k && i < table->N_group) { //... особый случай проверки исключаемой позиции;
				if (table->table[k].type < CHAR_TYPE_RECORD); else
				if (table->table[k].type <  INT_TYPE_RECORD) {
					char * str1 = (char *)(long long)GetTableParam(table, i, k),
						  * str2 = (char *)(long long)GetTableParam(table, i+1, k);
					if (str1 && ! str2 || ! str1 && str2 || str1 && str2 && ::strcmp(str1, str2)) first = 0;
				}
				else
				if (table->table[k].type < NUM_OF_TYPE_RECORD) {
					if (fabs(GetTableParam(table, i, k)-GetTableParam(table, i+1, k)) >= EE_ker) first = 0;
				}
			}
		}	
		if (! first) {
			del_table_group(table, 1, i, id_static_char); 
			j = 0;
		}
		if (id_sequent) break;
	}
}
void table_double_correct(Table * table, Table * tab_ank, int id_static_char) 
{
	for (int i = table ? table->N_group : 0; i > 1; i--)
	for (int j = i-1; j > 0; j--) {
		int first = 0;
		for (int k = 0; ! first && k < table->N; k++)
		if (fabs(GetTableParam(table, i, k)-GetTableParam(table, j, k)) >= EE_ker) first = 1;
		if (! first) {
			del_table_group(table, 1, i, id_static_char);
			for (int l = tab_ank ? tab_ank->N_group : 0; l > 1; l--)
			if ((int)GetTableParam(tab_ank, l, NUMCIRCUIT+8) == i) SetTableParam(tab_ank, l, NUMCIRCUIT+8, (double)j);
			j = 0;
		}
	}
}

//////////////////////////////////////////////////////////////////
//...переустанавливаем линию для анализа схем виброзащиты провода;
void ResetLineForDampher(void ** anker_list) 
{
//////////////////////////////
//...устанавливаем свои опоры;
	int N_group = GetNGroup(anker_list[NUMANKERTAB]), k;
	for (k = 0;  k < N_group; k++)
		if (id_anker_num((Table *)GetTable(anker_list[NUMANKERTAB]), k)) 
			 SetTableParam(anker_list[NUMANKERTAB], k+1, 1, (double)(long long)"У"); else
			 SetTableParam(anker_list[NUMANKERTAB], k+1, 1, (double)(long long)"П");

//////////////////////////////////////////////
//...устанавливаем таблицу режимов нагружений;
	SetTableParam(anker_list[NUMWIRE1MODETAB], 0, GetTableParam(anker_list[NUMANKERTAB], 1, 0));

///////////////////////////////////
//...стираем все остальные таблицы;
	for (k = NUMANKERTAB+1;  k < NUMLINEPROJECT; k++)
		if (k != NUMCLIMATETAB && k != NUMWIRE1MODETAB && k != NUMWIRE1MONTAGETAB && k != NUMWIRE1DAMPHERTAB)
			SetTable(anker_list[k], NULL);
}

// ======================================================================================================
// Line reading  ========================================================================================
// ======================================================================================================
char * LineReading(char * pchar, void * context, int head_reading)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);
	int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);

////////////////////////////////////
//...зачитываем имя и описание линии;
	if (((p_end = strstr(pchar, TABLESTR)) != NULL|| (p_end = p_eof) != NULL) && head_reading) {
		if ((pnext = strstr(pchar, STRSEPARATOR)) != NULL || (pnext = p_end) != NULL) {
			char Temp0 = pnext[0];
			pnext[0]   = '\x0';

			SetSampleComment(context, pchar);

			pnext[0] = Temp0;
			if ((pchar = pnext) < p_end) pchar += lenStrSeparator;
		}
		if ( pchar < (pnext = p_end)) pnext -= lenStrSeparator;
		char Temp2 =  pnext[0];
		pnext[0]   = '\x0';

		SetSampleDescription(context, pchar);

		pnext[0] = Temp2;
	}
	pchar = p_end;

////////////////////////////////
//...зачитываем параметры линии;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

///////////////////////////////////
//...зачитываем число опор в линии;
		if (((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

///////////////////////////////////
//...порождаем новую таблицу линии;
			Table * table = GetLineTable(N_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Non-standart towers reading  =========================================================================
// ======================================================================================================
char * LineTowerReading(char * pchar, void * context)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

/////////////////////////////////////////////
//...зачитываем параметры нестандартных опор;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

///////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице нестандартных опор;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

//////////////////////////////////////////////////////////
//...порождаем новую таблицу нестандартных опор для линии;
			Table * table = GetTowerTable(N_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Line hanging reading  ===============================================================================
// ======================================================================================================
char * LineHangReading(char * pchar, void * context, int convert)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

///////////////////////////////////////////
//...зачитываем параметры условий подвески;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

//////////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице условий подвески на линии;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

///////////////////////////////////////////////////////
//...порождаем новую таблицу условий подвески на линии;
			Table * table = GetLineHangTable(N_group);
			if (table) {
///////////////////////////////////////
//...заполняем таблицу данными с диска;
					if (convert) {
					Shablon records = {CHAR_TYPE_RECORD, NULL};
					add_table_param(table, 1, & records, 1);
				}
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
				if (convert) del_table_param(table, 1, 2, NULL_STATE);

				SetTable (context, table);
				SetTableIndex(context, 1);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Line intersections reading  =========================================================================
// ======================================================================================================
char * LineSectReading(char * pchar, void * context, int convert)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

///////////////////////////////////////////////
//...зачитываем параметры пересечений на линии;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

/////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице пересечений на линии;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

//////////////////////////////////////////////////
//...порождаем новую таблицу пересечений на линии;
			Table * table = GetLineSectTable(N_group);
			if (table) {
///////////////////////////////////////
//...заполняем таблицу данными с диска;
				if (convert) {
					Shablon records = {CHAR_TYPE_RECORD, NULL};
					add_table_param(table, 1, & records, 1);
				}
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
				if (convert) del_table_param(table, 1, 2, NULL_STATE);
				SetTable (context, table);
				SetTableIndex(context, 1);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Line profile reading  ================================================================================
// ======================================================================================================
char * LineProfileReading(char * pchar, void * context, int convert)
{
	char * p_end, * pnext, * pname, 
		  * p_eof = pchar+strlen(pchar);

////////////////////////////////////////
//...зачитываем параметры профиля линии;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

//////////////////////////////////////////////////////
//...зачитываем число строчек в таблице профиля линии;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

///////////////////////////////////////////
//...порождаем новую таблицу профиля линии;
			Table * table = GetLineProfileTable(N_group);
			if (table) {
///////////////////////////////////////
//...заполняем таблицу данными с диска;
				if (convert) {
					Shablon records = {CHAR_TYPE_RECORD, NULL};
					add_table_param(table, 1, & records, 1);
				}
				int i = 0, N_points = 0, m = convert ? 1 : 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
/////////////////////////////////
//...выравниваем строчки таблицы;
						Shablon records[] = {{DLENGTH_TYPE_RECORD, NULL},//(0) -- L_point/H_point
													{	  ERR_TYPE_RECORD, NULL},//(1) -- признак продолжения...
						};
						if (pchar) pname = strstr(pchar,   SEMICOLONSTR);
						if (pname && convert) pname = strstr(pname+1, SEMICOLONSTR);
						if (! pname || sscanf(pname+strlen(COLONSTR)+1, "%i", &N_points) != 1) N_points = 0;
						if (N_points >= 0 && N_points > table->N-2-m) 
							add_table_param(table, N_points-table->N+2+m, records);
						if (N_points < 0 && 2*abs(N_points) > table->N-2-m) 
							add_table_param(table, 2*abs(N_points)-table->N+2+m, records);
////////////////////////////////
//...зачитываем строчку таблицы;
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
				if (convert) del_table_param(table, 1, 2, NULL_STATE);
				SetTable (context, table);
				SetTableIndex(context, 1);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Climate parameters reading  ==========================================================================
// ======================================================================================================
char * LineClimateReading(char * pchar, void * context)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

/////////////////////////////////////////////////////////
//...зачитываем параметры климатических условий на линии;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

/////////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице климатических параметров;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

//////////////////////////////////////////////////////
//...порождаем новую таблицу климатических параметров;
			Table * table = GetLineClimateTable(N_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading mode parameters ==============================================================================
// ======================================================================================================
char * LineModeReading(char * pchar, void * context, int convert)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

/////////////////////////////////////////////
//...зачитываем параметры условий нагружения;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

///////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице нагрузочных параметров;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);

		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

////////////////////////////////////////////////////
//...порождаем новую таблицу нагрузочных параметров;
			Table * table = GetLineModeTable(N_group);
			if (table) {
///////////////////////////////////////
//...заполняем таблицу данными с диска;
				if (convert) {
					Shablon records = {CHAR_TYPE_RECORD, NULL};
					add_table_param(table, 1, & records, 0);
				}
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
				if (convert) {
					Record temp = table->table[0];
					table->table[0] = table->table[2];
					table->table[2] = temp;
					del_table_param(table, 1, 1, NULL_STATE);
				}
				SetTable (context, table);
				SetTableIndex(context, 1);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading montage table  ===============================================================================
// ======================================================================================================
char * LineMontageReading(char * pchar, void * context)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

////////////////////////////////////
//...зачитываем параметры монтажной;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

//////////////////////////////////////////////////
//...зачитываем число строчек в монтажной таблице;
		int  N_group = 0, N = 0, lenStrSeparator = (int)strlen(STRSEPARATOR), lenSemiColonStr = (int)strlen(SEMICOLONSTR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) && 
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

////////////////////////////////
//...вычисляем число параметров;
			char *  p_str = strstr(pchar, STRSEPARATOR);
			while ((pchar = strstr(pchar, SEMICOLONSTR)) != NULL && (pchar += lenSemiColonStr) <= p_str) N++;
			pchar = pnext;

///////////////////////////////////////
//...порождаем новую монтажную таблицу;
			Table * table = GetMontageTable(N_group, N);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext, " ");
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading scheme parameters ============================================================================
// ======================================================================================================
char * LineSchemeReading(char * pchar, void * context)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

/////////////////////////////////////////////
//...зачитываем параметры условий нагружения;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

///////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице нагрузочных параметров;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);

		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

////////////////////////////////////////////////////
//...порождаем новую таблицу нагрузочных параметров;
			Table * table = GetDampherSchemeTable(N_group, NULL_STATE);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

char * LineDampherReading(char * pchar, void * context, void * ankertab)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

///////////////////////////////////////////////////////////////
//...зачитываем параметры демпфирующих схем защиты от вибрации;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

////////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице демпфирующих параметров;
		Table * tab_ank = (Table *)GetTable(ankertab);
		int  N_anker_group = tab_ank ? tab_ank->N_group : 0, lenStrSeparator = (int)strlen(STRSEPARATOR), 
			  N_group = 0;

		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

/////////////////////////////////////////
//...порождаем новую таблицу демпфирующих параметров;
			Table * table = GetLineDampherTable(N_group ? N_group : N_anker_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading tension parameters (reading or setting table) ================================================
// ======================================================================================================
char * LineTensionReading(char * pchar, void * context, void * ankertab)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

///////////////////////////////////////////
//...зачитываем параметры тяжений на линии;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

////////////////////////////////////////////////
//...зачитываем число строчек в таблице тяжений;
		Table * tab_ank = (Table *)GetTable(ankertab);
		int  N_anker_group = tab_ank && tab_ank->N_group > 1 ? tab_ank->N_group-1 : 0, lenStrSeparator = (int)strlen(STRSEPARATOR), 
			  N_group = 0;

		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

/////////////////////////////////////
//...порождаем новую таблицу тяжений;
			Table * table = GetLineTensionTable(N_group ? N_group : N_anker_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading dampher parameters (reading or setting table) ================================================
// ======================================================================================================
char * LineDampherReading_Aeol(char * pchar, void * context, void * ankertab)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

///////////////////////////////////////////////////////////////
//...зачитываем параметры демпфирующих схем защиты от вибрации;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

////////////////////////////////////////////////////////////////
//...зачитываем число строчек в таблице демпфирующих параметров;
		Table * tab_ank = (Table *)GetTable(ankertab);
		int  N_anker_group = tab_ank ? tab_ank->N_group : 0, lenStrSeparator = (int)strlen(STRSEPARATOR), 
			  N_group = 0;

		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

/////////////////////////////////////////
//...порождаем новую таблицу демпфирующих параметров;
			Table * table = GetLineDampherTable_Aeol(N_group ? N_group : N_anker_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Loading dampher data base ============================================================================
// ======================================================================================================
char * DampherDataBaseReading(char * pchar, void * context, int head_reading)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);
	int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);

///////////////////////////////////////////////////////////////
//...зачитываем имя и описание элемента демпферной базы данных;
	if (((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) && head_reading) {
		if ((pnext = strstr(pchar, STRSEPARATOR)) != NULL || (pnext = p_end) != NULL) {
			char Temp0 = pnext[0];
			pnext[0]   = '\x0';

			SetSampleComment(context, pchar);

			pnext[0] = Temp0;
			if ((pchar = pnext) < p_end) pchar += lenStrSeparator;
		}
		if ( pchar < (pnext = p_end)) pnext -= lenStrSeparator;
		char Temp2 =  pnext[0];
		pnext[0]   = '\x0';

		SetSampleDescription(context, pchar);

		pnext[0] = Temp2;
		pchar = p_end;
	}

//////////////////////////////////////////////////////////
//...зачитываем параметры элемента демпферной базы данных;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

/////////////////////////////////
//...зачитываем число параметров;
		if (((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

/////////////////////////////////////////////////////////////
//...порождаем новую таблицу элемента демпферной базы данных;
			Table * table = GetDampherDataBaseTable(N_group);
			if (table) {
				SetTable (context, table);
				SetTableIndex(context, 1);

///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
						SetTableParamsAsString(table, i, pchar, pnext, " || ");
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
			}
			else p_end = pchar;
			if (head_reading) p_end = pchar;//...оставляем конец на последем зачитывании;
		}
	}
	return(p_end);
}

// ======================================================================================================
// Tower exploration reading  ===========================================================================
// ======================================================================================================
////////////////////////////////////////////////////////////////////
//...конвертация данных CSV таблиц поопорной ведомости в LPR формат;
void TowerCSVConvertToLPR(void * context, void ** anker_list, int mufta_num)
{
	Table * tab_csv = (Table *)GetTable(context);
	if (tab_csv) {
//////////////////////////////////////////////////////////
//...вычисляем число опор и распределяем анкерную таблицу;
		char * pchar, * pnext, * sloc = _strdup(setlocale(LC_COLLATE, NULL));
		setlocale(LC_COLLATE, "rus"); //...сравниваем без учета регистра;
		int N_group = 0, i;
		for ( i = 1; i <= tab_csv->N_group; i++)
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && ! _stricmp("шифр", pchar)) N_group++;
		
		SetTable(anker_list[NUMANKERTAB], GetLineTable(N_group));

////////////////////////////////////////////////////////////////////
//...вычисляем число пересечений и распределяем таблицу пересечений;
		for ( N_group = 0, i = 1; i <= tab_csv->N_group; i++)
			if ((pchar = ((char **)tab_csv->table[5].parm_values)[i]) != NULL && ! _stricmp("да", pchar)) {
				int i_value = 0;
				if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i]) != NULL && sscanf(pchar, "%i", &i_value) == 1)
					N_group += i_value;
			}

		SetTable(anker_list[NUMSECTIONTAB], GetLineSectTable(N_group));

//////////////////////////////////////////////////////
//...записываем название линии (и пустой комментарий);
		SetSampleComment    (anker_list[NUMANKERTAB], (char *)(long long)GetTableParam(context, 2, 0));
		SetSampleDescription(anker_list[NUMANKERTAB], "");

////////////////////////////////
//...заполняем анкерную таблицу;
		int i_value, j = 0, l = 0, m = mufta_num;
		for (i = 1; i <= tab_csv->N_group; i++) {
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && ! _stricmp("шифр", pchar)) {
				j++;
				SetTableParam(anker_list[NUMANKERTAB], j, 0, GetTableParam(context, i-3, 1));//...проектный номер;
				SetTableParam(anker_list[NUMANKERTAB], j, 2, GetTableParam(context, i-2, 1));//...эксплуатационный номер;
				SetTableParam(anker_list[NUMANKERTAB], j, 3, GetTableParam(context, i+2, 6));//...название фотографии;
				if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i-1]) != NULL && sscanf(pchar, "%i", &i_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT, i_value); else
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT, -1);

				double d_value = 0.;
				if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-5, d_value);//...длина пролета вперед;
				if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i-2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && j > 1)
				SetTableParam(anker_list[NUMANKERTAB], j-1, NUMCIRCUIT-3, -d_value);//...разность высот назад;

///////////////////////////////////
//...отмечаем анкер на месте муфты;
				if ( ( pchar = ((char **)tab_csv->table[10].parm_values)[i-3]) != NULL && (pnext = strstr(pchar, "анкер")) != NULL)
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-1, (double)(long long)"анкер");

/////////////////////////////////////////////////////////////////////////////////////////////////////
//...корректируем шифр опоры -- убираем все пробелы и нумеруем муфты;
				if ((pchar = ((char **)tab_csv->table[1].parm_values)[i]) != NULL)
				while ((pnext = strstr(pchar, " ")) != NULL) {
						memmove(pnext, pnext+1, strlen(pnext)+1);
						pchar = pnext;
				}
				if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i]) != NULL   && (pnext = strstr(pchar, "(z)")) != NULL ||
					  (pchar = ((char **)tab_csv->table[1].parm_values)[i-2]) != NULL && (pnext = strstr(pchar, "(z)")) != NULL) {
					char buff[200];
					if ((pchar = ((char **)tab_csv->table[10].parm_values)[i-3]) != NULL && (pnext = strstr(pchar, "анкер")) != NULL)
						sprintf(buff, "анкер М%i", m++); else
						sprintf(buff, "М%i", m++);
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-1, (double)(long long)buff);
					sprintf(buff, "№%i", m-1);
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT+2, (double)(long long)buff);//...шифр муфты;
				}
				SetTableParam(anker_list[NUMANKERTAB], j, 1, GetTableParam(context, i, 1));//...шифр опоры;


/////////////////////////////
//...разбираем угол поворота;
				if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i+1]) != NULL) {
						d_value = fabs(strtod(pchar, & pnext))*60.;
						if (pnext && '\xB0' == pnext[0]) { //...символ градусов -- °;
							 pchar = pnext+1;
							 double f = fabs(strtod(pchar, & pnext));
							 if (pnext && ('\x60' == pnext[0] || '\x27' == pnext[0])) { //...символ минут -- ` или ';
								  pchar = pnext+1;
								  d_value += f;
							 }
						}
						if (pnext && strstr(pnext, "право") ) d_value = -d_value;
						SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-6, d_value);//...закодированный угол поворота;
				}

///////////////////////////////////////////////////////////
//...разбираем GPS координаты (N угол и E угол независимо);
				char * pname; 
				if ( (pnext = ((char **)tab_csv->table[3].parm_values)[i+1]) != NULL && (pchar = strstr(pnext, "N")) != NULL) {
						d_value = fabs(strtod(pchar+1, & pnext))*3600.;
						if (pnext && '\xB0' == pnext[0]) { //...символ градусов -- °;
							 pchar = pnext+1;
							 double f = fabs(strtod(pchar, & pnext));
							 if (pnext && ('\x60' == pnext[0] || '\x27' == pnext[0])) { //...символ минут -- ` или ';
								  pchar = pnext+1;
								  d_value += f*60;
								  if (pnext && (pname = strstr(pnext, ",")) != NULL) { //...проверяем десятичную запятую;
										pname[0] = '.';
										f = fabs(strtod(pchar, & pnext));
										pname[0] = ',';
								  }
								  else f = fabs(strtod(pchar, & pnext));
								  if (pnext && ! strncmp("\x60\x60", pnext, 2)) { //...символ секунд -- ``;
										pchar = pnext+2;
										d_value += f;
								  }
							 }
						}
						SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS, d_value);//...закодированная N-GPS координата;
				}
				if ( (pnext = ((char **)tab_csv->table[3].parm_values)[i+1]) != NULL && (pchar = strstr(pnext, "E")) != NULL) {
						d_value = fabs(strtod(pchar+1, & pnext))*3600.;
						if (pnext && '\xB0' == pnext[0]) { //...символ градусов -- °;
							 pchar = pnext+1;
							 double f = fabs(strtod(pchar, & pnext));
							 if (pnext && ('\x60' == pnext[0] || '\x27' == pnext[0])) { //...символ минут -- ` или ';
								  pchar = pnext+1;
								  d_value += f*60;
								  if (pnext && (pname = strstr(pnext, ",")) != NULL) { //...проверяем десятичную запятую;
										pname[0] = '.';
										f = fabs(strtod(pchar, & pnext));
										pname[0] = ',';
								  }
								  else f = fabs(strtod(pchar, & pnext));
								  if (pnext && ! strncmp("\x60\x60", pnext, 2)) { //...символ секунд -- ``;
										pchar = pnext+2;
										d_value += f;
								  }
							 }
						}
						SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS+1, d_value);//...закодированная Е-GPS координата;
				}

///////////////////////////////////////////////////////////////
//...заполняем индекс местности и индекс климатических условий;
				if ((pchar = ((char **)tab_csv->table[9].parm_values)[i-3]) != NULL && strlen(pchar) > 0) {
					unsigned char T[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
					T[0] = strstr(pchar, "1-1") != NULL;
					T[1] = strstr(pchar, "1-2") != NULL;
					T[2] = strstr(pchar, "1-3") != NULL;
					T[3] = strstr(pchar, "2-1") != NULL;
					T[4] = strstr(pchar, "2-2") != NULL;
					T[5] = strstr(pchar, "2-3") != NULL;
					T[6] = strstr(pchar, "3-1") != NULL;
					T[7] = strstr(pchar, "3-2") != NULL;
					T[8] = strstr(pchar, "3-3") != NULL;
					T[9] = strstr(pchar, "4-1") != NULL;
					T[10] = strstr(pchar, "4-2") != NULL;
					T[11] = strstr(pchar, "4-3") != NULL;
					T[12] = strstr(pchar, "5-1") != NULL;
					T[13] = strstr(pchar, "5-2") != NULL;
					T[14] = strstr(pchar, "5-3") != NULL;

					unsigned long terra = T[0];
					for (int k = 1; k < 15; k++) 
						terra = 2*terra+T[k];
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT+9, (double)terra);
				}
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT+8, 1.);

////////////////////////////////////
//...заполняем комментарий по опоре;
				int first = 1, comment_length = 0, k;
				for (k = -3; k <= 2; k++)
				if ((pchar = ((char **)tab_csv->table[11].parm_values)[i+k]) != NULL && strlen(pchar) > 0)
					comment_length += (int)strlen(pchar)+2;

				char * comment = new_struct((comment_length+1)*sizeof(char));
				for ( k = -3; k <= 2; k++)
				if ((pchar = ((char **)tab_csv->table[11].parm_values)[i+k]) != NULL && strlen(pchar) > 0 && strstr(pchar, " ")) {
					if (! first)
					strcat(comment, ";\xA");
					strcat(comment, pchar);
					first = 0;
				}
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT+10, (double)(long long)comment);
				delete_struct(comment);
			}

/////////////////////////////////////////////////////////////////////////////////////////
//...заполняем позицию пересечений среди всех опор (пересечения берутся из правой опоры);
			if ((pchar = ((char **)tab_csv->table[5].parm_values)[i]) != NULL && ! _stricmp("да", pchar) && j > 1) {
				i_value = 0;
				if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i]) != NULL && sscanf(pchar, "%i", &i_value) == 1)
				for (int k = 0; k < i_value; k++) {
					l++;
					SetTableParam(anker_list[NUMSECTIONTAB], l, 0, j);//...шифр опоры берется из левой, предыдущей опоры;
					SetTableParam(anker_list[NUMSECTIONTAB], l, NUMSCTNAME, GetTableParam(context, i+2, 4));
					SetTableParam(anker_list[NUMSECTIONTAB], l, NUMSCTTYPE+8, (double)(long long)"");//...комментарий к пересечениям;
				}
			}
		}
		setlocale(LC_COLLATE, sloc); free(sloc);
	}
}

//////////////////////////////////////////////////////////////////
//...конвертация (коррекция) данных на основе поопорной ведомости;
void TowerCSVConvertToLPR(void * context, void ** anker_list, int beg_pos, int set_anker, int beg_pos_comment)
{
	Table * tab_csv = (Table *)GetTable(context);
	if (tab_csv) {

//////////////////////////
//...вычисляем число опор;
		char * pchar;
		int N_group = 0, i;
		for (i = beg_pos; i <= tab_csv->N_group; i++)
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && strlen(pchar) > 0) N_group++;

		if (set_anker) {
			SetTable(anker_list[NUMANKERTAB], GetLineTable(N_group));
			SetTable(anker_list[NUMSECTIONTAB], NULL);

//////////////////////////////////////////////////////
//...записываем название линии (и пустой комментарий);
			SetSampleComment    (anker_list[NUMANKERTAB], (char *)(long long)GetTableParam(context, 2, 0));
			SetSampleDescription(anker_list[NUMANKERTAB], "");
		}

////////////////////////////////
//...заполняем анкерную таблицу;
		int j = 0;
		if (N_group == GetNGroup(anker_list[NUMANKERTAB]))
		for (i = beg_pos; i <= tab_csv->N_group; i++) {
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && strlen(pchar) > 0) {
				j++;
				SetTableParam(anker_list[NUMANKERTAB], j, 0, GetTableParam(context, i, 0));//...проектный номер;
				SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS, 0.);//...N-GPS координата;
				SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS+1, 0.);//...E-GPS координата;

				double d_value = 0.;
				if ( i < tab_csv->N_group)
				if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i+1]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-5, d_value);//...длина пролета вперед;

				d_value = GetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-3);
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-3, d_value = 0.01*((int)(100.*d_value+.5)));

				if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-6, d_value*60.);//...угол поворота;

////////////////////
//...отмечаем анкер;
				if ( ( pchar = ((char **)tab_csv->table[1].parm_values)[i]) != NULL && ! strncmp(pchar, "У", 1))
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-1, (double)(long long)"анкер");
				else 
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-1, NULL);

////////////////////////////////////////////////////
//...корректируем шифр опоры -- убираем все пробелы;
				char * pnext;
				if (pchar)
				while ((pnext = strstr(pchar, " ")) != NULL) {
						memmove(pnext, pnext+1, strlen(pnext)+1);
						pchar = pnext;
				}
				SetTableParam(anker_list[NUMANKERTAB], j, 1, GetTableParam(context, i, 1));//...шифр опоры;

/////////////////////////////
//...накапливаем пересечения;
				if (i < tab_csv->N_group)
				if (1 && ( pchar = ((char **)tab_csv->table[5].parm_values)[i+1]) != NULL && strlen(pchar) > 0) {
					while ((pnext = strstr(pchar, ",")) != NULL) {
							char temp = pnext[0];
							pnext[0] = '\x0';
							if (! GetTable(anker_list[NUMSECTIONTAB])) {
								SetTable(anker_list[NUMSECTIONTAB], GetLineSectTable(1));
								SetTableParam(anker_list[NUMSECTIONTAB], 1, 0, j);
								SetTableParam(anker_list[NUMSECTIONTAB], 1, NUMSCTNAME, (double)(long long)pchar);
							}
							else {
								AddTableParameters (anker_list[NUMSECTIONTAB], 1);
								N_group = GetNGroup(anker_list[NUMSECTIONTAB]);
								SetTableParam(anker_list[NUMSECTIONTAB], N_group, 0, j);
								SetTableParam(anker_list[NUMSECTIONTAB], N_group, NUMSCTNAME, (double)(long long)pchar);
							}
							pnext[0] = temp;
							pchar = pnext+1;
					}
					if (! GetTable(anker_list[NUMSECTIONTAB])) {
						SetTable(anker_list[NUMSECTIONTAB], GetLineSectTable(1));
						SetTableParam(anker_list[NUMSECTIONTAB], 1, 0, j);
						SetTableParam(anker_list[NUMSECTIONTAB], 1, NUMSCTNAME, (double)(long long)pchar);
					}
					else {
						AddTableParameters (anker_list[NUMSECTIONTAB], 1);
						N_group = GetNGroup(anker_list[NUMSECTIONTAB]);
						SetTableParam(anker_list[NUMSECTIONTAB], N_group, 0, j);
						SetTableParam(anker_list[NUMSECTIONTAB], N_group, NUMSCTNAME, (double)(long long)pchar);
					}
				}
			}
		}
	}
}


////////////////////////////////////////////////////
//...конвертация (коррекция) данных по пересечениям;
void TowerCSVSetSectionsToLPR(void * context, void ** anker_list, int id_section, int id_inverse)
{
	Table * tab_csv = (Table *)GetTable(context);
	if (tab_csv) {

///////////////////////////
//...заполняем пересечение;
		char * pchar;
		for (int i = 1; i <= tab_csv->N_group; i++) {
			if ( ( pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && strlen(pchar) > 0) {
				double d_value = 0.;
				if (strstr(pchar, "Наименование")) 
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTNAME, GetTableParam(context, i, 2));//...наименование пересечения;
				else
				if (strstr(pchar, "Дата измерений")) {
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTNAME+1, GetTableParam(context, i, 3));//...дата измерений;
					if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTNAME+2, d_value);//...температура;
				}
				else
				if (strstr(pchar, "Продольный профиль")) {

/////////////////////////////////////
//...вставляем параметры пересечения;
					if (! id_inverse) {
						if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-5, d_value);//...H_L_Podves;

						if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-4, d_value);//...H_R_Podves;
					}
					else {
						if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-4, d_value);//...H_R_Podves;

						if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-5, d_value);//...H_L_Podves;
					}
					if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, 2, d_value); else 
					if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, 2, d_value); else 
					if ( (pchar = ((char **)tab_csv->table[5].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, 2, d_value);//...H_Obj;

					if ( (pchar = ((char **)tab_csv->table[6].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-2, d_value); else
					if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-2, d_value); else 
					if ( (pchar = ((char **)tab_csv->table[8].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE-2, d_value);//...H_Provod;
				}
				else
				if (strstr(pchar, "План")) {
//////////////////////////////////////////////////////////
//...вставляем неизвесный тип пересечения (SCTTYPE_STATE);
					SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE, SCTTYPE_STATE);//...Type
					if (! id_inverse) {
						if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, d_value); else
						if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, -d_value); else
						if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, -d_value);//...L_to_Obj;

						if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+1, d_value); else
						if ( (pchar = ((char **)tab_csv->table[8].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+1, d_value);//...UU

						if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+4, d_value); else
						if ( (pchar = ((char **)tab_csv->table[6].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+4, d_value);//...L1

						if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+5, d_value); else
						if ( (pchar = ((char **)tab_csv->table[5].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+5, d_value);//...L2
					}
					else { //...угол не меняется при переворачивании линии;
						if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, -d_value); else
						if ( (pchar = ((char **)tab_csv->table[1].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, d_value); else
						if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, 1, d_value);//...L_to_Obj;

						if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+1, d_value); else
						if ( (pchar = ((char **)tab_csv->table[8].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+1, d_value);//...UU;

						if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+4, d_value); else
						if ( (pchar = ((char **)tab_csv->table[5].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+4, d_value);//...L1;

						if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+5, d_value); else
						if ( (pchar = ((char **)tab_csv->table[6].parm_values)[i+2]) != NULL && sscanf(pchar, "%lg", &d_value) == 1 && d_value != 0.)
						SetTableParam(anker_list[NUMSECTIONTAB], id_section, NUMSCTTYPE+5, d_value);//...L2;
					}
				}
			}
		}
	}
}


// ======================================================================================================
// Simple table reading  ================================================================================
// ======================================================================================================
/////////////////////////////////////////////////////////////////////
//...конвертация данных CSV таблиц из "простого" списка в LPR формат;
void SimpleCSVConvertToLPR(void * context, void ** anker_list, int mufta_num)
{
	Table * tab_csv = (Table *)GetTable(context);
	if (tab_csv) {

//////////////////////////////////////////////////////////
//...вычисляем число опор и распределяем анкерную таблицу;
		char * pchar;
		int N_group = 0, i;
		for (i = 5; i <= tab_csv->N_group; i++)
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && strlen(pchar) > 0) N_group++;

		SetTable(anker_list[NUMANKERTAB], GetLineTable(N_group));
		SetTable(anker_list[NUMSECTIONTAB], NULL);

//////////////////////////////////////////////////////
//...записываем название линии (и пустой комментарий);
		SetSampleComment    (anker_list[NUMANKERTAB], (char *)(long long)GetTableParam(context, 2, 2));
		SetSampleDescription(anker_list[NUMANKERTAB], "");

////////////////////////////////
//...заполняем анкерную таблицу;
		int j = 0;
		for (i = 5; i <= tab_csv->N_group; i++) {
			if ( (pchar = ((char **)tab_csv->table[0].parm_values)[i]) != NULL && strlen(pchar) > 0) {
				j++;
				SetTableParam(anker_list[NUMANKERTAB], j, 0, GetTableParam(context, i, 0));//...проектный номер;
				SetTableParam(anker_list[NUMANKERTAB], j, 1, GetTableParam(context, i, 1));//...шифр опоры;
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT, 2);

				double d_value = 0.;
				if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-5, d_value);//...длина пролета вперед;

				if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &d_value) == 1)
				SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-6, d_value*60.);//...угол поворота;

////////////////////
//...отмечаем анкер;
				char * pnext;
				if ( ( pchar = ((char **)tab_csv->table[1].parm_values)[i]) != NULL && ! strncmp(pchar, "У", 1))
					SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-1, (double)(long long)"анкер");

/////////////////////////////
//...разбираем угол поворота;
				if ( (pchar = ((char **)tab_csv->table[5].parm_values)[i]) != NULL) {
						d_value = strtod(pchar, & pnext)*60.;
						SetTableParam(anker_list[NUMANKERTAB], j, NUMCIRCUIT-6, d_value);//...закодированный угол поворота;
				}

///////////////////////////////////////////////////////////
//...разбираем GPS координаты (N угол и E угол независимо);
				if ( (pchar = ((char **)tab_csv->table[2].parm_values)[i]) != NULL) {
						d_value = strtod(pchar, & pnext);
						SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS, d_value);//...закодированная N-GPS координата;
				}
				if ( (pchar = ((char **)tab_csv->table[3].parm_values)[i]) != NULL) {
						d_value = strtod(pchar, & pnext);
						SetTableParam(anker_list[NUMANKERTAB], j, NUMGPS_POS+1, d_value);//...закодированная Е-GPS координата;
				}
			}
		}

///////////////////////////////////////////
//...перепад уровня земли (разность высот);
		for (i = 1; i < N_group; i++)
			SetTableParam(anker_list[NUMANKERTAB], i, NUMCIRCUIT-3, 
			GetTableParam(anker_list[NUMANKERTAB], i, NUMCIRCUIT-3)-GetTableParam(anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3));
	}
}

//////////////////////////////////////////////////////////////////////
//...установка тяжения провода по данным CSV таблиц "простого" списка;
void SimpleCSVSetTensionToLPR(void * context, void ** anker_list, double t_norm)
{
	Table * tab_csv = (Table *)GetTable(context);
	if (tab_csv) {

//////////////////////////////////////////////////////////////////////
//...вычисляем тяжения по стрелам провеса и заполняет таблицу тяжений;
		char * pchar, * sss;
		for (int j = 1, i = 5; i <= tab_csv->N_group; i++) {
			if ((pchar = ((char **)tab_csv->table[1].parm_values)[i]) != NULL && strlen(pchar) > 0  &&
				 (sss = (char *)(long long)GetTableParam(anker_list[NUMANKERTAB], j, 1)) != NULL && ! ::strcmp(sss, pchar)) {
				double f_strela = 0.;
				if ( (pchar = ((char **)tab_csv->table[4].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &f_strela) == 1) {
					//...f_strela -- стрела провеса верхней траверсы;
					f_strela = f_strela;
				}
				if ( (pchar = ((char **)tab_csv->table[7].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &f_strela) == 1) {
					//...f_strela -- стрела провеса верхней траверсы;
					f_strela = f_strela;
				}
				if ( (pchar = ((char **)tab_csv->table[10].parm_values)[i]) != NULL && sscanf(pchar, "%lg", &f_strela) == 1) {
					//...f_strela -- стрела провеса верхней траверсы;
					f_strela = f_strela;
				}
				j++;
			}
		}
	}
}
#endif
