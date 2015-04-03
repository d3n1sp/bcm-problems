#include "stdafx.h"
#include "cpro_bpr.h"

#ifdef ___PRO_BPR_cpp___
///////////////////////////
//...готовые шаблоны таблиц;
Table * GetTopoTable(int N_group, int N_points)
{
	Shablon records[] = {{CHAR_TYPE_RECORD, NULL},//(0) -- TypeBlock
								{ INT_TYPE_RECORD, NULL},//(1) -- ExternalNumber
								{ INT_TYPE_RECORD, NULL},//(2) -- PropNum
								{ INT_TYPE_RECORD, NULL},//(3) -- GridNum
								{ ERR_TYPE_RECORD, NULL},//(4) -- признак продолжения...
	};
	Table * table = get_shablon_table(N_points+3, N_group, records, NULL);
	return table;
}

Table * GetGridTable(int N_group)
{
	Shablon records[] = {{	  INT_TYPE_RECORD, NULL},//(0) -- ExternalNumber
								{	  INT_TYPE_RECORD, NULL},//(1) -- TypeCoord
								{DLENGTH_TYPE_RECORD, NULL},//(2) -- XCoord 
								{DLENGTH_TYPE_RECORD, NULL},//(3) -- YCoord 
								{DLENGTH_TYPE_RECORD, NULL},//(4) -- ZCoord 
								{	  INT_TYPE_RECORD, NULL},//(5) -- PropNum 
	};
	Table * table = get_shablon_table(6, N_group, records, NULL);
	for (int k = 1; table && k <= table->N_group; k++)
		((int *)(table->table[0].parm_values))[k] = k;
	return table;
}

//////////////////////////////////////////////
//...число позиций в описании топологии блока;
int topo_num(char * element, int & N_points)
{
	if (! element || ! strlen(element)) return(0);

	char * sloc = _strdup(setlocale(LC_COLLATE, NULL));
	setlocale(LC_COLLATE, "rus"); //...сравниваем без учета регистра;

	if (! _stricmp("CROD",   element) || ! _stricmp("CBEAM",  element)) N_points = 2; else 
	if (! _stricmp("CTRIA3", element)) N_points = 3; else 
	if (! _stricmp("CQUAD4", element) || ! _stricmp("CTETRA", element)) N_points = 4; else 
	if (! _stricmp("CPENTA", element)) N_points = 5; else 
	if (! _stricmp("CHEXA",  element)) N_points = 6; else N_points = 0;

	setlocale(LC_COLLATE, sloc); free(sloc);
	return(1);
}

// ======================================================================================================
// Topo reading  ========================================================================================
// ======================================================================================================
char * BlockTopoReading(char * pchar, void * context, int head_reading)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);
	int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);

///////////////////////////////////////
//...зачитываем имя и описание проекта;
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
	}
	pchar = p_end;

//////////////////////////////////////////
//...зачитываем параметры тополгии блоков;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

/////////////////////////////
//...зачитываем число блоков;
		if (((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

//////////////////////////////////////////////
//...порождаем новую таблицу топологии блоков;
			Table * table = GetTopoTable(N_group);
			if (table) {
///////////////////////////////////////
//...заполняем таблицу данными с диска;
				int i = 0;
				while (++i <= table->N_group && 
						((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL)) {
/////////////////////////////////
//...выравниваем строчки таблицы;
						int id_align = 1, N_points = 0;
						if (id_align) {
							char * pname, Temp = 0;
							Shablon records[] = {{INT_TYPE_RECORD, NULL},//(0) -- GridNum
														{ERR_TYPE_RECORD, NULL},//(1) -- признак продолжения...
							};
							if ( (pname = strstr(pchar, SEMICOLONSTR)) != NULL) Temp = pname[0];
							if (! pname || topo_num(pchar+strlen(COLONSTR)+1, N_points) != 1) N_points = 0;
							if (  pname) pname[0] = Temp;
							if (N_points > table->N-3) 
								add_table_param(table, N_points-table->N+3, records);
						}
////////////////////////////////
//...зачитываем строчку таблицы;
						SetTableParamsAsString(table, i, pchar, pnext);
						pchar = pnext;
				}
				// pchar -- имя профиля (занести в таблицу ???);
				SetTable(context, table);
			}
		}
	}
	return(p_end);
}

// ======================================================================================================
// Grid reading  ========================================================================================
// ======================================================================================================
char * BlockGridReading(char * pchar, void * context)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);

//////////////////////////////////////
//...зачитываем параметры сетки узлов;
	if ( pchar < p_eof) pchar += strlen(TABLESTR);
	if ((p_end = strstr(pchar, TABLESTR)) != NULL || (p_end = p_eof) != NULL) {

////////////////////////////////////////////////////
//...зачитываем число строчек в таблице сетки узлов;
		int  N_group = 0, lenStrSeparator = (int)strlen(STRSEPARATOR);
		if ( ((pnext = strstr(pchar, STRSEPARATOR)) != NULL && (pnext += lenStrSeparator) <= p_end || (pnext = p_end) != NULL) &&
			(sscanf(pchar, "%i", &N_group) == 1 || pnext == p_end)) {
			pchar = pnext;

/////////////////////////////////////////
//...порождаем новую таблицу сетки узлов;
			Table * table = GetGridTable(N_group);
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
#endif
