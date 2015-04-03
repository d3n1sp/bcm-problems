/*===========================================*/
/*                  PRO_BPR                  */
/*===========================================*/
#ifndef ___PRO_BPR___
#define ___PRO_BPR___
#define ___PRO_BPR_cpp___
#include "cprofile.h"

////////////////////////////////////////////////////
//...нумерация таблиц в проекте электрической линии;
enum Num_table_BlocksProject {
		NUMTOPOTAB = 0,	//...для загрузки таблицы топологии блоков;
		NUMGRIDTAB,			//...для загрузки таблицы опорных узлов блоков; 
		NUMBLOCKSPROJECT, //...общее число таблиц в проекте блочного метода мультиполей;
};

///////////////////////////////////////////////////////////
//...максимальное число позиций в описании топологии блока;
const unsigned NUM_OF_TOPO_ELEMENTS = 6;

////////////////////////////
//...готовые шаблоны таблиц;
Table * GetTopoTable(int N_group, int N_points = NUM_OF_TOPO_ELEMENTS);
Table * GetGridTable(int N_group);

//////////////////////////////////////////////
//...число позиций в описании топологии блока;
int topo_num(char * element, int & N_points);

////////////////////////////////////////////////////
//...вспомогательные функции для зачитывания таблиц;
char * BlockTopoReading(char * pchar, void * context, int head_reading = 0);
char * BlockGridReading(char * pchar, void * context);

#endif
