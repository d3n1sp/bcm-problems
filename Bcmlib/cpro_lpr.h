/*===========================================*/
/*                  PRO_LPR                  */
/*===========================================*/
#ifndef ___PRO_LPR___
#define ___PRO_LPR___
#define ___PRO_LPR_cpp___
#include "cprofile.h"

////////////////////////////////////////////////////
//...нумерация таблиц в проекте электрической линии;
enum Num_table_LineProject {
		NUMANKERTAB = 0,    //...для загрузки таблицы проектируемой линии; 
		NUMTOWERTAB,        //...для загрузки таблицы параметров опор;
		NUMANKERHANGINGTAB, //...для загрузки таблицы условий подвески кабеля и провода на анкерной опоре;
		NUMSTANDHANGINGTAB, //...для загрузки таблицы условий подвески кабеля и провода на поддерживающей опоре;
		NUMSECTIONTAB,      //...для загрузки таблицы пересечений линии;
		NUMPROFILETAB,      //...для загрузки таблицы профиля линии;
		NUMCLIMATETAB,      //...для загрузки таблицы климатических условий линии;
		NUMCABLEMODETAB,    //...для загрузки таблицы условий нагружения кабеля;     
		NUMGROZOMODETAB,    //...для загрузки таблицы условий нагружения грозотроса;
		NUMWIRE1MODETAB,    //...для загрузки таблицы условий нагружения провода;
		NUMCABLEMONTAGETAB, //...для хранения монтажной таблицы кабеля;
		NUMGROZOMONTAGETAB, //...для хранения монтажной таблицы грозотроса;
		NUMWIRE1MONTAGETAB, //...для хранения монтажной таблицы провода;
		NUMCABLESCHEMETAB,  //...для загрузки таблицы схем виброзащиты для кабеля;
		NUMGROZOSCHEMETAB,  //...для загрузки таблицы схем виброзащиты для грозотроса;
		NUMWIRE1SCHEMETAB,  //...для загрузки таблицы схем виброзащиты для провода;
		NUMTENSIONTAB,      //...для загрузки таблицы тяжений для провода, кабеля и грозотроса;
		NUMCABLEDAMPHERTAB, //...для загрузки таблицы демпфирующих параметров для кабеля;
		NUMGROZODAMPHERTAB, //...для загрузки таблицы демпфирующих параметров для грозотроса;
		NUMWIRE1DAMPHERTAB, //...для загрузки таблицы демпфирующих параметров для провода;
		NUMLINEPROJECT,     //...общее число таблиц в проекте электрической линии;
};

////////////////////////////////////////////////////////
//...максимальное число кабелей, проводов и грозотросов;
const unsigned MAXCABLES  = 1;
const unsigned MAXGROZOS  = 1;
const unsigned MAXCIRCUIT = 6;

////////////////////////////////////////////
//...максимальное число гасителей на пролет;
const unsigned MAXDAMPHER = 4;

////////////////////////////////
//...позиции пунктов в таблицах;
const unsigned NUMGPS_POS = 4; //...позиция начала группы GPS координат (GPS_N) в таблице линии;
const unsigned NUMCIRCUIT = 13;//...позиция числа объектов на опоре (N_circuit -- "фазность") в таблице линии;
const unsigned NUMMODE_P1 = 5; //...позиция начала группы режимов нагружений (Mode_p1) в таблице режимов нагружений;
const unsigned NUMTENSION = 9; //...позиция начала группы тяжений (T_break) в таблице режимов нагружений;
const unsigned NUMEFMODUL = 16;//...позиция группы модулей (EF) в таблице режимов нагружений;
const unsigned NUMSCTNAME = 3; //...позиция наименования пересечения (Name) в таблице пересечений;
const unsigned NUMSCTTYPE = 11;//...позиция индекса пересечений (Type) в таблице пересечений;
const unsigned NUMAC70_72 = 9; //...позиция провода АС-70/72 (грозотроса!) в таблице проводов типа АС (_WIRE4);
const unsigned NUMSTENPOS = 5; //...число параметров для одного объекта в таблице тяжений;
const unsigned NUMMAXSPAN = 5; //...позиция максимальной длины защищаемого пролета в демпферной таблице;

//////////////////////////////////////////////////////////////////////////////////////////////////
//...максимальное число позиций в базе данных опто-волоконных кабелей и в описании геометрии опор;
const unsigned NUM_OF_CABLE_BASE_ELEMENTS = 8;
const unsigned NUM_OF_TOWER_BASE_ELEMENTS = 22;

////////////////////////////
//...готовые шаблоны таблиц;
Table * GetLineTable (int N_group);
Table * GetTowerTable(int N_group);
Table * GetLineHangTable(int N_group);
Table * GetLineSectTable(int N_group);
Table * GetLineProfileTable(int N_group, int N_points = 0);
Table * GetLineClimateTable(int N_group);
Table * GetLineModeTable(int N_group);
Table * GetLineTensionTable(int N_group);
Table * GetLineDampherTable(int N_group);
Table * GetLineDampherTable_Aeol(int N_group);
Table * GetMontageTable (int N_group, int N, double * temper = NULL);

///////////////////////////////////////////////////
//...шаблон таблиц для задания предельных нагрузок;
Table * GetTenLimitTable(int N_group = 1, char * table_names[] = NULL);

//////////////////////////////////////////
//...таблицы для представления баз данных;
Table * GetWireDataBaseTable   (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetCableDataBaseTable  (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetTowerDataBaseTable	 (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetDampherDataBaseTable(int N_group);
Table * GetDampherSchemeTable	 (int N_group, int id_name = OK_STATE);

////////////////////////////////////////////////////////////////
// ...эксплуатационный номер опоры и определение анкерной опоры;
char * operating_num(Table * table, int i, char * buff, int num_project = NULL_STATE);
int	 id_anker_num (Table * table, int i, int id_tower = 0);

////////////////////////////////////////
// ...корректировка пустых комментариев;
void comment_correct(Table * table, int id_static_char = NULL_STATE);
void comment_correct(void * context);

////////////////////////////////////
// ...корректировка таблиц подвески;
void hunging_correct(Table * table, Table * tab_ank, int * anker_num, int id_anker = 1);

//////////////////////////////////////////////////////////////
// ...корректировка таблицы: выбрасываем повторяющиеся строки;
void table_double_correct(Table * table, int excl_pos = ERR_STATE, int id_static_char = NULL_STATE, int id_sequent = OK_STATE);
void table_double_correct(Table * table, Table * tab_ank, int id_static_char = NULL_STATE);

//////////////////////////////////////////////////////////////////
//...переустанавка линии для анализа схем виброзащиты провода;
void ResetLineForDampher(void ** anker_list);

////////////////////////////////////////////////////
//...вспомогательные функции для зачитывания таблиц;
char * LineReading       (char * pchar, void * context, int head_reading = NULL_STATE);
char * LineTowerReading  (char * pchar, void * context);
char * LineHangReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineSectReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineProfileReading(char * pchar, void * context, int convert = NULL_STATE);
char * LineClimateReading(char * pchar, void * context);
char * LineModeReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineMontageReading(char * pchar, void * context);
char * LineSchemeReading (char * pchar, void * context);
char * LineDampherReading(char * pchar, void * context, void * ankertab = NULL);
char * LineTensionReading(char * pchar, void * context, void * ankertab = NULL);
char * LineDampherReading_Aeol(char * pchar, void * context, void * ankertab = NULL);
char * DampherDataBaseReading (char * pchar, void * context, int head_reading = NULL_STATE);

////////////////////////////////////////////////////////////////////
//...конвертация данных CSV таблиц поопорной ведомости в LPR формат;
void TowerCSVConvertToLPR		(void * context, void ** anker_list, int mufta_num = OK_STATE);
void TowerCSVConvertToLPR		(void * context, void ** anker_list, int beg_pos, int set_anker, int beg_pos_comment);
void TowerCSVSetSectionsToLPR (void * context, void ** anker_list, int id_section, int id_inverse = NULL_STATE);

void SimpleCSVConvertToLPR		(void * context, void ** anker_list, int mufta_num = OK_STATE);
void SimpleCSVSetTensionToLPR	(void * context, void ** anker_list, double t_norm = 18.);

#endif
