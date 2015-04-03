/*============================================*/
/*                  PRO_TOWER                 */
/*============================================*/
#ifndef ___PRO_TOWER___
#define ___PRO_TOWER___
#define ___PRO_TOWER3D_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

////////////////////////////////////////////////////
//...число статических образцов для опор и их имена;
int    GetTOWERSampleCount(void);
char * GetTOWERSampleName (int N_sm);

///////////////////////////////////////////////////////////////////////
//...функция образования контекста для статических образцов библиотеки;
void * CreateTOWERContext(int N_sm);

/////////////////////////////////////////////////////////////////
//...пpедваpительная инициализация статического обpазца для опор;
#ifndef ___ABRIDGE_PROFILE_MODE___
int tower3D_init (void * context, CGrid * block_nd = NULL);
#endif
enum Tower_Topo {
	_TOPO_ONE_CIRCLE_TREE = 1,  
	_TOPO_TWO_CIRCLE_TREE,  
	_TOPO_SPLIT_CIRCLE_TREE,  
	_TOPO_LEFT_GROZO_TREE,  
	_TOPO_HORIZONTAL_CROSS_TREE,  
	_TOPO_ONE_CIRCLE_RUMKA,

	_TOPO_ORDINARY_DUPLEX,
	_TOPO_ONE_CIRCUIT_DUPLEX,
	_TOPO_ONE_CIRCUIT_DUPLEX_A,
	_TOPO_ONE_CIRCUIT_DUPLEX_B, //...выяснить геометрию опор УС110-3, УС110-5, УС330-2, ПС330-6, ПС35-4;

	_TOPO_ORDINARY_TRIPLEX,  
	_TOPO_ONE_CIRCUIT_TRIPLEX,

	_TOPO_SIMPLE_TREE,    //...независимые столбы;
	_TOPO_SIMPLE_DUPLEX,  //...независимые столбы;
	_TOPO_SIMPLE_TRIPLEX, //...независимые столбы;
	_NUM_TOWER_TOPOS,
};
 
/////////////////////////////////////////////////////////////////////
// ...база данных стандартных опор (включая дополнительные проверки);
void TowerBase(void ** tower_list, char * tower, double * par, void * context = NULL);

/////////////////////////////////////
//...pешатель функций фоpмы для опор;
int tower3D_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////
//...результаты для 2D образца (задачи для опор);
char * GetFuncTOWERName(int num);
enum Func_Tower {
    _FTOWER1 = 0,  
    _FTOWER2,  
    _FTOWER3,  
    _FTOWER4,  
    _FTOWER5,  
    _FTOWER6,
    _NUM_TOWER_FUNCTIONS
};
#define NUM_TOWER_FUNCTIONS _NUM_TOWER_FUNCTIONS
inline int GetFuncTOWERCount  (void) { return(NUM_TOWER_FUNCTIONS);}
inline int GetFuncTOWERDefault(void) { return(_FTOWER1);}

/////////////////////////////////////////////////////
//...получение результатов по решению задач для опор;
#ifndef ___ABRIDGE_PROFILE_MODE___
void solver_Tower(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncTOWERDefault());
#endif

#endif
