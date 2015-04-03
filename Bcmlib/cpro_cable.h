/*============================================*/
/*                  PRO_CABLE                 */
/*============================================*/
#ifndef ___PRO_CABLE___
#define ___PRO_CABLE___
#define ___PRO_CABLE2D_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

////////////////////////////////////////////////////
//...число статических образцов для опор и их имена;
int    GetCABLESampleCount(void);
char * GetCABLESampleName (int N_sm);

///////////////////////////////////////////////////////////////////////
//...функция образования контекста для статических образцов библиотеки;
void * CreateCABLEContext(int N_sm);

/////////////////////////////////////////////////////////////////
//...пpедваpительная инициализация статического обpазца для опор;
#ifndef ___ABRIDGE_PROFILE_MODE___
int cable2D_init (void * context, CGrid * block_nd = NULL);
#endif

///////////////////////////////////////////////
// ...база данных имеющихся параметров кабелей;
void CableBase(void ** cable_list, char * cable, double * par);

/////////////////////////////////////
//...pешатель функций фоpмы для опор;
int cable3D_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////
//...результаты для 2D образца (задачи для опор);
char * GetFuncCABLEName(int num);
enum Func_Cable {
    _FCABLE1 = 0,  
    _FCABLE2,  
    _FCABLE3,  
    _NUM_CABLE_FUNCTIONS
};
#define NUM_CABLE_FUNCTIONS _NUM_CABLE_FUNCTIONS
inline int GetFuncCABLECount  (void) { return(NUM_CABLE_FUNCTIONS);}
inline int GetFuncCABLEDefault(void) { return(_FCABLE1);}

/////////////////////////////////////////////////////
//...получение результатов по решению задач для опор;
#ifndef ___ABRIDGE_PROFILE_MODE___
void solver_Cable(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncCABLEDefault());
#endif

#endif
