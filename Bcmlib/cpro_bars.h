/*============================================*/
/*                  PRO_BARS                  */
/*============================================*/
#ifndef ___PRO_BARS___
#define ___PRO_BARS___
#define ___PRO_BARS2D_cpp___
#include "cgrid.h"

///////////////////////////////////////////////////////////////
//...число 2D статических образцов-сечений из ГОСТА и их имена;
int    GetBARSGOSTSampleCount(void);
char * GetBARSGOSTSampleName (int N_sm);

/////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для статических образцов библиотеки ГОСТА;
void * CreateBARSGOSTContext(int N_sm);

/////////////////////////////////////////////////////////
//...число 2D статических образцов для Bars25 и их имена;
int    GetBARSSampleCount(void);
char * GetBARSSampleName (int N_sm);

///////////////////////////////////////////////////////////////////////
//...функция образования контекста для статических образцов библиотеки;
void * CreateBARSContext(int N_sm);

////////////////////////////////////////////////////////////////////////////////////
//...пpедваpительная инициализация  статического обpазца для Bars25 и ГОСТа сечений;
int bars2D_init(void * context, CGrid * block_nd = NULL);

///////////////////////////////////////////
//...pешатель функций фоpмы для скин-задач;
int bars2D_solv(void * context, int & sec, int & hund, int id_solv);

//////////////////////////////////////////////////////////
//...результаты для 2D образца (задачи пограничного слоя);
char * GetFuncBARSName(int num);
enum Func_Bars {
    _FBARS1 = 0,  
    _FBARS2,  
    _FBARS3,  
    _FBARS4,  
    _FBARS5,  
    _FBARS6,
    _NUM_BARS_FUNCTIONS
};
#define NUM_BARS_FUNCTIONS _NUM_BARS_FUNCTIONS
inline int GetFuncBARSCount  (void) { return(NUM_BARS_FUNCTIONS);}
inline int GetFuncBARSDefault(void) { return(_FBARS1);}

////////////////////////////////////////
//...рассчет результатов для скин-задач;
void solver_Bars(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncBARSDefault());

/////////////////////////////////////////////////////////////////
//...интерфейсные функции для взаимодействия с внешними модулями;
double get_square_bars(void * context);
double get_inertia_Y_bars(void * context);
double get_inertia_Z_bars(void * context);
double get_inertia_YZ_bars(void * context);
double get_static_Y_bars(void * context);
double get_static_Z_bars(void * context);
double get_torsion_J_bars(void * context);
double get_E_core_bars(void * context);
double get_nju_core_bars(void * context);

#endif
