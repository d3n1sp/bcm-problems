/*============================================*/
/*                  PRO_FITTING               */
/*============================================*/
#ifndef ___PRO_FITTING___
#define ___PRO_FITTING___
#define ___PRO_FITTING_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

///////////////////////////////////////////////////
//...число типов защитного оборудования и их имена;
int    GetFITTINGSampleCount(void);
char * GetFITTINGSampleName (int N_sm);
char * GetFITTINGGOSTName   (int N_sm);

///////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для базы данных защитного оборудования проводов;
  void * CreateFITTINGContext(int N_sm);

////////////////////////////////////////////////////////////////////////
//...пpедваpительная инициализация статического обpазца для оборудования;
#ifndef ___ABRIDGE_PROFILE_MODE___
  int fitting_init(void * context, CGrid * block_nd = NULL);
#endif

/////////////////////////////////////////////
//...pешатель функций фоpмы для оборудования;
int fitting_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////////////////////
//...интерфейсные функции для взаимодействия с внешними модулями;

#endif
