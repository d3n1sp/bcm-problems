/*========================================*/
/*                 CBASE                  */
/*========================================*/
#ifndef ___CBASE___
#define ___CBASE___

#include "utils.h"

typedef double Param;

/////////////////////////////////////////////////////////
//...bases class for serialization of multipoles objects;
class CBase {
protected:
		Param * param; //...data base parameters;
		static unsigned long COD; //...data base current code key;
public:
		virtual int size_of_param() {	return(0);}
		void set_param(int k, Param value) { 
			if (0 <= k && k < size_of_param()) param[k] = value;
		}
		void get_param(int k, Param & value) { 
			if (0 <= k && k < size_of_param()) value = param[k];
		}
		Param get_param(int k) { 
			Param value = 0.; get_param(k, value); 
			return(value);
		}
//...constructor/destructor;
		CBase() {
			param = NULL;
		}
      virtual ~CBase(void)	{
			delete_struct(param);
		}
};
#endif
