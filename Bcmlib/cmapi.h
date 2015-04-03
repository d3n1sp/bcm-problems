/*========================================*/
/*                 CMAPI                  */
/*========================================*/
#ifndef ___CMAPI___
#define ___CMAPI___

#include "cshapes.h"

//////////////////////////////////////////////
//														  //
//      SYSTEMS of SPHEROID MULTIPOLES      //
//														  //
//////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//...class of acoustic multipoles, based on three plane waves;
class CMapi3DSpheroid : public CShape<double>  {
public:
		int type()          { return MP3D_SPHEROID_SHAPE;}
		int size_of_param() { return(2);}
		int freedom (int m) { return 4;} //...вычисляем только четыре функции; 
		void set_shape(double R0, double eps = 1., double p1 = 1., double p2 = 0., double p3 = 0., double p4 = 0.);
public:
//...initialization and calculation of multipoles;
		void parametrization		 (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0) { parametrization(P, m);} 
		void parametrization_hess(double * P = NULL, int m = 0) { parametrization(P, m);} 
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CMapi3DSpheroid () {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

class CMapi3DSpheroidFull : public CShape<double>  {
public:
		int type()          { return MP3D_SPHEROID_FULL_SHAPE;}
		int size_of_param() { return(2);}
		int freedom(int m)  { return sqr(m+1);}
		void set_shape(double R0, double eps = 1., double p1 = 1., double p2 = 0., double p3 = 0., double p4 = 0.);
public:
//...initialization and calculation of multipoles;
		void parametrization		 (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0) { parametrization(P, m_dop);} 
		void parametrization_hess(double * P = NULL, int m_dop = 0) { parametrization(P, m_dop);} 
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CMapi3DSpheroidFull () {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_spheroid_eps					0. //...параметр сфероида -- соотношение полуосей;
#undef  SHAPE_spheroid_eps						//...param(0);
#define SHAPE_spheroid_A					1. //...большая полуось сфероида;
#undef  SHAPE_spheroid_A						//...param(1);
#endif
