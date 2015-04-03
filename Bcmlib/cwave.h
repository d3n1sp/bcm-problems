/*========================================*/
/*                 CWAVE                  */
/*========================================*/
#ifndef ___CWAVE___
#define ___CWAVE___

#include "cshapes.h"

////////////////////////////////////////////////////////
//   WAVE MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//...class of 2D wave-trasfer multipoles with polynomial behaviour;
class CWave2DPoly : public CShape<double> {
public:
		int type() {return WV2D_POLY_SHAPE;}
		int freedom(int  m) { return m*2+1; }
		int size_of_param() { return(3);};
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0., double r = 1., double L1 = 0., double L2 = 0., double L3 = 0.);
		void init2 (int N, int M, int dim);
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_t(double * deriv, double f = 1.);
//...constructor;
		CWave2DPoly () {
			delete_struct(param);
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
		};
};

///////////////////////////////////////////////////////////////////
//...class of 3D wave-trasfer multipoles with polynomial behaviour;
class CWave3DPoly : public CWave2DPoly {
public:
		int type()				{ return WV3D_POLY_SHAPE;}
		int freedom (int  m) { return (m+1)*(m+1); }
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0);
		void parametrization_hess(double * P = NULL, int m = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
		void deriv_t(double * deriv, double f = 1.);
//...constructor;
		CWave3DPoly () {
			delete_struct(param);
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_C                 0. //...normalized heat-transfer coefficient;
#undef  SHAPE_norm_C                    //...param(0);
#define SHAPE_norm_C_inv             1. //...inverse of normalized heat-transfer coefficient;
#undef  SHAPE_norm_C_inv                //...param(1);
#define SHAPE_heat_C                 0. //...non-normalized heat-transfer coefficient;
#undef  SHAPE_heat_C                    //...param(2);

#endif
