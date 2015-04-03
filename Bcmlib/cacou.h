/*=========================================*/
/*                 CACOU                   */
/*=========================================*/
#ifndef ___CACOU___
#define ___CACOU___

#include "cshapes.h"

//////////////////////////////////////////////
//														  //
//  BASIC SYSTEM of ELLIPSOIDAL MULTIPOLES  //
//														  //
//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles, based on plane wave in Z-direction;
class CAcou3DEll : public CShape<double> {
public:
		int type()          { return AU3D_ELLI_SHAPE;}
		int freedom (int m) { return (m+1)*(m+2);}
		int size_of_param() { return(5);};
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CAcou3DEll () {
			delete_struct(param);
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_kappa_shift               0. //...shift-normalized wave number;
#undef  SHAPE_kappa_shift                  //...param(4);


///////////////////////////////////////////////////////////
//																			//
//  ACOUSTIC MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
//																			//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles with polynomial behaviour along the line;
class CAcou3DPoly : public CShape<double> {
protected:
		double * au, * E1;
public:
		int type()          { return AU3D_POLY_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);}
		void release();
public:
//...calculation of multipoles;
		void parametrization(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CAcou3DPoly () {
			delete_struct(param);
			param   = (Param *)new_struct(size_of_param()*sizeof(Param));
			au = E1 = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);


////////////////////////////////////////////////////////////
//																			 //
//  ACOUSTIC MULTIPOLES for SPHERE INTERIOUR OR EXTERIOR  //
//																			 //
////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//...class of classic acoustic multipoles for interior or exterior;
class CAcou3DZoom : public CShape<double> {
protected:
		double * au;
public:
		int type()          { return AU3D_ZOOM_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);}
		void release();
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CAcou3DZoom (int inverse = 0) {
			delete_struct(param);
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			au = NULL;
			id_inverse = inverse;
//       id_cmpl	  = 1; //...we can install id_cmpl by void change()!!!
       };
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);


//////////////////////////////////////////
//													 //
//     SUPERPOSITION of PLANE WAVES     //
//													 //
//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//...class of exponential plane waves (separated as complex potential with real and imaginary part);
class CAcou3DWave : public CShape<double> {
protected:
       double * cz, * sz, * cy, * sy;
public:
       int type()          { return AU3D_WAVE_SHAPE;}
		 int freedom (int m) { return sqr(m+1);}
       int size_of_param() { return(4);}
public:
//...initialization and calculation of multipoles;
       void init1(int N, int dim) { init2(2, N, dim);};
       void init2(int id_flag, int N, int dim);
       void parametrization(double * P, int m_dop = 1);
       int  add_expo_power (double X, double Y, double Z, int flag = 0);
//...differentiation;
       void deriv_X(double * deriv, double f = 1.);
       void deriv_Y(double * deriv, double f = 1.);
       void deriv_Z(double * deriv, double f = 1.);
       void deriv_N();
//...constructor;
       CAcou3DWave () {
           delete_struct(param);
           param   = (Param *)new_struct(size_of_param()*sizeof(Param));
           cz = sz = cy = sy = NULL;
           //id_cmpl = 1;
       };
//...destructor;
      ~CAcou3DWave (void);
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);


////////////////////////////////////////////////////////////
//																			 //
//  BASIC SYSTEM of ACOUSTIC MULTIPOLES (acoustic beam)   //
//																			 //
////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles of beams type with two oscillation characteristics;
class CAcou3DBeam : public CShape<double> {
protected:
		double  * au, * az, * pim, * pxim, * pyim, * pzim;
		complex * E1;
public:
		int type()          { return AU3D_BEAM_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(6);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0., double kk_dop = 0., double L1 = 0., double L2 = 0.);
		void release();
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CAcou3DBeam () {
			delete_struct(param);
			param   = (Param *)new_struct(size_of_param()*sizeof(Param));
			au = az = pim = pxim = pyim = pzim = NULL;
			E1 = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized cylindrical wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_norm_kappa_dop            0. //...normalized directional wave number;
#undef  SHAPE_norm_kappa_dop               //...param(4);
#define SHAPE_kappa_dop                 0. //...non-normalized wave number;
#undef  SHAPE_kappa_dop                    //...param(5);

/////////////////////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles of beams type with oscillation characteristics along Z;
class CAcou3DBeamZ : public CShape<double> {
protected:
		double  * az, * pim, * pxim, * pyim, * pzim;
		complex * E1;
public:
		int type()          { return AU3D_BEAMZSHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);};
		void release();
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
//...constructor;
		CAcou3DBeamZ () {
			delete_struct(param);
			param   = (Param *)new_struct(size_of_param()*sizeof(Param));
			az = pim = pxim = pyim = pzim = NULL;
			E1 = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number (directional);
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#endif
