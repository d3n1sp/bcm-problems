/*===========================================*/
/*                  CMINDL2D                 */
/*===========================================*/
#ifndef ___CMindl2D___
#define ___CMindl2D___

#include "ccomput2d.h"

////////////////////////////////////////////////////////
//...class of blocks partition for double plane problem;
class CMindl2D : public CComput2D<double> {
public:
		Num_Draft type	  () { return MINDL2D_DRAFT;}
		int size_of_param() { return(14);}
//...constructor;
		CMindl2D (int num_phase = 7) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
		};
protected:
static int NUM_SHEAR, NUM_SHIFT, NUM_ADHES, MAX_PHASE, NUM_HESS, regul;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_common_x	(double * P, int i, int m);
		void jump1_compress_x(double * P, int i, int m);
		void jump1_common_y	(double * P, int i, int m);
		void jump1_compress_y(double * P, int i, int m);
		void jump2_common_x	(double * P, int i, int m);
		void jump2_compress_x(double * P, int i, int m);
		void jump2_common_y	(double * P, int i, int m);
		void jump2_compress_y(double * P, int i, int m);
		void jump3_common_x	(double * P, int i, int m);
		void jump3_compress_x(double * P, int i, int m);
		void jump3_common_y	(double * P, int i, int m);
		void jump3_compress_y(double * P, int i, int m);
		void jump3_common_z	(double * P, int i, int m);
		void jump4_common_x	(double * P, int i, int m);
		void jump4_compress_x(double * P, int i, int m);
		void jump4_common_y	(double * P, int i, int m);
		void jump4_compress_y(double * P, int i, int m);
		void hessian_deriv_N (int k, double * P, int i);
		void jump_admittance (int l, int i, int m, double adm_re, int k = -1, double adm_im = 0.);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local){ return gram4(nd, i, id_local);};
		Num_State gram2    (CGrid * nd, int i, int id_local){ return gram3(nd, i, id_local);};
		Num_State gram3    (CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local){ return transfer4(nd, i, k, id_local);};
		Num_State transfer2(CGrid * nd, int i, int k, int id_local){ return transfer3(nd, i, k, id_local);};
		Num_State transfer3(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State rigidy1  (CGrid * nd, int i, double * K);
		Num_State computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B);
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void GetEnergyValue	(int k, double * energy);
//...аналитические модели для дисперсных включений;
		double TakeEshelby_volm_two (double ff);
		double TakeEshelby_shear_two(double ff, double eps = EE, int max_iter = 100);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_A                     0.   //...boundary jump coefficient in matrix;
#undef  DRAFT_A                          //...param(3);
#define DRAFT_B 							0    //...shear adhesion parameter in matrix;
#undef  DRAFT_B                          //...param(4);
#define DRAFT_G1                    1.   //...normalized shear modulus in matrix;
#undef  DRAFT_G1                         //...param(5);
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(6);
#define DRAFT_kappa_kk1					0.   //...wave number in matrix for compression wave;
#undef  DRAFT_kappa_kk1						  //...param(7);
#define DRAFT_kappa_mu1             0.   //...wave number in matrix for shear wave;
#undef  DRAFT_kappa_mu1                  //...param(8);
#define DRAFT_G2                    1.   //...normalized shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(9);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(10);
#define DRAFT_kappa_kk2					0.   //...wave number in inclusion for compression wave;
#undef  DRAFT_kappa_kk2						  //...param(11);
#define DRAFT_kappa_mu2             0.   //...wave number in inclusion for shear wave;
#undef  DRAFT_kappa_mu2                  //...param(12);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(13);
#endif
