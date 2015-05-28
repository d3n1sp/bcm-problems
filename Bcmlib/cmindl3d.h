/*===========================================*/
/*                  CMINDL3D                 */
/*===========================================*/
#ifndef ___CMindl3D___
#define ___CMindl3D___

#include "ccomput3d.h"

////////////////////////////////////////////////////////
//...class of blocks partition for double plane problem;
class CMindl3D : public CComput3D<double> {
public:
		Num_Draft type	  () { return MINDL3D_DRAFT;}
		int size_of_param() { return(15);}
//...constructor;
		CMindl3D (int num_phase = 7) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
		};
protected:
static int NUM_SHEAR, NUM_SHIFT, NUM_ADHES, MAX_PHASE, NUM_HESS;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_common_x(double * P, int i, int m);
		void jump1_common_y(double * P, int i, int m);
		void jump1_common_z(double * P, int i, int m);
		void jump2_common_x(double * P, int i, int m);
		void jump2_common_y(double * P, int i, int m);
		void jump2_common_z(double * P, int i, int m);
		void jump3_common_x(double * P, int i, int m);
		void jump3_common_y(double * P, int i, int m);
		void jump3_common_z(double * P, int i, int m);
		void jump4_common_x(double * P, int i, int m);
		void jump4_common_y(double * P, int i, int m);
		void jump4_common_z(double * P, int i, int m);
		void hessian_deriv_N (int k, double * P, int i);
public:
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B);
//...аналитические модели;
		double TakeEshelby_volm_two (double ff);
		double TakeEshelby_shear_two(double ff, double eps = EE_fine, int max_iter = 200);
		double TakeEshelby_shear_det(double ff);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_elem                10   //...number of element in facet subdivision;
#undef  DRAFT_N_elem                     //...param(2);
#define DRAFT_Q_facet               1.   //...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet                    //...param(3);
#define DRAFT_A                     0.   //...boundary jump coefficient in matrix;
#undef  DRAFT_A                          //...param(4);
#define DRAFT_B 							0    //...shear adhesion parameter in matrix;
#undef  DRAFT_B                          //...param(5);
#define DRAFT_G1                    1.   //...normalized shear modulus in matrix;
#undef  DRAFT_G1                         //...param(6);
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(7);
#define DRAFT_kappa_kk1					0.   //...wave number in matrix for compression wave;
#undef  DRAFT_kappa_kk1						  //...param(8);
#define DRAFT_kappa_mu1             0.   //...wave number in matrix for shear wave;
#undef  DRAFT_kappa_mu1                  //...param(9);
#define DRAFT_G2                    1.   //...normalized shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(10);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(11);
#define DRAFT_kappa_kk2					0.   //...wave number in inclusion for compression wave;
#undef  DRAFT_kappa_kk2						  //...param(12);
#define DRAFT_kappa_mu2             0.   //...wave number in inclusion for shear wave;
#undef  DRAFT_kappa_mu2                  //...param(13);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(14);
#endif
