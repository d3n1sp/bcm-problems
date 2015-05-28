/*===========================================*/
/*                  CCOHES3D                 */
/*===========================================*/
#ifndef ___CCohes3D___
#define ___CCohes3D___

#include "clame3d.h"
#define real_T double

////////////////////////////////////////////////////////
//...class of blocks partition for double plane problem;
class CCohes3D : public CLame3D {
public:
		Num_Draft type() { return COHES3D_DRAFT;}
//...constructor;
		CCohes3D (int num_phase = 7) {
			NUM_PHASE = num_phase;
			MAX_PHASE = 2;
		};
protected:
		int  block_shape_init(Block<double> & B, Num_State id_free);
public:
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2);
		void set_adhesion(double AA, double BB) { set_param(NUM_ADHES, AA); set_param(NUM_ADHES+1, BB);}
//...одномерные аналитические модели с межфазным слоем и адгезией;
		double TakeLayer_kk(int N, double * ff, double * kk, double * ll);
//...аналитические модели;
		double TakeEshelby_volm_two (double ff);
		double TakeEshelby_volm_sym (double ff);
		double TakeEshelby_shear_two(double ff, double eps = EE, int max_iter = 100);
		double TakeEshelby_shear_sym(double ff, double eps = EE, int max_iter = 100);
		double TakeEshelby_shear	 (double ff, double nju1, double nju2, double E1, double E2, double l1, double l2);
		double TakeEshelby_shear_old(double ff, double nju1, double nju2, double E1, double E2, double l1, double l2);
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
#define DRAFT_AA                    0.   //...boundary jump coefficient in matrix;
#undef  DRAFT_AA                         //...param(4);
#define DRAFT_BB 							0    //...shear adhesion parameter in matrix;
#undef  DRAFT_BB                         //...param(5);
#define DRAFT_reserve_alpha			0.   //...зарезервированный параметр;
#undef  DRAFT_reserve_alpha              //...param(6);
#define DRAFT_C1                    0.   //...cohesion parameter in matrix;
#undef  DRAFT_C1                         //...param(7);
#define DRAFT_G1                    1.   //...normalized shear modulus in matrix;
#undef  DRAFT_G1                         //...param(8);
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(9);
#define DRAFT_alpha1                0.   //...coefficient of Neuber-Papkovich representation in matrix;
#undef  DRAFT_alpha1                     //...param(10);
#define DRAFT_kappa_2mu_lm          0.   //...wave number for compression wave in matrix;
#undef  DRAFT_kappa_2mu_lm               //...param(11);
#define DRAFT_kappa_mu              0.   //...wave number for shear wave in matrix;
#undef  DRAFT_kappa_mu                   //...param(12);
#define DRAFT_C2                    0.   //...cohesion parameter in inclusion;
#undef  DRAFT_C2                         //...param(13);
#define DRAFT_G2                    1.   //...normalized shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(14);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(15);
#define DRAFT_alpha2                0.   //...coefficient of Neuber-Papkovich representation in inclusion;
#undef  DRAFT_alpha2                     //...param(16);
#define DRAFT_kappa_2mu_lm          0.   //...wave number for compression wave in inclusion;
#undef  DRAFT_kappa_2mu_lm               //...param(17);
#define DRAFT_kappa_mu              0.   //...wave number for shear wave in inclusion;
#undef  DRAFT_kappa_mu                   //...param(18);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(19);
#endif
