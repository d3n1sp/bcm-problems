/*==========================================*/
/*                  CLAME2D                 */
/*==========================================*/
#ifndef ___CLame2D___
#define ___CLame2D___

#include "ccomput2d.h"

/////////////////////////////////////////////////
//...class of blocks partition for plane problem;
class CLame2D : public CComput2D<double> {
public:
		Num_Draft type   () { return LAME2D_DRAFT;}
		int size_of_param() { return(11);}
//...constructor;
		CLame2D(int num_phase = 8) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
			NUM_GEOMT = 3;
		}
protected:
static int NUM_SHEAR, NUM_SHIFT;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_classic_x (double * P, int i, int m);
		void jump1_classic_y (double * P, int i, int m);
		void jump2_classic_x (double * P, int i, int m);
		void jump2_classic_y (double * P, int i, int m);
		void jump3_classic_x (double * P, int i, int m);
		void jump3_classic_y (double * P, int i, int m);
		void jump4_classic_x (double * P, int i, int m);
		void jump4_classic_y (double * P, int i, int m);
		void jump5_classic_x (double * P, int i, int m);
		void jump5_classic_y (double * P, int i, int m);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local);
		Num_State gram2    (CGrid * nd, int i, int id_local);
		Num_State gram3    (CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer2(CGrid * nd, int i, int k, int id_local);
		Num_State transfer3(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4(CGrid * nd, int i, int k, int id_local){ return transfer3(nd, i, k, id_local);}
		Num_State rigidy1  (CGrid * nd, int i, double * K);
		Num_State computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2) { set_fasa_hmg (nju1, nju2, nju2, G1, G2, G2);}
		void set_fasa_hmg(double nju1, double nju2, double nju3, double G1, double G2, double G3);
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void GetEnergy			(double * energy, Num_Value _FMF = ENERGY_VALUE);
		void GetEnergyValue(int k, double * energy);
//...одномерные аналитические модели с межфазным слоем в осевой симметрии;
		double TakeLayer_kk(int N, double * ff, double * kk);
		double TakeLayer_kk(int N, double * ff, double * kp, double * mu);
		double TakeLayer_kk(int N, double * ff, double * lm, double * mu, double * nj);
		double TakeLayer_GG(int N, double * ff, double * mu, double * nj);
		double TakeLayer_GG(int N, double * ff, double * kp, double * mu, double * nj, double eps = EE, int max_iter = 100);
		double TakeLayer_sh(int N, double * ff, double * kp, double * mu, double * nj, double eps = EE, int max_iter = 100);
//...аналитические модели;
		double TakeLayer_E1(double ff);
		double TakeLayer_E2(double ff);
		double TakeLayer_G1(double ff);
		double TakeLayer_G2(double ff);
		double TakeEshelby_volm	   (double ff, double ff_l) {return 0.;};
		double TakeEshelby_volm_two(double ff) {return 0.;};
};

////////////////////////////////// 
//...parametrization of the draft;
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_requl						1.   //...normalization coefficient in reqularization;
#undef  DRAFT_requl							  //...param(3);
#define DRAFT_G1                    1.   //...shear modulus in matrix;
#undef  DRAFT_G1                         //...param(4);
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(5);
#define DRAFT_G2                    1.   //...shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(6);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(7);
#define DRAFT_G3                    1.   //...shear modulus in inclusion;
#undef  DRAFT_G3                         //...param(8);
#define DRAFT_nju3                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju3                       //...param(9);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(10);
#endif
