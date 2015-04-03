/*==========================================*/
/*                  CHeat2D                 */
/*==========================================*/
#ifndef ___CHeat2D___
#define ___CHeat2D___

#include "ccomput2d.h"

//////////////////////////////////////////////////////////////
//...class of blocks partition for space electrstatic problem;
class CHeat2D : public CComput2D<double> {
public:
		Num_Draft type   () { return HEAT2D_DRAFT;}
		int size_of_param() { return(9);}
//...constructor;
		CHeat2D (int num_phase = 8) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
		};
protected:
static int NUM_BASIC, NUM_SHIFT, NUM_GEOMT, MAX_PHASE;
		 int block_shape_init(Block<double> & B, int id_free);
//...auxilliary operations with block matrix;
		void jump1(double * P, int i, int m);
		void jump2(double * P, int i, int m);
		void jump2_compos(double * P, int i, int m);
//...forming block matrix elements;
		int  gram1    (CGrid * nd, int i, int id_local);
		int  gram2    (CGrid * nd, int i, int id_local);
		int  gram2peri(CGrid * nd, int i, int id_local);
		int  gram3    (CGrid * nd, int i, int id_local);
		int  gram4    (CGrid * nd, int i, int id_local);
		int  transfer1(CGrid * nd, int i, int k, int id_local);
		int  transfer2(CGrid * nd, int i, int k, int id_local);
		int  transfer3(CGrid * nd, int i, int k, int id_local) { return transfer4(nd, i, k, id_local);};
		int  transfer4(CGrid * nd, int i, int k, int id_local);
		int  rigidy1  (CGrid * nd, int i, double * K);
		int  rigidy2  (CGrid * nd, int i, double * K);
		int  rigidy5  (CGrid * nd, int i, double * K);
		int  computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double K1, double K2, double K3 = 0.) { set_fasa_hmg(0., 0., K1, K2, K3);}
		void set_fasa_hmg(double R1, double R2, double K3, double K1, double K2);
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
//...аналитические модели;
		double TakeEshelby	 (double ff, double ff_l);
		double TakeEshelby_two(double ff);
		double TakeEshelby(int nn, double * RR, double * KK);
//...функция пересчета поля температур в слоистой модели (схема установления);
		void TakeStabStep(double * T, int NN, double alpha);
		void TakeStabStep_layer(double * T, int N_SC, int N_CU, int N_cells, double alpha);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_R1							1.	  //...geometry of inclusion;
#undef  DRAFT_R1								  //...param(3);
#define DRAFT_R2							1.	  //...geometry of intermediate layer;
#undef  DRAFT_R2								  //...param(4);
#define DRAFT_K3							1.   //...heat conduction in matrix;
#undef  DRAFT_K3								  //...param(5);
#define DRAFT_K1							1.   //...heat conduction in inclusion;
#undef  DRAFT_K1								  //...param(6);
#define DRAFT_K2							1.   //...heat conduction in intermediate;
#undef  DRAFT_K2								  //...param(7);
#define DRAFT_lagrange              0.   //...Lagrange coefficient for LSM in block functional;
#undef  DRAFT_lagrange                   //...param(8);
#endif
