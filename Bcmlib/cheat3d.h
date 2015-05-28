/*==========================================*/
/*                  CHEAT3D                 */
/*==========================================*/
#ifndef ___CHeat3D___
#define ___CHeat3D___

#include "ccomput3d.h"

//////////////////////////////////////////////////////////////
//...class of blocks partition for space electrstatic problem;
class CHeat3D : public CComput3D<double> {
public:
		Num_Draft type   () { return HEAT3D_DRAFT;}
		int size_of_param() { return(11);}
//...constructor;
		CHeat3D (int num_phase = 8) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
		}
protected:
static int NUM_BASIC, NUM_SHIFT, NUM_GEOMT, MAX_PHASE;
		 int block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1(double * P, int i, int m);
		void jump2(double * P, int i, int m);
		void jump3(double * P, int i, int m);
		void jump4(double * P, int i, int m);
		void jump4_compos (double * P, int i, int m);
		void jump1_classic(double * P, int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local);
		Num_State gram2    (CGrid * nd, int i, int id_local);
		Num_State gram2_old(CGrid * nd, int i, int id_local);
		Num_State gram2peri(CGrid * nd, int i, int id_local);
		Num_State gram3    (CGrid * nd, int i, int id_local);
		Num_State gram3_old(CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer2(CGrid * nd, int i, int k, int id_local);
		Num_State trans_esh(CGrid * nd, int i, int k, int id_local);
		Num_State transfer3(CGrid * nd, int i, int k, int id_local) { return transfer4(nd, i, k, id_local);}
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State rigidy1  (CGrid * nd, int i, double * K);
		Num_State rigidy2  (CGrid * nd, int i, double * K);
		Num_State rigidy5  (CGrid * nd, int i, double * K);
		Num_State computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double K1, double K2, double K3) { set_fasa_hmg (0., 0., K1, K2, K3, 0.);}
		void set_fasa_hmg(double K1, double K2, double K3, double C) { set_fasa_hmg (0., 0., K1, K2, K3, C);}
		void set_fasa_hmg(double R1, double R2, double K3, double K1, double K2) { set_fasa_hmg (R1, R2, K3, K1, K2, 0.);}
		void set_fasa_hmg(double R1, double R2, double K3, double K1, double K2, double C);
		void set_geometry(double rad, double layer = 0.) { set_param(NUM_GEOMT, rad); set_param(NUM_GEOMT+1, rad+layer);}
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
//...аналитические модели;
		double TakeEshelby(double ff, double ff_l);
		double TakeEshelby_two (double ff);
		double TakeEshelby_grad(double ff);
		double TakeEllipsoidEshelby(double ff, double eps, double phi_stream = 0., int NX = 20, int NY = 20);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_Q_facet               0.   //...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet                    //...param(3);
#define DRAFT_C							0.   //...параметр когезионного поля;
#undef  DRAFT_C									//...param(4);
#define DRAFT_R1							1.	//...geometry of inclusion;
#undef  DRAFT_R1									//...param(5);
#define DRAFT_R2							1.	//...geometry of intermediate layer;
#undef  DRAFT_R2									//...param(6);
#define DRAFT_K3							1.   //...heat conduction in matrix;
#undef  DRAFT_K3									//...param(7);
#define DRAFT_K1							1.   //...heat conduction in inclusion;
#undef  DRAFT_K1									//...param(8);
#define DRAFT_K2							1.   //...heat conduction in intermediate;
#undef  DRAFT_K2									//...param(9);
#define DRAFT_lagrange              0.   //...Lagrange coefficient for LSM in block functional;
#undef  DRAFT_lagrange                   //...param(10);
#endif
