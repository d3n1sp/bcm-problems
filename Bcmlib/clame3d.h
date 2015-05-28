/*==========================================*/
/*                  CLAME3D                 */
/*==========================================*/
#ifndef ___CLame3D___
#define ___CLame3D___

#include "ccomput3d.h"
#define ___3BODY_ESHELBY___

////////////////////////////////////////////////////////
//...class of blocks partition for double plane problem;
class CLame3D : public CComput3D<double> {
public:
		Num_Draft type   () { return LAME3D_DRAFT;}
		int size_of_param() { return(26);}
//...constructor;
		CLame3D (int num_phase = 15) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
			TT = TH = NULL;
			C0 = C1 = C2 = A1 = B1 = A2 = B2 = NULL;
		}
      virtual ~CLame3D (void);
protected:
static int NUM_SHEAR, NUM_SHIFT, NUM_HESS;
		double *** TT,*** TH,  * C0, * C1, * C2, * A1, * B1, * A2, * B2;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_classic_x(double * P, int i, int m);
		void jump1_classic_y(double * P, int i, int m);
		void jump1_classic_z(double * P, int i, int m);
		void jump1_sphere_x_2body(double * P, int i, int m);
		void jump1_sphere_y_2body(double * P, int i, int m);
		void jump1_sphere_z_2body(double * P, int i, int m);
		void jump1_sphere_x_3body(double * P, int i, int m);
		void jump1_sphere_y_3body(double * P, int i, int m);
		void jump1_sphere_z_3body(double * P, int i, int m);
		void jump1_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump1_sphere_x_2body(P, i, m); else jump1_classic_x(P, i, m);
#else
				jump1_sphere_x_3body(P, i, m); else jump1_classic_x(P, i, m);
#endif
		}
		void jump1_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump1_sphere_y_2body(P, i, m); else jump1_classic_y(P, i, m);
#else
				jump1_sphere_y_3body(P, i, m); else jump1_classic_y(P, i, m);
#endif
		}
		void jump1_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump1_sphere_z_2body(P, i, m); else jump1_classic_z(P, i, m);
#else
				jump1_sphere_z_3body(P, i, m); else jump1_classic_z(P, i, m);
#endif
		}
		void jump2_classic_x(double * P, int i, int m);
		void jump2_classic_y(double * P, int i, int m);
		void jump2_classic_z(double * P, int i, int m);
		void jump2_sphere_x_2body(double * P, int i, int m);
		void jump2_sphere_y_2body(double * P, int i, int m);
		void jump2_sphere_z_2body(double * P, int i, int m);
		void jump2_sphere_x_3body(double * P, int i, int m);
		void jump2_sphere_y_3body(double * P, int i, int m);
		void jump2_sphere_z_3body(double * P, int i, int m);
		void jump2_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump2_sphere_x_2body(P, i, m); else jump2_classic_x(P, i, m);
#else
				jump2_sphere_x_3body(P, i, m); else jump2_classic_x(P, i, m);
#endif
		}
		void jump2_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump2_sphere_y_2body(P, i, m); else jump2_classic_y(P, i, m);
#else
				jump2_sphere_y_3body(P, i, m); else jump2_classic_y(P, i, m);
#endif
		}
		void jump2_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump2_sphere_z_2body(P, i, m); else jump2_classic_z(P, i, m);
#else
				jump2_sphere_z_3body(P, i, m); else jump2_classic_z(P, i, m);
#endif
		}
		void jump4_classic_x(double * P, int i, int m);
		void jump4_classic_y(double * P, int i, int m);
		void jump4_classic_z(double * P, int i, int m);
		void jump4_sphere_x_2body(double * P, int i, int m);
		void jump4_sphere_y_2body(double * P, int i, int m);
		void jump4_sphere_z_2body(double * P, int i, int m);
		void jump4_sphere_x_3body(double * P, int i, int m);
		void jump4_sphere_y_3body(double * P, int i, int m);
		void jump4_sphere_z_3body(double * P, int i, int m);
		void jump4_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump4_sphere_x_2body(P, i, m); else jump4_classic_x(P, i, m);
#else
				jump4_sphere_x_3body(P, i, m); else jump4_classic_x(P, i, m);
#endif
		}
		void jump4_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump4_sphere_y_2body(P, i, m); else jump4_classic_y(P, i, m);
#else
				jump4_sphere_y_3body(P, i, m); else jump4_classic_y(P, i, m);
#endif
		}
		void jump4_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
#ifdef ___2BODY_ESHELBY___
				jump4_sphere_z_2body(P, i, m); else jump4_classic_z(P, i, m);
#else
				jump4_sphere_z_3body(P, i, m); else jump4_classic_z(P, i, m);
#endif
		}
		void hessian_deriv_N (double * P, int i);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1	  (CGrid * nd, int i, int id_local);
		Num_State gram2	  (CGrid * nd, int i, int id_local);
		Num_State gram2peri(CGrid * nd, int i, int id_local);
		Num_State gram3	  (CGrid * nd, int i, int id_local);
		Num_State gram4	  (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer2(CGrid * nd, int i, int k, int id_local);
		Num_State trans_esh(CGrid * nd, int i, int k, int id_local);
		Num_State transfer3(CGrid * nd, int i, int k, int id_local) { return transfer4(nd, i, k, id_local);};
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State rigidy1  (CGrid * nd, int i, double * K);
		Num_State rigidy2  (CGrid * nd, int i, double * K);
		Num_State rigidy5  (CGrid * nd, int i, double * K);
		Num_State computing_header(Num_Comput Num);
public:
//...collocation rigidy matrix;
		int  rigidy1_collocat(CGrid * nd, int i, double * K);
//...analytical method for junction on spherical inclusion;
		void rgradf_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.);
		void rdivrf_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.);
		void unit_f_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.);
		void toreal_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0);
		void grdivf_transf(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.);
		void toreal_transf(double ** T, int n, int shift_m = 0, int shift_n = 0);
//...algorithm of junction on spherical inclusion;
		void set_eshelby_matrix(int N_mpl);
		void		eshelby_matrix_3body(double **& TT, double **& TX, int n, double **& TD, double **& TH, int num_bound);
		void		eshelby_matrix_2body(double **& TT, double **& TX, int n, double **& TD, double **& TH, int num_bound);
		void add_collocation(int i, int m, int j, int n, double * ff = NULL, double * ff_reg = NULL);
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double nju3, double G3);
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2) { set_fasa_hmg(nju1, nju2, G1, G2, 0.3, 1.);}
		void set_fasa_hmg(double R1, double R2, double nju1, double nju2, double nju3, double G1, double G2, double G3, double alpha = -1.);
		void set_geometry(double rad, double layer = 0.) { set_param(NUM_GEOMT, rad); set_param(NUM_GEOMT+1, rad+layer);}
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void block_descrap(char * OUT_FILE);
		void GetEnergy(double * energy, int id_variant = AXIS_Z, Num_Comput Num = BASIC_COMPUT);
		void GetEnergyValue(int k, double * energy);
//...одномерные аналитические модели;
		double TakeLayer_kk(int N, double * ff, double * kk);
//...одномерные аналитические модели с межфазным слоем в сферической симметрии;
		double TakeLayer_kk(int N, double * ff, double * kv, double * mu);
		double TakeLayer_GG(int N, double * ff, double * kv, double * mu, double * nj, double eps = EE, int max_iter = 100);
//...аналитические модели (самосогласованные);
		double TakeEshelby_volm (double ff, double ff_l);
		double TakeEshelby_shear(double ff, double ff_l, double eps = 1e-12, int max_iter = 100);
		double TakeEshelby_volm_two (double ff);
		double TakeEshelby_shear_two(double ff, double eps = EE, int max_iter = 100);
		double TakeEshelby_shear_sys(double ff, double eps = EE, int max_iter = 100);
		double TakeEshelby_shear_det(double ff, double alpha = -1);
		double TakeEshelby_shear    (double ff, double nju1, double nju2, double E1, double E2);
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
#define DRAFT_R1							1.	  //...geometry of inclusion;
#undef  DRAFT_R1								  //...param(4);
#define DRAFT_R2							1.	  //...geometry of intermediate layer;
#undef  DRAFT_R2								  //...param(5);
#define DRAFT_alpha					  -1.	  //...friction stiffnes (III type boundary condition);
#undef  DRAFT_alpha							  //...param(6);
#define DRAFT_reserveC1				   0.   //...зарезервированный параметр;
#undef  DRAFT_reserveC1                  //...param(7);
#define DRAFT_G1                    1.   //...shear modulus in matrix;
#undef  DRAFT_G1                         //...param(8); 
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(9);
#define DRAFT_alpha1                0.   //...coefficient of Neuber-Papkovich representation in matrix;
#undef  DRAFT_alpha1                     //...param(10);
#define DRAFT_reserveA1				   0.   //...зарезервированный параметр;
#undef  DRAFT_reserveA1						  //...param(11);
#define DRAFT_reserveB1             0.   //...зарезервированный параметр;
#undef  DRAFT_reserveB1                  //...param(12);
#define DRAFT_reserveC2             0.   //...зарезервированный параметр;
#undef  DRAFT_reserveC2                  //...param(13);
#define DRAFT_G2                    1.   //...shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(14);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(15);
#define DRAFT_alpha2                0.   //...coefficient of Neuber-Papkovich representation in matrix;
#undef  DRAFT_alpha2                     //...param(16);
#define DRAFT_reserveA2				   0.   //...зарезервированный параметр;
#undef  DRAFT_reserveA2						  //...param(17);
#define DRAFT_reserveB2             0.   //...зарезервированный параметр;
#undef  DRAFT_reserveB2                  //...param(18);
#define DRAFT_reserveC3             0.   //...зарезервированный параметр;
#undef  DRAFT_reserveC3                  //...param(19);
#define DRAFT_G3                    1.   //...normalized shear modulus in interphase layer;
#undef  DRAFT_G3                         //...param(20);
#define DRAFT_nju3                  0.3  //...Poisson coefficient in interphase layer;
#undef  DRAFT_nju3                       //...param(21);
#define DRAFT_alpha3                0.   //...coefficient of Neuber-Papkovich representation in matrix;
#undef  DRAFT_alpha3                     //...param(22);
#define DRAFT_reserveA3				 0.   //...зарезервированный параметр;
#undef  DRAFT_reserveA3						//...param(23);
#define DRAFT_reserveB3             0.   //...зарезервированный параметр;
#undef  DRAFT_reserveB3                  //...param(24);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(25);
#endif
