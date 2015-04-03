/*==========================================*/
/*                  CHYDRO3D                */
/*==========================================*/
#ifndef ___CHydro3D___
#define ___CHydro3D___

#include "ccomput3d.h"
#define real_T double

/////////////////////////////////////////////////
//...class of blocks partition for Stokes problem;
class CHydro3D : public CComput3D<real_T> {
public:
		Num_Draft type   () { return HYDRO3D_DRAFT;}
		int size_of_param() { return(10);}
//...constructor and destructor;
		CHydro3D (int num_phase = 8) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			param[NUM_MPLS] = 
			param[NUM_QUAD] = PackInts(2, 2);
			param[NUM_MPLS+1]  = 1.;
			param[NUM_GEOMT]	 = 1.;
			param[NUM_SHEAR]	 = 1.;
			param[NUM_SHEAR+2] = 1.;
			NUM_PHASE = num_phase;
			TT = TH = NULL; TN = NULL;
			BOX_LINK_PERIOD = OK_STATE;
		}
      virtual ~CHydro3D (void);
protected:
static int NUM_SHEAR, NUM_HESS, NUM_GEOMT, MAX_PHASE;
		real_T *** TT, *** TH, ** TN;
		int  block_shape_init(Block<real_T> & B, int id_free);
//...auxilliary operations with block matrix;
		void jump0_pressure(double * P, int i, int m);
		void jump0_sphere  (double * P, int i, int m);
		void jump0_current (double * P, int i, int m);
		void jump0_darcy	 (double * P, int i, int m);
		void jump0			 (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)){
				if (! get_param(NUM_SHEAR+1)) jump0_sphere(P, i, m); 
				else jump0_current(P, i, m);
			}
			else if ((B[i].type & ERR_CODE) == POLY_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump0_darcy	 (P, i, m);
			else jump0_pressure(P, i, m);
		}
		void jump1_fluid_x  (double * P, int i, int m); 
		void jump1_fluid_y  (double * P, int i, int m);
		void jump1_fluid_z  (double * P, int i, int m);
		void jump1_sphere_x (double * P, int i, int m);
		void jump1_sphere_y (double * P, int i, int m);
		void jump1_sphere_z (double * P, int i, int m);
		void jump1_brinkmn_x(double * P, int i, int m);
		void jump1_brinkmn_y(double * P, int i, int m);
		void jump1_brinkmn_z(double * P, int i, int m);
		void jump1_current_x(double * P, int i, int m);
		void jump1_current_y(double * P, int i, int m);
		void jump1_current_z(double * P, int i, int m);
		void jump1_darcy_x  (double * P, int i, int m);
		void jump1_darcy_y  (double * P, int i, int m);
		void jump1_darcy_z  (double * P, int i, int m);
		void jump1_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump1_sphere_x(P, i, m); 
				else jump1_current_x(P, i, m);
			}
			else if ((B[i].type & ERR_CODE) == POLY_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_x(P, i, m);
			else
			if (! get_param(NUM_SHEAR+1)) jump1_fluid_x(P, i, m);
			else								 jump1_brinkmn_x(P, i, m);
		}
		void jump1_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump1_sphere_y(P, i, m); 
				else jump1_current_y(P, i, m);
			}
			else if ((B[i].type & ERR_CODE) == POLY_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_y(P, i, m);
			else
			if (! get_param(NUM_SHEAR+1)) jump1_fluid_y(P, i, m);
			else								 jump1_brinkmn_y(P, i, m);
		}
		void jump1_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump1_sphere_z(P, i, m); 
				else jump1_current_z(P, i, m);
			}
			else if ((B[i].type & ERR_CODE) == POLY_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_z(P, i, m);
			else
			if (! get_param(NUM_SHEAR+1)) jump1_fluid_z(P, i, m);
			else								 jump1_brinkmn_z(P, i, m);
		}
		void jump2_fluid_x  (double * P, int i, int m);
		void jump2_fluid_y  (double * P, int i, int m);
		void jump2_fluid_z  (double * P, int i, int m);
		void jump2_sphere_x (double * P, int i, int m);
		void jump2_sphere_y (double * P, int i, int m);
		void jump2_sphere_z (double * P, int i, int m);
		void jump2_brinkmn_x(double * P, int i, int m);
		void jump2_brinkmn_y(double * P, int i, int m);
		void jump2_brinkmn_z(double * P, int i, int m);
		void jump2_current_x(double * P, int i, int m);
		void jump2_current_y(double * P, int i, int m);
		void jump2_current_z(double * P, int i, int m);
		void jump2_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump2_sphere_x(P, i, m); 
				else jump2_current_x(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump2_fluid_x(P, i, m);
			else								 jump2_brinkmn_x(P, i, m);
		}
		void jump2_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump2_sphere_y(P, i, m); 
				else jump2_current_y(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump2_fluid_y(P, i, m);
			else								 jump2_brinkmn_y(P, i, m);
		}
		void jump2_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump2_sphere_z(P, i, m); 
				else jump2_current_z(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump2_fluid_z(P, i, m);
			else								 jump2_brinkmn_z(P, i, m);
		}
		void jump3_fluid	(double * P, int i, int m);
		void jump3_sphere	(double * P, int i, int m);
		void jump3_brinkmn(double * P, int i, int m);
		void jump3_current(double * P, int i, int m);
		void jump3			(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump3_sphere(P, i, m); 
				else jump3_current(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump3_fluid(P, i, m);
			else								 jump3_brinkmn(P, i, m);
		}
		void jump4_fluid_x  (double * P, int i, int m);
		void jump4_fluid_y  (double * P, int i, int m);
		void jump4_fluid_z  (double * P, int i, int m);
		void jump4_sphere_x (double * P, int i, int m){} //...для граничной функции Стокса реализация напряжений отсутствует;
		void jump4_sphere_y (double * P, int i, int m){}
		void jump4_sphere_z (double * P, int i, int m){}
		void jump4_brinkmn_x(double * P, int i, int m);
		void jump4_brinkmn_y(double * P, int i, int m);
		void jump4_brinkmn_z(double * P, int i, int m);
		void jump4_current_x(double * P, int i, int m){} //...для граничной функции Бринкмана реализация напряжений отсутствует;
		void jump4_current_y(double * P, int i, int m){}
		void jump4_current_z(double * P, int i, int m){}
		void jump4_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump4_sphere_x(P, i, m); 
				else jump4_current_x(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump4_fluid_x(P, i, m);
			else								 jump4_brinkmn_x(P, i, m);
		}
		void jump4_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump4_sphere_y(P, i, m); 
				else jump4_current_y(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump4_fluid_y(P, i, m);
			else								 jump4_brinkmn_y(P, i, m);
		}
		void jump4_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump4_sphere_z(P, i, m); 
				else jump4_current_z(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump4_fluid_z(P, i, m);
			else								 jump4_brinkmn_z(P, i, m);
		}
		void jump5_fluid_x  (double * P, int i, int m);
		void jump5_fluid_y  (double * P, int i, int m);
		void jump5_fluid_z  (double * P, int i, int m);
		void jump5_sphere_x (double * P, int i, int m);
		void jump5_sphere_y (double * P, int i, int m);
		void jump5_sphere_z (double * P, int i, int m);
		void jump5_brinkmn_x(double * P, int i, int m);
		void jump5_brinkmn_y(double * P, int i, int m);
		void jump5_brinkmn_z(double * P, int i, int m);
		void jump5_current_x(double * P, int i, int m);
		void jump5_current_y(double * P, int i, int m);
		void jump5_current_z(double * P, int i, int m);
		void jump5_x		  (double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump5_sphere_x(P, i, m); 
				else jump5_current_x(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump5_fluid_x(P, i, m);
			else								 jump5_brinkmn_x(P, i, m);
		}
		void jump5_y(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump5_sphere_y(P, i, m); 
				else jump5_current_y(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump5_fluid_y(P, i, m);
			else								 jump5_brinkmn_y(P, i, m);
		}
		void jump5_z(double * P, int i, int m) {
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! get_param(NUM_SHEAR+1)) jump5_sphere_z(P, i, m); 
				else jump5_current_z(P, i, m);
			}
			else
			if (! get_param(NUM_SHEAR+1)) jump5_fluid_z(P, i, m);
			else								 jump5_brinkmn_z(P, i, m);
		}
		void hessian_deriv_N (double * P, int i);
		void hessian_deriv_N (int k, double * P, int i);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		int  gram1_phase1(CGrid * nd, int i, int id_local);
		int  gram1_phase2(CGrid * nd, int i, int id_local);
		int  gram1		  (CGrid * nd, int i, int id_local) {
			if (B[i].link[NUM_PHASE] == -1) return gram1_phase1(nd, i, id_local); else
			if (B[i].link[NUM_PHASE] == -2) return gram1_phase2(nd, i, id_local); else return(ERR_STATE);
		}
		int  gram2	  (CGrid * nd, int i, int id_local);
		int  gram2peri(CGrid * nd, int i, int id_local);
		int  gram3	  (CGrid * nd, int i, int id_local);
		int  gram4_phase1(CGrid * nd, int i, int id_local);
		int  gram4_phase2(CGrid * nd, int i, int id_local);
		int  gram4		  (CGrid * nd, int i, int id_local) {
			if (B[i].link[NUM_PHASE] == -1) return gram4_phase1(nd, i, id_local); else
			if (B[i].link[NUM_PHASE] == -2) return gram4_phase2(nd, i, id_local); else return(ERR_STATE);
		}
		int  transfer1_phase1(CGrid * nd, int i, int k, int id_local);
		int  transfer1_phase2(CGrid * nd, int i, int k, int id_local);
		int  transfer1			(CGrid * nd, int i, int k, int id_local) {
			if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -1) return transfer1_phase1(nd, i, k, id_local); else
			if (B[i].link[NUM_PHASE] == -2 && B[k].link[NUM_PHASE] == -2) return transfer1_phase2(nd, i, k, id_local); else return(ERR_STATE);
		}
		int  transfer2			(CGrid * nd, int i, int k, int id_local);
	 //int  transfer3			(CGrid * nd, int i, int k, int id_local);
		int  transfer4_phase1(CGrid * nd, int i, int k, int id_local);
		int  transfer4_phase2(CGrid * nd, int i, int k, int id_local);
		int  transfer4			(CGrid * nd, int i, int k, int id_local) {
			if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -1) return transfer4_phase1(nd, i, k, id_local); else
			if (B[i].link[NUM_PHASE] == -2 && B[k].link[NUM_PHASE] == -2) return transfer4_phase2(nd, i, k, id_local); else return(ERR_STATE);
		}
		int  rigidy1  (CGrid * nd, int i, double * K);
		int  rigidy2  (CGrid * nd, int i, double * K);
		int  computing_header (Num_Comput Num);
		int  computing_kernel5();//...счетная схема для явного аналитического решения;
public:
//...analytical method for boundary functions;
		void add_collocation(int i, int m, int j);
		void intern_matrix(real_T **& T, int N, real_T **& TH);
		void extern_matrix(real_T **& T, int N);
		void set_eshelby_matrix (int N_mpl);
		void set_normaliz_vector(int N_mpl, double rad, double kk);
		void set_normaliz_vector_old(int N_mpl, double rad, double kk);
//...параметры задачи;
		void set_fasa_hmg(double R0, double ll, double G_dynamic, double alpha, double kk);
		void set_fasa_hmg(double G_dynamic, double alpha) { set_fasa_hmg(get_param(NUM_GEOMT), get_param(NUM_GEOMT+1), G_dynamic, alpha, get_param(NUM_SHEAR+2));}
		void set_fasa_hmg(double alpha = 0.) {set_fasa_hmg(get_param(NUM_GEOMT), get_param(NUM_GEOMT+1), get_param(NUM_SHEAR), alpha, get_param(NUM_SHEAR+2));}
		void set_geometry(double rad, double layer = 0.) { set_param(NUM_GEOMT, rad); set_param(NUM_GEOMT+1, rad+layer);}
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
//...аналитические модели (течение Бринкмана);
		double TakeLayer	 (double ff);
		double TakeCylinder(double ff, double eps = EE);
		double CylinderVelocity(double rr, double RR, double eps = EE);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     1.	//...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                2    //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_Q_facet               1.   //...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet                    //...param(3);
#define DRAFT_R0                    1.   //...radius of inclusion;
#undef  DRAFT_R0                         //...param(4);
#define DRAFT_R1							1.	//...radius of intermediate layer;
#undef  DRAFT_R1									//...param(5);
#define DRAFT_G0                    1.   //...dynamic shear modulus in fluid;
#undef  DRAFT_G0                         //...param(6);
#define DRAFT_alpha                 0.   //...normalized Brinkman;
#undef  DRAFT_alpha                      //...param(7);
#define DRAFT_kk                    1.   //...permability of intermediate layer;
#undef  DRAFT_kk                         //...param(8);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(9);
#endif
