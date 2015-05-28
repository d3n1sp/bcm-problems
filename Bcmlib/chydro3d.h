/*==========================================*/
/*                  CHYDRO3D                */
/*==========================================*/
#ifndef ___CHydro3D___
#define ___CHydro3D___

#include "shapes.h"
#include "ccomput3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}
int regul = 1;

/////////////////////////////////////////////////
//...class of blocks partition for Stokes problem;
template <typename T> 
class CHydro3D : public CComput3D<T> {
public:
		Num_Draft type   () { return HYDRO3D_DRAFT;}
		int size_of_param() { return(10);}
//...constructor and destructor;
		CHydro3D (int num_phase = 8) {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			this->param[this->NUM_MPLS] = 
			this->param[this->NUM_QUAD] = PackInts(2, 2);
			this->param[this->NUM_MPLS+1]  = 1.;
			this->param[this->NUM_GEOMT]	 = 1.;
			this->param[this->NUM_SHEAR]	 = 1.;
			this->param[this->NUM_SHEAR+2] = 1.;
			this->NUM_PHASE = num_phase;
			TT = TH = NULL; TN = NULL;
			this->BOX_LINK_PERIOD = OK_STATE;
		}
      virtual ~CHydro3D (void);
protected:
static int NUM_SHEAR, NUM_HESS, NUM_GEOMT, MAX_PHASE;
		T *** TT, *** TH, ** TN;
		int  block_shape_init(Block<T> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump0_pressure(double * P, int i, int m);
		void jump0_sphere  (double * P, int i, int m);
		void jump0_current (double * P, int i, int m);
		void jump0_darcy	 (double * P, int i, int m);
		void jump0			 (double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)){
				if (! this->get_param(this->NUM_SHEAR+1)) jump0_sphere(P, i, m); 
				else jump0_current(P, i, m);
			}
			else if ((this->B[i].type & ERR_CODE) == POLY_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
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
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump1_sphere_x(P, i, m); 
				else jump1_current_x(P, i, m);
			}
			else if ((this->B[i].type & ERR_CODE) == POLY_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_x(P, i, m);
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump1_fluid_x(P, i, m);
			else								 jump1_brinkmn_x(P, i, m);
		}
		void jump1_y(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump1_sphere_y(P, i, m); 
				else jump1_current_y(P, i, m);
			}
			else if ((this->B[i].type & ERR_CODE) == POLY_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_y(P, i, m);
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump1_fluid_y(P, i, m);
			else								 jump1_brinkmn_y(P, i, m);
		}
		void jump1_z(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump1_sphere_z(P, i, m); 
				else jump1_current_z(P, i, m);
			}
			else if ((this->B[i].type & ERR_CODE) == POLY_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) 
				  jump1_darcy_z(P, i, m);
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump1_fluid_z(P, i, m);
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
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump2_sphere_x(P, i, m); 
				else jump2_current_x(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump2_fluid_x(P, i, m);
			else								 jump2_brinkmn_x(P, i, m);
		}
		void jump2_y(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump2_sphere_y(P, i, m); 
				else jump2_current_y(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump2_fluid_y(P, i, m);
			else								 jump2_brinkmn_y(P, i, m);
		}
		void jump2_z(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump2_sphere_z(P, i, m); 
				else jump2_current_z(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump2_fluid_z(P, i, m);
			else								 jump2_brinkmn_z(P, i, m);
		}
		void jump3_fluid	(double * P, int i, int m);
		void jump3_sphere	(double * P, int i, int m);
		void jump3_brinkmn(double * P, int i, int m);
		void jump3_current(double * P, int i, int m);
		void jump3			(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump3_sphere(P, i, m); 
				else jump3_current(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump3_fluid(P, i, m);
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
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump4_sphere_x(P, i, m); 
				else jump4_current_x(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump4_fluid_x(P, i, m);
			else								 jump4_brinkmn_x(P, i, m);
		}
		void jump4_y(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump4_sphere_y(P, i, m); 
				else jump4_current_y(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump4_fluid_y(P, i, m);
			else								 jump4_brinkmn_y(P, i, m);
		}
		void jump4_z(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump4_sphere_z(P, i, m); 
				else jump4_current_z(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump4_fluid_z(P, i, m);
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
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump5_sphere_x(P, i, m); 
				else jump5_current_x(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump5_fluid_x(P, i, m);
			else								 jump5_brinkmn_x(P, i, m);
		}
		void jump5_y(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump5_sphere_y(P, i, m); 
				else jump5_current_y(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump5_fluid_y(P, i, m);
			else								 jump5_brinkmn_y(P, i, m);
		}
		void jump5_z(double * P, int i, int m) {
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK && this->B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				if (! this->get_param(this->NUM_SHEAR+1)) jump5_sphere_z(P, i, m); 
				else jump5_current_z(P, i, m);
			}
			else
			if (! this->get_param(this->NUM_SHEAR+1)) jump5_fluid_z(P, i, m);
			else								 jump5_brinkmn_z(P, i, m);
		}
		void hessian_deriv_N (double * P, int i);
		void hessian_deriv_N (int k, double * P, int i);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1_phase1(CGrid * nd, int i, int id_local);
		Num_State gram1_phase2(CGrid * nd, int i, int id_local);
		Num_State gram1		  (CGrid * nd, int i, int id_local) {
			if (this->B[i].link[this->NUM_PHASE] == -1) return gram1_phase1(nd, i, id_local); else
			if (this->B[i].link[this->NUM_PHASE] == -2) return gram1_phase2(nd, i, id_local); else return(ERR_STATE);
		}
		Num_State gram2	  (CGrid * nd, int i, int id_local);
		Num_State gram2peri(CGrid * nd, int i, int id_local);
		Num_State gram3	  (CGrid * nd, int i, int id_local);
		Num_State gram4_phase1(CGrid * nd, int i, int id_local);
		Num_State gram4_phase2(CGrid * nd, int i, int id_local);
		Num_State gram4		  (CGrid * nd, int i, int id_local) {
			if (this->B[i].link[this->NUM_PHASE] == -1) return gram4_phase1(nd, i, id_local); else
			if (this->B[i].link[this->NUM_PHASE] == -2) return gram4_phase2(nd, i, id_local); else return(ERR_STATE);
		}
		Num_State transfer1_phase1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer1_phase2(CGrid * nd, int i, int k, int id_local);
		Num_State transfer1			(CGrid * nd, int i, int k, int id_local) {
			if (this->B[i].link[this->NUM_PHASE] == -1 && this->B[k].link[this->NUM_PHASE] == -1) return transfer1_phase1(nd, i, k, id_local); else
			if (this->B[i].link[this->NUM_PHASE] == -2 && this->B[k].link[this->NUM_PHASE] == -2) return transfer1_phase2(nd, i, k, id_local); else return(ERR_STATE);
		}
		Num_State transfer2			(CGrid * nd, int i, int k, int id_local);
	 //Num_State transfer3			(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4_phase1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4_phase2(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4			(CGrid * nd, int i, int k, int id_local) {
			if (this->B[i].link[this->NUM_PHASE] == -1 && this->B[k].link[this->NUM_PHASE] == -1) return transfer4_phase1(nd, i, k, id_local); else
			if (this->B[i].link[this->NUM_PHASE] == -2 && this->B[k].link[this->NUM_PHASE] == -2) return transfer4_phase2(nd, i, k, id_local); else return(ERR_STATE);
		}
		Num_State rigidy1  (CGrid * nd, int i, T * K);
		Num_State rigidy2  (CGrid * nd, int i, T * K);
		Num_State computing_header (Num_Comput Num);
		Num_State computing_kernel5();//...счетная схема для явного аналитического решения;
public:
//...analytical method for boundary functions;
		void add_collocation(int i, int m, int j);
		void intern_matrix(T **& TT, int N, T **& TH);
		void extern_matrix(T **& TT, int N);
		void set_eshelby_matrix (int N_mpl);
		void set_normaliz_vector(int N_mpl, double rad, double kk);
		void set_normaliz_vector_old(int N_mpl, double rad, double kk);
//...параметры задачи;
		void set_fasa_hmg(double R0, double ll, double G_dynamic, double alpha, double kk);
		void set_fasa_hmg(double G_dynamic, double alpha) { set_fasa_hmg(this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1), G_dynamic, alpha, this->get_param(this->NUM_SHEAR+2));}
		void set_fasa_hmg(double alpha = 0.) {set_fasa_hmg(this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1), this->get_param(this->NUM_SHEAR), alpha, this->get_param(this->NUM_SHEAR+2));}
		void set_geometry(double rad, double layer = 0.) { set_param(this->NUM_GEOMT, rad); set_param(this->NUM_GEOMT+1, rad+layer);}
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, T * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
//...аналитические модели (течение Бринкмана);
		double TakeLayer	 (double ff);
		double TakeCylinder(double ff, double eps = EE);
		double CylinderVelocity(double rr, double RR, double eps = EE);
};

/////////////////////////////////////////////////
//...this->parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...this->param(0);
#define DRAFT_Q                     1.	//...normalization coefficient;
#undef  DRAFT_Q                          //...this->param(1);
#define DRAFT_N_quad                2    //...this->parameters of quadrature;
#undef  DRAFT_N_quad                     //...this->param(2);
#define DRAFT_Q_facet               1.   //...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet                    //...this->param(3);
#define DRAFT_R0                    1.   //...radius of inclusion;
#undef  DRAFT_R0                         //...this->param(4);
#define DRAFT_R1							1.	//...radius of intermediate layer;
#undef  DRAFT_R1									//...this->param(5);
#define DRAFT_G0                    1.   //...dynamic shear modulus in fluid;
#undef  DRAFT_G0                         //...this->param(6);
#define DRAFT_alpha                 0.   //...normalized Brinkman;
#undef  DRAFT_alpha                      //...this->param(7);
#define DRAFT_kk                    1.   //...permability of intermediate layer;
#undef  DRAFT_kk                         //...this->param(8);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...this->param(9);

//////////////////////////////////////////////////////
//        TEMPLATE VARIANT OF CHydro3D CLASS        //
//////////////////////////////////////////////////////
template <typename T> int CHydro3D<T>::NUM_HESS = 14;
template <typename T> int CHydro3D<T>::NUM_SHEAR = 6;
template <typename T> int CHydro3D<T>::NUM_GEOMT = 4;
template <typename T> int CHydro3D<T>::MAX_PHASE = 2;

//////////////////////////////////
//...destructor of CHydro3D class;
template <typename T>
CHydro3D<T>::~CHydro3D (void)
{
	for (int k = 0; TT && TT[k]; k++) {
		delete_struct(TT[k]);
		delete_struct(TH[k]);
	}
   delete_struct(TT);
   delete_struct(TH);
   delete_struct(TN);
}

//////////////////////////////////
//...initialization of the blocks;
template <typename T>
int CHydro3D<T>::block_shape_init(Block<T> & B, Num_State id_free)
{
	int k;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<T>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<T>(MP3D_ZOOM_SHAPE));

			B.shape->init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, sqr(B.mp[8])/(this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7])));

			if (this->get_param(this->NUM_SHEAR+1)) { //...подключаем две дополнительные экспоненциальные системы функций;
				B.shape->add_shape(CreateShape<T>(SK3D_EXPP_SHAPE),	 NULL_STATE);
				B.shape->add_shape(CreateShape<T>(SK3D_EXPP_SHAPE, 1), NULL_STATE);

				B.shape->init1(2, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, 0);
				B.shape->set_shape(2, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)));

				B.shape->init1(3, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, 0);
				B.shape->set_shape(3, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)));
			}
		}
		else
		if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<T>(MP3D_ZOOM_SHAPE), NULL_STATE);

			B.shape->init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, 0);
			B.shape->set_shape(1, sqr(B.mp[8] = this->get_param(this->NUM_GEOMT))/(this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7])));

			if (this->get_param(this->NUM_SHEAR+1)) { //...подключаем две дополнительные экспоненциальные системы функций;
				B.shape->add_shape(CreateShape<T>(SK3D_EXPP_SHAPE),	  NULL_STATE);
				B.shape->add_shape(CreateShape<T>(SK3D_EXPP_SHAPE, 1), NULL_STATE);

				B.shape->init1(2, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, 0);
				B.shape->set_shape(2, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), fabs(B.mp[8])); 

				B.shape->init1(3, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, 0);
				B.shape->set_shape(3, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), fabs(B.mp[8])); 
			}
		}
		else
		if ((B.type & ERR_CODE) == POLY_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<T>(MP3D_SPHERE_SHAPE));

			B.shape->init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, 1/*draft_dim(type())*/);
 			B.shape->set_shape(0, fabs(B.mp[7]), fabs(B.mp[8]/B.mp[7]));
		}
		else {
			B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE));

			B.shape->init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));

			if (this->get_param(this->NUM_SHEAR+1)) { //...подключаем дополнительную экспоненциальную систему функций;
				B.shape->add_shape(CreateShape<T>(SK3D_ZOOM_SHAPE), NULL_STATE);

				B.shape->init1(1, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, 0);
				B.shape->set_shape(1, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
			}
		}

/////////////////////////////////////////////////////
//...local system of coordinates and this->parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release();
   }

///////////////////////////////////////
//...setting this->parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, sqr(B.mp[8])/(this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7])));
				B.shape->set_shape(2, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), fabs(B.mp[8])); 
				B.shape->set_shape(3, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), fabs(B.mp[8])); 
			}
			else 
			if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, sqr(B.mp[8])/(this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7])));
			}
 			else {
				B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]), sqrt(this->get_param(this->NUM_SHEAR+1)), this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
			}
		}
		else
			if (id_free == OK_STATE && this->solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < this->solver.id_norm; m++)
				B.shape->set_potential(this->solver.hh[k][0][m], m);
		else
		if (id_free == NO_STATE && this->solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < this->solver.id_norm; m++)
				B.shape->get_potential(this->solver.hh[k][0][this->solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

/////////////////////////////////////////////////////////////
//...realization pressure (p) for fluid and brinkman current;
template <typename T>
void CHydro3D<T>::jump0_pressure(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, -1.);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, 0., 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->p_cpy, -1.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, -1.);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
}

////////////////////////////////////////////////////////////
//...realization of fluid pressure (p) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump0_sphere(double * P, int i, int m)
{
	int j = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, -1.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));

	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, 0., 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->p_cpy, -1.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, -1.);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization pressure for brinkman current arround sphere (p);
template <typename T>
void CHydro3D<T>::jump0_current(double * P, int i, int m)
{
	int nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(0); m += this->solver.id_norm;
	double RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
/////////////////////////////////////////
//...главная лапласовская часть давления;
	this->B[i].shape->cpy_x(0, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), this->B[i].shape->deriv, 0., -1.);

	this->B[i].shape->cpy_y(0, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->deriv, 0., -1.);

	this->B[i].shape->cpy_z(0, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->deriv, 0., -1.);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm), 1., RR);

	this->B[i].shape->cpy_y(nm, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm), 1., RR);

	this->B[i].shape->cpy_z(nm, this->B[i].shape->deriv);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm), 1., RR);
}

////////////////////////////////////////////////
//...realization pressure (p) for darcy current;
template <typename T>
void CHydro3D<T>::jump0_darcy(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, 0., 0.);
	this->B[i].shape->adm		 (this->B[i].shape->p_cpy, -.5/this->get_param(this->NUM_SHEAR+2));
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vx);
template <typename T>
void CHydro3D<T>::jump1_fluid_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_x     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[0]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2], 0.);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[0]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[0]*G1*this->B[i].shape->get_R_inv();
	}
}

///////////////////////////////////////////////////////////// 
//...realization of fluid velocity (Vx) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump1_sphere_x(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_x     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[0]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2], 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[0]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[0]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xx(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	this->B[i].shape->cpy_xy(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	this->B[i].shape->cpy_xz(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vx);
template <typename T>
void CHydro3D<T>::jump1_brinkmn_x(double * P, int i, int m)
{                                                                                         
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy_xx(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_xy(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_xz(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm	(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), G1);
		this->B[i].shape->adm_xx(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_xy(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_xz(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vx);
template <typename T>
void CHydro3D<T>::jump1_current_x(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
///////////////////////////////// 
//...лапласовская часть скорости;
	this->B[i].shape->cpy_xx(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy_xy(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy_xz(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xx(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_xy(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_xz(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., G1);

	this->B[i].shape->cpy_xx(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_xy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_xz(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., G1);

	this->B[i].shape->cpy_xx(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_xy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_xz(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vx);
template <typename T>
void CHydro3D<T>::jump1_darcy_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_x		 (this->B[i].shape->deriv, .5/this->get_param(this->NUM_SHEAR));
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vy);
template <typename T>
void CHydro3D<T>::jump1_fluid_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_y     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[1]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2], 0.);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[1]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[1]*G1*this->B[i].shape->get_R_inv();
	}
}

/////////////////////////////////////////////////////////////
//...realization of fluid velocity (Vy) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump1_sphere_y(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_y     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[1]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2], 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[1]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[1]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xy(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	this->B[i].shape->cpy_yy(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	this->B[i].shape->cpy_yz(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vy);
template <typename T>
void CHydro3D<T>::jump1_brinkmn_y (double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy_xy(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_yy(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_yz(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm	(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), G1);
		this->B[i].shape->adm_xy(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_yy(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_yz(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vy);
template <typename T>
void CHydro3D<T>::jump1_current_y(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	this->B[i].shape->cpy_xy(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy_yy(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy_yz(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xy(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_yy(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_yz(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., G1);

	this->B[i].shape->cpy_xy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_yy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_yz(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., G1);

	this->B[i].shape->cpy_xy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_yy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_yz(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vy);
template <typename T>
void CHydro3D<T>::jump1_darcy_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_y		 (this->B[i].shape->deriv, .5/this->get_param(this->NUM_SHEAR));
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vz);
template <typename T>
void CHydro3D<T>::jump1_fluid_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_z     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[2]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[1], 0.);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[2]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[2]*G1*this->B[i].shape->get_R_inv();
	}
}

/////////////////////////////////////////////////////////////
//...realization of fluid velocity (Vz) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump1_sphere_z(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_z     (this->B[i].shape->deriv);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, this->B[i].shape->deriv, 1., -P[2]);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, G1, 0.);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0], 0.);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[1], 0.);
	
	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[2]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[2]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xz(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	this->B[i].shape->cpy_yz(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	this->B[i].shape->cpy_zz(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vz);
template <typename T>
void CHydro3D<T>::jump1_brinkmn_z (double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy_xz(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_yz(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_zz(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->deriv, 0., alpha);
	}

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm	(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), G1);
		this->B[i].shape->adm_xz(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_yz(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_zz(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vz);
template <typename T>
void CHydro3D<T>::jump1_current_z(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	this->B[i].shape->cpy_xz(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy_yz(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy_zz(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_xz(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_yz(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_zz(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., G1);

	this->B[i].shape->cpy_xz(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_yz(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_zz(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., G1);

	this->B[i].shape->cpy_xz(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_yz(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_zz(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vz);
template <typename T>
void CHydro3D<T>::jump1_darcy_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, 0., 0.);
	this->B[i].shape->adm_z		 (this->B[i].shape->deriv, .5/this->get_param(this->NUM_SHEAR));
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVx);
template <typename T>
void CHydro3D<T>::jump2_fluid_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_xx	 (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[5]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, P[4]*G1);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, P[5]*G1);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[1]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->p_cpy, -P[4]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, -P[5]*G1);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[3]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[3]*G1*this->B[i].shape->get_R_inv();
	}
}

//////////////////////////////////////////////////////////////////////////////////// 
//...realization of normal derivative of fluid velocity (dnVx) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump2_sphere_x(double * P, int i, int m)
{
	int num_hess = NUM_HESS+this->solver.id_norm, j = num_hess+2, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_xx	 (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[5]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, P[4]*G1);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, P[5]*G1);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[1]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->p_cpy, -P[4]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, -P[5]*G1);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[3]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[3]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0); //...deriv_N_xx;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0); //...deriv_N_xz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVx);
template <typename T>
void CHydro3D<T>::jump2_brinkmn_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1); T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[3]*G1);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[4]*G1);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[5]*G1);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),		this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l),	 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l),	 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVx);
template <typename T>
void CHydro3D<T>::jump2_current_x(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), num_hess = NUM_HESS+this->solver.id_norm, k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
///////////////////////////////// 
//...лапласовская часть скорости;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0),    this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0), 0., 1.); //...deriv_N_xx;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 1), 0., 1.); //...deriv_N_xy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0), 0., 1.); //...deriv_N_xz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0),    this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], 0), alpha, -alpha*RR); //...deriv_N_xx;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], 0, 1), alpha, -alpha*RR); //...deriv_N_xy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], 0), alpha, -alpha*RR); //...deriv_N_xz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy_x(nm+1, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+1, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+1, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+1, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], 0); //...deriv_N_xx;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm+2, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+2, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+2, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+2, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], 0); //...deriv_N_xx;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVy);
template <typename T>
void CHydro3D<T>::jump2_fluid_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_yy	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[4], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[5]); 

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.); 
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, P[3]*G1);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, P[5]*G1);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->p_cpy, -P[3]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, -P[5]*G1);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[4]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[4]*G1*this->B[i].shape->get_R_inv();
	}
}

////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVy) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump2_sphere_y(double * P, int i, int m)
{
	int num_hess = NUM_HESS+this->solver.id_norm, j = num_hess+2, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_yy	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[4], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[5]); 

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.); 
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, P[3]*G1);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, P[5]*G1);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->p_cpy, -P[3]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, -P[5]*G1);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[4]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[4]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 2); //...deriv_N_yy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVy);
template <typename T>
void CHydro3D<T>::jump2_brinkmn_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1); T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[3]*G1);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[4]*G1);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[5]*G1);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),		this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 2), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVy);
template <typename T>
void CHydro3D<T>::jump2_current_y(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), num_hess = NUM_HESS+this->solver.id_norm, k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0),    this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 1), 0., 1.); //...deriv_N_xy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 2), 0., 1.); //...deriv_N_yy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 1), 0., 1.); //...deriv_N_yz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], 0, 1), alpha, -alpha*RR); //...deriv_N_xy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], 0, 2), alpha, -alpha*RR); //...deriv_N_yy;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], 0, 1), alpha, -alpha*RR); //...deriv_N_yz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy_x(nm+1, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+1, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+1, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+1, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], 0, 2); //...deriv_N_yy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm+2, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+2, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+2, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+2, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], 0, 2); //...deriv_N_yy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVz);
template <typename T>
void CHydro3D<T>::jump2_fluid_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy_zz	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[5], 0.);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[4]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, P[3]*G1);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, P[4]*G1);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->p_cpy, -P[3]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, -P[4]*G1);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
	}
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[5]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[5]*G1*this->B[i].shape->get_R_inv();
	}
}

////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVz) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump2_sphere_z(double * P, int i, int m)
{
	int num_hess = NUM_HESS+this->solver.id_norm, j = num_hess+2, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy_zz	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[5], 0.);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[4]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv, P[3]*G1);
	this->B[i].shape->adm_y     (this->B[i].shape->deriv, P[4]*G1);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->p_cpy, -P[3]*G1);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, -P[4]*G1);
	
	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[5]*G1*this->B[i].shape->get_R_inv();
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[5]*G1*this->B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0); //...deriv_N_xz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 2); //...deriv_N_zz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVz);
template <typename T>
void CHydro3D<T>::jump2_brinkmn_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1); T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[3]*G1);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[4]*G1);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[5]*G1);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),		this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l),	 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 2), 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVz);
template <typename T>
void CHydro3D<T>::jump2_current_z(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), num_hess = NUM_HESS+this->solver.id_norm, k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0),    this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0), 0., 1.); //...deriv_N_xz;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 1), 0., 1.); //...deriv_N_yz;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 2), 0., 1.); //...deriv_N_zz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], 0), alpha, -alpha*RR); //...deriv_N_xz;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], 0, 1), alpha, -alpha*RR); //...deriv_N_yz;
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], 0, 2), alpha, -alpha*RR); //...deriv_N_zz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy_x(nm+1, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+1, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+1, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+1, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], 0, 2); //...deriv_N_zz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm+2, this->B[i].shape->deriv); this->B[i].shape->admittance(nm+2, this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_y(nm+2, this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_z(nm+2, this->B[i].shape->deriv, P[5]); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., G1);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], 0, 2); //...deriv_N_zz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////////////////////
//...realization biharmonic potential of fluid velocity;
template <typename T>
void CHydro3D<T>::jump3_fluid(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	this->B[i].shape->cpy       (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2]*G1, 0.);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
}

/////////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential of fluid velocity for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump3_sphere(double * P, int i, int m)
{
	int num_hess = NUM_HESS+this->solver.id_norm, j = num_hess+2, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = .5/this->get_param(this->NUM_SHEAR), RR = sqr(this->B[i].mp[8]); T * p;
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));

	this->B[i].shape->cpy       (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	this->B[i].shape->cpy       (this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[2]*G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(0, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));
	this->B[i].shape->cpy(1, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0), p, 1., -G1);

	this->B[i].shape->cpy_y(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), p, 1., -G1);

	this->B[i].shape->cpy_z(1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman velocity;
template <typename T>
void CHydro3D<T>::jump3_brinkmn(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy_x(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_y(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->deriv, 0., alpha);

		this->B[i].shape->cpy_z(l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman current arround sphere;
template <typename T>
void CHydro3D<T>::jump3_current(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1), 
			 RR = this->get_param(this->NUM_MPLS+1)*fabs(this->B[i].mp[7])/fabs(this->B[i].mp[8]);
	T * p;

///////////////////////////////// 
//...лапласовская часть скорости;
	this->B[i].shape->cpy_x(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy_y(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy_z(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_y(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	this->B[i].shape->cpy_z(nm, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm);
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->cpy_x(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_y(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_z(nm+1, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->cpy_x(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., -alpha);

	this->B[i].shape->cpy_y(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	this->B[i].shape->cpy_z(nm+2, this->B[i].shape->deriv); p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (px);
template <typename T>
void CHydro3D<T>::jump4_fluid_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->cpy_xx	 (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[3], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[4]);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[5]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[0], 0.);
	this->B[i].shape->adm_x     (this->B[i].shape->deriv,	P[3]);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[1], 0.);
	//this->B[i].shape->adm_y(this->B[i].shape->p_cpy, P[3]);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2], 0.);
	//this->B[i].shape->adm_z(this->B[i].shape->deriv, P[3]);
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->p_cpy[3] = 0.;
		this->B[i].shape->deriv[3] = 0.;

		this->B[i].shape->p_cpy[2] = -P[3]*.5*this->B[i].shape->get_R_inv();
		this->B[i].shape->deriv[1] =  P[3]*.5*this->B[i].shape->get_R_inv();
	}
	this->B[i].shape->adm_y(this->B[i].shape->p_cpy, P[3]);
	this->B[i].shape->adm_z(this->B[i].shape->deriv, P[3]);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (px);
template <typename T>
void CHydro3D<T>::jump4_brinkmn_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double alpha = 1./this->get_param(this->NUM_SHEAR+1);	T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);
		this->B[i].shape->adm_x		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), P[3]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);
		this->B[i].shape->adm_y		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), P[3]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
		this->B[i].shape->adm_z		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), P[3]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[3]);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[4]);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[5]);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),	 this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l),	 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l),	 1., -alpha);
	}
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (py);
template <typename T>
void CHydro3D<T>::jump4_fluid_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->cpy_yy	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[4], 0.);
	this->B[i].shape->adm_xy    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[5]); 

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1], 0.); 
	this->B[i].shape->adm_y     (this->B[i].shape->deriv,	P[4]);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0], 0.);
	//this->B[i].shape->adm_x(this->B[i].shape->p_cpy, P[4]);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2], 0.);
	//this->B[i].shape->adm_z(this->B[i].shape->deriv, P[4]);
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->p_cpy[2] = 0.;
		this->B[i].shape->deriv[2] = 0.;

		this->B[i].shape->p_cpy[3] =  P[4]*.5*this->B[i].shape->get_R_inv();
		this->B[i].shape->deriv[1] = -P[4]*.5*this->B[i].shape->get_R_inv();
	}
	this->B[i].shape->adm_x(this->B[i].shape->p_cpy, P[4]);
	this->B[i].shape->adm_z(this->B[i].shape->deriv, P[4]);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (py);
template <typename T>
void CHydro3D<T>::jump4_brinkmn_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1); T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);
		this->B[i].shape->adm_x		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), P[4]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);
		this->B[i].shape->adm_y		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), P[4]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
		this->B[i].shape->adm_z		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), P[4]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[3]*G1);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[4]*G1);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[5]*G1);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),	 this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], l, 2), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
	}
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (pz);
template <typename T>
void CHydro3D<T>::jump4_fluid_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->cpy_zz	(this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, P[5], 0.);
	this->B[i].shape->adm_xz    (this->B[i].shape->deriv, P[3]);
	this->B[i].shape->adm_yz    (this->B[i].shape->deriv, P[4]);

	this->B[i].shape->cpy       (this->B[i].shape->deriv, this->B[i].shape->p_cpy);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[2], 0.);
	this->B[i].shape->adm_z     (this->B[i].shape->deriv, P[5]);
	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));

	this->B[i].shape->cpy       (this->B[i].shape->p_cpy, this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->p_cpy, NULL, -P[0], 0.);
	//this->B[i].shape->adm_x(this->B[i].shape->p_cpy, P[5]);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, -P[1], 0.);
	//this->B[i].shape->adm_y(this->B[i].shape->deriv, P[5]);
	if (this->B[i].shape->get_N() && regul) {
		this->B[i].shape->p_cpy[1] = 0.;
		this->B[i].shape->deriv[1] = 0.;

		this->B[i].shape->p_cpy[3] = -P[5]*.5*this->B[i].shape->get_R_inv();
		this->B[i].shape->deriv[2] =  P[5]*.5*this->B[i].shape->get_R_inv();
	}
	this->B[i].shape->adm_x(this->B[i].shape->p_cpy, P[5]);
	this->B[i].shape->adm_y(this->B[i].shape->deriv, P[5]);
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->p_cpy, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (pz);
template <typename T>
void CHydro3D<T>::jump4_brinkmn_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR), alpha = G1/this->get_param(this->NUM_SHEAR+1); T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), p, 0., alpha);
		this->B[i].shape->adm_x		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l), P[5]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), p, 0., alpha);
		this->B[i].shape->adm_y		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), P[5]);

		p = this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), p, 0., alpha);
		this->B[i].shape->adm_z		 (l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), P[5]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->adm_x(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[3]*G1);
		this->B[i].shape->adm_y(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[4]*G1);
		this->B[i].shape->adm_z(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[5]*G1);

		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l),		this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l),	 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
		this->B[i].shape->admittance(l, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2), this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], l, 2), 1., -alpha);
	}
}

////////////////////////////////////////////////////
//...realization primitive of vector potential (fx);
template <typename T>
void CHydro3D<T>::jump5_fluid_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l));

	if (this->B[i].shape->get_N() && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[0]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[0]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump5_sphere_x(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);
	
	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, -1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[3] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[3] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] = -P[0]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] =  P[0]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for brinkman velocity;
template <typename T>
void CHydro3D<T>::jump5_brinkmn_x(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->primitive (nm+l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), this->B[i].shape->deriv, 0., G1);
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for brinkman current arround sphere;
template <typename T>
void CHydro3D<T>::jump5_current_x(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR); T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->primitive(nm+1, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->primitive(nm+2, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0), p, 1., G1);
}

////////////////////////////////////////////////////
//...realization primitive of vector potential (fy);
template <typename T>
void CHydro3D<T>::jump5_fluid_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 1));

	if (this->B[i].shape->get_N(0) && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[1]*.5*this->B[i].shape->get_R_inv(0)*P[2]*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[1]*.5*this->B[i].shape->get_R_inv(0)*P[2]*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump5_sphere_y(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[2] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[2] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] =  P[1]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2)[1] = -P[1]*.5*this->B[i].shape->get_R_inv()*P[2]*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for brinkman velocity;
template <typename T>
void CHydro3D<T>::jump5_brinkmn_y(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->primitive (nm+l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), this->B[i].shape->deriv, 0., G1);
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for brinkman current arround sphere;
template <typename T>
void CHydro3D<T>::jump5_current_y(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR); T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->primitive(nm+1, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->primitive(nm+2, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), p, 1., G1);
}

///////////////////////////////////////////////////
//...realization primitiveof vector potential (fz);
template <typename T>
void CHydro3D<T>::jump5_fluid_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], l, 2));

	if (this->B[i].shape->get_N() && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[2]*.5*this->B[i].shape->get_R_inv()*P[2]*.5*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[2]*.5*this->B[i].shape->get_R_inv()*P[2]*.5*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fz) for sphere inclusion;
template <typename T>
void CHydro3D<T>::jump5_sphere_z(double * P, int i, int m)
{
	int j = NUM_HESS+2+this->solver.id_norm, N_mpl = this->B[i].shape->get_N(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

	this->B[i].shape->primitive (this->B[i].shape->deriv);
	this->B[i].shape->admittance(this->B[i].shape->deriv, NULL, G1, 0.);

	this->B[i].shape->cpy(0, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2));
	this->B[i].shape->cpy(1, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[1] = 0.;
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[1] = 0.;

		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0)[3] = -P[2]*.5*this->B[i].shape->get_R_inv()*P[2]*.5*P[5];
		this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1)[2] =  P[2]*.5*this->B[i].shape->get_R_inv()*P[2]*.5*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fz) for brinkman velocity;
template <typename T>
void CHydro3D<T>::jump5_brinkmn_z(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		this->B[i].shape->primitive (nm+l, this->B[i].shape->deriv);
		this->B[i].shape->admittance(nm+l, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), this->B[i].shape->deriv, 0., G1);
	}
}

//////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman current arround sphere;
template <typename T>
void CHydro3D<T>::jump5_current_z(double * P, int i, int m)
{
	int    nm = this->B[i].shape->num_usual(), N_mpl = this->B[i].shape->get_N(), k, j; m += this->solver.id_norm;
	double G1 = 1./this->get_param(this->NUM_SHEAR); T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	this->B[i].shape->primitive(nm+1, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	this->B[i].shape->primitive(nm+2, this->B[i].shape->deriv);	p = this->B[i].shape->FULL(this->B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	this->B[i].shape->admittance(0, this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), p, 1., G1);
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
template <typename T>
void CHydro3D<T>::hessian_deriv_N(double * P, int i)
{
	int l, nm = this->B[i].shape->num_usual(), num_hess = NUM_HESS+this->solver.id_norm;
	this->B[i].shape->cpy_xx(); this->B[i].shape->deriv_N(); 
	this->B[i].shape->cpy_xx();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], nm+l), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], nm*2+l), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], nm*3+l), nm*3, 0, -nm*3+l));
	}
	this->B[i].shape->cpy_xy(); this->B[i].shape->deriv_N();
	this->B[i].shape->cpy_xy();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 1));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], nm+l, 1), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], nm*2+l, 1), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], nm*3+l, 1), nm*3, 0, -nm*3+l));
	}
	this->B[i].shape->cpy_yy(); this->B[i].shape->deriv_N();
	this->B[i].shape->cpy_yy();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], l, 2));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+2], nm+l, 2), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+4], nm*2+l, 2), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+6], nm*3+l, 2), nm*3, 0, -nm*3+l));
	}
	this->B[i].shape->cpy_xz(); this->B[i].shape->deriv_N();
	this->B[i].shape->cpy_xz();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], nm+l), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], nm*2+l), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], nm*3+l), nm*3, 0, -nm*3+l));
	}
	this->B[i].shape->cpy_yz(); this->B[i].shape->deriv_N();
	this->B[i].shape->cpy_yz();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 1));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], nm+l, 1), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], nm*2+l, 1), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], nm*3+l, 1), nm*3, 0, -nm*3+l));
	}
	this->B[i].shape->cpy_zz(); this->B[i].shape->deriv_N();
	this->B[i].shape->cpy_zz();
	for (l = 0; l < nm; l++) {
		this->B[i].shape->cpy(l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], l, 2));
		this->B[i].shape->cpy(nm+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+3], nm+l, 2), nm, 0, -nm+l));
		this->B[i].shape->cpy(nm*2+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+5], nm*2+l, 2), nm*2, 0, -nm*2+l));
		this->B[i].shape->cpy(nm*3+l, this->B[i].shape->deriv, this->B[i].shape->FULL(this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+7], nm*3+l, 2), nm*3, 0, -nm*3+l));
	}
}

///////////////////////////////////////////////////////////
//...realization normal derivatives of hessian for dim = 0;
template <typename T>
void CHydro3D<T>::hessian_deriv_N(int k, double * P, int i)
{
	int num_hess = NUM_HESS+this->solver.id_norm;
	this->B[i].shape->cpy_xx(); this->B[i].shape->deriv_N(k); 
	this->B[i].shape->cpy_xx();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, -k));

	this->B[i].shape->cpy_xy(); this->B[i].shape->deriv_N(k);
	this->B[i].shape->cpy_xy();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 1-k));

	this->B[i].shape->cpy_yy(); this->B[i].shape->deriv_N(k);
	this->B[i].shape->cpy_yy();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess], 0, 2-k));

	this->B[i].shape->cpy_xz(); this->B[i].shape->deriv_N(k);
	this->B[i].shape->cpy_xz();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, -k));

	this->B[i].shape->cpy_yz(); this->B[i].shape->deriv_N(k);
	this->B[i].shape->cpy_yz();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 1-k));

	this->B[i].shape->cpy_zz(); this->B[i].shape->deriv_N(k);
	this->B[i].shape->cpy_zz();
	this->B[i].shape->cpy(k, this->B[i].shape->deriv, this->B[i].shape->FULL(this->solver.hh[i][0][num_hess+1], 0, 2-k));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
template <typename T>
void CHydro3D<T>::jump_make_local(int i, int m)
{
	m  += this->solver.id_norm;
	for (int j = 0; j < this->solver.dim[i]; j++) {
		T P[] = { this->solver.hh[i][0][m  ][j], 
					 this->solver.hh[i][0][m+1][j], 
					 this->solver.hh[i][0][m+2][j] };
		this->B[i].shape->norm_local_T(P);
		this->solver.hh[i][0][m  ][j] = P[0];
		this->solver.hh[i][0][m+1][j] = P[1];
		this->solver.hh[i][0][m+2][j] = P[2];
	}
}

template <typename T>
void CHydro3D<T>::jump_make_common(int i, int m)
{
	m  += this->solver.id_norm;
	for (int j = 0; j < this->solver.dim[i]; j++) {
		T P[] = { this->solver.hh[i][0][m  ][j], 
					 this->solver.hh[i][0][m+1][j], 
					 this->solver.hh[i][0][m+2][j] };
		this->B[i].shape->norm_common_T(P);
		this->solver.hh[i][0][m  ][j] = P[0];
		this->solver.hh[i][0][m+1][j] = P[1];
		this->solver.hh[i][0][m+2][j] = P[2];
	}
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для фазы свободного течения жидкости;
template <typename T>
Num_State CHydro3D<T>::gram1_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);
		
////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
				hy = nd->nY[l]*hx;
				hz = nd->nZ[l]*hx;
				hx *= nd->nX[l];
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

/////////////////////////////////////////////////////////////////
//...jump of all displacement moments and composition functional;
			this->B[i].shape->parametrization_hess(P, 1);
			if (p4 == SKEWS_BND-SPECIAL_BND) { //...краевые условия для скорости поперек потока на ячейке;
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				jump2_x(P, i, 0); this->solver.admittance(i, 0, G1); this->solver.admittance (i, 3, 0., 0, P[3]);
				jump2_y(P, i, 1); this->solver.admittance(i, 1, G1); this->solver.admittance (i, 3, 1., 1, P[4]);
				jump2_z(P, i, 2); this->solver.admittance(i, 2, G1); this->solver.admittance (i, 3, 1., 2, P[5]);
				jump_make_common(i, 0);
				this->solver.admittance (i, 0, 1., 3, -nd->nX[l]);
				this->solver.admittance (i, 1, 1., 3, -nd->nY[l]);
				this->solver.admittance (i, 2, 1., 3, -nd->nZ[l]);
				jump1_x(P, i, 4); this->solver.admittance (i, 3, 0., 4, P[3]);
				jump1_y(P, i, 5); this->solver.admittance (i, 3, 1., 5, P[4]);
				jump1_z(P, i, 6); this->solver.admittance (i, 3, 1., 6, P[5]);

				this->solver.to_equationDD(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m],   f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+1], f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+2], f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
			}
			else
			if (p4 == MAX_HIT || p4 == NUMS_BND) { //...краевые условия для давления вдоль потока на ячейке;
				jump0(P, i, 3);
 				this->solver.to_equationDD(i, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
				
				if (fabs(hx) > EE) 
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+3], hx*f);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump1_x(P, i, 0); this->solver.admittance(i, 0, G1);
				jump1_y(P, i, 1); this->solver.admittance(i, 1, G1);
				jump1_z(P, i, 2); this->solver.admittance(i, 2, G1); jump_make_common(i, 0);

				this->solver.to_equationDD(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m],   f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+1], f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+2], f);
				
				if (fabs(hx) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m], hx*(f *= G1));
				if (fabs(hy) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+1], hy*f);
				if (fabs(hz) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+2], hz*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


//////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для фазы фильтрационного течения жидкости;
template <typename T>
Num_State CHydro3D<T>::gram1_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double hx, p4, f, P[6];
		int 	 m  = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);
		
////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == 0.) hx = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

/////////////////////////////////////////////////////////////////
//...jump of all displacement moments and composition functional;
			this->B[i].shape->parametrization_hess(P, 1);
			jump1_x(P, i, 0); this->solver.admittance (i, 3, 0., 0, P[3]);
			jump1_y(P, i, 1); this->solver.admittance (i, 3, 1., 1, P[4]);
			jump1_z(P, i, 2); this->solver.admittance (i, 3, 1., 2, P[5]);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
			if (fabs(hx) > EE) 
				this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+3], hx*f);

		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
template <typename T>
Num_State CHydro3D<T>::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz, 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1;
      int id_isolated = 1, m  = this->solver.id_norm, id_dir, k; int j;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hz = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hz = -AZ; break;
					case 6: TZ = -AZ; hz =  AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[3] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				
				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);

				this->solver.admittance(i, 0, g1, 4, g2); this->solver.admittance(i, 4, g0, 0, f1); 
				this->solver.admittance(i, 1, g1, 5, g2); this->solver.admittance(i, 5, g0, 1, f1); 
				this->solver.admittance(i, 2, g1, 6, g2); this->solver.admittance(i, 6, g0, 2, f1);

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);
				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);
				
				this->B[k].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6); jump_make_common(k, 4);

				this->solver.admittance(k, 0, g1, 4, g2); this->solver.admittance(k, 4, g0, 0, f1); 
				this->solver.admittance(k, 1, g1, 5, g2); this->solver.admittance(k, 5, g0, 1, f1); 
				this->solver.admittance(k, 2, g1, 6, g2); this->solver.admittance(k, 6, g0, 2, f1);

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m],   f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+4], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+4], this->solver.hh[k][0][m+4], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+5], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+5], this->solver.hh[k][0][m+5], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+6], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+6], this->solver.hh[k][0][m+6], f);

				if (abs(id_dir) < 5) { //...давление с учетом регуляризации;
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				}
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);

				if (fabs(hz) > EE) {
  				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+3], -hz*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m+3],  hz*f);
				}
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
			}
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для периодической задачи (один блок);
template <typename T>
Num_State CHydro3D<T>::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), AX, AY, AZ, f, P[6], 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1;
      int id_isolated = 1, m  = this->solver.id_norm, id_dir;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);
	           
////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l] && 
			((id_dir = (int)nd->get_param(3, l)) == 1 || id_dir == 3 || id_dir == 5)) {
			P[0] = nd->X[l]; P[3] = nd->nX[l];
			P[1] = nd->Y[l]; P[4] = nd->nY[l];
			P[2] = nd->Z[l]; P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

////////////////////////////////////////////////////
//...вычисляем функции для формирования функционала;
			this->B[i].shape->parametrization_hess(P, 1);
			if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
			
			jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
			jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);

			this->B[i].shape->make_common(P);
			this->B[i].shape->norm_common(P+3);

			if (id_dir == 1) this->B[i].mp[1] -= AX; else
			if (id_dir == 3) this->B[i].mp[2] -= AY; else
			if (id_dir == 5) this->B[i].mp[3] -= AZ; 
			this->B[i].shape->set_local_P0(this->B[i].mp+1);
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);
			
			this->B[i].shape->parametrization_hess(P, 1);
			if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			jump1_x(P, i,  7); jump1_y(P, i,  8); jump1_z(P, i,  9); jump_make_common(i,  7); jump0(P, i, 10);
			jump2_x(P, i, 11); jump2_y(P, i, 12); jump2_z(P, i, 13); jump_make_common(i, 11);

			this->solver.admittance(i, 0, 1.,  7, -1.); this->solver.admittance(i,  4, 1., 11, -1.);
			this->solver.admittance(i, 1, 1.,  8, -1.); this->solver.admittance(i,  5, 1., 12, -1.);
			this->solver.admittance(i, 2, 1.,  9, -1.); this->solver.admittance(i,  6, 1., 13, -1.); 
			this->solver.admittance(i, 3, 1., 10, -1.); this->solver.admittance(i, 10, 2.,  3,  1.);

			this->solver.admittance(i, 0, g1,  4, g2);  this->solver.admittance(i,  4, g0, 0, f1); 
			this->solver.admittance(i, 1, g1,  5, g2);  this->solver.admittance(i,  5, g0, 1, f1); 
			this->solver.admittance(i, 2, g1,  6, g2);  this->solver.admittance(i,  6, g0, 2, f1); 

/////////////////////////////////////////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов (скорость, давление, производная);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m],   f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+1], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+2], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+4], this->solver.hh[i][0][m+4], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+5], this->solver.hh[i][0][m+5], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+6], this->solver.hh[i][0][m+6], f);

			if (id_dir == 5) {//...правая часть -- скачок давления и регуляризация матрицы;
				this->solver.to_equationHH(i, 0,	this->solver.hh[i][0][m+3], AZ*f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+10], this->solver.hh[i][0][m+10], f);
			}
			if (id_dir == 1) this->B[i].mp[1] += AX; else
			if (id_dir == 3) this->B[i].mp[2] += AY; else
			if (id_dir == 5) this->B[i].mp[3] += AZ; 
			this->B[i].shape->set_local_P0(this->B[i].mp+1);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе свободное - свободное течение жидкости;
template <typename T>
Num_State CHydro3D<T>::transfer1_phase1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N) {
      double G1 = this->get_param(this->NUM_SHEAR), f, P[6], 
				 f0 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int id_isolated = 0, m = this->solver.id_norm;
		if (id_isolated) {
			g0 = g1 = G1;
			f0 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					//this->solver.hh[i][0].clear_row_matrix(num, this->solver.dim[i]);
					//this->solver.hh[k][0].clear_row_matrix(num, this->solver.dim[k]);
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				
				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);	

				this->solver.admittance (i, 0, g1, 4, g2); this->solver.admittance(i, 4, g0, 0, f0); 
				this->solver.admittance (i, 1, g1, 5, g2); this->solver.admittance(i, 5, g0, 1, f0); 
				this->solver.admittance (i, 2, g1, 6, g2); this->solver.admittance(i, 6, g0, 2, f0);

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);
				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);
				
				this->B[k].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0);	jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6); jump_make_common(k, 4);

				this->solver.admittance (k, 0, g1, 4, g2); this->solver.admittance(k, 4, g0, 0, f0); 
				this->solver.admittance (k, 1, g1, 5, g2); this->solver.admittance(k, 5, g0, 1, f0); 
				this->solver.admittance (k, 2, g1, 6, g2); this->solver.admittance(k, 6, g0, 2, f0); 

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m],   f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+4], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+4], this->solver.hh[k][0][m+4], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+5], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+5], this->solver.hh[k][0][m+5], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+6], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+6], this->solver.hh[k][0][m+6], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


//////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе фильтрационное - фильтрационное течение жидкости;
template <typename T>
Num_State CHydro3D<T>::transfer1_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N) {
      int m = this->solver.id_norm;
      double f, P[6];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				jump1_x(P, i, 0); this->solver.admittance (i, 3, 0., 0, P[3]);
				jump1_y(P, i, 1); this->solver.admittance (i, 3, 1., 1, P[4]);
				jump1_z(P, i, 2); this->solver.admittance (i, 3, 1., 2, P[5]);
				
				this->B[k].shape->parametrization_hess(P, 1);
				jump0(P, k, 3);

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе свободное - фильтрационное течение жидкости;
template <typename T>
Num_State CHydro3D<T>::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < this->solver.N && this->B && this->B[i].link && this->B[i].link[0] > this->NUM_PHASE) {
      double f, P[6], beta = 0.1;
      int m = this->solver.id_norm;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

				f = nd->get_param(0, l);

				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

///////////////////////////////////////////////
//...jump of all neaded moments for this block;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				if (this->B[i].link[this->NUM_PHASE] == -1) {
					jump1_x(P, i, 0); this->solver.admittance (i, 3, 0., 0, P[3]);
					jump1_y(P, i, 1); this->solver.admittance (i, 3, 1., 1, P[4]);
					jump1_z(P, i, 2); this->solver.admittance (i, 3, 1., 2, P[5]);
					this->solver.admittance (i, 0, 1., 3, -P[3]);
					this->solver.admittance (i, 1, 1., 3, -P[4]);
					this->solver.admittance (i, 2, 1., 3, -P[5]);

					jump4_x(P, i, 4); this->solver.admittance (i, 7, 0., 4, P[3]);
					jump4_y(P, i, 5); this->solver.admittance (i, 7, 1., 5, P[4]);
					jump4_z(P, i, 6); this->solver.admittance (i, 7, 1., 6, P[5]);
					this->solver.admittance (i, 4, 1., 7, -P[3]); this->solver.admittance (i, 4, 1., 0, beta);
					this->solver.admittance (i, 5, 1., 7, -P[4]); this->solver.admittance (i, 5, 1., 1, beta);
					this->solver.admittance (i, 6, 1., 7, -P[5]); this->solver.admittance (i, 6, 1., 2, beta);
					jump0(P, i, 7);
				}
				if (this->B[i].link[this->NUM_PHASE] == -2) {
					jump1_x(P, i, 0); this->solver.admittance (i, 3, 0., 0, P[3]);
					jump1_y(P, i, 1); this->solver.admittance (i, 3, 1., 1, P[4]);
					jump1_z(P, i, 2); this->solver.admittance (i, 3, 1., 2, P[5]);
					jump0  (P, i, 7); 
					this->solver.admittance (i, 0, beta, 3, -P[3]*beta);
					this->solver.admittance (i, 1, beta, 3, -P[4]*beta);
					this->solver.admittance (i, 2, beta, 3, -P[5]*beta);
				}

///////////////////////////////////////////////////////
//...jump of all neaded moments for neighbouring block;
				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);
				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);
	
				this->B[k].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				if (this->B[k].link[this->NUM_PHASE] == -1) {
					jump1_x(P, k, 0); this->solver.admittance (k, 3, 0., 0, P[3]);
					jump1_y(P, k, 1); this->solver.admittance (k, 3, 1., 1, P[4]);
					jump1_z(P, k, 2); this->solver.admittance (k, 3, 1., 2, P[5]);
					this->solver.admittance (k, 0, 1., 3, -P[3]);
					this->solver.admittance (k, 1, 1., 3, -P[4]);
					this->solver.admittance (k, 2, 1., 3, -P[5]);

					jump4_x(P, k, 4); this->solver.admittance (k, 7, 0., 4, P[3]);
					jump4_y(P, k, 5); this->solver.admittance (k, 7, 1., 5, P[4]);
					jump4_z(P, k, 6); this->solver.admittance (k, 7, 1., 6, P[5]);
					this->solver.admittance (k, 4, 1., 7, -P[3]); this->solver.admittance (k, 4, 1., 0, beta);
					this->solver.admittance (k, 5, 1., 7, -P[4]); this->solver.admittance (k, 5, 1., 1, beta);
					this->solver.admittance (k, 6, 1., 7, -P[5]); this->solver.admittance (k, 6, 1., 2, beta);
					jump0(P, k, 7);
				}
				if (this->B[k].link[this->NUM_PHASE] == -2) {
					jump1_x(P, k, 0); this->solver.admittance (k, 3, 0., 0, P[3]);
					jump1_y(P, k, 1); this->solver.admittance (k, 3, 1., 1, P[4]);
					jump1_z(P, k, 2); this->solver.admittance (k, 3, 1., 2, P[5]);
					jump0	 (P, k, 7); 
					this->solver.admittance (k, 0, beta, 3, -P[3]*beta);
					this->solver.admittance (k, 1, beta, 3, -P[4]*beta);
					this->solver.admittance (k, 2, beta, 3, -P[5]*beta);
				}

//////////////////////////////////////////////////////////////
//...сшивка скоростей и давлений методом наименьших квадратов;
				if (this->B[i].link[this->NUM_PHASE] == -1 && this->B[k].link[this->NUM_PHASE] == -2) {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+4], this->solver.hh[k][0][m],   f);
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+5], this->solver.hh[k][0][m+1], f);
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+6], this->solver.hh[k][0][m+2], f);
					this->solver.to_transferTT(i, j, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
					this->solver.to_transferTT(i, j, this->solver.hh[i][0][m+4], this->solver.hh[i][0][m+4], f);
					this->solver.to_transferTT(i, j, this->solver.hh[i][0][m+5], this->solver.hh[i][0][m+5], f);
					this->solver.to_transferTT(i, j, this->solver.hh[i][0][m+6], this->solver.hh[i][0][m+6], f);
					this->solver.to_transferTD(i, j, this->solver.hh[k][0][m+7], this->solver.hh[k][0][m+7], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+7], this->solver.hh[k][0][m+7], f);
				}
				if (this->B[k].link[this->NUM_PHASE] == -1 && this->B[i].link[this->NUM_PHASE] == -2) {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+7], this->solver.hh[k][0][m+7], f);
					this->solver.to_transferTT(i, j, this->solver.hh[i][0][m+7], this->solver.hh[i][0][m+7], f);
					this->solver.to_transferTD(i, j, this->solver.hh[k][0][m+3], this->solver.hh[k][0][m+3], f);
					this->solver.to_transferTD(i, j, this->solver.hh[k][0][m+4], this->solver.hh[k][0][m+4], f);
					this->solver.to_transferTD(i, j, this->solver.hh[k][0][m+5], this->solver.hh[k][0][m+5], f);
					this->solver.to_transferTD(i, j, this->solver.hh[k][0][m+6], this->solver.hh[k][0][m+6], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+4], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+5], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+6], f);
				}
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
template <typename T>
Num_State CHydro3D<T>::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz;
      int	 m  = this->solver.id_norm, id_dir, k; int j;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hz = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hz = -G1*AZ; break;
					case 6: TZ = -AZ; hz =  G1*AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[2] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
				this->solver.admittance(i, 0, G1); this->solver.admittance(i, 1, G1); this->solver.admittance(i, 2, G1); 
				this->solver.admittance(i, 7, 0., 0, nd->nX[l]); 
				this->solver.admittance(i, 7, 1., 1, nd->nY[l]); 
				this->solver.admittance(i, 7, 1., 2, nd->nZ[l]);

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);
				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6);	jump_make_common(k, 4);
				this->solver.admittance(k, 0, G1); this->solver.admittance(k, 1, G1); this->solver.admittance(k, 2, G1); 
				this->solver.admittance(k, 7, 0., 0, nd->nX[l]); 
				this->solver.admittance(k, 7, 1., 1, nd->nY[l]); 
				this->solver.admittance(k, 7, 1., 2, nd->nZ[l]);

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);

				if (fabs(hz) > EE) {
				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+3], -hz*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m+3],  hz*f);
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				this->solver.to_equationER(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m+4],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m], this->solver.hh[k][0][m+4], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+5],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+1], this->solver.hh[k][0][m+5], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+6],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+2], this->solver.hh[k][0][m+6], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+7], this->solver.hh[i][0][m+3], -f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+7], this->solver.hh[k][0][m+3],  f);

////////////////////////////////////
//...естественные граничные условия;
				if (fabs(hz) > EE) {
					this->solver.to_equationEH(i, 0, this->solver.hh[i][0][m+7],  hz*f);
					this->solver.to_equationEH(k, 0, this->solver.hh[k][0][m+7], -hz*f);
				}				
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
			}
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для фазы свободного течения жидкости;
template <typename T>
Num_State CHydro3D<T>::gram4_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

//////////////////////////////////////
//...jump of all displacement moments;
			this->B[i].shape->parametrization_hess(P, 1);
			if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
			jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
			this->solver.admittance(i, 0, G1); this->solver.admittance(i, 1, G1); this->solver.admittance(i, 2, G1); 
			this->solver.admittance(i, 7, 0., 0, nd->nX[l]); 
			this->solver.admittance(i, 7, 1., 1, nd->nY[l]); 
			this->solver.admittance(i, 7, 1., 2, nd->nZ[l]);

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p4 == MAX_HIT || p4 == NUMS_BND) {
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+3], this->solver.hh[i][0][m+3], f);
				
				if (fabs(hx) > EE) 
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+3], hx*f);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				this->solver.to_equationDD(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m],   f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+1], f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+2], f);

				if (fabs(hx) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m],   hx*G1*f);
				if (fabs(hy) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+1], hy*G1*f);
				if (fabs(hz) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+2], hz*G1*f);
			}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			this->solver.to_equationEE(i, this->solver.hh[i][0][m],	 this->solver.hh[i][0][m+4],  f);
			this->solver.to_equationEE(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+5],  f);
			this->solver.to_equationEE(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+6],  f);
			this->solver.to_equationEE(i, this->solver.hh[i][0][m+7], this->solver.hh[i][0][m+3], -f);

////////////////////////////////////
//...естественные граничные условия;
			if ((p4 == MAX_HIT || p4 == NUMS_BND) && fabs(hx) > EE)
				this->solver.to_equationEH(i, 0, this->solver.hh[i][0][m+7], -hx*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для фазы фильтрационного течения жидкости;
template <typename T>
Num_State CHydro3D<T>::gram4_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)this->solver.N && this->B[i].shape && this->B[i].mp) {
		double G1 = this->get_param(this->NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

//////////////////////////////////////
//...jump of all displacement moments;
			this->B[i].shape->parametrization_hess(P, 1);
			if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе свободное - свободное течение жидкости;
template <typename T>
Num_State CHydro3D<T>::transfer4_phase1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < (int)this->solver.N) {
      double G1 = this->get_param(this->NUM_SHEAR), f, P[6];
      int m = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
				this->solver.admittance(i, 0, G1); this->solver.admittance(i, 1, G1); this->solver.admittance(i, 2, G1); 
				this->solver.admittance(i, 7, 0., 0, nd->nX[l]); 
				this->solver.admittance(i, 7, 1., 1, nd->nY[l]); 
				this->solver.admittance(i, 7, 1., 2, nd->nZ[l]);

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);
				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6);	jump_make_common(k, 4);
				this->solver.admittance(k, 0, G1); this->solver.admittance(k, 1, G1); this->solver.admittance(k, 2, G1); 
				this->solver.admittance(k, 7, 0., 0, nd->nX[l]); 
				this->solver.admittance(k, 7, 1., 1, nd->nY[l]); 
				this->solver.admittance(k, 7, 1., 2, nd->nZ[l]);

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);

				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				this->solver.to_equationER(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m+4],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m],   this->solver.hh[k][0][m+4], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+5],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+1], this->solver.hh[k][0][m+5], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+6],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+2], this->solver.hh[k][0][m+6], -f);

				this->solver.to_equationER(i, this->solver.hh[i][0][m+7], this->solver.hh[i][0][m+3], -f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m+7], this->solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе фильтрационное - фильтрационное течение жидкости;
template <typename T>
Num_State CHydro3D<T>::transfer4_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < (int)this->solver.N) {
      double G1 = this->get_param(this->NUM_SHEAR), f, P[6];
      int m = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////
//...интегрирование скоростей фильтрации по границе блоков;
template <typename T>
Num_State CHydro3D<T>::rigidy1(CGrid * nd, int i, T * K)
{
	if (nd) {
      int l, m = this->solver.id_norm;
      double f, P[6]; T V[3];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);
////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

///////////////////////////////////////////////////////////////////////////////////
//...вычисляем скорости в сепарированном виде для формирования интеграла по объему;
			this->B[i].shape->parametrization_hess(P, 1);
			jump3  (P, i, 0); 
			jump5_x(P, i, 1); //...нужно ли приведение нормали к общей системе координат???
			jump5_y(P, i, 2); 
			jump5_z(P, i, 3); 

			V[0] = this->B[i].shape->potential(this->solver.hh[i][0][m], 0);
			K[0] += V[0]*nd->nX[l]*f;
			K[1] += V[0]*nd->nY[l]*f;
			K[2] += V[0]*nd->nZ[l]*f;

			V[0] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], 0);
			V[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], 0);
			V[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+3], 0);
			this->B[i].shape->norm_common_T(V);

			K[0] += V[0]*f;//...нормаль убрали внутрь функций;
			K[1] += V[1]*f;
			K[2] += V[2]*f;
			K[3+(-this->B[i].link[this->NUM_PHASE]-1)] += T(nd->Z[l]*nd->nZ[l]*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);}

//////////////////////////////////////////////////////////
//...интегрирование скоростей фильтрации по объему блоков;
template <typename T>
Num_State CHydro3D<T>::rigidy2(CGrid * nd, int i, T * K)
{
	if (nd) {
      int m = this->solver.id_norm, l;
      double f, P[6]; T V[3];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 0., 0., 0., AXIS_Z, 1);
/////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = -1.;
			P[1] = nd->Y[l];  P[4] = 0.;
			P[2] = nd->Z[l];  P[5] = 0.;
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

///////////////////////////////////
//...вычисляем скорости фильтрации;
			this->B[i].shape->parametrization_hess(P, 1);

			jump1_x(P, i, 0); 
			jump1_y(P, i, 1); 
			jump1_z(P, i, 2); 

			V[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   0);
			V[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], 0);
			V[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], 0);
			this->B[i].shape->norm_common_T(V);
			
			K[0] += V[0]*f;
			K[1] += V[1]*f;
			K[2] += V[2]*f;
			K[3] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////
//...istalling this->parameters of doubly connected media;
template <typename T>
void CHydro3D<T>::set_fasa_hmg(double R0, double ll, double G_dynamic, double alpha, double kk)
{
	if (size_of_param() > this->NUM_GEOMT+1) {
		this->param[this->NUM_GEOMT] = R0;
		this->param[this->NUM_GEOMT+1] = R0+ll;
	}
	if (size_of_param() > this->NUM_SHEAR+1) {
		this->param[this->NUM_SHEAR]	 = G_dynamic;
		this->param[this->NUM_SHEAR+1] = alpha;
		this->param[this->NUM_SHEAR+2] = kk;
	}
	return;
}

///////////////////////////////////////////////////////////
//...counting header for solving linear elasticity problem;
template <typename T>
Num_State CHydro3D<T>::computing_header(Num_Comput Num)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), id_dir, i, k, l, elem, n_rhs = 2;
	char msg[201];
	
	if (! this->solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CHydro3D_draft: N_sm = %d, N_mpl = %d, N_elem = %d, alpha = %g", this->N, UnPackInts(this->get_param(this->NUM_MPLS)), N_elem, this->get_param(this->NUM_SHEAR+1));
		Message(msg);

		Message(" ");
		switch (Num){
			case	 BASIC_COMPUT: Message("Junction Blocks...");	break;
			case MAPPING_COMPUT: Message("Mapping  Blocks..."); break;
			case  PERIOD_COMPUT: Message("Periodic Blocks..."); break;
		}
		Message(" ");
	}
	if (this->N == 1 && (this->B[0].type & ERR_CODE) == CLAYER_BLOCK) {
		set_eshelby_matrix (UnPackInts(this->get_param(this->NUM_MPLS)));
		set_normaliz_vector_old(UnPackInts(this->get_param(this->NUM_MPLS)), this->get_param(this->NUM_GEOMT), sqrt(this->get_param(this->NUM_SHEAR+1)));
	}

///////////////////////
//...блочную структуру;
	this->solver.set_blocks(this->N, n_rhs); //<==== number of saved potentials !!!
	this->solver.n += NUM_HESS+9;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < this->solver.N;  k++)
		  this->solver.set_links(k, this->B[k].link);

	this->shapes_init(INITIAL_STATE);
	this->shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (this->solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; this->SetGeomBounding(par);
		for (k = 0; k < this->N; k++) SkeletonBounding(this->B[k], par);
		for (k = 0; k < this->N; k++) if (this->B[k].link) {
			for (i = 0; i < this->B[k].link[0]; i++) if ((elem = geom_plink_3D(this->B[k], l = i, id_dir, par)) >= 0) 
			this->solver.add_link(k, elem);
			for (i = 0; i < this->B[k].link[0]; i++) if ((elem = block_plink_3D(this->B[k], l = i, id_dir, par)) >= 0) 
			this->solver.add_link(k, elem);
		}
	}
	this->LinkPhase3D(MAX_PHASE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! this->solver.struct_permutat(this->solver.id_change == EXTERN_STATE ? NULL_STATE : /*NULL*/OK_STATE) || 
		 ! this->solver.inverse_index()) {
		return ERR_STATE;
	}
	for (k = 0; k < (int)this->solver.N; k++)
		this->solver.set_dimension(k, this->freedom_block(k));

	this->solver.struct_init();

	if (this->solver.mode(FULLY_MODE)) { 
		this->solver.test_struct("1", 0);
		this->solver.test_struct("2", 1);
	}
   return OK_STATE;
}

////////////////////////////////////////////////////////////
//...счетная схема для явного аналитического решения задачи;
template <typename T>
Num_State CHydro3D<T>::computing_kernel5()
{
	this->set_param(this->NUM_MPLS, PackInts(2, 2)); //...multipoles degree;
	if (computing_header(SPECIAL3COMPUT) != OK_STATE) 
		return ERR_STATE;

///////////////////////////////////////////////////////////////
//...переносим результат в функции формы (задаем решение явно);
	this->shapes_init(OK_STATE);
	for (int i = 0; i < 1; i++) {
		memset(this->solver.hh[i][0][0], 0, this->solver.dim[i]*sizeof(T));
		this->B[i].shape->FULL(this->solver.hh[i][0][0], 0, 0)[6] = T(0.5*sqr(this->B[i].shape->get_R()));	
		this->B[i].shape->FULL(this->solver.hh[i][0][0], 0, 1)[5] = T(0.5*sqr(this->B[i].shape->get_R()));	
		this->B[i].shape->FULL(this->solver.hh[i][0][0], 0, 2)[0] = T(0.5);	
		this->B[i].shape->set_potential(this->solver.hh[i][0][0], 0);
	}
   return OK_STATE;
}

//////////////////////////////////////////////////////////////
//...сложение коллокационных векторов с учетом матрицы Эшелби;
template <typename T>
void CHydro3D<T>::add_collocation(int i, int m, int j)
{
	if (TT && TH) {
		T * px = this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 0), 
				  * py = this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 1), 
				  * pz = this->B[i].shape->FULL(this->solver.hh[i][0][j], 0, 2), * p, 
					ff = fabs(this->B[i].mp[8])/this->B[i].shape->get_R(0), ff_inv = 1./ff;
		int   N_mpl = this->B[i].shape->get_N(0), nn, mm, kk, k, l, ll;
		for (k = N_mpl; k >= 0; k--) {
			nn = 2*k+1; kk = k*k; mm = nn+4; ll = (k+2)*(k+2);
			if (TT[k]) //...вычисляем текущую степень;
			for (l = 0; l < nn; l++) {				
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l]+py[kk+j]*TT[k][j+nn][l]+pz[kk+j]*TT[k][j+nn*2][l]);
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l+nn]+py[kk+j]*TT[k][j+nn][l+nn]+pz[kk+j]*TT[k][j+nn*2][l+nn]);
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l+nn*2]+py[kk+j]*TT[k][j+nn][l+nn*2]+pz[kk+j]*TT[k][j+nn*2][l+nn*2]);
			}
			if (TH[k] && k+2 <= N_mpl) //...правим более высокую степень;
			for (l = 0; l < mm; l++) {
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 0), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l]+py[kk+j]*TH[k][j+nn][l]+pz[kk+j]*TH[k][j+nn*2][l]);
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 1), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l+mm]+py[kk+j]*TH[k][j+nn][l+mm]+pz[kk+j]*TH[k][j+nn*2][l+mm]);
				for ( p = this->B[i].shape->FULL(this->solver.hh[i][0][m], 0, 2), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l+mm*2]+py[kk+j]*TH[k][j+nn][l+mm*2]+pz[kk+j]*TH[k][j+nn*2][l+mm*2]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
//...формирование внутренних матриц в аналитическом методе для граничных функций;
template <typename T>
void CHydro3D<T>::intern_matrix(T **& TT, int N, T **& TH)
{
	int l, m, n = N, nn = 2*n+1, mm = nn-4; 
	if (! TT)   set_matrix(TT, 3*nn, 3*nn);
	else		clean_matrix(TT, 3*nn, 3*nn);
//...Re fw;
	double ff = (2.*n+3.)/(2.*n-1.)*.5;
	TT[nn][nn] = -(TT[0][0] = ff*((n-1.)*n*.5-(2.*n-1.))*.5);
	if (n >= 2) TT[0][4] = TT[0][3+nn] = TT[nn][4+nn] = -(TT[nn][3] = ff);
	if (n >= 1) TT[0][2+nn*2] = -(TT[nn][1+nn*2] = -ff*n*.5);
	for (m = 1; m <= N; m++) { 
		TT[m*2][m*2-1+nn] = TT[m*2-1][m*2-1] = TT[m*2-1][m*2+nn] = -(TT[m*2][m*2] = ff*((m-1.)*(2.*n-1.)+(n-m-1.)*(n-m)*.5)*.25);
		if (m+2 <= n) TT[m*2][m*2+4] = TT[m*2][m*2+3+nn] = TT[m*2-1][m*2+4+nn] = -(TT[m*2-1][m*2+3] = ff*(m+1.)*(m+2.)*.5);
		if (m+1 <= n) TT[m*2][m*2+2+nn*2] = -(TT[m*2-1][m*2+1+nn*2] = -ff*(m+1.)*(n+m)*.5);
	}
//...Im fw;
	for (m = 1; m <= N; m++) { 
		TT[m*2+nn][m*2] = TT[m*2+nn][m*2-1+nn] = TT[m*2-1+nn][m*2+nn] = -(TT[m*2-1+nn][m*2-1] = -ff*((n-m-1.)*(n-m)*.5-(2.*n-1.))*.25);
		if (m == 2) TT[m*2+nn][0] = -(TT[m*2-1+nn][nn] = ff*(2.*n-1.+(n-3.)*(n-2.)*.25)*(n-1.)*n*.125);
		if (m == 1)	{	
			double dd = ff*((n-2.)*(n-1.)*.5+(2.*n-1.))*.25;
			TT[m*2+nn][m*2] += dd; TT[m*2+nn][m*2-1+nn] += dd; TT[m*2-1+nn][m*2-1] += dd; TT[m*2-1+nn][m*2+nn] -= dd;
			TT[m*2+nn][nn*2] = ff*(1.-n*n)*n*.25;
		}
		if (m > 2) TT[m*2+nn][m*2-5+nn] = TT[m*2-1+nn][m*2-5] = TT[m*2-1+nn][m*2-4+nn] = 
					-(TT[m*2+nn][m*2-4] = -ff*(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*(n-m+1.)*(n-m+2.)*.0625/(m-1.));
		if (m > 1) TT[m*2+nn][m*2-2+nn*2] = -(TT[m*2-1+nn][m*2-3+nn*2] = ff*(n+m)*(n-m)*(n-m+1.)*.125/m);
	}
//...fz;
	TT[nn*2][nn*2] = ff*(n-1.)*(n-1.);
	if (n >= 1) TT[nn*2][2] = TT[nn*2][1+nn] = -ff*(n-1.);
	for (m = 1; m <= N; m++) { 
		TT[m*2+nn*2][m*2+nn*2] = -(TT[m*2-1+nn*2][m*2-1+nn*2] = -ff*(n+m-1.)*(n-m-1.)*.5);
		if (m == 1)	TT[m*2+nn*2][0] = -(TT[m*2-1+nn*2][nn] = -ff*n*(2.*n-1.+(n-2.)*(n-1.)*.5)*.5);
		if (n >= 1+m) TT[m*2+nn*2][m*2+2] = TT[m*2+nn*2][m*2+1+nn] = TT[m*2-1+nn*2][m*2+2+nn] = 
						-(TT[m*2-1+nn*2][m*2+1] = ff*(m+1.)*(n-m-1.)*.5);
		if (m > 1) TT[m*2+nn*2][m*2-3+nn] = TT[m*2-1+nn*2][m*2-3] = TT[m*2-1+nn*2][m*2-2+nn] = 
					-(TT[m*2+nn*2][m*2-2] = ff*(n-m+1.)*(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*.25);
	}
	if (N == 1 && regul) {
		double G1 = 1./this->get_param(this->NUM_SHEAR);
		TT[8][0] -= 1.25*G1;
		TT[1][1] += .625*G1; TT[4][1] -= .625*G1;
		TT[2][2] += .625*G1; TT[5][2] -= .625*G1; TT[6][2] += 2.5*G1;
		TT[7][3] += 1.25*G1;
		TT[2][4] += .625*G1; TT[5][4] += .625*G1; TT[6][4] -= 2.5*G1;
		TT[1][5] += .625*G1; TT[4][5] += .625*G1;
		TT[2][6] -= 1.25*G1;
		TT[3][7] += 1.25*G1;
		TT[0][8] -= 1.25*G1;
	}

/////////////////////////////////
//...раскомплексификация строчек;
	for (l = 3*nn-1; l >= 0; l--) {				
		TT[0][l] *= 2.; TT[nn][l] *= -2.;
		for (m = 1; m <= N; m++) {
			TT[m*2][l] = 2.*(TT[m*2+nn][l]+TT[m*2][l]); 
			TT[m*2-1][l] = -2.*(TT[m*2-1+nn][l]+TT[m*2-1][l]);
			TT[m*2+nn][l] = 4.*TT[m*2+nn][l]-TT[m*2][l];
			TT[m*2-1+nn][l] = 4.*TT[m*2-1+nn][l]+TT[m*2-1][l];
			TT[m*2+nn*2][l] *=  2.; TT[m*2-1+nn*2][l] *= -2.; 
			swap(TT[m*2+nn][l], TT[m*2-1+nn][l]);
		}
	}

//////////////////////////////////////////////
//...повторяем процедуру для матрицы перехода;
	if (n > 1) {
	  	if (! TH) set_matrix(TH, 3*mm, 3*nn);	
		else		clean_matrix(TH, 3*mm, 3*nn); N -= 2;
//...Re fw;
		TH[0][0] = -(TH[mm][nn] = n*(n-1.)*.125);
		TH[0][4] = TH[0][3+nn] = TH[mm][4+nn] = -(TH[mm][3] = -.5);
		TH[0][2+nn*2] = -(TH[mm][1+nn*2] = -(n-1.)*.25);
		for (m = 1; m <= N; m++) { 
			TH[m*2][m*2-1+nn] = TH[m*2-1][m*2-1] = TH[m*2-1][m*2+nn] = -(TH[m*2][m*2] = -(n-m-1.)*(n-m)*.0625);
			TH[m*2][m*2+4] = TH[m*2][m*2+3+nn] = TH[m*2-1][m*2+4+nn] = -(TH[m*2-1][m*2+3] = -(m+1.)*(m+2.)*.25);
			TH[m*2][m*2+2+nn*2] = -(TH[m*2-1][m*2+1+nn*2] = -(m+1.)*(n-m-1.)*.25);
		}
//...Im fw;
		for (m = 1; m <= N; m++) { 
			TH[m*2+mm][m*2] = TH[m*2+mm][m*2-1+nn] = TH[m*2-1+mm][m*2+nn] = -(TH[m*2-1+mm][m*2-1] = (n-m-1.)*(n-m)*.0625);
			if (m == 2) TH[m*2+mm][0] = -(TH[m*2-1+mm][nn] = -(n-3.)*(n-2.)*(n-1.)*n*.015625);
			if (m == 1)	{	
				double dd = -(n-2.)*(n-1.)*.0625;
				TH[m*2+mm][m*2] += dd; TH[m*2+mm][m*2-1+nn] += dd; TH[m*2-1+mm][m*2-1] += dd; TH[m*2-1+mm][m*2+nn] -= dd;
				TH[m*2+mm][nn*2] = -(n-2.)*(n-1.)*n*.125;
			}
			if (m > 2) TH[m*2+mm][m*2-5+nn] = TH[m*2-1+mm][m*2-4+nn] = TH[m*2-1+mm][m*2-5] = 
						-(TH[m*2+mm][m*2-4] = (n-m-1.)*(n-m)*(n-m+1.)*(n-m+2.)*.015625/(m*(m-1.)));
			if (m > 1) TH[m*2+mm][m*2-2+nn*2] = -(TH[m*2-1+mm][m*2-3+nn*2] = (n-m-1.)*(n-m)*(n-m+1.)*.0625/m);
		}
//...fz;
		TH[mm*2][nn*2] = n*(n-1.)*.5;
		TH[mm*2][2] = TH[mm*2][1+nn] = (n-1.)*.5;
		for (m = 1; m <= N; m++) { 
			TH[m*2+mm*2][m*2+nn*2] = -(TH[m*2-1+mm*2][m*2-1+nn*2] = -(n-m-1.)*(n-m)*.25);
			if (m == 1)	TH[m*2+mm*2][0] = -(TH[m*2-1+mm*2][nn] = (n-2.)*(n-1.)*n*.125 );
			TH[m*2+mm*2][m*2+2] = TH[m*2+mm*2][m*2+1+nn] = TH[m*2-1+mm*2][m*2+2+nn] = 
									  -(TH[m*2-1+mm*2][m*2+1] = -(m+1.)*(n-m-1.)*.25);
			if (m > 1) TH[m*2+mm*2][m*2-3+nn] = TH[m*2-1+mm*2][m*2-3] = TH[m*2-1+mm*2][m*2-2+nn] = 
						-(TH[m*2+mm*2][m*2-2] = -(n-m-1.)*(n-m)*(n-m+1.)*.0625/m);
		}

/////////////////////////////////
//...раскомплексификация строчек;
		for (l = 3*nn-1; l >= 0; l--) {				
			TH[0][l] *= 2.; TH[mm][l] *= -2.;
			for (m = 1; m <= N; m++) {
				TH[m*2][l] = 2.*(TH[m*2+mm][l]+TH[m*2][l]); 
				TH[m*2-1][l] = -2.*(TH[m*2-1+mm][l]+TH[m*2-1][l]);
				TH[m*2+mm][l] = 4.*TH[m*2+mm][l]-TH[m*2][l];
				TH[m*2-1+mm][l] = 4.*TH[m*2-1+mm][l]+TH[m*2-1][l];
				TH[m*2+mm*2][l] *=  2.; TH[m*2-1+mm*2][l] *= -2.; 
				swap(TH[m*2+mm][l], TH[m*2-1+mm][l]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//...формирование внешней матрицы в аналитическом методе для граничных функций;
template <typename T>
void CHydro3D<T>::extern_matrix(T **& TT, int N)
{
	int l, m, n = N, nn = 2*n+1; 
	if (! TT)  set_matrix(TT, 3*nn, 3*nn);	
	else		clean_matrix(TT, 3*nn, 3*nn);
//...Re fw;
	TT[nn][nn] = -(TT[0][0] = (n+1.)+(n-1.)*n*.125);
	if (n >= 2) TT[0][4] = TT[0][3+nn] = TT[nn][4+nn] = -(TT[nn][3] = .5);
	if (n >= 1) TT[0][2+nn*2] = -(TT[nn][1+nn*2] = (n+1.)*.25);
	for (m = 1; m <= N; m++) { 
		TT[m*2][m*2-1+nn] = TT[m*2-1][m*2-1] = TT[m*2-1][m*2+nn] = -(TT[m*2][m*2] = (n+1.)*.5-m*.25+(n-m-1.)*(n-m)*.0625);
		if (m+2 <= n) TT[m*2][m*2+4] = TT[m*2][m*2+3+nn] = TT[m*2-1][m*2+4+nn] = -(TT[m*2-1][m*2+3] = (m+1.)*(m+2.)*.25);
		if (m+1 <= n) TT[m*2][m*2+2+nn*2] = -(TT[m*2-1][m*2+1+nn*2] = (m+1.)*(n-m+1.)*.25);
	}
//...Im fw;
	for (m = 1; m <= N; m++) { 
		TT[m*2+nn][m*2] = TT[m*2+nn][m*2-1+nn] = TT[m*2-1+nn][m*2+nn] = -(TT[m*2-1+nn][m*2-1] = -((n+1.)*.5+m*(n+.5)*.25+(n-m-1.)*(n-m)*.0625));
		if (m == 2) TT[m*2+nn][0] = -(TT[m*2-1+nn][nn] = (2.*n-1.+(n-3.)*(n-2.)*.25)*(n-1.)*n*.0625);
		if (m == 1)	{	
			double dd = (n-.5)*.25+(n-2.)*(n-1.)*.0625;
			TT[m*2+nn][m*2] += dd; TT[m*2+nn][m*2-1+nn] += dd; TT[m*2-1+nn][m*2-1] += dd; TT[m*2-1+nn][m*2+nn] -= dd;
			TT[m*2+nn][nn*2] = (2.*n+1.+(n-1.)*n*.5)*n*.25;
		}
		if (m > 2) TT[m*2+nn][m*2-5+nn] = TT[m*2-1+nn][m*2-5] = TT[m*2-1+nn][m*2-4+nn] = 
					-(TT[m*2+nn][m*2-4] = -(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*(n-m+1.)*(n-m+2.)*.03125/(m-1.));
		if (m > 1) TT[m*2+nn][m*2-2+nn*2] = -(TT[m*2-1+nn][m*2-3+nn*2] = -(2.*n+1.+(n-m)*(n-m+1.)*.5/m)*(n-m+1.)*.125);
	}
//...fz;
	TT[nn*2][nn*2] = (n+2.)*(n+2.)*.5;
	if (n >= 1) TT[nn*2][2] = TT[nn*2][1+nn] = (n+2.)*.5;
	for (m = 1; m <= N; m++) { 
		TT[m*2+nn*2][m*2+nn*2] = -(TT[m*2-1+nn*2][m*2-1+nn*2] = -((n+1.)+(n-m)*(n+m)*.25));
		if (m == 1)	TT[m*2+nn*2][0] = -(TT[m*2-1+nn*2][nn] = n*((n-1.)*(n+3.)+4.)*.125);
		if (n >= 1+m) TT[m*2+nn*2][m*2+2] = TT[m*2+nn*2][m*2+1+nn] = TT[m*2-1+nn*2][m*2+2+nn] = 
						-(TT[m*2-1+nn*2][m*2+1] = -(m+1.)*(n+m+2.)*.25);
		if (m > 1) TT[m*2+nn*2][m*2-3+nn] = TT[m*2-1+nn*2][m*2-3] = TT[m*2-1+nn*2][m*2-2+nn] = 
					-(TT[m*2+nn*2][m*2-2] = -((n-m)*(n+m+2.)+4.*m)*(n-m+1.)*.0625/m);
	}

/////////////////////////////////
//...раскомплексификация строчек;
	for (l = 3*nn-1; l >= 0; l--) {				
		TT[0][l] *= 2.; TT[nn][l] *= -2.;
		for (m = 1; m <= N; m++) {
			TT[m*2][l] = 2.*(TT[m*2+nn][l]+TT[m*2][l]); 
			TT[m*2-1][l] = -2.*(TT[m*2-1+nn][l]+TT[m*2-1][l]);
			TT[m*2+nn][l] = 4.*TT[m*2+nn][l]-TT[m*2][l];
			TT[m*2-1+nn][l] = 4.*TT[m*2-1+nn][l]+TT[m*2-1][l];
			TT[m*2+nn*2][l] *=  2.; TT[m*2-1+nn*2][l] *= -2.; 
			swap(TT[m*2+nn][l], TT[m*2-1+nn][l]);
		}
	}
}

/////////////////////////////////////////////////////
//...формироване матриц Эшелби для граничных функций;
template <typename T>
void CHydro3D<T>::set_eshelby_matrix(int N_mpl)
{
	T ** TX = NULL, * H = (T *)new_struct((2*N_mpl+1)*3*sizeof(T));
	TT = (T ***)new_struct((N_mpl+2)*sizeof(T **));
	TH = (T ***)new_struct((N_mpl+1)*sizeof(T **));

/////////////////////////////////////////////// 
//...формирование решения граничного уравнения;
	if (TT && TH && H) {
		this->solver.pivot_init((2*N_mpl+1)*3);

		int  k, j, l, num, nn;
		for (k = N_mpl; k >= 0; k--) { 
			intern_matrix(TT[k], k, TH[k-2]);
			extern_matrix(TX, k);
			if (this->solver.GaussJ(TX,  nn = (2*k+1)*3)) {
				for (num = 0; num < nn && TT[k]; num++) {
					for (           l = 0; l < nn; l++) 
					for (H[l] = 0., j = 0; j < nn; j++) H[l] += TX[l][j]*TT[k][j][num]; 
					for (           l = 0; l < nn; l++) TT[k][l][num] = H[l]/*TT[k][l][num]/(2.*k+3.)*/; 
				}
				for (num = 0; num < nn+12 && TH[k]; num++) {
					for (           l = 0; l < nn; l++) 
					for (H[l] = 0., j = 0; j < nn; j++) H[l] += TX[l][j]*TH[k][j][num]; 
					for (           l = 0; l < nn; l++) TH[k][l][num] = H[l]/*TH[k][l][num]/(2.*k+3.)*/; 
				}
			}
		}
	}
	delete_struct(H);
	delete_struct(TX);
}

///////////////////////////////////////////////////////////////////////////////
//...вычисление нормирующих множителей в представлении для уравнения Бринкмана;
template <typename T>
void CHydro3D<T>::set_normaliz_vector(int N_mpl, double rad, double kk)
{
	if (rad > 0. && kk > 0.) {
		T h1 = 1., h2 = -1./(kk*rad), h0 = -h2, hh = sqr(kk*rad);
		//T he = exp(kk*rad), hi = 1./he, h0 = 1./(kk*rad), h1 = (he+hi)*.5, h2, hh = sqr(kk*rad);
		int m;
		set_matrix(TN, 2, N_mpl+2);

///////////////////////////////////
//...prepare of the skin functions;
		for ( TN[0][m = 0] = h1/*, h2 = -(he-hi)*.5*h0*/; m <= N_mpl; TN[0][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}
		for ( TN[1][m = 0] = (h1 = 1.), h2 = h0; m <= N_mpl; TN[1][++m] = h1) {
		//for ( TN[1][m = 0] = (h1 = (he-hi)*.5), h2 = -(he+hi)*.5*h0; m <= N_mpl; TN[1][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}

////////////////////////////////                                                                                                                                              
//...пересчитываем коэффициенты;
		for ( m = 0, h1 = sqr(hh = -kk); m <= N_mpl; m++, hh *= -h1/(2.*m+1.)) { 
			h2 = TN[1][m];
			TN[1][m] =  0.5*TN[0][m]/hh;
			TN[0][m] = -0.5*h2/hh;
			//TN[1][m] =  TN[0][m]/hh;
			//TN[0][m] = -h2/hh;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//...прежний способ вычисления нормирующих множителей в представлении для уравнения Бринкмана;
template <typename T>
void CHydro3D<T>::set_normaliz_vector_old(int N_mpl, double rad, double kk)
{
	if (rad > 0. && kk > 0.) {
		T h1 = 1., h2 = -1./(kk*rad), h0 = -h2, hh = sqr(kk*rad);
		int m;
		set_matrix(TN, 2, N_mpl+2);

///////////////////////////////////
//...prepare of the skin functions;
		for ( TN[0][m = 0] = h1; m <= N_mpl; TN[0][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}
		for ( TN[1][m = 0] = (h1 = 1.), h2 = h0; m <= N_mpl; TN[1][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1 ; 
			  swap(h1, h2);
		}

////////////////////////////////
//...пересчитываем коэффициенты;
		for ( m = 0, h1 = (hh = rad)*rad; m <= N_mpl; m++, hh *= h1/(2.*m-1.)) { 
			TN[1][m] = -TN[0][m]/TN[1][m]; h2 = (TN[0][m+1]+TN[1][m]*TN[1][m+1]);
			if (to_double(h2)) TN[0][m] = hh/h2; else TN[0][m] = 0.; 
			TN[1][m] *= TN[0][m]; 
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
template <typename T>
void CHydro3D<T>::GetFuncAllValues(double X, double Y, double Z, T * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6] = { X, Y, Z, 0., 0., 1.};

/////////////////////////////////////
//...operation with all input points;
	if (0 <= i && i < this->N && this->B[i].shape && this->B[i].mp) {
		int m = this->solver.id_norm;
//////////////////////////////////////////////////////////////////
//memset(this->solver.hh[i][0][0], 0, this->solver.dim[i]*sizeof(double));
//this->B[i].shape->FULL(this->solver.hh[i][0][0], 0, 2)[4] = this->B[i].shape->get_R();	
//this->B[i].shape->set_potential(this->solver.hh[i][0][0], 0);
////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < this->solver.n; num++)
			memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(double));

		this->B[i].shape->make_local(P);
		this->B[i].shape->norm_local(P+3);
		switch (id_F) {
		case PRESSURE_VALUE: {
/////////////////////////////////////
//...calculation filtration pressure;
				this->B[i].shape->parametrization_hess(P, 1);
				jump0(P, i, 0); 

				F[0] = F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);
				F[0] += Z;
		}		break;
		case VELOCITY_VALUE: {
/////////////////////////////////////
//...calculation filtration velocity;
				this->B[i].shape->parametrization_hess(P, 1);
				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   id_variant);
				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], id_variant);
////////////////////////////
//...вычисляем сглаживатель;
//double rr = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])), r0 = this->get_param(this->NUM_GEOMT), r1 = r0 + 0.03, 
//		smooth = 1.;
//if (rr < r1) smooth = 1.-exp((rr-r0)/(rr-r1));
////////////////////////////
				F[2] = /*smooth**/this->B[i].shape->potential(this->solver.hh[i][0][m+2], id_variant);
				this->B[i].shape->norm_common_T(F);
		}		break;
		case NORMAL_R_VALUE: {
/////////////////////////////////////////////
//...normal component of filtration velocity;
			   P[5] = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
				P[3] = P[0]/P[5]; P[4] = P[1]/P[5]; P[5] = P[2]/P[5];

				this->B[i].shape->parametrization_hess(P, 1);
				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   id_variant);
				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], id_variant);
				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], id_variant);
				this->B[i].shape->norm_common_T(F); 
				F[0] = F[0]*P[3]+F[1]*P[4]+F[2]*P[5];
		}		break;
		case NORMAL_X_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[3] = 1.; P[4] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   id_variant);
				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], id_variant);
				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], id_variant);
				this->B[i].shape->norm_common_T(F);
		}   break;
		case NORMAL_Y_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[4] = 1.; P[3] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   id_variant);
				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], id_variant);
				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], id_variant);
				this->B[i].shape->norm_common_T(F);
		}   break;
		case NORMAL_Z_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[5] = 1.; P[3] = P[4] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_hess(P, 1);
				if (this->get_param(this->NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m],   id_variant);
				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m+1], id_variant);
				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m+2], id_variant);
				this->B[i].shape->norm_common_T(F);
		}   break;
		case POTENTIAL_VALUE: {
//////////////////////////////////////////
//...calculation potentials fX, fY and fZ;
				this->B[i].shape->parametrization_hess(P, 1);

				F[0]  = this->B[i].shape->potential(0, this->B[i].shape->p_cpy, id_variant);
				F[1]  = this->B[i].shape->potential(1, this->B[i].shape->p_cpy-this->B[i].shape->get_NN(), id_variant);
				F[2]  = this->B[i].shape->potential(2, this->B[i].shape->p_cpy-this->B[i].shape->get_NN()*2, id_variant);
				this->B[i].shape->norm_common_T(F);
		}		break;
		case ANALYT_VALUE: {
/////////////////////////////////////
//...flow in the cylindrical channel;
				F[0] = F[1] = 0.;
				F[2] = CylinderVelocity(sqr(P[0])+sqr(P[1]), sqr(this->get_param(this->NUM_GEOMT)), this->get_param(this->NUM_SHEAR+1));
				this->B[i].shape->norm_common_T(F);
		}		break;
		case TESTI_VALUE: {
/////////////////////////////////
//...test calculation potentials;
				this->B[i].shape->parametrization_hess(P, 1);

				F[0]  = this->B[i].shape->potential(0, this->B[i].shape->FULL(this->B[i].shape->p_cpy, 0, iparam),	id_variant);
				F[1]  = this->B[i].shape->potential(0, this->B[i].shape->FULL(this->B[i].shape->p_cpy, 0, iparam+1), id_variant);
				F[2]  = this->B[i].shape->potential(0, this->B[i].shape->FULL(this->B[i].shape->p_cpy, 0, iparam+2), id_variant);
		}		break;
		default: F[0] = i; F[1] = F[2] = 0.;
		}
	}
}

//////////////////////////////////////////////////
//...проницаемость Бринкмана в трещиноватой среде;
template <typename T>
double CHydro3D<T>::TakeLayer(double ff)
{
	double ll = ff*.5, alpha = this->get_param(this->NUM_SHEAR+1);
	if (alpha == 0.) return(2.*sqr(ll)*ll/3.);
	else {
		double kappa = sqrt(fabs(alpha))*ll, EE = exp(-2.*kappa);
		return 2.*ll*(1.-(1.-EE)/((1.+EE)*kappa))/alpha;
	}
}

////////////////////////////////////////////////////
//...проницаемость Бринкмана в цилиндрических порах;
template <typename T>
double CHydro3D<T>::TakeCylinder(double ff, double eps)
{ //...rr - квадрат радиуса;
	double rr = fabs(ff/M_PI), alpha = this->get_param(this->NUM_SHEAR+1);
	if (alpha == 0.) return(M_PI*sqr(rr)*.125);
	else {
		double kappa = fabs(alpha)*rr*.25, f1 = .5, f2 = 1., sum1 = f1, sum2 = f2;
		int k = 0, k_max = 100;
		while ((f1 > eps || f2 > eps) && k++ < k_max ) {
			sum1 += (f1 *= kappa/(k*(k+2.)));
			sum2 += (f2 *= kappa/(k*(k+0.)));
		}
		return M_PI*sqr(rr)*sum1/(4.*sum2);
	}
}

//////////////////////////////////////////////
//...течение Бринкмана в цилиндрических порах;
template <typename T>
double CHydro3D<T>::CylinderVelocity(double rr, double RR, double eps)
{ //...rr, RR - квадраты радиусов;
	double alpha = this->get_param(this->NUM_SHEAR+1);
	if (alpha == 0.) return((RR-rr)*.25);
	else {
		double kappa1 = fabs(alpha)*rr*.25, kappa2 = fabs(alpha)*RR*.25, f1 = 1., f2 = 1., sum1 = f1, sum2 = f2;
		int k = 0, k_max = 100;
		while ((f1 > eps || f2 > eps) && k++ < k_max ) {
			sum1 += (f1 *= kappa1/(k*k));
			sum2 += (f2 *= kappa2/(k*k));
		}
		return((1.-sum1/sum2)/fabs(alpha));
	}
}
#undef  Message
#endif
