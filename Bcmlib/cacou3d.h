/*==========================================*/
/*                  CACOU3D                 */
/*==========================================*/
#ifndef ___CAcou3D___
#define ___CAcou3D___

#include "ccomput3d.h"

/////////////////////////////////////////////////////
//...class of blocks partition for acoustics problem;
class CAcou3D : public CComput3D<double> {
public:
		Num_Draft type   () { return ACOU3D_DRAFT;}
		int size_of_param() { return(19);};
//...constructor;
		CAcou3D (int num_phase = 7) {
			param = (Param *)new_struct(size_of_param()*sizeof(Param));
			NUM_PHASE = num_phase;
			NUM_QUAD  = 3;
			NUM_VIBRO = 5;
		};
      virtual ~CAcou3D (void);
protected:
static int NUM_KAPPA, NUM_SHEAR, NUM_HESS;
		int  block_shape_init(Block<double> & B, int id_free);
//...auxilliary operations with block matrix;
		void jump1				(double * P, int i, int m);
		void jump2				(double * P, int i, int m);
		void jump3				(double * P, int i, int m);
		void jump_cmpl_cpy	(double * p_cpy, double * p_cmp, complex * hh, int NN, int id_memset = OK_STATE);
		void jump_conj_cpy	(int i, int m, int k, complex adm = comp(1.));
		void jump1_common_x	(double * P, int i, int m);
		void jump1_compress_x(double * P, int i, int m);
		void jump1_common_y	(double * P, int i, int m);
		void jump1_compress_y(double * P, int i, int m);
		void jump1_common_z	(double * P, int i, int m);
		void jump1_compress_z(double * P, int i, int m);
		void jump2_common_x	(double * P, int i, int m);
		void jump2_compress_x(double * P, int i, int m);
		void jump2_common_y	(double * P, int i, int m);
		void jump2_compress_y(double * P, int i, int m);
		void jump2_common_z	(double * P, int i, int m);
		void jump2_compress_z(double * P, int i, int m);
		void jump4_common_x	(double * P, int i, int m);
		void jump4_compress_x(double * P, int i, int m);
		void jump4_common_y	(double * P, int i, int m);
		void jump4_compress_y(double * P, int i, int m);
		void jump4_common_z	(double * P, int i, int m);
		void jump4_compress_z(double * P, int i, int m);
		void hessian_deriv_N (double * P, int i);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		int  gram1_phase1(CGrid * nd, int i, int id_local);
		int  gram1_phase2(CGrid * nd, int i, int id_local);
		int  gram1		  (CGrid * nd, int i, int id_local) {
			if (B[i].link[NUM_PHASE] == -1) return gram1_phase1(nd, i, id_local); else
			if (B[i].link[NUM_PHASE] == -2) return gram1_phase2(nd, i, id_local); else return(NULL_STATE);
		}
		int  gram4_phase1(CGrid * nd, int i, int id_local);
		int  gram4_phase2(CGrid * nd, int i, int id_local);
		int  gram4		  (CGrid * nd, int i, int id_local) {
			if (B[i].link[NUM_PHASE] == -1) return gram4_phase1(nd, i, id_local); else
			if (B[i].link[NUM_PHASE] == -2) return gram4_phase2(nd, i, id_local); else return(NULL_STATE);
		}
		int  transfer1_phase1(CGrid * nd, int i, int k, int id_local);
		int  transfer1_phase2(CGrid * nd, int i, int k, int id_local);
		int  transfer1			(CGrid * nd, int i, int k, int id_local) {
			if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -1) return transfer1_phase1(nd, i, k, id_local); else
			if (B[i].link[NUM_PHASE] == -2 && B[k].link[NUM_PHASE] == -2) return transfer1_phase2(nd, i, k, id_local); else return(NULL_STATE);
		}
		int  transfer2			(CGrid * nd, int i, int k, int id_local);
		int  transfer3			(CGrid * nd, int i, int k, int id_local);
		int  transfer4_phase1(CGrid * nd, int i, int k, int id_local);
		int  transfer4_phase2(CGrid * nd, int i, int k, int id_local);
		int  transfer4			(CGrid * nd, int i, int k, int id_local) {
			if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -1) return transfer4_phase1(nd, i, k, id_local); else
			if (B[i].link[NUM_PHASE] == -2 && B[k].link[NUM_PHASE] == -2) return transfer4_phase2(nd, i, k, id_local); else return(NULL_STATE);
		}
		int  counting_header (int Num);
public:
//...operations with structure and block matrix;
		void set_fasa_hmg (double kHz, double Ro0, double C0, double Ro1, double nju1, double G1);
		void set_fasa_hmg (double kHz, double Ro0, double C0);
		void set_fasa_hmg (double kHz, double Ro1, double nju1, double G1);
		void block_descrap(char * OUT_FILE);
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void GetSurferFormat(FILE * SURF, FILE * SURF1, FILE * SURF2, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2		//...degree of multipoles;
#undef  DRAFT_N									//...param(0);
#define DRAFT_Q                     0.92	//...normalization coefficient;
#undef  DRAFT_Q									//...param(1);
#define DRAFT_acselerator           50		//...degree of acselerator;
#undef  DRAFT_acselerator						//...param(2);
#define DRAFT_N_elem                10		//...number of element in facet subdivision;
#undef  DRAFT_N_elem								//...param(3);
#define DRAFT_Q_facet               1.		//...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet							//...param(4);
#define DRAFT_Hz                    0.		//...frequency (hz);
#undef  DRAFT_Hz									//...param(5);
#define DRAFT_Ro							1.225	//...media density (air);
#undef  DRAFT_Ro									//...param(6);
#define DRAFT_Velocity					340.	//...media velocity (air);
#undef  DRAFT_Velocity							//...param(7);
#define DRAFT_kappa                 0.		//...wave number;
#undef  DRAFT_kappa								//...param(8);
#define DRAFT_Ro_omega              0.		//...normalized corfficient for velocity;
#undef  DRAFT_Ro_omega							//...param(9);
#define DRAFT_C0                    0.		//...wave parameter in air;
#undef  DRAFT_C0									//...param(10);
#define DRAFT_Ro1                   9.8	//...media density (body);
#undef  DRAFT_Ro1									//...param(11);
#define DRAFT_C1                    0.		//...wave parameter in body;
#undef  DRAFT_C1									//...param(12);
#define DRAFT_G1                    1.		//...shear modulus in body;
#undef  DRAFT_G1									//...param(13);
#define DRAFT_nju1                  0.3	//...Poisson coefficient in body;
#undef  DRAFT_nju1								//...param(14);
#define DRAFT_alpha1                0.		//...coefficient of Neuber-Papkovich representation;
#undef  DRAFT_alpha1								//...param(15);
#define DRAFT_C1_2mu_lm				   0.		//...compression wave number;
#undef  DRAFT_C1_2mu_lm							//...param(16);
#define DRAFT_C1_mu						 0.	//...shear wave number;
#undef  DRAFT_C1_mu								//...param(17);
#define DRAFT_lagrange            10000.	//...Lagrange coefficient for LSM in block functional;
#undef  DRAFT_lagrange							//...param(18);
#endif
