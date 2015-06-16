/*=========================================*/
/*                 CDRAFT                  */
/*=========================================*/
#ifndef ___CDRAFT___
#define ___CDRAFT___

#include <typeinfo>
#include "dd_real.h"
#include "qd_real.h"

#include "ccells.h"
#include "csmixer.h"
#include "cgrid_el.h"
#include "cblocked.h"

////////////////////////////////////////////////////////
//...description of existing types of blocks partitions;
enum Num_Draft {
			  ERR_DRAFT = -1,
			 NULL_DRAFT =  0,
		  LAME1D_DRAFT,
		  LAME2D_DRAFT,
		  LAME3D_DRAFT,
		  HEAT2D_DRAFT,
		  HEAT3D_DRAFT,
		HEATTR2D_DRAFT,
		HEATTR3D_DRAFT,
		 HYDRO3D_DRAFT,
		 COHES2D_DRAFT,
		 COHES3D_DRAFT,
		 MINDL2D_DRAFT,
		 MINDL3D_DRAFT,
		 VISCO2D_DRAFT,
		 VISCO3D_DRAFT,
		 FRACT2D_DRAFT,
  VISCO2D_GRAD_DRAFT,
/////////////
//...стержни;
		  BARS25_DRAFT,
		  BEAM3D_DRAFT,
/////////////
//...аустика;
		  ACOU3D_DRAFT,
		 VIBRO3D_DRAFT,
///////////////////////////////////////
//...оборудование вместе со структурой;
		  WIRE2D_DRAFT,
		 CABLE3D_DRAFT,
		 TOWER3D_DRAFT,
		 FITTING_DRAFT,
//////////////////////////////////
//...отображение исходной области;
		 IMAPI2D_DRAFT,
		  MAPI2D_DRAFT,
			 NUMS_DRAFT
};
//////////////////////////////////////////////////
//...разменость задачи (метрическое пространство);
inline int draft_dim(Num_Draft type) {
	if (type ==			  NULL_DRAFT)	return(0); else
	if (type ==		   LAME1D_DRAFT ||
		 type ==		   HEAT2D_DRAFT ||
		 type ==	    HEATTR2D_DRAFT ||
		 type ==		   HEAT3D_DRAFT ||
		 type ==	    HEATTR3D_DRAFT ||
		 type ==		  TOWER3D_DRAFT ||
		 type ==		  CABLE3D_DRAFT ||
		 type ==		  FITTING_DRAFT)	return(1); else
	if (type ==		   LAME2D_DRAFT ||
		 type ==		  VISCO2D_DRAFT ||
		 type ==		  COHES2D_DRAFT ||
		 type ==		  MINDL2D_DRAFT ||
		 type ==		  FRACT2D_DRAFT ||
		 type ==		  IMAPI2D_DRAFT ||
		 type ==		   MAPI2D_DRAFT ||
		 type ==		   BARS25_DRAFT ||
		 type ==	VISCO2D_GRAD_DRAFT ||
		 type ==		   WIRE2D_DRAFT)	return(2); else
	if (type ==		   ACOU3D_DRAFT ||
		 type ==		  VIBRO3D_DRAFT ||
		 type ==		   LAME3D_DRAFT ||
		 type ==		  VISCO3D_DRAFT ||
		 type ==		  HYDRO3D_DRAFT ||
		 type ==		  COHES3D_DRAFT ||
		 type ==		  MINDL3D_DRAFT ||
		 type ==		   BEAM3D_DRAFT)	return(3); else
	return(ERR_DIM);
}
///////////////////////////////////////////////////////
//...число компонент в полях (при выводе через SURFER);
inline int surfer_dim(Num_Draft type) {
	if (type ==			  NULL_DRAFT)	return(0); else
	if (type ==		   LAME1D_DRAFT ||
		 type ==		  FITTING_DRAFT)	return(1); else
	if (type ==		   LAME2D_DRAFT ||
		 type ==		  VISCO2D_DRAFT ||
		 type ==		  COHES2D_DRAFT ||
		 type ==		  MINDL2D_DRAFT ||
		 type ==		  FRACT2D_DRAFT ||
		 type ==		   HEAT2D_DRAFT ||
		 type ==		 HEATTR2D_DRAFT ||
		 type ==		   BARS25_DRAFT ||
		 type ==		   WIRE2D_DRAFT ||
		 type ==		  IMAPI2D_DRAFT ||
		 type == VISCO2D_GRAD_DRAFT ||
		 type ==		   MAPI2D_DRAFT)	return(2); else
	if (type ==		   ACOU3D_DRAFT ||
		 type ==		  VIBRO3D_DRAFT ||
		 type ==		   LAME3D_DRAFT ||
		 type ==		  VISCO3D_DRAFT ||
		 type ==		  HYDRO3D_DRAFT ||
		 type ==		  COHES3D_DRAFT ||
		 type ==		  MINDL3D_DRAFT ||
		 type ==		   HEAT3D_DRAFT ||
		 type ==		 HEATTR3D_DRAFT ||
		 type ==		   BEAM3D_DRAFT ||
		 type ==		  CABLE3D_DRAFT ||
		 type ==		  TOWER3D_DRAFT)	return(3); else
	return(ERR_DIM);
}

///////////////////////////////////////////////
//...кодировка используемых в программе блоков;
enum Num_Block {
			 NULL_BLOCK = 0x00000000,
			 POLY_BLOCK = 0x00000010,
		  POLY1D_BLOCK = 0x00001010,
		  POLY2D_BLOCK = 0x00002010,
		  POLY3D_BLOCK = 0x00003010,
			 IPOL_BLOCK = 0x00000011,
		  IPOL1D_BLOCK = 0x00001011,
		  IPOL2D_BLOCK = 0x00002011,
		  IPOL3D_BLOCK = 0x00003011,
			 ELLI_BLOCK = 0x00000020,
		  ELLI1D_BLOCK = 0x00001020,
		  ELLI2D_BLOCK = 0x00002020,
		  ELLI3D_BLOCK = 0x00003020,
			 ZOOM_BLOCK = 0x00000030,
		  STREEP_BLOCK = 0x00000040,
		  CIRCLE_BLOCK = 0x00000050,
		  SPHERE_BLOCK = 0x00000060,
		  CLAYER_BLOCK = 0x00000070,
			 BEAM_BLOCK = 0x00000080,
			 WAVE_BLOCK = 0x00000090,
			 ESHE_BLOCK = 0x000000A0,
	  ESHE_ZOOM_BLOCK = 0x000000B0,
//////////////////////////////////////////
//...список кодов основных блоков для BARs;
		  RECTAN_BLOCK = 0x00010000,
			 ROLL_BLOCK = 0x00020000,
		  UGOLOK_BLOCK = 0x00030000,
		  SREZKA_BLOCK = 0x00040000,
			BULBA_BLOCK = 0x00050000,
		  BORDER_BLOCK = 0x00060000,
//////////////////////////////////////////
//...список дополнительных кодов для BARs;
			 TRIA_BLOCK = 0x00070000,
			QUADR_BLOCK = 0x00080000,
			 RING_BLOCK = 0x00090000,
			 GEAR_BLOCK = 0x00100000,
			SHAFT_BLOCK = 0x00110000,
			 HEXA_BLOCK = 0x00120000,
////////////////////////////////////////////////////////////
//...список кодов для элементов усложненной библиотеки BARs;
			 BASE_BLOCK = 0x00130000,
			  LEG_BLOCK = 0x00140000,
//////////////////////////////////////////////////////
//...коды блоков для 2D-, осесимметричных и 3D- задач;
			 DIM3_BLOCK = 0x00150000,
			 AXIS_BLOCK = 0x00160000,
			 DIM2_BLOCK = 0x00170000,
/////////////////////////////////////////////////////
//...коды блоков для моделирования блочной структуры;
	 SUB_UGOLOK_BLOCK = 0x00180000,
	 SUB_UGOLOK2BLOCK = 0x00190000,
	 SUB_SREZKA_BLOCK = 0x00200000,
		SUB_QUAD_BLOCK = 0x00210000,
///////////////////////////////////////////////////
//...коды блоков для конкретных акустических задач;
	 L_TYPE_BOX_BLOCK = 0x00220000,
	L_TYPE_ROOM_BLOCK = 0x00230000,
////////////////////////////////////////////
//...коды блоков для задач теории упругости;
	 INTR_CONUS_BLOCK = 0x00240000,
/////////////////////////////////////////////////////////////////
//...нулевой блок для работы и оптимизации и аналитические блоки;
			 ZERO_BLOCK = 0x00250000,
	 RECTAN_NEW_BLOCK = 0x00260000
};

//////////////////////////////////////////////////
//...description of the state of computing kernel;
enum Num_Comput {
			 BASIC_COMPUT = 1000,
		  E_BASIC_COMPUT,
			 BASICtCOMPUT,
			PERIOD_COMPUT,
		 E_PERIOD_COMPUT,
			VOLUME_COMPUT,
		  MAPPING_COMPUT,
		  MAPPINGtCOMPUT,
		  ESHELBY_COMPUT,
		  SPECIAL_COMPUT,
		  SPECIAL2COMPUT,
		  SPECIAL3COMPUT,
		 COVERING_COMPUT,
};

/////////////////////////////
//...type of solving problem;
enum Num_Solving {
	 SQUARE_SOLVING =  0, //...решение общей задачи методом наименьших квадратов;
  PERIODIC_SOLVING =  1, //...решение периодической задачи методом наименьших квадратов;
   SPECIAL_SOLVING =  2, //...специальный метод решения задачи методом наименьших квадратов;
	 ENERGY_SOLVING = 10, //...решение общей задачи энергетческим методом;
E_PERIODIC_SOLVING = 11, //...решение периодической задачи энергетческим методом;
 E_SPECIAL_SOLVING = 12, //...специальный метод решения задачи энергетическим методом;
};

///////////////////////////////
//...type of block integration;
enum Num_Volume {
		  ERR_VOLUME = -1, //...интегрирование отсутствует;
		BOUND_VOLUME =  0, //...интегрирование только по границе блоков;
		BLOCK_VOLUME =  1, //...интегрирование только по объему блоков;
		MIXED_VOLUME =  2, //...cмешанное интегрирование по границе и по объему блоков;
};

//////////////////////////////
//...type of resulting values;
enum Num_Value {	
										ERR_VALUE = -1, //...результаты отсутствуют;
										SPL_VALUE,
									  HEAT_VALUE,
									  FLUX_VALUE,
							 FLUX_COMPOS_VALUE,
							 HEAT_ANALYT_VALUE,
							 FLUX_ANALYT_VALUE,
									  LINK_VALUE,
									 DISPL_VALUE,
									 DILAT_VALUE,
									 CROSS_VALUE,
									 FRACT_VALUE,
									 TESTI_VALUE,
									ENERGY_VALUE,
									MOMENT_VALUE,
								 PRESSURE_VALUE,
								 VELOCITY_VALUE,
								 ROTATION_VALUE,
								POTENTIAL_VALUE,
								PROCESSOR_VALUE,
								PUREFRACT_VALUE,
									FLUX_R_VALUE,
									FLUX_X_VALUE,
									FLUX_Y_VALUE,
									FLUX_Z_VALUE,
									THIN_R_VALUE,
									THIN_X_VALUE,
									THIN_Y_VALUE,
									THIN_Z_VALUE,
								 MOMENT_R_VALUE,
								 MOMENT_X_VALUE,
								 MOMENT_Y_VALUE,
								 MOMENT_Z_VALUE,
								 STRESS_R_VALUE,
								 STRESS_X_VALUE,
								 STRESS_Y_VALUE,
								 STRESS_Z_VALUE,
								 NORMAL_R_VALUE,
								 NORMAL_X_VALUE,
								 NORMAL_Y_VALUE,
								 NORMAL_Z_VALUE,
								 STRAIN_R_VALUE,
								 STRAIN_X_VALUE,
								 STRAIN_Y_VALUE,
								 STRAIN_Z_VALUE,
								HEAT_ESHE_VALUE,
								FLUX_ESHE_VALUE,
							  DISPL_JUMP_VALUE,
							  DISPL_UNIX_VALUE,
							  DISPL_INDI_VALUE,
							  DISPL_ESHE_VALUE,
						  DISPL_CLASSIC_VALUE,
						 DISPL_COHESION_VALUE,
						  STRESS_R_JUMP_VALUE,
						  STRESS_X_JUMP_VALUE,
						  STRESS_Y_JUMP_VALUE,
						  STRESS_Z_JUMP_VALUE,
						  STRESS_R_UNIX_VALUE,
						  STRESS_X_UNIX_VALUE,
						  STRESS_Y_UNIX_VALUE,
						  STRESS_Z_UNIX_VALUE,
						  STRESS_X_ESHE_VALUE,
						  STRESS_Y_ESHE_VALUE,
						  STRESS_Z_ESHE_VALUE,
						  STRESS_R_INDI_VALUE,
						  NORMAL_R_JUMP_VALUE,
						  NORMAL_X_JUMP_VALUE,
						  NORMAL_Y_JUMP_VALUE,
						  NORMAL_Z_JUMP_VALUE,
					  STRESS_R_CLASSIC_VALUE,
					  STRESS_X_CLASSIC_VALUE,
					  STRESS_Y_CLASSIC_VALUE,
					  STRESS_Z_CLASSIC_VALUE,
					 STRESS_R_COHESION_VALUE,
					 STRESS_X_COHESION_VALUE,
					 STRESS_Y_COHESION_VALUE,
					 STRESS_Z_COHESION_VALUE,
									ANALYT_VALUE = 1000,
									ANALYT2VALUE,
									ANALYT3VALUE,
									ANALYT4VALUE,
									ANALYT5VALUE,
									ANALYT6VALUE,
									ANALYT7VALUE,
									ANALYT8VALUE,
									ANALYT9VALUE,
							TESTI_ANALYT_VALUE,
					  SPL_ANALYT_RIGID_VALUE,
					PRESS_ANALYT_RIGID_VALUE,
				  FLUX_X_ANALYT_RIGID_VALUE,
				  FLUX_Z_ANALYT_RIGID_VALUE,
					 SPL_ANALYT_ABSORB_VALUE,
				  PRESS_ANALYT_ABSORB_VALUE,
				 FLUX_Z_ANALYT_ABSORB_VALUE,
	 HEAT_HOMOG_ALONG_LAYER_ANALYT_VALUE,
	 HEAT_HOMOG_ACROSSLAYER_ANALYT_VALUE,
	 FLUX_HOMOG_ALONG_LAYER_ANALYT_VALUE,
	 FLUX_HOMOG_ACROSSLAYER_ANALYT_VALUE,
};

////////////////////////////
//...structure of one block;
template <typename T>
struct Block {
		Num_Block		type;		//...type of element of blocks structure;
		Block				* B;     //...description of all blocks structure;
		CCells			* bar;   //...geometry of blocks element;
		Topo				* link;  //...links of the blocks element;
		CShapeMixer<T> * shape; //...blocks multipoles;
		CMap				* mp;    //...geometrical mapping of blocks multipoles;
};

////////////////////////////////////////////////////////////////////////////////////
//...циклы для перебора блоков при заполнения всей матрицы или одной блочной строки;
#define LOOP_OPT(loop, opt, k)\
for (this->solver.IL ? loop = 0 : loop = -(opt = -this->N);\
	 (opt >= 0 && loop < this->solver.IL[this->solver.p[opt]][0] && (k = this->solver.IL[this->solver.p[opt]][loop+this->solver.JR_SHIFT]) < this->N || \
	  opt  < 0 && loop < this->N && (k = this->solver.p[loop]) < this->N); loop++)

#define LOOP_MASK(opt, k)\
if (opt < 0)				this->solver.clean_mode (NO_TR| NO_TL); else \
if (k == this->solver.p[opt]) this->solver.change_mode(NO_TR, NO_TL); else \
								this->solver.change_mode(NO_TL, NO_TR)
///////////////////////////////////////////
//...циклы для перебора фасет одного блока;
#define IF_FACET(bar, num, m, ANY)\
if  (bar->ce[num]->cells_dim() == 2 && bar->ce[num]->mp && bar->ce[num]->graph && \
	  bar->ce[num]->mp[m = size_of_map (bar->ce[num]->mp)] == FACET_CELL && \
	  bar->ce[num]->graph[1] > 2 && (ANY || bar->ce[num]->mp[m+1] > 0.))

#define IF_LINK_FACET(bar, num, m, k_ext)\
if  (bar->ce[num]->cells_dim() == 2 && bar->ce[num]->mp && bar->ce[num]->graph && \
	  bar->ce[num]->mp[m = size_of_map (bar->ce[num]->mp)] == FACET_CELL && \
	  bar->ce[num]->mp[m+1] == -k_ext-1. && bar->ce[num]->graph[1] > 2)

#define IF_ANY_FACET(bar, num)\
if  (bar->ce[num]->cells_dim() == 2 && bar->ce[num]->mp && bar->ce[num]->graph && \
	  bar->ce[num]->mp[size_of_map(bar->ce[num]->mp)] == (CMap)FACET_CELL && \
	  bar->ce[num]->graph[1] > 2)

#define LOOP_FACET(bar, num, m)\
for (num = 0; num < bar->graph[0]; num++) \
IF_FACET(bar, num, m, 0)

#define LOOP_LINK_FACET(bar, num, m, k_ext)\
for (num = 0; num < bar->graph[0]; num++) \
IF_LINK_FACET(bar, num, m, k_ext)

/////////////////////////////////////
//...basic class of blocks partition;
template <typename T>
class CDraft : public CBase {
public:
		int		  N;	  //...blocks number;
		Block<T> * B;    //...description block structure;
		CCells	* bar;  //...description boundary surface;
		Topo		* link; //...description types of the boundary components;
public:
static int NUM_MPLS, NUM_QUAD, NUM_TIME, NUM_ADHES, NUM_VIBRO, NUM_GEOMT, NUM_PHASE, MAX_PHASE, BOX_LINK_PERIOD;
		CBlockedSolver<T> solver;
public:
		CGrid_el stru;
		int	 * id_prop; //...types of boundary conditions;
		double * pp_cond; //...description of the boundary parameters;
      double * src;	   //...description of the point sources;
		void bnd_marker2D(Topo j, double * pp) {
			if (pp) switch (j) { //...задаем граничные условия для помеченных границ;
				case BOUND1_STATE: pp[0] = 1.; pp[1] = 0.; pp[2] = 1.; break; //...растяжение единичной силой;
				case BOUND2_STATE: pp[0] = 0.; pp[1] = 0.; pp[2] = 1.; break; //...свободная граница;
				case BOUND3_STATE: pp[0] = 0.; pp[1] = 0.; pp[2] = MAX_HIT; break; //...жесткое закрепление;
				case BOUND4_STATE: pp[0] = 1.; pp[1] = 0.; pp[2] = MAX_HIT; break; //...жесткое растяжение;
			}
		}
		void bnd_marker3D(Topo j, double * pp) {
			if (pp) switch (j) { //...задаем граничные условия для помеченных границ;
				case BOUND1_STATE: pp[0] = -1.; pp[1] = 0.; pp[2] =  0.; pp[3] = 1; break; //...растяжение единичной силой;
				case BOUND2_STATE: pp[0] =  0.; pp[1] = 0.; pp[2] =  0.; pp[3] = 1; break; //...свободная граница;
				case BOUND3_STATE: pp[0] =  0.; pp[1] = 0.; pp[2] =  0.; pp[3] = MAX_HIT;  break; //...жесткое закрепление;
				case BOUND4_STATE: pp[0] =  0.; pp[1] = 0.; pp[2] =  0.; pp[3] = MAX_HIT;  break; //...жесткое растяжение;
//...для задач акустики;
				case BOUND5_STATE: pp[0] =  0.; pp[1] = 0.; pp[2] = -1.; pp[3] = 0.; break; //...излучение;
				case BOUND6_STATE: pp[0] =  0.; pp[1] = 0.; pp[2] =  0.; pp[3] = 0.; break; //...отражение;
				case BOUND7_STATE: pp[0] =  1.; pp[1] = 0.; pp[2] =  0.; pp[3] = 0.; break; //...поглощение;
			}
		}
protected:
		int solv, volm, buf_count, N_buf; //...type of solving problem, type of block integration and bufferization;
		virtual void comput1(int opt){}
		virtual void comput2(int opt){}
		virtual void comput3(int opt){}
		virtual void comput4(int opt){}
		virtual void comput5(int opt, T * K, Num_Value _FMF, int id_variant){}
		virtual void comput6(int opt, T * K, Num_Value _FMF, int id_variant){}
		virtual void comput7(int opt, T * K){}
		virtual Num_State comput_kernel1(Num_Comput Num) { return NULL_STATE;}
		virtual Num_State comput_kernel2(Num_Comput Num) { return NULL_STATE;}
		virtual Num_State comput_kernel3()	{ return NULL_STATE;}
		virtual Num_State comput_kernel4()	{ return NULL_STATE;}
		virtual Num_State comput_kernel5()	{ return NULL_STATE;}
		virtual Num_State computing_header(Num_Comput Num){ return NULL_STATE;}
		virtual int block_shape_init(Block<T> & B, Num_State id_free) { return 0;}
//...инициализация и линкование блоков;
		virtual void blocks_release();
		virtual int  block_geom_link (Block<T> & B, int * geom, int * geom_ptr, int max_env, int * block_env, int N_buf);
		virtual int	 block_iddir_2D  (Block<T> & B, int i, int j, double * par, double eps = 0.0005);
		virtual int  block_iddir_3D  (Block<T> & B, int i, double * par, int id_fast = OK_STATE, double eps = 0.0005);
		virtual int	 block_plink_2D  (Block<T> & B, int & i, int & j, int & id_dir, double * par, double eps = 0.0005);
		virtual int	 block_plink_3D  (Block<T> & B, int & i, int & id_dir, double * par, int id_fast = OK_STATE, double eps = 0.0005);
		virtual int  geom_iddir_2D	  (Block<T> & B, int i, double * par, double eps = 0.0005);
		virtual int  geom_iddir_3D	  (Block<T> & B, int i, double * par, double eps = 0.0005);
		virtual int	 geom_plink_2D   (Block<T> & B, int & i, int & id_dir, double * par, double eps = 0.0005);
		virtual int	 geom_plink_3D   (Block<T> & B, int & i, int & id_dir, double * par, double eps = 0.0005);
		virtual void bar_activate	  (Block<T> & B, int id_status = OK_STATE, int id_fast = OK_STATE);
//...заполнение блочной матрицы;
		virtual Num_State gram1		(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State gram2		(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State gram2peri(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State gram3		(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State gram4		(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State gram5		(CGrid * nd, int i, int id_local) { return NULL_STATE;}
		virtual Num_State rigidy1	(CGrid * nd, int i, T * K) { return NULL_STATE;}
		virtual Num_State rigidy2	(CGrid * nd, int i, T * K) { return NULL_STATE;}
		virtual Num_State rigidy3	(CGrid * nd, int i, T * K) { return NULL_STATE;}
		virtual Num_State rigidy4	(CGrid * nd, int i, T * K) { return NULL_STATE;}
		virtual Num_State rigidy5	(CGrid * nd, int i, T * K) { return NULL_STATE;}
		virtual void energy1			(CGrid * nd, int i, T * energy, int id_variant){}
		virtual void energy2			(CGrid * nd, int i, T * energy, int id_variant){}
		virtual void energy3			(CGrid * nd, int i, T * energy, int id_variant){}
		virtual void energy4			(CGrid * nd, int i, T * energy, int id_variant){}
		virtual Num_State transfer1(CGrid * nd, int i, int k, int id_local) { return NULL_STATE;}
		virtual Num_State transfer2(CGrid * nd, int i, int k, int id_local) { return NULL_STATE;}
		virtual Num_State trans_esh(CGrid * nd, int i, int k, int id_local) { return NULL_STATE;}
		virtual Num_State transfer3(CGrid * nd, int i, int k, int id_local) { return NULL_STATE;}
		virtual Num_State transfer4(CGrid * nd, int i, int k, int id_local) { return NULL_STATE;}
public:
//...интерфейс вычислительных схем;
		void computing(int opt = -1, Num_Comput Num = BASIC_COMPUT) {
				if (Num ==	 BASIC_COMPUT || Num ==   BASICtCOMPUT) comput1(opt); else
				if (Num == MAPPING_COMPUT || Num == MAPPINGtCOMPUT) comput2(opt); else
				if (Num ==	PERIOD_COMPUT) comput3(opt); else
				if (Num == ESHELBY_COMPUT) comput4(opt);
		}
		void GetRigidy(T * K, int opt = -1, Num_Comput Num = BASIC_COMPUT, Num_Value _FMF = ERR_VALUE, int id_variant = 0) {
				if (Num ==    BASIC_COMPUT) comput5(opt, K, _FMF, id_variant); else
				if (Num ==   VOLUME_COMPUT) comput6(opt, K, _FMF, id_variant);
				if (Num == COVERING_COMPUT) comput7(opt, K);
		}
		Num_State computing_kernel(Num_Comput Num = BASIC_COMPUT) {
				if (Num ==	 BASIC_COMPUT || Num == MAPPING_COMPUT || Num == PERIOD_COMPUT || Num == ESHELBY_COMPUT) return comput_kernel1(Num); else
				if (Num ==	 BASICtCOMPUT || Num == MAPPINGtCOMPUT) return comput_kernel2(Num); else
				if (Num == SPECIAL_COMPUT) return comput_kernel3(); else//...дополнительные схемы;
				if (Num == SPECIAL2COMPUT) return comput_kernel4(); else
				if (Num == SPECIAL3COMPUT) return comput_kernel5(); else return NULL_STATE;
		}
public:
//...интерфейс заполнения блочной матрицы;
		Num_State GramAll(CGrid * nd, int i, int id_local = 0, Num_Comput Num = BASIC_COMPUT) {
				if (Num ==    BASIC_COMPUT) return gram1(nd, i, id_local); else
				if (Num ==  E_BASIC_COMPUT) return gram4(nd, i, id_local); else
				if (Num ==	 PERIOD_COMPUT) return gram2(nd, i, id_local); else
				if (Num == E_PERIOD_COMPUT) return gram3(nd, i, id_local); else
				if (Num == 	  BASICtCOMPUT) return gram5(nd, i, id_local); else return NULL_STATE;
		}
		Num_State TransferBB(CGrid * nd, int i, int k, int id_local = 0, Num_Comput Num = BASIC_COMPUT) {
				if (Num == 	  BASIC_COMPUT) return transfer1(nd, i, k, id_local); else
				if (Num ==  E_BASIC_COMPUT) return transfer4(nd, i, k, id_local); else
				if (Num ==	 PERIOD_COMPUT) return transfer2(nd, i, k, id_local); else
				if (Num == E_PERIOD_COMPUT) return transfer3(nd, i, k, id_local); else return NULL_STATE;
		}
		Num_State RigidyAll(CGrid * nd, int i, T * K, Num_Comput Num = BASIC_COMPUT) {
				if (Num ==	  BASIC_COMPUT) return rigidy1(nd, i, K); else
				if (Num ==   VOLUME_COMPUT) return rigidy2(nd, i, K); else
				if (Num ==  SPECIAL_COMPUT) return rigidy3(nd, i, K); else
				if (Num ==  SPECIAL2COMPUT) return rigidy4(nd, i, K); else
				if (Num == COVERING_COMPUT) return rigidy5(nd, i, K); else return NULL_STATE;
		}
		void EnergyAll(CGrid * nd, int i, T * energy, int id_variant = 0, Num_Comput Num = BASIC_COMPUT) {
				if (Num ==	 BASIC_COMPUT) energy1(nd, i, energy, id_variant); else
				if (Num ==  VOLUME_COMPUT) energy2(nd, i, energy, id_variant); else
				if (Num == SPECIAL_COMPUT) energy3(nd, i, energy, id_variant); else
				if (Num == SPECIAL2COMPUT) energy4(nd, i, energy, id_variant);
		}
public:
//...заполнение параметров блочной структуры;
		virtual void set_fasa_hmg(double CC){}
		virtual void set_fasa_hmg(double CC, double CC2){}
		virtual void set_fasa_hmg(double Hz, double Ro1, double C1){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2){}
		virtual void set_fasa_hmg(double Hz, double Ro1, double Ro2, double C1, double C2){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B){}
		virtual void set_fasa_hmg(double R1, double R2, double nju1, double nju2, double nju3, double G1, double G2, double G3, double alpha = -1.){}
		virtual void set_fasa_hmg(double K1, double K2, double G1, double G2, double l1, double l2, double d1, double d2, double A, double B, double d0){}
		virtual void set_time	 (double d_time = 0.1, double end_time = 1.) {}
		//	set_param(NUM_TIME+1, 0.); set_param(NUM_TIME+2, d_time); set_param(NUM_TIME+3, end_time);
		//}
		virtual void set_adhesion(double AA, double BB){}
		virtual void set_geometry(double rad, double layer = 0.){}
		virtual void set_mpls    (double mpls) {set_param(NUM_MPLS, mpls);}
		virtual void set_quad	 (double quad) {set_param(NUM_QUAD, quad);}
		virtual void set_vibro   (double Hz)   {set_param(NUM_VIBRO,  Hz);}
		virtual void set_lattice (double lattice)  {set_param(NUM_TIME, lattice);}
		virtual void set_normaliz(double normaliz) {set_param(NUM_MPLS+1,		normaliz);}
		virtual void set_lagrange(double lagrange) {set_param(size_of_param()-1, lagrange);}
//////////////////////////////
//...параметризация в образце;
		virtual Param * ID_EPS()        { return(param);}
		virtual Param * ID_CUSTOMER()   { return(param);}
		virtual Param * ID_MAIN()       { return(param);}
		virtual Param * ID_YOUNG_MODUL(){ return(param);}
		virtual Param * ID_INERTIA_MOM(){ return(param);}
		virtual Param * ID_TORSION_MOM(){ return(param);}
		virtual Param * ID_ENERGY_MOM() { return(param);}
public:
//...интерфейс результатов решения задачи;
		virtual void GetFuncAllValues(double X, double Y, double Z, T * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0){}
		virtual void GetFuncAllValues(double * P, T * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0) { GetFuncAllValues(P[0], P[1], P[2], F, id_block, id_F, id_variant, iparam);}
//... формат выдачи Golden Surfer;
		void GetSurferFormat(char * SURF_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0);
		void GetSurferFormat(const char * SURF_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0);
		void GetDataFormat(char * DATA_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0);
		void GetDataFormat(const char * DATA_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0);
//... определение функций;
		virtual void GetSurferFormat(FILE * SURF, FILE * SURF1, FILE * SURF2, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam);
		virtual void GetDataFormat  (FILE * DATA, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam);
//... формат выдачи результатов на блочной структуре в виде CSV таблиц;
		void GetCsvFormat(char * CSV_FILE, CGrid * nd, int id_variant = 0, int id_centroid = 0, CGrid * bnd = NULL);
		void GetCsvFormat(const char * CSV_FILE, CGrid * nd, int id_variant = 0, int id_centroid = 0, CGrid * bnd = NULL);
//... определение функций;
		virtual void GetCsvFormat (FILE * CSV, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd);
//...аналитические модели;
		virtual double TakeLayer				(double ff)						  { return 0.;}
		virtual double TakeLayer_E1			(double ff)						  { return 0.;}
		virtual double TakeLayer_E2			(double ff)						  { return 0.;}
		virtual double TakeLayer_G1			(double ff)						  { return 0.;}
		virtual double TakeLayer_G2			(double ff)						  { return 0.;}
		virtual double TakeCylinder			(double ff, double eps = EE) { return 0.;}
		virtual double TakeEshelby				(double ff, double ff_l)	  { return 0.;}
		virtual double TakeEshelby_two		(double ff)						  { return 0.;}
		virtual double TakeEshelby_grad		(double ff)						  { return 0.;}
		virtual double TakeEshelby_volm		(double ff, double ff_l)	  { return 0.;}
		virtual double TakeEshelby_volm_two (double ff)						  { return 0.;}
		virtual double TakeEshelby_volm_sym (double ff)						  { return 0.;}
		virtual double TakeEshelby_shear		(double ff, double ff_l, double eps = EE, int max_iter = 100) { return 0.;}
		virtual double TakeEshelby_shear_two(double ff, double eps = EE, int max_iter = 100)				  { return 0.;}
		virtual double TakeEshelby_shear_sym(double ff, double eps = EE, int max_iter = 100)				  { return 0.;}
		virtual double TakeEshelby_shear_det(double ff, double alpha = -1)										  { return 0.;}
		virtual double TakeLayer_kk(int N, double * ff, double * kk)												  { return 0.;}
		virtual double TakeLayer_kk(int N, double * ff, double * kk, double * mu)								  { return 0.;}
		virtual double TakeLayer_kk(int N, double * ff, double * kk, double * mu, double * nj)				  { return 0.;}
		virtual double TakeLayer_GG(int N, double * ff, double * mu, double * nj)								  { return 0.;}
		virtual double TakeLayer_GG(int N, double * ff, double * kk, double * mu, double eps = EE, int max_iter = 100) { return 0.;}
		virtual double TakeLayer_GG(int N, double * ff, double * kp, double * mu, double * nj, double eps = EE, int max_iter = 100) { return 0.;}
//...функция пересчета полей (температур) согласно схеме установления;
		virtual void TakeStabStep(double * Temp, int NN, double alpha){};
		virtual void TakeStabStep_layer(double * Temp, int N_SC, int N_CU, int N_cells, double alpha){};
public:
		virtual Num_Draft type() { return NULL_DRAFT;}
public:
//...constructor/destructor;
	   double Acou3D_homotat, Acou3D_L_inf; int Acou3D_N_elem2;
		CDraft() {
			N		= 0;
			B     = NULL;
			bar   = NULL;
			link  = NULL;
			id_prop = NULL;
			pp_cond = src = NULL;
			volm	= BOUND_VOLUME;
			solv  = SQUARE_SOLVING;
			N_buf = 20; buf_count = 0;
			Acou3D_homotat = Acou3D_L_inf = 0.;
			Acou3D_N_elem2 = 0;
		};
		virtual ~CDraft (void) {
			blocks_release ();
			if (bar) {
				 bar->zero_cells(); 
				 delete bar; bar = NULL;
			}
			delete_struct(link);
			delete_struct(id_prop);
			delete_struct(pp_cond);
			delete_struct(src);
		}
public:
//...линкование блоков (и степени свободы);
		int set_link(Block<T>& B, int m)	{
			delete_struct(B.link);
			if ( NULL != (B.link = (Topo *)new_struct((m+1)*sizeof(Topo)))) {
				memset(B.link, ERR_STATE, (m+1)*sizeof(Topo)); B.link[0] = m;
				return(1);
			}
			return(0);
		}
		int link_add(Block<T> & B, Topo new_element, int id_element = OK_STATE) {
			if (! id_element || ! link_id(B, new_element)) {
			  if (link_pad(B, 1)) B.link[B.link[0]] = new_element;
			  else                return(0);
			}
			return(1);
		}
		int link_id(Block<T> & B, int some_element) {
			if (B.link) {
				int j = B.link[0];

				while (0 < j && some_element != B.link[j]) j--;
				return(0 < j);
			}
			return(0);
		}
		int link_pad(Block<T> & B, int N_pad) {
			if (N_pad > 0) {
				int m = (int)(B.link ? B.link[0] : 0)+1;
				Topo * new_link = (Topo *)new_struct((m+N_pad)*sizeof(Topo));
				if  (! new_link) return(0);

				memset (new_link, ERR_STATE, (m+N_pad)*sizeof(Topo));
				if (B.link) memcpy(new_link, B.link, m*sizeof(Topo)); else new_link[0] = 0;

				delete_struct(B.link); B.link = new_link; B.link[0] += N_pad;
			}
			return(1);
		}
		int freedom_block(int k) { return(k < N && B[k].shape ? B[k].shape->NN : 0);}
public:
//...инициализация блоков;
		void change_solv(int set_solv = SQUARE_SOLVING, int * save_solv = NULL) { if (save_solv) *save_solv = volm; solv = set_solv;}
		void change_volm(int set_volm = BOUND_VOLUME,   int * save_volm = NULL) { if (save_volm) *save_volm = volm; volm = set_volm;}
		void shapes_init  (Num_State id_free);
		void SetBUniStruct(Num_Block type = POLY_BLOCK, int genus = SPHERE_GENUS);
		void SetBUniStruct(int k, Num_Block type, int genus) {
			if (0 <= k && k < N) {
				B[k].type = type;
				if (genus > ERR_GENUS) set_block3D(B[k], genus, 1.);
			}
		}
		void RotateUniStruct(double fi = 0., double theta = 0., double f0 = 0.);
		void RotateUniStruct(int k, double fi, double theta, double f0) {
			if (0 <= k && k < N && B[k].mp) {
				B[k].mp[4] = fi;
				B[k].mp[5] = theta;
				B[k].mp[6] = f0;
			}
		}
		void BlockActivate(int id_status = OK_STATE, int id_fast = OK_STATE);
		void LinkUniStruct();
		void LinkPhase2D	(int max_phase = 1) {
			if (NUM_PHASE > 0)
				for (int i, k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE) {
					for (int j = B[k].link[0]; j > 0; j--) 
					if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE]) {
						B[k].link[j] = -i+SRF_STATE;
						link_add(B[k], i);
					}
					if (B[k].link[NUM_PHASE] < -max_phase) B[k].link[NUM_PHASE] = -max_phase;
				}
		}
		void LinkPhase3D	(int max_phase = 1) {
			if (NUM_PHASE > 0) 
				for (int i, k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE) {
					for (int j = B[k].link[0]; j > 0; j--) 
					if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE]) {
						B[k].link[j] = -i+SRF_STATE;
						link_add(B[k], i);
						link_add(B[k], SRF_STATE, NULL_STATE);
					}
					if (B[k].link[NUM_PHASE] < -max_phase) B[k].link[NUM_PHASE] = -max_phase;
				}
		}
public:
//...bounding blocks;
		double MaxSkeletonRadius (Block<T> & B, int id_set_origine = NULL_STATE, int id_fast = OK_STATE);
		double MinSkeletonRadius (Block<T> & B, int id_set_origine = NULL_STATE, int id_fast = OK_STATE);
		double MaxDirectionRadius(Block<T> & B, int axis = AXIS_X, int id_fast = OK_STATE);
		void   SkeletonBounding  (Block<T> & B, double * par,    int id_fast = OK_STATE);
		void   SkeletonRadius  (int id_set_origine = NULL_STATE, int id_fast = OK_STATE);
		void	 SetGeomBounding (double * par);
		void   SetBlockBounding(double * par, int id_fast = OK_STATE, int id_reset = NULL_STATE);
		void	 SetBounding	  (double * par, int id_fast = OK_STATE) {
			SetGeomBounding (par);
			SetBlockBounding(par, id_fast);
		}
public:
//...установки блоков;
		int  add_block  (Num_Block type);
		int  set_block2D(Block<T> & B, int genus, double R, int id_dop = 0,  double x0 = 0., double y0 = 0., int id_local = NULL_STATE);
		int  set_block3D(Block<T> & B, int genus, double R, int id_dop = 0,  double x0 = 0., double y0 = 0., double z0 = 0., 
																			 double  fi = 0., double th = 0., double f0 = 0., int id_local = NULL_STATE);
		void init_blocks(int N_sm = 0, CCells * bar = NULL);
		void add_sph_surface(double R, double X0, double Y0, double Z0); 
		void add_cyl_surface(double R, double X0, double Y0, double Z0, int axis = AXIS_Z, double fi = 0., double theta = 0.); 
//...reading boundary condition from ABAQUS/NASTRAN file;
		int  bar_condit_in(char	 * id_CONDIT, unsigned long & comput, unsigned long upper_limit, int id_print_stat = NULL_STATE);
		void bar_condit_in(char	* ch_CONDIT, int id_print_stat = NULL_STATE);
		void bar_condit_in(const char * ch_CONDIT, int id_print_stat = NULL_STATE);
//...конвертация плоской структуры в призматическую или цилиндрически симметричную;
		int  ConvertToBeamStruct(int * NN, double * LL, int m_cell = NULL_STATE, int id_phase = NULL_STATE);
		void ConvertToCylStruct	(int NN, double axis_dist, double fi0 = 360., int m_close = 1);
		void CurveCorrect	(double X0, double Y0, double rad, double eps = EE_usl);
//...идентификация точек внутри системы блоков;
		int  Poly_struc_in2D (int & k, double X, double Y, int id_fast = ZERO_STATE, double eps = EE_ker);
		int  Poly_struc_in3D (int & k, double X, double Y, double Z, int id_local = ZERO_STATE, double eps = EE_ker);
		int  Sph_struc_in3D  (int & k, double X, double Y, double Z, int id_map1into = OK_STATE);
		void StructEllCorrect(int & k, double X, double Y);
		void hit_beam_struct	(CGrid * nd, int N0, int * NN, double * LL, int axis = AXIS_Z); 
//...библиотека готовых структур;
		int GetBoxStruct				 (double AA,	double BB, double CC, int NX = 1, int NY = 1, int NZ = 1);
		int GetSphBoxStruct			 (double AA, double BB, double CC, double rad = 0., double ll = 0., int id_bound = NULL_STATE);
		int GetSpheroidBoxStruct	 (double AA, double BB, double CC, double rad = 0., double delt = 0.5, double fi = 0., double theta = 0., int id_bound = NULL_STATE);
		int GetLatticeBox3DStruct	 (double * X, double * Y, double * Z, int NX, int NY, int NZ);
		int GetLatticeBox3DStruct	 (double * X, double * Y, double * Z, int NX, int NY, int NZ, Num_Block id_block);
		int GetOnoBoxStruct			 (double AA,	double BB, double aa, double bb, double rad, int * NN, double * LL, int id_phase, int id_curve, double eps);
		int GetOnoBoxStruct			 (double AA,	double BB, double aa, double bb, double rad, int * NN, double * LL, char * name, int id_phase = NULL_STATE, int id_curve = OK_STATE, double eps = EE_usl);
		int GetOnoBoxStruct			 (double AA,	double BB, double aa, double bb, double rad, int * NN, double * LL, const char * name, int id_phase = NULL_STATE, int id_curve = OK_STATE, double eps = EE_usl);
		int GetPenetrateSphere		 (double rad, double L);
		int GetCylinderStruct		 (double rad, double R, double length, int NX = 3, int NY = 1, int NZ = 1);
		int GetCylinderChannelStruct(double rad, double R, double length, int NZ, int id_curve);
		int GetCylinderChannelStruct(double rad, double R, double length, char * name, int NZ = 1, int id_curve = OK_STATE);
		int GetCylinderChannelStruct(double rad, double R, double length, const char * name, int NZ = 1, int id_curve = OK_STATE);
		int GetCylinderChannelStruct(double rad, double R, double length, double r0,  int NX = 3, int NY = 1, int N2 = 1, int NZ = 1);
		int GetCylSphStruct			 (double ff_vol, double rad, int id_bound = NULL_STATE);
		int GetCylSphStruct			 (double ff_vol, double rad, double L, int id_bound = NULL_STATE);
		int GetSph2BodyStruct		 (double ff_vol, double rad, int id_bound = NULL_STATE);
		int GetBars25Struct			 (CDraft<T> * bars25, double length, int NZ = 1);
};
/////////////////////////////////////////////////////////////////////////////////////////////////////
//...global functions for construction all existing drafts types of in all existing types of numbers;
CDraft<double>  * CreateDraftR(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8);
CDraft<complex> * CreateDraftC(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8);
CDraft<dd_real> * CreateDraftD(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8);
CDraft<qd_real> * CreateDraftQ(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8);

//////////////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR EXTERN BLOCK MATRIX SOLVER          //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////
//...structures for extern solver;
template <typename T> struct BCM_draft {
	CDraft<T> * sm; Num_Comput Num_computing;
};

//////////////////////////////////////
//...число блоков в системе уравнений;
template <typename T> 
void Number_of_Blocks(void * _pBCM, int & _NBlks)
{
	_NBlks = ((BCM_draft<T> *)_pBCM)->sm->solver.N;
}

////////////////////////
//...размерность блоков;
template <typename T> 
void Blocks_Partitioning(void * _pBCM, int * _Blks)
{
	_Blks[0] = 0;
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int * p = sm->solver.p;
	for (int i = 0; i < sm->solver.N; i++) 
		_Blks[i+1] = _Blks[i]+sm->solver.dim[p[i]];
}

/////////////////////////////////////////////////////////////////////////
//...структура разреженности блочной матрицы для несимметричного солвера;
template <typename T> 
void Blocks_Sparsity(void * _pBCM, int *& _IA, int *& _JA)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int  *  p = sm->solver.p;
	int ** JR = sm->solver.JR, i, j;

	for (j = i = 0; i < sm->solver.N; i++)
		j += JR[i][0]*2-1;

	_IA = (int *)new_struct((sm->solver.N+1)*sizeof(int));
	_JA = (int *)new_struct(j*sizeof(int));

	for (int m = 0, k = 0; k < sm->solver.N;  k++) {
		for ( i = 0; i < k; i++) {
			for (j = JR[p[i]][0]; j > 0; j--) 
			if (JR[p[i]][j+sm->solver.JR_SHIFT-1] == p[k]) break;
			if (j) _JA[m++] = i;
		}
		for (j = 0; j < JR[p[k]][0]; j++) 
		_JA[m++] = p[sm->solver.N+JR[p[k]][j+sm->solver.JR_SHIFT]];
		_IA[k+1] = m;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...заполнение строки блочной вещественной симметричной матрицы для несимметричного солвера;
template <typename T> 
void Blocks_Row(void * _pBCM, int _BlkRow, double * _defc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	if (! sm->solver.mode(PROCESSOR_ID)) {
		
///////////////////////////
//...расчет блочной строки;
		sm->solver.set_mode(NO_MESSAGE);
		sm->computing(_BlkRow, ((BCM_draft<T> *)_pBCM)->Num_computing);
		sm->solver.clean_mode(NO_MESSAGE);

////////////////////////////////
//...сборка диагональных матриц;
		int  *  p = sm->solver.p;
		int ** JR = sm->solver.JR;
		int ** IL = sm->solver.IL, i, j, k, l, m, ii;
		sm->solver.diagonal(p[_BlkRow], 1.);

///////////////////////////////
//...заполнение блочной строки;
		double alpha = max(1., sm->get_param(sm->size_of_param()-1)), beta = 1.-1./alpha;
		if (sm->solver.mode(TESTI_GRAM))   alpha = 1., beta = 0.; else
		if (sm->solver.mode(TESTI_ENERGY)) alpha = 0., beta = 1.;

		for (m = 0, i = 0; i < IL[p[_BlkRow]][0]; i++) {
			for (j = JR[ii = IL[p[_BlkRow]][i+sm->solver.JR_SHIFT]][0]; j > 0; j--)
			if ( ii != p[_BlkRow] && JR[ii][j+sm->solver.JR_SHIFT-1] == p[_BlkRow]) break;

			if (j)
			for (k = 0; k < sm->solver.dim[ii]; k++)
			for (l = 0; l < sm->solver.dim[p[_BlkRow]]; l++) 
				_defc[m++] = to_double(sm->solver.TL[ii][j-1][k][l]*alpha);
		}
		for (j = 0; j < JR[p[_BlkRow]][0]; j++) {
			if (j+sm->solver.JR_SHIFT != JR[p[_BlkRow]][sm->solver.JR_DIAG]) {
				for (l = 0; l < sm->solver.dim[JR[p[_BlkRow]][j+sm->solver.JR_SHIFT]]; l++)
				for (k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
					_defc[m++] = to_double(sm->solver.TR[p[_BlkRow]][j][k][l]*alpha);
			}
			else {
				for (l = 0; l < sm->solver.dim[p[_BlkRow]]; l++)
				for (k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
					_defc[m++] = to_double(sm->solver.TR[p[_BlkRow]][j][k][l]*alpha+
												  sm->solver.TL[p[_BlkRow]][j][k][l]*beta);
			}
		}
	}
	else { //...тестовая запись распределения по процессорам;
#ifdef ___MPI_INIT___
		extern CMPIComm comm_mpi;
		sm->solver.hh[sm->solver.p[_BlkRow]][0][sm->solver.id_norm][0] = T(comm_mpi.GetMyid()+1.);
#endif
	}
}

//////////////////////////////////////////////////
//...заполнение правой части вещественной матрицы;
template <typename T> 
void Right_Handside(void * _pBCM, int _BlkRow, double * _refc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	if (! sm->solver.mode(PROCESSOR_ID)) {

		double alpha = max(1., sm->get_param(sm->size_of_param()-1)), beta = 1.-1./alpha;
		if (sm->solver.mode(TESTI_GRAM))   alpha = 1., beta = 0.; else
		if (sm->solver.mode(TESTI_ENERGY)) alpha = 0., beta = 1.;

		int * p = sm->solver.p;
		for (int k = 0; k < sm->solver.dim[p[_BlkRow]]; k++)
			_refc[k] = to_double(sm->solver.hh[p[_BlkRow]][0][0][k]*alpha+
										sm->solver.hh[p[_BlkRow]][2][0][k]*beta);
	}
	else {//...тестовая запись распределения по процессорам;
#ifdef ___MPI_INIT___
		extern CMPIComm comm_mpi;
		_refc[0] = to_double(comm_mpi.GetMyid()+1.);

		if (sm->solver.mode(PRINT_MODE)) {
			FILE *  TST = fopen("proceccor_id.dat", "a");
			fprintf(TST, "%i = %i\n", sm->solver.p[_BlkRow], abs(comm_mpi.GetMyid())+1);
			fclose (TST);
		}
#endif
	}
}

/////////////////////////////////////
//...начальное приближение к решению;
template <typename T> 
void Initial_Guess(void * _pBCM, int _BlkRow, double * _refc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int * p = sm->solver.p;
	for (int k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
		_refc[k] = to_double(sm->solver.hh[p[_BlkRow]][0][sm->solver.id_norm][k]);
}

///////////////////////////////////////
//...перенос решения в блочную матрицу;
template <typename T> 
void Store_Solution(void * _pBCM, int _BlkRow, double * _refc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int * p = sm->solver.p;
	for (int k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
		sm->solver.hh[p[_BlkRow]][0][0][k] = T(_refc[k]);
}

//////////////////////////////////////////////////////////
//...структура разреженности блочной симметричной матрицы;
template <typename T> 
void Blocks_SparsitySym(void * _pBCM, int *& _IA, int *& _JA)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int  *  p = sm->solver.p;
	int ** JR = sm->solver.JR, i, j;

	for (j = i = 0; i < sm->solver.N; i++)
		j += JR[i][0];

	_IA = (int *)new_struct((sm->solver.N+1)*sizeof(int));
	_JA = (int *)new_struct(j*sizeof(int));

	for (int m = 0, k = 0; k < sm->solver.N; k++) {//...fill arrays, mapped rows;

		for (int j = 0; j < JR[p[k]][0]; j++) 
		_JA[m++] = p[sm->solver.N+JR[p[k]][j+sm->solver.JR_SHIFT]];
		_IA[k+1] = m;
	}
}

/////////////////////////////////////////////////////////////////
//...заполнение строки блочной вещественной симметричной матрицы;
template <typename T> 
void Blocks_RowSym(void * _pBCM, int _BlkRow, double * _defc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	int  *  p = sm->solver.p;
	int ** JR = sm->solver.JR;
	for (int m = 0, j = 0; j < JR[p[_BlkRow]][0]; j++)
	for (int l = 0; l < sm->solver.dim[JR[p[_BlkRow]][j+sm->solver.JR_SHIFT]]; l++) 
	for (int k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
		_defc[m++] = to_double(sm->solver.TR[p[_BlkRow]][j][k][l]);
}

///////////////////////////////////////////////////
//...заполнение блочной строки комплексной матрицы;
template <typename T> 
void Blocks_RowC(void * _pBCM, int _BlkRow, void * _defc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	if (! sm->solver.mode(PROCESSOR_ID)) {
		if (typeid(T) != typeid(complex)) return;

///////////////////////////
//...расчет блочной строки;
		sm->solver.set_mode(NO_MESSAGE);
		sm->computing(_BlkRow, ((BCM_draft<T> *)_pBCM)->Num_computing);
		sm->solver.clean_mode(NO_MESSAGE);

////////////////////////////////
//...сборка диагональных матриц;
		int  *  p = sm->solver.p;
		int ** JR = sm->solver.JR;
		int ** IL = sm->solver.IL, i, j, k, l, m, ii;
		sm->solver.diagonal(p[_BlkRow], 1.);

///////////////////////////////
//...заполнение блочной строки;
		double alpha = max(1., sm->get_param(sm->size_of_param()-1)), beta = 1.-1./alpha;
		if (sm->solver.mode(TESTI_GRAM))   alpha = 1., beta = 0.; else
		if (sm->solver.mode(TESTI_ENERGY)) alpha = 0., beta = 1.;

		for (m = 0, i = 0; i < IL[p[_BlkRow]][0]; i++) {
			for (j = JR[ii = IL[p[_BlkRow]][i+sm->solver.JR_SHIFT]][0]; j > 0; j--)
			if ( ii != p[_BlkRow] && JR[ii][j+sm->solver.JR_SHIFT-1] == p[_BlkRow]) break;

			if (j)
			for (k = 0; k < sm->solver.dim[ii]; k++)
			for (l = 0; l < sm->solver.dim[p[_BlkRow]]; l++) 
				((T *)_defc)[m++] = sm->solver.TL[ii][j-1][k][l]*alpha;
		}
		for (j = 0; j < JR[p[_BlkRow]][0]; j++) {
			if (j+sm->solver.JR_SHIFT != JR[p[_BlkRow]][sm->solver.JR_DIAG]) {
				for (l = 0; l < sm->solver.dim[JR[p[_BlkRow]][j+sm->solver.JR_SHIFT]]; l++)
				for (k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
					((T *)_defc)[m++] = sm->solver.TR[p[_BlkRow]][j][k][l]*alpha;
			}
			else {
				for (l = 0; l < sm->solver.dim[JR[p[_BlkRow]][j+sm->solver.JR_SHIFT]]; l++)
				for (k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
					((T *)_defc)[m++] = sm->solver.TR[p[_BlkRow]][j][k][l]*alpha+
											  sm->solver.TL[p[_BlkRow]][j][k][l]*beta;
			}
		}
	}
	else {//...тестовая запись распределения по процессорам;
#ifdef ___MPI_INIT___
		extern CMPIComm comm_mpi;
		sm->solver.hh[sm->solver.p[_BlkRow]][0][sm->solver.id_norm][0] = T(comm_mpi.GetMyid()+1.);
#endif
	}
}

/////////////////////////////////////////////////
//...заполнение правой части комплексной матрицы;
template <typename T> 
void Right_HandsideC(void * _pBCM, int _BlkRow, void * _refc)
{
	CDraft<T> * sm = ((BCM_draft<T> *)_pBCM)->sm;
	if (! sm->solver.mode(PROCESSOR_ID)) {
		if (typeid(T) != typeid(complex)) return;

		double alpha = max(1., sm->get_param(sm->size_of_param()-1)), beta = 1.-1./alpha;
		if (sm->solver.mode(TESTI_GRAM))   alpha = 1., beta = 0.; else
		if (sm->solver.mode(TESTI_ENERGY)) alpha = 0., beta = 1.;

		int * p = sm->solver.p;
		for (int k = 0; k < sm->solver.dim[p[_BlkRow]]; k++) 
			((T *)_refc)[k] = sm->solver.hh[p[_BlkRow]][0][0][k]*alpha+
									sm->solver.hh[p[_BlkRow]][2][0][k]*beta;
	}
	else {//...тестовая запись распределения по процессорам;
#ifdef ___MPI_INIT___
		extern CMPIComm comm_mpi;
		((T *)_refc)[0] = T(comm_mpi.GetMyid()+1.);

		if (sm->solver.mode(PRINT_MODE)) {
			FILE *  TST = fopen("proceccor_id.dat", "a");
			fprintf(TST, "%i = %i\n", sm->solver.p[_BlkRow], abs(comm_mpi.GetMyid())+1);
			fclose (TST);
		}
#endif
	}
}

////////////////////////////////////////////////////
//          TEMPLATE VARIANT OF CDRAFT            //
////////////////////////////////////////////////////
template <typename T> int CDraft<T>::NUM_MPLS  = 0;
template <typename T> int CDraft<T>::NUM_QUAD  = 2;
template <typename T> int CDraft<T>::NUM_TIME  = 4;
template <typename T> int CDraft<T>::NUM_ADHES = 4;
template <typename T> int CDraft<T>::NUM_VIBRO = 4;
template <typename T> int CDraft<T>::NUM_GEOMT = 4;
template <typename T> int CDraft<T>::NUM_PHASE = 0;
template <typename T> int CDraft<T>::MAX_PHASE = 3;
template <typename T> int CDraft<T>::BOX_LINK_PERIOD = 0;

////////////////////
//...blocks release;
template <typename T> 
void CDraft<T>::blocks_release()
{
   if (B) 
	while (--N >= 0) {
      delete B[N].shape;
      if (B[N].bar) {
          B[N].bar->zero_cells(); 
			 delete B[N].bar;
      }
      delete_struct(B[N].link);
      delete_struct(B[N].mp);
   }
   delete[] B; B = NULL; 
	N = buf_count = 0;
}

/////////////////////////////////
//...инициализация функций формы;
template <typename T> 
void CDraft<T>::shapes_init(Num_State id_free)
{
	for (int k = 0; B && k < N; k++)
		block_shape_init(B[k], id_free);
}

/////////////////////////////
//...defining type of blocks;
template <typename T> 
void CDraft<T>::SetBUniStruct(Num_Block type, int genus)
{
	for (int k = 0; k < N; k++) {
		B[k].type = type;
		if (genus > ERR_GENUS) set_block3D(B[k], genus, 1.);
	}
}

//////////////////////////////////////
//...rotation of the blocks structure;
template <typename T> 
void CDraft<T>::RotateUniStruct(double fi, double theta, double f0)
{
	for (int k = 0; k < N; k++) if (B[k].mp) { //...разворачиваем блоки;
		B[k].mp[4] = fi;
		B[k].mp[5] = theta;
		B[k].mp[6] = f0;
	}
}

///////////////////////////////////////////////
//...activization geonetry of blocks structure;
template <typename T> 
void CDraft<T>::BlockActivate(int id_status, int id_fast)
{
	for (int k = 0; k < N; k++)
		bar_activate(B[k], id_status, id_fast);
}

////////////////////////////////////////////////////////////////////////////
//...линкование блочной структуры на базе целочисленного описания геометрии;
template <typename T> 
void CDraft<T>::LinkUniStruct()
{
	if (stru.geom && stru.geom[0] && N_buf > 1) {
		int i, j, k, l, m, max_env, * env_ptr;
		init_blocks(stru.geom[0], NULL);

		int ** node_env = NULL; set_matrix(node_env, stru.N, N_buf+1);
		for (l = 0; l < stru.N; l++) node_env[l][0] = N_buf;

		for (i = 1, k = 0; k < stru.geom[0]; k++, i += stru.geom[++i]+1) 
		for (j = 3; j < stru.geom[i+1]; j++) if ((l  = stru.geom[i+j+2]) >= 0 && l < stru.N) {

//////////////////////////////////////////
//...буфферизация вспомогательной матрицы;
			env_ptr = node_env[l]; m = 0;
			while ((m += N_buf) && ! env_ptr[0]) env_ptr = (int *)(long long)env_ptr[1];

			if (env_ptr[0] == 1) {
				 env_ptr[1] = (int)(long long)new_struct((N_buf+1)*sizeof(int));
				 env_ptr[0] = 0;
				((int *)(long long)env_ptr[1])[0] = N_buf;
			}
			if (i == 1 || max_env < m) max_env = m;

////////////////////////////
//...заносим индексы блоков;
			env_ptr = node_env[l];
			while (! env_ptr[0]) env_ptr = (int *)(long long)env_ptr[1];
			env_ptr[ env_ptr[0]--] = k;
		}
		max_env *= 2;

//////////////////////////////////
//...формируем связи между блокми;
		int * block_env = (int *)new_struct((max_env+1)*sizeof(int)); block_env[0] = max_env;
		for (i = 1, k = 0; k < stru.geom[0]; k++, i += stru.geom[++i]+1) {

////////////////////////////
//...заносим индексы блоков;
			for (j = 3; j < stru.geom[i+1]; j++) if ((l  = stru.geom[i+j+2]) >= 0 && l < stru.N) {
				env_ptr = node_env[l];
				do {
///////////////////////////////////////////
//...буфферизация вспомогательного массива;
					if (N_buf-max(env_ptr[0], 1) > block_env[0]) {
						int * new_env = (int *)new_struct((max_env*2+1)*sizeof(int));
						memcpy(new_env+max_env+1, block_env+1, max_env*sizeof(int));
						new_env[0] = block_env[0]+max_env;
						delete_struct(block_env); block_env = new_env;
						max_env *= 2;
					}
/////////////////////////////////////////////////////
//...занесение и переход к следующей порции индексов;
					for (m = max(env_ptr[0], 1)+1; m <= N_buf; m++) if (env_ptr[m] != k) block_env[block_env[0]--] = env_ptr[m];
					if (! env_ptr[0]) env_ptr = (int *)(long long)env_ptr[1]; 
					else  env_ptr = NULL;
				}
				while (env_ptr);
			}

////////////////////////////////////////////////////////////
//...упорядочивание индексов и линкование блочной структуры;
			stru.quick_sort_inverse(block_env[0]+1, max_env, block_env);
			block_geom_link (B[k], stru.geom, stru.geom_ptr, max_env, block_env, NUM_PHASE);
			block_env[0] = max_env;
		}
	
/////////////////////////////////////
//...удаляем вспомогательные массивы;
		for (l = 0; l < stru.N; l++) 
		if (! (env_ptr = node_env[l])[0]) {
			env_ptr = (int *)(long long)env_ptr[1];
			do {
				int * temp = env_ptr; 
				if (! temp[0]) env_ptr = (int *)(long long)temp[1]; else env_ptr = NULL;
				delete_struct(temp);
			}
			while (env_ptr);
		}
		delete_struct(node_env);
		delete_struct(block_env);
	}
}

////////////////////////////////////////////////////////////////
//...links block structure (and store phase as property number);
template <typename T> 
int CDraft<T>::block_geom_link(Block<T> & B, int * geom, int * geom_ptr, int max_env, int * block_env, int N_buf)
{
	if (! geom || ! geom_ptr) return(0);
	int i, m, m1, m2, m3, m4, m5, m6, m7, m8;
   
/////////////////////////////////////
//...forming links of block elements;
	switch (geom[i = geom_ptr[(int)(&B-B.B)]]) {
		case GL_LINE_STRIP: if (geom[i+1] == 5) {
			set_link(B, N_buf = max(N_buf, max(NUM_PHASE, 2)));	
			m1 = geom[i+5];
			m2 = geom[i+6];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 2;
/////////////////////
//...all environment;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1])) {
				if (geom_link_id1s(m1, geom_ptr[block_env[m]], geom)) {
					if (B.link[0] < N_buf) B.link[++B.link[0]] = block_env[m]; 
					else  link_add(B, block_env[m]);
					B.link[1] = block_env[m];
				}
				if (geom_link_id1s(m2, geom_ptr[block_env[m]], geom)) {
					if (B.link[0] < N_buf) B.link[++B.link[0]] = block_env[m]; 
					else  link_add(B, block_env[m]);
					B.link[2] = block_env[m];
				}
			}
		}  break;
		case GL_TRIANGLES: if (geom[i+1] == 6) {
			set_link(B, N_buf = max(N_buf, 3));	
			m1 = geom[i+5];
			m2 = geom[i+6];
			m3 = geom[i+7];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 3;
////////////
//...side 1;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m1, m2, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
////////////
//...side 2;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m2, m3, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
////////////
//...side 3;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m3, m1, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
		}  break;
		case GL_QUADS:
		case GL_QUAD_STRIP: if (geom[i+1] == 7) {
			set_link(B, N_buf = max(N_buf, 4));	
			m1 = geom[i+5];
			m2 = geom[i+6];
			m3 = geom[i+7];
			m4 = geom[i+8];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 4;
////////////
//...side 1;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m1, m2, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
////////////
//...side 2;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m2, m3, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
////////////
//...side 3;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m3, m4, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
////////////
//...side 4;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id2s(m4, m1, geom_ptr[block_env[m]], geom)) {
				B.link[4] = block_env[m];
				break;
			}
		}  break;
		case GL_BOXS: if (geom[i+1] == 11) {
			set_link(B, N_buf = max(N_buf, 6));	
			m1   = geom[i+5];
			m2   = geom[i+6];
			m3   = geom[i+7];
			m4   = geom[i+8];
			m5   = geom[i+9];
			m6   = geom[i+10];
			m7   = geom[i+11];
			m8   = geom[i+12];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 6;
/////////////
//...face X';
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m1, m5, m8, m4, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
/////////////
//...face X;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m2, m3, m7, m6, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
/////////////
//...face Y';
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m1, m2, m6, m5, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
/////////////
//...face Y;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m4, m8, m7, m3, geom_ptr[block_env[m]], geom)) {
				B.link[4] = block_env[m];
				break;
			}
/////////////
//...face Z';
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m1, m4, m3, m2, geom_ptr[block_env[m]], geom)) {
				B.link[5] = block_env[m];
				break;
			}
/////////////
//...face Z;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m5, m6, m7, m8, geom_ptr[block_env[m]], geom)) {
				B.link[6] = block_env[m];
				break;
			}
		}  break;
		case GL_PENTA: if (geom[i+1] == 9) {
			set_link(B, N_buf = max(N_buf, 5));	
			m1   = geom[i+5];
			m2   = geom[i+6];
			m3   = geom[i+7];
			m4   = geom[i+8];
			m5   = geom[i+9];
			m6   = geom[i+10];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 5;
////////////////////
//...face lateral 1;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m1, m2, m5, m4, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
////////////////////
//...face lateral 2;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m2, m3, m6, m5, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
////////////////////
//...face lateral 3;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m3, m1, m4, m6, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
/////////////
//...face Z';
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m1, m3, m2, geom_ptr[block_env[m]], geom)) {
				B.link[4] = block_env[m];
				break;
			}
/////////////
//...face Z;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m4, m5, m6, geom_ptr[block_env[m]], geom)) {
				B.link[5] = block_env[m];
				break;
			}
		}  break;
		case GL_TETRA: if (geom[i+1] == 7 || geom[i+1] == 13) {
			set_link(B, N_buf = max(N_buf, 4));	
			m1   = geom[i+5];
			m2   = geom[i+6];
			m3   = geom[i+7];
			m4   = geom[i+8];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 4;
///////////////
//...lateral 1;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m1, m3, m2, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
///////////////
//...lateral 2;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m2, m3, m4, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
///////////////
//...lateral 3;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m3, m1, m4, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
///////////////
//...lateral 4;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m1, m2, m4, geom_ptr[block_env[m]], geom)) {
				B.link[4] = block_env[m];
				break;
			}
		}  break;
		case GL_PYRAMID: if (geom[i+1] == 8) {
			set_link(B, N_buf = max(N_buf, 5));	
			m1   = geom[i+5];
			m2   = geom[i+6];
			m3   = geom[i+7];
			m4   = geom[i+8];
			m5   = geom[i+9];
			if (NUM_PHASE && N_buf >= NUM_PHASE) B.link[B.link[0] = NUM_PHASE] = geom[i+4]; 
			else											 B.link[0] = 5;
///////////////
//...lateral 1;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m1, m2, m5, geom_ptr[block_env[m]], geom)) {
				B.link[1] = block_env[m];
				break;
			}
///////////////
//...lateral 2;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m2, m3, m5, geom_ptr[block_env[m]], geom)) {
				B.link[2] = block_env[m];
				break;
			}
///////////////
//...lateral 3;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m3, m4, m5, geom_ptr[block_env[m]], geom)) {
				B.link[3] = block_env[m];
				break;
			}
///////////////
//...lateral 4;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id3(m4, m1, m5, geom_ptr[block_env[m]], geom)) {
				B.link[4] = block_env[m];
				break;
			}
////////////
//...base 5;
			for (m = max_env; m > block_env[0]; m--) if ((m == max_env || block_env[m] != block_env[m+1]) &&
				geom_link_id4(m1, m4, m3, m2, geom_ptr[block_env[m]], geom)) {
				B.link[5] = block_env[m];
				break;
			}
		}  break;
	}
	return(1);
}

///////////////////////////////////////////////////////////////////////
//...определение стороны прямоугольной ячейки в формате (X:1,2, Y:3,4);
template <typename T> 
int CDraft<T>::block_iddir_2D(Block<T> & B, int i, int j, double * par, double eps)
{
   int k = (int)(&B-B.B), id_dir = 0;

	if (0 <= k && k < N && B.bar && 0 <= i && i < B.bar->graph[0]) IF_ANY_FACET(B.bar, i)
	if (0 <= j && j < B.bar->ce[i]->graph[1]) {
		  CCells * ce = B.bar->ce[B.bar->ce[i]->graph[j+2]];
		  ce->line_correct();
////////////////////////////////////////////////////////////////////
//...определяем граничную сторону ячейки и разбиваем ее на элементы;
		  double  X0 = ce->ce[get_num(ce->graph, 0)]->mp[1],
					 Y0 = ce->ce[get_num(ce->graph, 0)]->mp[2], 
					 X1 = ce->ce[get_num(ce->graph, 1)]->mp[1],
					 Y1 = ce->ce[get_num(ce->graph, 1)]->mp[2], 
					 nX =  sin(ce->mp[4])*sin(ce->mp[5]), 
					 nY = -cos(ce->mp[4])*sin(ce->mp[5]);
		  if (fabs(nX+1.) < eps && fabs(X0-par[0]) < eps && fabs(X1-par[0]) < eps) id_dir = 1; else
		  if (fabs(nX-1.) < eps && fabs(X0-par[1]) < eps && fabs(X1-par[1]) < eps) id_dir = 2; else
		  if (fabs(nY+1.) < eps && fabs(Y0-par[2]) < eps && fabs(Y1-par[2]) < eps) id_dir = 3; else
		  if (fabs(nY-1.) < eps && fabs(Y0-par[3]) < eps && fabs(Y1-par[3]) < eps) id_dir = 4; else id_dir = 0;
	}
	return(id_dir);
}

//////////////////////////////////////////////////////////////////////////////////
//...определение стороны параллелепипедной ячейки в формате (X:1,2, Y:3,4, Z:5,6);
template <typename T> 
int CDraft<T>::block_iddir_3D(Block<T> & B, int i, double * par, int id_fast, double eps)
{
	int l, k = (int)(&B-B.B), N_ini = 4, id_dir = 0;
	double pm[12];

	if (0 <= k && k < N && B.bar && 0 <= i && i < B.bar->graph[0]) IF_ANY_FACET(B.bar, i) {
//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = B.bar->ce[i]->graph[(N_arc = B.bar->ce[i]->graph[1])+1], 
			 N_min = min(N_arc, N_ini), m1, m2;
		for (l = 1; l <= N_min; l++, prev = arc) {
			arc = B.bar->ce[i]->graph[l+1];
			if (id_fast == OK_STATE)  m1 = arc;	else {
				m1  = get_num(B.bar->ce[i]->ce[arc]->graph, 0),
				m2  = get_num(B.bar->ce[i]->ce[arc]->graph, 1);
				if (! B.bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			pm[3*(l-1)]   = B.bar->ce[i]->ce[m1]->mp[1];
			pm[3*(l-1)+1] = B.bar->ce[i]->ce[m1]->mp[2];
			pm[3*(l-1)+2] = B.bar->ce[i]->ce[m1]->mp[3];
		}

/////////////////////////////////////////
//...определяем граничную сторону ячейки;
		double nX = B.bar->ce[i]->mp[4],
				 nY = B.bar->ce[i]->mp[5],
				 nZ = B.bar->ce[i]->mp[6];

		if (! id_dir) {
			if (fabs(nX+1.) < eps) {
				for (id_dir = 1, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[0]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nX-1.) < eps) {
				for (id_dir = 2, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[1]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (fabs(nY+1.) < eps) {
				for (id_dir = 3, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[2]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nY-1.) < eps) {
				for (id_dir = 4, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[3]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (fabs(nZ+1.) < eps) {
				for (id_dir = 5, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[4]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nZ-1.) < eps) {
				for (id_dir = 6, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[5]) >= eps) id_dir = 0;
			}
		}
	}
	return(id_dir);
}

//////////////////////////////////////////////////////////////////
//...определение соответственных блоков в 2D периодической задаче;
template <typename T> 
int CDraft<T>::block_plink_2D(Block<T> & B, int & i, int & j, int & id_dir, double * par, double eps)
{
   int m, l, s, r, k = (int)(&B-B.B), elem = -1; id_dir = 0;

	if (i >= 0 && j >= 0)
	if (0 <= k && k < N && B.bar && B.link)
	for (; elem < 0 && i < B.bar->graph[0]; i++, j = 0) IF_ANY_FACET(B.bar, i)
	for (; elem < 0 && j < B.bar->ce[i]->graph[1]; j++) if (B.link[j+1] < 0) {
		  CCells * ce = B.bar->ce[B.bar->ce[i]->graph[j+2]], * cell;
		  ce->line_correct();
////////////////////////////////////////////////////////////////////
//...определяем граничную сторону ячейки и разбиваем ее на элементы;
		  double  X0 = ce->ce[get_num(ce->graph, 0)]->mp[1],
					 Y0 = ce->ce[get_num(ce->graph, 0)]->mp[2], 
					 X1 = ce->ce[get_num(ce->graph, 1)]->mp[1],
					 Y1 = ce->ce[get_num(ce->graph, 1)]->mp[2], 
					 nX =  sin(ce->mp[4])*sin(ce->mp[5]), 
					 nY = -cos(ce->mp[4])*sin(ce->mp[5]);
		  if (fabs(nX+1.) < eps && fabs(X0-par[0]) < eps && fabs(X1-par[0]) < eps) id_dir = 1; else
		  if (fabs(nX-1.) < eps && fabs(X0-par[1]) < eps && fabs(X1-par[1]) < eps) id_dir = 2; else
		  if (fabs(nY+1.) < eps && fabs(Y0-par[2]) < eps && fabs(Y1-par[2]) < eps) id_dir = 3; else
		  if (fabs(nY-1.) < eps && fabs(Y0-par[3]) < eps && fabs(Y1-par[3]) < eps) id_dir = 4; else id_dir = 0;
		  if (solv%ENERGY_SOLVING != PERIODIC_SOLVING) elem = 0;

/////////////////////////////////////////////
//...определяем номер соответственной ячейки;
		  if (id_dir && (solv%ENERGY_SOLVING == PERIODIC_SOLVING))
		  for (s = 0; elem < 0 && s < N; s++) if (s != k && B.B[s].bar && B.B[s].link && B.B[s].link[NUM_PHASE] == B.link[NUM_PHASE])
		  for (m = 0; elem < 0 && m < B.B[s].bar->graph[0]; m++) IF_ANY_FACET(B.B[s].bar, m)
		  for (l = 0; elem < 0 && l < B.B[s].bar->ce[m]->graph[1]; l++) if (B.B[s].link[l+1] < 0) {
				cell = B.B[s].bar->ce[B.B[s].bar->ce[m]->graph[l+2]];
				cell->line_correct();
				double XX0 = cell->ce[get_num(cell->graph, 0)]->mp[1],
						 YY0 = cell->ce[get_num(cell->graph, 0)]->mp[2], 
						 XX1 = cell->ce[get_num(cell->graph, 1)]->mp[1],
						 YY1 = cell->ce[get_num(cell->graph, 1)]->mp[2], 
						 nXX =  sin(cell->mp[4])*sin(cell->mp[5]), 
						 nYY = -cos(cell->mp[4])*sin(cell->mp[5]);

				if (id_dir == 1 && fabs(nXX-1.) < eps && fabs(XX0-par[1]) < eps && fabs(XX1-par[1]) < eps && (fabs(Y1-YY0) < eps || fabs(Y0-YY0) < eps) && (fabs(Y0-YY1) < eps || fabs(Y1-YY1) < eps) ||
					 id_dir == 2 && fabs(nXX+1.) < eps && fabs(XX0-par[0]) < eps && fabs(XX1-par[0]) < eps && (fabs(Y1-YY0) < eps || fabs(Y0-YY0) < eps) && (fabs(Y0-YY1) < eps || fabs(Y1-YY1) < eps) ||
					 id_dir == 3 && fabs(nYY-1.) < eps && fabs(YY0-par[3]) < eps && fabs(YY1-par[3]) < eps && (fabs(X1-XX0) < eps || fabs(X0-XX0) < eps) && (fabs(X0-XX1) < eps || fabs(X1-XX1) < eps) ||
					 id_dir == 4 && fabs(nYY+1.) < eps && fabs(YY0-par[2]) < eps && fabs(YY1-par[2]) < eps && (fabs(X1-XX0) < eps || fabs(X0-XX0) < eps) && (fabs(X0-XX1) < eps || fabs(X1-XX1) < eps)) elem = s; 
				if (elem >= 0) 
            for (r = NUM_PHASE; r < B.link[0]; r++) 
            if  (B.link[j+1] == -B.link[r+1]+SRF_STATE) elem = -1;
		  }
		  if (elem >= 0 && (solv%ENERGY_SOLVING == PERIODIC_SOLVING)) {
			  j++;
			  return(elem);
		  }
	}
	return(-1);
}

//////////////////////////////////////////////////////////////////
//...определение соответственных блоков в 3D периодической задаче;
template <typename T> 
int CDraft<T>::block_plink_3D(Block<T> & B, int & i, int & id_dir, double * par, int id_fast, double eps)
{
	int m, l, s, k = (int)(&B-B.B), j, mm, nn, elem = -1, N_ini = 4;
	double pm[12], pp_cond[12];

	if (i >= 0)
	if (0 <= k && k < N && B.bar && B.link)
	for (; elem < 0 && i < B.bar->graph[0]; i++) if (B.link[i+1] < 0)
	IF_ANY_FACET(B.bar, i) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = B.bar->ce[i]->graph[(N_arc = B.bar->ce[i]->graph[1])+1], 
			 N_min = min(N_arc, N_ini), m1, m2, NN_min; id_dir = 0;
		for (l = 1; l <= N_min; l++, prev = arc) {
			arc = B.bar->ce[i]->graph[l+1];
			if (id_fast == OK_STATE) m1 = arc; else {
				m1  = get_num(B.bar->ce[i]->ce[arc]->graph, 0),
				m2  = get_num(B.bar->ce[i]->ce[arc]->graph, 1);
				if (! B.bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			pm[3*(l-1)]   = B.bar->ce[i]->ce[m1]->mp[1];
			pm[3*(l-1)+1] = B.bar->ce[i]->ce[m1]->mp[2];
			pm[3*(l-1)+2] = B.bar->ce[i]->ce[m1]->mp[3];
		}

/////////////////////////////////////////
//...определяем граничную сторону ячейки;
		double nX = B.bar->ce[i]->mp[4],
				 nY = B.bar->ce[i]->mp[5],
				 nZ = B.bar->ce[i]->mp[6];

		if (! id_dir) {
			if (fabs(nX+1.) < eps) {
				for (id_dir = 1, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[0]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nX-1.) < eps) {
				for (id_dir = 2, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[1]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (fabs(nY+1.) < eps) {
				for (id_dir = 3, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[2]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nY-1.) < eps) {
				for (id_dir = 4, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[3]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (fabs(nZ+1.) < eps) {
				for (id_dir = 5, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[4]) >= eps) id_dir = 0;
			}
			else
			if (fabs(nZ-1.) < eps) {
				for (id_dir = 6, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[5]) >= eps) id_dir = 0;
			}
		}
		if (solv%ENERGY_SOLVING != PERIODIC_SOLVING) elem = 0;

/////////////////////////////////////////////
//...определяем номер соответственной ячейки;
		if (id_dir && (solv%ENERGY_SOLVING == PERIODIC_SOLVING))
		for (s = 0; elem < 0 && s < N; s++) if (s != k && B.B[s].bar && B.B[s].link && B.B[s].link[NUM_PHASE] == B.link[NUM_PHASE])
		for (m = 0; elem < 0 && m < B.B[s].bar->graph[0]; m++) IF_ANY_FACET(B.B[s].bar, m) 
		if (B.B[s].link[m+1] < 0) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
			prev = B.B[s].bar->ce[m]->graph[(N_arc = B.B[s].bar->ce[m]->graph[1])+1]; 
			NN_min = min(N_arc, N_ini);
			for (j = 1; j <= NN_min; j++, prev = arc) {
				arc = B.B[s].bar->ce[m]->graph[j+1];
				if (id_fast == OK_STATE)  m1 = arc;	else {
					m1  = get_num(B.B[s].bar->ce[m]->ce[arc]->graph, 0),
					m2  = get_num(B.B[s].bar->ce[m]->ce[arc]->graph, 1);
					if (! B.B[s].bar->ce[m]->ce[prev]->topo_id(m1)) swap(m1, m2);
				}
				pp_cond[3*(j-1)]   = B.B[s].bar->ce[m]->ce[m1]->mp[1];
				pp_cond[3*(j-1)+1] = B.B[s].bar->ce[m]->ce[m1]->mp[2];
				pp_cond[3*(j-1)+2] = B.B[s].bar->ce[m]->ce[m1]->mp[3];
			}

///////////////////////////////////////
//...определяем соответственную ячейку;
			double nXX = B.B[s].bar->ce[m]->mp[4],
					 nYY = B.B[s].bar->ce[m]->mp[5],
					 nZZ = B.B[s].bar->ce[m]->mp[6];

			mm = nn = 0;
			if (NN_min == N_min) {
				if (id_dir == 1 && fabs(nXX-1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)]-par[1]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)+1]-pm[3*(l-1)+1]) < eps && 
							 fabs(pp_cond[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 2 && fabs(nXX+1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)]-par[0]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)+1]-pm[3*(l-1)+1]) < eps && 
							 fabs(pp_cond[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 3 && fabs(nYY-1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)+1]-par[3]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp_cond[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 4 && fabs(nYY+1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)+1]-par[2]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp_cond[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 5 && fabs(nZZ-1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)+2]-par[5]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp_cond[3*(j-1)+1]-pm[3*(l-1)+1]) < eps) nn = 1;
				}
				else
				if (id_dir == 6 && fabs(nZZ+1.) < eps) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp_cond[3*(j-1)+2]-par[4]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp_cond[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp_cond[3*(j-1)+1]-pm[3*(l-1)+1]) < eps) nn = 1;
				}
			}
			if (mm && nn) elem = s;
			if (elem >= 0) 
			for (j = NUM_PHASE; j < B.link[0]; j++) 
			if  (B.link[i+1] == -B.link[j+1]+SRF_STATE) elem = -1;
		}
		if (elem >= 0 && (solv%ENERGY_SOLVING == PERIODIC_SOLVING)) { //...убрать i++ !!!
			i++; return(elem);
		}
	}
	return(-1);
}

///////////////////////////////////////////////////////////////////////
//...определение стороны прямоугольной ячейки в формате (X:1,2, Y:3,4);
template <typename T> 
int CDraft<T>::geom_iddir_2D(Block<T> & B, int i, double * par, double eps)
{
   int m, m1, m2, j, k, id_dir = 0;
	if (i >= 0 && stru.geom && stru.geom_ptr && stru.X && stru.Y)
	for (m = stru.geom[(j = stru.geom_ptr[k = (int)(&B-B.B)])+1]-3; i < m; i++) {
		m1 = stru.geom[j+5+i%m];
		m2 = stru.geom[j+5+(i+1)%m];
		double X0 = stru.X[m1], Y0 = stru.Y[m1], 
				 X1 = stru.X[m2], Y1 = stru.Y[m2];
		if (Y1 < Y0 && fabs(X0-par[0]) < eps && fabs(X1-par[0]) < eps) id_dir = 1; else
		if (Y0 < Y1 && fabs(X0-par[1]) < eps && fabs(X1-par[1]) < eps) id_dir = 2; else
		if (X0 < X1 && fabs(Y0-par[2]) < eps && fabs(Y1-par[2]) < eps) id_dir = 3; else
		if (X1 < X0 && fabs(Y0-par[3]) < eps && fabs(Y1-par[3]) < eps) id_dir = 4; else id_dir = 0;
	}
	return(id_dir);
}

//////////////////////////////////////////////////////////////////////////////////
//...определение стороны параллелепипедной ячейки в формате (X:1,2, Y:3,4, Z:5,6);
template <typename T> 
int CDraft<T>::geom_iddir_3D(Block<T> & B, int i, double * par, double eps)
{
	int l, mf, k, N_ini = 4, id_dir = 0, mm_facet[6][5];
	double pm[12];

	if ((mf = stru.stru_install(k = (int)(&B-B.B), mm_facet)) != 0 && 
		  0 <= i && i < mf && stru.X && stru.Y) {
//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_min = mm_facet[i][0];
		for (l = 1;	l <= N_min; l++) {
			pm[3*(l-1)]   = stru.X[mm_facet[i][l]];
			pm[3*(l-1)+1] = stru.Y[mm_facet[i][l]];
			pm[3*(l-1)+2] = stru.Z[mm_facet[i][l]];
		}
      double nX = (pm[4]-pm[1])*(pm[N_min*3-1]-pm[2])-(pm[5]-pm[2])*(pm[N_min*3-2]-pm[1]),
             nY = (pm[5]-pm[2])*(pm[N_min*3-3]-pm[0])-(pm[3]-pm[0])*(pm[N_min*3-1]-pm[2]),
             nZ = (pm[3]-pm[0])*(pm[N_min*3-2]-pm[1])-(pm[4]-pm[1])*(pm[N_min*3-3]-pm[0]);

/////////////////////////////////////////
//...определяем граничную сторону ячейки;
		if (! id_dir) {
			if (nX < 0.) {
				for (id_dir = 1, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[0]) >= eps) id_dir = 0;
			}
			else
			if (nX > 0.) {
				for (id_dir = 2, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[1]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (nY < 0.) {
				for (id_dir = 3, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[2]) >= eps) id_dir = 0;
			}
			else
			if (nY > 0.) {
				for (id_dir = 4, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[3]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (nZ < 0.) {
				for (id_dir = 5, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[4]) >= eps) id_dir = 0;
			}
			else
			if (nZ > 0.) {
				for (id_dir = 6, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[5]) >= eps) id_dir = 0;
			}
		}
	}
	return(id_dir);
}

//////////////////////////////////////////////////////////////////////////////////////
//...определение соответственных блоков в 2D периодической задаче на основе геометрии;
template <typename T> 
int CDraft<T>::geom_plink_2D(Block<T> & B, int & i, int & id_dir, double * par, double eps)
{
   int m, m1, m2, j, k, l, s, r, elem = -1; id_dir = 0;
	if (i >= 0 && B.link && stru.geom && stru.geom_ptr && stru.X && stru.Y) {
		for (m = stru.geom[(j = stru.geom_ptr[k = (int)(&B-B.B)])+1]-3; elem < 0 && i < m; i++) 
		if (B.link[i+1] < 0) {
			m1 = stru.geom[j+5+i%m];
			m2 = stru.geom[j+5+(i+1)%m];
			double X0 = stru.X[m1], Y0 = stru.Y[m1], 
					 X1 = stru.X[m2], Y1 = stru.Y[m2];
			if (Y1 < Y0 && fabs(X0-par[0]) < eps && fabs(X1-par[0]) < eps) id_dir = 1; else
			if (Y0 < Y1 && fabs(X0-par[1]) < eps && fabs(X1-par[1]) < eps) id_dir = 2; else
			if (X0 < X1 && fabs(Y0-par[2]) < eps && fabs(Y1-par[2]) < eps) id_dir = 3; else
			if (X1 < X0 && fabs(Y0-par[3]) < eps && fabs(Y1-par[3]) < eps) id_dir = 4; else id_dir = 0;
			if (solv%ENERGY_SOLVING != PERIODIC_SOLVING) elem = 0;

/////////////////////////////////////////////
//...определяем номер соответственной ячейки;
			if (id_dir && (solv%ENERGY_SOLVING == PERIODIC_SOLVING))
			for (s = 0; elem < 0 && s < stru.geom[0]; s++) if (s != k && B.B[s].link && B.B[s].link[NUM_PHASE] == B.link[NUM_PHASE])
			for (r = stru.geom[(j = stru.geom_ptr[s])+1]-3, l = 0; elem < 0 && l < r; l++) 
			if (B.B[s].link[l+1] < 0) {
				m1 = stru.geom[j+5+l%r];
				m2 = stru.geom[j+5+(l+1)%r];
				double XX0 = stru.X[m1], YY0 = stru.Y[m1], 
						 XX1 = stru.X[m2], YY1 = stru.Y[m2];
				if (id_dir == 1 && YY0 < YY1 && fabs(XX0-par[1]) < eps && fabs(XX1-par[1]) < eps && (fabs(Y1-YY0) < eps || fabs(Y0-YY0) < eps) && (fabs(Y0-YY1) < eps || fabs(Y1-YY1) < eps) ||
					 id_dir == 2 && YY1 < YY0 && fabs(XX0-par[0]) < eps && fabs(XX1-par[0]) < eps && (fabs(Y1-YY0) < eps || fabs(Y0-YY0) < eps) && (fabs(Y0-YY1) < eps || fabs(Y1-YY1) < eps) ||
					 id_dir == 3 && XX1 < XX0 && fabs(YY0-par[3]) < eps && fabs(YY1-par[3]) < eps && (fabs(X1-XX0) < eps || fabs(X0-XX0) < eps) && (fabs(X0-XX1) < eps || fabs(X1-XX1) < eps) ||
					 id_dir == 4 && XX0 < XX1 && fabs(YY0-par[2]) < eps && fabs(YY1-par[2]) < eps && (fabs(X1-XX0) < eps || fabs(X0-XX0) < eps) && (fabs(X0-XX1) < eps || fabs(X1-XX1) < eps)) elem = s; 
				if (elem >= 0) 
				for (j = NUM_PHASE; j < B.link[0]; j++) 
				if  (B.link[i+1] == -B.link[j+1]+SRF_STATE) elem = -1;
			}
			if (elem >= 0 && (solv%ENERGY_SOLVING == PERIODIC_SOLVING)) {
				i++;
				return(elem);
			}
		}
	}
	return(-1);
}

//////////////////////////////////////////////////////////////////////////////////////
//...определение соответственных блоков в 3D периодической задаче на основе геометрии;
template <typename T> 
int CDraft<T>::geom_plink_3D(Block<T> & B, int & i, int & id_dir, double * par, double eps)
{
	int m, j, l, s, mf, nf, mm, nn, k, elem = -1, N_ini = 4, mm_facet[6][5], nn_facet[6][5];
	double pm[12], pp[12];

	if (0 <= i && B.link && stru.X && stru.Y && stru.Z && (mf = stru.stru_install(k = (int)(&B-B.B), mm_facet)) != 0)
	for (; elem < 0 && i < mf; i++) 
	if (B.link[i+1] < 0) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_min = mm_facet[i][0]; id_dir = 0;
		for (l = 1;	l <= N_min; l++) {
			pm[3*(l-1)]   = stru.X[mm_facet[i][l]];
			pm[3*(l-1)+1] = stru.Y[mm_facet[i][l]];
			pm[3*(l-1)+2] = stru.Z[mm_facet[i][l]];
		}
      double nX = (pm[4]-pm[1])*(pm[N_min*3-1]-pm[2])-(pm[5]-pm[2])*(pm[N_min*3-2]-pm[1]),
             nY = (pm[5]-pm[2])*(pm[N_min*3-3]-pm[0])-(pm[3]-pm[0])*(pm[N_min*3-1]-pm[2]),
             nZ = (pm[3]-pm[0])*(pm[N_min*3-2]-pm[1])-(pm[4]-pm[1])*(pm[N_min*3-3]-pm[0]);

/////////////////////////////////////////
//...определяем граничную сторону ячейки;
		if (! id_dir) {
			if (nX < 0.) {
				for (id_dir = 1, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[0]) >= eps) id_dir = 0;
			}
			else
			if (nX > 0.) {
				for (id_dir = 2, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)]-par[1]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (nY < 0.) {
				for (id_dir = 3, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[2]) >= eps) id_dir = 0;
			}
			else
			if (nY > 0.) {
				for (id_dir = 4, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+1]-par[3]) >= eps) id_dir = 0;
			}
		}
		if (! id_dir) {
			if (nZ < 0.) {
				for (id_dir = 5, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[4]) >= eps) id_dir = 0;
			}
			else
			if (nZ > 0.) {
				for (id_dir = 6, l = 1; id_dir && l <= N_min; l++) 
					if (fabs(pm[3*(l-1)+2]-par[5]) >= eps) id_dir = 0;
			}
		}
		if (solv%ENERGY_SOLVING != PERIODIC_SOLVING) elem = 0;

/////////////////////////////////////////////
//...определяем номер соответственной ячейки;
		if (id_dir && (solv%ENERGY_SOLVING == PERIODIC_SOLVING))
		for (s = 0; elem < 0 && s < N; s++) if ((N == 1 || s != k) && B.B[s].link && B.B[s].link[NUM_PHASE] == B.link[NUM_PHASE] &&
			(nf = stru.stru_install(s, nn_facet)) != 0)
		for (m = 0; elem < 0 && m < nf; m++) if (B.B[s].link[m+1] < 0) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
			int NN_min = nn_facet[m][0];
			for (j = 1; j <= NN_min; j++) {
				pp[3*(j-1)]   = stru.X[nn_facet[m][j]];
				pp[3*(j-1)+1] = stru.Y[nn_facet[m][j]];
				pp[3*(j-1)+2] = stru.Z[nn_facet[m][j]];
			}
			double nXX = (pp[4]-pp[1])*(pp[NN_min*3-1]-pp[2])-(pp[5]-pp[2])*(pp[NN_min*3-2]-pp[1]),
					 nYY = (pp[5]-pp[2])*(pp[NN_min*3-3]-pp[0])-(pp[3]-pp[0])*(pp[NN_min*3-1]-pp[2]),
					 nZZ = (pp[3]-pp[0])*(pp[NN_min*3-2]-pp[1])-(pp[4]-pp[1])*(pp[NN_min*3-3]-pp[0]);

///////////////////////////////////////
//...определяем соответственную ячейку;
			mm = nn = 0;
			if (NN_min == N_min) {
				if (id_dir == 1 && nXX > 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)]-par[1]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)+1]-pm[3*(l-1)+1]) < eps && 
							 fabs(pp[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 2 && nXX < 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)]-par[0]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)+1]-pm[3*(l-1)+1]) < eps && 
							 fabs(pp[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 3 && nYY > 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)+1]-par[3]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 4 && nYY < 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)+1]-par[2]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp[3*(j-1)+2]-pm[3*(l-1)+2]) < eps) nn = 1;
				}
				else
				if (id_dir == 5 && nZZ > 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)+2]-par[5]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp[3*(j-1)+1]-pm[3*(l-1)+1]) < eps) nn = 1;
				}
				else
				if (id_dir == 6 && nZZ < 0.) {
					for (mm = nn = 1, j = 1; mm && j <= N_min; j++) 
						if (fabs(pp[3*(j-1)+2]-par[4]) >= eps) mm = 0;

					for (    j = 1; mm && nn && j <= N_min; j++) 
					for (nn = 0, l = 1; ! nn && l <= N_min; l++) 
						if (fabs(pp[3*(j-1)]  -pm[3*(l-1)])   < eps && 
							 fabs(pp[3*(j-1)+1]-pm[3*(l-1)+1]) < eps) nn = 1;
				}
			}
			if (mm && nn) elem = s;
			if (elem >= 0) 
			for (j = NUM_PHASE; j < B.link[0]; j++) 
			if  (B.link[i+1] == -B.link[j+1]+SRF_STATE) elem = -1;
		}
		if (elem >= 0 && (solv%ENERGY_SOLVING == PERIODIC_SOLVING)) { //...убрать i++ !!!
			i++; return(elem);
		}
		mf = stru.stru_install(k, mm_facet);
	}
	return(-1);
}

////////////////////////////////////////
//...активизация геометрии одного блока;
template <typename T> 
void CDraft<T>::bar_activate(Block<T> & B, int id_status, int id_fast)
{
	if (! B.bar) {
			B.bar = new CCells;
		B.bar->get_nd_bar_directly(& stru, (int)(&B-B.B));
		
//////////////////////////////////////////
//...переустанавливаем нормирующий радиус;
		if (B.mp && size_of_map(B.mp) > 7)
			 B.mp[7] = MaxSkeletonRadius(B, OK_STATE, id_fast);
		block_shape_init(B, SPECIAL_STATE);

		if (id_status) {
///////////////////////////////////////////
//...устанавливаем линки в геометрию фасет;
			if (B.link) {
				if (NUM_PHASE) {
					int  i, l;
					for (i = 1;	i < NUM_PHASE; i++) if (B.link[i] <= SRF_STATE)
					for (l = NUM_PHASE; l < B.link[0]; l++)
					if  (B.link[l+1] >= 0 && B.link[i] == -B.link[l+1]+SRF_STATE)
					B.bar->SetFacetParam(i-1, -B.link[l+1]-1.);
					
					for (i = 1; i < NUM_PHASE; i++) if (B.link[i] >= 0)
					B.bar->SetFacetParam(i-1, -B.link[i]-1.);
				}
				else
					for (int i = 0; i < B.link[0]; i++) if (B.link[i+1] >= 0)
					B.bar->SetFacetParam (i, -B.link[i+1]-1.);
			}

/////////////////////////////////////////////////////////////
//...заносим граничные условия в фасетные элементы геометрии;
			if (stru.cond && stru.cond_ptr && id_prop && pp_cond) {
				double square;
				int l, m, k, i/*, num*/;
				for (i = 0;  i < stru.cond[0]; i++)
					if (stru.cond[l = stru.cond_ptr[i]] == GL_PLOAD4 && stru.cond[l+3] == (int)(&B-B.B)) {
						for (m = id_prop[0]; m > 0; m--) if (id_prop[m*2-1] == stru.cond[l+2]) break;
						if (m && 
							(k = B.bar->set_nd_bar_condit(stru.cond[l+3], & stru, stru.cond[l+4], stru.cond[l+5], m, pp_cond, id_prop, square)) < 0)
//						if (B.B[-k-1].bar)
						B.B[-k-1].bar->set_nd_bar_condit(stru.cond[l+3], m, pp_cond, id_prop, square);
				}

#ifdef  ___CONSTRUCTION___
//???????????????????????????????????????????????????????????????????//
//...заносим граничные условия в геометрию из ABAQUS'а (переделать!!!);
				for (i = 0;  i < stru.cond[0]; i++) 
				if (stru.cond[l = stru.cond_ptr[i]] == GL_ELEMENTS)
				for (int j = stru.cond[l+1]; j > 1; j--) if (stru.cond[l+1+j] == (int)(&B-B.B)) {
					for (m = id_prop[0]; m > 0; m--) if (id_prop[m*2-1] == -stru.cond[l+2]) break;
					if (m) {
						for (k = B.bar->arcs_number(); k > 0; k--) if (B.link[k] <= ERR_STATE) {
							for (num = NUM_PHASE; num < B.link[0]; num++) 
							if (B.link[k] == -B.link[num+1]+SRF_STATE) break;
							if (num == B.link[0])
								((CCells *)B.bar)->set_nd_bar_condit(stru.cond[l+3], & stru, -1, -1, m, pp_cond, id_prop, square, k-1);
						}
					}
				}
				for (i = 0;  i < stru.cond[0]; i++) 
				if (stru.cond[l = stru.cond_ptr[i]] == GL_ELEMENTS_GENERATE)
				for (int j = stru.cond[l+3]; j <= stru.cond[l+4]; j += stru.cond[l+5]) if (j == (int)(&B-B.B)) {
					for (m = id_prop[0]; m > 0; m--) if (id_prop[m*2-1] == -stru.cond[l+2]) break;
					if (m) {
						for (k = B.bar->arcs_number(); k > 0; k--) if (B.link[k] <= ERR_STATE) {
							for (num = NUM_PHASE; num < B.link[0]; num++) 
							if (B.link[k] == -B.link[num+1]+SRF_STATE) break;
							if (num == B.link[0])
								((CCells *)B.bar)->set_nd_bar_condit(stru.cond[l+3], & stru, -1, -1, m, pp_cond, id_prop, square, k-1);
						}
					}
				}
//???????????????????????????????????????????????????????????????????//
#endif
			}

/////////////////////////////////////////
//...устанавливаем криволинейную границу;
			int num, m, j;
			if (B.bar && B.link) LOOP_FACET(B.bar, num, m)
			if (B.link[0] >= num) {
				for (j = MARKER1_BND; j < NUMS_BND; j++)
				if (B.bar->ce[num]->mp[m+5] == j-SPECIAL_BND+1.) B.link[num+1] = SRF_STATE-j+SPECIAL_BND+1;
			}
		}
	}
}


/////////////////////////////////////////////////////
//...defining of the maximal radius of block element;
template <typename T> 
double CDraft<T>::MaxSkeletonRadius(Block<T> & B, int id_set_origine, int id_fast)
{
  double f, rad = 0., p[3] = {0., 0., 0.};
  int    m, l, k, i; 

//////////////////////
//...centroid setting;
  if (id_set_origine && B.bar && B.mp && B.bar->graph) {
      for (i =  k = 0; k < B.bar->graph[0]; k++) 
		if  (2 == (l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
      for (l = 0;  l < B.bar->ce[k]->graph[1];   l++) {
           p[0] += B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1];
           p[1] += B.bar->ce[m                           ]->mp[2];
           p[2] += B.bar->ce[m                           ]->mp[3];
           i    += 1;
      }
      if (i) {
          B.mp[1] = p[0]*(f = 1./i);
          B.mp[2] = p[1]* f;
          B.mp[3] = p[2]* f;
      }
  }

////////////////////////
//...radius calculation;
	if (B.bar && B.mp && B.bar->graph)
	for (int m, l, k = 0; k < B.bar->graph[0]; k++) 
	if  (2 == (l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
	for (l = 0;  l < B.bar->ce[k]->graph[1];   l++) {
       p[0] = B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1]-B.mp[1];
       p[1] = B.bar->ce[m                           ]->mp[2]-B.mp[2];
       p[2] = B.bar->ce[m                           ]->mp[3]-B.mp[3];
       rad  = max(f = sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2])), rad);
	}
	return(rad);
}

/////////////////////////////////////////////////////
//...defining of the minimal radius of block element;
template <typename T> 
double CDraft<T>::MinSkeletonRadius(Block<T> & B, int id_set_origine, int id_fast)
{
  double f, rad = MAX_HIT, p[3] = {0., 0., 0.};
  int    m, l, k, i; 

//////////////////////
//...centroid setting;
  if (id_set_origine && B.bar && B.mp && B.bar->graph) {
      for (i =  k = 0; k < B.bar->graph[0]; k++) 
		if  (2 == (l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
      for (l = 0;  l < B.bar->ce[k]->graph[1];   l++) {
           p[0] += B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1];
           p[1] += B.bar->ce[m                           ]->mp[2];
           p[2] += B.bar->ce[m                           ]->mp[3];
           i    += 1;
      }
      if (i) {
          B.mp[1] = p[0]*(f = 1./i);
          B.mp[2] = p[1]* f;
          B.mp[3] = p[2]* f;
      }
  }

////////////////////////
//...radius calculation;
	if (B.bar && B.mp && B.bar->graph)
	for (int m, l, k = 0; k < B.bar->graph[0]; k++) 
	if  (2 == (l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
	for (l = 0;  l < B.bar->ce[k]->graph[1];   l++) {
       p[0] = B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1]-B.mp[1];
       p[1] = B.bar->ce[m                           ]->mp[2]-B.mp[2];
       p[2] = B.bar->ce[m                           ]->mp[3]-B.mp[3];
       rad  = min(f = sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2])), rad);
  }
/*
  if (B.bar && B.mp && B.bar->graph)
  for (int m, k = 0; k < B.bar->graph[0]; k++) 
  if  (B.bar->ce[k]->cells_dim() == 2 && 
       B.bar->ce[k]->mp[m = size_of_map (B.bar->ce[k]->mp)] == FACET_CELL) {

       p[2] = (B.bar->ce[k]->mp[1]-B.mp[1])*B.bar->ce[k]->mp[4]+
              (B.bar->ce[k]->mp[2]-B.mp[2])*B.bar->ce[k]->mp[5]+
              (B.bar->ce[k]->mp[3]-B.mp[3])*B.bar->ce[k]->mp[6];
       rad  = min(fabs(p[2]), rad);
  }
*/
	return(rad);
}

//////////////////////////////////////////////////////////////////
//...defining of the maximal radius of block element in direction;
template <typename T> 
double CDraft<T>::MaxDirectionRadius(Block<T> & B, int axis, int id_fast)
{
  double rad = 0., p[3] = {0., 0., 0.};

////////////////////////
//...radius calculation;
	if (B.bar && B.mp && B.bar->graph)
	for (int m, l, k = 0; k < B.bar->graph[0]; k++) 
	if  (2 == (l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
	for (l = 0;  l < B.bar->ce[k]->graph[1];   l++) {
		p[0] = B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1];
		p[1] = B.bar->ce[m                           ]->mp[2];
		p[2] = B.bar->ce[m                           ]->mp[3];
		make_local(p, B.mp, OK_STATE);
		switch (axis) {
			case AXIS_X: rad = max(fabs(p[0]), rad); break;
			case AXIS_Y: rad = max(fabs(p[1]), rad); break;
			case AXIS_Z: rad = max(fabs(p[2]), rad); break;
		}
	}
	return(rad);
}

///////////////////////////////////////////////////
//...defining of the bounding box of block element;
template <typename T> 
void CDraft<T>::SkeletonBounding(Block<T> & B, double * par, int id_fast)
{
//////////////////////////
//...bounding calculation;
	if (B.bar && B.mp && B.bar->graph && par)
	for (int m, l, k = 0; k < B.bar->graph[0]; k++)
	if  (2 == ( l = B.bar->ce[k]->cells_dim()) && id_fast == OK_STATE || 1 == l && id_fast != OK_STATE)
	for (l = 0; l < B.bar->ce[k]->graph[1]; l++) {
		par[0] = min(par[0], B.bar->ce[m = B.bar->ce[k]->graph[l+2]]->mp[1]);
		par[1] = max(par[1], B.bar->ce[m                           ]->mp[1]);
		par[2] = min(par[2], B.bar->ce[m                           ]->mp[2]);
		par[3] = max(par[3], B.bar->ce[m                           ]->mp[2]);
		par[4] = min(par[4], B.bar->ce[m                           ]->mp[3]);
		par[5] = max(par[5], B.bar->ce[m                           ]->mp[3]);
	}
}

////////////////////////////////////////////////////////
//...defining of the bounding box for the all structure;
template <typename T> 
void CDraft<T>::SetBlockBounding(double * par, int id_fast, int id_reset)
{
	if (N && B[0].bar) {
		if (par && id_reset) {
			par[0] = par[2] = par[4] = MAX_HIT;
			par[1] = par[3] = par[5] = MIN_HIT;
		}
		for (int k = 0; k < N; k++) 
			SkeletonBounding(B[k], par, id_fast);
	}
}

/////////////////////////////////////////////////////
//...defining of the skeleton radius of the geometry;
template <typename T> 
void CDraft<T>::SkeletonRadius(int id_set_origine, int id_fast)
{
	for (int k = 0; k < N; k++)
		if (B[k].bar && B[k].mp && size_of_map(B[k].mp) > 7)
			B[k].mp[7] = MaxSkeletonRadius(B[k], id_set_origine, id_fast);
}

///////////////////////////////////////////////////
//...defining of the bounding box on grid geometry;
template <typename T> 
void CDraft<T>::SetGeomBounding(double * par)
{
	if (par) {
		par[0] = par[2] = par[4] = MAX_HIT;
		par[1] = par[3] = par[5] = MIN_HIT;
	}
	if (stru.geom && stru.X && stru.Y && stru.Z && stru.N > 0) {
		int i, j, k, l;
		for (i = 1, k = 0; k < stru.geom[0]; k++, i += stru.geom[++i]+1) 
		for (j = 2; j < stru.geom[i+1]; j++) if ((l  = stru.geom[i+j+2]) >= 0 && l < stru.N) {
			if (i == 1 && j == 2) {
				 if (stru.X) par[0] = par[1] = stru.X[l];
				 if (stru.Y) par[2] = par[3] = stru.Y[l];
				 if (stru.Z) par[4] = par[5] = stru.Z[l];
			}
			if (stru.X) {
				par[0] = min(par[0], stru.X[l]);
				par[1] = max(par[1], stru.X[l]);
			}
			if (stru.Y) {
				par[2] = min(par[2], stru.Y[l]);
				par[3] = max(par[3], stru.Y[l]);
			}
			if (stru.Z) {
				par[4] = min(par[4], stru.Z[l]);
				par[5] = max(par[5], stru.Z[l]);
			}
		}
	}
}

//////////////////////////////////
//...blocks increase in structure;
template <typename T> 
int CDraft<T>::add_block(Num_Block type)
{
  if (! buf_count) {
      Block<T> * new_B = (Block<T> *)new_struct((N+N_buf)*sizeof(Block<T>));
      int  i = N;
      while (--i >= 0) {
             new_B[i] = B[i];
             new_B[i].B = new_B;
      }
      delete[] B; B = new_B; buf_count = N_buf;
  }
  B[N].B = B;
  B[N++].type = type; buf_count--;

  return(1);
}

///////////////////////////////////////////////////////////
//...setting type and location of the multipoles (2D case);
template <typename T> 
int CDraft<T>::set_block2D(Block<T> & B, int genus, double R, int id_dop, double x0, double y0, int id_local)
{
  int k = size_of_map(1, genus); delete_struct (B.mp);
  if (k < 9 && fabs(R) != 0. || 
      k > 8 && fabs(R) <  EE_ker || id_dop < 0) return(0);
  if ((B.mp = (CMap *)new_struct((k+1+id_dop)*sizeof(CMap))) != NULL) {
       B.mp[0] = ID_MAP(1, genus);
       B.mp[1] = x0;
       B.mp[2] = y0;
       if (k > 8) {
           B.mp[7] = R;
           B.mp[8] = M_PI;
       }
       B.mp[k] = (CMap)NULL_CELL;

       if (id_local) B.bar->cells_local(B.mp);
       return(1);
  }
  return(0);
}

///////////////////////////////////////////////////////////
//...setting type and location of the multipoles (3D case);
template <typename T> 
int CDraft<T>::set_block3D(Block<T> & B, int genus, double R, int id_dop, double x0, double y0, double z0,
                                  double fi, double th, double f0, int id_local)
{
  int k = size_of_map(2, genus); delete_struct (B.mp);
  if (k < 8 || fabs(R) < EE_ker || id_dop < 0) return(0);
  if ((B.mp = (CMap *)new_struct((k+1+id_dop)*sizeof(CMap))) != NULL) {
       B.mp[0] = ID_MAP(2, genus);
       B.mp[1] = x0;
       B.mp[2] = y0;
       B.mp[3] = z0;
       B.mp[4] = fi;
       B.mp[5] = th;
       B.mp[6] = f0;
       B.mp[7] = R;
       B.mp[k] = (CMap)NULL_CELL;
       if (id_local) B.bar->cells_local(B.mp);
       return(1);
  }
  return(0);
}

/////////////////////////////////////////////
//...initialization of the sample of N parts;
template <typename T> 
void CDraft<T>::init_blocks(int N_sm, CCells * bar)
{
   blocks_release();
   if (CDraft<T>::bar) {
       CDraft<T>::bar->zero_cells(); 
       delete CDraft<T>::bar;
   }
   CDraft<T>::bar = bar;
   N              = (N_sm = abs(N_sm));
   B              = (Block<T> *)new_struct(N_sm*sizeof(Block<T>));
   while (--N_sm >= 0) B[N_sm].B = B;
}

////////////////////////////////////////////////////////////////
//...добавление сферической поверхности для коррекции квадратур;
template <typename T> 
void CDraft<T>::add_sph_surface(double R, double X0, double Y0, double Z0)
{
	if (! bar) bar = new(CCells);
	CCells	* ce  = new(CCells);
	int l = size_of_map(2, SPHERE_GENUS);
	ce->cells_new(1, 2, l+1);
	ce->mp[0] = (CMap)ID_MAP(2, SPHERE_GENUS);
	ce->mp[1] = (CMap)X0;
	ce->mp[2] = (CMap)Y0;
	ce->mp[3] = (CMap)Z0;
	ce->mp[7] = (CMap)R;
	ce->mp[l] = (CMap)NULL_CELL;
	bar->bar_add(ce);
	bar->bar_ord();
}

///////////////////////////////////////////////////////////////////
//...добавление цилиндрической поверхности для коррекции квадратур;
template <typename T> 
void CDraft<T>::add_cyl_surface(double R, double X0, double Y0, double Z0, int axis, double fi, double theta)
{
	if (! bar) bar = new(CCells);
	CCells	* ce  = new(CCells);
	int l = size_of_map(2, CYL_GENUS);
	ce->cells_new(1, 2, l+1);
	ce->mp[0] = (CMap)ID_MAP(2, CYL_GENUS);
	ce->mp[4] = (CMap)(axis == AXIS_Y ? M_PI_2 : 0.);
	ce->mp[5] = (CMap)(axis == AXIS_Y || axis == AXIS_X ? M_PI_2 : 0.);
	ce->cells_iso(NULL, fi, theta);
	ce->mp[1] = (CMap)X0;
	ce->mp[2] = (CMap)Y0;
	ce->mp[3] = (CMap)Z0;
	ce->mp[7] = (CMap)R;
	ce->mp[l] = (CMap)NULL_CELL;
	bar->bar_add(ce);
	bar->bar_ord();
}

/////////////////////////////////////////////////////////
//...подготовка данных для визуализации в формате Surfer;
template <typename T> 
void CDraft<T>::GetSurferFormat(FILE * SURF, FILE * SURF1, FILE * SURF2, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	size_t res;
	if (nd && nd->N > 0 && nd->N1 > 0 && SURF) {
		short int i0 = (short int)nd->N, 
                j0 = (short int)nd->N1;
		T *	 out_F = (T *)new_struct(surfer_dim(type())*sizeof(T)); 
		double X, Y, Z, 
				 min1F = 0., max1F = 1., 
				 min2F = 0., max2F = 1., 
				 min3F = 0., max3F = 1.;
		int hit;
		if (SURF) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF);
			res = fwrite(& i0,   sizeof(short int), 1, SURF); res = fwrite(& j0,         sizeof(short int), 1, SURF);
			res = fwrite(nd->X,  sizeof(double),  1, SURF); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF);
			res = fwrite(& min1F, sizeof(double), 1, SURF); res = fwrite(& max1F,        sizeof(double),  1, SURF);
		}
		if (SURF1) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF1);
			res = fwrite(& i0,   sizeof(short int), 1, SURF1); res = fwrite(& j0,         sizeof(short int), 1, SURF1);
			res = fwrite(nd->X,  sizeof(double),  1, SURF1); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF1);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF1); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF1);
			res = fwrite(& min2F, sizeof(double), 1, SURF1); res = fwrite(& max2F,        sizeof(double),  1, SURF1);
		}
		if (SURF2) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF2);
			res = fwrite(& i0,   sizeof(short int), 1, SURF2); res = fwrite(& j0,         sizeof(short int), 1, SURF2);
			res = fwrite(nd->X,  sizeof(double),  1, SURF2); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF2);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF2); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF2);
			res = fwrite(& min3F, sizeof(double), 1, SURF2); res = fwrite(& max3F,        sizeof(double),  1, SURF2);
		}
		min1F = min2F = min3F = MAX_HIT; 
		max1F = max2F = max3F = MIN_HIT;
      for (int j = 0; j < nd->N1; j++)
      for (int i = 0; i < nd->N;  i++) {
			if (id_axis == AXIS_X) {
				Y = nd->X[i];
				Z = nd->Y[j];
				X = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Y) {
            Z = nd->X[i];
            X = nd->Y[j];
            Y = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Z) {
            X = nd->X[i];
            Y = nd->Y[j];
            Z = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_SPH) {
            Z = nd->N2 == 1 ? nd->Z[0] : 1.;
            X = cos(nd->X[i])*sin(nd->Y[j])*Z*1.;
            Y = sin(nd->X[i])*sin(nd->Y[j])*Z*2.;
            Z = cos(nd->Y[j])*Z*3.;
			}
			else
			if (id_axis == AXIS_TIME) {
            X = nd->X[i];
            Y = 0.;
            Z = 0.;
				set_param(6, nd->Y[j]);
			}

			float ff = NOT_HIT; 
			if (nd->hit) hit = nd->hit[i+j*nd->N];
			if (hit != -1)  {
				if (SURF) {
					GetFuncAllValues(X, Y, Z, out_F, hit, _FMF, id_variant, iparam);
               res = fwrite(& (ff = (float)to_double(out_F[0])), sizeof(float), 1, SURF);
               if (min1F > to_double(out_F[0])) min1F = to_double(out_F[0]);
               if (max1F < to_double(out_F[0])) max1F = to_double(out_F[0]);
				}
				if (SURF1) {
               res = fwrite(& (ff = (float)to_double(out_F[1])), sizeof(float), 1, SURF1);
               if (min2F > to_double(out_F[1])) min2F = to_double(out_F[1]);
               if (max2F < to_double(out_F[1])) max2F = to_double(out_F[1]);
				}
				if (SURF2) {
               res = fwrite(& (ff = (float)to_double(out_F[2])), sizeof(float), 1, SURF2);
               if (min3F > to_double(out_F[2])) min3F = to_double(out_F[2]);
               if (max3F < to_double(out_F[2])) max3F = to_double(out_F[2]);
				}
			}
			else {
				if (SURF ) res = fwrite(& ff, sizeof(float), 1, SURF);
				if (SURF1) res = fwrite(& ff, sizeof(float), 1, SURF1);
				if (SURF2) res = fwrite(& ff, sizeof(float), 1, SURF2);
			}
      }
      delete_struct(out_F);

//////////////////////////////////////////////////////////////
//...перезапись максимального и минимального значения функции;
		if (SURF) {
			res = fseek(SURF, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min1F, sizeof(double), 1, SURF);
			res = fwrite(& max1F, sizeof(double), 1, SURF);
			res = fseek(SURF, 0L, SEEK_END);
		}
		if (SURF1) {
			res = fseek(SURF1, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min2F, sizeof(double), 1, SURF1);
			res = fwrite(& max2F, sizeof(double), 1, SURF1);
			res = fseek(SURF1, 0L, SEEK_END);
		}
		if (SURF2) {
			res = fseek(SURF2, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min3F, sizeof(double), 1, SURF2);
			res = fwrite(& max3F, sizeof(double), 1, SURF2);
			res = fseek(SURF2, 0L, SEEK_END);
		}
	}
}


////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetSurferFormat(char * SURF_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * SURF = NULL, * SURF1 = NULL, * SURF2 = NULL;
	::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	if (surfer_dim(type()) > 0) SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	if (surfer_dim(type()) > 1) SURF1 = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_2.grd");
	if (surfer_dim(type()) > 2) SURF2 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, SURF2, nd, _FMF, id_variant, id_axis, iparam);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
	if (SURF2) fclose(SURF2);
}
template <typename T> 
void CDraft<T>::GetSurferFormat(const char * SURF_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * SURF = NULL, * SURF1 = NULL, * SURF2 = NULL;
	::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	if (surfer_dim(type()) > 0) SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	if (surfer_dim(type()) > 1) SURF1 = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_2.grd");
	if (surfer_dim(type()) > 2) SURF2 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, SURF2, nd, _FMF, id_variant, id_axis, iparam);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
	if (SURF2) fclose(SURF2);
}

//////////////////////////////////////////////////////////
//...подготовка данных для визуализации в формате таблицы;
template <typename T> 
void CDraft<T>::GetDataFormat(FILE * DATA, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	int dim = surfer_dim(type());
	if (nd && nd->N > 0 && nd->N1 > 0 && DATA) {
		T *	 out_F = (T *)new_struct(surfer_dim(type())*sizeof(T));
		double X, Y, Z;
		int hit;
      for (int j = 0; j < nd->N1; j++)
      for (int i = 0; i < nd->N;  i++) {
			if (id_axis == AXIS_X) {
				Y = nd->X[i];
				Z = nd->Y[j];
				X = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Y) {
            Z = nd->X[i];
            X = nd->Y[j];
            Y = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Z) {
            X = nd->X[i];
            Y = nd->Y[j];
            Z = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_SPH) {
            Z = nd->N2 == 1 ? nd->Z[0] : 1.;
            X = cos(nd->X[i])*sin(nd->Y[j])*Z*1.;
            Y = sin(nd->X[i])*sin(nd->Y[j])*Z*2.;
            Z = cos(nd->Y[j])*Z*3.;
			}
			else
			if (id_axis == AXIS_TIME) {
            X = nd->X[i];
            Y = 0.;
            Z = 0.;
				set_param(6, nd->Y[j]);
			}

			float ff = NOT_HIT, f1, f2, f3; 
			if (nd->hit) hit = nd->hit[i+j*nd->N];
			if (hit != -1)  {
				if (DATA) {
					GetFuncAllValues(X, Y, Z, out_F, hit, _FMF, id_variant, iparam);
               f1 = (float)to_double(out_F[0]);
               f2 = (float)to_double(out_F[1]);
               f3 = (float)to_double(out_F[2]);
					if (dim > 2) fprintf(DATA, "%g   %g   %g   %g   %g   %g\n", X, Y, Z, f1, f2, f3); else
					if (dim > 1) fprintf(DATA, "%g   %g   %g   %g   %g\n", X, Y, Z, f1, f2); else
					if (dim > 0) fprintf(DATA, "%g   %g   %g   %g\n", X, Y, Z, f1);
				}
			}
			else {
				if (DATA) {
					if (dim > 2) fprintf(DATA, "%g   %g   %g   %g   %g   %g\n", X, Y, Z, ff, ff, ff); else
					if (dim > 1) fprintf(DATA, "%g   %g   %g   %g   %g\n", X, Y, Z, ff, ff); else
					if (dim > 0) fprintf(DATA, "%g   %g   %g   %g\n", X, Y, Z, ff);
				}
			}
      }
      delete_struct(out_F);
	}
	if (DATA) fclose(DATA);
}

////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetDataFormat(char * DATA_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * DATA = NULL;
	::strcpy(buff, DATA_FILE); strcat(buff, ".dat");
	int dim = surfer_dim(type());
	if (dim > 0)   DATA = fopen(buff, "w");
	GetDataFormat (DATA, nd, _FMF, id_variant, id_axis, iparam);
	if (DATA) fclose(DATA);
}
template <typename T> 
void CDraft<T>::GetDataFormat(const char * DATA_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * DATA = NULL;
	::strcpy(buff, DATA_FILE); strcat(buff, ".dat");
	int dim = surfer_dim(type());
	if (dim > 0)   DATA = fopen(buff, "w");
	GetDataFormat (DATA, nd, _FMF, id_variant, id_axis, iparam);
	if (DATA) fclose(DATA);
}

/////////////////////////////////////////////////////////////////
//...подготовка данных в формате CSV для пространственной задачи;
template <typename T> 
void CDraft<T>::GetCsvFormat(FILE * CSV, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	float ff[9] = {NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT}; 
	if (CSV && nd) {
      T * out_F = (T *)new_struct(9*sizeof(T));
      if (out_F && nd->X && nd->Y && nd->Z && nd->geom) {
			switch (id_variant) {
			case 0: fprintf(CSV, ",\"Ux\",\"Uy\",\"Uz\"\n"); break;
			case 1: fprintf(CSV, ",\"txx\",\"tyy\",\"tzz\",\"txy\",\"txz\",\"tyz\"\n"); break;
			}
			Num_Value _FU0 = SPL_VALUE, _FU1 = FLUX_X_VALUE, _FU2 = FLUX_Y_VALUE, _FU3 = FLUX_Z_VALUE;

//////////////////////////////////////////////////////////////////
//...записываем результаты на диск в узлах визуализационной сетки;
          if (! id_centroid && nd->hit)
          for (int k = 0;  k < nd->N; k++) {
				 if (id_variant == 0) {
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F, nd->hit[k], _FU0, id_variant);
					fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg\n", k+1/*(nd->hit ? nd->hit[k] : (k+1))*/, 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
				 }
				 if (id_variant == 1) {
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F,   nd->hit[k], _FU1, id_variant);
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F+3, nd->hit[k], _FU2, id_variant);
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F+6, nd->hit[k], _FU3, id_variant);
               
					fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg,%0.15lg,%0.15lg,%0.15lg\n", k+1/*(nd->hit ? nd->hit[k] : (k+1))*/, 
						ff[0] = (float)to_double(out_F[0]), ff[4] = (float)to_double(out_F[4]), ff[8] = (float)to_double(out_F[8]),
						ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]), ff[5] = (float)to_double(out_F[5]));
				 }
          }
          else

/////////////////////////////////////////////////////////
//...записываем результаты на диск в цетроидах элементов;
          if (! bnd)
          for (int j = 1, k = 0;  k < N; k++) {
               double X0 = 0., Y0 = 0., Z0 = 0., f;
               int id_set_centroid = (B[k].type & ERR_CODE) == ZOOM_BLOCK, l, m, i, num;
//////////////////////
//...centroid setting;
               if (id_set_centroid && B[k].bar->graph) {
                   for (i = num = 0; num < B[k].bar->graph[0]; num++) 
                   if  (1 == B[k].bar->ce[num]->cells_dim()) 
                   for (l = 0;  l < B[k].bar->ce[num]->graph[1];   l++) {
                        X0 += B[k].bar->ce[m = B[k].bar->ce[num]->graph[l+2]]->mp[1];
								Y0 += B[k].bar->ce[m                                ]->mp[2];
                        Z0 += B[k].bar->ce[m                                ]->mp[3];
                        i  += 1;
                   }
                   if (i) {
                       X0 *= (f = 1./i);
                       Y0 *=  f;
                       Z0 *=  f;
                   }
               }
               else {
                   X0 = B[k].mp[1];
                   Y0 = B[k].mp[2];
                   Z0 = B[k].mp[3];
               }
               GetFuncAllValues(X0, Y0, Z0, out_F, k, _FU0, id_variant);
               fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg\n", -nd->geom[j+2], 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
               j += nd->geom[++j]+1;
          }
          else

///////////////////////////////////////////////////////////////////////////////
//...записываем результаты на диск в цетроидах элементов более подробной сетки;
          if (bnd->geom && 0 < id_centroid && id_centroid <= bnd->geom[0])
			 for (int k = bnd->geom[4*(id_centroid-1)+3]; k < bnd->geom[4*(id_centroid-1)+4]; k++) {
               double X = bnd->X[k], Y = bnd->Y[k], Z = bnd->Z[k];
               int id_block = (int)bnd->hit[k], 
                   ex_block = (int)bnd->nY[k];

               GetFuncAllValues(X, Y, Z, out_F, id_block, _FU0, id_variant);

               fprintf(CSV, "%d,%0.15lg, %0.15lg,%0.15lg\n", ex_block, 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
          }
		}

///////////////////////
//...addditional point;
      double X0, Y0, Z0;
		int hit = -1;
		X0 = .5;
		Y0 =  1.;
		Z0 = .5;
      Poly_struc_in3D (hit, X0, Y0, Z0);
		GetFuncAllValues(X0, Y0, Z0, out_F, hit, SPL_VALUE, id_variant);
		//fprintf(CSV, "additional point:%0.15lg,%0.15lg,%0.15lg\n", 
		//	ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
      delete_struct(out_F);
	}
}

////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetCsvFormat(char * CSV_FILE, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	FILE * CSV = fopen(CSV_FILE, "w");
	GetCsvFormat(CSV, nd, id_variant, id_centroid, bnd);
	if (CSV) fclose(CSV);
}
template <typename T> 
void CDraft<T>::GetCsvFormat(const char * CSV_FILE, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	FILE * CSV = fopen(CSV_FILE, "w");
	GetCsvFormat(CSV, nd, id_variant, id_centroid, bnd);
	if (CSV) fclose(CSV);
}

////////////////////////////////////////////////////////////
//...reading boundary condition from NASTRAN or ABAQUS file;
template <typename T> 
int  CDraft<T>::bar_condit_in(char * id_CONDIT, unsigned long & count, unsigned long upper_limit, int id_print_stat)
{
  if (id_CONDIT &&  count < upper_limit) {
		char   buf[BUF_SIZE+1];
		int    i, j, k, buf_count = 0;
		double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12;

		delete_struct (id_prop); id_prop = (int *)new_struct((1+N_buf*2)*sizeof(int));
		delete_struct (pp_cond); pp_cond = (double *)new_struct(N_buf*6*sizeof(double));
		if (id_prop && pp_cond) buf_count = N_buf;

///////////////////////////////////////
//...зачитываем типы граничных условий;
		unsigned long count1;
		while (user_Read (buf, id_CONDIT, count, upper_limit) && ::strcmp(buf, "FEMAP") && ::strcmp(buf, "GEOSTATIC"));
		while (! ::strcmp(buf, "FEMAP") || ! ::strcmp(buf, "GEOSTATIC")) {
			user_Read(buf, id_CONDIT, count, upper_limit, 8); i = DEFAULT_BND;
			if (! ::strcmp(buf, "LOAD SET")) {
                 user_Count(     id_CONDIT, count, count1, ':');
             if (user_Read (buf, id_CONDIT, count, upper_limit, count1-count)) j = atoi(buf); else j = -1;
                 user_Read (buf, id_CONDIT, count, upper_limit);
				 if (user_Read (buf, id_CONDIT, count, upper_limit) && ! ::strcmp(buf, "RIGID")) {
                 p1 = p2 = p3 = p4 = p5 = p6 = 0.;
                 i  = RIGID_BND;
             }
				 else if (! ::strcmp(buf, "PN")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  p2 = MIN_HIT;
                  p1 = p5 = p6 = 0.;
                  i  = PN_BND;
             }
				 else if (! ::strcmp(buf, "VEL")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p5 = user_strtod(buf); else p5 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p6 = user_strtod(buf); else p6 = 0.; 
                  i  = VEL_BND;
             }
				 else if (! ::strcmp(buf, "VN")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  p1 = p2 = p5 = p6 = 0.;
                  i  = VN_BND;
             }
				 else if (! ::strcmp(buf, "VNTABLE")) {
                  p1 = p2 = p3 = p5 = p6 = 0.; p4 = MIN_HIT;
                  i  = VNTABLE_BND;
             }
				 else if (! ::strcmp(buf, "ABSORB")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  p3 = p4 = p5 = p6 = 0.;
                  i  = ABSORB_BND;
             }
				 else if (! ::strcmp(buf, "ABSORBTABLE")) {
                  p1 = 1.; p2 = p3 = 0.; p4 = MAX_HIT; p5 = p6 = 0.;
                  i  = ABSORBTABLE_BND;
             }
				 else if (! ::strcmp(buf, "SOURCE")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
						p5 = 0.;
                  p6 = -1.;
                  i  = SOURCE_BND;

                  int    N_sou = src ? (int)src[0] : 0;
                  double * new_src = (double *)new_struct((6*(N_sou+1)+1)*sizeof(double));
                  if (new_src) {
                      if (src) memcpy(new_src, src, (6*N_sou+1)*sizeof(double));
                      delete_struct(src); src = new_src;
                      src[6*N_sou+1] = p1; 
                      src[6*N_sou+2] = p2; 
                      src[6*N_sou+3] = p3; 
                      src[6*N_sou+4] = p4; 
                      src[6*N_sou+5] = p5; 
                      src[6*N_sou+6] = p6; 
                      src[0] += 1.;
                  }
             }
				 else if (! ::strcmp(buf, "SOURCE2")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  p5 = 0.;
                  p6 = -1.;
						i  = SOURCE2_BND;

                  int    N_sou = src ? (int)src[0] : 0;
                  double * new_src = (double *)new_struct((6*(N_sou+2)+1)*sizeof(double));
                  if (new_src) {
                      if (src) memcpy(new_src, src, (6*N_sou+1)*sizeof(double));
                      delete_struct(src); src = new_src;
                      src[6*N_sou+1] = p1; 
                      src[6*N_sou+2] = p2; 
                      src[6*N_sou+3] = p3; 
                      src[6*N_sou+4] = p4; 
                      src[6*N_sou+5] = p5; 
                      src[6*N_sou+6] = p6; 
                      src[0] += 1.;
                  }
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p7 = user_strtod(buf); else p7 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p8 = user_strtod(buf); else p8 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p9 = user_strtod(buf); else p9 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p10 = user_strtod(buf); else p10 = 0.; 
                  p11 = 0.;
                  p12 = -1.;

                  if (new_src) {
                      src[6*N_sou+7 ] = p7; 
                      src[6*N_sou+8 ] = p8; 
							 src[6*N_sou+9 ] = p9;
                      src[6*N_sou+10] = p10; 
                      src[6*N_sou+11] = p11; 
                      src[6*N_sou+12] = p12; 
                      src[0] += 1.;
                      src[6*N_sou+5 ] = -src[0]; //...initialization of double source for GetFunctionValues;
                  }
						p4 = p7; 
						p5 = p8; 
						p6 = p9;
             }
				 else if (! ::strcmp(buf, "BSOURCE")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p5 = user_strtod(buf); else p5= 0.; 
                  p6 = -1.-p5;
                  i  = BSOURCE_BND;
                  int    N_sou = src ? (int)src[0] : 0;
                  double * new_src = (double *)new_struct((6*(N_sou+1)+1)*sizeof(double));
                  if (new_src) {
                      if (src) memcpy(new_src, src, (6*N_sou+1)*sizeof(double));
                      delete_struct(src); src = new_src;
							 src[6*N_sou+1] = p1;
                      src[6*N_sou+2] = p2; 
                      src[6*N_sou+3] = p3; 
                      src[6*N_sou+4] = p4; 
                      src[6*N_sou+5] = p5; 
                      src[6*N_sou+6] = p6; 
                      src[0] += 1.;
                  }
             }
				 else if (! ::strcmp(buf, "BSOURCE2")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p5 = user_strtod(buf); else p5= 0.; 
                  p6 = -1.-p5;
                  i  = BSOURCE2_BND;
                  int    N_sou = src ? (int)src[0] : 0;
                  double * new_src = (double *)new_struct((6*(N_sou+2)+1)*sizeof(double));
                  if (new_src) {
                      if (src) memcpy(new_src, src, (6*N_sou+1)*sizeof(double));
                      delete_struct(src); src = new_src;
                      src[6*N_sou+1] = p1; 
                      src[6*N_sou+2] = p2; 
							 src[6*N_sou+3] = p3;
                      src[6*N_sou+4] = p4; 
                      src[6*N_sou+5] = p5; 
                      src[6*N_sou+6] = p6; 
                      src[0] += 1.;
                  }
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p7 = user_strtod(buf); else p7 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p8 = user_strtod(buf); else p8 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p9 = user_strtod(buf); else p9 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p10 = user_strtod(buf); else p10 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p11 = user_strtod(buf); else p11= 0.; 
                  p12 = -1.-p11;
                  if (new_src) {
                      src[6*N_sou+7 ] = p7; 
                      src[6*N_sou+8 ] = p8; 
                      src[6*N_sou+9 ] = p9; 
                      src[6*N_sou+10] = p10; 
                      src[6*N_sou+11] = p11; 
                      src[6*N_sou+12] = p12; 
                      src[0] += 1.;
                      src[6*N_sou+5 ] = -src[0]; //...initialization of double source for GetFunctionValues;
                  }
						p4 = p7; 
						p5 = p8; 
						p6 = p9;
             }
				 else if (! ::strcmp(buf, "EXTERN")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p2 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p3 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p4 = user_strtod(buf); else p4 = 0.; 
                  p5 = p6 = 0.;
                  i  = EXTERN_BND;
             }
				 else if (! ::strcmp(buf, "DISPLN")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  p2 = p3 = p5 = p6 = 0.;
						p4 = MIN_HIT;
                  i  = DISPLN_BND;
             }
				 else if (! ::strcmp(buf, "FORCEN")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  p2 = p3 = p5 = p6 = 0.;
						p4 = 1.;
                  i  = FORCEN_BND;
             }
				 else if (! ::strcmp(buf, "ADHESION")) {
                  p1 = p2 = p3 = p5 = p6 = 0.;
						p4 = 0;
                  i  = ADHESION_BND;
             }
				 else if (! ::strcmp(buf, "PRESSURE")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  p2 = p3 = p5 = p6 = 0.;
						p4 = MAX_HIT;
                  i  = PRESSURE_BND;
             }
				 else if (! ::strcmp(buf, "VELOCITY")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
						if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.;
                  p5 = p6 = 0.;
						p4 = 1.;
                  i  = VELOCITY_BND;
             }
				 else if (! ::strcmp(buf, "VELOCITY_N")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  p2 = p3 = p5 = p6 = 0.;
						p4 = 2.;
                  i  = VELOCITY_N_BND;
             }
				 else if (! ::strcmp(buf, "RIGID_DISPL3D")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
						if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.;
                  p5 = p6 = 0.;
						p4 = MIN_HIT;
                  i  = RIGID_DISPL3D_BND;
             }
				 else if (! ::strcmp(buf, "DISPL3D")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
						if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.;
                  p5 = p6 = 0.;
						p4 = MAX_HIT;
                  i  = DISPL3D_BND;
             }
				 else if (! ::strcmp(buf, "FORCE3D")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.; 
                  p5 = p6 = 0.;
						p4 = MARKER1_BND-SPECIAL_BND;
                  i  = FORCE3D_BND;
             }
				 else if (! ::strcmp(buf, "MOMENT3D")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.; 
                  p5 = p6 = 0.;
						p4 = MOMENTS_BND-SPECIAL_BND;
                  i  = MOMENT3D_BND;
             }
				 else if (! ::strcmp(buf, "NORMS")) {
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p2 = user_strtod(buf); else p1 = 0.; 
                  if (user_Read (buf, id_CONDIT, count, upper_limit)) p3 = user_strtod(buf); else p1 = 0.; 
                  p5 = p6 = 0.;
						p4 = NORMS_BND-SPECIAL_BND;
						i  = FNORMS_BND;
             }
				 else if (! ::strcmp(buf, "SKEWS")) {
						if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
                  p2 = p3 = p5 = p6 = 0.;
						p4 = SKEWS_BND-SPECIAL_BND;
                  i  = FSKEWS_BND;
             }
				 else if (! ::strcmp(buf, "UPTAKE")) {
						if (user_Read (buf, id_CONDIT, count, upper_limit)) p1 = user_strtod(buf); else p1 = 0.; 
					   p1 = 0.;
                  p2 = p3 = p5 = p6 = 0.;
						p4 = SKEWS_BND-SPECIAL_BND;
                  i  = UPTAKES_BND;
             }
				 else if (! ::strcmp(buf, "SPECIAL")) {
						p1 = p2 = p3 = p5 = p6 = 0.;
						p4 = NUMS_BND;
						i  = SPECIAL_BND;
             }
             else {
					 for (k = MARKER1_BND; k < NUMS_BND; k++) {
						 char str[200];
						 sprintf(str, "MARKER%i", k-SPECIAL_BND);
						 if (! ::strcmp(buf, str)) {
							p1 = p2 = p3 = p5 = p6 = 0.;
							p4 = k-SPECIAL_BND+1.;
							i  = k;
							break;
      				 }
					 }
					 if (k == NUMS_BND) p1 = p2 = p3 = p4 = p5 = p6 = 0.;
				 }

////////////////////////////////////////////////////////////
//...буфферизация и запись граничных условий в общий массив;
             if (! buf_count) {
                 int    * new_id_prop = (int *)new_struct(((N_buf+id_prop[0])*2+1)*sizeof(int));
					  double * new_pp_cond = (double *)new_struct((N_buf+id_prop[0])*6*sizeof(double));
                 if (new_id_prop) {
                     memcpy   (new_id_prop,  id_prop, (id_prop[0]*2+1)*sizeof(int));
                     delete_struct(id_prop); id_prop = new_id_prop;
                 }
                 if (new_pp_cond && id_prop) {
                     memcpy	(new_pp_cond,  pp_cond, id_prop[0]*4*sizeof(double));
                     delete_struct(pp_cond); pp_cond = new_pp_cond;
                 }
                 else delete_struct(pp_cond);
                 if (id_prop && pp_cond) buf_count = N_buf;
             }
             if (buf_count--) {
                 pp_cond[id_prop[0]*6  ] = p1;
                 pp_cond[id_prop[0]*6+1] = p2;
                 pp_cond[id_prop[0]*6+2] = p3;
                 pp_cond[id_prop[0]*6+3] = p4;
                 pp_cond[id_prop[0]*6+4] = p5;
                 pp_cond[id_prop[0]*6+5] = p6;
                 id_prop[id_prop[0]*2+1] = j;
                 id_prop[id_prop[0]*2+2] = i;
                 id_prop[0]++;
             }
         }
			while (user_Read (buf, id_CONDIT, count, upper_limit) && ::strcmp(buf, "FEMAP") && ::strcmp(buf, "GEOSTATIC"));
      }

///////////////////////////////////////////
//...распечатываем файл граничных значений;
      if (id_prop && pp_cond && id_print_stat) {
          FILE *  TST = fopen("BoundCondit.sta", "w");
          fprintf(TST, "Boundary condition, N_loading_sets = %d: \n", id_prop[0]);
          for (j = 0; j < id_prop[0]; j++) {
              fprintf(TST, " %d  ", id_prop[j*2+1]);

              if (id_prop[j*2+2] == RIGID_BND)
                  fprintf(TST, "%s\n", " RIGID_BND");
              else
              if (id_prop[j*2+2] == PN_BND)
                  fprintf(TST, "%s    %g  %g\n", "    PN_BND", pp_cond[j*6+2], pp_cond[j*6+3]);
              else
              if (id_prop[j*2+2] == VEL_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g  %g  %g\n", "VEL_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3], pp_cond[j*6+4], pp_cond[j*6+5]);
              else
              if (id_prop[j*2+2] == VN_BND)
                  fprintf(TST, "%s    %g  %g\n", "    VN_BND", pp_cond[j*6+2], pp_cond[j*6+3]);
              else
              if (id_prop[j*2+2] == VNTABLE_BND)
                  fprintf(TST, "%s\n", "VNTABLE_BND");
				  else
              if (id_prop[j*2+2] == ABSORB_BND)
                  fprintf(TST, "%s    %g  %g\n", "ABSORB_BND",  pp_cond[j*6], pp_cond[j*6+1]);
              else
              if (id_prop[j*2+2] == ABSORBTABLE_BND)
                  fprintf(TST, "%s\n", "ABSORBTABLE_BND");
              else
              if (id_prop[j*2+2] == SOURCE_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g\n", "SOURCE_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3]);
              else
              if (id_prop[j*2+2] == SOURCE2_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g  %g  %g\n", "SOURCE2_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3], pp_cond[j*6+4], pp_cond[j*6+5]);
              else
              if (id_prop[j*2+2] == BSOURCE_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g  %g\n", "BSOURCE_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3], pp_cond[j*6+4]);
              else
              if (id_prop[j*2+2] == BSOURCE2_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g  %g  %g\n", "BSOURCE2_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3], pp_cond[j*6+4], pp_cond[j*6+5]);
              else
              if (id_prop[j*2+2] == EXTERN_BND)
                  fprintf(TST, "%s    %g  %g  %g  %g\n", "EXTERN_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2], pp_cond[j*6+3]);
              else
              if (id_prop[j*2+2] == DISPLN_BND)
                  fprintf(TST, "%s    %g\n", "DISPLN_BND", pp_cond[j*6]);
              else
				  if (id_prop[j*2+2] == FORCEN_BND)
                  fprintf(TST, "%s    %g\n", "FORCEN_BND", pp_cond[j*6]);
              else
              if (id_prop[j*2+2] == ADHESION_BND)
                  fprintf(TST, "%s\n", "ADHESION_BND");
              else
				  if (id_prop[j*2+2] == PRESSURE_BND)
                  fprintf(TST, "%s    %g\n", "PRESSURE_BND", pp_cond[j*6]);
              else
              if (id_prop[j*2+2] == VELOCITY_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "VELOCITY_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
				  if (id_prop[j*2+2] == VELOCITY_N_BND)
                  fprintf(TST, "%s    %g\n", "VELOCITY_N_BND", pp_cond[j*6]);
				  else
              if (id_prop[j*2+2] == RIGID_DISPL3D_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "RIGID_DISPL3D_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
              if (id_prop[j*2+2] == DISPL3D_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "DISPL3D_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
              if (id_prop[j*2+2] == FORCE3D_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "FORCE3D_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
              if (id_prop[j*2+2] == MOMENT3D_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "MOMENT3D_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
              if (id_prop[j*2+2] == FNORMS_BND)
                  fprintf(TST, "%s    %g  %g  %g\n", "FNORMS_BND", pp_cond[j*6], pp_cond[j*6+1], pp_cond[j*6+2]);
              else
              if (id_prop[j*2+2] == FSKEWS_BND)
                  fprintf(TST, "%s    %g\n", "FSKEWS_BND", pp_cond[j*6]);
              else
              if (id_prop[j*2+2] == UPTAKES_BND)
                  fprintf(TST, "%s    %g\n", "UPTAKES_BND", pp_cond[j*6]);
              else
              if (id_prop[j*2+2] == SPECIAL_BND)
                  fprintf(TST, "%s\n", "SPECIAL_BND");
              else {
  					  for (k = MARKER1_BND; k < NUMS_BND; k++) {
						  if (id_prop[j*2+2] == k) {
								fprintf(TST, "MARKER%i_BND\n", k-SPECIAL_BND); 
								break;
						  }
					  }
					  if (k == NUMS_BND) fprintf(TST, "%s\n", "UNKNOWN_BND");
				  }
			 }
			 fclose (TST);
		}
  }
  return(1);
}

template <typename T> 
void CDraft<T>::bar_condit_in(char * ch_CONDIT, int id_print_stat)
{
	unsigned long count, upper_limit;
	char        * id_CONDIT = read_struct_ascii(ch_CONDIT);
	if         (! id_CONDIT) return;
	user_Count   (id_CONDIT, 0, upper_limit, '\x0');
	bar_condit_in(id_CONDIT, count = 0, upper_limit, id_print_stat);
	delete_struct(id_CONDIT);
}

template <typename T> 
void CDraft<T>::bar_condit_in(const char * ch_CONDIT, int id_print_stat)
{
	unsigned long count, upper_limit;
	char        * id_CONDIT = read_struct_ascii(ch_CONDIT);
	if         (! id_CONDIT) return;
	user_Count   (id_CONDIT, 0, upper_limit, '\x0');
	bar_condit_in(id_CONDIT, count = 0, upper_limit, id_print_stat);
	delete_struct(id_CONDIT);
}

/////////////////////////////////////////////////////
//...block structure on the base of the plane sample;
template <typename T> 
int CDraft<T>::ConvertToBeamStruct(int * NN, double * LL, int m_cell, int id_phase)
{
	if (! NN || ! LL || NN[0] <= 0) return(0);

	double pp_cond[] = {0., 0., 0.}, dd, leng;
	int i, j, k, l, m, n, j_shift, N_arcs, Max_phase, shift = bar ? 2 : 0, N0 = N;

//////////////////////////////////////////
//...construction of the blocks structure;
//	for (                           k = 0; k < N0; k++) B[k].bar->bar_invers();
	for (Max_phase = 0,             k = 0; k < N0; k++) if (B[k].link[NUM_PHASE] < Max_phase) Max_phase = B[k].link[NUM_PHASE];
	for (leng = 0., j_shift = 0,    n = 0; n < NN[0]; leng = LL[n], n++)
	for (dd = (LL[n]-leng)/NN[n+1], j = 0; j < NN[n+1]; j_shift++,  j++) if (j || n)
	for (pp_cond[2] = leng+dd*(j+.5),	  k = 0; k < N0;						 k++) {

/////////////////////////////////
//...образование геометрии блока;
		CCells * beam = new CCells;

		beam->get_beam(m_cell ? B[k].bar->bar_cpy(m_cell) : B[k].bar->bar_cpy(), dd);
		//beam->cells_out("bar_prev");

		beam->segms_bar();
		//beam->cells_iso(pp_cond);
		beam->bar_iso(pp_cond);
		//beam->cells_out("bar");

////////////////////////////////////
//...достраиваем торцевые сегменты;
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[7];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { swap(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_SREZKA_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[5];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { swap(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.75);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point1]->mp[1]+B1*cos(f0+M_PI*.75)-A1*sin(f0+M_PI*.75);
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point1]->mp[2]+B1*cos(f0+M_PI*.75)+A1*sin(f0+M_PI*.75);

			if (B[k].type & ERR_CODE) swap(A1, B1);
			if (add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)90.;
			}
			if (add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)90.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK2BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[6];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { swap(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}

/////////////////////////////////
//...образование блока в образце;
		add_block(B[k].type);
		B[k+j_shift*N0].bar = beam;
		set_block3D(B[k+j_shift*N0], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		set_link (B[k+j_shift*N0], max(NUM_PHASE, B[k].link[0]));
		if (j_shift > 1)
		B[k+(j_shift-1)*N0].link[2] = k+(j_shift+0)*N0; 
		B[k+(j_shift+0)*N0].link[1] = k+(j_shift-1)*N0; 
		if (m_cell) {
			for (m = min(NUM_PHASE-2, B[k].link[0]), i = 0; i < m; i++)
				if (B[k].link[i+1] >= 0) B[k+j_shift*N0].link[i+3] = B[k].link[i+1]+j_shift*N0; else
				if (B[k].link[i+1] <= SRF_STATE) {
					B[k+j_shift*N0].link[i+3] = B[k].link[i+1]-shift;	
					for (l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
						if (B[k].link[l+1] == SRF_STATE-B[k].link[i+1]) {
							B[k+j_shift*N0].link[i+3] += shift-j_shift*N0; break;
						}
				}
				else B[k+j_shift*N0].link[i+3] = B[k].link[i+1];
		}
		else {
			for (m = min(NUM_PHASE-2, B[k].bar->arcs_number()), i = 0; i < m; i++)
				if (B[k].link[i+1] >= 0) B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1]+j_shift*N0; else
				if (B[k].link[i+1] <= SRF_STATE) {
					B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1]-shift;
					for (l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
						if (B[k].link[l+1] == SRF_STATE-B[k].link[i+1]) {
							B[k+j_shift*N0].link[m-i+2] += shift-j_shift*N0; break;
						}
				}
				else B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1];
		}
		for (i = NUM_PHASE; i < B[k+j_shift*N0].link[0]; i++) {
			if (B[k].link[i+1] >= 0)	B[k+j_shift*N0].link[i+1] = B[k].link[i+1]+j_shift*N0;
			else								B[k+j_shift*N0].link[i+1] = B[k].link[i+1];
		}
		if (B[k].link[0] < NUM_PHASE)	B[k+j_shift*N0].link[0] = 2+B[k].link[0];
		if (id_phase)						B[k+j_shift*N0].link[NUM_PHASE] = B[k].link[NUM_PHASE]+n*(Max_phase-1);
		else									B[k+j_shift*N0].link[NUM_PHASE] = B[k].link[NUM_PHASE];
	}

///////////////////////////////////////////////////////////////
//...преобразование блоков в нулевом слое и окаймляющего блока;
	for (dd = LL[0]/NN[1], pp_cond[2] = dd*.5, k = 0; k < N0; k++) {
		CCells * beam = new CCells;

		beam->get_beam(m_cell ? B[k].bar->bar_cpy(m_cell) : B[k].bar->bar_cpy(), dd);
		beam->segms_bar();
		beam->bar_iso(pp_cond); /*beam->cells_iso(pp_cond); */N_arcs = m_cell ? B[k].link[0] : B[k].bar->arcs_number(); delete B[k].bar;

////////////////////////////////////
//...достраиваем торцевые сегменты;
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[7];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { swap(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_SREZKA_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[5];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { swap(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.75);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point1]->mp[1]+B1*cos(f0+M_PI*.75)-A1*sin(f0+M_PI*.75);
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point1]->mp[2]+B1*cos(f0+M_PI*.75)+A1*sin(f0+M_PI*.75);

			if (B[k].type & ERR_CODE) swap(A1, B1);
			if (add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)90.;
			}
			if (add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)90.;
			}
		}

//////////////////////////////////////////
//...преобразуем начальный блок в образце;
		B[k].bar = beam;
		set_block3D(B[k], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		int * link = B[k].link; B[k].link = NULL;
		set_link (B[k], max(NUM_PHASE, link[0]));

		if (m_cell) {
			for (m = min(NUM_PHASE-2, N_arcs), i = 0; i < m; i++)
				if (link[i+1] <= SRF_STATE) {
					B[k].link[i+3] = link[i+1]-shift;	
					for (l = NUM_PHASE; l < link[0]; l++) //...коррекция многофазных сред;
						if (link[l+1] == SRF_STATE-link[i+1]) {
							B[k].link[i+3] += shift; break;
						}
				}
				else B[k].link[i+3] = link[i+1];
		}
		else {
			for (m = min(NUM_PHASE-2, N_arcs), i = 0; i < m; i++)
				if (link[i+1] <= SRF_STATE) {
					B[k].link[m-i+2] = link[i+1]-shift;
					for (l = NUM_PHASE; l < link[0]; l++) //...коррекция многофазных сред;
						if (link[l+1] == SRF_STATE-link[i+1]) {
							B[k].link[m-i+2] += shift; break;
						}
				}
				else B[k].link[m-i+2] = link[i+1];
		}
		for (i = NUM_PHASE; i < link[0]; i++)	B[k].link[i+1] = link[i+1];
		if  (		NUM_PHASE > link[0])				B[k].link[0] = 2+link[0];
															B[k].link[NUM_PHASE] = link[NUM_PHASE];
		if (NN[0] > 0 && NN[1] > 1 || NN[0] > 1) 
		B[k].link[2] = k+N0; else B[k].link[2] = ERR_STATE;
		B[k].link[1] = ERR_STATE;
		delete_struct(link);
	}
	if (bar) {
		CCells * beam = new CCells;
		beam->get_beam(m_cell ? bar->bar_cpy(m_cell) : bar->bar_cpy(), leng);
		beam->segms_bar(); 

		delete bar; bar = beam;
	}
	SkeletonRadius(OK_STATE); return(Max_phase);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//...cylindrical symmetry block structure on the base of the plane sample (not complete code);
template <typename T> 
void CDraft<T>::ConvertToCylStruct(int NN, double axis_dist, double fi0, int m_close)
{
	double pp_cond[] = {axis_dist, 0., 0.}, delta_fi = (fi0 = m_close ? 360. : fi0)/(NN = m_close ? max(NN, 2) : NN);
	int k, i, j, N0 = N;

//////////////////////////////////////////
//...construction of the blocks structure;
	for (k = 0; k < N0; k++) { B[k].bar->cells_iso(pp_cond); B[k].bar->bar_invers();}
	for (j = 1; j < NN; j++)
	for (k = 0; k < N0; k++) {

/////////////////////////////////
//...образование геометрии блока;
		CCells * beam = new CCells;

		beam->get_cyl_beam(B[k].bar->bar_cpy(), delta_fi);
		beam->segms_bar();
		beam->cells_iso(NULL, 0., delta_fi*(j+.5)*M_PI/180.);

////////////////////////////////////
//...достраиваем торцевые сегменты;

/////////////////////////////////
//...образование блока в образце;
		add_block(B[k].type);
		B[k+j*N0].bar = beam;
		set_block(B[k+j*N0], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		set_link (B[k+j*N0], max(NUM_PHASE, B[k].link[0]));
		if (j > 1)
		B[k+(j-1)*N0].link[2] = k+(j+0)*N0; 
		B[k+(j+0)*N0].link[1] = k+(j-1)*N0; 
		for (i = 0; i < B[k].link[0]; i++) { 
			if (B[k].link[i+1] >= 0)			B[k+j*N0].link[i+3] = B[k].link[i+1]+j*N0; else
			if (B[k].link[i+1] <= SRF_STATE) B[k+j*N0].link[i+3] = B[k].link[i+1]-2;	
			else										B[k+j*N0].link[i+3] = B[k].link[i+1];
		}
	}

///////////////////////////////////////////////////////////////
//...преобразование блоков в нулевом слое и окаймляющего блока;
	for (k = 0; k < N0; k++) {
		CCells * beam = new CCells;

		beam->get_cyl_beam(B[k].bar->bar_cpy(), delta_fi);
		beam->segms_bar(); 
		beam->cells_iso(NULL, 0., delta_fi*.5*M_PI/180.); delete B[k].bar;

		B[k].bar = beam;
		set_block(B[k], SPHERE_GENUS, 1.);
		for (i = B[k].link[0]; i > 0 ; i--) 
		B[k].link[i+2] = B[k].link[i];
		if (NN > 1) {
			B[k].link[2] = k+N0; 
			if (m_close) {
				B[k].link[1] = k+(NN-1)*N0; 
				B[k+(NN-1)*N0].link[2] = k; 
			}
		}
	}
	CCells * beam = new CCells;
	beam->get_cyl_beam(bar->bar_cpy(), fi0);
	beam->segms_bar(); 
	beam->cells_iso(NULL, 0., fi0*.5*M_PI/180.);

	delete bar; bar = beam;
	SkeletonRadius();
}

//////////////////////////////////////////////////////////////
//...корректировка криволинейных дуг блоков в плоском образце;
template <typename T> 
void CDraft<T>::CurveCorrect(double X0, double Y0, double rad, double eps)
{
	for (int j, m, l, k = 0; k < N; k++) 
	for (j = B[k].link[0]; j > 0; j--) 
	if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
		for (l = 0; l < NUM_PHASE; l++) 
			if (B[k].link[l+1] == -m+SRF_STATE) break;
			if (l < B[k].bar->ce[0]->graph[1]) {
				int m1 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 0),
					 m2 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 1);
				double X1 = B[k].bar->ce[m1]->mp[1], Y1 = B[k].bar->ce[m1]->mp[2],
						 X2 = B[k].bar->ce[m2]->mp[1], Y2 = B[k].bar->ce[m2]->mp[2];
				if (fabs(sqr(X1-X0)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-Y0)-sqr(rad)) < eps) {
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
				}
			}
	}
}

///////////////////////////////////////////////////////
//...идентификация точек внутри системы плоских блоков;
template <typename T> 
int CDraft<T>::Poly_struc_in2D(int & k, double X, double Y, int id_fast, double eps)
{
	if (N > 0) {
////////////////////////////////////
//...проверка ближайшей окрестности;
		if (0 <= k && k < N && B[k].link) {
			if (B[k].bar->ce[0]->in_poly_facet(X, Y, 0., id_fast, eps)) return(1);
			for (int j = 0; j < B[k].link[0]; j++) 
				if (	B[k].link[j+1] >= 0 && B[B[k].link[j+1]].bar->ce[0]->in_poly_facet(X, Y, 0., id_fast, eps)) {
				  k = B[k].link[j+1]; return(1);
				}
		}
///////////////////////////////////
//...проверка всего массива блоков;
		for (int i = 0;  i < N; i++) {
			//B[i].bar->cells_out("CCells");
			if  (B[i].bar->ce[0]->in_poly_facet(X, Y, 0., id_fast, eps)) {
				  k = i; return(1);
			}
		}
	}
	k = ERR_STATE; return(0);
}
////////////////////////////////////////////////////////////////
//...идентификация точек внутри системы пространственных блоков;
template <typename T> 
int CDraft<T>::Poly_struc_in3D(int & k, double X, double Y, double Z, int id_local, double eps)
{
  if (N > 0) {
      k = -1;
      for (int i = 0;  i < N; i++)
      if  (B[i].bar->in_poly_body(X, Y, Z, id_local, NULL, eps)) {
           k = i; return(1);
      }
      return(1);
  }
  k = -1; return(0);
}

////////////////////////////////////////////////////////////////
//...проверка попадания точки внутрь системы сферических блоков;
template <typename T> 
int CDraft<T>::Sph_struc_in3D(int & k, double X, double Y, double Z, int id_map1into)
{
  if (N > 0) {
      int    i  = /*N-1*/0, l = N-1;
      double fx = X-B[i].mp[1], fy = Y-B[i].mp[2],
             fz = Z-B[i].mp[3], f, f0;


      if (B[l].type == ZOOM_BLOCK && B[l].mp[7] < 0. && sqr(X-B[l].mp[1])+sqr(Y-B[l].mp[2])+sqr(Z-B[l].mp[3]) > sqr(B[l].mp[7])) {
          k = N-1;
          return(1);
      }
      if (id_map1into && ((int)B[i].mp[l = size_of_map(B[i].mp)] == SPECIAL_BOX ||
                          (int)B[i].mp[l                       ] ==   BASIS_BOX ||
                          (int)B[i].mp[l                       ] == SPECIAL_SHEET)) {
           if (map1into(B[i].mp, fx, fy, fz)) {
                k = i; return(1);
           }
           for (k = i, i = 1; i < N; i++) {
               fx = X-B[i].mp[1];
               fy = Y-B[i].mp[2];
               fz = Z-B[i].mp[3];
               if (map1into(B[i].mp, fx, fy, fz)) {
                    k = i; return(1);
               }
           }
      }
      else
          for (f = B[i].mp[7] > 0 ? (fx*fx+fy*fy+fz*fz)/sqr(B[i].mp[7]) : 
                                  ((f0 = fx*fx+fy*fy+fz*fz) > 0. ? sqr(B[i].mp[7])/f0 : MAX_HIT), 
               k = i, i = 1;  i < N; i++) {
          fx = X-B[i].mp[1];
          fy = Y-B[i].mp[2];
          fz = Z-B[i].mp[3];
          if (B[i].mp[7] > 0 && (f0 = (fx*fx+fy*fy+fz*fz)/sqr(B[i].mp[7])) < f ||
              B[i].mp[7] < 0 && (f0 = sqr(B[i].mp[7])/(fx*fx+fy*fy+fz*fz)) < f) {
              k = i;
              f = f0;
          }
      }
      return(1);
  }
  k = -1; return(0);
}

////////////////////////////////////////////////
//...корректировка номеров криволинейных блоков;
template <typename T> 
void CDraft<T>::StructEllCorrect(int & k, double X, double Y)
{
	double P[3] = {X, Y, -1.};
	if (k >= 0 && B[k].link[NUM_PHASE] == -1 && bar && bar->ce[0] && stru.maps_in(bar->ce[0]->mp, P)) {
		if (B[k].link[0] > NUM_PHASE && B[k].link[NUM_PHASE+1] >= 0) k = B[k].link[NUM_PHASE+1];
		else { //..проверяем еще один уровень связей;
			int l, i;
			for (l = 0; l < B[k].link[0]; l++)
			if ((i = B[k].link[l+1]) >= 0 && B[i].link[NUM_PHASE] == -1 && 
						B[i].link[0] > NUM_PHASE && B[i].link[NUM_PHASE+1] >= 0) break;
			if ( l < B[k].link[0]) k = B[i].link[NUM_PHASE+1];	
			else						  k = -1;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//...распределение множества точек по номерам блоков для призматической структуры образца;
template <typename T> 
void CDraft<T>::hit_beam_struct(CGrid * nd, int N0, int * NN, double * LL, int axis)
{
	if (! NN || ! LL || NN[0] <= 0) return;

	double dd, leng;

	if (axis == AXIS_X) {
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			int k, m, n, l;
			for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
			for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
			if (leng+m*dd <= nd->Y[j] && nd->Y[j] < leng+(m+1.)*dd+EE_ker) {
				n = NN[0]; break;
			}
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->Z[0], nd->X[i])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
	else
	if (axis == AXIS_Y) {
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			int k, m, n, l;
			for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
			for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
			if (leng+m*dd <= nd->X[i] && nd->X[i] < leng+(m+1.)*dd+EE_ker) {
				n = NN[0]; break;
			}
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->Y[j], nd->Z[0])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
	else
	if (axis == AXIS_Z) {
		int k, m, n, l;
		for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
		for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
		if (leng+m*dd <= nd->Z[0] && nd->Z[0] < leng+(m+1.)*dd+EE_ker) {
			n = NN[0]; break;
		}
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->X[i], nd->Y[j])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
}

/////////////////////////////////////////////
//														 //
//   БИБЛИОТЕКА ГОТОВЫХ БЛОЧНЫХ СТРУКТУР   //
//														 //
/////////////////////////////////////////////
////////////////////////////////////////////////
//...блочная структура прямоугольного гексаэдра;
template <typename T> 
int CDraft<T>::GetBoxStruct(double AA, double BB, double CC, int NX, int NY, int NZ)
{
	double delta_A = AA/(NX = max(NX, 1)), delta_B = BB/(NY = max(NY, 1)), P0[3];
	CCells * ce;

//////////////////////
//...обнуляем образец;
	init_blocks(NULL);

//////////////////////////////////////
//...добавлям описание плоского торца;
	ce = new CCells;
	ce->get_sheet(AA, BB);
	bar = new CCells;	
	bar->bar_add(ce);

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int j = 0; j < NY;  j++)
	for (int i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCells;
		ce->get_sheet(delta_A, delta_B);
		P0[0] = (i+.5)*delta_A-AA*.5;
		P0[1] = (j+.5)*delta_B-BB*.5;
		P0[2] = 0.;
		ce->cells_iso(P0);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block3D(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link   (B[i+j*NX], 4);
		B[i+j*NX].link[1] = j ? i+(j-1)*NX : BOUND4_STATE;
		B[i+j*NX].link[2] = (i+1)%NX ? (i+1)+j*NX : BOUND5_STATE;
		B[i+j*NX].link[3] = (j+1)%NY ? i+(j+1)*NX : BOUND6_STATE;
		B[i+j*NX].link[4] = i ? (i-1)+j*NX : BOUND7_STATE;
	}
	SkeletonRadius(OK_STATE);
///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &CC, OK_STATE);
	return(N0);
}

////////////////////////////////////////////////////////////////////////
//...блочная структура для сферического включения с промежуточным слоем;
template <typename T> 
int CDraft<T>::GetSphBoxStruct(double AA, double BB, double CC, double rad, double ll, int id_bound)
{
	if (rad > AA*.5 || rad > BB*.5 || rad > CC*.5) return(0);
	double P1[3] = {-.5*AA, -.5*BB, -.5*CC},
			 P2[3] = {0.5*AA, -.5*BB, -.5*CC},
			 P3[3] = {0.5*AA, 0.5*BB, -.5*CC},
			 P4[3] = {-.5*AA, 0.5*BB, -.5*CC},
			 P5[3] = {-.5*AA, -.5*BB, 0.5*CC},
			 P6[3] = {0.5*AA, -.5*BB, 0.5*CC},
			 P7[3] = {0.5*AA, 0.5*BB, 0.5*CC},
			 P8[3] = {-.5*AA, 0.5*BB, 0.5*CC};
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks();

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCells;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block3D(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = ll > 0. ? k+2 : k+1;
		B[k].link[NUM_PHASE] = -2;
	}

////////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- параллелепипедная матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCells;
		ce->get_sphere(-rad-ll);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}
	ce = new CCells;
	ce->get_quad_facet(P1, P5, P8, P4, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P2, P3, P7, P6, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P1, P2, P6, P5, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P4, P8, P7, P3, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P1, P4, P3, P2, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P5, P6, P7, P8, OK_STATE);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block3D(B[k], SPHEROID_GENUS, sqrt(AA*AA+BB*BB+CC*CC)*.5);
	B[k].mp[8] = rad+ll;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? (ll > 0. ? k+1 : k-1) : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

////////////////////////////////////////////////
//...добавляем при необходимости межфазный слой;
	if (rad > 0 && ll > 0.) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type = ZOOM_BLOCK;

		ce = new CCells;
		ce->get_sphere(rad+ll);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block3D(B[k], SPHEROID_GENUS, rad+ll);
		B[k].mp[8] = rad;
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k-1;
		B[k].link[2] = k-2;
		B[k].link[NUM_PHASE] = id_bound == OK_STATE ? -2 : -3;
	}
	return(N);
}

//////////////////////////////////////////////////////////////
//...block structure for box sample with spheroidal inclusion;
template <typename T> 
int CDraft<T>::GetSpheroidBoxStruct(double AA, double BB, double CC, double RR, double rr, double fi, double theta, int id_bound)
{
	if (rr > AA*.5 || rr > BB*.5 || RR > CC*.5) return(0);
	double P1[3] = {-.5*AA, -.5*BB, -.5*CC},
			 P2[3] = {0.5*AA, -.5*BB, -.5*CC},
			 P3[3] = {0.5*AA, 0.5*BB, -.5*CC},
			 P4[3] = {-.5*AA, 0.5*BB, -.5*CC},
			 P5[3] = {-.5*AA, -.5*BB, 0.5*CC},
			 P6[3] = {0.5*AA, -.5*BB, 0.5*CC},
			 P7[3] = {0.5*AA, 0.5*BB, 0.5*CC},
			 P8[3] = {-.5*AA, 0.5*BB, 0.5*CC};
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks();

////////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сфероидальное включение;
	if (RR > 0 && rr > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCells;
		ce->get_spheroid(RR, rr);
		ce->cells_iso(NULL, fi, theta);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block3D(B[k], SPHEROID_GENUS, RR, 0, 0., 0., 0., fi, theta);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

//////////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- параллелепипедная матрица со сфероидальным включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;
	if (RR > 0 && rr > 0.) {
		ce = new CCells;
		ce->get_spheroid(-RR, -rr);
		ce->cells_iso(NULL, fi, theta);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}
	ce = new CCells;
	ce->get_quad_facet(P1, P5, P8, P4, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P2, P3, P7, P6, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P1, P2, P6, P5, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P4, P8, P7, P3, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P1, P4, P3, P2, OK_STATE);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_quad_facet(P5, P6, P7, P8, OK_STATE);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block3D(B[k], SPHEROID_GENUS, sqrt(AA*AA+BB*BB+CC*CC)*.5*get_param(NUM_MPLS+1));
	B[k].mp[8] = RR;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = RR > 0. && rr > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////////////////////////////////////////
//...блочная структура для 3D прямоугольной решетки (прямое моделирование);
template <typename T> 
int CDraft<T>::GetLatticeBox3DStruct(double * X, double * Y, double * Z, int NX, int NY, int NZ)
{
	if (! X || ! Y || ! Z || NX < 2 || NY < 2 || NZ < 2) return(0);

//////////////////////
//...обнуляем образец;
	int mX, mY, mZ, k = -1;
	init_blocks();

///////////////////////////////////////////////////////////////////
//...образуем описание блочной структуры параллелепипедной решетки;
	for (mZ = 1; mZ < NZ; mZ++)
	for (mY = 1; mY < NY; mY++)
	for (mX = 1; mX < NX; mX++) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

///////////////////////////////////////////////////////
//...образуем описание граничных точек параллелепипеда;
		double P[24];
		P[0]  = X[mX-1]; P[1]  = Y[mY-1]; P[2]  = Z[mZ-1];
		P[3]  = X[mX];   P[4]  = Y[mY-1]; P[5]  = Z[mZ-1];
		P[6]  = X[mX];   P[7]  = Y[mY];   P[8]  = Z[mZ-1];
		P[9]  = X[mX-1]; P[10] = Y[mY];   P[11] = Z[mZ-1];
		P[12] = X[mX-1]; P[13] = Y[mY-1]; P[14] = Z[mZ];
		P[15] = X[mX];   P[16] = Y[mY-1]; P[17] = Z[mZ];
		P[18] = X[mX];   P[19] = Y[mY];   P[20] = Z[mZ];
		P[21] = X[mX-1]; P[22] = Y[mY];   P[23] = Z[mZ];

//////////////////////////////////////////////
//...образуем описание граней параллелепипеда;
		CCells * ce = new CCells;
		ce->get_quad_facet(P, P+12, P+21, P+9, OK_STATE);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_quad_facet(P+3, P+6, P+18, P+15, OK_STATE);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_quad_facet(P, P+3, P+15, P+12, OK_STATE);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_quad_facet(P+9, P+21, P+18, P+6, OK_STATE);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_quad_facet(P, P+9, P+6, P+3, OK_STATE);
		B[k].bar->bar_add(ce);

		ce = new CCells;
		ce->get_quad_facet(P+12, P+15, P+18, P+21, OK_STATE);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();

//////////////////////////////////////////
//...образуем описание блока и его связей;
		set_block3D(B[k], SPHEROID_GENUS, sqrt(sqr(X[mX]-X[mX-1])+sqr(Y[mY]-Y[mY-1])+sqr(Z[mZ]-Z[mZ-1]))*.5, 0,
					  (X[mX]+X[mX-1])*.5, (Y[mY]+Y[mY-1])*.5, (Z[mZ]+Z[mZ-1])*.5);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = mX > 1    ?				  -1+k : /*NX-2+k*/BOUND2_STATE;
		B[k].link[2] = mX < NX-1 ?					1+k :/*-NX+2+k*/BOUND2_STATE;
		B[k].link[3] = mY > 1    ?			  -NX+1+k : /*(NY-2)*(NX-1)+k*/BOUND2_STATE;
		B[k].link[4] = mY < NY-1 ?				NX-1+k :/*-(NY-2)*(NX-1)+k*/BOUND2_STATE;
		B[k].link[5] = mZ > 1    ?-(NY-1)*(NX-1)+k : /*(NZ-2)*(NY-1)*(NX-1)+k*/BOUND3_STATE;
		B[k].link[6] = mZ < NZ-1 ? (NY-1)*(NX-1)+k :/*-(NZ-2)*(NY-1)*(NX-1)+k*/BOUND1_STATE;
		B[k].link[NUM_PHASE] = -1;
	}
	return(N);
}

////////////////////////////////////////////////////////////////////////////////////////
//...блочная структура для 3D прямоугольной решетки (моделирование на уровне геометрии);
template <typename T> 
int CDraft<T>::GetLatticeBox3DStruct(double * X, double * Y, double * Z, int NX, int NY, int NZ, Num_Block id_block)
{
	if (! X || ! Y || ! Z || NX < 2 || NY < 2 || NZ < 2) return(0);
	int mX, mY, mZ, k, l, N_geom = 0, m = 0, j = 1;

///////////////////////////////////////////////////////
//...инициализируем и просчитываем геометрию элементов;
	stru.zero_grid();
	for (mZ = 1; mZ < NZ; mZ++)
	for (mY = 1; mY < NY; mY++)
	for (mX = 1; mX < NX; mX++) stru.geom_ptr_add(11, N_geom);


//////////////////////////////////////////
//...заполнякм массив геометрии элементов;
	if ((stru.geom = (int *)new_struct((stru.geom_ptr[N_geom]+1)*sizeof(int))) != NULL) {
		for (mZ = 0; mZ < NZ-1; mZ++)
		for (mY = 0; mY < NY-1; mY++)
		for (mX = 0; mX < NX-1; mX++) {
			l = stru.geom_ptr[stru.geom[0]++];
			stru.geom[l]   = GL_BOXS;
			stru.geom[l+1] = 11;
			stru.geom[l+2]  = -1;
			stru.geom[l+3]  = -(++m);
			stru.geom[l+4]  = -j;
			stru.geom[l+5]  = (k = mX+NX*(mY+NY*mZ));
			stru.geom[l+6]  = k+1;
			stru.geom[l+7]  = k+1+NX;
			stru.geom[l+8]  = k+NX;
			stru.geom[l+9]  = (k += NX*NY);
			stru.geom[l+10]  = k+1;
			stru.geom[l+11] = k+1+NX;
			stru.geom[l+12] = k+NX;
		}

////////////////////////////////////////////////////
//...образуем формальное описание граничных условий;
		stru.cond = (int *)new_struct(2*sizeof(int));
		stru.cond_ptr = (int *)new_struct(2*sizeof(int));
		stru.cond_ptr[0] = 1;

////////////////////////////////////////////////////////////////////////
//...заполняем массив граничных точек - вершин параллелепипедных блоков;
		k = 0;
		for (mZ = 0; mZ < NZ; mZ++)
		for (mY = 0; mY < NY; mY++)
		for (mX = 0; mX < NX; mX++) {
			stru.add_new_point(X[mX], Y[mY], Z[mZ]); 
			stru.hit[k] = k+1; ++k;
		}

/////////////////////////////////////////////
//...устанавливаем блоки и линкуем структуру;
		LinkUniStruct();
		SetBUniStruct(id_block, SPHEROID_GENUS);
	}
	return(N_geom);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...block structure for Ono Box sample on the base of extern description block structure;
template <typename T> 
int CDraft<T>::GetOnoBoxStruct(double AA, double BB, double aa, double bb, double rad, int * NN, double * LL, int id_phase, int id_curve, double eps)
{
////////////////////////////////////////////
//...корректируем криволинейные дуги блоков;
	if (id_curve) {
		for (int j, m, l, k = 0; k < N; k++) 
		for (j = B[k].link[0]; j > 0; j--) 
		if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
			for (l = 0; l < NUM_PHASE; l++) 
				if (B[k].link[l+1] == -m+SRF_STATE) break;
				if (l < B[k].bar->ce[0]->graph[1]) {
					int m1 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 0),
						 m2 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 1);
					double X1 = B[k].bar->ce[m1]->mp[1], Y1 = B[k].bar->ce[m1]->mp[2],
							 X2 = B[k].bar->ce[m2]->mp[1], Y2 = B[k].bar->ce[m2]->mp[2], 
							 X0 = (AA-aa)*.5+rad, Y0 = (BB-bb)*.5+rad, 
							 XX = (AA+aa)*.5-rad, YY = (BB+bb)*.5-rad;
					if (fabs(sqr(X1-X0)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-Y0)-sqr(rad)) < eps ||
						 fabs(sqr(X1-X0)+sqr(Y1-YY)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-YY)-sqr(rad)) < eps ||
						 fabs(sqr(X1-XX)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-XX)+sqr(Y2-Y0)-sqr(rad)) < eps ||
						 fabs(sqr(X1-XX)+sqr(Y1-YY)-sqr(rad)) < eps && fabs(sqr(X2-XX)+sqr(Y2-YY)-sqr(rad)) < eps) {
						B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
						B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
					}
				}
		}
	}

////////////////////////////////////////////////////////////////////////////
//...преобразуем геометрию образца и добавлям описание структуры излучателя;
	int N0 = N, 
		 N_phase = ConvertToBeamStruct(NN, LL, OK_STATE, id_phase);
	if (id_curve) add_cyl_surface(rad, 0., 0., 0.);

///////////////////////////////////////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем торце и корректируем связи (убираем включение);
	int  k;
	for (k = 0; k < N; k++) {
		if (B[k].link[1] <= ERR_STATE && B[k].link[NUM_PHASE] == -2) B[k].link[1] = BOUND1_STATE;
		if (B[k].link[2] <= ERR_STATE)	B[k].link[2] = BOUND3_STATE;

		for (int m = min(NUM_PHASE, B[k].link[0]), i = 0; i < m; i++) if (B[k].link[i+1] <= SRF_STATE)
		for (int l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
		if (B[k].link[i+1] == SRF_STATE-B[k].link[l+1]) {
			 B[k].link[i+1] = B[k].link[l+1]; break;
		}
		B[k].link[0] = NUM_PHASE;
		B[k].link[NUM_PHASE] = B[k].link[NUM_PHASE]/(1-N_phase)-1;
	}

////////////////////////////////////////////////////////////////////////////////
//...перемещение линков в конец списка (криволинейную границу не устанавливаем);
	for (k = 0; k < N; k++)
	for (int i, j = B[k].link[0]; j > 0; j--) 
	if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE]) {
		B[k].link[j] =-i+SRF_STATE;
		link_add(B[k], i);
		link_add(B[k], ERR_STATE, NULL_STATE);
	}

////////////////////////////////////
//...коррекция нескольких включений;
	if (id_phase == OK_STATE)
	for (k = 0; k < N; k++)
	  if (B[k].link[NUM_PHASE] <= -3)
			B[k].link[NUM_PHASE]  = -1;

	return(N0);
}

/////////////////////////////////////////////////////////////
//...строковые варианты функции (зачитывание внешнего файла);
template <typename T> 
int CDraft<T>::GetOnoBoxStruct(double AA, double BB, double aa, double bb, double rad, int * NN, double * LL, char * name, int id_phase, int id_curve, double eps)
{
////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");
	int	 N0 = GetOnoBoxStruct(AA, BB, aa, bb, rad, NN, LL, id_phase, id_curve, eps);
	return(N0);
}
template <typename T> 
int CDraft<T>::GetOnoBoxStruct(double AA, double BB, double aa, double bb, double rad, int * NN, double * LL, const char * name, int id_phase, int id_curve, double eps)
{
////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");
	int	 N0 = GetOnoBoxStruct(AA, BB, aa, bb, rad, NN, LL, id_phase, id_curve, eps);
	return(N0);
}

///////////////////////////////////////////////////////////////////
//...блочная структура для взаимопроникающих сферических включений;
template <typename T> 
int CDraft<T>::GetPenetrateSphere(double rad, double L)
{
	if (rad < L || rad > M_SQRT2*L) return(0);
	double P1[3] = {-L, 0., 0.},
			 P2[3] = {0., -L, 0.},
			 P3[3] = {0., 0., -L}, RR = sqrt(rad*rad-L*L);
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks();

/////////////////////////////////////////////////////////////////////////////////////
//...образуем описание одного блока -- сферическое включение с проникновением в кубе;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = ZOOM_BLOCK;

////////////////////////////////////////////////////////////////
//...образуем геометрию сферического включения с проникновением;
	ce = new CCells;
	ce->get_sph_intrusion(-rad, L);
	B[k].bar->bar_add(ce);

/////////////////////////////////////////////////
//...образуем геометрию боковых граней с дырками;
	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P1, 0., -M_PI_2); P1[0] = -P1[0];
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P1, 0., M_PI_2);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P2, M_PI_2, -M_PI_2); P2[1] = -P2[1];
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P2, M_PI_2, M_PI_2);
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P3, 0., M_PI); P3[2] = -P3[2];
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P3, 0., 0.);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block3D(B[k], SPHEROID_GENUS, sqrt(3.)*L*get_param(1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

/////////////////////////////////////////
//...block structure for cylinder sample;
template <typename T> 
int CDraft<T>::GetCylinderStruct(double rad, double R, double length, int NX, int NY, int NZ)
{
	double delta_fi = M_PI*2./(NX = max(NX, 2)), delta_rr = (R-rad)/(NY = max(NY, 1));
	CCells * ce; 

//////////////////////
//...обнуляем образец;
   init_blocks(NULL);

//////////////////////////////////////
//...добавлям описание плоского торца;
	ce = new CCells;
	ce->get_ring(rad, R);
	bar = new CCells;	
	bar->bar_add(ce);

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int j = 0; j < NY;  j++)
	for (int i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCells;
		ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, rad+j*delta_rr, rad+(j+1)*delta_rr);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link (B[i+j*NX], 4);
		B[i+j*NX].link[1] = (i-1+NX)%NX+j*NX;
		B[i+j*NX].link[2] = (j+1)%NY ? i+(j+1)*NX : SRF_STATE-1;
		B[i+j*NX].link[3] = (i+1)%NX+j*NX;
		B[i+j*NX].link[4] = j ? i+(j-1)*NX : SRF_STATE;
	}
	SkeletonRadius(OK_STATE);
///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	return(N0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//...block structure for cylinder sample on the base of the extern description of the plane structure;
template <typename T> 
int CDraft<T>::GetCylinderChannelStruct(double rad, double R, double length, int NZ, int id_curve)
{
////////////////////////////////////////////
//...корректируем криволинейные дуги блоков;
	if (id_curve) {
		int k, j, m, l;
		for (k = 0; k < N; k++)
		for (j = B[k].link[0]; j > 0; j--)
		if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
			for (l = 0; l < NUM_PHASE; l++)
				if (B[k].link[l+1] == -m+SRF_STATE) break;
				if (l < B[k].bar->ce[0]->graph[1]) {
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
				}
		}
		for (k = 0; k < N; k++) if (B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph)
		for (j = min(B[k].bar->ce[0]->graph[1], NUM_PHASE); j > 0; j--)
		if (B[k].link[j] == ERR_STATE) {
			B[k].bar->ce[B[k].bar->ce[0]->graph[j+1]]->cells_to(SPHERE_GENUS);
			B[k].bar->ce[B[k].bar->ce[0]->graph[j+1]]->circ_correct(R, 1, NULL_STATE);
		}
	}

////////////////////////////////////////////////////////////////////////////
//...преобразуем геометрию образца и добавлям описание структуры излучателя;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	add_cyl_surface(rad, 0., 0., 0.);

///////////////////////////////////////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем торце и корректируем связи (убираем включение);
	for (int k = 0; k < N; k++) {
		if (B[k].link[NUM_PHASE] == -1 && B[k].link[1] <= ERR_STATE) B[k].link[1] = BOUND2_STATE;
		if (B[k].link[NUM_PHASE] == -2 && B[k].link[1] <= ERR_STATE) B[k].link[1] = BOUND1_STATE;

		for (int m = min(NUM_PHASE, B[k].link[0]), i = 0; i < m; i++) if (B[k].link[i+1] <= SRF_STATE)
		for (int l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
		if (B[k].link[i+1] == SRF_STATE-B[k].link[l+1]) {
			 B[k].link[i+1] = B[k].link[l+1]; break;
		}
		B[k].link[0] = B[k].bar->arcs_number();
	}
	return(N0);
}

/////////////////////////////////////////////////////////////
//...строковые варианты функции (зачитывание внешнего файла);
template <typename T> 
int CDraft<T>::GetCylinderChannelStruct(double rad, double R, double length, char * name, int NZ, int id_curve)
{
	if (rad >= R) return(0); 

////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");
	int	 N0 = GetCylinderChannelStruct(rad, R, length, name, NZ, id_curve);
	return(N0);
}
template <typename T> 
int CDraft<T>::GetCylinderChannelStruct(double rad, double R, double length, const char * name, int NZ, int id_curve)
{
	if (rad >= R) return(0); 

////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");
	int	 N0 = GetCylinderChannelStruct(rad, R, length, name, NZ, id_curve);
	return(N0);
}

/////////////////////////////////////////
//...block structure for cylinder sample;
template <typename T> 
int CDraft<T>::GetCylinderChannelStruct(double rad, double R, double length, double r0, int NX, int NY, int N2, int NZ)
{
	if (rad >= R) return(0); 

	double delta_rr = (R-rad )/(NY = max(NY, 1)), delta_fi = M_PI*2./(NX = max(NX, 2)), 
			 delta_r0 = (rad-r0)/(N2 = max(N2, 1)), RR = rad;
	CCells * ce;

//////////////////////
//...обнуляем образец;
	init_blocks(NULL);

////////////////////////////////////////////
//...добавлям описание структуры излучателя;
	ce = new CCells;
	ce->get_ring(rad, R);
	bar = new CCells;
	bar->bar_add(ce);

//	ce = new CCells;
//	ce->get_disk(rad);
//	bar = new CCells;
//	bar->bar_add(ce);

//////////////////////////////////////////////////////
//...образуем блочную структуру стенки плоского торца;
	int i, j;
	for (j = 0; j < NY;  j++)
	for (i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCells;
		ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, rad+j*delta_rr, rad+(j+1)*delta_rr);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link (B[i+j*NX], 4);
		B[i+j*NX].link[1] = (i-1+NX)%NX+j*NX;
		B[i+j*NX].link[2] = (j+1)%NY ? i+(j+1)*NX : BOUND3_STATE;
		B[i+j*NX].link[3] = (i+1)%NX+j*NX;
		B[i+j*NX].link[4] = j ? i+(j-1)*NX : SRF_STATE;
	}

//////////////////////////////////////////////////////////
//...образуем блочную структуру излучателя плоского торца;
	int shift = 0, m0 = NX*NY;
	if (r0 < rad) {
		shift = m0; m0 += NX*N2; RR = r0;
		for (j = 0; j < N2;  j++)
		for (i = 0; i < NX;  i++) {
			add_block (NULL_BLOCK);
			ce = new CCells;
			ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, r0+j*delta_r0, r0+(j+1)*delta_r0);
			B[shift+i+j*NX].bar = new CCells;
			B[shift+i+j*NX].bar->bar_add(ce);
			set_block(B[shift+i+j*NX], SPHERE_GENUS, 1.);
			set_link (B[shift+i+j*NX], 4);
			B[shift+i+j*NX].link[1] = shift+(i-1+NX)%NX+j*NX;
			if ((j+1)%N2) B[shift+i+j*NX].link[2] = shift+i+(j+1)*NX;
			else {
				B[shift+i+j*NX].link[2] = i;
				B[i].link[4] = shift+i+j*NX;
			}
			B[shift+i+j*NX].link[3] = shift+(i+1)%NX+j*NX;
			B[shift+i+j*NX].link[4] = j ? shift+i+(j-1)*NX : SRF_STATE;
		}
	}

////////////////////////////////////////////////
//...образуем центральный блок на плоском торце;
	add_block (NULL_BLOCK);
	ce = new CCells;
	for (i = 0; i < NX;  i++) {
		CCells * arc = new CCells;
		arc->get_arc(RR, i*delta_fi, (i+1)*delta_fi);
		ce->bar_add(arc);
	}
	ce->bar_invers();
	ce->bar_generate();
	CMap * mp = get_map(2, NULL_GENUS, RING_SEGMENT);
	if (mp) {
		mp[++(i = size_of_map(2, NULL_GENUS))] = (CMap)M_PI*2.;
		mp[++i] = (CMap)0.;
		mp[++i] = (CMap)RR;
	}
	ce->bar_span(mp);

	B[m0].bar = new CCells;
	B[m0].bar->bar_add(ce);
	set_block(B[m0], SPHERE_GENUS, 1.);
	set_link (B[m0], NX);

	for (i = 0; i < NX; i++) {
		B[m0].link[i+1] = shift+i;
		B[shift+i].link[4] = m0;
	}
	SkeletonRadius(OK_STATE);

///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);

//////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем плоском торце;
	for (j = 0; j < NY;  j++)
	for (i = 0; i < NX;  i++) 
		B[i+j*NX].link[1] = BOUND2_STATE;

	if (r0 < rad)
	for (j = 0; j < N2;  j++)
	for (i = 0; i < NX;  i++) 
		B[shift+i+j*NX].link[1] = BOUND1_STATE;

	B[m0].link[1] = BOUND1_STATE;

	return(N0);
}

///////////////////////////////////////////////////////////////////////////////////
//...block structure for cylinder sample with spherical inclusion (equal distance);
template <typename T> 
int CDraft<T>::GetCylSphStruct(double ff_vol, double rad, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = pow(4./(3.*(1.-ff_vol)), 1./3.)*rad;
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCells;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

/////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- цилиндрическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCells;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}

///////////////////////////////////////////////
//...добавляем два торца и боковую поверхность;
	ce = new CCells;
	ce->get_ring_segment(-M_PI, M_PI, R1, 0.);
	ce->mp[3] = -R1*.5;
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_ring_segment(-M_PI, M_PI, 0., R1);
	ce->mp[3] = R1*.5;
	B[k].bar->bar_add(ce);

	ce = new CCells;
	//ce->get_cyl_segment(R1, -M_PI, M_PI, -R1, R1);
	ce->get_cylinder(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, M_SQRT2*R1*get_param(NUM_MPLS+1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////////////////////////////////
//...block structure for cylinder sample with spherical inclusion;
template <typename T> 
int CDraft<T>::GetCylSphStruct(double ff_vol, double rad, double L, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = sqrt(rad/(3.*L*(1.-ff_vol)))*2.*rad;
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCells;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

/////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- цилиндрическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCells;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}

///////////////////////////////////////////////
//...добавляем два торца и боковую поверхность;
	ce = new CCells;
	ce->get_ring_segment(-M_PI, M_PI, R1, 0.);
	ce->mp[3] = -L*.5;
	B[k].bar->bar_add(ce);

	ce = new CCells;
	ce->get_ring_segment(-M_PI, M_PI, 0., R1);
	ce->mp[3] = L*.5;
	B[k].bar->bar_add(ce);

	ce = new CCells;
	//ce->get_cyl_segment(R1, -M_PI, M_PI, -L*.5, L*.5);
	ce->get_cylinder(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, sqrt(R1*R1+L*L*.25)*get_param(NUM_MPLS+1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////////////////////////////////
//...block structure for spherical sample with spherical inclusion;
template <typename T> 
int CDraft<T>::GetSph2BodyStruct(double ff_vol, double rad, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = pow(1./(1.-ff_vol), 1./3.)*rad;
	CCells * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCells;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

//////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- сферическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCells;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}
	ce = new CCells;
	ce->get_sphere(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, R1);
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////
//...block structure for Bars25 sample;
template <typename T> 
int CDraft<T>::GetBars25Struct(CDraft<T> * bars25, double length, int NZ)
{
//////////////////////
//...обнуляем образец;
   init_blocks();

//////////////////////////////////////
//...добавлям описание плоского торца;
	bar = bars25->bar; bars25->bar = NULL;

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int k = 0; k < bars25->N;  k++) {
		add_block (bars25->B[k].type);
		B[k].bar = bars25->B[k].bar; bars25->B[k].bar = NULL;
		B[k].bar->cells_common(bars25->B[k].mp);
		set_block(B[k], SPHERE_GENUS, 1.);
		set_link (B[k], bars25->B[k].link[0]);
		memcpy(B[k].link, bars25->B[k].link, (bars25->B[k].link[0]+1)*sizeof(Topo)); 
	}
	SkeletonRadius(OK_STATE);

///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	return(N0);
}
#endif
