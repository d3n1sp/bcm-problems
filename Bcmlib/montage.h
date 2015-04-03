/*============================================*/
/*                  MONTAGE                   */
/*============================================*/
#ifndef ___MONTAGE___
#define ___MONTAGE___

#include "cprofile.h"

#define CIRCUIT(NWires) (NWires <= 3 ? 1 :\
								(NWires <= 6 ? 2 :\
								(NWires <= 9 ? 3 : -1)))

///////////////////////////////////////////
//...description of the current anker mode;
enum Num_Anker_mode {
			 NULL_MODE = 0,
			GROZO_MODE,
	  OPERATING_MODE,
				P1_MODE,
	 MIN_TEMPER_MODE,
	 MAX_TEMPER_MODE,
				P3_MODE,
 P6_RIGHT_RUSH_MODE,
  P6_LEFT_RUSH_MODE,
 P7_RIGHT_RUSH_MODE,
  P7_LEFT_RUSH_MODE,
		 MONTAGE_MODE,
	  NUM_ANKER_MODE,
};

//////////////////////////////////////////
//...математика для расчета стрел провеса;
void transit     (double H0, double p0, double t0, double EF0, double & H, double p, double t, double EF, double alpha, double L, double creep = 0., double creep0 = 0.);
void transit_p	  (double H0, double p0, double t0, double EF0, double H, double & p, double t, double EF, double alpha, double L, double creep = 0., double creep0 = 0.);
void transit_d   (double H0, double p0, double EF0, double & H, double L, double delta_L);
void transit_H	  (double H0, double p0, double t0, double EF0, double & H, double p, double t, double EF, double alpha, double L, double delta_h);
void strela (double p, double L, double DA, double H, double & f);
void strela0(double p, double L, double DA, double H, double & f);
void tension(double p, double L, double DA, double & H, double & T, int m = 1);
double H_tension(double T,  double p,  double L,  double delta_h, int hyperb = 1);
double L_transit(double H0, double p0, double t0, double EF0, double p, double t, double alpha);
double L_average(int i1, int i2, int id_set_table = 1);
double H_average(Table * tab_mod, Table * tab_clm, int id_modetab, int i1, int i2, int mode = MAX_TEMPER_MODE, double * Te = NULL, double * leng = NULL, double * hh = NULL, double * ff = NULL);
double H_average(int num_modetab, int i1, int i2, int mode = MAX_TEMPER_MODE, double * Te = NULL, double * leng = NULL, double * hh = NULL, double * ff = NULL);
double grozo_strela(int num_modetab, int i, int i_grozo, int * i_wire, double del_izol, double fi0 = M_PI/9.);
double interpolation4(const double Kp_table[][4], int M, double parm, int m_index, int round = OK_STATE);
double interpolatchk2(const double Kp_table[][2], int M, double parm, int check, int round = OK_STATE);
double interpolation2(const double Kp_table[][2], int M, double parm, int round = OK_STATE);

/////////////////////////////////////////////////////
//...уточненная математика для расчета стрел провеса;
void   bottom (double L, double del_h, double p, double H, double * Pf);
double strela0(double L, double del_h, double p, double H, double * Pf);
double sag    (double X, double p, double H, double * Pf);
double tension(double X, double p, double H, double * Pf);

////////////////////////////////////////////////
//...дифференциальная для расчета стрел провеса;
void save_anker (int i1, int i2, int id_wire, double * H, int m = NULL_STATE);
void set_loading(int i1, int i2, double * p, double value);
void transit_anker     (int i1, int i2, int id_wire, double * p0, double t0, double EF0, 
												double * p, double t, double EF, double alpha, double eps = 1e-15);
void transit_anker_test(int i1, int i2, int id_wire, double * p0, double t0, double EF0, 
												double * p, double t, double EF, double alpha, double eps = 1e-15);

////////////////////////////////////
//...аварийный режим обрыва провода;
void transit_damage_anker(int i1, int i2, int id, int id_wire, double p, double EF, double eps = 1e-15);

////////////////////////////////////////////////////////////////////////////////////
// ...установка параметров подвески и режимов нагружения провода и кабеля в пролете;
void set_hanging_index(int i1, int i2, int num_hangingtab, int id_second = 1, int id_second2 = 1);
void set_mode_index   (int i1, int i2, int num_modetab);
void set_mode_loading (int num_modetab);
void set_pue_loading  (int i_ank = 1);
void set_pue_wire		 (int i_ank = 1);
void set_pue_grozo	 (int i_ank = 1);
void set_pue_cable	 (int i_ank = 1);

////////////////////////////////////////////////////////////////////////
// ...извлечение температурных и нагрузочных режимов (провода и кабеля);
void get_temper_mode(double & t_glaze, double & t_rush, double & t_aver, double & t_abs_max, double & t_abs_min, double & t_mtg);
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, 
											 double & T_abs_max, double & T_opt, double & T_max, double & T_min, 
											 double & EF_creep,  double & EF, double & EF_fin, double & alpha, double & creep, char *& mark);
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, 
											 double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, 
											 double & EF_creep, double & EF, double & EF_fin, double & alpha, double & creep);
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, double & EF, double & EF_fin, double & alpha, double & creep, char *& mark);
void get_loading(int num_modetab, double & p1, double & p3, double & p7, double & T_abs_max, double & T_opt, double & EF, double & EF_fin, double & alpha, double & creep);
void get_loading(int num_modetab, double & p1, double & p7, double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, double & EF, double & EF_fin, double & alpha, double & creep);
void get_loading(int num_modetab, double & p1, double & EF, double & alpha, char *& mark);
void get_loading(int num_modetab, char *& mark);
void set_montage_table(int num_montagetab, int i_ank, double H, int mode, double * t_op = NULL);
void set_mode_tension (int num_montagetab, int i_ank, double & pV, double & pH, int mode, double * t_op = NULL);

///////////////////////////////
// ...монтаж анкерного участка;
void AnkerMontage(int i_ank, double S_stand, double S_anker, double S_sect);
void AnkerEarth  (int i_ank, int i, double & H_cab, double & H_wr1, double & S_cable, double & S_wire, int mode);
void AnkerSection(int i_ank, int i, double & H_cab, double & H_wr1, double & S_cable, double & S_wire, int N_sect, int mode);
void AnkerStrela (int i_ank, int i, double & H_cab, double & H_gro, double & H_wr, double & f_cable, double & f_grozo, double & f_wire, int mode);
void AnkerPhasa  (int i_ank, int i, double & H_cab, double & H_wr1, double & fw_cable, double & fw_wire, double & Sw_cable, int mode);
void AnkerExtremeTension(int i_ank, double T_cable = 0., double T_grozo = 0., double T_wire = 0.);
void SetMontageState(int i1, int i2, int i_ank, int mode = 0, double D_phase = 1., int * anker_num = NULL);

////////////////////////////////////////////////////////////////////////////////////////
// ...визуализация линии в первый момент после загрузки (и заполнение монтажных таблиц);
void FirstVisualization(int id_montage_table = 0);

///////////////////////////////////////////////////////////////
// ...сообщение об установленном монтажном состоянии в пролете;
void MessageMontageState(int mode);

///////////////////////////////////////////////////////////////////
// ...рассчет монтажных таблиц по классическому уравнению перехода;
void CalculateCableMontageTable(int i_ank);
void CalculateGrozoMontageTable(int i_ank);
void CalculateWireMontageTable (int i_ank);

//////////////////////////////////////////////////////
//...установка предельно возможного тяжения в пролете;
double MaxTension (Table * tab_mod, Table * tab_clm, int i1, int i2, int mode = MONTAGE_MODE, int max_tension = 1, Table * tab_montage = NULL);
double MaxTension (int num_modetab, int i1, int i2, int mode = MONTAGE_MODE, int max_tension = 1);

double OptTension (int num_modetab, int i1, int i2, int mode = MONTAGE_MODE);
double L_Transit  (int num_modetab);
double ModeTension(Table * tab_mod, Table * tab_clm, int i1, int i2, int mode_ini, int mode_fin, double Te);

void creep_transit(Table * tab_mod, Table * tab_clm, double H, double p, double t, double & H_norm, double t_norm, double L, int m);
void creep_transit(int num_modetab, double H, double p, double t, double & H_norm, double t_norm, double L, int m = 1);

//////////////////////////////////////////////////////////
//...выравнивание опор согласно данным линейных измерений;
void TowerJump(int i1, int i2, int id_cross = 0);
void TowerJump2(int i1, int i2, int id_cross = 0);

//////////////////////////////////////////////////
//...вычисление запаса до конкретного пересечения;
double SectCableVGap(int i, int mode = 0);
double SectWireVGap (int i, int mode = 0, int id_cross = 0);

//////////////////////////////////////////////////////
//...округление данных до нужного числа значащих цифр;
void Rounding (double & S1_cab, double & S3_cab, double & S6_cab, double & S7_cab,
				   double & T1_cab, double & T3_cab, double & T6_cab, double & T7_cab);
void RoundingS(double & S1_cab, double & S3_cab, double & S6_cab, double & S7_cab);
void RoundingT(double & T1_cab, double & T3_cab, double & T6_cab, double & T7_cab);

///////////////////////////////////////////////////////////////
//...основной пакет документов для монтажа линии в CSV формате;
void ConvertToTowerSpanTableCSV    (const char * montage_file = NULL, int id_context = 0);

void ConvertToCableMontageTableCSV (const char * montage_file = NULL, int id_context = 1);
void ConvertToGrozoMontageTableCSV (const char * montage_file = NULL, int id_context = 2);
void ConvertToWireMontageTableCSV  (const char * montage_file = NULL, int id_context = 3);

void ConvertToDampherTableCSV		  (const char * montage_file = NULL, int id_context = 4);

void ConvertToCableSpanLoadTableCSV(const char * montage_file = NULL, int i = -1, int id_context = 5);
void ConvertToGrozoSpanLoadTableCSV(const char * montage_file = NULL, int i = -1, int id_context = 6);
void ConvertToWireSpanLoadTableCSV (const char * montage_file = NULL, int i = -1, int id_context = 7);

/////////////////////////////////////////////////////////////////////
//...дополнительный пакет документов для монтажа линии в CSV формате;
void GetBMMBendingTableCSV			  (int i_ank, int i1, const char * montage_file = NULL);
void ConvertToBendingTableCSV		  (const char * montage_file = NULL, int id_context = 5);

void ConvertToCableLoadingTableCSV (const char * montage_file = NULL);
void ConvertToGrozoLoadingTableCSV (const char * montage_file = NULL);
void ConvertToWireLoadingTableCSV  (const char * montage_file = NULL);

void ConvertToCableHangingTableCSV (const char * montage_file = NULL);
void ConvertToGrozoHangingTableCSV (const char * montage_file = NULL);
void ConvertToWireHangingTableCSV  (const char * montage_file = NULL);

void ConvertToCableSizeTableCSV	  (const char * montage_file = NULL);
void ConvertToGrozoSizeTableCSV    (const char * montage_file = NULL);
void ConvertToWireSizeTableCSV     (const char * montage_file = NULL);

void ConvertToSectCableSizeTableCSV(const char * montage_file = NULL);
void ConvertToSectGrozoSizeTableCSV(const char * montage_file = NULL);
void ConvertToSectWireSizeTableCSV (const char * montage_file = NULL);
 
void ConvertToGrozoDistTableCSV    (const char * montage_file = NULL);


/////////////////////////////////////////
//...Convert to Elements from Anker Span;
void AnkerConstructionGlobal(int i1, int i2, int id_center = OK_STATE, CMap * ext_mp = NULL);
void AnkerSectConstructionGlobal(CMap * mp, int i1, int i2);
void AnkerProfileConstructionGlobal(CMap * mp, int i1, int i2);

////////////////////////////////////////
//...раздельное моделирование элементов;
void TowerToElementsCalculation(double * par, int & NumberOfAddedElements, int & NumberOfAddedNodes);
void TowerToElements(double * par, CMap * mp, int & PrevNumberOfElements, int & PrevNumberOfNodes, int _ExtPropNum);
void HungingPoints(double * par, double * PWires, int & NWires_max, double * PGrozo, int & NGrozo_max, CMap * mp);
double CrossHeight(double * par, int id_cross);
int    WiresCross (double * par, int id_wire);

void WireToElementsCalculation(int N, int & NumberOfAddedElements, int & NumberOfAddedNodes);
void WireToElements(int N, double * par, int & PrevNumberOfElements, int & PrevNumberOfNodes, int _ExtPropNum);
void InsulatorsToElementsCalculation(int N, int & NumberOfAddedElements, int & NumberOfAddedNodes);
void InsulatorsToElements(int N, double * par, int & PrevNumberOfElements, int & PrevNumberOfNodes, int _ExtPropNum);
void SectToElementsCalculation(int & NumberOfAddedElements, int & NumberOfAddedNodes);
void SectToElements(CMap * mp, double * par, int & PrevNumberOfElements, int & PrevNumberOfNodes, int _ExtPropNum);
void ProfileToElementsCalculation(int N, int & NumberOfAddedElements, int & NumberOfAddedNodes);
void ProfileToElements(int N, CMap * mp, double * par, int Tower1Node, int Tower2Node, int & PrevNumberOfElements, int & PrevNumberOfNodes, int _ExtPropNum);

/////////////////////////////////////
//...моделирование анкерного участка;
int  AnkerSubdivision (int *& anker_num);
void AnkerMuftaSubdivision(int i1, int i2, int *& mufta_num);

void   SetVertMontageSag(CMap * mp, double p_vert, double p, double H, double L1, double L2, int _ExtPropNum);
void	 SetHorzMontageSag(CMap * mp, double p_horz, double p, double H, double L1, double L2, int _ExtPropNum, int sg = -1);
void   ZeroMontageSag     (CMap * mp, double L1, double L2, int _ExtPropNum);
double GetMontageDiversity(CMap * mp, double L1, double L2, int _ExtPropNum);
double GetMontageDistance (CMap * mp, double L1, double L2, int _ExtPropNum1, int _ExtPropNum2);
double GetMontageDistance (CMap * mp, double L1, double L2, int _ExtPropNum, double * pp);
void   GetMontageSize(double & Sv, double & Sh, CMap * mp, double L1, double L2, int _ExtPropNum1, int _ExtPropNum2, double sg1 = 1., double sg2 = 1.);
double GetHorzMontageSize (CMap * mp, double L1, double L2, int _ExtPropNum1, int _ExtPropNum2, double sg = 1);
double GetVertMontageSize (CMap * mp, double L1, double L2, int _ExtPropNum1, int _ExtPropNum2, double sg = 1);
double GetVertMontageSize (CMap * mp, double L1, double L2, int _ExtPropNum, double del_earth, double * pp);
double GetVertMontageSize (CMap * mp, double L1, double L2, int _ExtPropNum, double del_earth, double sg = 1);

#endif
