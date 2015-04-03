// ======================================================================================================
// Montage.cpp : mathematical foundation of cable montage
//
#include "stdafx.h"
#include "montage.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifdef ___WINDOWS_LOG_MESSAGE___
#define  Message(Msg)              theView11 ? theView11->Message(Msg) : printf("%s\n", Msg)
#else
#define  Message(Msg)              printf("%s\n", Msg)
#endif

#define DECREASE_STEP_ROUGH  0.98 
#define INCREASE_STEP_ROUGH  1.02 
#define DECREASE_STEP  0.998 
#define INCREASE_STEP  1.002 

// ======================================================================================================
// Расчет стрел провеса и монтажных таблиц (Convert from Anker Span) ====================================
// ======================================================================================================
////////////////////////////////////////////////////////////////
//...уравнение перехода из одного монтажного состояния в другое;
void transit(double H0, double p0, double t0, double EF0, double & H, double p, double t, double EF, double alpha, double L, double creep, double creep0)
{
/////////////////////////////////////
//...определение конечного состояния;
	int i, M = 1;
	double p_delta = (p-p0)/M, p_curr = p0, eps = 1e-12,
			 t_delta = (t-t0)/M, t_curr = t0,
			 EF_delta = (EF-EF0)/M, EF_curr = EF0, A, Y, ggg;

	for ( i = 0; i < M; i++) {
			p_curr += p_delta;
			EF_curr += EF_delta;

			A = H0/EF0+creep0-sqr(p0*L/(2.*H0))*(1./6)-alpha*t_delta;
			Y = H0;

			while ((ggg = Y/EF_curr+creep-sqr(p_curr*L/(2.*Y))*(1./6)-A) < 0.) 
				Y += H0;

			H0 = 0.;
			H = Y;
			do {
				Y = (H+H0)*.5;
				if ((ggg = Y/EF_curr+creep-sqr(p_curr*L/(2.*Y))*(1./6)-A) < 0.) H0 = Y; else H = Y;
			}
			while (fabs(H-H0) > H*eps);

			H0 = (H+H0)*.5;
			p0 = p_curr;
			EF0 = EF_curr;
	}
	H = H0;
}

////////////////////////////////////
//...уравнение перехода для тяжения;
void transit_p(double H0, double p0, double t0, double EF0, double H, double & p, double t, double EF, double alpha, double L, double creep, double creep0)
{
/////////////////////////////////////
//...определение конечного состояния;
	int i, M = 1;
	double H_delta = (H-H0)/M, H_curr = H0,
			 t_delta = (t-t0)/M, t_curr = t0,
			 EF_delta = (EF-EF0)/M, EF_curr = EF0, A;

	for ( i = 0; i < M; i++) {
			H_curr += H_delta;
			EF_curr += EF_delta;

			A = H0/EF0+creep0-sqr(p0*L/(2.*H0))*(1./6)-alpha*t_delta;
			p = 2.*H_curr/L*sqrt(6.*fabs(H_curr/EF_curr+creep-A));

			p0 = p;
			H0 = H_curr;
			EF0 = EF_curr;
	}
}

//////////////////////////////////////////
//...уравнение перехода в аварийный режим;
void transit_d(double H0, double p0, double EF0, double & H, double L, double delta_L)
{
/////////////////////////////////////
//...определение конечного состояния;
	int i, M = 1;
	double p_curr = p0, eps = 1e-12,
			 EF_curr = EF0, l_delta = delta_L/L, A, Y, ggg;

	for ( i = 0; i < M; i++) {
			A = H0/EF0-sqr(p0*L/(2.*H0))*(1./6)-l_delta;
			Y = H0;

			while ((ggg = Y/EF_curr-sqr(p_curr*L/(2.*Y))*(1./6)-A) < 0.) 
				Y += H0;

			H0 = 0.;
			H = Y;
			do {
				Y = (H+H0)*.5;
				if ((ggg = Y/EF_curr-sqr(p_curr*L/(2.*Y))*(1./6)-A) < 0.) H0 = Y; else H = Y;
			}
			while (fabs(H-H0) > H*eps);

			H0 = (H+H0)*.5;
	}
	H = H0;
}

////////////////////////////////////////////////////////////////////////////////
//...гиперболическое уравнение перехода из одного монтажного состояния в другое;
void transit_H(double H0, double p0, double t0, double EF0, double & H, double p, double t, double EF, double alpha, double L, double delta_h)
{
/////////////////////////////////////
//...определение конечного состояния;
	int i, M = 1;
	double p_delta = (p-p0)/M, p_curr = p0, eps = 1e-12,
			 t_delta = (t-t0)/M, t_curr = t0,
			 EF_delta = (EF-EF0)/M, EF_curr = EF0, 
			 del_h = delta_h/L, kap, x_f, LL, A, C, D, Y, ggg;

	for ( i = 0; i < M; i++) {
			p_curr += p_delta;
			EF_curr += EF_delta;

			C = ((C = exp(kap = p0*L/(2.*H0)))-1./C)*.5/kap;
			x_f = del_h/C;
			x_f = log(x_f+sqrt(sqr(x_f)+1.))/kap;

			D = LL = ((D = exp(kap*x_f))+1./D)*.5*C;
			A = (1.-H0/EF0*D+alpha*t_delta)*D;
			Y = H0;

			while ((ggg = (1.-Y/EF_curr*D)*LL-A) > 0.) {
				Y += H0;
				C = ((C = exp(kap = p_curr*L/(2.*Y)))-1./C)*.5/kap;

				x_f = del_h/C;
				x_f = log(x_f+sqrt(sqr(x_f)+1.))/kap;
				LL = ((LL = exp(kap*x_f))+1./LL)*.5*C;
			}

			H0 = 0.;
			H = Y;
			do {
				Y = (H+H0)*.5;
				C = ((C = exp(kap = p_curr*L/(2.*Y)))-1./C)*.5/kap;

				x_f = del_h/C;
				x_f = log(x_f+sqrt(sqr(x_f)+1.))/kap;
				LL = ((LL = exp(kap*x_f))+1./LL)*.5*C;

  				if ((ggg = (1.-Y/EF_curr*D)*LL-A) > 0.) H0 = Y; else H = Y;
			}
			while (fabs(H-H0) > H*eps);

			H0 = (H+H0)*.5;
			p0 = p_curr;
			EF0 = EF_curr;
	}
	H = H0;
}

////////////////////////////////////////////////////////////
//...определение стрел провеса (математической и монтажной);
void strela(double p, double L, double DA, double H, double & f)
{
	double LA = L*.5+(DA*H)/(L*p);
	f = p*LA*LA/H*.5;
}
void strela0(double p, double L, double DA, double H, double & f)
{
	double LA = L*.5+(DA*H)/(L*p);
	f = (p/H*L*(LA-L*.25)-DA)*.5;
}

/////////////////////////////////////////////////////
//...определение тяжения и напряжения в нижней точке;
void tension(double p, double L, double DA, double & H, double & T, int m)
{//...m = 0; -- напряжение H по тяжению T, иначе тяжение T по напряжению H;
	if (m) T = H*(1.+sqr(DA/L)*.5)+p*p*L*L/H*.125+DA*p*.5;
	else   H = (T-DA*p*.5+sqrt(fabs(T*(T-DA*p)-p*p*L*L*.5)))/(2.+sqr(DA/L)); 
}

////////////////////////////////////////////////////////////////////////////////////////////
//...определение редуцированного тяжения аварийного режима (по уравнению нерастяжимой нити);
double reduce_tension(double p, double L, double H, double delta_L)
{
	double kappa  = 24./(L*sqr(L))*sqr(H/p),
			 reduce = H/sqrt(1.+kappa*delta_L);
	return(reduce);
}

////////////////////////////////////////////////////////////////////////////
//...вычисление приращения пролета в зависимости от редуцированного тяжения;
double reduce_I(double H0, double p0, double EF0, double H, double L)
{
	double delta_L = L*((H0-H)/EF0-(sqr(p0*L/(2.*H0))-sqr(p0*L/(2.*H)))*(1./6));
	return delta_L;
}

////////////////////////////////////////////////////////////////////////////
//...вычисление редуцированного тяжения в зависимости от приращения пролета;
double reduce_H(double p0, double delta_L, double IL, double IW, double L)
{
	double H = delta_L*(p0*L+IW)*.5/sqrt(fabs(sqr(IL)-sqr(delta_L)));
	return H;
}

///////////////////////////////////////////////////////////////////////////
//...вычисление горизонтального тяжения по максимальному тяжению в проводе;
double H_tension(double T, double p, double L, double delta_h, int hyperb)
{
	if (! L) return T;

	double A = 2.+sqr(delta_h/L),	B = T+fabs(delta_h)*p*.5, C = sqr(B)-sqr(p*L*.5)*A, 
			eps = 1e-8, H1, H2, H, Y, ff;

	if ( C < eps) return T;
	if ((H = (B+sqrt(C))/A) < p*L+eps) return H;

	if (hyperb) {
		H1 = H2 = H;
		Y  = exp(p*L/H*.5);
		Y -= 1./Y;
		Y  = fabs(delta_h)/Y;
		Y *= p/H;
		Y  = 1.-log(Y+sqrt(sqr(Y)+1.));
		Y  = exp(p*L/H*.5*Y);
		Y += 1./Y;
		ff = (Y *= H*.5);
		while (Y > T) {
			H1 *= .5;
			Y  = exp(p*L/H1*.5);
			Y -= 1./Y;
			Y  = fabs(delta_h)/Y;
			Y *= p/H1;
			Y  = 1.-log(Y+sqrt(sqr(Y)+1.));
			Y  = exp(p*L/H1*.5*Y);
			Y += 1./Y;
			Y *= H1*.5;
		}
		Y = ff;
		while (Y < T) {
			H2 *= 2.;
			Y  = exp(p*L/H2*.5);
			Y -= 1./Y;
			Y  = fabs(delta_h)/Y;
			Y *= p/H2;
			Y  = 1.-log(Y+sqrt(sqr(Y)+1.));
			Y  = exp(p*L/H2*.5*Y);
			Y += 1./Y;
			Y *= H2*.5;
		}
		do {
			H  = (H1+H2)*.5;
			Y  = exp(p*L/H*.5);
			Y -= 1./Y;
			Y  = fabs(delta_h)/Y;
			Y *= p/H;
			Y  = 1.-log(Y+sqrt(sqr(Y)+1.));
			Y  = exp(p*L/H*.5*Y);
			Y += 1./Y;
			Y *= H*.5;
			if (Y < T) H1 = H; else H2 = H;
		}
		while (fabs(H2-H1) > H*eps);
		H = (H1+H2)*.5;
	}
	return(H);
}

///////////////////////////////////////////////////////////////////////
// ...расчет критического пролета для классического уравнения перехода;
double L_transit(double H0, double p0, double t0, double EF0, double H, double p, double t, double alpha)
{
	double A = (sqr(p/(2.*H))-sqr(p0/(2.*H0)))*(1./6), L = A ? ((H-H0)/EF0+alpha*(t-t0))/A : 0.;
	return(L > 0. ? sqrt(L) : 0.);
}

///////////////////////////////////////////////////////////////////////
// ...расчет приведенного пролета для классического уравнения перехода;
double L_average(int i1, int i2, int id_set_table)
{
////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double L, del_h, L_av, f0, f1;
	int    i;

#ifndef ___OLD_LENGTH_OF_AVERAGE_SPAN___
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ...вычисляем длину приведенного пролета и осредненный модуль упругости по точной асимптотической формуле; 
	double f2, f_alpha;
	for (f0 = f1 = f2 = 0., i = i1; i <= i2; i++)
	{
		L  = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, (MAXCABLES+MAXGROZOS)*NUMSTENPOS);//...первый провод;

		f0 += L*sqr(L);
		f1 += L*(1.+sqr(del_h/L)*.5);
		f2 += L*sqr(1.+sqr(del_h/L)*.5);
	}
	L_av = f1 != 0. ? sqrt(f0/f1) : 0.;
	f_alpha = f2 != 0. ? f1/f2 : 0.;
#else
/////////////////////
// ...старая формула; 
	for (f0 = f1 = 0., i = i1; i <= i2; i++)
	{
		L  = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, (MAXCABLES+MAXGROZOS)*NUMSTENPOS);//...первый провод;

		f0 += sqr(L*L)/sqrt(sqr(L)+sqr(del_h));
		f1 += sqrt(sqr(L)+sqr(del_h));
	}
	L_av = f1 != 0. ? sqrt(f0/f1) : 0.;
#endif
////////////////////////////////////////////////////////////////////
// ...запоминаем полученную информацию о длине приведенного пролета;
	for (i = i1; id_set_table && i <= i2; i++)
		SetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-2, L_av);

	return L_av;
}

//////////////////////////////////////////////////////////////////////////////////////
// ...расчет высоты приведенного центра тяжести пролета для задания ветровой нагрузки;
double H_average(Table * tab_mod, Table * tab_clm, int id_modetab, int i1, int i2, int mode, double * Te, double * leng, double * hh, double * ff)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table) return(MAX_HIT);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double par[NUM_OF_TOWER_BASE_ELEMENTS], h_pr = 0., del_earth = 0., Leng = 0.,
			 p1 = GetTableParam(tab_mod, NUMMODE_P1),
			 p3 = GetTableParam(tab_mod, NUMMODE_P1+1), H0, Pf[3], h0, p0;
	if (mode == P3_MODE) p0 = p3; else p0 = p1;
	if (! Te) H0 = MaxTension(tab_mod, tab_clm, i1, i2, mode, 1); else H0 = Te[0];

///////////////////////////////////////////////////////////
// ...вычисляем высоту приведенного центра тяжести пролета;
	SetTableIndex(theApp.anker_list[NUMANKERHANGINGTAB], 0);
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);

	for (int ll = 0, i = i1; i <= i2; i++) {
		double h_sr = MAX_HIT, f = MAX_HIT, L = GetTableParam(table, i+1, NUMCIRCUIT-5), 
				 h1, jump;
		del_earth += GetTableParam(table, i+1, NUMCIRCUIT-3);
		Leng += L;

//////////////////////////////////////////////////////////////////////
//...вычисляем среднюю высоту подвески и стрелу в максимальном режиме;
		if (i == i1) {
			h0 = MAX_HIT;
			int num_hangingtab = id_anker_num(table, i) ? NUMANKERHANGINGTAB : NUMSTANDHANGINGTAB, nn, l, j;
			set_hanging_index(0, i, NUMANKERHANGINGTAB, 1);
			set_hanging_index(0, i, NUMSTANDHANGINGTAB, 1);

			TowerBase(theApp.tower_list, (char *)(unsigned long)GetTableParam(table, i+1, 1), par, theApp.anker_list[NUMTOWERTAB]);
			jump = GetTableParam(table, i+1, NUMCIRCUIT-4);

			if (id_modetab == 3 && tab_mod)
			for (h0 = 0., nn = j = 0; j < MAXCABLES; j++) 
			if ((l = (int)GetTableParam(theApp.anker_list[num_hangingtab], 1+j*3)) >= 0) {
				h0 = CrossHeight(par, l)-GetTableParam(theApp.anker_list[num_hangingtab], 2+j*3)+jump;
				nn++;
			}
			if (id_modetab == 2 || id_modetab == 1) {
				int	 NWires, NGrozo;
				double PWires[MAXCIRCUIT*3], PGrozo[MAXGROZOS*3];
				HungingPoints(par, PWires, NWires = MAXCIRCUIT, PGrozo, NGrozo = MAXGROZOS, NULL);

				if (id_modetab == 2 && tab_mod)
				for (h0 = 0., nn = j = 0; j < NGrozo; j++) 
				if ((l = min((int)GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+j)*3+1), NGrozo)) > 0) {
					h0 += PGrozo[l*3-2]-GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+j)*3+2)+jump;
					nn++;
				}
				if (id_modetab == 1 && tab_mod)
				for (h0 = 0., nn = j = 0; j < NWires; j++) 
				if ((l = min((int)GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+MAXGROZOS+j)*3+1), NWires)) > 0) {
					h0 += PWires[l*3-2]-GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+MAXGROZOS+j)*3+2)+jump;
					nn++;
				}
			}
			if (nn) h0 /= nn; else h0 = MAX_HIT;
		}
		h1 = MAX_HIT;
		int num_hangingtab = id_anker_num(table, i+1) ? NUMANKERHANGINGTAB : NUMSTANDHANGINGTAB, nn, l, j;
		set_hanging_index(i+1, i+1, NUMANKERHANGINGTAB, 0);
		set_hanging_index(i+1, i+1, NUMSTANDHANGINGTAB, 0);

		TowerBase(theApp.tower_list, (char *)(unsigned long)GetTableParam(table, i+2, 1), par, theApp.anker_list[NUMTOWERTAB]);
		jump = GetTableParam(table, i+2, NUMCIRCUIT-4);

		if (id_modetab == 3 && tab_mod)
		for (h1 = 0., nn = j = 0; j < MAXCABLES; j++) 
		if ((l = (int)GetTableParam(theApp.anker_list[num_hangingtab], 1+j*3)) >= 0) {
			h1 = CrossHeight(par, l)-GetTableParam(theApp.anker_list[num_hangingtab], 2+j*3)+jump+del_earth;
			nn++;
		}
		if (id_modetab == 2 || id_modetab == 1) {
			int	 NWires, NGrozo;
			double PWires[MAXCIRCUIT*3], PGrozo[MAXGROZOS*3];
			HungingPoints(par, PWires, NWires = MAXCIRCUIT, PGrozo, NGrozo = MAXGROZOS, NULL);

			if (id_modetab == 2 && tab_mod)
			for (h1 = 0., nn = j = 0; j < NGrozo; j++) 
			if ((l = min((int)GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+j)*3+1), NGrozo)) > 0) {
				h1 += PGrozo[l*3-2]-GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+j)*3+2)+jump+del_earth;
				nn++;
			}
			if (id_modetab == 1 && tab_mod)
			for (h1 = 0., nn = j = 0; j < NWires; j++) 
			if ((l = min((int)GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+MAXGROZOS+j)*3+1), NWires)) > 0) {
			  h1 += PWires[l*3-2]-GetTableParam(theApp.anker_list[num_hangingtab], (MAXCABLES+MAXGROZOS+j)*3+2)+jump+del_earth;
			  nn++;
			}
		}
		if (nn) h1 /= nn; else h1 = MAX_HIT;
		if (h0 != MAX_HIT && h1 != MAX_HIT) h_sr = (h0+h1)*.5;
		if (p0 && H0) { 
			bottom (L, 0., p0, H0, Pf); 
			f = strela0(L, 0., p0, H0, Pf);
		}
		set_hanging_index(i+1, i+1, NUMANKERHANGINGTAB, 1);
		set_hanging_index(i+1, i+1, NUMSTANDHANGINGTAB, 1);
		h0 = h1;

		if (h_sr != MAX_HIT && f != MAX_HIT) {
			h_pr += (h_sr-(2./3.)*f)*L;
			ll++;
		}
		if (ff) ff[i-i1] = f;
		if (hh) hh[i-i1] = h_sr;
		if (leng) leng[i-i1] = L;
	}
	if (Leng != 0. && ll) h_pr /= Leng; else h_pr = MAX_HIT;
	return h_pr;
}

double H_average(int num_modetab, int i1, int i2, int mode, double * Te, double * leng, double * hh, double * ff)
{
	int id_modetab = (num_modetab == NUMCABLEMODETAB ? 3 :
						  (num_modetab == NUMGROZOMODETAB ? 2 :
						  (num_modetab == NUMWIRE1MODETAB ? 1 : 0)));
	return H_average((Table *)GetTable(theApp.anker_list[num_modetab]), 
						  (Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), id_modetab, i1, i2, mode, Te, leng, hh, ff);
}

//////////////////////////////////////////////////////////////////////////////
// ...расчет предельной стрелы провеса для грозотроса по условиям грозозащиты;
double grozo_strela(int num_modetab, int i, int i_grozo, int * i_wire, double del_izol, double fi0)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table) return(0.);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double par[NUM_OF_TOWER_BASE_ELEMENTS], f_grozo = 0., f_wire,
			 p_wire  = GetTableParam(theApp.anker_list[num_modetab], NUMMODE_P1),
			 H_wire  = MaxTension(num_modetab, i, i, GROZO_MODE), Pf[3];

///////////////////////////////////////////////////////////
// ...вычисляем высоту приведенного центра тяжести пролета;
	if (0 <= i && i < table->N_group-1) {
		double f = 0., L = GetTableParam(table, i+1, NUMCIRCUIT-5);

//////////////////////////////////////////////////////////////
//...вычисляем точку подвески грозотроса относительно провода;
		int num_hangingtab = id_anker_num(table, i) ? NUMANKERHANGINGTAB : NUMSTANDHANGINGTAB, l, j;

		SetTableIndex(theApp.anker_list[num_hangingtab], 0);
		set_hanging_index(0, i, num_hangingtab, 1);

		TowerBase(theApp.tower_list, (char *)(unsigned long)GetTableParam(table, i+1, 1), par, theApp.anker_list[NUMTOWERTAB]);

		int	 NWires, NGrozo, first = 0;
		double PWires[MAXCIRCUIT*3], PGrozo[MAXGROZOS*3], fff;
		HungingPoints(par, PWires, NWires = MAXCIRCUIT, PGrozo, NGrozo = MAXGROZOS, NULL);
		i_grozo = min(i_grozo, NGrozo);

		if (0 < i_grozo && i_grozo <= NGrozo && p_wire && H_wire) {
			bottom  (L, 0., p_wire, H_wire, Pf); 
			f_wire = strela0(L, 0., p_wire, H_wire, Pf);

			if  (i_wire) 
			for (j = 0; j < i_wire[0]; j++) if ((l = i_wire[j+1]) <= NWires) {
				fff = f_wire+PGrozo[i_grozo*3-2]-PWires[l*3-2]-fabs(PGrozo[i_grozo*3-1]-PWires[l*3-1])*cos(fi0)/sin(fi0)+del_izol;
				if (! first) {
					f_grozo = fff;
					first = 1;
				}
				else f = min(f, fff);
			}
		}
	}
	return f_grozo;
}

//////////////////////////////////////////////////////////////////////
// ...линейная интерполяция коэффициента табличной величины (для ПУЭ);
double interpolation4(const double Kp_table[][4], int M, double parm, int m_index, int round)
{
	double Kp = MAX_HIT; 
	if (parm != MAX_HIT && m_index < 3) { 
      if (parm < Kp_table[0][0]) Kp = Kp_table[0][m_index+1]; else
      if (parm > Kp_table[M-1][0]) Kp = Kp_table[M-1][m_index+1]; else
      {
         int i = 0; 
         while (i < M-1 && parm >= Kp_table[i+1][0]) i++;
         if (i < M-1) {
            Kp = Kp_table[i][m_index+1]+(Kp_table[i+1][m_index+1]-Kp_table[i][m_index+1])/
                  (Kp_table[i+1][0]-Kp_table[i][0])*(parm-Kp_table[i][0]);
         }
         else
            Kp = Kp_table[M-1][m_index+1];
      }
      if (round) Kp = 0.01*((int)(100.*Kp+.5));
	}
	return Kp;
}
double interpolatchk2(const double Kp_table[][2], int M, double parm, int check, int round)
{
	double Kp = MAX_HIT; 
	if (parm != MAX_HIT && check) Kp = 1.;
	if (parm != MAX_HIT && ! check) { 
      if (parm < Kp_table[0][0]) Kp = Kp_table[0][1]; else
      if (parm > Kp_table[M-1][0]) Kp = Kp_table[M-1][1]; else
      {
         int i = 0; 
         while (i < M-1 && parm >= Kp_table[i+1][0]) i++;
         if (i < M-1) {
            Kp = Kp_table[i][1]+(Kp_table[i+1][1]-Kp_table[i][1])/
                  (Kp_table[i+1][0]-Kp_table[i][0])*(parm-Kp_table[i][0]);
         }
         else
            Kp = Kp_table[M-1][1];
      }
      if (round) Kp = 0.01*((int)(100.*Kp+.5));
	}
	return Kp;
}
double interpolation2(const double Kp_table[][2], int M, double parm, int round)
{
	double Kp = MAX_HIT; 
	if (parm != 0.) { 
      if (parm < Kp_table[0][0]) Kp = Kp_table[0][1]; else
      if (parm > Kp_table[M-1][0]) Kp = Kp_table[M-1][1]; else
      {
         int i = 0; 
         while (i < M-1 && parm >= Kp_table[i+1][0]) i++;
         if (i < M-1) {
            Kp = Kp_table[i][1]+(Kp_table[i+1][1]-Kp_table[i][1])/
                  (Kp_table[i+1][0]-Kp_table[i][0])*(parm-Kp_table[i][0]);
         }
         else
            Kp = Kp_table[M-1][1];
      }
		if (round) Kp = 0.01*((int)(100.*Kp+.5));
	}
	return Kp;
}

// ======================================================================================================
// Уточненная математика для расчета стрел провеса (уравнение тяжелой нити) =============================
// ======================================================================================================
/////////////////////////////////////////////////////////////////////////
//...определение положения нижней точки провода относительно левой опоры;
void bottom(double L, double del_h, double p, double H, double * Pf)
{
	Pf[0]  = exp(p*L/H*.5);
	Pf[0] -= 1./Pf[0];
	Pf[0]  = del_h/Pf[0];
	Pf[0] *= p/H;
	Pf[0]  = log(Pf[0]+sqrt(sqr(Pf[0])+1.));
	Pf[0]  = L*.5+H/p*Pf[0];
	Pf[1]  = exp(p*Pf[0]/H*.5); 
	Pf[1] -= 1./Pf[1];
	Pf[1] *= Pf[1];
	Pf[1]  = -H/p*.5*Pf[1];
}

////////////////////////////////////////////////////////
//...уравнение прогиба провода относительно левой опоры;
double sag(double X, double p, double H, double * Pf)
{
	double Y = exp(p*(X-Pf[0])/H*.5); 
	Y -= 1./Y;
	Y *= Y;
	Y  = H/p*.5*Y+Pf[1];
	return(Y);
}

//////////////////////////////////////////
//...определение монтажной стрелы провеса;
double strela0(double L, double del_h, double p, double H, double * Pf)
{
	double d = -del_h/L, X = log(d+sqrt(sqr(d)+1.));
	X = Pf[0]+H/p*X;
	return (d*X-sag(X, p, H, Pf));
}

///////////////////////////////////
//...определение тяжения в проводе;
double tension(double X, double p, double H, double * Pf)
{
	double T = exp(p/H*(X-Pf[0]));
	T += 1./T;
	T *= H*.5;
	return(T);
}

/////////////////////////////////////////////////////////////////////////////
//...расчет горизонтального смещения точек подвески за счет разности тяжения;
double DeltaLength(int i, int id_wire, double p_prev, double p_next)
{
	double del_L = 0.;
	if (i > 0) {
		double L = GetTableParam(theApp.anker_list[NUMANKERTAB], i, NUMCIRCUIT-5),
				 del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i,id_wire*NUMSTENPOS),
				 H1 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i, 1+id_wire*NUMSTENPOS), 
				 IL = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i, 3+id_wire*NUMSTENPOS), 
				 IW = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i, 4+id_wire*NUMSTENPOS), 
				 H2, T1, T2, V1, V2, Pf[3];

		bottom(L, del_h, p_prev, H1, Pf);
		T1 = tension (L, p_prev, H1, Pf);
		V1 = sqrt(fabs(sqr(T1)-sqr(H1)));

		L  = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1,id_wire*NUMSTENPOS);
		H2 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS);

		bottom (L, del_h, p_next, H2, Pf);
		T2 = tension (0., p_next, H2, Pf);
		V2 = sqrt(fabs(sqr(T2)-sqr(H2)));

		del_L = (H2-H1)*IL/sqrt(sqr(V1+V2+IW*.5)+sqr(H2-H1));
	}
	return del_L;
}

///////////////////////////////////////////////////////////////////////
//...расчет горизонтального смещения точек подвески в аварийном режиме;
double DamageLength(int i, int id_wire, double p, double H)
{
	double L		 = GetTableParam(theApp.anker_list[NUMANKERTAB],   i+1, NUMCIRCUIT-5),
			 IL	 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS), 
			 IW	 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS),
			 del_L = IL*H/sqrt(sqr(H)+sqr((p*L+IW)*.5));
	return del_L;
}

///////////////////////////////////////////////////////////////////////////////
//...расчет горизонтальных смещений точек подвески провода (для кабеля -- нет);
void DeltaLength(int i1, int i2, int id_wire, double p_lambda, double * p0, double * p)
{
	if (id_wire > 1)
	for (int i = i1+1; i < i2; i++)
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 2+id_wire*NUMSTENPOS,
				DeltaLength(i, id_wire, 
					p0[i-i1-1]+(p[i-i1-1]-p0[i-i1-1])*p_lambda, 
					p0[i-i1]+(p[i-i1]-p0[i-i1])*p_lambda)
		);
}

///////////////////////////////////////////////////////////////////////////////
//...расчет горизонтальных смещений точек подвески провода (для кабеля -- нет);
void DamageLength(int i1, int i2, int id_wire, double p, double * H, int m = OK_STATE)
{
	if  (m == OK_STATE)
	for (int i = i1; i <= i2; i++)
	SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 2+id_wire*NUMSTENPOS,
			DamageLength(i, id_wire, p, H[i-i1]-H[i-i1-1])
	);
	else
	for (int i = i2; i > i1; i--)
	SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 2+id_wire*NUMSTENPOS,
			DamageLength(i, id_wire, p, H[i-i1]-H[i-i1-1])
	);
}

//////////////////////////////////////////
//...сохранение тяжений анкерного участка;
void save_anker(int i1, int i2, int id_wire, double * H, int m)
{
	if  (m == OK_STATE)
	for (int i = i1; i <= i2; i++) 
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H[i-i1]);
	else
	for (int i = i1; i <= i2; i++) 
		H[i-i1] = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS);
}

///////////////////////////
//...установление нагрузки;
void set_loading(int i1, int i2, double * p, double value)
{
	for (int i = i1; i < i2; i++) 
		p[i-i1] = value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ...разрешение одного шага по тяжению анкерного участка (один шаг нижнего уровня, гиперболический обезразмеренный случай);
void anker_low_level(int i1, int i2, int id_wire, double & step, double & p_lambda, 
							double * H, double rk, double * kk, double * p0, double * p, 
							double & EF0, double & del_EF, double & del_t)
{
	DeltaLength(i1, i2, id_wire, p_lambda, p0, p);
	double EF_cur = 1./EF0+del_EF*p_lambda;
	for (int i = i1; i < i2; i++)
	{
		double  L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5), del_h, 
				kap_cur, del_kap, f0, f1, f2, f3, f4, f5, A, B, C, D, x_f, 
				LL, TT, dp_LL, dx_LL, dp_TT, dx_TT, HH;
		
		L    += GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+2, 2+(id_wire+MAXCABLES+MAXGROZOS)*NUMSTENPOS)-
				  GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 2+(id_wire+MAXCABLES+MAXGROZOS)*NUMSTENPOS);
		del_h = GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, id_wire*NUMSTENPOS)/L;
		HH    = GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS);

		del_kap = (p[i-i1]-p0[i-i1])*L*(f0 = 1./HH);
		kap_cur = (p0[i-i1]+(p[i-i1]-p0[i-i1])*p_lambda)*L*f0;
		
		x_f  = exp(kap_cur*.5);
		x_f -= 1./x_f;
		x_f  = del_h/x_f;
		x_f *= kap_cur;
		x_f  = log(x_f+sqrt(sqr(x_f)+1.));
		x_f  = .5+x_f/kap_cur;

		f2  = f4 = exp(kap_cur*.5); //...1/2 coefficient for sh, ch !!!
		f2 -= 1./f2;
		f4 += 1./f4;

		f5  = exp(kap_cur*(.5-x_f));
		f5 += 1./f5;

		D = 0.5+2.*del_h*(1.-kap_cur*.5*f4/f2)/(f2*f5);

		f1  = f3 = exp(kap_cur*x_f);
		f1 -= 1./f1;
		f3 += 1./f3;

		f2  = f4 = exp(kap_cur*(1.-x_f));
		f2 -= 1./f2;
		f4 += 1./f4;

		LL = (f1+f2)*(f0 = 1./kap_cur)*.5;
		dp_LL = ((x_f*f3+(1.-x_f)*f4)*.5-LL)*f0;
		dx_LL = (f3-f4)*f0*.5;

		f1  = f3 = exp(2.*kap_cur*x_f);
		f1 -= 1./f1;
		f3 += 1./f3;

		f2  = f4 = exp(2.*kap_cur*(1.-x_f));
		f2 -= 1./f2;
		f4 += 1./f4;

		TT = (f1+f2)*(f0 = .5/kap_cur)*.5;
		dp_TT = 2.*((x_f*f3+(1.-x_f)*f4)*.5-TT)*f0;
		dx_TT = 2.*(f3-f4)*f0*.5;

		C = (dp_LL+(D-x_f)*dx_LL)-EF_cur*HH*.5*(dp_TT+(D-x_f)*dx_TT);

		B = EF_cur*HH*(TT+0.5)+kap_cur*C;
		A = del_kap*C*HH-(del_t*LL+del_EF*HH*(TT+0.5))*HH;

		H[i-i1] += rk*(kk[i-i1] = step*A/B);
	}
}

//////////////////////////////////////////////////////////////////////////
// Явный метод Рунге-Кутты 4-го порядка с переменным шагом интегрирования 
// и экстраполяцией по Ричардсону (один шаг нижнего уровня)
void RunKut4e_Richardson_low_level(int i1, int i2, int id_wire, double & step, double & p_lambda, double & eps, 
											  double * H, double * H0, double * H2, double * kk, double * p0, double * p, 
											  double & EF0, double & del_EF, double & del_t, 
                                   int Richardson)
{
	int i, istep;
	double step2, HH, err, norm, f,
			 pp_lambda, p0_lambda, p1_lambda, p2_lambda;

// --------------------------------------------------------------------------------------
// 1.1. Запоминание вектора решения предыдущего шага
	save_anker(i1, i2, id_wire, H0);
	p0_lambda = p_lambda;
	 
// 1.2. Вычисление вектора k1 и накопление составляющих решения в H и p_lambda
	memset (H, 0, (i2-i1)*sizeof(double));
	anker_low_level(i1, i2, id_wire, step, p_lambda, H, 1., kk, p0, p, EF0, del_EF, del_t);
	p1_lambda = pp_lambda = step;

// --------------------------------------------------------------------------------------
// 2.1. Вычисление второго аргумента вектора k2
	for (p_lambda = p0_lambda+0.5*pp_lambda, i = i1; i < i2; i++) 
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+0.5*kk[i-i1]);

// 2.2. Вычисление вектора k2 и накопление составляющих решения
	anker_low_level(i1, i2, id_wire, step, p_lambda, H, 2., kk, p0, p, EF0, del_EF, del_t);
	p1_lambda += 2.*(pp_lambda = step);

// --------------------------------------------------------------------------------------
// 3.1. Вычисление второго аргумента вектора k3
	for (p_lambda = p0_lambda+0.5*pp_lambda, i = i1; i < i2; i++) 
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+0.5*kk[i-i1]);
      
// 3.2. Вычисление вектора k3 и накопление составляющих решения
	anker_low_level(i1, i2, id_wire, step, p_lambda, H, 2., kk, p0, p, EF0, del_EF, del_t);
	p1_lambda += 2.*(pp_lambda = step);

// --------------------------------------------------------------------------------------
// 4.1. Вычисление второго аргумента вектора k4
	for (p_lambda = p0_lambda+pp_lambda, i = i1; i < i2; i++) 
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+kk[i-i1]);

// 4.2. Вычисление вектора k4 и накопление составляющих решения
	anker_low_level(i1, i2, id_wire, step, p_lambda, H, 1., kk, p0, p, EF0, del_EF, del_t);
	p1_lambda += (pp_lambda = step);

// --------------------------------------------------------------------------------------
// 5.1. Вычисление вектора решения
	for (p_lambda = p0_lambda+p1_lambda*(1./6.), i = i1; i < i2; i++) 
		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+H[i-i1]*(1./6.));

// --------------------------------------------------------------------------------------
// 5.2. Повторение вычислений с половинным шагом
	if (Richardson)
	{
   // --------------------------------------------------------------------------------------
   // 5.3. Запоминаем полученное решение и восстанавливаем решение с предыдущего шага
		save_anker(i1, i2, id_wire, H2);
		p2_lambda = p_lambda;
	 
		save_anker(i1, i2, id_wire, H0, OK_STATE);
		p_lambda = p0_lambda;

		step2 = .5*step;
		istep = 2;

   // --------------------------------------------------------------------------------------
   // 5.4. Два шага цикла вычислений Рунге-Кутты с половинным шагом
		while(--istep >= 0)
		{
			save_anker(i1, i2, id_wire, H0);
			p0_lambda = p_lambda;
			 
			memset (H, 0, (i2-i1)*sizeof(double));
			anker_low_level(i1, i2, id_wire, step2, p_lambda, H, 1., kk, p0, p, EF0, del_EF, del_t);
			p1_lambda = pp_lambda = step2;

			for (p_lambda = p0_lambda+0.5*pp_lambda, i = i1; i < i2; i++) 
				SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+0.5*kk[i-i1]);

			anker_low_level(i1, i2, id_wire, step2, p_lambda, H, 2., kk, p0, p, EF0, del_EF, del_t);
			p1_lambda += 2.*(pp_lambda = step2);

			for (p_lambda = p0_lambda+0.5*pp_lambda, i = i1; i < i2; i++) 
				SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+0.5*kk[i-i1]);
      
			anker_low_level(i1, i2, id_wire, step2, p_lambda, H, 2., kk, p0, p, EF0, del_EF, del_t);
			p1_lambda += 2.*(pp_lambda = step2);

			for (p_lambda = p0_lambda+pp_lambda, i = i1; i < i2; i++) 
				SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+kk[i-i1]);

			anker_low_level(i1, i2, id_wire, step2, p_lambda, H, 1., kk, p0, p, EF0, del_EF, del_t);
			DeltaLength(i1, i2, id_wire, p_lambda, p0, p);
			p1_lambda += (pp_lambda = step2);

			for (p_lambda = p0_lambda+p1_lambda*(1./6.), i = i1; i < i2; i++) 
				SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, H0[i-i1]+H[i-i1]*(1./6.));
		}

   // --------------------------------------------------------------------------------------
   // 5.4. Уточнение решения и выбор шага интегрирования
		double fac_max = 2., fac_min = 0.5, fac = 0.9;
		for (err = 0., i = i1; i < i2; i++)
		{
			HH = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS);

			SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS, HH += (f = (HH-H2[i-i1])/15.));
			norm = max(fabs(HH), 1.); 
			err  = max(err, fabs(f)/norm);
		}
		p_lambda += (p_lambda-p2_lambda)/15.;

		if (err > 0.) f = min(fac_max, max(fac_min, fac*pow(eps/err, .2))); 
		else          f = fac_max;    

		step *= f;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...уравнение перехода из одного монтажного состояния в другое для всего анкерного участка (с учетом отклонения гирлянд);
void transit_anker(int i1, int i2, int id_wire, double * p0, double t0, double EF0, 
												double * p, double t, double EF, double alpha, double eps)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || i2 <= i1) return;

//////////////////////////////////////
//...устанавливаем начальные значения;
	double step = eps*.01, del_EF = 1./EF-1./EF0, del_t = (t-t0)*alpha, p_lambda = 0.;
	double * H  = new double[i2-i1]; 
	double * H0 = new double[i2-i1]; 
	double * H2 = new double[i2-i1];
	double * kk = new double[i2-i1];

/////////////////////////////////////////////////////
//...извлекаем параметры гирлянды в анкерном участке;
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);
	for (int i = i1+1; id_wire > 1 && i < i2; i++) {
		set_hanging_index(i == i1 ? 0 : i, i, NUMSTANDHANGINGTAB);

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 3+id_wire*3));

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4+id_wire*3));
	}

///////////////////////////////////////////////////////////////////////
//...интегрируем систему уравнений для вычисления распределения тяжений;
	while (p_lambda < 1.-step)
	RunKut4e_Richardson_low_level(i1, i2, id_wire, step, p_lambda, eps, 
											H, H0, H2, kk, p0, p, EF0, del_EF, del_t, OK_STATE);
	RunKut4e_Richardson_low_level(i1, i2, id_wire, step = 1.-p_lambda, p_lambda, eps, 
											H, H0, H2, kk, p0, p, EF0, del_EF, del_t, NULL_STATE);
   delete [] H;
   delete [] H0;
   delete [] H2;
   delete [] kk;
}

void transit_anker_test(int i1, int i2, int id_wire, double * p0, double t0, double EF0, 
												double * p, double t, double EF, double alpha, double eps)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || i2 <= i1) return;

//////////////////////////////////////
//...устанавливаем начальные значения;
	double step = eps*.01, del_EF = 1./EF-1./EF0, del_t = (t-t0)*alpha, p_lambda = 0.;
	double * H  = new double[i2-i1]; 
	double * H0 = new double[i2-i1]; 
	double * H2 = new double[i2-i1];
	double * kk = new double[i2-i1];

///////////////////////////////////////////////////////////////////
//...извлекаем параметры гирлянды в анкерном участке (для провода);
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);
	for (int i = i1+1; id_wire > 1 && i < i2; i++) {
		set_hanging_index(i == i1 ? 0 : i, i, NUMSTANDHANGINGTAB);

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 3+id_wire*3));

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4+id_wire*3));
	}

///////////////////////////////////////////////////////////////////////
//...интегрируем систему уравнений для вычисления распределения тяжений;
	FILE * TST = fopen("transit_H.dat", "w");
	double HH;
	fprintf(TST, "%g", t0+(t-t0)*p_lambda);
	for (i = i1; i < i2; i++) 
	fprintf(TST, "\t%g", HH = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS));
	fprintf(TST, "\n");

	while (p_lambda < 1.-step)
	{
	RunKut4e_Richardson_low_level(i1, i2, id_wire, step, p_lambda, eps, 
											H, H0, H2, kk, p0, p, EF0, del_EF, del_t, OK_STATE);
	fprintf(TST, "%g", t0+(t-t0)*p_lambda);
	for (i = i1; i < i2; i++) 
	fprintf(TST, "\t%g", HH = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS));
	fprintf(TST, "\n");
	}
	RunKut4e_Richardson_low_level(i1, i2, id_wire, step = 1.-p_lambda, p_lambda, eps, 
											H, H0, H2, kk, p0, p, EF0, del_EF, del_t, NULL_STATE);
	fprintf(TST, "%g", t0+(t-t0)*p_lambda);
	for (i = i1; i < i2; i++) 
	fprintf(TST, "\t%g", HH = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS));
	fprintf(TST, "\n");
	fclose (TST);
   delete [] H;
   delete [] H0;
   delete [] H2;
   delete [] kk;
}

//////////////////////////////////////////////////////////////////////////
//...уравнение перехода в авапийное состояния для всего анкерного участка;
void transit_damage_anker(int i1, int i2, int id, int id_wire, double p, double EF, double eps)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || i2 <= i1 || id < i1 || i2 < id) return;

////////////////////////////////////////////////////////////////////////////////////////
//...извлекаем параметры гирлянды в анкерном участке и устанавливаем начальные значения;
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);
	for (int i = i1+1; id_wire > 0 && i <= i2; i++) {
		set_hanging_index(i == i1+1 ? 0 : i, i, NUMSTANDHANGINGTAB);

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 3+id_wire*3));

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4+id_wire*3));
	}

///////////////////////////////////////////////////////////////////////////////////////
//...устанавливаем предельное редуцированное тяжение (максимальное отклонение гирлянд);
	double * H = new double[i2-i1+1]; 
	save_anker(i1, i2, id_wire, H); H[id-i1] = 0.;
	save_anker(i1, i2, id_wire, H, OK_STATE);

///////////////////////////////////////////////////////////
//...устрамваем итерационный процесс по отклонению гирлянд;
	for (int k = 0; k < 40; k++) {
//		DamageLength(id+1, i2, id_wire, p, H+id+1-i1);
		for (i = id+1; i < i2; i++) {
			DamageLength(i, i, id_wire, p, H+i-i1);
			transit_d(GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS), 
						 p, EF, H[i-i1], GetTableParam(table, i+1, NUMCIRCUIT-5),
						 GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 2+id_wire*NUMSTENPOS)-
						 GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+2, 2+id_wire*NUMSTENPOS));
		}
		DamageLength(i, i, id_wire, p, H+i-i1);
		transit_d(GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS), 
					 p, EF, H[i-i1], GetTableParam(table, i+1, NUMCIRCUIT-5),
					 GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 2+id_wire*NUMSTENPOS));
		save_anker(i1, i2, id_wire, H, OK_STATE);
	}
	delete [] H;
}

///////////////////////////////////
//...обратный итерационный процесс;
void transit_damage_anker_inverse(int i1, int i2, int id, int id_wire, double p, double EF, double eps = 1E-15)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || i2 <= i1 || id < i1 || i2 < id) return;

////////////////////////////////////////////////////////////////////////////////////////
//...извлекаем параметры гирлянды в анкерном участке и устанавливаем начальные значения;
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);
	for (int i = i1+1; id_wire > 0 && i <= i2; i++) {
		set_hanging_index(i == i1+1 ? 0 : i, i, NUMSTANDHANGINGTAB);

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 3+id_wire*3));

		SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS, 
		GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4+id_wire*3));
	}

///////////////////////////////////////////////////////////////////////////////////////
//...устанавливаем предельное редуцированное тяжение (максимальное отклонение гирлянд);
	double * H = new double[i2-i1+1], delta_L; 
	save_anker(i1, i2, id_wire, H); H[id-i1] = 0.;
	save_anker(i1, i2, id_wire, H, OK_STATE);

///////////////////////////////////////////////////////////
//...устрамваем итерационный процесс по отклонению гирлянд;
	for (int k = 0; k < 40; k++) {
		for (i = id+1; i < i2; i++) {
			if (k == 0)
			delta_L  = GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS)*.98;
			H[i-i1] += reduce_H(p, delta_L, 
									GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS),
									GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS), GetTableParam(table, i+1, NUMCIRCUIT-5));
			delta_L = reduce_I(GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS), 
									p, EF, H[i-i1], GetTableParam(table, i+1, NUMCIRCUIT-5));
		}
		if (k == 0)
		delta_L = GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS)*.98;
		H[i-i1] = reduce_H(p, delta_L, 
								GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 3+id_wire*NUMSTENPOS),
								GetTableParam (theApp.anker_list[NUMTENSIONTAB], i+1, 4+id_wire*NUMSTENPOS), GetTableParam(table, i+1, NUMCIRCUIT-5));
		delta_L = reduce_I(GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+id_wire*NUMSTENPOS), 
								p, EF, H[i-i1], GetTableParam(table, i+1, NUMCIRCUIT-5));
	}
	save_anker(i1, i2, id_wire, H, OK_STATE);
	delete [] H;
}

/////////////////////////////////////////////////////////////////////
// ...устанавливаем индекс параметров подвески для заданного пролета;
void set_hanging_index(int i1, int i2, int num_hangingtab, int id_second, int id_second2)
{
	Table * table = (Table *)GetTable(theApp.anker_list[num_hangingtab]);
	int k, id_first = 1;

	if (table)
	for (int i = i2; id_first && i >= i1; i--)
	for (    k = 0;  id_first && k < table->N_group; k++) {
		if (((int *)table->table[0].parm_values)[k+1] == i+1) id_first = 0;
		if (! id_first && i != i2) id_second = id_second2; 

		while (! id_first && id_second && k+1 < table->N_group &&
				((int *)table->table[0].parm_values)[k+2] == i+1) {
			id_second--;
			k++;
		}
	}
	if (! id_first) SetTableIndex(theApp.anker_list[num_hangingtab], k); //...k -- т.к.цикл сдвигает k на 1 !!!
	//...индекс сохраняет свое старое значение, если нет вхождений;
}

////////////////////////////////////////////////////////////////////
// ...устанавливаем индекс режимов нагружения для заданного пролета;
void set_mode_index(int i1, int i2, int num_modetab)
{
	Table * table = (Table *)GetTable(theApp.anker_list[num_modetab]);
	int k, id_first = 1;

	if (table)
	for (int i = i2; id_first && i >= i1; i--)
	for (    k = 0;  id_first && k < table->N_group; k++)
	if (((int *)table->table[0].parm_values)[k+1] == i+1) id_first = 0;
	if (! id_first) SetTableIndex(theApp.anker_list[num_modetab], k); //...k -- т.к.цикл сдвигает k на 1 !!!
	//...индекс сохраняет свое старое значение, если нет вхождений;
}

///////////////////////////////////////////////////////////////////
// ...устанавливаем режимы нагружения провода или кабеля в пролете;
void set_mode_loading(int num_modetab)
{
///////////////////////////////////////////////////////////////
//...вычисление нагрузок по найденным климатическим параметрам;
	for ( int k = GetNGroup(theApp.anker_list[num_modetab]); k >= 1; k--) {
		double D = GetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1-1),
				p1 = GetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1), 
				p2 = 0.,
				p3 = p1, 
				p4 = 0.,
				p5 = 0.,
				p6 = p1,
				p7 = p1,
				C_glaze      = GetTableParam(theApp.anker_list[NUMCLIMATETAB], k, 2),
				C_rush_glaze = GetTableParam(theApp.anker_list[NUMCLIMATETAB], k, 3),
				Q_rush       = GetTableParam(theApp.anker_list[NUMCLIMATETAB], k, 9),
				Q_rush_glaze = GetTableParam(theApp.anker_list[NUMCLIMATETAB], k, 10),
				Ro_glaze = GetTableParam(theApp.anker_list[NUMCLIMATETAB], k, 13), 
				Cx_rush	= D < 20 ? 1.2 : 1.1,
				Cx_glaze = 1.2;

		int id_set_ice_and_wind_mode = 1;
		if (id_set_ice_and_wind_mode) { 
			p2 = Ro_glaze*M_PI*C_glaze*(D+C_glaze)*M_G*1e-3;
			p4 = Q_rush*D*Cx_rush;
			p5 = Q_rush_glaze*(D+C_rush_glaze*2.)*Cx_glaze;
			p3 = p1+p2;
			p6 = sqrt(sqr(p1)+sqr(p4));
			p7 = sqrt(sqr(p3)+sqr(p5));
		}
		SetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1,	p1);
		SetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1+1, p3);
		SetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1+2, p6);
		SetTableParam(theApp.anker_list[num_modetab], k, NUMMODE_P1+3, p7);
	}
}

////////////////////////////////////////////////////
// ...устанавливаем все режимы нагружения через ПУЕ;
void set_pue_loading(int i_ank)
{
/////////////////////////////////////////////
//...запускаем справочник по ПУЭ-6 или ПУЭ-7;
	unsigned long cond = (unsigned long)GetTableParam(theApp.anker_list[NUMCLIMATETAB], 14);
	if ((cond & 0xF0000000) == 0x10000000) {
		CPue7Prp * property = new CPue7Prp( (Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMWIRE1MODETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMGROZOMODETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMCABLEMODETAB]), 
														(theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
														(theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else if ((cond & 0xF0000000) == 0x00000000) {
		CPue6Prp * property = new CPue6Prp( (Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMWIRE1MODETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMGROZOMODETAB]), 
														(Table *)GetTable(theApp.anker_list[NUMCABLEMODETAB]), 
														(theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
														(theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else {
		set_mode_loading(NUMWIRE1MODETAB);
		set_mode_loading(NUMGROZOMODETAB);
		set_mode_loading(NUMCABLEMODETAB);
	}
}
void set_pue_wire(int i_ank)
{
/////////////////////////////////////////////
//...запускаем справочник по ПУЭ-6 или ПУЭ-7;
	unsigned long cond = (unsigned long)GetTableParam(theApp.anker_list[NUMCLIMATETAB], 14);
	if ((cond & 0xF0000000) == 0x10000000) {
		CPue7Prp * property = new CPue7Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 (Table *)GetTable(theApp.anker_list[NUMWIRE1MODETAB]), 
																 NULL, 
																 NULL, 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else if ((cond & 0xF0000000) == 0x00000000) {
		CPue6Prp * property = new CPue6Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 (Table *)GetTable(theApp.anker_list[NUMWIRE1MODETAB]), 
																 NULL, 
																 NULL, 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else set_mode_loading(NUMWIRE1MODETAB);
}
void set_pue_grozo(int i_ank)
{
/////////////////////////////////////////////
//...запускаем справочник по ПУЭ-6 или ПУЭ-7;
	unsigned long cond = (unsigned long)GetTableParam(theApp.anker_list[NUMCLIMATETAB], 14);
	if ((cond & 0xF0000000) == 0x10000000) {
		CPue7Prp * property = new CPue7Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 NULL, 
																 (Table *)GetTable(theApp.anker_list[NUMGROZOMODETAB]), 
																 NULL, 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else if ((cond & 0xF0000000) == 0x00000000) {
		CPue6Prp * property = new CPue6Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 NULL, 
																 (Table *)GetTable(theApp.anker_list[NUMGROZOMODETAB]), 
																 NULL, 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else set_mode_loading(NUMGROZOMODETAB);
}
void set_pue_cable(int i_ank)
{
/////////////////////////////////////////////
//...запускаем справочник по ПУЭ-6 или ПУЭ-7;
	unsigned long cond = (unsigned long)GetTableParam(theApp.anker_list[NUMCLIMATETAB], 14);
	if ((cond & 0xF0000000) == 0x10000000) {
		CPue7Prp * property = new CPue7Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 NULL, 
																 NULL, 
																 (Table *)GetTable(theApp.anker_list[NUMCABLEMODETAB]), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else if ((cond & 0xF0000000) == 0x00000000) {
		CPue6Prp * property = new CPue6Prp((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), 
																 NULL, 
																 NULL, 
																 (Table *)GetTable(theApp.anker_list[NUMCABLEMODETAB]), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank] : 0), 
																 (theApp.anker_num && i_ank <= theApp.anker_num[0] ? theApp.anker_num[i_ank+1]-1 : 0));
		delete property;
	}
	else set_mode_loading(NUMCABLEMODETAB);
}

/////////////////////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения провода и кабеля; 
void get_temper_mode(Table * tab_clm, double & t_glaze, double & t_rush, double & t_aver, double & t_abs_max, double & t_abs_min, double & t_mtg, Table * tab_montage = NULL)
{
//////////////////////////////////////
//...извлечение температурных режимов;
	t_glaze   = GetTableParam(tab_clm, 4);
	t_rush    = GetTableParam(tab_clm, 5); //...температура ветрового режима;
	t_aver    = GetTableParam(tab_clm, 6);
	t_abs_max = GetTableParam(tab_clm, 7);
	t_abs_min = GetTableParam(tab_clm, 8);

//////////////////////////////////////
//...извлечение монтажной температуры;
   if (! tab_montage) t_mtg = 0.;
	else {
		theApp.Opt.montage_temper = max(theApp.Opt.montage_temper, 0);
		theApp.Opt.montage_temper = min(theApp.Opt.montage_temper, tab_montage->N-1);
		t_mtg = user_strtod(tab_montage->table[theApp.Opt.montage_temper].parm_name);
	}
}

void get_temper_mode(double & t_glaze, double & t_rush, double & t_aver, double & t_abs_max, double & t_abs_min, double & t_mtg)
{
	Table * table; 
   if ( ! (table = (Table *)GetTable(theApp.anker_list[NUMCABLEMONTAGETAB])) &&
		  ! (table = (Table *)GetTable(theApp.anker_list[NUMGROZOMONTAGETAB])) &&
		  ! (table = (Table *)GetTable(theApp.anker_list[NUMWIRE1MONTAGETAB]))) t_mtg = 0.;
	get_temper_mode((Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_mtg, table);
}

/////////////////////////////////////////////////////////////////////
// ...извлекаем все режимы нагружения провода, кабеля или грозотроса; 
void get_loading(Table * tab_mode, double & p1, double & p3, double & p6, double & p7, 
											  double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, 
											  double & EF_creep,  double & EF, double & EF_fin, double & alpha, double & creep, char *& mark)
{
/////////////////////////
//...извлечение нагрузок;
	p1 = GetTableParam(tab_mode, NUMMODE_P1);
	p3 = GetTableParam(tab_mode, NUMMODE_P1+1);
	p6 = GetTableParam(tab_mode, NUMMODE_P1+2);
	p7 = GetTableParam(tab_mode, NUMMODE_P1+3);

///////////////////////////////////
//...извлечение предельных тяжений;
	double T_break = GetTableParam(tab_mode, NUMTENSION);
	T_max		  = GetTableParam(tab_mode, NUMTENSION+1); //...максимальная горизонтальная нагрузка на опору;
	T_opt		  = GetTableParam(tab_mode, NUMTENSION+2);
	T_abs_max  = GetTableParam(tab_mode, NUMTENSION+3)*T_break*.01; //...максимальное напряжение (в терминах тяжения);
	T_abs_tmin = GetTableParam(tab_mode, NUMTENSION+5)*T_break*.01;

///////////////////////////////////
//...извлечение параметров тяжения;
	EF_creep = GetTableParam(tab_mode, NUMEFMODUL-1);
	EF			= GetTableParam(tab_mode, NUMEFMODUL);
	EF_fin	= GetTableParam(tab_mode, NUMEFMODUL+1);
	alpha		= GetTableParam(tab_mode, NUMEFMODUL+2);
//	creep		= GetTableParam(tab_mode, NUMEFMODUL+4);
	creep		= T_abs_max*(1./EF_creep-1./EF_fin); //...модель нормальной вытяжки по касательному модулю;

//////////////////////////////////////////
//...марка провода, кабеля или грозотроса; 
	mark = (char *)((unsigned long)GetTableParam(tab_mode, 1));
}
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, 
											 double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, 
											 double & EF_creep,  double & EF, double & EF_fin, double & alpha, double & creep, char *& mark)
{	get_loading((Table *)GetTable(theApp.anker_list[num_modetab]), p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, 
											EF_creep,  EF, EF_fin, alpha, creep, mark);
}
void get_loading(Table * tab_mode, double & p1, double & p3, double & p6, double & p7, double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, double & EF_creep, double & EF, double & EF_fin, double & alpha, double & creep)
{
	char * mark;
	get_loading(tab_mode, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, double & EF_creep, double & EF, double & EF_fin, double & alpha, double & creep)
{	
	get_loading((Table *)GetTable(theApp.anker_list[num_modetab]), p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep);
}
void get_loading(int num_modetab, double & p1, double & p3, double & p6, double & p7, double & EF, double & EF_fin, double & alpha, double & creep, char *& mark)
{
	double T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep;
	get_loading(num_modetab, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}
void get_loading(Table * tab_mode, double & p1, double & p3, double & p7, double & T_abs_max, double & T_opt, double & EF, double & EF_fin, double & alpha, double & creep)
{
	char * mark;
	double p6, T_max, T_abs_tmin, EF_creep;
	get_loading(tab_mode, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}
void get_loading(int num_modetab, double & p1, double & p3, double & p7, double & T_abs_max, double & T_opt, double & EF, double & EF_fin, double & alpha, double & creep)
{
	get_loading((Table *)GetTable(theApp.anker_list[num_modetab]), p1, p3, p7, T_abs_max, T_opt, EF, EF_fin, alpha, creep);
}
void get_loading(int num_modetab, double & p1, double & p7, double & T_abs_max, double & T_abs_tmin, double & T_opt, double & T_max, double & EF, double & EF_fin, double & alpha, double & creep)
{
	char * mark;
	double p3, p6, EF_creep;
	get_loading(num_modetab, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}
void get_loading(int num_modetab, double & p1, double & EF, double & alpha, char *& mark)
{
	double p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_fin, creep, EF_creep;
	get_loading(num_modetab, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}
void get_loading(int num_modetab, char *& mark)
{
	double p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep;
	get_loading(num_modetab, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep, mark);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ...устанавливаем монтажную таблицу по тяжению конкретного режима нагружения провода, кабеля или грозотроса;
void set_montage_table(int num_montagetab, int i_ank, double H, int mode, double * t_op)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;
	int i1 = theApp.anker_num[i_ank],
		 i2 = theApp.anker_num[i_ank+1]-1, m = 1;

/////////////////////////////////////////////////////////////////////////////////
// ...извлекаем параметры провода, кабеля и грозотроса (и все режимы нагружения); 
	char * mark;
	double p1, p3, p6, p7, EF, EF_fin, alpha, creep, p;
	get_loading(num_montagetab-3, p1, p3, p6, p7, EF, EF_fin, alpha, creep, mark);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, t_norm, t_grozo = 15.; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

///////////////////////////
//...нагрузочное состояние;
	double L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT-2), H1, f0;
	switch (mode) {
	case			  NULL_MODE:  p = 0.; t_norm = t_montage; break;
	case			 GROZO_MODE:  p = p1; t_norm = t_grozo; break;
	case		OPERATING_MODE:  p = p1; t_norm = t_op ? t_op[0] : t_montage; break;
	case				 P1_MODE:  p = p1; t_norm = t_aver;    break;
	case	  MIN_TEMPER_MODE:  p = p1; t_norm = t_abs_min; break;
	case	  MAX_TEMPER_MODE:  p = p1; t_norm = t_abs_max; break;
	case				 P3_MODE:  p = p3; t_norm = t_glaze;   break;
	case P6_RIGHT_RUSH_MODE:  p = p6; t_norm = t_rush; break;
	case  P6_LEFT_RUSH_MODE:  p = p6; t_norm = t_rush; break;
	case P7_RIGHT_RUSH_MODE:  p = p7; t_norm = t_glaze; break;
	case  P7_LEFT_RUSH_MODE:  p = p7; t_norm = t_glaze; break;
	case		  MONTAGE_MODE:  p = p1; t_norm = t_montage; m = 0; break;
	default:						  p = p1; t_norm = 0.; break;
	}
	if (L_av == 0.) L_av = L_average(i1, i2, 1);
	if (mode == NULL_MODE || H < ___EE___) H1 = 0.; else {
		if (m) {
			transit(H, p, t_norm, EF_fin, f0, p1, t_aver, EF_fin, alpha, L_av);
			creep_transit(num_montagetab-3, f0, p1, t_aver, H1, t_montage, L_av, 0);
		}
		else
			transit(H, p, t_norm, EF, H1, p1, t_montage, EF, alpha, L_av);
	}

//////////////////////////////////////////////////////////////////////////////////////////////////
// ...вносим полученную информацию о состоянии провода, кабеля или грозотроса в монтажную таблицу;
	SetTableParam(theApp.anker_list[num_montagetab], i_ank, theApp.Opt.montage_temper, H1);

	if (num_montagetab == NUMCABLEMONTAGETAB)
		CalculateCableMontageTable(i_ank);

	if (num_montagetab == NUMGROZOMONTAGETAB)
		CalculateGrozoMontageTable(i_ank);

	if (num_montagetab == NUMWIRE1MONTAGETAB)
		CalculateWireMontageTable(i_ank);
}

/////////////////////////////////////////////////////////////////////////////////////////
// ...устанавливаем тяжение конкретного режима нагружения провода, кабеля или грозотроса;
void set_mode_tension(int num_montagetab, int i_ank, double & pV, double & pH, int mode, double * t_op)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;
	int i1 = theApp.anker_num[i_ank],
		 i2 = theApp.anker_num[i_ank+1]-1, i, j, m = 1;

/////////////////////////////////////////////////////////////////////////////////
// ...извлекаем параметры провода, кабеля и грозотроса (и все режимы нагружения); 
	char * mark;
	double p1, p3, p6, p7, EF, EF_fin, alpha, creep, p;
	get_loading(num_montagetab-3, p1, p3, p6, p7, EF, EF_fin, alpha, creep, mark);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, t_norm, t_grozo = 15.; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT-2), 
			 H1 = GetTableParam(theApp.anker_list[num_montagetab], i_ank, theApp.Opt.montage_temper), H, f0;

///////////////////////////
//...нагрузочное состояние;
	switch (mode) {
	case			  NULL_MODE: pV = 0.; pH = 0.; p = 0.; t_norm = t_montage; break;
	case			 GROZO_MODE: pV = p1; pH = 0.; p = p1; t_norm = t_grozo; break;
	case		OPERATING_MODE: pV = p1; pH = 0.; p = p1; t_norm = t_op ? t_op[0] : t_montage; break;
	case				 P1_MODE: pV = p1; pH = 0.; p = p1; t_norm = t_aver;    break;
	case	  MIN_TEMPER_MODE: pV = p1; pH = 0.; p = p1; t_norm = t_abs_min; break;
	case	  MAX_TEMPER_MODE: pV = p1; pH = 0.; p = p1; t_norm = t_abs_max; break;
	case				 P3_MODE: pV = p3; pH = 0.; p = p3; t_norm = t_glaze;   break;
	case P6_RIGHT_RUSH_MODE: pV = p1; pH = -sqrt(sqr(p6)-sqr(p1)); p = p6; t_norm = t_rush; break;
	case  P6_LEFT_RUSH_MODE: pV = p1; pH =  sqrt(sqr(p6)-sqr(p1)); p = p6; t_norm = t_rush; break;
	case P7_RIGHT_RUSH_MODE: pV = p3; pH = -sqrt(sqr(p7)-sqr(p3)); p = p7; t_norm = t_glaze; break;
	case  P7_LEFT_RUSH_MODE: pV = p3; pH =  sqrt(sqr(p7)-sqr(p3)); p = p7; t_norm = t_glaze; break;
	case		  MONTAGE_MODE: pV = p1; pH = 0.; p = p1; H = H1; m = 0; break;
	default:						  p = p1; t_norm = 0.; break;
	}
	if (L_av == 0.) L_av = L_average(i1, i2, 1);
	if (mode == NULL_MODE || H1 < ___EE___) H = 0.; else {
		if (m) {
			creep_transit(num_montagetab-3, H1, p1, t_montage, f0, t_aver, L_av);
			transit(f0, p1, t_aver, EF_fin, H, p, t_norm, EF_fin, alpha, L_av);
		}
	}

////////////////////////////////////////////////////////////////////////////////////////////////
// ...вносим полученную информацию о состоянии провода, кабеля или грозотроса в таблицу тяжений;
	if (num_montagetab == NUMCABLEMONTAGETAB)
		for (j = 0;  j < MAXCABLES; j++)
		for (i = i1; i <= i2; i++)
			SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+j*NUMSTENPOS, H);

	if (num_montagetab == NUMGROZOMONTAGETAB)
		for (j = 0;  j < MAXGROZOS; j++)
		for (i = i1; i <= i2; i++)
			SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+j)*NUMSTENPOS, H);

	if (num_montagetab == NUMWIRE1MONTAGETAB)
		for (i = i1; i <= i2; i++) {
			int  N_wires = min((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT), MAXCIRCUIT);
			for (j = 0;  j < N_wires; j++)
				SetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+MAXGROZOS+j)*NUMSTENPOS, H);
		}
}

///////////////////////////////////////////////////////////////
// ...сообщение об установленном монтажном состоянии в пролете;
void MessageMontageState(int mode)
{
	switch (mode) {
	case			  NULL_MODE: Message("\"NULL\" -- mode"); break;
	case			 GROZO_MODE: Message("\"Grozo\" -- mode"); break;
	case		OPERATING_MODE: Message("\"Operating\" -- mode"); break;
	case				 P1_MODE: Message("\"P1\" -- mode"); break;
	case	  MIN_TEMPER_MODE: Message("\"Min temperature\" -- mode"); break;
	case	  MAX_TEMPER_MODE: Message("\"Max temperature\" -- mode"); break;
	case				 P3_MODE: Message("\"P3\" -- mode"); break;
	case P6_RIGHT_RUSH_MODE: Message("\"P6\" -- mode, right rush"); break;
	case  P6_LEFT_RUSH_MODE: Message("\"P6\" -- mode, left rush"); break;
	case P7_RIGHT_RUSH_MODE: Message("\"P7\" -- mode, right rush"); break;
	case  P7_LEFT_RUSH_MODE: Message("\"P7\" -- mode, left rush"); break;
	case		  MONTAGE_MODE: Message("\"Montage state\" -- mode"); break;
	default:						 Message("Not defined mode"); break;
	}
}

// ======================================================================================================
// Anker Montage (Convert to Elements from Anker Span and cycle of montage ) ============================
// ======================================================================================================
//////////////////////////////////////////////////////////////////////////
// ...монтаж одного анкерного участка (с проверкой максимального тяжения);
void AnkerMontage(int i_ank, double S_stand, double S_anker, double S_sect)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || ! theApp.anker_num) return;
	int i1 = theApp.anker_num[i_ank],
		 i2 = theApp.anker_num[i_ank+1]-1, i, j;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i1+1, NUMCIRCUIT+8)));

	set_mode_index(0, i1, NUMCABLEMODETAB);
	set_mode_index(0, i1, NUMGROZOMODETAB);
	set_mode_index(0, i1, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////////
// ...конструируем пролет в виде FEM модели;
	Nodes.SetNodes(0, NULL, 0);
	Elements.SetFElements(-1, & Nodes, & Properties);
	Elements.Numbering();
	Elements.NumGroups();
	Elements.Initials();

	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1.}; 
	AnkerConstructionGlobal(i1, i2, OK_STATE, (double *)(& mp));

///////////////////////////////////////////
// ...вычисляем длину критического пролета;
/*	double L_tr;
	L_tr = L_Transit(NUMWIRE1MODETAB);
	L_tr = L_Transit(NUMCABLEMODETAB);
	L_tr = L_Transit(NUMGROZOMODETAB);

///////////////////////////////////////////
// ...вычисляем длину критического пролета;
	int i_grozo = 1, i_wire[] = {1, 1};
	double f_grozo = grozo_strela(NUMWIRE1MODETAB, i1, i_grozo, i_wire, 4.5);
*/
////////////////////////////////////////////////////////////////////////
// ...определяем необходимые параметры для провода, кабеля и грозотроса; 
	double p1_wr,  p3_wr,  p7_wr,  T0_opt, T0_abs_max, EF_wr,  EF0_fin, alpha_wr,  creep_wr,
			 p1_cab, p3_cab, p7_cab, T1_opt, T1_abs_max, EF_cab, EF1_fin, alpha_cab, creep_cab,
			 p1_gro, p3_gro, p7_gro, T2_opt, T2_abs_max, EF_gro, EF2_fin, alpha_gro, creep_gro;

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double L_av, L1, L2, pp[3], H_cab, H1_cab, H3_cab, H7_cab, H_wr, H1_wr, H3_wr, H_gro, H1_gro, H3_gro,	f0, f1, 
			 S_wind = 1.5, pV_wr, pH_wr, pV_cab, pH_cab, S;
	int N_wires = min((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT), MAXCIRCUIT), 
		 mode = P7_RIGHT_RUSH_MODE, id_Sw = NULL_STATE,//...параметры боковой нагрузки;
		 N1wires = min(2, N_wires), N2wires = min(3, N_wires); 

///////////////////////////////////////////////////////////////////////////
// ...вычисляем предельное тяжение провода и монтируем провод по габаритам;
	L_av  = L_average(i1, i2);
	H1_wr = MaxTension(NUMWIRE1MODETAB, i1, i2);
	f0    = OptTension(NUMWIRE1MODETAB, i1, i2);

	get_loading(NUMWIRE1MODETAB, p1_wr,  p3_wr,  p7_wr,  T0_abs_max, T0_opt, EF_wr,  EF0_fin, alpha_wr, creep_wr);
	double Size[3], Y0 = mp[2]; //...S_stand минимальное расстояние до земли в промежуточном пролете;
	do {								 //...S_anker минимальное расстояние до земли в анкерном (начальном и конечном) пролете;
		H1_wr *= DECREASE_STEP_ROUGH;
		transit(H1_wr, p1_wr, t_aver, EF0_fin, H3_wr, p3_wr, t_glaze, EF0_fin, alpha_wr, L_av);

		if (S_stand != 0. || S_anker != 0.)
		for (mp[2] = Y0, L2 = 0., Size[0] = Size[1] = MAX_HIT, i = i1; i <= i2; i++)
		{
			L1 = L2;
			L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

			SetVertMontageSag(mp, p3_wr, p3_wr, H3_wr, L1, L2, WIRES_PROP+1);
			f1 = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3));

			mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

			if (i1 < i && i < i2) {
				if (f1 < Size[0]) Size[0] = f1; 
			}
			else {
				if (f1 < Size[1]) Size[1] = f1;
			}
		}
	}
	while ((S_stand == 0. || S_stand < Size[0]) && (S_anker == 0. || S_anker < Size[1]) && T0_opt < H1_wr);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
	creep_transit(NUMWIRE1MODETAB, H1_wr, p1_wr, t_aver, H_wr, t_montage, L_av, 0);
	SetTableParam(theApp.anker_list[NUMWIRE1MONTAGETAB], i_ank, theApp.Opt.montage_temper, H_wr);
	CalculateWireMontageTable(i_ank);

//////////////////////////////////////////////////////////////////////
// ...монтируем кабель до заданного габарита или до заданного тяжения;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB])) {

/////////////////////////////////////////////////////////////////////
// ...вычисляем предельное тяжение кабеля (и боковую стрелу провеса);
		H1_cab = MaxTension(NUMCABLEMODETAB, i1, i2);
		f0     = OptTension(NUMCABLEMODETAB, i1, i2);
		if (id_Sw == OK_STATE) {
			set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr,  pH_wr,  mode);
			set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode);
		}

////////////////////////////////////////////////////////////////
// ...делаем расчет (c учетом стрелы бокового ветра для кабеля);
		get_loading(NUMCABLEMODETAB, p1_cab, p3_cab, p7_cab, T1_abs_max, T1_opt, EF_cab, EF1_fin, alpha_cab, creep_cab);
		do {
			H1_cab *= DECREASE_STEP_ROUGH;
			transit(H1_cab, p1_cab, t_aver, EF1_fin, H7_cab, p7_cab, t_glaze, EF1_fin, alpha_cab, L_av);
			transit(H1_cab, p1_cab, t_aver, EF1_fin, H3_cab, p3_cab, t_glaze, EF1_fin, alpha_cab, L_av);
			
			if (S_stand != 0. || S_anker != 0. || S_wind != 0. && id_Sw == OK_STATE)
			for (mp[2] = Y0, L2 = 0., Size[0] = Size[1] = MAX_HIT, Size[2] = S_wind, i = i1; i <= i2; i++)
			{
				L1 = L2;
				L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

				SetVertMontageSag(mp, p3_cab, p3_cab, H3_cab, L1, L2, CABLE_PROP+1);
				f1 = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3));

				if (i1 < i && i < i2) {
					if (f1 < Size[0]) Size[0] = f1; 
				}
				else {
					if (f1 < Size[1]) Size[1] = f1;
				}

///////////////////////////////////////////////////////////////
// ...устанавливаем стрелу бокового ветра для кабеля и провода;
				if (id_Sw == OK_STATE) {
					SetHorzMontageSag(mp, pH_cab, p7_cab, H7_cab, L1, L2, CABLE_PROP+1);
					for (j = N1wires; j <= N2wires; j++)
						SetHorzMontageSag(mp, pH_wr, p7_wr, GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(j-1+MAXCABLES+MAXGROZOS)*NUMSTENPOS), L1, L2, WIRES_PROP+j);

					for (j = N1wires; j <= N2wires; j++)
					if ((S = GetHorzMontageSize(mp, L1, L2, CABLE_PROP+1, WIRES_PROP+j, -1.)) < f1 || j == 1) f1 = S; 
					if (i1 == i || Size[2] < f1) Size[2] = f1; 

					SetHorzMontageSag(mp, -pH_cab, p7_cab, H7_cab, L1, L2, CABLE_PROP+1);
					for (j = N1wires; j <= N2wires; j++)
						SetHorzMontageSag(mp, -pH_wr, p7_wr, GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(j-1+MAXCABLES+MAXGROZOS)*NUMSTENPOS), L1, L2, WIRES_PROP+j);

					for (j = N1wires; j <= N2wires; j++)
					if ((S = GetHorzMontageSize(mp, L1, L2, CABLE_PROP+1, WIRES_PROP+j, 1.)) < f1) f1 = S; 
					if (Size[2] < f1) Size[2] = f1; 
				}
				mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);
			}
		}
		while ((S_stand == 0. || S_stand < Size[0]) && (S_anker == 0. || S_anker < Size[1]) && 
				 (S_wind  == 0. || S_wind <= Size[2] || id_Sw != OK_STATE) && T1_opt < H1_cab);

/////////////////////////////////////////////////////////////////////////////////////////
// ...натягиваем кабель до заданного габарита в пересечении или до максимального тяжения;
		if (tab_sct)
		do { //...S_sect минимальное расстояние до пересечения в пролете;    
			H1_cab *= INCREASE_STEP_ROUGH;
			transit(H1_cab, p1_cab, t_aver, EF1_fin, H3_cab, p3_cab, t_glaze, EF1_fin, alpha_cab, L_av);
			transit(H1_cab, p1_cab, t_aver, EF1_fin, H7_cab, p7_cab, t_glaze, EF1_fin, alpha_cab, L_av);

			if (S_sect != 0.)
			for (mp[2] = Y0, L2 = 0., Size[0] = MAX_HIT, i = i1; i <= i2; i++)
			{
				L1 = L2;
				L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

				SetVertMontageSag (mp, p3_cab, p3_cab, H3_cab, L1, L2, CABLE_PROP+1); //...p3 mode -- for maximal vertical sag;
				for (j = 0; j < tab_sct->N_group; j++)
				if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1) {
					pp[0] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1); 
					pp[1] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
					pp[2] = 0.; 

					f0 = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3), pp);
					if (f0 < Size[0]) Size[0] = f0; 
				}
				mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);
			}
		}
		while (S_sect && S_sect > Size[0] && H7_cab < T1_abs_max && H1_cab < T1_abs_max);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		creep_transit(NUMCABLEMODETAB, H1_cab, p1_cab, t_aver, H_cab, t_montage, L_av, 0);
//		creep_transit(NUMCABLEMODETAB, H_cab,  p1_cab, t_montage, f0, t_aver, L_av, 1);

		SetTableParam(theApp.anker_list[NUMCABLEMONTAGETAB], i_ank, theApp.Opt.montage_temper, H_cab);
		CalculateCableMontageTable(i_ank);
	}

//////////////////////////////////////////////////////////////////////////////////
// ...ослабляем тяжение грозотроса до заданного габарита или до заданного тяжения;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 4) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMGROZOMODETAB])) {

//////////////////////////////////////////////
// ...вычисляем предельное тяжение грозотроса;
		double S_grozo = S_sect;
		H1_gro = MaxTension(NUMGROZOMODETAB, i1, i2);
		f0     = OptTension(NUMGROZOMODETAB, i1, i2);

		get_loading(NUMGROZOMODETAB, p1_gro, p3_gro, p7_gro, T2_abs_max, T2_opt, EF_gro, EF2_fin, alpha_gro, creep_gro);
		do { //...S_grozo минимальное расстояние до верхнего траверса в пролете;
			H1_gro *= DECREASE_STEP_ROUGH;
			transit(H1_gro, p1_gro, t_aver, EF_gro, H3_gro, p3_gro, t_glaze, EF2_fin, alpha_gro, L_av);

			if (S_grozo != 0.)
			for (mp[2] = Y0, L2 = 0., Size[0] = MAX_HIT, i = i1; i <= i2; i++)
			{
				L1 = L2;
				L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

				SetVertMontageSag(mp, p3_gro, p3_gro, H3_gro, L1, L2, GROZO_PROP+1);
				SetVertMontageSag(mp, p3_wr,  p3_wr,  H3_wr,  L1, L2, WIRES_PROP+(int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT));

				Size[0] = GetMontageDistance(mp, L1, L2, GROZO_PROP+1, WIRES_PROP+(int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT));

				mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);
			}
		}
		while ((S_grozo == 0. || S_grozo < Size[0]) && T2_opt < H1_gro);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		creep_transit(NUMGROZOMODETAB, H1_gro, p1_gro, t_aver, H_gro, t_montage, L_av, 0);
		SetTableParam(theApp.anker_list[NUMGROZOMONTAGETAB], i_ank, theApp.Opt.montage_temper, H_gro);
		CalculateGrozoMontageTable(i_ank);
	}
}

////////////////////////////////////////////////////////////////////////////////////
// ...контроль габарита до земли в одном участке без проверки максимального тяжения;
void AnkerEarth(int i_ank, int i, double & H_cab, double & H_wr1, double & S_cable, double & S_wire, int mode)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i+1, NUMCIRCUIT+8)));

	set_mode_index(0, i, NUMCABLEMODETAB);
	set_mode_index(0, i, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////////
// ...конструируем пролет в виде FEM модели;
	Nodes.SetNodes(0, NULL, 0);
	Elements.SetFElements(-1, & Nodes, & Properties);
	Elements.Numbering();
	Elements.NumGroups();
	Elements.Initials();

	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1.}; 
	AnkerConstructionGlobal(i, i, OK_STATE, (double *)(& mp));

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double pV_cab, pH_cab, pV_wr1, pH_wr1, L1 = 0., L2, del_earth,
			 Size; //...S_cable, S_wire -- минимальное расстояние до земли (для кабеля и провода);

/////////////////////
// ...расчета кабеля;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB])) {

//////////////////////////////////////////////
//...устанавливаем начальное состояния кабеля;
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode);
		H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1);

		L2 = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_earth = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

///////////////////////////////////////////////////////////////////////////
// ...натягиваем или ослабляем кабель до заданного габарита в i-ом пролете;
		if (fabs(S_cable) > ___EE___)
		do {
			H_cab *= INCREASE_STEP;

			SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);
			Size = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth);
		}
		while (S_cable > Size);
		do {
			H_cab *= DECREASE_STEP;

			SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);
			Size = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth);
		}
		while (fabs(S_cable) > ___EE___ && S_cable <= Size);
		H_cab /= DECREASE_STEP;

///////////////////////////////////////////////
// ...передаем назад реальный габарит до земли;
		SetVertMontageSag (mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);
		S_cable = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMCABLEMONTAGETAB, i_ank, H_cab, mode);
	}

////////////////////////////////////////////////////////////////
//...устанавливаем начальное состояния провода (расчет провода);
	set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr1, pH_wr1, mode);
	H_wr1 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+MAXGROZOS)*NUMSTENPOS);

	L2 = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
	del_earth = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

///////////////////////////////////////////////////////////////////////////
// ...натягиваем или ослабляем провод до заданного габарита в i-ом пролете;
	if (fabs(S_wire) > ___EE___)
	do {
		H_wr1 *= INCREASE_STEP;

		SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);
		Size = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth);
	}
	while (S_wire > Size);
	do {
		H_wr1 *= DECREASE_STEP;

		SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);
		Size = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth);
	}
	while (fabs(S_wire) > ___EE___ && S_wire <= Size);
	H_wr1 /= DECREASE_STEP;

///////////////////////////////////////////////
// ...передаем назад реальный габарит до земли;
	SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);
	S_wire = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
	set_montage_table(NUMWIRE1MONTAGETAB, i_ank, H_wr1, mode);
}

/////////////////////////////////////////////////////////////////////////////////////////
// ...подтягивание одного анкерного участка в сечении без проверки максимального тяжения;
void AnkerSection(int i_ank, int i, double & H_cab, double & H_wr1, double & S_cable, double & S_wire, int N_sect, int mode)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || ! tab_sct || ! theApp.anker_num) return;
	N_sect = max(1, N_sect);

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i+1, NUMCIRCUIT+8)));

	set_mode_index(0, i, NUMCABLEMODETAB);
	set_mode_index(0, i, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////////
// ...конструируем пролет в виде FEM модели;
	Nodes.SetNodes(0, NULL, 0);
	Elements.SetFElements(-1, & Nodes, & Properties);
	Elements.Numbering();
	Elements.NumGroups();
	Elements.Initials();

	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1.}; 
	AnkerConstructionGlobal(i, i, OK_STATE, (double *)(& mp));

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double pV_cab, pH_cab, pV_wr1, pH_wr1, L1 = 0., L2, del_earth, pp[3],
			 Size; //...S_cable, S_wire -- минимальное расстояние до пересечения (для кабеля и провода);

/////////////////////////////////////////////
// ...ищем пересечение в таблице пересечений;
	int id_first = N_sect, j;
	for (j = 0; id_first && j < tab_sct->N_group; j++)
		if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1)
			id_first--;
	j--;

/////////////////////
// ...расчета кабеля;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB]) && ! id_first) {

//////////////////////////////////////////////
//...устанавливаем начальное состояния кабеля;
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode);
		H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1);

		L2 = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_earth = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

///////////////////////////////////////////////////////////////////////////
// ...натягиваем или ослабляем кабель до заданного габарита в i-ом пролете;
		if (fabs(S_cable) > ___EE___)
		do {
			H_cab *= INCREASE_STEP;
			SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);

			pp[0] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1); if (pp[0] < 0.) pp[0] += L2;
			pp[1] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
			pp[2] = 0.; 
			Size  = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth, pp);
		}
		while (S_cable > Size);
		do {
			H_cab *= DECREASE_STEP;
			SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);

			pp[0] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1);  if (pp[0] < 0.) pp[0] += L2;
			pp[1] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
			pp[2] = 0.; 
			Size  = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth, pp);
		}
		while (fabs(S_cable) > ___EE___ && S_cable <= Size);
		H_cab /= DECREASE_STEP;

////////////////////////////////////////////////////////
// ...передаем назад реальное расстояние до пересечения;
		SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);
		pp[0]   = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1);  if (pp[0] < 0.) pp[0] += L2;
		pp[1]   = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
		pp[2]   = 0.; 
		S_cable = GetVertMontageSize(mp, L1, L2, CABLE_PROP+1, del_earth, pp);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMCABLEMONTAGETAB, i_ank, H_cab, mode);
	}

////////////////////
//...расчет провода;
	if (! id_first) {
///////////////////////////////////////////////
//...устанавливаем начальное состояния провода;
		set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr1, pH_wr1, mode);
		H_wr1 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+MAXGROZOS)*NUMSTENPOS);

		L2 = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_earth = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

///////////////////////////////////////////////////////////////////////////
// ...натягиваем или ослабляем провод до заданного габарита в i-ом пролете;
		if (fabs(S_wire) > ___EE___)
		do {
			H_wr1 *= INCREASE_STEP;
			SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);

			pp[0] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1);  if (pp[0] < 0.) pp[0] += L2;
			pp[1] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
			pp[2] = 0.; 
			Size  = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth, pp);
		}
		while (S_wire > Size);
		do {
			H_wr1 *= DECREASE_STEP;
			SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);

			pp[0] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1);  if (pp[0] < 0.) pp[0] += L2;
			pp[1] = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
			pp[2] = 0.; 
			Size  = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth, pp);
		}
		while (fabs(S_wire) > ___EE___ && S_wire <= Size);
		H_wr1 /= DECREASE_STEP;

////////////////////////////////////////////////////////
// ...передаем назад реальное расстояние до пересечения;
		SetVertMontageSag(mp, pV_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, L1, L2, WIRES_PROP+1);
		pp[0]  = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 1);  if (pp[0] < 0.) pp[0] += L2;
		pp[1]  = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, 2); 
		pp[2]  = 0.; 
		S_wire = GetVertMontageSize(mp, L1, L2, WIRES_PROP+1, del_earth, pp);

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMWIRE1MONTAGETAB, i_ank, H_wr1, mode);
	}
}

///////////////////////////////////////////////////////
// ...установка тяжения провода, кабеля или грозотроса;
void AnkerExtremeTension(int i_ank, double T_cable, double T_grozo, double T_wire)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (  ! table || ! theApp.anker_num) return;
	int i1 = theApp.anker_num[i_ank],
		 i2 = theApp.anker_num[i_ank+1]-1;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i1+1, NUMCIRCUIT+8)));

	set_mode_index(0, i1, NUMCABLEMODETAB);
	set_mode_index(0, i1, NUMGROZOMODETAB);
	set_mode_index(0, i1, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////////////////////////////////////
// ...определяем необходимые параметры для провода, кабеля и грозотроса; 
	double p1_wr,  p3_wr,  p6_wr,  p7_wr,  T0_opt, T0_max, T0_abs_max, T0_abs_tmin, EF_wr,  EF0_fin, EF0_creep, alpha_wr,  creep_wr,
			 p1_cab, p3_cab, p6_cab, p7_cab, T1_opt, T1_max, T1_abs_max, T1_abs_tmin, EF_cab, EF1_fin, EF1_creep, alpha_cab, creep_cab,
			 p1_gro, p3_gro, p6_gro, p7_gro, T2_opt, T2_max, T2_abs_max, T2_abs_tmin, EF_gro, EF2_fin, EF2_creep, alpha_gro, creep_gro;

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double L_av, H1_cab, H1_wr, H1_gro, f0;

///////////////////////////////////////////
// ...вычисляем предельное тяжение провода;
	get_loading(NUMWIRE1MODETAB, p1_wr, p3_wr, p6_wr, p7_wr, T0_abs_max, T0_abs_tmin, T0_opt, T0_max, EF0_creep, EF_wr, EF0_fin, alpha_wr, creep_wr);
	T0_max = GetTableParam(theApp.anker_list[NUMWIRE1MODETAB], NUMTENSION)*T_wire;
	L_av = L_average(i1, i2);

	if (T_wire > 0. && p1_wr != 0. && p6_wr != 0. && p7_wr != 0.) {
		transit(T0_opt = T0_max, p7_wr, t_glaze,   EF0_fin, f0, p1_wr, t_aver, EF0_fin, alpha_wr, L_av); H1_wr = f0;
		transit(T0_opt,		    p6_wr, t_rush,    EF0_fin, f0, p1_wr, t_aver, EF0_fin, alpha_wr, L_av); if (f0 < H1_wr) H1_wr = f0;
		transit(T0_opt,			 p1_wr, t_abs_min, EF0_fin, f0, p1_wr, t_aver, EF0_fin, alpha_wr, L_av); if (f0 < H1_wr) H1_wr = f0;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		creep_transit(NUMWIRE1MODETAB, H1_wr, p1_wr, t_aver, f0, t_montage, L_av, 0);
		SetTableParam(theApp.anker_list[NUMWIRE1MONTAGETAB], i_ank, theApp.Opt.montage_temper, f0);
		CalculateWireMontageTable(i_ank);
	}

//////////////////////////////////////////
// ...вычисляем предельное тяжение кабеля;
	get_loading(NUMCABLEMODETAB, p1_cab, p3_cab, p6_cab, p7_cab, T1_abs_max, T1_abs_tmin, T1_opt, T1_max, EF1_creep, EF_cab, EF1_fin, alpha_cab, creep_cab);
	T1_max = GetTableParam(theApp.anker_list[NUMCABLEMODETAB], NUMTENSION)*T_cable;

	if (T_cable > 0. && p1_cab != 0. && p6_cab != 0. && p7_cab != 0.) {
		transit(T1_opt = T1_max, p7_cab, t_glaze,   EF1_fin, f0, p1_cab, t_aver, EF1_fin, alpha_cab, L_av, 0., creep_cab); H1_cab = f0;
		transit(T1_opt,		    p6_cab, t_rush,    EF1_fin, f0, p1_cab, t_aver, EF1_fin, alpha_cab, L_av, 0., creep_cab); if (f0 < H1_cab) H1_cab = f0;
		transit(T1_opt,			 p1_cab, t_abs_min, EF1_fin, f0, p1_cab, t_aver, EF1_fin, alpha_cab, L_av, 0., creep_cab); if (f0 < H1_cab) H1_cab = f0;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		creep_transit(NUMCABLEMODETAB, H1_cab, p1_cab, t_aver, f0, t_montage, L_av, 0);
		SetTableParam(theApp.anker_list[NUMCABLEMONTAGETAB], i_ank, theApp.Opt.montage_temper, f0);
		CalculateCableMontageTable(i_ank);
	}

//////////////////////////////////////////////
// ...вычисляем предельное тяжение грозотроса;
	get_loading(NUMGROZOMODETAB, p1_gro, p3_gro, p6_gro, p7_gro, T2_abs_max, T2_abs_tmin, T2_opt, T2_max, EF2_creep, EF_gro, EF2_fin, alpha_gro, creep_gro);
	T2_max = GetTableParam(theApp.anker_list[NUMGROZOMODETAB], NUMTENSION)*T_grozo;

	if (T_grozo > 0. && p1_gro != 0. && p6_gro != 0. && p7_gro != 0.) {
		transit(T2_opt = T2_max, p7_gro, t_glaze,   EF2_fin, f0, p1_gro, t_montage, EF_gro, alpha_gro, L_av, 0., creep_gro); H1_gro = f0;
		transit(T2_opt,			 p6_gro, t_rush,    EF2_fin, f0, p1_gro, t_montage, EF_gro, alpha_gro, L_av, 0., creep_gro); if (f0 < H1_gro) H1_gro = f0;
		transit(T2_opt,			 p1_gro, t_abs_min, EF2_fin, f0, p1_gro, t_montage, EF_gro, alpha_gro, L_av, 0., creep_gro); if (f0 < H1_gro) H1_gro = f0;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		creep_transit(NUMGROZOMODETAB, H1_gro, p1_gro, t_aver, f0, t_montage, L_av, 0);
		SetTableParam(theApp.anker_list[NUMGROZOMONTAGETAB], i_ank, theApp.Opt.montage_temper, f0);
		CalculateGrozoMontageTable(i_ank);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// ...установка монтажной стрелы провеса в одном анкерном участке для  провода, кабеля и грозотроса;
void AnkerStrela(int i_ank, int i, double & H_cab, double & H_gro, double & H_wr1, 
					  double & f_cable, double & f_grozo, double & f_wire, int mode)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i+1, NUMCIRCUIT+8)));

	set_mode_index(0, i, NUMCABLEMODETAB);
	set_mode_index(0, i, NUMGROZOMODETAB);
	set_mode_index(0, i, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double pV_cab, pH_cab, pV_gro, pH_gro, pV_wr1, pH_wr1, L, del_h, Pf[3],
			 Size; //...f_cable, f_grozo, f_wire -- максимальная стрела для кабеля, грозотроса и провода;

/////////////////////
// ...расчета кабеля;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB])) {

//////////////////////////////////////////////
//...устанавливаем начальное состояния кабеля;
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode);
		H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1);

		L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 0);

//////////////////////////////////////////////////////////////////////
// ...ослабляем начальное состояния кабеля до заданной стрелы провеса;
		if (H_cab && fabs(f_cable) > ___EE___)
		do {
			H_cab *= DECREASE_STEP;

			bottom (L, del_h, pV_cab, H_cab, Pf);
			Size = strela0(L, del_h, pV_cab, H_cab, Pf);
		}
		while (f_cable > Size);
		do {
			H_cab *= INCREASE_STEP;

			if (H_cab) {
				bottom (L, del_h, pV_cab, H_cab, Pf);
				Size = strela0(L, del_h, pV_cab, H_cab, Pf);
			}
		}
		while (H_cab && fabs(f_cable) > ___EE___ && f_cable <= Size);
		H_cab /= INCREASE_STEP;

////////////////////////////////////////////////////
// ...передаем назад реальную стрелу провеса кабеля;
		if (H_cab) {
			bottom (L, del_h, pV_cab, H_cab, Pf);
			f_cable = strela0(L, del_h, pV_cab, H_cab, Pf);
		}
		else f_cable = 0.;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMCABLEMONTAGETAB, i_ank, H_cab, mode);
	}

/////////////////////////
// ...расчета грозотроса;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 4) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMGROZOMODETAB])) {

//////////////////////////////////////////////////
//...устанавливаем начальное состояния грозотроса;
		set_mode_tension(NUMGROZOMONTAGETAB, i_ank, pV_gro, pH_gro, mode);
		H_gro = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+MAXCABLES*NUMSTENPOS);

		L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, MAXCABLES*NUMSTENPOS);

//////////////////////////////////////////////////////////////
// ...ослабляем тяжение грозотроса до заданной стрелы провеса;
		if (H_gro && fabs(f_grozo) > ___EE___)
		do {
			H_gro *= DECREASE_STEP;

			bottom (L, del_h, pV_gro, H_gro, Pf);
			Size = strela0(L, del_h, pV_gro, H_gro, Pf);
		}
		while (f_grozo > Size);
		do {
			H_gro *= INCREASE_STEP;

			if (H_gro) {
				bottom (L, del_h, pV_gro, H_gro, Pf);
				Size = strela0(L, del_h, pV_gro, H_gro, Pf);
			}
		}
		while (H_gro && fabs(f_grozo) > ___EE___ && f_grozo <= Size);
		H_gro /= INCREASE_STEP;

////////////////////////////////////////////////////////
// ...передаем назад реальную стрелу провеса грозотроса;
		if (H_gro) {
			bottom (L, del_h, pV_gro, H_gro, Pf);
			f_grozo = strela0(L, del_h, pV_gro, H_gro, Pf);
		}
		else f_grozo = 0.;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMGROZOMONTAGETAB, i_ank, H_gro, mode);
	}

////////////////////////////////////////////////////////////////
//...устанавливаем начальное состояния провода (расчет провода);
	set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr1, pH_wr1, mode);
	H_wr1 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+MAXGROZOS)*NUMSTENPOS);

	L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
	del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, (MAXCABLES+MAXGROZOS)*NUMSTENPOS);

///////////////////////////////////////////////////////////
// ...ослабляем тяжение провода до заданной стрелы провеса;
	if (H_wr1 && fabs(f_wire) > ___EE___)
	do {
		H_wr1 *= DECREASE_STEP;

		bottom (L, del_h, pV_wr1, H_wr1, Pf);
		Size = strela0(L, del_h, pV_wr1, H_wr1, Pf);
	}
	while (f_wire > Size);
	do {
		H_wr1 *= INCREASE_STEP;

		if (H_wr1) {
			bottom (L, del_h, pV_wr1, H_wr1, Pf);
			Size = strela0(L, del_h, pV_wr1, H_wr1, Pf);
		}
	}
	while (H_wr1 && fabs(f_wire) > ___EE___ && f_wire <= Size);
	H_wr1 /= INCREASE_STEP;

/////////////////////////////////////////////////////
// ...передаем назад реальную стрелу провеса провода;
	if (H_wr1) {
		bottom (L, del_h, pV_wr1, H_wr1, Pf);
		f_wire = strela0(L, del_h, pV_wr1, H_wr1, Pf);
	}
	else f_wire = 0.;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
	set_montage_table(NUMWIRE1MONTAGETAB, i_ank, H_wr1, mode);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// ...установка монтажной стрелы провеса в одном анкерном участке для  провода, кабеля и грозотроса;
void AnkerPhasa(int i_ank, int i, double & H_cab, double & H_wr1, 
					  double & fw_cable, double & fw_wire, double & Sw_cable, int mode)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i+1, NUMCIRCUIT+8)));

	set_mode_index(0, i, NUMCABLEMODETAB);
	set_mode_index(0, i, NUMWIRE1MODETAB);

	set_pue_loading(i_ank);

////////////////////////////////////////////
// ...конструируем пролет в виде FEM модели;
	Nodes.SetNodes(0, NULL, 0);
	Elements.SetFElements(-1, & Nodes, & Properties);
	Elements.Numbering();
	Elements.NumGroups();
	Elements.Initials();

	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1.}; 
	AnkerConstructionGlobal(i, i, OK_STATE, (double *)(& mp));

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double pV_cab, pH_cab, pV_wr1, pH_wr1, L, del_h, Pf[3], S,
			 Size; //...fw_cable, fw_wire, Sw_cable -- стрела ветровой нагрузки и габарит для кабеля и провода;

////////////////////////////////////////////////////////////////
//...устанавливаем начальное состояния провода (расчет провода);
	set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr1, pH_wr1, mode);
	H_wr1 = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(MAXCABLES+MAXGROZOS)*NUMSTENPOS);

	L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
	del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, (MAXCABLES+MAXGROZOS)*NUMSTENPOS);

///////////////////////////////////////////////////////////
// ...ослабляем тяжение провода до заданной стрелы провеса;
	if (H_wr1 && fabs(pH_wr1) && fabs(fw_wire) > ___EE___)
	do {
		H_wr1 *= DECREASE_STEP;

		bottom (L, del_h, fabs(pH_wr1), H_wr1, Pf);
		Size = strela0(L, del_h, fabs(pH_wr1), H_wr1, Pf);
	}
	while (fw_wire > Size);
	do {
		H_wr1 *= INCREASE_STEP;

		if (H_wr1 && fabs(pH_wr1)) {
			bottom (L, del_h, fabs(pH_wr1), H_wr1, Pf);
			Size = strela0(L, del_h, fabs(pH_wr1), H_wr1, Pf);
		}
	}
	while (H_wr1 && fabs(pH_wr1) && fabs(fw_wire) > ___EE___ && fw_wire <= Size);
	H_wr1 /= INCREASE_STEP;

/////////////////////////////////////////////////////
// ...передаем назад реальную стрелу провеса провода;
	if (H_wr1 && fabs(pH_wr1)) {
		bottom (L, del_h, fabs(pH_wr1), H_wr1, Pf);
		fw_wire = strela0(L, del_h, fabs(pH_wr1), H_wr1, Pf);
	}
	else fw_wire = 0.;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
	set_montage_table(NUMWIRE1MONTAGETAB, i_ank, H_wr1, mode);

/////////////////////
// ...расчета кабеля;
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB])) {

//////////////////////////////////////////////
//...устанавливаем начальное состояния кабеля;
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode);
		H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1);

		L = GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);
		del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 0);

//////////////////////////////////////////////////////////////////////
// ...ослабляем начальное состояния кабеля до заданной стрелы провеса;
		if (H_cab && fabs(pH_cab) && fabs(fw_cable) > ___EE___)
		do {
			H_cab *= DECREASE_STEP;

			bottom (L, del_h, fabs(pH_cab), H_cab, Pf);
			Size = strela0(L, del_h, fabs(pH_cab), H_cab, Pf);
		}
		while (fw_cable > Size);
		do {
			H_cab *= INCREASE_STEP;

			if (H_cab && fabs(pH_cab)) {
				bottom (L, del_h, fabs(pH_cab), H_cab, Pf);
				Size = strela0(L, del_h, fabs(pH_cab), H_cab, Pf);
			}
		}
		while (H_cab && fabs(pH_cab) && fabs(fw_cable) > ___EE___ && fw_cable <= Size);
		H_cab /= INCREASE_STEP;

///////////////////////////////////////////////////////////////////////////
// ...натягиваем или ослабляем кабель до заданного габарита в i-ом пролете;
		int N_wires = min((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT), MAXCIRCUIT), j;
		for (j = 1; j <= N_wires; j++)
			SetHorzMontageSag(mp, pH_wr1, sqrt(sqr(pV_wr1)+sqr(pH_wr1)), H_wr1, 0., L, WIRES_PROP+j);

		if (H_cab && fabs(Sw_cable) > ___EE___)
		do {
			H_cab *= INCREASE_STEP;

			SetHorzMontageSag(mp, pH_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, 0., L, CABLE_PROP+1);

			for (Size = 0., j = 1; j <= N_wires; j++)
//			if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)*((j % 2)*2-1.)*(1.-(j == N_wires && N_wires % 2)*2.)) < Size || j == 1) Size = S; 
			if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)) < Size || j == 1) Size = S; 
		}
		while (Sw_cable > Size);
		do {
			H_cab *= DECREASE_STEP;

			if (H_cab)
			SetHorzMontageSag(mp, pH_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, 0., L, CABLE_PROP+1);

			for (Size = 0., j = 1; j <= N_wires; j++)
//			if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)*((j % 2)*2-1.)*(1.-(j == N_wires && N_wires % 2)*2.)) < Size || j == 1) Size = S; 
			if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)) < Size || j == 1) Size = S; 
		}
		while (H_cab && fabs(Sw_cable) > ___EE___ && Sw_cable <= Size);
		H_cab /= DECREASE_STEP;

////////////////////////////////////////////////////////////////////////
// ...передаем назад реальный габарит бокового ветра до фазных проводов;
		if (H_cab)
		SetHorzMontageSag(mp, pH_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, 0., L, CABLE_PROP+1);

		for (Size = 0., j = 1; j <= N_wires; j++)
//		if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)*((j % 2)*2-1.)*(1.-(j == N_wires && N_wires % 2)*2.)) < Size || j == 1) Size = S; 
		if ((S = GetHorzMontageSize(mp, 0, L, CABLE_PROP+1, WIRES_PROP+j, pH_cab < 0. ? -1. : 1.)) < Size || j == 1) Size = S; 
		Sw_cable = Size;

////////////////////////////////////////////////////
// ...передаем назад реальную стрелу провеса кабеля;
		if (H_cab && fabs(pH_cab)) {
			bottom (L, del_h, fabs(pH_cab), H_cab, Pf);
			fw_cable = strela0(L, del_h, fabs(pH_cab), H_cab, Pf);
		}
		else fw_cable = 0.;

//////////////////////////////////////////////////
// ...запоминаем полученную информацию о тяжениях;
		set_montage_table(NUMCABLEMONTAGETAB, i_ank, H_cab, mode);
	}
}
///////////////////////////////////////////////
// ...установка монтажного состояния в пролете;
void SetMontageState(int i1, int i2, int i_ank, int mode, double D_phase, int * anker_num)
{
	Table * table = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]);
	if (! table || ! theApp.anker_num) return;
	int i, j;

///////////////////////////////////////////////////////////////////////
//...устанавливаем климатический, нагрузочный индекс и нагрузку по ПУЭ;
	SetTableIndex(theApp.anker_list[NUMCLIMATETAB], max(1, (int)GetTableParam(table, i1+1, NUMCIRCUIT+8)));

	set_mode_index(0, i1, NUMCABLEMODETAB);
	set_mode_index(0, i1, NUMGROZOMODETAB);
	set_mode_index(0, i1, NUMWIRE1MODETAB);

	set_pue_loading(i_ank); //...делаем только одну итерацию для экономии времени;

////////////////////////////////////////////
// ...конструируем пролет в виде FEM модели;
	Nodes.SetNodes(0, NULL, 0);
	Elements.SetFElements(-1, & Nodes, & Properties);
	Elements.Numbering();
	Elements.NumGroups();
	Elements.Initials();

	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1.}; 
	AnkerConstructionGlobal(i1, i2, OK_STATE, (double *)(& mp));
	AnkerSectConstructionGlobal(mp, i1, i2);
	AnkerProfileConstructionGlobal(mp, i1, i2);

///////////////////////////////////////
//...вспомогательные параметры расчета;
	double H_wr, H_cab, H_gro, pV_wr, pH_wr, pV_cab, pH_cab, pV_gro, pH_gro, L1, L2;

//////////////////////////////////////////////////////////////////////
//...перенос состояния провода, грозотроса и кабеля в таблицу тяжений;
	set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr, pH_wr, mode);

	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB]))
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, mode); else
		set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, NULL_MODE);

	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 4) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 4) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMGROZOMODETAB]))
		set_mode_tension(NUMGROZOMONTAGETAB, i_ank, pV_gro, pH_gro, mode); else
		set_mode_tension(NUMGROZOMONTAGETAB, i_ank, pV_gro, pH_gro, NULL_MODE);

////////////////////////////
//...тест аварийного режима;
{/*
char * mark;
double p1, EF, alpha;
get_loading(NUMWIRE1MODETAB, p1, EF, alpha, mark);
transit_damage_anker_inverse(i1, i2, i1, 2, p1 = sqrt(sqr(pV_wr)+sqr(pH_wr)), EF);
/**/}
//...тест аварийного режима;
////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//...внесение провиса провода, грозотроса и кабеля в конечно-элементную модель;
	for (L2 = 0., i = i1; i <= i2; i++)
	{
		L1 = L2;
		L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

		if ((H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1)) != 0.) {
			SetVertMontageSag(mp, pV_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);
			SetHorzMontageSag(mp, pH_cab, sqrt(sqr(pV_cab)+sqr(pH_cab)), H_cab, L1, L2, CABLE_PROP+1);

////////////////////////////
//...тест аварийного режима;
{/*
char * mark;
double H_reduce = reduce_tension(sqrt(sqr(pV_cab)+sqr(pH_cab)), L2-L1, H_cab, 0.4), p1, EF, alpha, H;
get_loading(NUMCABLEMODETAB, p1, EF, alpha, mark);
transit_d(H_cab, p1, EF, H, L2-L1, 0.4);
/**/}
//...тест аварийного режима;
////////////////////////////
		}
		int N_wires = min((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT), MAXCIRCUIT);
		for (j = 0; j < N_wires; j++)
			if ((H_wr = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+(j+MAXCABLES+MAXGROZOS)*NUMSTENPOS)) != 0.) {
				SetVertMontageSag(mp, pV_wr, sqrt(sqr(pV_wr)+sqr(pH_wr)), H_wr, L1, L2, WIRES_PROP+1+j);
				SetHorzMontageSag(mp, pH_wr, sqrt(sqr(pV_wr)+sqr(pH_wr)), H_wr, L1, L2, WIRES_PROP+1+j);
		}
		if ((H_gro = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, 1+MAXCABLES*NUMSTENPOS)) != 0.) {
			SetVertMontageSag(mp, pV_gro, sqrt(sqr(pV_gro)+sqr(pH_gro)), H_gro, L1, L2, GROZO_PROP+1);
			SetHorzMontageSag(mp, pH_gro, sqrt(sqr(pV_gro)+sqr(pH_gro)), H_gro, L1, L2, GROZO_PROP+1);

////////////////////////////
//...тест аварийного режима;
{/*
char * mark;
double H_reduce = reduce_tension(sqrt(sqr(pV_gro)+sqr(pH_gro)), L2-L1, H_gro, 0.4), p1, EF, alpha, H;
get_loading(NUMGROZOMODETAB, p1, EF, alpha, mark);
transit_d(H_gro, p1, EF, H, L2-L1, 0.4);
/**/}
//...тест аварийного режима;
////////////////////////////
		}
	}

////////////////////////////////////////////////////////////////////
// ...вычисляем получающиеся габариты и расстояния до всех проводов;
	double dist[MAXCIRCUIT]; //...D_phase -- минимальное расстояние до ближайшего фазного провода;
	if (pV_cab)
	for (L2 = 0., i = i1; i <= i2; i++)
	{
		int N_wires = min((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT), MAXCIRCUIT);
		char * num1 = ((char **)table->table[0].parm_values)[i+1], 
			  * num2 = ((char **)table->table[0].parm_values)[i+2];
		L1 = L2;
		L2 = L1 + GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-5);

		memset(dist, 0, MAXCIRCUIT*sizeof(double));
		for (j = 1; j <= N_wires; j++)
			dist[j-1] = GetMontageDistance(mp, L1, L2, CABLE_PROP+1, WIRES_PROP+j);

		mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

		for (int j = 0; j < N_wires; j++) 
			if (0. < dist[j] && dist[j] <= D_phase && 0) {
				char buff[2000];
				sprintf(buff, "\ndist = %g", dist[j]);
				AfxMessageBox(CString("Не соблюдается габаритное расстояние до фазных проводов в пролете ") + 
								  CString(num1) + CString(" -- ") + CString(num2) + CString(buff));
			}

//////////////////////////////
// ...контроль бокового ветра;
		mp[2] -= GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

		double f = pH_wr < 0. ? -1. : 1.;
		memset(dist, 0, MAXCIRCUIT*sizeof(double));
		for (j = 1; j <= N_wires; j++)
//			dist[j-1] = GetHorzMontageSize(mp, L1, L2, CABLE_PROP+1, WIRES_PROP+j, f)*(g = (j % 2)*2-1.)*(h = 1.-(j == N_wires && N_wires % 2)*2.);
			dist[j-1] = GetHorzMontageSize(mp, L1, L2, CABLE_PROP+1, WIRES_PROP+j, f);

		mp[2] += GetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-3);

		for (j = 0; j < N_wires; j++) 
			if (dist[j] < ___EE_dop___ && 0) {
				char buff[2000];
				sprintf(buff, "\nSw = %g", dist[j]);
				AfxMessageBox(CString("Кабель и провода сближаются в пролете ") + 
								  CString(num1) + CString(" -- ") + CString(num2) + CString(buff));
			}
	}
}

////////////////////////////////////////////////////////
// ...визуализация линии в первый момент после загрузки;
void FirstVisualization(int id_montage_table)
{
	if (theView12 && theApp.anker_num && theApp.anker_num[0] > 0 && ! theApp.wrongVersion) {
/////////////////////////////////////////////////////////////////////////////////////
//...создание КЭ-модели электрической линии и визуализация первого анкерного участка;
		theApp.current_anker = 1;
		theApp.current_anker_span = 0;
		theApp.current_anker_span = min(theApp.current_anker_span, theApp.anker_num[0]);
		theApp.current_anker_span = max(theApp.current_anker_span, 0);

		SetMontageState(theView12->first_span = theApp.anker_num[theApp.current_anker], theView12->last_span = theApp.anker_num[theApp.current_anker+1]-1, 
							 theApp.current_anker, theApp.Opt.current_anker_mode);
		Message("Finish!");
		MessageMontageState(theApp.Opt.current_anker_mode);

		double H_cab = GetTableParam(theApp.anker_list[NUMCABLEMONTAGETAB], theApp.current_anker, theApp.Opt.montage_temper),
				 H_gro = GetTableParam(theApp.anker_list[NUMGROZOMONTAGETAB], theApp.current_anker, theApp.Opt.montage_temper),
				 H_wr1 = GetTableParam(theApp.anker_list[NUMWIRE1MONTAGETAB], theApp.current_anker, theApp.Opt.montage_temper);

		H_cab = ((int)(H_cab/M_G+.5));
		H_gro = ((int)(H_gro/M_G+.5));
		H_wr1 = ((int)(H_wr1/M_G+.5));

		char buff[2000];
		sprintf(buff, "Anker #%i: H_cab = %g, H_gro = %g, H_wr = %g", theApp.current_anker, H_cab, H_gro, H_wr1);
		Message(buff);
		Message("");
	}
	if (id_montage_table && theApp.anker_num) 
		for (int i = 1; i <= theApp.anker_num[0]; i++) { //...рабочее заполнение монтажных таблиц;
			CalculateCableMontageTable(i);
			CalculateGrozoMontageTable(i);
			CalculateWireMontageTable (i);
	}
}

//========================================================//
//                 РАСЧЕТ МОНТАЖНЫХ ТАБЛИЦ                //
//========================================================//
/////////////////////////////////////////
// ...расчет монтажных таблиц для кабеля;
void CalculateCableMontageTable(int i_ank)
{
//////////////////////////////////
// ...извлекаем монтажную таблицу;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMCABLEMONTAGETAB]);
	char  * cable_mark;
	if (! table || ! theApp.anker_num) return;

///////////////////////////////////////////////////////////////////////
// ...извлекаем параметры кабеля (и все необходимые режимы нагружения); 
	double p1_cab, EF_cab, alpha_cab;
	get_loading(NUMCABLEMODETAB, p1_cab, EF_cab, alpha_cab, cable_mark);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double H1, L_av, H, t_mtg, t_montage;
	theApp.Opt.montage_temper = max (theApp.Opt.montage_temper, 0);
	theApp.Opt.montage_temper = min (theApp.Opt.montage_temper, table->N-1);
	t_mtg = user_strtod(table->table[theApp.Opt.montage_temper].parm_name);

	if ((H1 = GetTableParam(theApp.anker_list[NUMCABLEMONTAGETAB], i_ank, theApp.Opt.montage_temper)) == 0.) return;
	L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], theApp.anker_num[i_ank]+1, NUMCIRCUIT-2); 

//////////////////////////////////////////////////////////
//...монтажное состояние кабеля при различной температуре;
	for (int k = 0; k < table->N; k++) {
		t_montage = user_strtod(table->table[k].parm_name);
		transit(H1, p1_cab, t_mtg, EF_cab, H, p1_cab, t_montage, EF_cab, alpha_cab, L_av);
		SetTableParam(theApp.anker_list[NUMCABLEMONTAGETAB], i_ank, k, H);
	}
}

/////////////////////////////////////////////
// ...расчет монтажных таблиц для грозотроса;
void CalculateGrozoMontageTable(int i_ank)
{
//////////////////////////////////
// ...извлекаем монтажную таблицу;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMGROZOMONTAGETAB]);
	char  * grozo_mark;
	if (! table || ! theApp.anker_num) return;

///////////////////////////////////////////////////////////////////////////
// ...извлекаем параметры грозотроса (и все необходимые режимы нагружения); 
	double p1_gro, EF_gro, alpha_gro;
	get_loading(NUMGROZOMODETAB, p1_gro, EF_gro, alpha_gro, grozo_mark);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double H2, L_av, H, t_mtg, t_montage;
	theApp.Opt.montage_temper =  max(theApp.Opt.montage_temper, 0);
	theApp.Opt.montage_temper =  min(theApp.Opt.montage_temper, table->N-1);
	t_mtg = user_strtod(table->table[theApp.Opt.montage_temper].parm_name);

	if ((H2 = GetTableParam(theApp.anker_list[NUMGROZOMONTAGETAB], i_ank, theApp.Opt.montage_temper)) == 0.) return;
	L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], theApp.anker_num[i_ank]+1, NUMCIRCUIT-2); 

//////////////////////////////////////////////////////////////
//...монтажное состояние грозотроса при различной температуре;
	for (int k = 0; k < table->N; k++) {
		t_montage = user_strtod(table->table[k].parm_name);
		transit(H2, p1_gro, t_mtg, EF_gro, H, p1_gro, t_montage, EF_gro, alpha_gro, L_av);
		SetTableParam(theApp.anker_list[NUMGROZOMONTAGETAB], i_ank, k, H);
	}
}

//////////////////////////////////////////
// ...расчет монтажных таблиц для провода;
void CalculateWireMontageTable(int i_ank)
{
//////////////////////////////////
// ...извлекаем монтажную таблицу;
	Table * table = (Table *)GetTable(theApp.anker_list[NUMWIRE1MONTAGETAB]);
	char  * wire_mark;
	if (! table || ! theApp.anker_num) return;

////////////////////////////////////////////////////////////////////////
// ...извлекаем параметры провода (и все необходимые режимы нагружения); 
	double p1_wr, EF_wr, alpha_wr;
	get_loading(NUMWIRE1MODETAB, p1_wr, EF_wr, alpha_wr, wire_mark);

////////////////////////////////////////
// ...вспомогательные параметры расчета;
	double H0, L_av, H, t_mtg, t_montage;
	theApp.Opt.montage_temper =  max(theApp.Opt.montage_temper, 0);
	theApp.Opt.montage_temper =  min(theApp.Opt.montage_temper, table->N-1);
	t_mtg = user_strtod(table->table[theApp.Opt.montage_temper].parm_name);

	if ((H0 = GetTableParam(theApp.anker_list[NUMWIRE1MONTAGETAB], i_ank, theApp.Opt.montage_temper)) == 0.) return;
	L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], theApp.anker_num[i_ank]+1, NUMCIRCUIT-2); 

///////////////////////////////////////////////////////////
//...монтажное состояние провода при различной температуре;
	for (int k = 0; k < table->N; k++) {
		t_montage = user_strtod(table->table[k].parm_name);
		transit(H0, p1_wr, t_mtg, EF_wr, H, p1_wr, t_montage, EF_wr, alpha_wr, L_av);
		SetTableParam(theApp.anker_list[NUMWIRE1MONTAGETAB], i_ank, k, H);
	}
}

// ==============================================================================================

// ==============================================================================================
// Новый подход к моделированию линии ===========================================================
// ==============================================================================================
//////////////////////////////////////////////////////
//...установка предельно возможного тяжения в пролете;
double MaxTension(Table * tab_mod, Table * tab_clm, int i1, int i2, int mode, int max_tension, Table * tab_montage)
{
//////////////////////////////////////////////////////////////
// ...извлекаем нагрузочные параметры и все режимы нагружения; 
	double p1, p3, p6, p7, T_opt, T_abs_max, T_abs_tmin, EF_creep, EF, EF_fin, alpha, creep, T_max, f0 = 0., H = 0., H1 = 0., H6 = 0., H7 = 0., ff;
	get_loading (tab_mod, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, L_av, L, del_h, t_grozo = 15.; 
	get_temper_mode(tab_clm, t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, tab_montage);

//////////////////////////////////////////////////////////////////////////////////////////////////
//...устанавливаем максимальную горизонтальную нагрузку (с учетом максимальной нагрузки на опору);
	int hyperb = 1;
	if (max_tension && p1 != 0. && p6 != 0. && p7 != 0. && T_abs_tmin != 0. && T_abs_max != 0.) {
		for (int i = i1; i <= i2; i++) {
			L	   = GetTableParam(theApp.anker_list[NUMANKERTAB],   i+1, NUMCIRCUIT-5);
			del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], i+1, (MAXCABLES+MAXGROZOS)*NUMSTENPOS);//...первый провод;
			ff = H_tension(T_abs_tmin, p1, L, del_h, hyperb); if (i == i1 || ff < H1) H1 = ff;
			ff = H_tension(T_abs_max,  p6, L, del_h, hyperb); if (i == i1 || ff < H6) H6 = ff;
			ff = H_tension(T_abs_max,  p7, L, del_h, hyperb); if (i == i1 || ff < H7) H7 = ff;
		}
		H1 = min(H1, T_max);
		H6 = min(H6, T_max);
		H7 = min(H7, T_max);
	}
	else {
		H1 = min(T_max, T_abs_tmin);
		H6 = 
		H7 = min(T_max, T_abs_max);
	}

///////////////////////////////////////////////////////////////////////////////////////////
// ...вычисляем тяжение провода, соответствующее разрывному и приводим его к нулю градусов;
	if ( (L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT-2)) == 0.)
			L_av = L_average(i1, i2, 1);

	if (p1 != 0. && p6 != 0. && p7 != 0.) {
		if (H1 > 0.) transit(H1, p1, t_abs_min, EF_fin, f0, p1, t_aver, EF_fin, alpha, L_av); H = f0;
		if (H6 > 0.) transit(H6, p6, t_rush,    EF_fin, f0, p1, t_aver, EF_fin, alpha, L_av); if (f0 < H) H = f0;
		if (H7 > 0.) transit(H7, p7, t_glaze,   EF_fin, f0, p1, t_aver, EF_fin, alpha, L_av); if (f0 < H) H = f0;
		if (H  > 0.) 
		switch (mode) {
			case			 P1_MODE: transit(f0 = H, p1, t_aver, EF_fin, H, p1, t_aver, EF_fin, alpha, L_av); break;
			case MIN_TEMPER_MODE: transit(f0 = H, p1, t_aver, EF_fin, H, p1, t_abs_min, EF_fin, alpha, L_av); break;
			case MAX_TEMPER_MODE: transit(f0 = H, p1, t_aver, EF_fin, H, p1, t_abs_max, EF_fin, alpha, L_av); break;
			case			 P3_MODE: transit(f0 = H, p1, t_aver, EF_fin, H, p3, t_glaze, EF_fin, alpha, L_av); break;
			case		 GROZO_MODE: transit(f0 = H, p1, t_aver, EF_fin, H, p1, t_grozo, EF_fin, alpha, L_av); break;
			case	  MONTAGE_MODE: creep_transit(tab_mod, tab_clm, f0 = H, p1, t_aver, H, t_montage, L_av, 0); break;
		}
	}
	return H;
}

double MaxTension(int num_modetab, int i1, int i2, int mode, int max_tension)
{
	Table * table; 
   if (! (table = (Table *)GetTable(theApp.anker_list[NUMCABLEMONTAGETAB])) &&
		 ! (table = (Table *)GetTable(theApp.anker_list[NUMGROZOMONTAGETAB])) &&
		 ! (table = (Table *)GetTable(theApp.anker_list[NUMWIRE1MONTAGETAB]))) table = table;
	return MaxTension((Table *)GetTable(theApp.anker_list[num_modetab]), 
							(Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), i1, i2, mode, max_tension, table);
}

//////////////////////////////////////////////////
//...установка эксплутационного тяжения в пролете;
double OptTension(int num_modetab, int i1, int i2, int mode)
{
//////////////////////////////////////////////////////////////
// ...извлекаем нагрузочные параметры и все режимы нагружения; 
	double p1, p3, p7, T_opt, T_abs_max, EF, EF_fin, alpha, creep, H = 0., f0;
	get_loading (num_modetab, p1, p3, p7, T_abs_max, T_opt, EF, EF_fin, alpha, creep);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, L_av; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

///////////////////////////////////////////////////////////////////////////////////////////
// ...вычисляем тяжение провода, соответствующее разрывному и приводим его к нулю градусов;
	if ( (L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT-2)) == 0.)
			L_av = L_average(i1, i2, 1);

	if (T_opt > 0.) 
	switch (mode) {
		case			 P1_MODE: transit(f0 = T_opt, p1, t_aver, EF_fin, H, p1, t_aver, EF_fin, alpha, L_av); break;
		case MIN_TEMPER_MODE: transit(f0 = T_opt, p1, t_aver, EF_fin, H, p1, t_abs_min, EF_fin, alpha, L_av); break;
		case MAX_TEMPER_MODE: transit(f0 = T_opt, p1, t_aver, EF_fin, H, p1, t_abs_max, EF_fin, alpha, L_av); break;
		case			 P3_MODE: transit(f0 = T_opt, p1, t_aver, EF_fin, H, p3, t_glaze, EF_fin, alpha, L_av); break;
		case	  MONTAGE_MODE: creep_transit(num_modetab, f0 = T_opt, p1, t_aver, H, t_montage, L_av, 0); break;
	}
	return H;
}

////////////////////////////////////////////////////////////////////
//...перевод тяжения из одного эксплуатационного состояния в другое;
double ModeTension(Table * tab_mod, Table * tab_clm, int i1, int i2, int mode_ini, int mode_fin, double Te)
{
//////////////////////////////////////////////////////////////
// ...извлекаем нагрузочные параметры и все режимы нагружения; 
	double p1, p3, p7, T_opt, T_abs_max, EF, EF_fin, alpha, creep, H = 0., f0, p, temper;
	get_loading (tab_mod, p1, p3, p7, T_abs_max, T_opt, EF, EF_fin, alpha, creep);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage, L_av; 
	get_temper_mode(tab_clm, t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);//...t_montage = 0;

///////////////////////////////////////////////////////////////////////////////////////////
// ...вычисляем тяжение провода, соответствующее разрывному и приводим его к нулю градусов;
	if ( (L_av = GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT-2)) == 0.)
			L_av = L_average(i1, i2, 1);

	switch (mode_ini) {
		case			 P1_MODE: p = p1; temper = t_aver; break;
		case MIN_TEMPER_MODE: p = p1; temper = t_abs_min; break;
		case MAX_TEMPER_MODE: p = p1; temper = t_abs_max; break;
		case			 P3_MODE: p = p3; temper = t_glaze; break;
	}
	if (Te > 0.) 
	switch (mode_fin) {
		case			 P1_MODE: transit(f0 = Te, p, temper, EF_fin, H, p1, t_aver, EF_fin, alpha, L_av); break;
		case MIN_TEMPER_MODE: transit(f0 = Te, p, temper, EF_fin, H, p1, t_abs_min, EF_fin, alpha, L_av); break;
		case MAX_TEMPER_MODE: transit(f0 = Te, p, temper, EF_fin, H, p1, t_abs_max, EF_fin, alpha, L_av); break;
		case			 P3_MODE: transit(f0 = Te, p, temper, EF_fin, H, p3, t_glaze, EF_fin, alpha, L_av); break;
	}
	return H;
}

///////////////////////////////////////////////////////////////////////////////////////////
//...функция вытяжки провода в пролете (переход из состояния до вытяжки в состояние после);
void creep_transit(Table * tab_mod, Table * tab_clm, double H, double p, double t, double & H_norm, double t_norm, double L, int m)
{
//////////////////////////////////////////////////////////////
// ...извлекаем нагрузочные параметры и все режимы нагружения; 
	double p1, p3, p6, p7, T_opt, T_abs_max, T_abs_tmin, EF_creep, EF, EF_fin, alpha, creep, T_max, H_max, p_max;
	get_loading (tab_mod, p1, p3, p6, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF_creep, EF, EF_fin, alpha, creep);
	T_max = max(T_max, T_abs_max);
	T_max = max(T_max, T_abs_tmin);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, ttt; 
	get_temper_mode(tab_clm, t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, ttt); //...ttt - не используется;

//////////////////////////////////////////////////////////
// ...осужествляем вытяжку провода, кабеля или грозотроса;
#ifndef ___CREEP_D_MODULE___
	if (p1 != 0. && p6 != 0. && p7 != 0.)
	if (m) {
//		double t_max;
//		transit(H, p, t, EF_creep, f0, p_max = p7, t_max = t_glaze, EF_creep, alpha, L);		H_max = f0;
//		transit(H, p, t, EF_creep, f0, p6, t_rush,    EF_creep, alpha, L); if (f0 > H_max) {H_max = f0; p_max = p6; t_max = t_rush;}
//		transit(H, p, t, EF_creep, f0, p1, t_abs_min, EF_creep, alpha, L); if (f0 > H_max) {H_max = f0; p_max = p1; t_max = t_abs_min;}
//		transit(H_max, p_max, t_max, EF_fin, H_norm, p, t_norm, EF_fin, alpha, L);

		transit_p(H, p, t, EF_creep, T_max, p_max, t_aver, EF_creep, alpha, L);
		transit  (H, p, t, EF_creep, H_max, p_max, t_aver, EF_creep, alpha, L);
		transit(H_max, p_max, t_aver, EF_fin, H_norm, p, t_norm, EF_fin, alpha, L);
	}
	else {
//		transit(H, p, t, EF_fin, f0, p_max = p7, t_max = t_glaze, EF_fin, alpha, L); H_max = f0;
//		transit(H, p, t, EF_fin, f0, p6, t_rush,    EF_fin, alpha, L); if (f0 > H_max) {H_max = f0; p_max = p6; t_max = t_rush;}
//		transit(H, p, t, EF_fin, f0, p1, t_abs_min, EF_fin, alpha, L); if (f0 > H_max) {H_max = f0; p_max = p1; t_max = t_abs_min;}
//		transit(H_max, p_max, t_max, EF_creep, H_norm, p, t_norm, EF_creep, alpha, L);

		transit_p(H, p, t, EF_fin, T_max, p_max, t_aver, EF_fin, alpha, L);
		transit  (H, p, t, EF_fin, H_max, p_max, t_aver, EF_fin, alpha, L);
		transit(H_max, p_max, t_aver, EF_creep, H_norm, p, t_norm, EF_creep, alpha, L);
	}
#else
////////////////////////////////////////////////////////
// ...другая модель вытяжки -- по остаточной деформации;
	if (p1 != 0. && p6 != 0. && p7 != 0.)
	if (m) {
//		transit_p(H, p, t, EF, T_max, p_max, t_aver, EF, alpha, L);
//		transit  (H, p, t, EF, H_max, p_max, t_aver, EF_fin, alpha, L, creep);
//		transit(H_max, p_max, t_aver, EF_fin, H_norm, p, t_norm, EF_fin, alpha, L);

		transit(H, p, t, EF, H_norm, p, t_norm, EF_fin, alpha, L, creep);
	}
	else {
//		transit_p(H, p, t, EF_fin, T_max, p_max, t_aver, EF, alpha, L, 0., creep);
//		transit  (H, p, t, EF_fin, f0, p_max, t_aver, EF_fin, alpha, L);
//		transit(f0, p_max, t_aver, EF_fin, H_max, p_max, t_aver, EF, alpha, L, 0., creep);
//		transit(H_max, p_max, t_aver, EF, H_norm, p, t_norm, EF, alpha, L);

		transit(H, p, t, EF_fin, H_norm, p, t_norm, EF, alpha, L, 0., creep);
	}
#endif
}
void creep_transit(int num_modetab, double H, double p, double t, double & H_norm, double t_norm, double L, int m)
{
	creep_transit ((Table *)GetTable(theApp.anker_list[num_modetab]), 
						(Table *)GetTable(theApp.anker_list[NUMCLIMATETAB]), H, p, t, H_norm, t_norm, L, m);
}

///////////////////////////////////////////
//...определение критической длины пролета;
double L_Transit(int num_modetab)
{
//////////////////////////////////////////////////////////////
// ...извлекаем нагрузочные параметры и все режимы нагружения; 
	double p1, p7, T_opt, T_abs_max, T_abs_tmin, EF, EF_fin, alpha, creep, T_max, H = 0.;
	get_loading (num_modetab, p1, p7, T_abs_max, T_abs_tmin, T_opt, T_max, EF, EF_fin, alpha, creep);

////////////////////////////////////////////////
// ...извлекаем температурные режимы нагружения; 
	double t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage; 
	get_temper_mode(t_glaze, t_rush, t_aver, t_abs_max, t_abs_min, t_montage);

///////////////////////////////////
 // ...расчет критического пролета;
	double L_tr = L_transit(T_abs_tmin, p1, t_abs_min, EF_fin, T_abs_max, p7, t_glaze, alpha);
	return L_tr;
}

//////////////////////////////////////////////////////////
//...выравнивание опор согласно данным линейных измерений;
void TowerJump(int i1, int i2, int id_cross)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || i1 < 0 || i1 > table->N_group-2 ||
						i2 < 0 || i2 > table->N_group-2 || ! tab_sct) return;
	double par     [NUM_OF_TOWER_BASE_ELEMENTS+19],
			 par_prev[NUM_OF_TOWER_BASE_ELEMENTS+19];

////////////////////////////////////////
//...образуем модель избранных пролетов;
	int i, j, l, id_first = 1, id_i1 = 0, id_i2 = 0,
		 N_circuit = CIRCUIT((int)GetTableParam(theApp.anker_list[NUMANKERTAB], i1+1, NUMCIRCUIT));

/////////////////////////////////////////////////////////
//...устанавливаем начальные значения параметров подвески;
	memset(par+NUM_OF_TOWER_BASE_ELEMENTS, 0, 18*sizeof(double));
	par[NUM_OF_TOWER_BASE_ELEMENTS+3]  = 1+N_circuit;
	par[NUM_OF_TOWER_BASE_ELEMENTS+10] = 1;
	par[NUM_OF_TOWER_BASE_ELEMENTS+14] = 2;
	memcpy(par_prev+NUM_OF_TOWER_BASE_ELEMENTS, par+NUM_OF_TOWER_BASE_ELEMENTS, 18*sizeof(double));

/////////////////////////
//...конструируем пролет;
	SetTableIndex(theApp.anker_list[NUMANKERHANGINGTAB], 0);
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);

	for (i = i1; i <= i2; i++, memcpy(par_prev, par, (NUM_OF_TOWER_BASE_ELEMENTS+19)*sizeof(double)))
	for (j = 0; j < tab_sct->N_group; j++)
	if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1) {

		if (theApp.anker_num)
		for (id_i1 = id_i2 = l = 1; (id_i1 || id_i2) && l <= theApp.anker_num[0]; l++) {
			if (i == theApp.anker_num[l]) id_i1 = 0;
			if (i+1 == theApp.anker_num[l+1]) id_i2 = 0;
		}
		int num_hangingtab = id_i1 ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
		set_hanging_index(0, i, NUMANKERHANGINGTAB);
		set_hanging_index(0, i, NUMSTANDHANGINGTAB);

		if (GetTable(theApp.anker_list[num_hangingtab]))
		for (l = 0; l < 18; l++) //...извлекаем параметры подвески;
			par_prev[NUM_OF_TOWER_BASE_ELEMENTS+l] = GetTableParam(theApp.anker_list[num_hangingtab], l+1);

		TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+1], par_prev, theApp.anker_list[NUMTOWERTAB]);
		double hunging  = par_prev[(int)par_prev[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4]*3+1]-par_prev[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4+2],
				 jump = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, NUMSCTTYPE-5);
		if (jump)
		SetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-4, jump-hunging);

		num_hangingtab  = id_i2 ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
		set_hanging_index(i+1, i+1, NUMANKERHANGINGTAB);
		set_hanging_index(i+1, i+1, NUMSTANDHANGINGTAB);

		if (GetTable(theApp.anker_list[num_hangingtab]))
		for (l = 0; l < 18; l++) //...извлекаем параметры подвески;
			par[NUM_OF_TOWER_BASE_ELEMENTS+l] = GetTableParam(theApp.anker_list[num_hangingtab], l+1);

		TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+2], par, theApp.anker_list[NUMTOWERTAB]);
		hunging = par[(int)par[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4]*3+1]-par[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4+2],
		jump = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, NUMSCTTYPE-4);
		if (jump)
		SetTableParam(theApp.anker_list[NUMANKERTAB], i+2, NUMCIRCUIT-4, jump-hunging);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...другой алгоритм выравнивания на основе точек подвески -- работает как-то странно (пока не используем);
void TowerJump2(int i1, int i2, int id_cross)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || i1 < 0 || i1 > table->N_group-2 ||
						i2 < 0 || i2 > table->N_group-2 || ! tab_sct) return;
	double par     [NUM_OF_TOWER_BASE_ELEMENTS+19],
			 par_prev[NUM_OF_TOWER_BASE_ELEMENTS+19];
	CMap mp[] = {1., 0., 0., 0., 0., 0., 0., -1. };

////////////////////////////////////////
//...образуем модель избранных пролетов;
	int i, j, l, id_first = 1, id_i1 = 0, id_i2 = 0;

/////////////////////////////////////////////////////////
//...устанавливаем начальные значения параметров подвески;
	memset(par+NUM_OF_TOWER_BASE_ELEMENTS, 0, 18*sizeof(double));
	memcpy(par_prev+NUM_OF_TOWER_BASE_ELEMENTS, par+NUM_OF_TOWER_BASE_ELEMENTS, 18*sizeof(double));

/////////////////////////
//...конструируем пролет;
	SetTableIndex(theApp.anker_list[NUMANKERHANGINGTAB], 0);
	SetTableIndex(theApp.anker_list[NUMSTANDHANGINGTAB], 0);

	for (i = i1; i <= i2; i++, memcpy(par_prev, par, (NUM_OF_TOWER_BASE_ELEMENTS+19)*sizeof(double)))
	for (j = 0; j < tab_sct->N_group; j++)
	if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1) {

		if (theApp.anker_num)
		for (id_i1 = id_i2 = l = 1; (id_i1 || id_i2) && l <= theApp.anker_num[0]; l++) {
			if (i == theApp.anker_num[l]) id_i1 = 0;
			if (i+1 == theApp.anker_num[l+1]) id_i2 = 0;
		}
		int num_hangingtab = id_i1 ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
		set_hanging_index(0, i, NUMANKERHANGINGTAB, 1);
		set_hanging_index(0, i, NUMSTANDHANGINGTAB, 1);

		if (GetTable(theApp.anker_list[num_hangingtab]))
		for (l = 0; l < 18; l++) //...извлекаем параметры подвески;
			par_prev[NUM_OF_TOWER_BASE_ELEMENTS+l] = GetTableParam(theApp.anker_list[num_hangingtab], l+1);

		TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+1], par_prev, theApp.anker_list[NUMTOWERTAB]);

///////////////////////////////////////////////////
//...вычисляем точки подвески провода и грозотроса;
		int N1Wires_max, N1Grozo_max;
		double P1Wires[18], P1Grozo[3], hunging, jump;
		HungingPoints(par_prev, P1Wires, N1Wires_max = 6, P1Grozo, N1Grozo_max = 1, mp);
		hunging = (id_cross*2 < N1Wires_max ? P1Wires[id_cross*6+1] : P1Grozo[1])-par_prev[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4+2];
		if ((jump = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, NUMSCTTYPE-5)) != 0.)
		SetTableParam(theApp.anker_list[NUMANKERTAB], i+1, NUMCIRCUIT-4, jump-hunging);

		num_hangingtab  = id_i2 ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
		set_hanging_index(i+1, i+1, NUMANKERHANGINGTAB);
		set_hanging_index(i+1, i+1, NUMSTANDHANGINGTAB);

		if (GetTable(theApp.anker_list[num_hangingtab]))
		for (l = 0; l < 18; l++) //...извлекаем параметры подвески;
			par[NUM_OF_TOWER_BASE_ELEMENTS+l] = GetTableParam(theApp.anker_list[num_hangingtab], l+1);

		TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+2], par, theApp.anker_list[NUMTOWERTAB]);

///////////////////////////////////////////////////
//...вычисляем точки подвески провода и грозотроса;
		int N2Wires_max, N2Grozo_max;
		double P2Wires[18], P2Grozo[3];
		HungingPoints(par, P2Wires, N2Wires_max = 6, P2Grozo, N2Grozo_max = 1, mp);
		hunging = (id_cross*2 < N2Wires_max ? P1Wires[id_cross*6+1] : P2Grozo[1])-par[NUM_OF_TOWER_BASE_ELEMENTS+6+id_cross*4+2];
		if ((jump = GetTableParam(theApp.anker_list[NUMSECTIONTAB], j+1, NUMSCTTYPE-4)) != 0.)
		SetTableParam(theApp.anker_list[NUMANKERTAB], i+2, NUMCIRCUIT-4, jump-hunging);
	}
}

//////////////////////////////////////////////////
//...вычисление запаса до конкретного пересечения;
double SectCableVGap(int i, int mode)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || ! tab_sct || ! theApp.anker_num 
					|| i < 0 || i >= tab_sct->N_group) return(0.);

/////////////////////
// ...расчета кабеля;
	double pV_cab, pH_cab, H_cab, t_norm, L, del_h, Pf[3], gap = 0.;

	SetTableIndex(theApp.anker_list[NUMCABLEMODETAB], 0);
	if ((! GetTable(theApp.anker_list[NUMANKERHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMANKERHANGINGTAB], 1) != ERR_STATE) &&
		 (! GetTable(theApp.anker_list[NUMSTANDHANGINGTAB]) || (int)GetTableParam(theApp.anker_list[NUMSTANDHANGINGTAB], 1) != ERR_STATE) &&
			 GetTable(theApp.anker_list[NUMCABLEMODETAB]))
		for (int i_ank = 0, j = 0; j < table->N_group; j++) {
			set_mode_index(j, j, NUMCABLEMODETAB);

			if (theApp.anker_num[0] > i_ank && j >= theApp.anker_num[i_ank+1]) i_ank++;

			if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1) {
				L = GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-5); 
				del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], j+1, 0); 

////////////////////////////////////
//...устанавливаем состояния кабеля;
				set_mode_tension(NUMCABLEMONTAGETAB, i_ank, pV_cab, pH_cab, OPERATING_MODE, & (t_norm = GetTableParam(theApp.anker_list[NUMCLIMATETAB], i+1, 6)));
				if ((H_cab = GetTableParam(theApp.anker_list[NUMTENSIONTAB], j+1, 1)) == 0.) return gap;

////////////////////////////////////////////////
//...устанавливаем значения параметров подвески;
				double par[NUM_OF_TOWER_BASE_ELEMENTS], cable_cross, cable_point, hanging;
				int num_hangingtab = j != theApp.anker_num[i_ank] ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
				set_hanging_index(0, j, NUMSTANDHANGINGTAB);
				set_hanging_index(0, j, NUMANKERHANGINGTAB);

				TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+1], par, theApp.anker_list[NUMTOWERTAB]);

				cable_cross = GetTableParam(theApp.anker_list[num_hangingtab], 1);
				cable_point = GetTableParam(theApp.anker_list[num_hangingtab], 2);
				hanging = par[(int)cable_cross*3+1]-cable_point+GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-4);

/////////////////////////////////////
// ...вычисляем запас до пересечения;
				double d = GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-3)/L, f;
				if ((f = GetTableParam(theApp.anker_list[NUMSECTIONTAB], i+1, 1)) <  ___EE___)  f += L;
				bottom(L, del_h, pV_cab, H_cab, Pf);
				gap = hanging+sag(f, pV_cab, H_cab, Pf)-GetTableParam(theApp.anker_list[NUMSECTIONTAB], i+1, 2)-f*d;
			}
	}
	return gap;
}

//////////////////////////////////////////////////
//...вычисление запаса до конкретного пересечения;
double SectWireVGap(int i, int mode, int id_cross)
{
///////////////////////////////////
// ...извлекаем таблицу всей линии;
	Table * table   = (Table *)GetTable(theApp.anker_list[NUMANKERTAB]),
			* tab_sct = (Table *)GetTable(theApp.anker_list[NUMSECTIONTAB]);
	if (! table || ! tab_sct || ! theApp.anker_num 
					|| i < 0 || i >= tab_sct->N_group) return(0.);

/////////////////////
// ...расчета провода;
	double pV_wr, pH_wr, H_wr, t_norm, L, del_h, Pf[3], gap = 0.;

	SetTableIndex(theApp.anker_list[NUMWIRE1MODETAB], 0);
	for (int i_ank = 0, j = 0; j < table->N_group; j++) {
		set_mode_index(j, j, NUMWIRE1MODETAB);

		if (theApp.anker_num[0] > i_ank && j >= theApp.anker_num[i_ank+1]) i_ank++;

		if (((int *)tab_sct->table[0].parm_values)[j+1] == i+1) {
			L = GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-5); 
			del_h = GetTableParam(theApp.anker_list[NUMTENSIONTAB], j+1, (id_cross*2+MAXCABLES+MAXGROZOS)*NUMSTENPOS); 

////////////////////////////////////
//...устанавливаем состояния кабеля;
			set_mode_tension(NUMWIRE1MONTAGETAB, i_ank, pV_wr, pH_wr, OPERATING_MODE, & (t_norm = GetTableParam(theApp.anker_list[NUMCLIMATETAB], i+1, 6)));
			if ((H_wr = GetTableParam(theApp.anker_list[NUMTENSIONTAB], j+1, 1+(id_cross*2+MAXCABLES+MAXGROZOS)*NUMSTENPOS)) == 0.) return gap;

////////////////////////////////////////////////
//...устанавливаем значения параметров подвески;
			double par[NUM_OF_TOWER_BASE_ELEMENTS], wire_cross, wire_point, hanging, h_sect, b_sect, f_sect;
			int num_hangingtab = j != theApp.anker_num[i_ank] ? NUMSTANDHANGINGTAB : NUMANKERHANGINGTAB;
			set_hanging_index(0, j, NUMSTANDHANGINGTAB);
			set_hanging_index(0, j, NUMANKERHANGINGTAB);

			TowerBase(theApp.tower_list, ((char **)table->table[1].parm_values)[i+1], par, theApp.anker_list[NUMTOWERTAB]);

			wire_cross = GetTableParam(theApp.anker_list[num_hangingtab], 7+id_cross*4);
			wire_point = GetTableParam(theApp.anker_list[num_hangingtab], 9+id_cross*4);
			hanging = par[(int)wire_cross*3+1]-wire_point+GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-4);

/////////////////////////////////////
// ...вычисляем запас до пересечения;
			double d = -GetTableParam(theApp.anker_list[NUMANKERTAB], j+1, NUMCIRCUIT-3)/L, f;
			if ((f = GetTableParam(theApp.anker_list[NUMSECTIONTAB], i+1, 1)) <  ___EE___)  f += L;
			bottom(L, del_h, pV_wr, H_wr, Pf);
			f_sect = -sag(f, pV_wr, H_wr, Pf)- del_h/L*f;
			h_sect = hanging+sag(f, pV_wr, H_wr, Pf)+d*f;
			gap = h_sect-(b_sect = GetTableParam(theApp.anker_list[NUMSECTIONTAB], i+1, 2));
		}
	}
	return gap;
}

