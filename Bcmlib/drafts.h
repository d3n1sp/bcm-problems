/*=========================================*/
/*                  DRAFTS                 */
/*=========================================*/
#ifndef ___DRAFTS___
#define ___DRAFTS___

#include <typeinfo>

#include "clame2d.h"
#include "cheat2d.h"
#include "cheat3d.h"
#include "ccohes2d.h"
#include "ccohes3d.h"
#include "cmindl2d.h"
#include "cmindl3d.h"
#include "chydro3d.h"
#include "cvisco2d_grad.h"

///////////////////////////////////////////////////////////////////////
//              ALL EXISTING TYPES OF DRAFT PROBLEMS                 //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...global functions for construction of all existing types of drafts;
template <typename T> 
CDraft<T> * CreateDraft(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8)
{
	switch (id_DRAFT) {
      case			  ERR_DRAFT: return NULL;
      case		  LAME2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CLame2D(id_dop));
      case		  LAME3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CLame3D(id_dop));
      case		  HEAT2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CHeat2D(id_dop));
      case		  HEAT3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CHeat3D(id_dop));
      case		 HYDRO3D_DRAFT: if (typeid(T) != typeid(complex)) return (CDraft<T> *)(new CHydro3D<T>(id_dop));
      case		 COHES2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CCohes2D(id_dop));
      case		 COHES3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CCohes3D(id_dop));
      case		 MINDL2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CMindl2D(id_dop));
      case		 MINDL3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CMindl3D(id_dop));
      case VISCO2D_GRAD_DRAFT: if (typeid(T) == typeid(complex)) return (CDraft<T> *)(new CVisco2D_grad(id_dop));
	}
   return new CDraft<T>;
}

////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR ESHELBY PROBLEMS          //
////////////////////////////////////////////////////////////////
double TakeLayer_EH(int N, double * ff, double * kk, double * ll)
{
	CCohes3D draft;
	CLame3D ldraft; double EH;
	if (ll == NULL) EH = ldraft.TakeLayer_kk(N, ff, kk);
	else				 EH =  draft.TakeLayer_kk(N, ff, kk, ll);
	return(EH);
}

double TakeCylinder_EH(int N, double * ff, double * kk)
{
	CLame2D ldraft; double EH;
	EH = ldraft.TakeLayer_kk(N, ff, kk);
	return(EH);
}

double TakeCylinder_KH(int N, double * ff, double * kp, double * mu)
{
	CLame2D ldraft; double KH;
	KH = ldraft.TakeLayer_kk(N, ff, kp, mu);
	return(KH);
}

double TakeCylinder_GH(int N, double * ff, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, mu, nj);
	return(GH);
}

double TakeCylinder_GH(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, kp, mu, nj);
	return(GH);
}

double TakeCylinder_SH(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_sh(N, ff, kp, mu, nj);
	return(GH);
}

double TakeCylinder_LH(int N, double * ff, double * lm, double * mu, double * nj)
{
	CLame2D ldraft; double LH;
	LH = ldraft.TakeLayer_kk(N, ff, lm, mu, nj);
	return(LH);
}

double TakeSphere_KH(int N, double * ff, double * kv, double * mu)
{
	CLame3D ldraft; double KH;
	KH = ldraft.TakeLayer_kk(N, ff, kv, mu);
	return(KH);
}

double TakeSphere_GH(int N, double * ff, double * kv, double * mu, double * nj)
{
	CLame3D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, kv, mu, nj);
	return(GH);
}

double TakeSphere_GH_det(double c0, double nju1, double nju2, double E1, double E2, double alpha)
{
	CLame3D ldraft; double GH;
	ldraft.set_fasa_hmg(nju2, nju1, E2/(1.+nju2)*.5, E1/(1.+nju1)*.5);
	GH = ldraft.TakeEshelby_shear_det(c0, alpha);
	return(GH);
}

double TakeSphere_volm_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2, double AA)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		draft.set_adhesion(AA, 0.);
		KH = draft.TakeEshelby_volm_two(c0);
	}
	return(KH);
}

double TakeSphere_shear_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2, double AA, double BB)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		GH = ldraft.TakeEshelby_shear_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		draft.set_adhesion(AA, BB);
		GH = draft.TakeEshelby_shear_two(c0);
	}
	return(GH);
}

double TakeSphere_volm_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		KH = draft.TakeEshelby_volm_sym(c0);
	}
	return(KH);
}

double TakeSphere_shear_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 1e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) <= 0.01 && fabs(l2) <= 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear_old(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}

double TakeSphere_shear(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 2e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) < 0.01 && fabs(l2) < 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция аппроксимации межфазного слоя ll = 1./sqrt(C0+sqr(t)*(C1+t*(C2+...+t*CN)...), coef[N+1], data[2*N+2];
int strata_approx(int N, double * data, double * coef)
{
	int  * ii = (int *)new_struct((N+1)*sizeof(int)), i, k, l, k0, l0;
	double f, ** matrix = 0; set_matrix(matrix, N+1, N+1);

/////////////////////////////////////////////////
//...переносим заданные значения межфазного слоя;
	for (i = 0; i <= N; i++) coef[i] = 1./sqr(data[2*i+1]); 

/////////////////////////////////////////////////////////////////
//...заполняем систему линейных уравнений с матрицей Вандермонда;
	for (i = 0; i <= N; i++)
	for (f = 1., k = 0; k <= N; k++, f *= data[2*i], f *= (k == 1 ? data[2*i] : 1.)) matrix[i][k] = f; 

/////////////////////////////////////////////////////////////////////////////////////////////
//...определяем коэффициенты путем решения системы линейных уравнений с матрицей Вандермонда;
	for (i = 0; i <= N; i++) {
		for (f = 0., k = 0; k <= N; k++) //...look for position maximal element;
			if (ii[k] != 1) 
				for (l = 0; l <= N; l++) 
					if (! ii[l]) {
						if (fabs(matrix[k][l]) >= f) f = fabs(matrix[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) {
						delete_struct(ii); delete_struct(matrix);	return(-1);
					}
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0)
			for (swap(coef[k0], coef[l0]), l = 0; l <= N; l++) swap(matrix[k0][l], matrix[l0][l]);
		if (matrix[l0][l0] == 0.) {
			delete_struct(ii); delete_struct(matrix);	return(0);
		}
////////////////////////////////
//...diagonal row normalization;
		f = 1./matrix[l0][l0]; matrix[l0][l0] = 1.;
		for (coef[l0] *= f, l = 0; l <= N; l++) matrix[l0][l] *= f;
////////////////////////////////
//...elimination all other rows;
		for (k = 0; k <= N; k++)
			if ( k != l0) {
				f = matrix[k][l0]; matrix[k][l0] = 0.;
				for (coef[k] -= coef[l0]*f, l = 0; l <= N; l++) matrix[k][l] -= matrix[l0][l]*f;
			}
	}
	delete_struct(ii); delete_struct(matrix);
	return(1);
}

////////////////////////////////////////////////////////////
//...аппроксимации радиуса агломерации разложением Лагранжа;
double strata_aglom(double t, int N, double * coef)
{
	double sum = 0;
	for (int i = N; i > 0; i--) sum = coef[i]+t*sum; 
	return 1./sqrt(sum = coef[0]+sqr(t)*sum);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция аппроксимации межфазного слоя наименьшими квадратами ll = 1./sqrt(ll0**(-2)+B*(t**2-t0**2)), data[2*N+2];
double strata_approx(int N, double * data, double ll0)
{
	double sum = 0., fnorm = 0.;

///////////////////////////////////////////////
//...определяем скалярное произведение и норму;
	for (int i = 0; i <= N; i++) {
		sum += (sqr(data[2*i])-sqr(data[0]))*(1./sqr(data[2*i+1])-1./sqr(ll0));
		fnorm += sqr(sqr(data[2*i])-sqr(data[0]));
	}
	return(sum = sum/fnorm);
}

////////////////////////////////////////////////////////////////////////
//...аппроксимации радиуса агломерации разложением наименьших квадратов;
double strata_aglom(double t, double t0, double B, double ll0)
{
	return 1./sqrt(1./sqr(ll0)+B*(sqr(t)-sqr(t0)));
}

///////////////////////////////////////////////////
//...аппроксимации радиуса агломерации экспонентой;
double aglom_approx(double t, double alpha, double r_max, double r_min, double t_min)
{
	return max(r_min, r_max-(r_max-r_min)*exp(-alpha*(t-t_min)));
}

///////////////////////////////////////////////////
//...аппроксимации радиуса агломерации экспонентой;
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double t_min, double E1, double E2, double E3, double t1, double t2, double t3)
{
	double a1 = -log((E1_max-E1)/(E1_max-E_min)), b1 = -log((E1_max-E2)/(E1_max-E_min)), c1 = -log((E2_max-E3)/(E2_max-E_min)),
			 alpha = (b1*sqr(t1)*sqr(t1)-a1*sqr(t2)*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))),
			 alpha2 = -(b1*sqr(t1)-a1*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))), beta = c1/sqr(t3);
	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-sqr(t/t_min-1.)*(alpha+alpha2*sqr(t/t_min-1.))) : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
}
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double t_min, double E1, double E2, double E3, double E4, double t1, double t2, double t3, double t4)
{
	double a1 = -log((E1_max-E1)/(E1_max-E_min)), b1 = -log((E1_max-E2)/(E1_max-E_min)), 
			 c1 = -log((E2_max-E3)/(E2_max-E_min)), d1 = -log((E2_max-E4)/(E2_max-E_min)),
			 alpha = (b1*sqr(t1)*sqr(t1)-a1*sqr(t2)*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))),
			 alpha2 = -(b1*sqr(t1)-a1*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))), 
			 beta = (d1*sqr(t3)*sqr(t3)-c1*sqr(t4)*sqr(t4))/(sqr(t3)*sqr(t4)*(sqr(t3)-sqr(t4))),
			 beta2 = -(d1*sqr(t3)-c1*sqr(t4))/(sqr(t3)*sqr(t4)*(sqr(t3)-sqr(t4)));
	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-sqr(t/t_min-1.)*(alpha+alpha2*sqr(t/t_min-1.))) : 
							  E2_max-(E2_max-E_min)*exp(-sqr(t/t_min-1.)*(beta+beta2*sqr(t/t_min-1.))));
}

//////////////////////////////////////////////////////////////
//...аппроксимации радиуса агломерации экспонентой и степенью;
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double alpha, double beta, double t_min)
{
	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-alpha*sqr(t/t_min-1.)) : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
	//return (t < t_min ? (61.138+27.877*(t/t_min-1.))*sqr(t/t_min-1.)+60. : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
}
#endif
