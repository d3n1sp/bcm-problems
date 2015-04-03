#include "stdafx.h"
#include "cmapi.h"

///////////////////////////////
//...setting of the parameters;
void CMapi3DSpheroid::set_shape(double R0, double eps, double p1, double p2, double p3, double p4)
{
    CMapi3DSpheroid::R0 = R0;
    R0_inv     = R0 > EE ? 1./R0 : 1.;
    int   k    = -1;
    if (++k < size_of_param()) param[k] = (Param)eps; else return;
    if (++k < size_of_param()) param[k] = (Param)p1; else return;
}

void CMapi3DSpheroidFull::set_shape(double R0, double eps, double p1, double p2, double p3, double p4)
{
    CMapi3DSpheroidFull::R0 = R0;
    R0_inv     = R0 > EE ? 1./R0 : 1.;
    int   k    = -1;
    if (++k < size_of_param()) param[k] = (Param)eps; else return;
    if (++k < size_of_param()) param[k] = (Param)p1; else return;
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CMapi3DSpheroid::parametrization(double * P, int m)
{
	if (! p) {
		p = (double *)new_struct((NN_dop = freedom(N = 1))*sizeof(double));
		px = (double *)new_struct(4*sizeof(double));
		py = (double *)new_struct(4*sizeof(double));
		pz = (double *)new_struct(4*sizeof(double));
	}
	if (P) { 
		double  X  = P[0]*R0_inv,
				  Y  = P[1]*R0_inv,
				  Z  = P[2]*R0_inv;
		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

///////////////////////////////////////////////
//...параметры сфероидальной системы координат;
		double eps = param[0], R0 = param[1], f0 = R0, z0 = 1., ch_alpha = 1., sh_alpha = 0., 
			h0, h1, h2, th, cth, d0, d1, d2, rr, f, ff, f1, f2, Co = 1., Si = 0., FN2;
		if (0. < eps && eps < 1.) {
			ch_alpha = 1./sqrt(1.-eps*eps); 
			sh_alpha = eps*ch_alpha;
			FN2 = 0.; 
			h0 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.));
			h1 =  1./ch_alpha+h0; 
			h2 = 1./(eps*sh_alpha)+h0;
			z0 = ch_alpha;
			f0 = R0/z0;
			d0 =  1./sh_alpha;
			d1 =  1./(sqr(ch_alpha)*sh_alpha);
			d2 = -2./(sqr(sh_alpha)*sh_alpha);
			th = eps;
			cth = 1./th;
		}

/////////////////////////////////////
//...фактические координаты сфероида;
		rr = sqrt(sqr(P[0])+sqr(P[1]));
		if (rr > EE) {
			Co = P[0]/rr; Si = P[1]/rr; 
		}
		f = sqr(rr)+sqr(f0);	f1 = sqr(P[2]);
		if ((ff = sqr(P[2])+sqr(rr/eps)) >= sqr(R0)) {
			f2 = (f1+2.*f-4.*sqr(f0))/sqr(f);
			if (f1 < 1e-5) 
				ch_alpha = (1.-f2*f*.5*(1.-f1*f2*.25))*.5;
			else
				ch_alpha = (1.-f/f1*(sqrt(1.+f1*f2)-1.))*.5;
			ch_alpha = 1./sqrt(ch_alpha); 
			sh_alpha = sqrt(sqr(ch_alpha)-1.);
			FN2 = 1./(sqr(f0*ch_alpha)-sqr(P[2]/ch_alpha)); 
			h0 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.));
			h1 = 1./ch_alpha+h0; 
			h2 = ch_alpha/(sqr(ch_alpha)-1.)+h0;
			d0 =  1./sh_alpha;
			d1 =  1./(sqr(ch_alpha)*sh_alpha);
			d2 = -2./(sqr(sh_alpha)*sh_alpha);
			th = sh_alpha/ch_alpha;
			cth = 1./th;
		}

///////////////////////////////////////////////////////////
//...первые четыре функции сфероидальной системы координат;
		p[0] = h0;
		p[1] = Z*h1; 
		p[2] = Y*h2; 
		p[3] = X*h2; 

/////////////////////////////////
//...сразу вычисл€ем производные;
		px[0] = d0*cth*P[0]*FN2;
		py[0] = d0*cth*P[1]*FN2;
		pz[0] = d0* th*P[2]*FN2;
		px[1] = d1*cth*P[0]*FN2*Z;
		py[1] = d1*cth*P[1]*FN2*Z;
		pz[1] = d1* th*P[2]*FN2*Z+h1*R0_inv;
		px[2] = d2*cth*P[0]*FN2*Y;
		py[2] = d2*cth*P[1]*FN2*Y+h2*R0_inv;
		pz[2] = d2* th*P[2]*FN2*Y;
		px[3] = d2*cth*P[0]*FN2*X+h2*R0_inv;
		py[3] = d2*cth*P[1]*FN2*X;
		pz[3] = d2* th*P[2]*FN2*X;
	}
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CMapi3DSpheroidFull::parametrization(double * P, int m_dop)
{
	if (! p) {
		p   = (double *)new_struct((NN_dop = 16/*min(freedom(N), 16)*/)*sizeof(double));
		px  = (double *)new_struct (NN_dop*sizeof(double));
		py  = (double *)new_struct (NN_dop*sizeof(double));
		pz  = (double *)new_struct (NN_dop*sizeof(double));
		pxx = (double *)new_struct (NN_dop*sizeof(double));
		pxy = (double *)new_struct (NN_dop*sizeof(double));
		pyy = (double *)new_struct (NN_dop*sizeof(double));
		pxz = (double *)new_struct (NN_dop*sizeof(double));
		pyz = (double *)new_struct (NN_dop*sizeof(double));
		pzz = (double *)new_struct (NN_dop*sizeof(double));
	}
	if (P && p) {
		double   Z = P[2]/**R0_inv*/, rr, h0, f0 = param[1], f2, ff, f, f1 = Z*Z, f3, f4, f5;
		complex  z = comp(P[0], P[1])/**R0_inv*/, zm, F1, F2;
		int      i, m; 
		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

///////////////////////////////////////////////
//...параметры сфероидальной системы координат;
		double delt = param[0], R0 = param[1], ch_alpha = 1., sh_alpha = 0., th, cth, cos_beta, Co = 1., Si = 0., FN2 = 0.,
			alpha_x, alpha_y, alpha_z, alpha_xx, alpha_xy, alpha_yy, alpha_xz, alpha_yz, alpha_zz, deriv_h, dd_h, hh, 
			dQ, dP, dQ_prev, dP_prev, ddP, ddQ;
		if (0. < delt && delt < 1.) f0 = R0*sqrt(1.-sqr(delt));
		if (0. < delt && delt < 1.) {
			ch_alpha = 1./sqrt(1.-sqr(delt)); 
			sh_alpha = delt*ch_alpha;
			h0 = log((ch_alpha-1.)/(ch_alpha+1.));
			f0 = R0/ch_alpha; /*f0 = R0;*/
			th = delt;
			cth = 1./th;
		}

/////////////////////////////////////
//...фактические координаты сфероида;
		if ((rr = sqrt(sqr(P[0])+sqr(P[1]))) > EE) {
			Co = P[0]/rr; Si = P[1]/rr; 
		}
		f = sqr(rr)+sqr(f0);
		if ((ff = sqr(P[2])+sqr(rr/delt)) >= sqr(R0)) {
			f2 = (f1+2.*f-4.*sqr(f0))/sqr(f);
			if (f1 < 1e-5) 
				ch_alpha = (1.-f2*f*.5*(1.-f1*f2*.25))*.5;
			else
				ch_alpha = (1.-f/f1*(sqrt(1.+f1*f2)-1.))*.5;
			ch_alpha = 1./sqrt(ch_alpha); 
			sh_alpha = sqrt(sqr(ch_alpha)-1.);
			FN2 = f0/(sqr(f0*ch_alpha)-sqr(P[2]/ch_alpha)); 
			h0  = log((ch_alpha-1.)/(ch_alpha+1.));
			th  = sh_alpha/ch_alpha;
			cth = 1./th;
		}
		cos_beta = Z/(f0*ch_alpha); z /= f0*sh_alpha;
///////////////////////////////
//...тестова€ проверка функций;
#ifdef ___CONSTRUCTION___
//////////////////////////////////////////
//...вычисление угловой части мультиполей;
		for (F1 = (zm = comp(1.)), F2 = comp(0.), i = 0; i < N+N_dop; i++) {
			  F2 = ((2.*i+1.)*cos_beta*(p[i*i] = real(F1))-i*real(F2))/(i+1.); 
			  cSWAP(F1, F2);
		}
		for (p[i*i] = real(F1), m = 1; m <= N+N_dop; m++) {
			 for (F1 = (zm *= z), F2 = comp(0.), i = m; i < N+N_dop; i++) {
					p[i*i+m*2-1] = imag(F1);
					p[i*i+m*2]   = real(F1);
					F2 = ((2.*i+1.)*cos_beta*F1-(i-m)*F2)/(i+m+1.); 
					cSWAP(F1, F2);
			 }
			 p[i*i+m*2-1] = imag(F1);
			 p[i*i+m*2]   = real(F1);
		}
/////////////////////////////////////////////
//...вычисление радиальной части мультиполей;
		for (f1 = f3 = 1., f2 = 0., i = 0; i < N+N_dop; i++) {
			p[i*i] *= f1;
			f2 = ((2.*i+1.)*ch_alpha*f1-i*f2)/(i+1.); swap(f1, f2);
		}
		for (p[i*i] *= f1, m = 1; m <= N+N_dop; m++) {
			for (f1 = (f3 *= sh_alpha), f2 = 0., i = m; i < N+N_dop; i++) {
				p[i*i+m*2-1] *= f1;
				p[i*i+m*2]   *= f1;
				f2 = ((2.*i+1.)*ch_alpha*f1-(i-m)*f2)/(i+m+1.);	swap(f1, f2);
			}
			p[i*i+m*2-1] *= f1;
			p[i*i+m*2]   *= f1;
		}
//////////////////////
//...проверка функций;
		p[0] -= 1.;
		p[1] -= P[2]/f0;
		p[2] -= P[1]/f0;
		p[3] -= P[0]/f0;
		p[4] -= (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/2.)*1.5/sqr(f0)-1./2.;
		p[5] -= P[2]*P[1]/sqr(f0);
		p[6] -= P[2]*P[0]/sqr(f0);
		p[7] -= 2.*P[0]*P[1]/sqr(f0);
		p[8] -= (sqr(P[0])-sqr(P[1]))/sqr(f0);
		p[9] -= P[2]*(sqr(P[2])-3.*(sqr(P[0])+sqr(P[1]))/2.)*2.5/(f0*sqr(f0))-1.5*P[2]/f0;
		p[10] -= (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/4.)*P[1]*1.25/(f0*sqr(f0))-0.25*P[1]/f0;
		p[11] -= (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/4.)*P[0]*1.25/(f0*sqr(f0))-0.25*P[0]/f0;
		p[12] -= P[2]*(2.*P[0]*P[1])/(f0*sqr(f0));
		p[13] -= P[2]*(sqr(P[0])-sqr(P[1]))/(f0*sqr(f0));
		p[14] -= (3.*sqr(P[0])*P[1]-sqr(P[1])*P[1])/(f0*sqr(f0));
		p[15] -= (sqr(P[0])*P[0]-3.*P[0]*sqr(P[1]))/(f0*sqr(f0));
		f3 = f3;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление радиальной части мультиполей второго рода и производные по z (только до третьего пор€дка);
		f2 = h0; f1 = ch_alpha*f2+2.; f3 = 2.;
		pz[0] = f3*cos_beta*p[0]*FN2; f3 = p[0]; p[0] *= f2; 
		if (N+N_dop > 0) {
			for (i = 1; i < N+N_dop; i++) {
				pz[i*i] = i*(ch_alpha*f1*f3-cos_beta*f2*p[i*i])*FN2; f3 = p[i*i];	p[i*i] *= f1;
				f2 = ((2.*i+1.)*ch_alpha*f1-i*f2)/(i+1.); swap(f1, f2);
			}
			pz[i*i] = i*(ch_alpha*f1*f3-cos_beta*f2*p[i*i])*FN2; p[i*i] *= f1;
		}
		for (m = 1; m <= N+N_dop; m++) {
			if (m == 1) { f2 = sh_alpha*h0+2.*cth; f1 = ch_alpha*f2-4./(3.*sh_alpha); f3 = ff = -4./sh_alpha;} else 
			if (m == 2) { f2 = sqr(sh_alpha)*h0+2.*ch_alpha-4.*ch_alpha/(3.*sqr(sh_alpha)); f1 = ch_alpha*f2+16./(15.*sqr(sh_alpha)); f3 = ff = 16./(3.*sqr(sh_alpha));} else
			if (m == 3) { f2 = sh_alpha*(sqr(sh_alpha)*h0+2.*ch_alpha)+.8*cth*(4./3.*sqr(cth)-3.); f1 = 0.; f3 = ff = -32./(5.*sh_alpha*sqr(sh_alpha));}
			pz[m*m+m*2-1] = f3*cos_beta*p[m*m+m*2-1]*FN2; f3 = p[m*m+m*2-1]; p[m*m+m*2-1] *= f2;
			pz[m*m+m*2]   = ff*cos_beta*p[m*m+m*2]*FN2;   ff = p[m*m+m*2];   p[m*m+m*2]   *= f2;
			if (N+N_dop > m) {
				for (i = m+1; i < N+N_dop; i++) {
					pz[i*i+m*2-1] =(i-m)*(ch_alpha*f1*f3-cos_beta*f2*p[i*i+m*2-1])*FN2; f3 = p[i*i+m*2-1]; p[i*i+m*2-1] *= f1;
					pz[i*i+m*2]   =(i-m)*(ch_alpha*f1*ff-cos_beta*f2*p[i*i+m*2])*FN2;   ff = p[i*i+m*2];   p[i*i+m*2]   *= f1;
					f2 = ((2.*i+1.)*ch_alpha*f1-(i-m)*f2)/(i+m+1.);	swap(f1, f2);
				}
				pz[i*i+m*2-1] =(i-m)*(ch_alpha*f1*f3-cos_beta*f2*p[i*i+m*2-1])*FN2; p[i*i+m*2-1] *= f1;
				pz[i*i+m*2]   =(i-m)*(ch_alpha*f1*ff-cos_beta*f2*p[i*i+m*2])*FN2;   p[i*i+m*2]   *= f1;
			}
		}
#endif
////////////////////////////////////////////
//...регул€рна€ часть функций и производные;
		p[0] = 1.;
		p[1] = P[2]/f0;
		p[2] = P[1]/f0;
		p[3] = P[0]/f0;
		p[4] = (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/2.)*1.5/sqr(f0)-1./2.;
		p[5] = P[2]*P[1]/sqr(f0);
		p[6] = P[2]*P[0]/sqr(f0);
		p[7] = 2.*P[0]*P[1]/sqr(f0);
		p[8] = (sqr(P[0])-sqr(P[1]))/sqr(f0);
		p[9] = P[2]*(sqr(P[2])-3.*(sqr(P[0])+sqr(P[1]))/2.)*2.5/(f0*sqr(f0))-1.5*P[2]/f0;
		p[10] = (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/4.)*P[1]*1.25/(f0*sqr(f0))-0.25*P[1]/f0;
		p[11] = (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/4.)*P[0]*1.25/(f0*sqr(f0))-0.25*P[0]/f0;
		p[12] = P[2]*(2.*P[0]*P[1])/(f0*sqr(f0));
		p[13] = P[2]*(sqr(P[0])-sqr(P[1]))/(f0*sqr(f0));
		p[14] = (3.*sqr(P[0])*P[1]-sqr(P[1])*P[1])/(f0*sqr(f0));
		p[15] = (sqr(P[0])*P[0]-3.*P[0]*sqr(P[1]))/(f0*sqr(f0));
//...по z;
		pz[0] = 0.;
		pz[1] = 1./f0;
		pz[2] = 0.;
		pz[3] = 0.;
		pz[4] = P[2]*3./sqr(f0);
		pz[5] = P[1]/sqr(f0);
		pz[6] = P[0]/sqr(f0);
		pz[7] = 0.;
		pz[8] = 0.;
		pz[9] = (sqr(P[2])-(sqr(P[0])+sqr(P[1]))/2.)*7.5/(f0*sqr(f0))-1.5/f0;
		pz[10] = P[2]*P[1]*2.5/(f0*sqr(f0));
		pz[11] = P[2]*P[0]*2.5/(f0*sqr(f0));
		pz[12] = P[0]*P[1]*2.0/(f0*sqr(f0));
		pz[13] = (sqr(P[0])-sqr(P[1]))/(f0*sqr(f0));
		pz[14] = 0.;
		pz[15] = 0.;
//...по y;
		py[0] = 0.;
		py[1] = 0.;
		py[2] = 1./f0;
		py[3] = 0.;
		py[4] = -P[1]*1.5/sqr(f0);
		py[5] = P[2]/sqr(f0);
		py[6] = 0.;
		py[7] =  2.*P[0]/sqr(f0);
		py[8] = -2.*P[1]/sqr(f0);
		py[9] = -3.*P[1]*P[2]*2.5/(f0*sqr(f0));
		py[10] = (sqr(P[2])-(sqr(P[0])+3.*sqr(P[1]))/4.)*1.25/(f0*sqr(f0))-0.25/f0;
		py[11] = -P[1]*P[0]*0.625/(f0*sqr(f0));
		py[12] =  2.*P[0]*P[2]/(f0*sqr(f0));
		py[13] = -2.*P[1]*P[2]/(f0*sqr(f0));
		py[14] =  3.*(sqr(P[0])-sqr(P[1]))/(f0*sqr(f0));
		py[15] = -6.*P[0]*P[1]/(f0*sqr(f0));
//...по x;
		px[0] = 0.;
		px[1] = 0.;
		px[2] = 0.;
		px[3] = 1./f0;
		px[4] = -P[0]*1.5/sqr(f0);
		px[5] = 0.;
		px[6] = P[2]/sqr(f0);
		px[7] =  2.*P[1]/sqr(f0);
		px[8] =  2.*P[0]/sqr(f0);
		px[9] = -3.*P[0]*P[2]*2.5/(f0*sqr(f0));
		px[10] = -P[0]*P[1]*0.625/(f0*sqr(f0));
		px[11] = (sqr(P[2])-(3.*sqr(P[0])+sqr(P[1]))/4.)*1.25/(f0*sqr(f0))-0.25/f0;
		px[12] = 2.*P[1]*P[2]/(f0*sqr(f0));
		px[13] = 2.*P[0]*P[2]/(f0*sqr(f0));
		px[14] = 6.*P[0]*P[1]/(f0*sqr(f0));
		px[15] = 3.*(sqr(P[0])-sqr(P[1]))/(f0*sqr(f0));
//...по xx;
		pxx[0] = 0.;
		pxx[1] = 0.;
		pxx[2] = 0.;
		pxx[3] = 0.;
		pxx[4] = -1.5/sqr(f0);
		pxx[5] = 0.;
		pxx[6] = 0.;
		pxx[7] =  0.;
		pxx[8] =  2./sqr(f0);
		pxx[9] = -3.*P[2]*2.5/(f0*sqr(f0));
		pxx[10] = -P[1]*0.625/(f0*sqr(f0));
		pxx[11] = -P[0]*1.875/(f0*sqr(f0));
		pxx[12] = 0.;
		pxx[13] = 2.*P[2]/(f0*sqr(f0));
		pxx[14] = 6.*P[1]/(f0*sqr(f0));
		pxx[15] = 6.*P[0]/(f0*sqr(f0));
//...по xy;
		pxy[0] = 0.;
		pxy[1] = 0.;
		pxy[2] = 0.;
		pxy[3] = 0.;
		pxy[4] = 0.;
		pxy[5] = 0.;
		pxy[6] = 0.;
		pxy[7] =  2./sqr(f0);
		pxy[8] =  0.;
		pxy[9] = 0.;
		pxy[10] = -P[0]*0.625/(f0*sqr(f0));
		pxy[11] = -P[1]*0.625/(f0*sqr(f0));
		pxy[12] = 2.*P[2]/(f0*sqr(f0));
		pxy[13] = 0.;
		pxy[14] =  6.*P[0]/(f0*sqr(f0));
		pxy[15] = -6.*P[1]/(f0*sqr(f0));
//...по y;
		pyy[0] = 0.;
		pyy[1] = 0.;
		pyy[2] = 0.;
		pyy[3] = 0.;
		pyy[4] = -1.5/sqr(f0);
		pyy[5] = 0.;
		pyy[6] = 0.;
		pyy[7] =  0.;
		pyy[8] = -2./sqr(f0);
		pyy[9] = -3.*P[2]*2.5/(f0*sqr(f0));
		pyy[10] = -P[1]*1.875/(f0*sqr(f0));
		pyy[11] = -P[0]*0.625/(f0*sqr(f0));
		pyy[12] =  0.;
		pyy[13] = -2.*P[2]/(f0*sqr(f0));
		pyy[14] = -6.*P[1]/(f0*sqr(f0));
		pyy[15] = -6.*P[0]/(f0*sqr(f0));
//...по xz;
		pxz[0] = 0.;
		pxz[1] = 0.;
		pxz[2] = 0.;
		pxz[3] = 0.;
		pxz[4] = 0.;
		pxz[5] = 0.;
		pxz[6] = 1./sqr(f0);
		pxz[7] =  0.;
		pxz[8] =  0.;
		pxz[9] = -P[0]*7.5/(f0*sqr(f0));
		pxz[10] = 0.;
		pxz[11] = P[2]*2.5/(f0*sqr(f0));
		pxz[12] = 2.*P[1]/(f0*sqr(f0));
		pxz[13] = 2.*P[0]/(f0*sqr(f0));
		pxz[14] = 0.;
		pxz[15] = 0.;
//...по yz;
		pyz[0] = 0.;
		pyz[1] = 0.;
		pyz[2] = 0.;
		pyz[3] = 0.;
		pyz[4] = 0.;
		pyz[5] = 1./sqr(f0);
		pyz[6] = 0.;
		pyz[7] = 0.;
		pyz[8] = 0.;
		pyz[9] = -P[1]*7.5/(f0*sqr(f0));
		pyz[10] = P[2]*2.5/(f0*sqr(f0));
		pyz[11] = 0.;
		pyz[12] =  2.*P[0]/(f0*sqr(f0));
		pyz[13] = -2.*P[1]/(f0*sqr(f0));
		pyz[14] = 0.;
		pyz[15] = 0.;
//...по zz;
		pzz[0] = 0.;
		pzz[1] = 0.;
		pzz[2] = 0.;
		pzz[3] = 0.;
		pzz[4] = 3./sqr(f0);
		pzz[5] = 0.;
		pzz[6] = 0.;
		pzz[7] = 0.;
		pzz[8] = 0.;
		pzz[9] = P[2]*15./(f0*sqr(f0));
		pzz[10] = P[1]*2.5/(f0*sqr(f0));
		pzz[11] = P[0]*2.5/(f0*sqr(f0));
		pzz[12] = 0.;
		pzz[13] = 0.;
		pzz[14] = 0.;
		pzz[15] = 0.;

//////////////////////////////////////////////////////////////////////
//...вычисление прозводных через дифференцирование радиальной функции;
		f2 = h0; f1 = ch_alpha*f2+2.; f3 = 2.;	f5 = ff = 1.; f4 = ch_alpha*f5;
		alpha_x  = ch_alpha*FN2*real(z);
		alpha_y  = ch_alpha*FN2*imag(z);
		alpha_z  = sh_alpha*FN2*cos_beta;
		alpha_xx = cth*FN2*(1./f0-sqr(real(z))*FN2*(1.+2.*sqr(sh_alpha)*(sqr(ch_alpha)+sqr(cos_beta))/(sqr(ch_alpha)-sqr(cos_beta))));
		alpha_xy = -cth*FN2*real(z)*imag(z)*FN2*(1.+2.*sqr(sh_alpha)*(sqr(ch_alpha)+sqr(cos_beta))/(sqr(ch_alpha)-sqr(cos_beta)));
		alpha_yy = cth*FN2*(1./f0-sqr(imag(z))*FN2*(1.+2.*sqr(sh_alpha)*(sqr(ch_alpha)+sqr(cos_beta))/(sqr(ch_alpha)-sqr(cos_beta))));
		alpha_xz = -FN2*real(z)*cos_beta*FN2*(1.+2.*sqr(ch_alpha)*(sqr(sh_alpha)+sqr(cos_beta)-1.)/(sqr(ch_alpha)-sqr(cos_beta)));
		alpha_yz = -FN2*imag(z)*cos_beta*FN2*(1.+2.*sqr(ch_alpha)*(sqr(sh_alpha)+sqr(cos_beta)-1.)/(sqr(ch_alpha)-sqr(cos_beta)));
		alpha_zz = th*FN2*(1./f0+sqr(cos_beta)*FN2*(1.-2.*sqr(ch_alpha)*(sqr(sh_alpha)+sqr(cos_beta)-1.)/(sqr(ch_alpha)-sqr(cos_beta))));
		deriv_h  = f3/(f5*sh_alpha); hh = f2/f5;
		dd_h  = -f3/f5*ch_alpha/sqr(sh_alpha); 
		pxx[0] = hh*pxx[0]+2.*deriv_h*alpha_x*px[0]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[0];
		pxy[0] = hh*pxy[0]+deriv_h*(alpha_x*py[0]+alpha_y*px[0])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[0];
		pyy[0] = hh*pyy[0]+2.*deriv_h*alpha_y*py[0]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[0];
		pxz[0] = hh*pxz[0]+deriv_h*(alpha_x*pz[0]+alpha_z*px[0])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[0];
		pyz[0] = hh*pyz[0]+deriv_h*(alpha_y*pz[0]+alpha_z*py[0])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[0];
		pzz[0] = hh*pzz[0]+2.*deriv_h*alpha_z*pz[0]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[0];
		pz[0] = (f2*pz[0]+f3*p[0]*FN2*cos_beta)/f5; 
		py[0] = (f2*py[0]+f3*p[0]*cth*FN2*imag(z))/f5; 
		px[0] = (f2*px[0]+f3*p[0]*cth*FN2*real(z))/f5; 
		p[0] *= hh; dQ_prev = f3/sh_alpha; dP_prev = 0.;
		if (N > 0) {
			for (i = 1; i < N; i++) {
				deriv_h  = i*((hh = f1/f4)*f5-f2)/(f4*sh_alpha);
				dP = i*(ch_alpha*f4-f5)/sh_alpha;
				dQ = i*(ch_alpha*f1-f2)/sh_alpha;
				ddP = ((i-1.)*ch_alpha*dP-i*dP_prev)/sh_alpha+i*f4;
				ddQ = ((i-1.)*ch_alpha*dQ-i*dQ_prev)/sh_alpha+i*f1;
				dd_h = (ddQ-hh*ddP)/f4-2.*dP*(dQ-hh*dP)/sqr(f4); dQ_prev = dQ; dP_prev= dP;
				pxx[i*i] = hh*pxx[i*i]+2.*deriv_h*alpha_x*px[i*i]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i];
				pxy[i*i] = hh*pxy[i*i]+deriv_h*(alpha_x*py[i*i]+alpha_y*px[i*i])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i];
				pyy[i*i] = hh*pyy[i*i]+2.*deriv_h*alpha_y*py[i*i]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i];
				pxz[i*i] = hh*pxz[i*i]+deriv_h*(alpha_x*pz[i*i]+alpha_z*px[i*i])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i];
				pyz[i*i] = hh*pyz[i*i]+deriv_h*(alpha_y*pz[i*i]+alpha_z*py[i*i])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i];
				pzz[i*i] = hh*pzz[i*i]+2.*deriv_h*alpha_z*pz[i*i]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i];
				pz[i*i] = (f1*pz[i*i]+i*(hh*f5-f2)*p[i*i]*FN2*cos_beta)/f4; 
				py[i*i] = (f1*py[i*i]+i*(hh*f5-f2)*p[i*i]*cth*FN2*imag(z))/f4; 
				px[i*i] = (f1*px[i*i]+i*(hh*f5-f2)*p[i*i]*cth*FN2*real(z))/f4; 
				p[i*i] *= hh;
				f2 = ((2.*i+1.)*ch_alpha*f1-i*f2)/(i+1.); swap(f1, f2);
				f5 = ((2.*i+1.)*ch_alpha*f4-i*f5)/(i+1.); swap(f4, f5);
			}
			deriv_h  = i*((hh = f1/f4)*f5-f2)/(f4*sh_alpha);
			dP = i*(ch_alpha*f4-f5)/sh_alpha;
			dQ = i*(ch_alpha*f1-f2)/sh_alpha;
			ddP = ((i-1.)*ch_alpha*dP-i*dP_prev)/sh_alpha+i*f4;
			ddQ = ((i-1.)*ch_alpha*dQ-i*dQ_prev)/sh_alpha+i*f1;
			dd_h = (ddQ-hh*ddP)/f4-2.*dP*(dQ-hh*dP)/sqr(f4);
			pxx[i*i] = hh*pxx[i*i]+2.*deriv_h*alpha_x*px[i*i]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i];
			pxy[i*i] = hh*pxy[i*i]+deriv_h*(alpha_x*py[i*i]+alpha_y*px[i*i])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i];
			pyy[i*i] = hh*pyy[i*i]+2.*deriv_h*alpha_y*py[i*i]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i];
			pxz[i*i] = hh*pxz[i*i]+deriv_h*(alpha_x*pz[i*i]+alpha_z*px[i*i])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i];
			pyz[i*i] = hh*pyz[i*i]+deriv_h*(alpha_y*pz[i*i]+alpha_z*py[i*i])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i];
			pzz[i*i] = hh*pzz[i*i]+2.*deriv_h*alpha_z*pz[i*i]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i];
			pz[i*i] = (f1*pz[i*i]+i*(f1/f4*f5-f2)*p[i*i]*FN2*cos_beta)/f4; 
			py[i*i] = (f1*py[i*i]+i*(f1/f4*f5-f2)*p[i*i]*cth*FN2*imag(z))/f4; 
			px[i*i] = (f1*px[i*i]+i*(f1/f4*f5-f2)*p[i*i]*cth*FN2*real(z))/f4; 
			p[i*i] *= hh;
		}
		for (m = 1; m <= N; m++) {
			f5 = (ff *= sh_alpha); f4 = ch_alpha*f5; 
			if (m == 1) { f2 = sh_alpha*h0+2.*cth; f1 = ch_alpha*f2-4./(3.*sh_alpha); f3 = -4./sh_alpha;} else 
			if (m == 2) { f2 = sqr(sh_alpha)*h0+2.*ch_alpha-4.*ch_alpha/(3.*sqr(sh_alpha)); f1 = ch_alpha*f2+16./(15.*sqr(sh_alpha)); f3 = 16./(3.*sqr(sh_alpha));} else
			if (m == 3) { f2 = sh_alpha*(sqr(sh_alpha)*h0+2.*ch_alpha)+.8*cth*(4./3.*sqr(cth)-3.); f1 = 0.; f3 = -32./(5.*sh_alpha*sqr(sh_alpha));}
			deriv_h  = f3/(f5*sh_alpha); hh = f2/f5;
			dd_h  = -(2.*m+1.)*f3/f5*ch_alpha/sqr(sh_alpha); 
			pxx[m*m+m*2-1] = hh*pxx[m*m+m*2-1]+2.*deriv_h*alpha_x*px[m*m+m*2-1]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[m*m+m*2-1];
			pxy[m*m+m*2-1] = hh*pxy[m*m+m*2-1]+deriv_h*(alpha_x*py[m*m+m*2-1]+alpha_y*px[m*m+m*2-1])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[m*m+m*2-1];
			pyy[m*m+m*2-1] = hh*pyy[m*m+m*2-1]+2.*deriv_h*alpha_y*py[m*m+m*2-1]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[m*m+m*2-1];
			pxz[m*m+m*2-1] = hh*pxz[m*m+m*2-1]+deriv_h*(alpha_x*pz[m*m+m*2-1]+alpha_z*px[m*m+m*2-1])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[m*m+m*2-1];
			pyz[m*m+m*2-1] = hh*pyz[m*m+m*2-1]+deriv_h*(alpha_y*pz[m*m+m*2-1]+alpha_z*py[m*m+m*2-1])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[m*m+m*2-1];
			pzz[m*m+m*2-1] = hh*pzz[m*m+m*2-1]+2.*deriv_h*alpha_z*pz[m*m+m*2-1]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[m*m+m*2-1];
			pxx[m*m+m*2] = hh*pxx[m*m+m*2]+2.*deriv_h*alpha_x*px[m*m+m*2]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[m*m+m*2];
			pxy[m*m+m*2] = hh*pxy[m*m+m*2]+deriv_h*(alpha_x*py[m*m+m*2]+alpha_y*px[m*m+m*2])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[m*m+m*2];
			pyy[m*m+m*2] = hh*pyy[m*m+m*2]+2.*deriv_h*alpha_y*py[m*m+m*2]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[m*m+m*2];
			pxz[m*m+m*2] = hh*pxz[m*m+m*2]+deriv_h*(alpha_x*pz[m*m+m*2]+alpha_z*px[m*m+m*2])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[m*m+m*2];
			pyz[m*m+m*2] = hh*pyz[m*m+m*2]+deriv_h*(alpha_y*pz[m*m+m*2]+alpha_z*py[m*m+m*2])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[m*m+m*2];
			pzz[m*m+m*2] = hh*pzz[m*m+m*2]+2.*deriv_h*alpha_z*pz[m*m+m*2]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[m*m+m*2];
			pz[m*m+m*2-1] = (f2*pz[m*m+m*2-1]+f3*p[m*m+m*2-1]*FN2*cos_beta)/f5; 
			py[m*m+m*2-1] = (f2*py[m*m+m*2-1]+f3*p[m*m+m*2-1]*cth*FN2*imag(z))/f5; 
			px[m*m+m*2-1] = (f2*px[m*m+m*2-1]+f3*p[m*m+m*2-1]*cth*FN2*real(z))/f5; 
			pz[m*m+m*2] = (f2*pz[m*m+m*2]+f3*p[m*m+m*2]*FN2*cos_beta)/f5; 
			py[m*m+m*2] = (f2*py[m*m+m*2]+f3*p[m*m+m*2]*cth*FN2*imag(z))/f5; 
			px[m*m+m*2] = (f2*px[m*m+m*2]+f3*p[m*m+m*2]*cth*FN2*real(z))/f5; 
			p[m*m+m*2-1] *= hh;
			p[m*m+m*2]   *= hh; dQ_prev = f3/sh_alpha+m*f2*cth; dP_prev = m*f5*cth;
			if (N > m) {
				for (i = m+1; i < N; i++) {
					deriv_h = (i-m)*((hh = f1/f4)*f5-f2)/(f4*sh_alpha);
					dP = (i*ch_alpha*f4-(i-m)*f5)/sh_alpha;
					dQ = (i*ch_alpha*f1-(i-m)*f2)/sh_alpha;
					ddP = ((i-1.)*ch_alpha*dP-(i-m)*dP_prev)/sh_alpha+i*f4;
					ddQ = ((i-1.)*ch_alpha*dQ-(i-m)*dQ_prev)/sh_alpha+i*f1;
					dd_h = (ddQ-hh*ddP)/f4-2.*dP*(dQ-hh*dP)/sqr(f4); dQ_prev = dQ; dP_prev= dP;
					pxx[i*i+m*2-1] = hh*pxx[i*i+m*2-1]+2.*deriv_h*alpha_x*px[i*i+m*2-1]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i+m*2-1];
					pxy[i*i+m*2-1] = hh*pxy[i*i+m*2-1]+deriv_h*(alpha_x*py[i*i+m*2-1]+alpha_y*px[i*i+m*2-1])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i+m*2-1];
					pyy[i*i+m*2-1] = hh*pyy[i*i+m*2-1]+2.*deriv_h*alpha_y*py[i*i+m*2-1]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i+m*2-1];
					pxz[i*i+m*2-1] = hh*pxz[i*i+m*2-1]+deriv_h*(alpha_x*pz[i*i+m*2-1]+alpha_z*px[i*i+m*2-1])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i+m*2-1];
					pyz[i*i+m*2-1] = hh*pyz[i*i+m*2-1]+deriv_h*(alpha_y*pz[i*i+m*2-1]+alpha_z*py[i*i+m*2-1])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i+m*2-1];
					pzz[i*i+m*2-1] = hh*pzz[i*i+m*2-1]+2.*deriv_h*alpha_z*pz[i*i+m*2-1]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i+m*2-1];
					pxx[i*i+m*2] = hh*pxx[i*i+m*2]+2.*deriv_h*alpha_x*px[i*i+m*2]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i+m*2];
					pxy[i*i+m*2] = hh*pxy[i*i+m*2]+deriv_h*(alpha_x*py[i*i+m*2]+alpha_y*px[i*i+m*2])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i+m*2];
					pyy[i*i+m*2] = hh*pyy[i*i+m*2]+2.*deriv_h*alpha_y*py[i*i+m*2]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i+m*2];
					pxz[i*i+m*2] = hh*pxz[i*i+m*2]+deriv_h*(alpha_x*pz[i*i+m*2]+alpha_z*px[i*i+m*2])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i+m*2];
					pyz[i*i+m*2] = hh*pyz[i*i+m*2]+deriv_h*(alpha_y*pz[i*i+m*2]+alpha_z*py[i*i+m*2])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i+m*2];
					pzz[i*i+m*2] = hh*pzz[i*i+m*2]+2.*deriv_h*alpha_z*pz[i*i+m*2]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i+m*2];
					pz[i*i+m*2-1] = (f1*pz[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*FN2*cos_beta)/f4; 
					py[i*i+m*2-1] = (f1*py[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*cth*FN2*imag(z))/f4; 
					px[i*i+m*2-1] = (f1*px[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*cth*FN2*real(z))/f4; 
					pz[i*i+m*2] = (f1*pz[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*FN2*cos_beta)/f4; 
					py[i*i+m*2] = (f1*py[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*cth*FN2*imag(z))/f4; 
					px[i*i+m*2] = (f1*px[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*cth*FN2*real(z))/f4; 
					p[i*i+m*2-1] *= hh;
					p[i*i+m*2]   *= hh;
					f2 = ((2.*i+1.)*ch_alpha*f1-(i-m)*f2)/(i+m+1.);	swap(f1, f2);
					f5 = ((2.*i+1.)*ch_alpha*f4-(i-m)*f5)/(i+m+1.);	swap(f4, f5);
				}
				deriv_h = (i-m)*((hh = f1/f4)*f5-f2)/(f4*sh_alpha);
				dP = (i*ch_alpha*f4-(i-m)*f5)/sh_alpha;
				dQ = (i*ch_alpha*f1-(i-m)*f2)/sh_alpha;
				ddP = ((i-1.)*ch_alpha*dP-(i-m)*dP_prev)/sh_alpha+i*f4;
				ddQ = ((i-1.)*ch_alpha*dQ-(i-m)*dQ_prev)/sh_alpha+i*f1;
				dd_h = (ddQ-hh*ddP)/f4-2.*dP*(dQ-hh*dP)/sqr(f4);
				pxx[i*i+m*2-1] = hh*pxx[i*i+m*2-1]+2.*deriv_h*alpha_x*px[i*i+m*2-1]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i+m*2-1];
				pxy[i*i+m*2-1] = hh*pxy[i*i+m*2-1]+deriv_h*(alpha_x*py[i*i+m*2-1]+alpha_y*px[i*i+m*2-1])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i+m*2-1];
				pyy[i*i+m*2-1] = hh*pyy[i*i+m*2-1]+2.*deriv_h*alpha_y*py[i*i+m*2-1]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i+m*2-1];
				pxz[i*i+m*2-1] = hh*pxz[i*i+m*2-1]+deriv_h*(alpha_x*pz[i*i+m*2-1]+alpha_z*px[i*i+m*2-1])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i+m*2-1];
				pyz[i*i+m*2-1] = hh*pyz[i*i+m*2-1]+deriv_h*(alpha_y*pz[i*i+m*2-1]+alpha_z*py[i*i+m*2-1])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i+m*2-1];
				pzz[i*i+m*2-1] = hh*pzz[i*i+m*2-1]+2.*deriv_h*alpha_z*pz[i*i+m*2-1]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i+m*2-1];
				pxx[i*i+m*2] = hh*pxx[i*i+m*2]+2.*deriv_h*alpha_x*px[i*i+m*2]+(dd_h*sqr(alpha_x)+deriv_h*alpha_xx)*p[i*i+m*2];
				pxy[i*i+m*2] = hh*pxy[i*i+m*2]+deriv_h*(alpha_x*py[i*i+m*2]+alpha_y*px[i*i+m*2])+(dd_h*alpha_x*alpha_y+deriv_h*alpha_xy)*p[i*i+m*2];
				pyy[i*i+m*2] = hh*pyy[i*i+m*2]+2.*deriv_h*alpha_y*py[i*i+m*2]+(dd_h*sqr(alpha_y)+deriv_h*alpha_yy)*p[i*i+m*2];
				pxz[i*i+m*2] = hh*pxz[i*i+m*2]+deriv_h*(alpha_x*pz[i*i+m*2]+alpha_z*px[i*i+m*2])+(dd_h*alpha_x*alpha_z+deriv_h*alpha_xz)*p[i*i+m*2];
				pyz[i*i+m*2] = hh*pyz[i*i+m*2]+deriv_h*(alpha_y*pz[i*i+m*2]+alpha_z*py[i*i+m*2])+(dd_h*alpha_y*alpha_z+deriv_h*alpha_yz)*p[i*i+m*2];
				pzz[i*i+m*2] = hh*pzz[i*i+m*2]+2.*deriv_h*alpha_z*pz[i*i+m*2]+(dd_h*sqr(alpha_z)+deriv_h*alpha_zz)*p[i*i+m*2];
				pz[i*i+m*2-1] = (f1*pz[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*FN2*cos_beta)/f4; 
				py[i*i+m*2-1] = (f1*py[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*cth*FN2*imag(z))/f4; 
				px[i*i+m*2-1] = (f1*px[i*i+m*2-1]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2-1]*cth*FN2*real(z))/f4; 
				pz[i*i+m*2] = (f1*pz[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*FN2*cos_beta)/f4; 
				py[i*i+m*2] = (f1*py[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*cth*FN2*imag(z))/f4; 
				px[i*i+m*2] = (f1*px[i*i+m*2]+(i-m)*(f1/f4*f5-f2)*p[i*i+m*2]*cth*FN2*real(z))/f4; 
				p[i*i+m*2-1] *= hh;
				p[i*i+m*2]   *= hh;
			}
		}
		P[0] = P[0];
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CMapi3DSpheroid::deriv_X(double * deriv, double f)
{
	if (! deriv || ! px) return;
	deriv[0] += px[0]*f;
	deriv[1] += px[1]*f;
	deriv[2] += px[2]*f;
	deriv[3] += px[3]*f;
}

void CMapi3DSpheroid::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! py) return;
	deriv[0] += py[0]*f;
	deriv[1] += py[1]*f;
	deriv[2] += py[2]*f;
	deriv[3] += py[3]*f;
}

void CMapi3DSpheroid::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! pz) return;
	deriv[0] += pz[0]*f;
	deriv[1] += pz[1]*f;
	deriv[2] += pz[2]*f;
	deriv[3] += pz[3]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CMapi3DSpheroidFull::deriv_X(double * deriv, double f)
{
	if (! deriv || ! px) return;
	for (int k = 0; k < NN; k++) deriv[k] += px[k]*f;
}

void CMapi3DSpheroidFull::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! py) return;
	for (int k = 0; k < NN; k++) deriv[k] += py[k]*f;
}

void CMapi3DSpheroidFull::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! pz) return;
	for (int k = 0; k < NN; k++) deriv[k] += pz[k]*f;
}
