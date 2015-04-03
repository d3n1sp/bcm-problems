#include "stdafx.h"
#include "cacou.h"

/////////////////////////////////////////////////////////////////
//																					//
//   BASIC SYSTEM of ELLIPSOIDAL MULTIPOLES (acoustic waves)   //
//																					//
/////////////////////////////////////////////////////////////////
///////////////////////////////////////
//...parametrization of the multipoles;
void CAcou3DEll::parametrization(double * P, int m_dop)
{
	if (! p) {
		p = (double *)new_struct((NN_dop = freedom(N))*sizeof(double));
	}
	if (P) {
		int m0(0), m1(2), m2, i, i0, i1, i2, m;
		double  X  = P[0]*R0_inv,
				  Y  = P[1]*R0_inv,
				  Z  = P[2]*R0_inv;
		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

//////////////////////////////////////////////
//...calculation direction of the integration;
		if ( param[0] > EE_dop) { //...индексация -- [n*(n+1)+m*2], [n*(n+1)+m*2+1];

//////////////////////////////////////////////////
//...симметричное и антисимметричное по Z решение;
			p[0] = cos(param[0]*(Z-param[4]*R0_inv)); 
			p[1] = sin(param[0]*(Z-param[4]*R0_inv))*param[1]; 
			for (m = 0; m < N+m_dop; m++, m0 = m1, m1 = m2) {
				for (i0 = (i1 = m0)-m*2, i = m; i < N; i++, i0 = i1, i1 = i2) {
					i2 = i1+2*(i+1);
					p[i2+m*2] = X*p[i1+m*2]-(i-m)*Z*p[i0+m*2+1];
					p[i2+m*2+1] = X*p[i1+m*2+1]+sqr(param[1])*(
									(i-m)*(Z*p[i0+m*2]-(i-m)*p[i0+m*2+1]+(i-m-1)*X*p[i0+(m-i)*2+3])-
									(m-1)*m*(p[i0+m*2-3]-X*p[i0+(m-i)*2-1]));
				}
				m2 = m1+2*(m+2);
				p[m2-2] = Y*p[m1-2]-m*Z*p[m0-1];
				p[m2-1] = Y*p[m1-1]+m*sqr(param[1])*(Z*p[m0-2]-m*p[m0-1]+(m-1)*Y*p[m0-2*m-1]);
			}
		}
		else {
///////////////////////////////////////////////////////////////
//...гармоничекое симметричное и антисимметричное по Z решение;
			p[0] = 1; 
			p[1] = Z; 
			for (m = 0; m < N+m_dop; m++, m0 = m1, m1 = m2) {
				for (i0 = (i1 = m0)-m*2, i = m; i < N; i++, i0 = i1, i1 = i2) {
					i2 = i1+2*(i+1);
					p[i2+m*2+1] = ((Z*(p[i2+m*2] = X*p[i1+m*2]-(i-m)*Z*p[i0+m*2+1])+(i-m+1)*X*p[i1+m*2+1])-
									(m-1)*m*(p[i2+m*2-3]-X*p[i1+m*2-3])/(i-m+2))/(i-m+2);
				}
				m2 = m1+2*(m+2);
				p[m2-1] = (Z*(p[m2-2] = Y*p[m1-2]-m*Z*p[m0-1])+(m+1)*Y*p[m1-1])/(m+2);
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
void CAcou3DEll::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! px && ! py && ! pz) {
		px = (double *)new_struct((NN_grad = NN_dop)*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
		pz = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py && pz) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		memset (pz, 0, NN_grad*sizeof(double));
		deriv_X(px);
		deriv_Y(py);
		deriv_Z(pz);
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
void CAcou3DEll::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (double *)new_struct((NN_hess = NN_grad)*sizeof(double));
		pxy = (double *)new_struct( NN_hess*sizeof(double));
		pyy = (double *)new_struct( NN_hess*sizeof(double));
		pxz = (double *)new_struct( NN_hess*sizeof(double));
		pyz = (double *)new_struct( NN_hess*sizeof(double));
		pzz = (double *)new_struct( NN_hess*sizeof(double));
	}
	if (P && pxx && pxy && pyy && pxz && pyz && pzz) {
		memset (pxx, 0, NN_hess*sizeof(double));
		memset (pxy, 0, NN_hess*sizeof(double));
		memset (pyy, 0, NN_hess*sizeof(double));
		memset (pxz, 0, NN_hess*sizeof(double));
		memset (pyz, 0, NN_hess*sizeof(double));
		memset (pzz, 0, NN_hess*sizeof(double));
		swap(p, px);
			deriv_X(pxx);
		swap(p, px);
		swap(p, py);
			deriv_X(pxy);
			deriv_Y(pyy);
		swap(p, py);
		swap(p, pz);
			deriv_X(pxz);
			deriv_Y(pyz);
			deriv_Z(pzz);
		swap(p, pz);
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DEll::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int m, jj, ii = N*(N+1), i = N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*i, m = i-1; m >= 0; m--) {
			deriv[ii+m*2+1] += (i-m)*p[jj+m*2+1]*f;
			deriv[ii+m*2]   += (i-m)*p[jj+m*2]*f;
		}
	}
}

void CAcou3DEll::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int m, jj, ii = N*(N+1), i = N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*i, m = i; m > 0; m--) {
			deriv[ii+m*2+1] += m*p[jj+m*2-1]*f;
			deriv[ii+m*2]   += m*p[jj+m*2-2]*f;
		}
	}
}

void CAcou3DEll::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int m, jj, ii = N*(N+1), i = N; i >= 0; i--, ii = i*(i+1)) {
		for (jj = (i-2)*(i-1), m = i; m >= 0; m--) {
			deriv[ii+m*2+1] += p[ii+m*2]*f;
			deriv[ii+m*2]	 -= p[ii+m*2+1]*param[2]*f;
			if (m < i-1) deriv[ii+m*2] -= (i-m)*(i-m-1)*p[jj+m*2+1]*f;
			if (m > 1)	 deriv[ii+m*2] -= m*(m-1)*p[jj+m*2-3]*f;
		}
	}
}

///////////////////////////////////////////////////////////
//																			//
//  ACOUSTIC MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
//																			//
///////////////////////////////////////////////////////////
///////////////////////////////////
//...realization of the multipoles;
void CAcou3DPoly::release()
{
	delete_struct(au);
	delete_struct(E1);
	CShape<double>::release();
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CAcou3DPoly::parametrization(double * P, int m_dop)
{
	//if (id_cmpl == 2) return;
	if (! p) {
			p = (double *)new_struct((NN_dop = freedom(N+m_dop))*sizeof(double));
///////////////////////
//...technical support;
			delete_struct(au); au = (double *)new_struct(( N+1+m_dop)*sizeof(double));        
			delete_struct(E1); E1 = (double *)new_struct(((N+m_dop)/2+1)*sizeof(double));        
	}
	if (P && p) {
		int  i, l, nm, m;
		complex z = comp(P[0], P[1])*R0_inv, zm = comp(1.);
		double rr = abs(z), r_inv,
				r0  = rr*.25,
				Z   = P[2]*R0_inv,
				Z0  = Z*Z, f1, f2, d, F;
		if (rr > EE)
			z *= (r_inv = 1./rr); else 
			z  = comp (r_inv = 1.);

		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		if (au) {
			const int descrpt[51] = {30, 36, 30, 36, 42, 42, 42, 42, 49, 56, 49,
												 56, 56, 64, 64, 64, 81, 81, 90, 90,110,
												100,100,121,132,144,156,182,196,196,196,
												225,240,272,272,289,289,272,342,342,361,
												441,380,462,462,462,484,506,552,552,600
			};
			double z = param[0]*rr, z_inv = param[1]*r_inv, f = 2.*param[1], h0, h1, rr0 = .25*z*z;
			int    m;
			for (m = 0, h0 = 1., h1 = sqrt(M_2_PI*z_inv); m <= N+m_dop; m++, h0 *= rr, h1 *= m*f)
			if ((m <= 50 && (int)rr0 < descrpt[m] || m >= 51 && (int)rr0 < (int)m*m*.25)) {
			  double f = 1., P = f;
			  int    k = 1;
			  do {
					 P += (f *= -rr0/((double)k*(m+k))); ++k;
			  }
			  while (fabs(f) > EE);
			  au[m] = h0*P;
			}
			else {
				double P = 1., Q = 0., f = 1., rr_inv = z_inv*.5, Co = cos(z-(m+.5)*M_PI_2),
						Si = sin(z-(m+.5)*M_PI_2), h;
				int    k = 1;
				do {
					  h  =  f;
					  Q += (f *=  (m-k+.5)*(m+k-.5)*rr_inv/(double)k); ++k;
					  P += (f *= -(m-k+.5)*(m+k-.5)*rr_inv/(double)k); ++k;
				}
				while ((k <= m-.5 || fabs(f) > EE && fabs(f) < fabs(h)));
				au[m] = h1*(P*Co-Q*Si);
			}
		}
		for (f1 = 1., f2 = Z, i = 0; i <= N+m_dop; i++) {
			for ( E1[l = 0] = 1., nm = i/2; l < nm; l++)
					E1[l+1] = -r0*(i-2*l-1)*(i-2*l)/sqr(l+1.)*E1[l];
			for ( F  = E1[nm]*au[nm]*(d  = f1), l = nm-1; l >= 0; l--)
					F += E1[l ]*au[l ]*(d *= Z0);
			p[i*i] = F; swap(f1, f2);
		}
		for (m = 1; m <= N+m_dop; m++)
		for (zm *= z, f1 = 1., f2 = Z, i = m; i <= N+m_dop; i++) {
			for ( E1[l = 0] = 1., nm = (i-m)/2; l < nm; l++)
					E1[l+1] = -r0*(i-m-2*l-1)*(i-m-2*l)/((l+1.)*(m+l+1.))*E1[l];
			for ( F  = E1[nm]*au[m+nm]*(d  = f1), l = nm-1; l >= 0; l--)
					F += E1[l ]*au[m+l ]*(d *= Z0);
			p[i*i+m*2-1] = F*imag(zm);
			p[i*i+m*2]   = F*real(zm); swap(f1, f2);
		}
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DPoly::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	int ll = (N+1)*(N+1);
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-3]-.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+1]+param[2]*p[ll+m*2+1]))*f;
			deriv[ii+m*2]   += (m*p[jj+m*2-2]-.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+2]+param[2]*p[ll+m*2+2]))*f;
		}
		if (m > 0) {
			deriv[ii+1] -= (.125*((i-1)*(i-2)*p[jj+3]+param[2]*p[ll+3]))*f;
			deriv[ii+2] -= (.125*((i-1)*(i-2)*p[jj+4]+param[2]*p[ll+4])-p[jj])*f;
		}
		deriv[ii] -= .5*(i*(i-1)*p[jj+2]+param[2]*p[ll+2])*f;
	}
	deriv[0] -= .5*param[2]*p[ll+2]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DPoly::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	int ll = (N+1)*(N+1);
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-2]+.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+2]+param[2]*p[ll+m*2+2]))*f;
			deriv[ii+m*2]   -= (m*p[jj+m*2-3]+.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+1]+param[2]*p[ll+m*2+1]))*f;
		}
		if (m > 0) {
			deriv[ii+m*2-1] += (.125*((i-1)*(i-2)*p[jj+4]+param[2]*p[ll+4])+p[jj])*f;
			deriv[ii+m*2]   -= (.125*((i-1)*(i-2)*p[jj+3]+param[2]*p[ll+3]))*f;
		}
		deriv[ii] -= .5*(i*(i-1.)*p[jj+1]+param[2]*p[ll+1])*f;
	}
	deriv[0] -= .5*param[2]*p[ll+1]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DPoly::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int m, jj, ii = N*N, i = N; i >= 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
			deriv[ii+m*2-1] += (i-m)*p[jj+m*2-1]*f;
			deriv[ii+m*2]   += (i-m)*p[jj+m*2]*f;
		}
		deriv[ii] += i*p[jj]*f;
	}
}

////////////////////////////////////////////////////////////
//																			 //
//  ACOUSTIC MULTIPOLES for SPHERE INTERIOUR OR EXTERIOR  //
//																			 //
////////////////////////////////////////////////////////////
///////////////////////////////////
//...realization of the multipoles;
void CAcou3DZoom::release()
{
    delete_struct(au);
    CShape<double>::release();
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CAcou3DZoom::parametrization(double * P, int m_dop)
{
	//if (id_cmpl == 2) return;
	if (! p) {
			p = (double *)new_struct((NN_dop = freedom(N+m_dop))*sizeof(double));
///////////////////////
//...technical support;
			delete_struct(au); au = (double *)new_struct((N+1+m_dop)*sizeof(double));        
	}
	if (P && p) {
		int i, m;
		complex zz, zm, F1, F2;
		double R = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]),	
				rr = R*param[3], 
				X  = P[0], Y = P[1], Z = P[2], 
				d  = R > EE ? 1./R : 1.;
		X *= d;
		Y *= d;
		Z *= d; zz = comp(X, Y);

/////////////////////////////////////////////////////////////////////
//...preliminary transformation of the point along the normal vector;
		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

//////////////////////////////////////////////////////
//...calculation of the multipoles (and its rotation);
		if (au) {
			double z_inv, h0 = 1., h1 = 1., f = 0., rr0 = .25*rr*rr, mu = .5;
			int m, asympt;
			if (! id_inverse) {
			 if ((asympt = rr0 < .25) != 0) f = .5*(h0 = rr); else z_inv = 1./rr;
				for (m = 0; m <= N+m_dop; m++, h0 *= f/(mu += 1))
				if (asympt) {
					double f = 1., P = f;
					int    k = 1 ;
					do {
						 P += (f *= -rr0/(double)(k*(mu+k))); ++k;
					}
					while (fabs(f) > EE);
					au[m] = h0*P;
				}
				else {
					double P = 1., Q = 0., f = 1., rr0_inv = z_inv*.5,
							Co = cos(rr-(mu+.5)*M_PI_2),  Si = sin(rr-(mu+.5)*M_PI_2), h;
					int   k  = 1;
					do {
						  h  =  f;
						  Q += (f *=  (mu-k+.5)*(mu+k-.5)*rr0_inv/(double)k); ++k;
						  P += (f *= -(mu-k+.5)*(mu+k-.5)*rr0_inv/(double)k); ++k;
					}
					while (k <= fabs(mu-.5) || fabs(f) > EE && fabs(f) < fabs(h));
					au[m] = h1*(P*Co-Q*Si);
				}
			}
			else {
				z_inv = 1./rr; h1 = 1.; mu = -.5;
				for (m = 0; m <= N+m_dop; m++, h1 *= -1., mu -= 1.) {
					double P = 1., Q = 0., f = 1., rr_inv = z_inv*.5,
							Co = cos(rr-(mu+.5)*M_PI_2),  Si = sin(rr-(mu+.5)*M_PI_2), h;
					int   k  = 1;
					do {
						  h  =  f;
						  Q += (f *=  (mu-k+.5)*(mu+k-.5)*rr_inv/(double)k); ++k;
						  P += (f *= -(mu-k+.5)*(mu+k-.5)*rr_inv/(double)k); ++k;
					}
					while (k <= fabs(mu-.5) || fabs(f) > EE && fabs(f) < fabs(h));
					au[m] = h1*(P*Co-Q*Si);
				}
			}
		}
		for (F1 = (zm = comp(d)), F2 = comp(0.), i = 0; i < N+m_dop; i++) { //...Z-axis;
			p[i*i] = real(F1)*au[i];
			F2 = ((2.*i+1.)*Z*F1-i*F2)/(i+1.); swap(F1, F2);
		}    
		p[i*i] = real(F1)*au[i];
		for (m = 1; m <= N+m_dop; m++) {
			for (F1 = (zm *= zz), F2 = comp(0.), i = m; i < N+m_dop; i++) {
				p[i*i+m*2-1] = imag(F1)*au[i];
				p[i*i+m*2]   = real(F1)*au[i];
				F2 = ((2.*i+1.)*Z*F1-(i-m)*F2)/(i+m+1.); swap(F1, F2);
			}
			p[i*i+m*2-1] = imag(F1)*au[i];
			p[i*i+m*2]   = real(F1)*au[i];
		}
//////////////////////
/*///...for debugging;
	 double mu;
	 p[0] = cos(rr)*d;
//*
	 double mu;
	 p[0] = sin(rr)*d;
/**///...end of debugging;
//////////////////////////
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
void CAcou3DZoom::parametrization_grad(double * P, int m)
{
	parametrization(P, m+1);
	if (! px && ! py && ! pz) {
		px = (double *)new_struct((NN_grad = freedom(N+m))*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
		pz = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py && pz) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		memset (pz, 0, NN_grad*sizeof(double));
		N += m;
		deriv_X(px);
		deriv_Y(py);
		deriv_Z(pz);
		N -= m;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
void CAcou3DZoom::parametrization_hess(double * P, int m)
{
	parametrization_grad(P, m+1);
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (double *)new_struct((NN_hess = freedom(N+m))*sizeof(double));
		pxy = (double *)new_struct( NN_hess*sizeof(double));
		pyy = (double *)new_struct( NN_hess*sizeof(double));
		pxz = (double *)new_struct( NN_hess*sizeof(double));
		pyz = (double *)new_struct( NN_hess*sizeof(double));
		pzz = (double *)new_struct( NN_hess*sizeof(double));
	}
	if (P && pxx && pxy && pyy && pxz && pyz && pzz) {
		memset (pxx, 0, NN_hess*sizeof(double));
		memset (pxy, 0, NN_hess*sizeof(double));
		memset (pyy, 0, NN_hess*sizeof(double));
		memset (pxz, 0, NN_hess*sizeof(double));
		memset (pyz, 0, NN_hess*sizeof(double));
		memset (pzz, 0, NN_hess*sizeof(double));
		N += m;
		swap(p, px);
			deriv_X(pxx);
			deriv_Y(pxy);
			deriv_Z(pxz);
		swap(p, px);
		swap(p, py);
			deriv_Y(pyy);
			deriv_Z(pyz);
		swap(p, py);
		swap(p, pz);
			deriv_Z(pzz);
		swap(p, pz);
		N -= m;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DZoom::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv*param[0];
////////////////////////////////////
//...recurrent derivation (about X);
	deriv[0] -= p[3]*f;
	for (int ii = 1, jj = 0, ll = 4, i = 1; i <= N; ll += 2, i++) {
		deriv[ii++] -= p[ll+2]*f; ll++;
		if (i > 1) {
			double ff = .5*i*(i-1.)/(2.*i+1.);
			deriv[ii-1] -= ff*(p[jj+2]+p[ll+1])*f; jj++;
			if (i > 2) {
				double ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] -= (.5*p[ll+2]+ff*(p[jj+2]+p[ll+2])/(2.*i+1.))*f; ll++;
				deriv[ii++] -= (.5*p[ll+2]-((p[jj-1]+p[ll-2])-ff*(p[jj+3]+p[ll+2]))/(2.*i+1.))*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] -= (.5*p[ll+2]-(m*(p[jj]+p[ll-2])-ff*(p[jj+4]+p[ll+2]))/(2.*i+1.))*f; jj++; ll++;
					deriv[ii++] -= (.5*p[ll+2]-(m*(p[jj]+p[ll-2])-ff*(p[jj+4]+p[ll+2]))/(2.*i+1.))*f; jj++; ll++;
				}
				ff = (i-1.)/(2.*i+1.);
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
				ff = i/(2.*i+1.);
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
			}
			else {
				ff = (i-1.)/(2.*i+1.);
				deriv[ii++] -= 0.5*p[ll+2]*f; ll++;
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj-1]+p[ll-2]))*f; ll++;
				ff = i/(2.*i+1.);
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+2]-ff*(p[jj++]+p[ll-2]))*f; ll++;
			}
		}
		else {
			deriv[ii++] -= 0.5*p[ll+2]*f; ll++;
			deriv[ii++] -= (.5*p[ll+2]-(p[jj++]+p[ll-2])/(2.*i+1.))*f; ll++;
		}
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DZoom::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv*param[0];
////////////////////////////////////
//...recurrent derivation (about Y);
	deriv[0] -= p[2]*f;
	for (int ii = 1, jj = 0, ll = 4, i = 1; i <= N; ll += 2, i++) {
		deriv[ii++] -= p[ll+1]*f; ll++;
		if (i > 1) {
			double ff = .5*i*(i-1.)/(2.*i+1.);
			deriv[ii-1] -= ff*(p[jj+1]+p[ll])*f; jj++;
			if (i > 2) {
				double ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] += (.5*p[ll+3]+((p[jj-1]+p[ll-1])+ff*(p[jj+3]+p[ll+3]))/(2.*i+1.))*f; ll++;
				deriv[ii++] -= (.5*p[ll+1]+ff*(p[jj+2]+p[ll+1])/(2.*i+1.))*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] += (.5*p[ll+3]+(m*(p[jj+1]+p[ll-1])+ff*(p[jj+5]+p[ll+3]))/(2.*i+1.))*f; jj++; ll++;
					deriv[ii++] -= (.5*p[ll+1]+(m*(p[jj-1]+p[ll-3])+ff*(p[jj+3]+p[ll+1]))/(2.*i+1.))*f; jj++; ll++;
				}
				ff = (i-1.)/(2.*i+1.);
				deriv[ii++] += (.5*p[ll+3]+ff*(p[(++jj)  ]+p[ll-1]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+1]+ff*(p[(++jj)-2]+p[ll-3]))*f; ll++;
				ff = i/(2.*i+1.);
				deriv[ii++] += (.5*p[ll+3]+ff*(p[(++jj)  ]+p[ll-1]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+1]+ff*(p[(++jj)-2]+p[ll-3]))*f; ll++;
			}
			else {
				ff = (i-1.)/(2.*i+1.);
				deriv[ii++] += (.5*p[ll+3]+ff*(p[jj-1]+p[ll-1]))*f; ll++;
				deriv[ii++] -= 0.5*p[ll+1]*f; ll++;
				ff = i/(2.*i+1.);
				deriv[ii++] += (.5*p[ll+3]+ff*(p[(++jj)  ]+p[ll-1]))*f; ll++;
				deriv[ii++] -= (.5*p[ll+1]+ff*(p[(++jj)-2]+p[ll-3]))*f; ll++;
			}
		}
		else {
			deriv[ii++] += (.5*p[ll+3]+(p[jj++]+p[ll-1])/(2.*i+1.)	)*f; ll++;
			deriv[ii++] -= 0.5*p[ll+1]*f; ll++;
		}
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DZoom::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv*param[0];
////////////////////////////////////
//...recurrent derivation (about Z);
	deriv[0] -= p[1]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= N; ll += 2, i++) {
		double ff = i/(2.*i+1.);
		deriv[ii++] -= (p[ll]-ff*(p[jj++]+p[ll]))*f; ll++;
		for (int m = 1; m < i; m++) {
			ff = (i-m)/(2.*i+1.);
			deriv[ii++] -= (p[ll]-ff*(p[jj++]+p[ll]))*f; ll++;
			deriv[ii++] -= (p[ll]-ff*(p[jj++]+p[ll]))*f; ll++;
		}
		deriv[ii++] -= p[ll++]*f;
		deriv[ii++] -= p[ll++]*f;
	}
}

////////////////////////////////////////////////////
//																  //
//      SUPERPOSITION of COMPLEX PLANE WAVES      //
//																  //
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//...destructor of class of exponential plane waves;
CAcou3DWave::~CAcou3DWave (void)
{
    delete_struct(cz);
    delete_struct(sz);
    delete_struct(cy);
    delete_struct(sy);
}

//////////////////////////////////////////////////
//...calculation of directions of the plane waves;
int CAcou3DWave::add_expo_power(double X, double Y, double Z, int flag)
{
    double g, d;
    int    buf_size = 30;

/////////////////////////////////////////////////////////////
//...calculation new number of directions of the plane waves;
    while (M >= N) {
        double * new_cz = (double *)new_struct((N+buf_size)*sizeof(double)), 
               * new_sz = (double *)new_struct((N+buf_size)*sizeof(double)), 
               * new_cy = (double *)new_struct((N+buf_size)*sizeof(double)), 
               * new_sy = (double *)new_struct((N+buf_size)*sizeof(double));

        memcpy(new_cz, cz, N*sizeof(double));
        memcpy(new_sz, cz, N*sizeof(double));
        memcpy(new_cy, cz, N*sizeof(double));
        memcpy(new_sy, cz, N*sizeof(double));

        delete_struct(cz); cz = new_cz;
        delete_struct(sz); sz = new_sz;
        delete_struct(cy); cy = new_cy;
        delete_struct(sy); sy = new_sy; if (! cz || ! sz || ! cy || ! sy)  return(0);
        N += buf_size;
    }

/////////////////////////////////
//...adding new nedded direction;
    g     = sqrt (d = X*X+Y*Y);
    cy[M] = Z*(d = 1./sqrt(d+Z*Z));
    sy[M] = g* d;
    cz[M] = 1.;
    sz[M] = 0.;
    if (g > EE) {
        cz[M] = X*(g = 1./g);
        sz[M] = Y* g;
    }
    M++; return(1);
}

//////////////////////////////////////
//...initialization of the multipoles;
void CAcou3DWave::init2(int id_flag, int N, int dim)
{
	if (id_flag != 1) {
		double * new_cz = (double *)new_struct(N*sizeof(double)), 
				 * new_sz = (double *)new_struct(N*sizeof(double)), 
				 * new_cy = (double *)new_struct(N*sizeof(double)), 
				 * new_sy = (double *)new_struct(N*sizeof(double));

		memcpy(new_cz, cz, CAcou3DWave::N*sizeof(double));
		memcpy(new_sz, cz, CAcou3DWave::N*sizeof(double));
		memcpy(new_cy, cz, CAcou3DWave::N*sizeof(double));
		memcpy(new_sy, cz, CAcou3DWave::N*sizeof(double));

		delete_struct(cz); cz = new_cz;
		delete_struct(sz); sz = new_sz;
		delete_struct(cy); cy = new_cy;
		delete_struct(sy); sy = new_sy;
		CAcou3DWave::N = N;
	}
	if (id_flag != 0) {
		CAcou3DWave::id_dim = dim;
		NN = M;
		release();
	}
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CAcou3DWave::parametrization(double * P, int m_dop)
{
    //if (id_cmpl == 2) return;
    if (! p) {
          p = (double *)new_struct(2*NN*sizeof(double));
    }
    cs[0] = P[3];
    cs[1] = P[4];
    cs[2] = P[5];

////////////////////////////////////
//...calculation of the plane waves;
	if (! id_inverse)
	for ( int k = 0; k < M; k++) {
		double Z = P[2]*cy[k]+(P[0]*cz[k]+P[1]*sz[k])*sy[k];
		p[k]     = cos(param[3]*Z);
		p[k+NN]  = sin(param[3]*Z);
	}
//	if (id_inverse) cpy(p+NN, p_cpy); 
//	else            cpy(p,    p_cpy); 
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DWave::deriv_X(double * deriv, double f)
{
     if (! deriv) return;
     if (id_inverse)
     for (int k = 0; k < M; k++) {
          double f = param[3]*cz[k]*sy[k];
          deriv[k] = p[k]*f;
     }
     else
     for (int k = 0; k < M; k++) {
          double f = param[3]*cz[k]*sy[k];
          deriv[k] = -p[k+NN]*f;
     }
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DWave::deriv_Y(double * deriv, double f)
{
     if (! deriv) return;
     if (id_inverse)
     for (int k = 0; k < M; k++) {
          double f = param[3]*sz[k]*sy[k];
          deriv[k] = p[k]*f;
     }
     else
     for (int k = 0; k < M; k++) {
          double f = param[3]*sz[k]*sy[k];
          deriv[k] = -p[k+NN]*f;
     }
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DWave::deriv_Z(double * deriv, double f)
{
     if (! deriv) return;
     if (id_inverse)
     for (int k = 0; k < M; k++) {
          double f = param[3]*cy[k];
          deriv[k] = p[k]*f;
     }
     else
     for (int k = 0; k < M; k++) {
          double f = param[3]*cy[k];
          deriv[k] = -p[k+NN]*f;
     }
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DWave::deriv_N()
{
     //if (id_cmpl == 2) return;
     if (id_inverse)
     for (int k = 0; k < M; k++) {
			double f = param[3]*((cs[0]*cz[k]+cs[1]*sz[k])*sy[k]+cs[2]*cy[k]);
//			deriv[k] = p[k]*f;
     }
     else
     for (int k = 0; k < M; k++) {
			double f = param[3]*((cs[0]*cz[k]+cs[1]*sz[k])*sy[k]+cs[2]*cy[k]);
//			deriv[k] = -p[k+NN]*f;
     }
}

////////////////////////////////////////////////////////////
//																			 //
//  BASIC SYSTEM of ACOUSTIC MULTIPOLES (acoustic beam)   //
//																			 //
////////////////////////////////////////////////////////////
///////////////////////////////////
//...realization of the multipoles;
void CAcou3DBeam::release()
{
	delete_struct(pim);
	delete_struct(pxim);
	delete_struct(pyim);
	delete_struct(pzim);
	delete_struct(au);
	delete_struct(az);
	delete_struct(E1);
	CShape<double>::release();
}

void CAcou3DBeamZ::release()
{
	delete_struct(pim);
	delete_struct(pxim);
	delete_struct(pyim);
	delete_struct(pzim);
	delete_struct(az);
	delete_struct(E1);
	CShape<double>::release();
}

///////////////////////////////
//...setting of the parameters;
void CAcou3DBeam::set_shape(double R0, double kk, double kk_dop, double L1, double L2)
{
	CAcou3DBeam::R0 = R0;
	R0_inv     = R0 > EE ? 1./R0 : 1.;
	int   k    = -1;
	if (++k < size_of_param()) param[k] = (Param)((kk = sqrt(fabs(kk*kk-kk_dop*kk_dop)))*R0);
	if (++k < size_of_param()) param[k] = (Param)(param[0] > EE ? 1./param[0] : 1.);
	if (++k < size_of_param()) param[k] = (Param)(param[0]*param[0]);
	if (++k < size_of_param()) param[k] = (Param)kk;
	if (++k < size_of_param()) param[k] = (Param)((kk_dop = fabs(kk_dop))*R0);
	if (++k < size_of_param()) param[k] = (Param)kk_dop;
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CAcou3DBeam::parametrization(double * P, int m_dop)
{
	//if (id_cmpl == 2) return;
	if (! p) {
			p   = (double *)new_struct((NN_dop = freedom(N+m_dop))*sizeof(double));
			pim = (double *)new_struct( NN_dop*sizeof(double));
///////////////////////
//...technical support;
			delete_struct(au); au = (double  *)new_struct((N+1+m_dop)*sizeof(double));
			delete_struct(az); az = (double  *)new_struct((N+1+m_dop)*sizeof(double));
			delete_struct(E1); E1 = (complex *)new_struct((N+1+m_dop)*sizeof(complex));
	}
	if (P && p && pim) {
		int  i, j, l, ii, nm, m;
		if (id_inverse) { //...мнимая часть акустических мультиполей;
			for (l = NN_dop-1; l >= 0; l--) 
				pim[l] = -pim[l];
			swap(p, pim);
			return;
		}
		complex z = comp(P[0], P[1])*R0_inv, zm = comp(1.), F, A;
		double rr = abs(z),	r_inv,
            Z   = P[2]*R0_inv,
				Co  = cos(param[4]*Z),
            Si  = sin(param[4]*Z), d, f, ff;
		if (rr > EE)
			z *= (r_inv = 1./rr); else 
			z  = comp (r_inv = 1.);

		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		for (az[i = 0] = 1.; i < N+m_dop; i++) az[i+1] = az[i]*Z;
		if (au) {
			const int descrpt[51] = {30, 36, 30, 36, 42, 42, 42, 42, 49, 56, 49,
												 56, 56, 64, 64, 64, 81, 81, 90, 90,110,
												100,100,121,132,144,156,182,196,196,196,
												225,240,272,272,289,289,272,342,342,361,
												441,380,462,462,462,484,506,552,552,600
			};
			double z = param[0]*rr, z_inv = param[1]*r_inv, f = 2.*param[1], h0, h1, rr0 = .25*z*z;
			int    m;
			for (m = 0, h0 = 1., h1 = sqrt(M_2_PI*z_inv); m <= N+m_dop; m++, h0 *= rr, h1 *= m*f)
			if ((m <= 50 && (int)rr0 < descrpt[m] || m >= 51 && (int)rr0 < (int)m*m*.25)) {
			  double f = 1., P = f;
			  int    k = 1;
			  do {
					 P += (f *= -rr0/((double)k*(m+k))); ++k;
			  }
			  while (fabs(f) > EE);
			  au[m] = h0*P;
			}
			else {
				double P = 1., Q = 0., f = 1., rr_inv = z_inv*.5, Co = cos(z-(m+.5)*M_PI_2),
						Si = sin(z-(m+.5)*M_PI_2), h;
				int    k = 1;
				do {
					  h  =  f;
					  Q += (f *=  (m-k+.5)*(m+k-.5)*rr_inv/(double)k); ++k;
					  P += (f *= -(m-k+.5)*(m+k-.5)*rr_inv/(double)k); ++k;
				}
				while ((k <= m-.5 || fabs(f) > EE && fabs(f) < fabs(h)));
				au[m] = h1*(P*Co-Q*Si);
			}
		}
		for (i = 0; i <= N+m_dop; i++) {
			memset(E1, 0, (N+m_dop+1)*sizeof(complex)); 
			E1[nm = i] = 1.;
			F = (ff = 1.)*E1[nm]*az[nm]*au[0];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/sqr((double)j);
				A = comp(0.);
				for (d = az[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(2.*param[4]*E1[l+1]*comp(0., 1.)-(l+2.)*E1[l+2])*f;
					A += E1[l]*d;
					d *= Z;
				}
				E1[l] = (l+1.)*(2.*param[4]*E1[l+1]*comp(0., 1.))*f;
				A += E1[l]*d;
				F += (ff *= rr)*A*au[j];
			}
			p  [ii] = real(F)*Co+imag(F)*Si;
			pim[ii] = imag(F)*Co-real(F)*Si;
		}
		for (m = 1; m <= N+m_dop; m++)
		for (zm *= z, i = m; i <= N+m_dop; i++) {
			memset(E1, 0, (N+m_dop+1)*sizeof(complex));
			E1[nm = i-m] = 1.;
			F = (ff = 1.)*E1[nm]*az[nm]*au[m];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/((double)j*(m+j));
				A = comp(0.);
				for (d = az[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(2.*param[4]*E1[l+1]*comp(0., 1.)-(l+2.)*E1[l+2])*f;
					A    += E1[l]*d;
					d    *= Z;
				}
				E1[l] = (l+1.)*(2.*param[4]*E1[l+1]*comp(0., 1.))*f;
				A += E1[l]*d;
				F += (ff *= rr)*A*au[m+j];
			}
			p  [ii+m*2-1] = (real(F)*Co+imag(F)*Si)*imag(zm);
			p  [ii+m*2]   = (real(F)*Co+imag(F)*Si)*real(zm);
			pim[ii+m*2-1] = (imag(F)*Co-real(F)*Si)*imag(zm);
			pim[ii+m*2]   = (imag(F)*Co-real(F)*Si)*real(zm);
		}
		if (id_inverse) { //...мнимая часть акустических мультиполей (repetition);
			for (l = NN_dop-1; l >= 0; l--) 
				pim[l] = -pim[l];
			swap(p, pim);
		}
	}
}

void CAcou3DBeamZ::parametrization(double * P, int m_dop)
{
	//if (id_cmpl == 2) return;
	if (! p) {
			p   = (double *)new_struct((NN_dop = freedom(N))*sizeof(double));
			pim = (double *)new_struct( NN_dop*sizeof(double));
///////////////////////
//...technical support;
			delete_struct(az); az = (double  *)new_struct((N+1)*sizeof(double));
			delete_struct(E1); E1 = (complex *)new_struct((N+1)*sizeof(complex));
	}
	if (P && p && pim) {
		int  i, j, l, ii, nm, m;
		if (id_inverse) { //...мнимая часть акустических мультиполей;
			for (l = NN_dop-1; l >= 0; l--)
				pim[l] = -pim[l];
			swap(p, pim);
			return;
		}
		complex z = comp(P[0], P[1])*R0_inv, zm = comp(1.), F, A;
		double rr = abs(z),
				r0  = rr*rr,
				Z   = P[2]*R0_inv,
				Co  = cos(param[0]*Z),
				Si  = sin(param[0]*Z), d, f, ff;
		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		for (az[i = 0] = 1.; i < N; i++) az[i+1] = az[i]*Z;
		for (i = 0; i <= N; i++) {
			memset(E1, 0, (N+1)*sizeof(complex));
			E1[nm = i] = 1.;
			F = (ff = 1.)*E1[nm]*az[nm];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/sqr((double)j);
				A = comp(0.);
				for (d = az[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(2.*param[0]*E1[l+1]*comp(0., 1.)-(l+2.)*E1[l+2])*f;
					A    += E1[l]*d;
					d    *= Z;
				}
				E1[l] = (l+1.)*(2.*param[0]*E1[l+1]*comp(0., 1.))*f;
				A    += E1[l]*d;
				F    += (ff *= r0)*A;
			}
			p  [ii] = real(F)*Co+imag(F)*Si;
			pim[ii] = imag(F)*Co-real(F)*Si;
		}
		for (m = 1; m <= N; m++)
		for (zm *= z, i = m; i <= N; i++) {
			memset(E1, 0, (N+1)*sizeof(complex));
			E1[nm = i-m] = 1.;
			F = (ff = 1.)*E1[nm]*az[nm];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/((double)j*(m+j));
				A = comp(0.);
				for (d = az[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(2.*param[0]*E1[l+1]*comp(0., 1.)-(l+2.)*E1[l+2])*f;
					A    += E1[l]*d;
					d    *= Z;
				}
				E1[l] = (l+1.)*(2.*param[0]*E1[l+1]*comp(0., 1.))*f;
				A    += E1[l]*d;
				F    += (ff *= r0)*A;
			}
			p  [ii+m*2-1] = (real(F)*Co+imag(F)*Si)*imag(zm);
			p  [ii+m*2]   = (real(F)*Co+imag(F)*Si)*real(zm);
			pim[ii+m*2-1] = (imag(F)*Co-real(F)*Si)*imag(zm);
			pim[ii+m*2]   = (imag(F)*Co-real(F)*Si)*real(zm);
		}
		if (id_inverse) { //...мнимая часть акустических мультиполей (repetition);
			for (l = NN_dop-1; l >= 0; l--)
				pim[l] = -pim[l];
			swap(p, pim);
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
void CAcou3DBeam::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! px && ! py && ! pz && ! pxim && ! pyim && ! pzim) {
		px = (double *)new_struct((NN_grad = freedom(N+m_dop))*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
		pz = (double *)new_struct( NN_grad*sizeof(double));
		pxim = (double *)new_struct( NN_grad*sizeof(double));
		pyim = (double *)new_struct( NN_grad*sizeof(double));
		pzim = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py && pz && pxim && pyim && pzim) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		memset (pz, 0, NN_grad*sizeof(double));
		memset (pxim, 0, NN_grad*sizeof(double));
		memset (pyim, 0, NN_grad*sizeof(double));
		memset (pzim, 0, NN_grad*sizeof(double));
		N += m_dop;
		swap(p, pim); param[0] = -param[0];
			deriv_X(pxim);
			deriv_Y(pyim);
			deriv_Z(pzim);
		swap(p, pim); param[0] = -param[0];
			deriv_X(px);
			deriv_Y(py);
			deriv_Z(pz);
		N -= m_dop;
	}
}

void CAcou3DBeamZ::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! px && ! py && ! pz && ! pxim && ! pyim && ! pzim) {
		px = (double *)new_struct((NN_grad = NN_dop)*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
		pz = (double *)new_struct( NN_grad*sizeof(double));
		pxim = (double *)new_struct( NN_grad*sizeof(double));
		pyim = (double *)new_struct( NN_grad*sizeof(double));
		pzim = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py && pz && pxim && pyim && pzim) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		memset (pz, 0, NN_grad*sizeof(double));
		memset (pxim, 0, NN_grad*sizeof(double));
		memset (pyim, 0, NN_grad*sizeof(double));
		memset (pzim, 0, NN_grad*sizeof(double));
		swap(p, pim); param[0] = -param[0];
			deriv_X(pxim);
			deriv_Y(pyim);
			deriv_Z(pzim);
		swap(p, pim); param[0] = -param[0];
			deriv_X(px);
			deriv_Y(py);
			deriv_Z(pz);
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
void CAcou3DBeam::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (double *)new_struct((NN_hess = freedom(N+m_dop))*sizeof(double));
		pxy = (double *)new_struct( NN_hess*sizeof(double));
		pyy = (double *)new_struct( NN_hess*sizeof(double));
		pxz = (double *)new_struct( NN_hess*sizeof(double));
		pyz = (double *)new_struct( NN_hess*sizeof(double));
		pzz = (double *)new_struct( NN_hess*sizeof(double));
	}
	if (P && pxx && pxy && pyy && pxz && pyz && pzz) {
		memset (pxx, 0, NN_hess*sizeof(double));
		memset (pxy, 0, NN_hess*sizeof(double));
		memset (pyy, 0, NN_hess*sizeof(double));
		memset (pxz, 0, NN_hess*sizeof(double));
		memset (pyz, 0, NN_hess*sizeof(double));
		memset (pzz, 0, NN_hess*sizeof(double));
		N += m_dop;
		swap(p, px); swap(pim, pxim);
			deriv_X(pxx);
		swap(p, px); swap(pim, pxim);
		swap(p, py); swap(pim, pyim);
			deriv_X(pxy);
			deriv_Y(pyy);
		swap(p, py); swap(pim, pyim);
		swap(p, pz); swap(pim, pzim);
			deriv_X(pxz);
			deriv_Y(pyz);
			deriv_Z(pzz);
		swap(p, pz); swap(pim, pzim);
		N -= m_dop;
	}
}

void CAcou3DBeamZ::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (double *)new_struct((NN_hess = NN_grad)*sizeof(double));
		pxy = (double *)new_struct( NN_hess*sizeof(double));
		pyy = (double *)new_struct( NN_hess*sizeof(double));
		pxz = (double *)new_struct( NN_hess*sizeof(double));
		pyz = (double *)new_struct( NN_hess*sizeof(double));
		pzz = (double *)new_struct( NN_hess*sizeof(double));
	}
	if (P && pxx && pxy && pyy && pxz && pyz && pzz) {
		memset (pxx, 0, NN_hess*sizeof(double));
		memset (pxy, 0, NN_hess*sizeof(double));
		memset (pyy, 0, NN_hess*sizeof(double));
		memset (pxz, 0, NN_hess*sizeof(double));
		memset (pyz, 0, NN_hess*sizeof(double));
		memset (pzz, 0, NN_hess*sizeof(double));
		swap(p, px); swap(pim, pxim);
			deriv_X(pxx);
		swap(p, px); swap(pim, pxim);
		swap(p, py); swap(pim, pyim);
			deriv_X(pxy);
			deriv_Y(pyy);
		swap(p, py); swap(pim, pyim);
		swap(p, pz); swap(pim, pzim);
			deriv_X(pxz);
			deriv_Y(pyz);
			deriv_Z(pzz);
		swap(p, pz); swap(pim, pzim);
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DBeam::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	int ll = (N+1)*(N+1);
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-3]-.25/(m+1.)*((i-m)*((i-m-1)*p[jj+m*2+1]+param[4]*pim[ii+m*2+1]*2.)+param[2]*p[ll+m*2+1]))*f;
			deriv[ii+m*2]   += (m*p[jj+m*2-2]-.25/(m+1.)*((i-m)*((i-m-1)*p[jj+m*2+2]+param[4]*pim[ii+m*2+2]*2.)+param[2]*p[ll+m*2+2]))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1.)*p[jj+i*2-5]-.25/(double)i*(param[4]*pim[ii+i*2-1]*2.+param[2]*p[ll+i*2-1]))*f;
			deriv[ii+i*2-2] += ((i-1.)*p[jj+i*2-4]-.25/(double)i*(param[4]*pim[ii+i*2]*2.+param[2]*p[ll+i*2]))*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += (i*p[jj+i*2-3]-.25/(i+1.)*param[2]*p[ll+i*2+1])*f;
			deriv[ii+i*2]   += (i*p[jj+i*2-2]-.25/(i+1.)*param[2]*p[ll+i*2+2])*f;

			deriv[ii+1] -= (.125*((i-1)*((i-2)*p[jj+3]+param[4]*pim[ii+3]*2.)+param[2]*p[ll+3]))*f;
			deriv[ii+2] -= (.125*((i-1)*((i-2)*p[jj+4]+param[4]*pim[ii+4]*2.)+param[2]*p[ll+4])-p[jj])*f;

			deriv[ii] -= .5*(i*((i-1)*p[jj+2]+param[4]*pim[ii+2]*2.)+param[2]*p[ll+2])*f;
		}
		if (i == 1) {
			deriv[ii+1] -= (.125*param[2]*p[ll+3])*f;
			deriv[ii+2] -= (.125*param[2]*p[ll+4]-p[jj])*f;

			deriv[ii] -= .5*(param[4]*pim[ii+2]*2.+param[2]*p[ll+2])*f;
		}
	}
	deriv[0] -= .5*param[2]*p[ll+2]*f;
}

void CAcou3DBeamZ::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-3]-.25/(m+1.)*(i-m)*((i-m-1)*p[jj+m*2+1]+param[0]*pim[ii+m*2+1]*2.))*f;
			deriv[ii+m*2]   += (m*p[jj+m*2-2]-.25/(m+1.)*(i-m)*((i-m-1)*p[jj+m*2+2]+param[0]*pim[ii+m*2+2]*2.))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1.)*p[jj+i*2-5]-.25/(double)i*param[0]*pim[ii+i*2-1]*2.)*f;
			deriv[ii+i*2-2] += ((i-1.)*p[jj+i*2-4]-.25/(double)i*param[0]*pim[ii+i*2]*2.)*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += i*p[jj+i*2-3]*f;
			deriv[ii+i*2]   += i*p[jj+i*2-2]*f;

			deriv[ii+1] -= (.125*(i-1)*((i-2)*p[jj+3]+param[0]*pim[ii+3]*2.))*f;
			deriv[ii+2] -= (.125*(i-1)*((i-2)*p[jj+4]+param[0]*pim[ii+4]*2.)-p[jj])*f;

			deriv[ii] -= .5*i*((i-1)*p[jj+2]+param[0]*pim[ii+2]*2.)*f;
		}
		if (i == 1) {
			deriv[ii+2] += p[jj]*f;
			deriv[ii]   -= param[0]*pim[ii+2]*f;
		}
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DBeam::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	int ll = (N+1)*(N+1);
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-2]+.25/(m+1.)*((i-m)*((i-m-1.)*p[jj+m*2+2]+param[4]*pim[ii+m*2+2]*2.)+param[2]*p[ll+m*2+2]))*f;
			deriv[ii+m*2]   -= (m*p[jj+m*2-3]+.25/(m+1.)*((i-m)*((i-m-1.)*p[jj+m*2+1]+param[4]*pim[ii+m*2+1]*2.)+param[2]*p[ll+m*2+1]))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1)*p[jj+i*2-4]+.25/(double)i*(param[4]*pim[ii+i*2]*2.+param[2]*p[ll+i*2]))*f;
			deriv[ii+i*2-2] -= ((i-1)*p[jj+i*2-5]+.25/(double)i*(param[4]*pim[ii+i*2-1]*2.+param[2]*p[ll+i*2-1]))*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += (i*p[jj+i*2-2]+.25/(i+1.)*param[2]*p[ll+i*2+2])*f;
			deriv[ii+i*2]   -= (i*p[jj+i*2-3]+.25/(i+1.)*param[2]*p[ll+i*2+1])*f;

			deriv[ii+1] += (.125*((i-1.)*((i-2.)*p[jj+4]+param[4]*pim[ii+4]*2.)+param[2]*p[ll+4])+p[jj])*f;
			deriv[ii+2] -= (.125*((i-1.)*((i-2.)*p[jj+3]+param[4]*pim[ii+3]*2.)+param[2]*p[ll+3]))*f;

			deriv[ii] -= .5*(i*((i-1.)*p[jj+1]+param[4]*pim[ii+1]*2.)+param[2]*p[ll+1])*f;
		}
		if (i == 1) {
			deriv[ii+1] += (.125*param[2]*p[ll+4]+p[jj])*f;
			deriv[ii+2] -= (.125*param[2]*p[ll+3])*f;

			deriv[ii] -= .5*(param[4]*pim[ii+1]*2.+param[2]*p[ll+1])*f;
		}
	}
	deriv[0] -= .5*param[2]*p[ll+1]*f;
}

void CAcou3DBeamZ::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	for (int m, jj, ii = N*N, i = N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*p[jj+m*2-2]+.25/(m+1.)*(i-m)*((i-m-1.)*p[jj+m*2+2]+param[0]*pim[ii+m*2+2]*2.))*f;
			deriv[ii+m*2]   -= (m*p[jj+m*2-3]+.25/(m+1.)*(i-m)*((i-m-1.)*p[jj+m*2+1]+param[0]*pim[ii+m*2+1]*2.))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1)*p[jj+i*2-4]+.25/(double)i*param[0]*pim[ii+i*2]*2.)*f;
			deriv[ii+i*2-2] -= ((i-1)*p[jj+i*2-5]+.25/(double)i*param[0]*pim[ii+i*2-1]*2.)*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += i*p[jj+i*2-2]*f;
			deriv[ii+i*2]   -= i*p[jj+i*2-3]*f;

			deriv[ii+1] += (.125*(i-1.)*((i-2.)*p[jj+4]+param[0]*pim[ii+4]*2.)+p[jj])*f;
			deriv[ii+2] -= (.125*(i-1.)*((i-2.)*p[jj+3]+param[0]*pim[ii+3]*2.))*f;

			deriv[ii] -= .5*i*((i-1.)*p[jj+1]+param[0]*pim[ii+1]*2.)*f;
		}
		if (i == 1) {
			deriv[ii+1] += p[jj]*f;
			deriv[ii]   -= param[0]*pim[ii+1]*f;
		}
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CAcou3DBeam::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	for (int m, jj, ii = N*N, i = N; i >= 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 0; m--) {
			deriv[ii+m*2-1] += ((i-m)*p[jj+m*2-1]+param[4]*pim[ii+m*2-1])*f;
			deriv[ii+m*2]   += ((i-m)*p[jj+m*2]+param[4]*pim[ii+m*2])*f;
		}
		deriv[ii] += (i*p[jj]+param[4]*pim[ii])*f;
	}
}

void CAcou3DBeamZ::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p || ! pim) return;
	f *= R0_inv;
	for (int m, jj, ii = N*N, i = N; i >= 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 0; m--) {
			deriv[ii+m*2-1] += ((i-m)*p[jj+m*2-1]+param[0]*pim[ii+m*2-1])*f;
			deriv[ii+m*2]   += ((i-m)*p[jj+m*2]+param[0]*pim[ii+m*2])*f;
		}
		deriv[ii] += (i*p[jj]+param[0]*pim[ii])*f;
	}
}
