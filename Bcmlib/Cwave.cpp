#include "stdafx.h"
#include "cwave.h"

////////////////////////////////////////////////////////
//   HEAT MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
////////////////////////////////////////////////////////
//////////////////////////////////////
//...initialization of the multipoles;
void CWave2DPoly::init2(int N, int M, int dim)
{
	CShape<double>::N = N;
	CShape<double>::M = max(M, 1);
	CShape<double>::id_dim = dim;
	NN = freedom(N)*(M+1);
	release();
}

///////////////////////////////
//...setting of the parameters;
void CWave2DPoly::set_shape(double R0, double kk, double r, double L1, double L2, double L3)
{
	CWave2DPoly::R0 = R0;
	R0_inv     = R0 > EE ? 1./R0 : 1.;
	int   k    = -1;
	if (++k < size_of_param()) param[k] = (Param)(kk*sqr(R0));
	if (++k < size_of_param()) param[k] = (Param)(param[0] > EE ? 1./param[0] : 1.);
	if (++k < size_of_param()) param[k] = (Param)(kk);
}

///////////////////////////////////////
//...parametrization of the multipoles;
void CWave2DPoly::parametrization(double * P, int m_dop)
{
	if (! P) {
		delete_struct(p);
	}
	if (! p) {
			p  = (double *)new_struct((NN_dop = freedom(N)*(M+1)/*+1000*/)*sizeof(double));
	}
	if (P && p) {
		double   t = P[6];
		complex  z = comp(P[0], P[1])*R0_inv, zm = 1., w = conj(z)*param[0]*.5, F3, F4;
		int      l, ii, mn, mm = NN_dop, m; 

		cs[0] = P[3];
		cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		for ( p[0] = 1., p[1] = t, m = 1;  m <= N+M; m++) {

				p[ii = (m*2-1)*(M+1)] = imag(F4 = (zm *= z));
				p[ii+1] = real(zm);

				for ( mn = min(m+1, M), l = 1; l <= mn; l++) {
						if ( l <= m) F3 = comp(p[ii+l*2-1], p[ii+l*2-2]);
						else			 F3 = p[l-1];
						F4 = F4*w*(l-1.)/(m-l+2.)+F3*t;

						if (l == m+1) p[l] = real(F4);
						else {
							p[ii+l*2]   = imag(F4);
							p[ii+l*2+1] = real(F4);
						}
						F4 = F3;	ii = max(ii-2*(M+1), 0);
				}
		}
	}
}

void CWave3DPoly::parametrization(double * P, int m_dop)
{
	if (! P) delete_struct(p);
	if (! p) {
			p  = (double *)new_struct(((NN_dop = freedom(N+m_dop)*(M+1))+2*M*(N+m_dop+1))*sizeof(double));
	}
	if (P && p) {
		double   Z = P[2]*R0_inv, f3 = Z*Z, t = P[6];
		complex  z = comp(P[0], P[1])*R0_inv, zm, 
					w = conj(z)*param[0]*.25, F1, F2, F3;
		int      i, l, ii, jj, mn, m0, mm = NN_dop, m; 

		cs[0] = P[3];
		cs[1] = P[4];
		cs[2] = P[5];
		f3   += real(z)*real(z)+imag(z)*imag(z);

///////////////////////////////////
//...calculation of the multipoles;
		for ( F1 = (zm = comp(1.)),  F2 = comp(0.), i = 0; i < N+m_dop; i++) {
				F2 = ((2.*i+1.)*Z*(p[i*i*(M+1)] = real(F1))-i*f3*F2)/(i+1.); 
				swap(F1, F2);
		}
		for ( p[i*i*(M+1)] = real(F1), m = 1; m <= N+m_dop+M; m++) {
			for ( F1 = (zm *= z), F2 = comp(0.), mn = min(m, M), i = m; i <= N+m_dop; i++) {
					p[ii = i*i*(M+1)+2*m-1] = imag(F3 = F1);
					p[ii+1] = real(F1);
					F2 = ((2.*i+1.)*Z*F1-(i-m)*f3*F2)/(i+m+1.); 
					swap(F1, F2);
					for (l = 1; l <= mn; l++) {
						F3 = F3*w/(m-l+1.);
						jj = (ii -= (2*(i-l)+1)*M+2*l-(l == m ? 1 : 0))-(2*(i-l)+1);
						if (l == m) p[ii] = real(F3 += p[jj]*t);
						else {
							p[ii]   = imag(F3 += comp(p[jj+1], p[jj])*t);
							p[ii+1] = real(F3);
						}
					}
			}
			for (m0 = N+m_dop+mn; i <= m0; i++, mm += 2) {
				p[mm]	  = imag(F3 = F1);
				p[mm+1] = real(F1);
				if (i < m0) {
					F2 = ((2.*i+1.)*Z*F1-(i-m)*f3*F2)/(i+m+1.);
					swap(F1, F2);
				}
				for (ii = i*i*(M+1)+2*m-1, jj = mm, l = 1; l < i-N-m_dop; l++) {
					F3 = F3*w/(m-l+1.); ii -= 2*l+(2*(i-l)+1)*M; jj -= 2*(N+m_dop+min(m-l, M)-max(m-l, N+m_dop)+1);
					p[jj]   = imag(F3 += comp(p[jj+1], p[jj])*t);
					p[jj+1] = real(F3);
				}
				for (; l <= mn; l++) {
					F3 = F3*w/(m-l+1.);
					jj = (ii -= (2*(i-l)+1)*M+2*l-(l == m ? 1 : 0))-(2*(i-l)+1);
					if (l == m) p[ii] = real(F3 += p[jj]*t);
					else {
						p[ii]   = imag(F3 += comp(p[jj+1], p[jj])*t);
						p[ii+1] = real(F3);
					}
				}
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
void CWave2DPoly::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! P) {
		delete_struct(px);
		delete_struct(py);
	}
	if (! px && ! py) {
		px = (double *)new_struct((NN_grad = NN_dop)*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		deriv_X(px);
		deriv_Y(py);
	}
}

void CWave3DPoly::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! P) {
		delete_struct(px);
		delete_struct(py);
		delete_struct(pz);
	}
	if (! px && ! py && ! pz) {
		px = (double *)new_struct((NN_grad = freedom(N+m_dop)*(M+1))*sizeof(double));
		py = (double *)new_struct( NN_grad*sizeof(double));
		pz = (double *)new_struct( NN_grad*sizeof(double));
	}
	if (P && px && py && pz) {
		memset (px, 0, NN_grad*sizeof(double));
		memset (py, 0, NN_grad*sizeof(double));
		memset (pz, 0, NN_grad*sizeof(double));
		N += m_dop;
		deriv_X(px);
		deriv_Y(py);
		deriv_Z(pz);
		N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
void CWave2DPoly::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! P) {
		delete_struct(pxx);
		delete_struct(pxy);
		delete_struct(pyy);
	}
	if (! pxx && ! pxy && ! pyy) {
		pxx = (double *)new_struct((NN_hess = NN_grad)*sizeof(double));
		pxy = (double *)new_struct( NN_hess*sizeof(double));
		pyy = (double *)new_struct( NN_hess*sizeof(double));
	}
	if (P && pxx && pxy && pyy) {
		memset (pxx, 0, NN_hess*sizeof(double));
		memset (pxy, 0, NN_hess*sizeof(double));
		memset (pyy, 0, NN_hess*sizeof(double));
		swap(p, px);
			deriv_X(pxx);
		swap(p, px);
		swap(p, py);
			deriv_X(pxy);
			deriv_Y(pyy);
		swap(p, py);
	}
}

void CWave3DPoly::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! P) {
		delete_struct(pxx);
		delete_struct(pxy);
		delete_struct(pyy);
		delete_struct(pxz);
		delete_struct(pyz);
		delete_struct(pzz);
	}
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (double *)new_struct((NN_hess = freedom(N+m_dop)*(M+1))*sizeof(double));
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
		N -= m_dop;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CWave2DPoly::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv; 
	
	double ff;
	int m, jj, mm, ii, i;
	for (ii = (2*N-1)*(M+1), i = N; i > 1; i--, ii = jj) {
		for (jj = ii-(M+1)*2, mm = ii+(M+1)*2, m = M; m > 0; m--) {
			deriv[ii+m*2]   += (i*p[jj+m*2]  +(ff = .25*param[0]/(i+1.)*m)*p[mm+m*2-2])*f;
			deriv[ii+m*2+1] += (i*p[jj+m*2+1]+ ff*p[mm+m*2-1])*f;
		}
		deriv[ii]   += i*p[jj]*f;
		deriv[ii+1] += i*p[jj+1]*f;
	}
	if (i > 0) {
		for (mm = ii+(M+1)*2, m = M; m > 0; m--) {
			deriv[ii+m*2]   += (ff = .125*param[0]*m)*p[mm+m*2-2]*f;
			deriv[ii+m*2+1] += (p[m]+ff*p[mm+m*2-1])*f;
		}
		deriv[ii+1] += p[0]*f;
	}
	for (mm = (M+1), m = M; m > 0; m--)
		deriv[m] += .5*param[0]*m*p[mm+m*2-1]*f;
}

void CWave3DPoly::deriv_X(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int l = 0; l <= M; l++) {
		int ll = (N+1)*(N+1)*(M+1)+(2*N+3)*l;
		for (int m, jj, ii = N*N*(M+1)+(2*N+1)*l, i = N; i > 0; i--, ll = ii, ii = jj) {
			for (jj = (i-1)*(i-1)*(M+1)+(2*i-1)*l, m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*p[jj+m*2-3]-.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+1]-param[0]*l*p[ll+(m-i)*2-2]))*f;
				deriv[ii+m*2]   += (m*p[jj+m*2-2]-.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+2]-param[0]*l*p[ll+(m-i)*2-1]))*f;
			}
			if (m > 0) {
				deriv[ii+1] -= (.125*((i-1)*(i-2)*p[jj+3]-param[0]*l*p[ll-i*2]))*f;
				deriv[ii+2] -= (.125*((i-1)*(i-2)*p[jj+4]-param[0]*l*p[ll-i*2+1])-p[jj])*f;
			}
			deriv[ii] -= .5*(i*(i-1)*p[jj+2]-param[0]*l*p[ll-i*2-1])*f;
		}
		deriv[l] += .5*param[0]*l*p[ll-1]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CWave2DPoly::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv; 
	
	double ff;
	int m, jj, mm, ii, i;
	for (ii = (2*N-1)*(M+1), i = N; i > 1; i--, ii = jj) {
		for (jj = ii-(M+1)*2, mm = ii+(M+1)*2, m = M; m > 0; m--) {
			deriv[ii+m*2]   += (i*p[jj+m*2+1]-(ff = .25*param[0]/(i+1.)*m)*p[mm+m*2-1])*f;
			deriv[ii+m*2+1] -= (i*p[jj+m*2]  - ff*p[mm+m*2-2])*f;
		}
		deriv[ii]   += i*p[jj+1]*f;
		deriv[ii+1] -= i*p[jj]*f;
	}
	if (i > 0) {
		for (mm = ii+(M+1)*2, m = M; m > 0; m--) {
			deriv[ii+m*2]   += (p[m]-(ff = .125*param[0]*m)*p[mm+m*2-1])*f;
			deriv[ii+m*2+1] += ff*p[mm+m*2-2]*f;
		}
		deriv[ii] += p[0]*f;
	}
	for (mm = (M+1), m = M; m > 0; m--)
		deriv[m] += .5*param[0]*m*p[mm+m*2-2]*f;
}

void CWave3DPoly::deriv_Y(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int l = 0; l <= M; l++) {
		int ll = (N+1)*(N+1)*(M+1)+(2*N+3)*l;
		for (int m, jj, ii = N*N*(M+1)+(2*N+1)*l, i = N; i > 0; i--, ll = ii, ii = jj) {
			for (jj = (i-1)*(i-1)*(M+1)+(2*i-1)*l, m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*p[jj+m*2-2]+.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+2]-param[0]*l*p[ll+(m-i)*2-1]))*f;
				deriv[ii+m*2]   -= (m*p[jj+m*2-3]+.25/(m+1.)*((i-m)*(i-m-1)*p[jj+m*2+1]-param[0]*l*p[ll+(m-i)*2-2]))*f;
			}
			if (m > 0) {
				deriv[ii+1] += (.125*((i-1)*(i-2)*p[jj+4]-param[0]*l*p[ll-i*2+1])+p[jj])*f;
				deriv[ii+2] -= (.125*((i-1)*(i-2)*p[jj+3]-param[0]*l*p[ll-i*2]))*f;
			}
			deriv[ii] -= .5*(i*(i-1)*p[jj+1]-param[0]*l*p[ll-i*2-2])*f;
		}
		deriv[l] += .5*param[0]*l*p[ll-2]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CWave3DPoly::deriv_Z(double * deriv, double f)
{
	if (! deriv || ! p) return;
	f *= R0_inv;
	for (int l = 0; l <= M; l++)
	for (int m, jj, ii = N*N*(M+1)+(2*N+1)*l, i = N;   i > 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1)*(M+1)+(2*i-1)*l, m = i-1; m > 0; m--) {
			deriv[ii+m*2-1] += (i-m)*p[jj+m*2-1]*f;
			deriv[ii+m*2]   += (i-m)*p[jj+m*2]*f;
		}
		deriv[ii] += i*p[jj]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
void CWave2DPoly::deriv_t(double * deriv, double f)
{
	if (! deriv || ! p) return;

	int m, ii, i;
	for (ii = (2*N-1)*(M+1), i = N; i > 0; i--, ii -= (M+1)*2)
	for (m = M; m > 0; m--) {
		deriv[ii+m*2]   += m*p[ii+m*2-2]*f;
		deriv[ii+m*2+1] += m*p[ii+m*2-1]*f;
	}
	for (m = M; m > 0; m--) 
		deriv[m] += m*p[m-1]*f;
}

void CWave3DPoly::deriv_t(double * deriv, double f)
{
	if (! deriv || ! p) return;

	int m, l, ii, i;
	for (ii = (N+1)*(N+1)*(M+1)-1, i = N; i > 0; ii -= i*2+1, i--)
		for ( l = M; l > 0; l--) {
		for ( m = i; m > 0; m--) {
			deriv[ii] += l*p[ii-i*2-1]*f; ii--;
			deriv[ii] += l*p[ii-i*2-1]*f; ii--;
		}
		deriv[ii] += l*p[ii-i*2-1]*f; ii--;
	}
	for (l = M; l > 0; l--) {
		deriv[ii] += l*p[ii-i*2-1]*f; ii--;
	}
}
