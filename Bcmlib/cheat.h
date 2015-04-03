/*========================================*/
/*                 CHEAT                  */
/*========================================*/
#ifndef ___CHEAT___
#define ___CHEAT___

#include "cshapes.h"

////////////////////////////////////////////////////////
//																		//
//   HEAT MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
//																		//
////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//...class of 1D heat-trasfer multipoles with polynomial behaviour;
template <typename T>
class CHeat1DPoly : public CShape<T> {
public:
		int type() { return HT1D_POLY_SHAPE;}
		int freedom(int  m) { return(2);}
		int size_of_param() { return(3);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double AA = 0., double kx = 0., double fx = 0., double p3 = 0.) {
			CHeat1DPoly::R0 = R0;
			this->R0_inv  	= R0 > EE ? 1./R0 : 1.;
			int   k = -1;
			if (++k < size_of_param()) this->param[k] = (Param)(AA*sqr(this->R0_inv));
			if (++k < size_of_param()) this->param[k] = (Param)(kx*R0);
			if (++k < size_of_param()) this->param[k] = (Param)(fx);
		}
		void init2(int N, int M, int dim) {
			CShape<T>::N  = N;
			CShape<T>::M  = M;
			CShape<T>::id_dim = dim;
			this->NN = freedom(N)*(M+1);
			CShape<T>::release();
		}
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat1DPoly () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		};
};

///////////////////////////////////////////////////////////////////
//...class of 2D heat-trasfer multipoles with polynomial behaviour;
template <typename T>
class CHeat2DPoly : public CShape<T> {
public:
		int type() { return HT2D_POLY_SHAPE;}
		int freedom (int m) { return m*2+1; }
		int size_of_param() { return(2);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk, double fo, double p2, double p3) {
			CHeat2DPoly::R0 = R0;
			this->R0_inv  	= R0 > EE ? 1./R0 : 1.;
			int   k = -1;
//////////////////////////////////////////////////////////////////////////////////////
//...R0 убираем из kk*R0 для получения характерного значения J_0(kk) на границе блока;
			if (++k < size_of_param()) this->param[k] = (Param)(kk);
			if (++k < size_of_param()) this->param[k] = (Param)(fo*sqr(kk*this->R0_inv));
		}
		void init2(int N, int M, int dim) {
			CShape<T>::N  = N;
			CShape<T>::M  = M;
			CShape<T>::id_dim = dim;
			this->NN = freedom(N)*(M+1);
			CShape<T>::release();
		}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0);
		void parametrization_hess(double * P = NULL, int m = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat2DPoly () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		};
};

///////////////////////////////////////////////////////////////////
//...class of 3D heat-trasfer multipoles with polynomial behaviour;
template <typename T>
class CHeat3DPoly : public CHeat1DPoly<T> {
public:
		int type()  { return HT3D_POLY_SHAPE;}
		int freedom(int m) { return sqr(m+1);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0);
		void parametrization_hess(double * P = NULL, int m = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat3DPoly () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(this->size_of_param()*sizeof(Param));
		};
};

//////////////////////////////////////////////////
//...parametrization of the polynomial multipoles;
#define SHAPE_heat_AA               1. //...normalized temperature-transfer coefficient;
#undef  SHAPE_heat_AA                  //...param(0);
#define SHAPE_heat_KX               0. //...normalized temperature-oscillaton coefficient on direction X;
#undef  SHAPE_heat_KX                  //...param(1);
#define SHAPE_heat_fx               0. //...shift on direction X;
#undef  SHAPE_heat_fx                  //...param(2);
#define SHAPE_heat_kappa				0. //...parameter of cooling rate;
#undef  SHAPE_heat_kappa					//...param(3);
#define SHAPE_heat_fo					0. //...parameter of Fourier time;
#undef  SHAPE_heat_fo						//...param(4);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CHeat1DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! P) delete_struct(this->p);
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N)*(this->M+1))*sizeof(T));
	}
	if (P && this->p) {
		double fo = P[6]*this->param[0]*sqr(this->R0), X = P[0]*this->R0_inv, f3 = X*sqr(this->R0)*.5;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		int m;
		for (this->p[0] = 1., this->p[m = 1] = X; m <= this->M; m++) {
			this->p[2*m] = fo*this->p[m*2-2]+f3*this->p[m*2-1]; 
			this->p[2*m+1] = (2*m*fo*this->p[m*2-1]+X*this->p[m*2])/(2.*m+1.); 
		}
	}
}

template <typename T>
void CHeat2DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! P) delete_struct(this->p);
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N)*(this->M+1)+2*this->M)*sizeof(T));
	}
	if (P && this->p) {
		complex z = comp(P[0], P[1])*this->R0_inv;
		T		 zm, zm_i, f, sum, hh;
		double rr = norm(z), q = -rr*.25*sqr(this->param[0]), dd = exp(-P[6]*this->param[1]);
		int i, m;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		for (zm = dd, zm_i = 0., m = 0; m <= this->N; zm = (hh = zm)*real(z)-zm_i*imag(z), zm_i = hh*imag(z)+zm_i*real(z), m++) {
			sum  = (f = 1.); i = 1; 
			do {
				sum += (f *= q/(i*(m+i))); ++i;
			}
			while (fabs(f) > EE);

			if (m) this->p[m*2-1] = zm_i*sum;
			this->p[m*2] = zm*sum;

		}
	}
}

template <typename T>
void CHeat3DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! P) delete_struct(this->p);
	if (! this->p) {
			this->p  = (T *)new_struct(((this->NN_dop = freedom(this->N)*(this->M+1))+2*this->M*(this->N+1))*sizeof(T));
	}
	if (P && this->p) {
		double   Z = P[2]*this->R0_inv, f3 = Z*Z, fo = P[6]*this->param[0]*sqr(this->R0);
		T			zm, zm_i, F1, F2, F3, F1_i, F2_i, F3_i, hh, h1, h2;
		complex  z = comp(P[0], P[1])*this->R0_inv, w = conj(z)*sqr(this->R0)*.25;
		int      i, l, ii, jj, mn, m0, mm = this->NN_dop, m; 

		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		f3   += real(z)*real(z)+imag(z)*imag(z);

///////////////////////////////////
//...calculation of the multipoles;
		for ( F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N; i++) {
				F2 = ((2.*i+1.)*Z*(this->p[i*i*(this->M+1)] = F1)-i*f3*F2)/(i+1.); 
				swap(F1, F2);
		}
		for ( this->p[i*i*(this->M+1)] = F1, m = 1; m <= this->N+this->M; m++) {
			for ( F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), 
					F2 = 0., F2_i = 0., mn = min(m, this->M), i = m; i <= this->N; i++) {
					this->p[ii = i*i*(this->M+1)+2*m-1] = (F3_i = F1_i);
					this->p[ii+1] = (F3 = F1);
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2); swap(F1_i, F2_i);
					for (l = 1; l <= mn; l++) {
						F3 = ((h1 = F3)*real(w)-F3_i*imag(w))*(hh = 1./(m-l+1.)); F3_i = (h1*imag(w)+F3_i*real(w))*hh;
						jj = (ii -= (2*(i-l)+1)*this->M+2*l-(l == m ? 1 : 0))-(2*(i-l)+1);
						if (l == m) this->p[ii] = (F3 += this->p[jj]*fo);
						else {
							this->p[ii]   = (F3_i += this->p[jj]*fo);
							this->p[ii+1] = (F3 += this->p[jj+1]*fo);
						}
					}
			}
			for (m0 = this->N+mn; i <= m0; i++, mm += 2) {
				this->p[mm]	  = (F3_i = F1_i);
				this->p[mm+1] = (F3 = F1);
				if (i < m0) {
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); F2_i= (h1*F1_i-h2*F2_i)*hh;
					swap(F1, F2); swap(F1_i, F2_i);
				}
				for (ii = i*i*(this->M+1)+2*m-1, jj = mm, l = 1; l < i-this->N; l++) {
					F3 = ((h1 = F3)*real(w)-F3_i*imag(w))*(hh = 1./(m-l+1.)); F3_i = (h1*imag(w)+F3_i*real(w))*hh;
					ii -= 2*l+(2*(i-l)+1)*this->M; jj -= 2*(this->N+min(m-l, this->M)-max(m-l, this->N)+1);
					this->p[jj]   = (F3_i += this->p[jj]*fo);
					this->p[jj+1] = (F3 += this->p[jj+1]*fo);
				}
				for (; l <= mn; l++) {
					F3 = ((h1 = F3)*real(w)-F3_i*imag(w))*(hh = 1./(m-l+1.)); F3_i = (h1*imag(w)+F3_i*real(w))*hh;
					jj = (ii -= (2*(i-l)+1)*this->M+2*l-(l == m ? 1 : 0))-(2*(i-l)+1);
					if (l == m) this->p[ii] = (F3 += this->p[jj]*fo);
					else {
						this->p[ii]   = (F3_i += this->p[jj]*fo);
						this->p[ii+1] = (F3 += this->p[jj+1]*fo);
					}
				}
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CHeat1DPoly<T>::parametrization_grad(double * P, int m)
{
	parametrization(P, m+1);
	if (! P) delete_struct(this->px);
	if (! this->px) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m)*(this->M+1))*sizeof(T));
	}
	if (P && this->px) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		this->N += m;
		deriv_X(this->px);
		this->N -= m;
	}
}

template <typename T>
void CHeat2DPoly<T>::parametrization_grad(double * P, int m)
{
	parametrization(P, m+1);
	if (! P) {
		delete_struct(this->px);
		delete_struct(this->py);
	}
	if (! this->px && ! this->py) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m)*(this->M+1))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		this->N += m;
		deriv_X(this->px);
		deriv_Y(this->py);
		this->N -= m;
	}
}

template <typename T>
void CHeat3DPoly<T>::parametrization_grad(double * P, int m)
{
	parametrization(P, m+1);
	if (! P) {
		delete_struct(this->px);
		delete_struct(this->py);
		delete_struct(this->pz);
	}
	if (! this->px && ! this->py && ! this->pz) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m)*(this->M+1))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
		this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py && this->pz) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		this->N += m;
		deriv_X(this->px);
		deriv_Y(this->py);
		deriv_Z(this->pz);
		this->N -= m;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CHeat1DPoly<T>::parametrization_hess(double * P, int m)
{
	parametrization_grad(P, m+1);
	if (! P) delete_struct(this->pxx);
	if (! this->pxx) {
		this->pxx = (T *)new_struct((this->NN_hess = freedom(this->N+m)*(this->M+1))*sizeof(T));
	}
	if (P && this->pxx) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		this->N += m;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		this->N -= m;
	}
}

template <typename T>
void CHeat2DPoly<T>::parametrization_hess(double * P, int m)
{
	parametrization_grad(P, m+1);
	if (! P) {
		delete_struct(this->pxx);
		delete_struct(this->pxy);
		delete_struct(this->pyy);
	}
	if (! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = (T *)new_struct((this->NN_hess = freedom(this->N+m)*(this->M+1))*sizeof(T));
		this->pxy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyy = (T *)new_struct( this->NN_hess*sizeof(T));
	}
	if (P && this->pxx && this->pxy && this->pyy) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		this->N += m;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		this->N -= m;
	}
}

template <typename T>
void CHeat3DPoly<T>::parametrization_hess(double * P, int m)
{
	parametrization_grad(P, m+1);
	if (! P) {
		delete_struct(this->pxx);
		delete_struct(this->pxy);
		delete_struct(this->pyy);
		delete_struct(this->pxz);
		delete_struct(this->pyz);
		delete_struct(this->pzz);
	}
	if (! this->pxx && ! this->pxy && ! this->pyy && ! this->pxz && ! this->pyz && ! this->pzz) {
		this->pxx = (T *)new_struct((this->NN_hess = freedom(this->N+m)*(this->M+1))*sizeof(T));
		this->pxy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pxz = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyz = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pzz = (T *)new_struct( this->NN_hess*sizeof(T));
	}
	if (P && this->pxx && this->pxy && this->pyy && this->pxz && this->pyz && this->pzz) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		memset (this->pxz, 0, this->NN_hess*sizeof(T));
		memset (this->pyz, 0, this->NN_hess*sizeof(T));
		memset (this->pzz, 0, this->NN_hess*sizeof(T));
		this->N += m;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		swap(this->p, this->pz);
			deriv_X(this->pxz);
			deriv_Y(this->pyz);
			deriv_Z(this->pzz);
		swap(this->p, this->pz);
		this->N -= m;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat1DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv; 

	for (int m = this->M; m > 0; m--) {
		deriv[m*2]   += m*this->p[m*2-1]*sqr(this->R0)*f;
		deriv[m*2+1] += this->p[m*2]*f;
	}
	deriv[1] += this->p[0]*f;
}

template <typename T>
void CHeat2DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv; 

	T  ff; 
	int m, jj, mm, ii, i;
	for (ii = (2*this->N-1)*(this->M+1), i = this->N; i > 1; i--, ii = jj) {
		for (jj = ii-(this->M+1)*2, mm = ii+(this->M+1)*2, m = this->M; m > 0; m--) {
			deriv[ii+m*2]   += (i*this->p[jj+m*2]  +(ff = .25*sqr(this->R0)/(i+1.)*m)*this->p[mm+m*2-2])*f;
			deriv[ii+m*2+1] += (i*this->p[jj+m*2+1]+ ff*this->p[mm+m*2-1])*f;
		}
		deriv[ii]   += i*this->p[jj]*f;
		deriv[ii+1] += i*this->p[jj+1]*f;
	}
	if (i > 0) {
		for (mm = ii+(this->M+1)*2, m = this->M; m > 0; m--) {
			deriv[ii+m*2]   += (ff = .125*sqr(this->R0)*m)*this->p[mm+m*2-2]*f;
			deriv[ii+m*2+1] += (this->p[m]+ff*this->p[mm+m*2-1])*f;
		}
		deriv[ii+1] += this->p[0]*f;
	}
	for (mm = (this->M+1), m = this->M; m > 0; m--)
		deriv[m] += .5*sqr(this->R0)*m*this->p[mm+m*2-1]*f;
}

template <typename T>
void CHeat3DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;

	for (int l = 0; l <= this->M; l++) {
		int ll = (this->N+1)*(this->N+1)*(this->M+1)+(2*this->N+3)*l;
		for (int m, jj, ii = sqr(this->N)*(this->M+1)+(2*this->N+1)*l, i = this->N; i > 0; i--, ll = ii, ii = jj) {
			for (jj = (i-1)*(i-1)*(this->M+1)+(2*i-1)*l, m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+1]-l*sqr(this->R0)*this->p[ll+(m-i)*2-2]))*f;
				deriv[ii+m*2]   += (m*this->p[jj+m*2-2]-.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+2]-l*sqr(this->R0)*this->p[ll+(m-i)*2-1]))*f;
			}
			if (m > 0) {
				deriv[ii+1] -= (.125*((i-1)*(i-2)*this->p[jj+3]-l*sqr(this->R0)*this->p[ll-i*2]))*f;
				deriv[ii+2] -= (.125*((i-1)*(i-2)*this->p[jj+4]-l*sqr(this->R0)*this->p[ll-i*2+1])-this->p[jj])*f;
			}
			deriv[ii] -= .5*(i*(i-1)*this->p[jj+2]-l*sqr(this->R0)*this->p[ll-i*2-1])*f;
		}
		deriv[l] += .5*sqr(this->R0)*l*this->p[ll-1]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat2DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv; 
	
	T  ff;
	int m, jj, mm, ii, i;
	for (ii = (2*this->N-1)*(this->M+1), i = this->N; i > 1; i--, ii = jj) {
		for (jj = ii-(this->M+1)*2, mm = ii+(this->M+1)*2, m = this->M; m > 0; m--) {
			deriv[ii+m*2]   += (i*this->p[jj+m*2+1]-(ff = .25*sqr(this->R0)/(i+1.)*m)*this->p[mm+m*2-1])*f;
			deriv[ii+m*2+1] -= (i*this->p[jj+m*2]  - ff*this->p[mm+m*2-2])*f;
		}
		deriv[ii]   += i*this->p[jj+1]*f;
		deriv[ii+1] -= i*this->p[jj]*f;
	}
	if (i > 0) {
		for (mm = ii+(this->M+1)*2, m = this->M; m > 0; m--) {
			deriv[ii+m*2]   += (this->p[m]-(ff = .125*sqr(this->R0)*m)*this->p[mm+m*2-1])*f;
			deriv[ii+m*2+1] += ff*this->p[mm+m*2-2]*f;
		}
		deriv[ii] += this->p[0]*f;
	}
	for (mm = (this->M+1), m = this->M; m > 0; m--)
		deriv[m] += .5*sqr(this->R0)*m*this->p[mm+m*2-2]*f;
}

template <typename T>
void CHeat3DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;

	for (int l = 0; l <= this->M; l++) {
		int ll = sqr(this->N+1)*(this->M+1)+(2*this->N+3)*l;
		for (int m, jj, ii = sqr(this->N)*(this->M+1)+(2*this->N+1)*l, i = this->N; i > 0; i--, ll = ii, ii = jj) {
			for (jj = (i-1)*(i-1)*(this->M+1)+(2*i-1)*l, m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+2]-l*sqr(this->R0)*this->p[ll+(m-i)*2-1]))*f;
				deriv[ii+m*2]   -= (m*this->p[jj+m*2-3]+.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+1]-l*sqr(this->R0)*this->p[ll+(m-i)*2-2]))*f;
			}
			if (m > 0) {
				deriv[ii+1] += (.125*((i-1)*(i-2)*this->p[jj+4]-l*sqr(this->R0)*this->p[ll-i*2+1])+this->p[jj])*f;
				deriv[ii+2] -= (.125*((i-1)*(i-2)*this->p[jj+3]-l*sqr(this->R0)*this->p[ll-i*2]))*f;
			}
			deriv[ii] -= .5*(i*(i-1)*this->p[jj+1]-l*sqr(this->R0)*this->p[ll-i*2-2])*f;
		}
		deriv[l] += .5*sqr(this->R0)*l*this->p[ll-2]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat3DPoly<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;

	for (int l = 0; l <= this->M; l++)
	for (int m, jj, ii = sqr(this->N)*(this->M+1)+(2*this->N+1)*l, i = this->N;   i > 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1)*(this->M+1)+(2*i-1)*l, m = i-1; m > 0; m--) {
			deriv[ii+m*2-1] += (i-m)*this->p[jj+m*2-1]*f;
			deriv[ii+m*2]   += (i-m)*this->p[jj+m*2]*f;
		}
		deriv[ii] += i*this->p[jj]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat1DPoly<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->param[0]*sqr(this->R0);

	for (int m = this->M; m > 0; m--) {
		deriv[m*2]   += m*this->p[m*2-2]*f;
		deriv[m*2+1] += m*this->p[m*2-1]*f;
	}
}

template <typename T>
void CHeat2DPoly<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->param[0]*sqr(this->R0);

	int m, ii, i;
	for (ii = (2*this->N-1)*(this->M+1), i = this->N; i > 0; i--, ii -= (this->M+1)*2)
	for (m = this->M; m > 0; m--) {
		deriv[ii+m*2]   += m*this->p[ii+m*2-2]*f;
		deriv[ii+m*2+1] += m*this->p[ii+m*2-1]*f;
	}
	for (m = this->M; m > 0; m--) 
		deriv[m] += m*this->p[m-1]*f;
}

template <typename T>
void CHeat3DPoly<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->param[0]*sqr(this->R0);

	int m, l, ii, i;
	for (ii = sqr(this->N+1)*(this->M+1)-1, i = this->N; i > 0; ii -= i*2+1, i--)
		for ( l = this->M; l > 0; l--) {
		for ( m = i; m > 0; m--) {
			deriv[ii] += l*this->p[ii-i*2-1]*f; ii--;
			deriv[ii] += l*this->p[ii-i*2-1]*f; ii--;
		}
		deriv[ii] += l*this->p[ii-i*2-1]*f; ii--;
	}
	for (l = this->M; l > 0; l--) {
		deriv[ii] += l*this->p[ii-i*2-1]*f; ii--;
	}
}

//////////////////////////////////////////////////////////////
//																				//
//        HEAT MULTIPOLES FOR CELL WITH INCLUSIONS          //
//																				//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//...class of heat multipoles for one streep inclusion;
template <typename T>
class CHeat2DStreep : public CShape<T> {
protected:
		T * pim;
		complex * E1;
public:
		int type()          { return HT2D_STREEP_SHAPE;}
		int size_of_param() { return(5);};
		int freedom(int m)  { return (m+1)*2; }
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0, double b = 0., double kap = 0, double q = 0, double fo = 0) {
			CHeat2DStreep::R0 = R0;
			this->R0_inv = R0 > EE ? 1./R0 : 1.;
			int   k = -1;
//////////////////////////////////////////////////////////////////////////////////////
//...R0 убираем из kk*R0 для получения характерного значения exp(kk) на границе блока;
			if (++k < size_of_param()) this->param[k] = (Param)(kk);
			if (++k < size_of_param()) this->param[k] = (Param)(b*this->R0_inv);
			if (++k < size_of_param()) this->param[k] = (Param)(kap);
			if (++k < size_of_param()) this->param[k] = (Param)(q);
			if (++k < size_of_param()) this->param[k] = (Param)(fo*sqr(kk*this->R0_inv));
		}
		void release() {
			delete_struct(pim);
			delete_struct(E1);
			CShape<T>::release();
		}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat2DStreep () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			pim = NULL;
			E1  = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_streep_kappa          0. //...parameter of cooling rate in inclusion;
#undef  SHAPE_streep_kappa             //...param(0);
#define SHAPE_streep_b					0. //...semi-width of inclusion;
#undef  SHAPE_streep_b						//...param(1);
#define SHAPE_streep_kap            0. //...parameter of cooling rate in inclusion kap_I/kap_M;
#undef  SHAPE_streep_kap               //...param(2);
#define SHAPE_streep_q					0. //...parameter of conjunction coefficients lamb_I/lamb_M;
#undef  SHAPE_streep_q                 //...param(3);
#define SHAPE_streep_fo					0. //...parameter of Fourier time;
#undef  SHAPE_streep_fo                //...param(4);


////////////////////////////////////////////////////
//...class of heat multipoles for circle inclusions;
template <typename T>
class CHeat2DCircle : public CShape<T> {
protected:
		T * pim;
		complex * E1;
public:
		int type()          { return HT2D_CIRCLE_SHAPE;}
		int size_of_param() { return(5);};
		int freedom(int m)  { return m*2+1; }
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0, double rad = 0., double kap = 0, double q = 0, double fo = 0) {
			CHeat2DCircle::R0 = R0;
			this->R0_inv = R0 > EE ? 1./R0 : 1.;
  			int   k = -1;
//////////////////////////////////////////////////////////////////////////////////////
//...R0 убираем из kk*R0 для получения характерного значения J_0(kk) на границе блока;
			if (++k < size_of_param()) this->param[k] = (Param)(kk);
			if (++k < size_of_param()) this->param[k] = (Param)(rad*this->R0_inv);
			if (++k < size_of_param()) this->param[k] = (Param)(kap);
			if (++k < size_of_param()) this->param[k] = (Param)(q);
			if (++k < size_of_param()) this->param[k] = (Param)(fo*sqr(kk*this->R0_inv));
		}
		void release() {
			delete_struct(pim);
			delete_struct(E1);
			CShape<T>::release();
		}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0);
		void parametrization_hess(double * P = NULL, int m = 0);
		void cpy(T * p_ext);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat2DCircle () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			pim = NULL;
			E1  = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_circle_kappa          0. //...parameter of cooling rate in inclusion;
#undef  SHAPE_circle_kappa             //...param(0);
#define SHAPE_circle_rad				0. //...radius of inclusion;
#undef  SHAPE_circle_rad					//...param(1);
#define SHAPE_circle_kap            0. //...parameter of cooling rate in inclusion;
#undef  SHAPE_circle_kap               //...param(2);
#define SHAPE_circle_q					0. //...parameter of conjunction coefficients;
#undef  SHAPE_circle_q                 //...param(3);
#define SHAPE_circle_fo					0. //...parameter of Fourier time;
#undef  SHAPE_circle_fo                //...param(4);


////////////////////////////////////////////////////
//...class of heat multipoles for circle inclusions;
template <typename T>
class CHeat3DSphere : public CShape<T> {
protected:
		T * pim;
		complex * E1;
public:
		int type() { return HT3D_SPHERE_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(5);};
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0, double rad = 0., double kap = 0, double q = 0, double fo = 0) {
			CHeat3DSphere::R0 = R0;
			this->R0_inv = R0 > EE ? 1./R0 : 1.;
			int   k = -1;
//////////////////////////////////////////////////////////////////////////////////////
//...R0 убираем из kk*R0 для получения характерного значения J_0(kk) на границе блока;
			if (++k < size_of_param()) this->param[k] = (Param)(kk);
			if (++k < size_of_param()) this->param[k] = (Param)(rad*this->R0_inv);
			if (++k < size_of_param()) this->param[k] = (Param)(kap);
			if (++k < size_of_param()) this->param[k] = (Param)(q);
			if (++k < size_of_param()) this->param[k] = (Param)(fo*sqr(kk*this->R0_inv));
		}
		void release() {
			delete_struct(pim);
			delete_struct(E1);
			CShape<T>::release();
		}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
		void cpy(T * p_ext);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
		void deriv_t(T * deriv, double f = 1.);
//...constructor;
		CHeat3DSphere () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			pim = NULL;
			E1  = NULL;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_sphere_kappa          0. //...parameter of cooling rate in inclusion;
#undef  SHAPE_sphere_kappa             //...param(0);
#define SHAPE_sphere_rad				0. //...radius of inclusion;
#undef  SHAPE_sphere_rad					//...param(1);
#define SHAPE_sphere_kap            0. //...parameter of cooling rate in inclusion;
#undef  SHAPE_sphere_kap               //...param(2);
#define SHAPE_sphere_q					0. //...parameter of conjunction coefficients;
#undef  SHAPE_sphere_q                 //...param(3);
#define SHAPE_sphere_fo					0. //...parameter of Fourier time;
#undef  SHAPE_sphere_fo                //...param(4);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CHeat2DStreep<T>::parametrization(double * P, int m_dop)
{
	int mmm = (this->N+m_dop)/2;
	if (! this->p) {
			this->p	 = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
			pim = (T *)new_struct( this->NN_dop*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(E1); E1 = (complex *)new_struct(2*(mmm+1)*sizeof(complex));

			if (this->param[0]) { //...реализация сшивочных коэффициентов в негармоническом случае;
				double q = this->param[0]*this->param[1], C2 = cos(q), S2 = sin(q), kap_M = 1./this->param[0], sg,
						 d = this->param[2]*q, C1 = cos(d), S1 = sin(d), kap_I = 1./(this->param[0]*this->param[2]),
						 f = this->param[2]*this->param[3], g = this->param[2]*f, Q = 1./(C2*C1+f*S2*S1), S = 1./(S2*S1+f*C2*C1);
				complex * UM = (complex *)new_struct(2*(mmm+1)*sizeof(complex));
				int i, i2, m2, m;

				UM[0]  = comp(C2, C1);
				UM[1]  = comp(S2*kap_M, S1*kap_I);
				kap_M *= kap_M;
				kap_I *= kap_I;

				for (m = 0; m < mmm; m++) {
					UM[m*2+2] = this->param[1]*UM[m*2+1]/(m*2+2);
					UM[m*2+3] = ((m*2+1)*UM[m*2+1]-this->param[1]*UM[m*2])/(m*2+2);
					UM[m*2+3] = comp(real(UM[m*2+3])*kap_M, imag(UM[m*2+3])*kap_I);
				}

/////////////////////////////////////////////////////////////////////////////////
//...расчет коэффициентов -- E1[m*2] (четная), E1[m*2+1] (нечетная), E1 = A + iB;
				E1[0] = Q*comp(1., (S2*C1-f*C2*S1)*this->param[0]);
				E1[1] = S*comp(this->param[2], (C2*S1-f*S2*C1)*kap_M*this->param[0]);

				for (m2 = (m = 1)*2; m <= mmm; m++, m2 = m*2) {
					for ( sg = -1., i2 = (i = 1)*2; i <= m; i++, i2 = i*2) {
						E1[m2] += (sg = -sg)*comp(real(UM[i2+1])*imag(E1[m2-i2])-imag(UM[i2])*real(E1[m2-i2]),
									kap_M*real(UM[i2])*imag(E1[m2-i2])+
									 (g*imag(UM[i2+1])-kap_M*this->param[3]*imag(UM[i2-1]))*real(E1[m2-i2]));

						E1[m2+1] += sg*comp(real(UM[i2])*imag(E1[m2-i2+1])-imag(UM[i2+1])*real(E1[m2-i2+1]),
										(kap_M*real(UM[i2-1])-real(UM[i2+1]))*imag(E1[m2-i2+1])
										-kap_M*this->param[3]*imag(UM[i2])*real(E1[m2-i2+1]));
					}
					E1[m2] += sg*comp(real(UM[m2]), kap_M*real(UM[m2-1])-real(UM[m2+1]));
					E1[m2] = Q*comp(this->param[0]*S2*imag(E1[m2])-C2*real(E1[m2]),
							this->param[0]*(this->param[0]*C1*imag(E1[m2])+S1*real(E1[m2])*f));

					E1[m2+1] += sg*comp(real(UM[m2+1]), kap_M*real(UM[m2]));
					E1[m2+1] = -S*comp(this->param[0]*this->param[2]*(this->param[0]*C2*imag(E1[m2+1])+S2*real(E1[m2+1])),
											 this->param[0]*S1*imag(E1[m2+1])- C1*real(E1[m2+1])*f);
				}
				delete_struct(UM);
			}
			else { //...реализация сшивочных коэффициентов в гармоническом случае;
				double f = 1./this->param[3], g = this->param[1], sg;
				double * UM = (double *)new_struct(2*(mmm+1)*sizeof(double));
				int i, i2, m2, m;

				for (UM[0] = 1., UM[1] = g, m = 0; m < mmm; m++) {
					UM[m*2+2] = g*UM[m*2+1]/(m*2+2.);
					UM[m*2+3] = g*UM[m*2+2]/(m*2+3.);
				}

/////////////////////////////////////////////////////////////////////////////////
//...расчет коэффициентов -- E1[m*2] (четная), E1[m*2+1] (нечетная), E1 = A + iB;
				for (E1[0] = comp(1., 0.), E1[1] = comp(f, (f-1.)*g), m2 = (m = 1)*2; m <= mmm; m++, m2 = m*2) {
					for ( sg = -1., i2 = (i = 1)*2; i <= m; i++, i2 = i*2) {
						E1[m2] += (sg = -sg)*comp(UM[i2+1]*imag(E1[m2-i2])-UM[i2]*real(E1[m2-i2]),
														  UM[i2]*imag(E1[m2-i2])-this->param[3]*UM[i2-1]*real(E1[m2-i2]));

						E1[m2+1] += sg*comp(UM[i2]*imag(E1[m2-i2+1])-UM[i2+1]*real(E1[m2-i2+1]),
											 -f*UM[i2-1]*imag(E1[m2-i2+1])+UM[i2]*real(E1[m2-i2+1]));
					}
					E1[m2] += sg*comp(UM[m2], UM[m2-1]);
					E1[m2] = comp(g*imag(E1[m2])-real(E1[m2]), imag(E1[m2]));

					E1[m2+1] += sg*comp(UM[m2+1],  -f*UM[m2]);
					E1[m2+1] = comp(imag(E1[m2+1]), g*imag(E1[m2+1])+real(E1[m2+1]));
				}
				delete_struct(UM);
			}
	}
	if (P && this->p) {
		int i, i2, m2, m;
		double X = P[0]*this->R0_inv,
				 Y = P[1]*this->R0_inv, kap, q, fact, dd = exp(-P[6]*this->param[4]);
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		if (! (this->id_inverse = (fabs(Y) > this->param[1]))) {//...internal system;
			if (this->param[0]) { //...реализация функций в негармоническом случае;
				kap  = 1./(q = this->param[0]*this->param[2]);
				pim[0] = cos(q*Y)*dd;     this->p[0] = 0.; 
				pim[1] = sin(q*Y)*kap*dd; this->p[1] = 0.;
				kap *= kap;

				for (i2 = i = 0; i < this->N; i++, i2 = i*2) {
					pim[i2+2] = X*pim[i2]-i*Y*pim[i2-1]; this->p[i2+2] = 0.;
					pim[i2+3] = X*pim[i2+1]+i*kap*(Y*pim[i2-2]-i*pim[i2-1]+(i-1)*X*pim[i2-3]); this->p[i2+3] = 0.;
				}
			}
			else { //...реализация функций в гармоническом случае;
				pim[0] = dd;   this->p[0] = 0.; 
				pim[1] = Y*dd; this->p[1] = 0.;

				for (i2 = i = 0; i < this->N; i++, i2 = i*2) {
					pim[i2+2] =  X*pim[i2]-i*Y*pim[i2-1]; this->p[i2+2] = 0.;
					pim[i2+3] = (Y*pim[i2+2]+(i+1)*X*pim[i2+1])/(i+2.); this->p[i2+3] = 0.;
				}
			}
			for (fact = 1., m = 0; m <= mmm; m++, fact *= (m2+1)*(m2+2))
			for (q = fact, i2 = (i = (m2 = m*2))*2;  i <= this->N; i++, i2 = i*2) {
				this->p[i2] += q*pim[i2-2*m2]*real(E1[m2]);
				this->p[i2+1] += q*pim[i2-2*m2+1]*real(E1[m2+1]);
				q *= (i+1.)/(i-m2+1.);
			}
		}
		else {//...external system;
			if (Y < 0.) this->id_inverse = -this->id_inverse;
			if (this->param[0]) { //...реализация функций в негармоническом случае;
				kap  = 1./( q = this->param[0]);
				pim[0] = this->p[0] = cos(q*Y)*dd; 
				pim[1] = this->p[1] = sin(q*Y)*kap*dd; 
				kap *= kap;

				for (i2 = i = 0; i < this->N; i++, i2 = i*2) {
					pim[i2+2] = this->p[i2+2] = X*pim[i2]-i*Y*pim[i2-1];
					pim[i2+3] = this->p[i2+3] = X*pim[i2+1]+i*kap*(Y*pim[i2-2]-i*pim[i2-1]+(i-1)*X*pim[i2-3]);
				}
			}
			else { //...реализация функций в гармоническом случае;
				pim[0] = this->p[0] = dd; 
				pim[1] = this->p[1] = Y*dd; 

				for (i2 = i = 0; i < this->N; i++, i2 = i*2) {
					pim[i2+2] = this->p[i2+2] =  X*pim[i2]-i*Y*pim[i2-1];
					pim[i2+3] = this->p[i2+3] = (Y*pim[i2+2]+(i+1)*X*pim[i2+1])/(i+2.);
				}
			}
			for (fact = 1., m = 0; m <= mmm; m++, fact *= (m2+1)*(m2+2))
			for (q = fact*this->id_inverse, i2 = (i = (m2 = m*2))*2;  i <= this->N; i++, i2 = i*2) {
				this->p[i2] += q*pim[i2-2*m2+1]*imag(E1[m2]);
				this->p[i2+1] += q*pim[i2-2*m2]*imag(E1[m2+1]);
				q *= (i+1.)/(i-m2+1.);
			}
		}
	}
}

template <typename T>
void CHeat2DCircle<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p	 = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
			pim = (T *)new_struct( this->NN_dop*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(E1); E1 = (complex *)new_struct((this->N+1+m_dop)*sizeof(complex));
			double q =  -sqr(this->param[0]*this->param[1])*.25, 
					 g = q*sqr(this->param[2]), f, d, hh, A, B, C, D, E, F;
			int i, m;
			for (m = 0; m <= this->N+m_dop; m++) {
				B = (A = (f = 1.))*m; i = 1; 
				do {
					A += (f *= q/(i*(m+i)));
					B +=  f *(m+2*i); ++i;
				}
				while (fabs(f) > EE);

				F = (E = (f = 1.))*m; i = 1; 
				do {
					E += (f *= g/(i*(m+i)));
					F +=  f *(m+2*i); ++i;
				}
				while (fabs(f) > EE);

				for (D = -(C = (d = 1.))*m, hh = 0., i = 1; i < m; i++) {
					C += (d *= q/(i*(-m+i))); 
					D +=  d *(-m+2*i); hh  -= 1./i;
				}
				if (! m) C = (D = 2.)*0.; 
				else	{
					C += (d *= q/i)*(hh -= 1./i);
					D +=  d *(2.+(-m+2*i)*hh); ++i;
				}
				do {
					C += (d *= q/(i*(-m+i)))*(hh -= 1./(i-m)+1./i); 
					D +=  d *(2.+(-m+2*i)*hh); ++i;
				}
				while (fabs(d) > EE);

				if ((hh = E*D-this->param[3]*F*C) != 0.)
				E1[m] = comp(A*D-B*C, this->param[3]*F*A-E*B)/hh;
			}
	}
	if (P && this->p) {
		complex z = comp(P[0], P[1])*this->R0_inv, w = conj(this->param[1]/z);
		T		 f, sum, zm, zm_i, hh;
		double rr = norm(z), q = -rr*.25*sqr(this->param[0]), dd = exp(-P[6]*this->param[4]);
		int i, m;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		if (! (this->id_inverse = (rr <= sqr(this->param[1])))) //...internal system;
		for (zm = dd, zm_i = 0., q *= sqr(this->param[2]), m = 0; m <= this->N+m_dop; 
			  zm = (hh = zm)*real(z)-zm_i*imag(z), zm_i = hh*imag(z)+zm_i*real(z), m++) {
			sum  = (f = 1.); i = 1; 
			do {
				sum += (f *= q/(i*(m+i))); ++i;
			}
			while (fabs(f) > EE);

			if (m) this->p[m*2-1] = zm_i*sum;
			this->p[m*2] = zm*sum;
		}
		else {//...external system;
			for (zm = dd, zm_i = 0., m = 0; m <= this->N+m_dop; 
				  zm = (hh = zm)*real(z)-zm_i*imag(z), zm_i = hh*imag(z)+zm_i*real(z), m++) {//...asymptotic;
				sum  = (f = 1.); i = 1; 
				do {
					sum += (f *= q/(i*(m+i))); ++i;
				}
				while (fabs(f) > EE);

				if (m) this->p[m*2-1] = zm_i*sum;
				this->p[m*2] = zm*sum;
			}
			for (zm = dd, zm_i = 0., rr = -log(norm(w)), w *= this->param[1], m = 0; m <= this->N+m_dop; 
				  zm = (hh = zm)*real(z)-zm_i*imag(z), zm_i = hh*imag(z)+zm_i*real(z), m++) {//...extern system;
				for (sum = (f = 1.), hh = rr, i = 1; i < m; i++) {
					sum += (f *= q/(i*(-m+i))); hh  -= 1./i;
				}
				if (! m) sum = rr; 
				else	{
					sum += (f *= q/i)*(hh -= 1./i); ++i;
				}
				do {
					sum += (f *= q/(i*(-m+i)))*(hh -= 1./(i-m)+1./i); ++i;
				}
				while (fabs(f) > EE);

				if (m) pim[m*2-1] = zm_i*sum;
				pim[m*2] = zm*sum;
			}
		}
	}
}

template <typename T>
void CHeat3DSphere<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p	 = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
			pim = (T *)new_struct( this->NN_dop*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(E1); E1 = (complex *)new_struct((this->N+1+m_dop)*sizeof(complex));
			double q =  -sqr(this->param[0]*this->param[1])*.25, 
					 g = q*sqr(this->param[2]), f, d, hh, A, B, C, D, E, F;
			int i, m;
			for (m = 0; m <= this->N+m_dop; m++) {
				B = (A = (f = 1.))*m; i = 1; 
				do {
					A += (f *= q/(i*(m+.5+i)));
					B +=  f *(m+2*i); ++i;
				}
				while (fabs(f) > EE);

				F = (E = (f = 1.))*m; i = 1; 
				do {
					E += (f *= g/(i*(m+.5+i)));
					F +=  f *(m+2*i); ++i;
				}
				while (fabs(f) > EE);

				D = -(C = (d = 1.))*(m+1.); i = 1;
				do {
					C += (d *= q/(i*(-m-.5+i))); 
					D +=  d *(-m-1.+2*i); ++i;
				}
				while (i <= m || fabs(d) > EE);

				if ((hh = E*D-this->param[3]*F*C) != 0.)
				E1[m] = comp(A*D-B*C, this->param[3]*F*A-E*B)/hh;
			}
	}
	if (P && this->p) {
		complex z = comp(P[0], P[1])*this->R0_inv;
		T		  f, sum, FF, FF_i, hh;
		double  Z = P[2]*this->R0_inv, rr = norm(z)+Z*Z, q = -rr*.25*sqr(this->param[0]), rad, dd = exp(-P[6]*this->param[4]);
		int i, ii, i1, k, m;
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		for (this->p[0] = dd, i = 0; i < this->N+m_dop; i++) {//...насчитываем полиномы;
			this->p[i1 = (i+1)*(i+1)] = ((2.*i+1.)*Z*this->p[ii = i*i]-i*rr*this->p[(i-1)*(i-1)])/(i+1.); 
			this->p[i1+1] = (FF_i = (i*Z*this->p[ii+1]+2.*this->p[ii]*imag(z))*(hh = 1./(i+2.)));
			this->p[i1+2] = (FF = (i*Z*this->p[ii+2]+2.*this->p[ii]*real(z))*hh);
			for (m = 1; m <= i; m++) {
				this->p[i1+1+m*2] = (FF_i = ((i-m)*Z*this->p[ii+1+m*2]+2.*(m+1.)*(this->p[ii+m*2]*imag(z)+this->p[ii-1+m*2]*real(z)))*(hh = 1./(i+m+2.)));
				this->p[i1+2+m*2] = (FF = ((i-m)*Z*this->p[ii+2+m*2]+2.*(m+1.)*(this->p[ii+m*2]*real(z)-this->p[ii-1+m*2]*imag(z)))*hh);
			}
		}
		if (! (this->id_inverse = (rr <= sqr(this->param[1])))) //...internal system;
		for (q *= sqr(this->param[2]), i = 0; i <= this->N+m_dop; i++) {//...насчитываем множители;
			sum  = (f = 1.); k = 1; i1 = 2*i; ii = i*i;
			do {
				sum += (f *= q/(k*(i+.5+k))); ++k;
			}
			while (fabs(f) > EE);
			for (m = 0; m <= i1; m++) this->p[ii+m] *= sum;
		}
		else {//...external system;
			for (i = 0; i <= this->N+m_dop; i++) {//...asymptotic;
				sum  = (f = 1.); k = 1; i1 = 2*i; ii = i*i;
				do {
					sum += (f *= q/(k*(i+.5+k))); ++k;
				}
				while (fabs(f) > EE);
				for (m = 0; m <= i1; m++) this->p[ii+m] *= sum;
			}

			rad = this->param[1]/sqrt(rr);
			Z  *=(rr = rad*rad);
			z  *= rr; 
			rr *= sqr(this->param[1]);
			for (pim[0] = dd*rad, i = 0; i < this->N+m_dop; i++) {//...насчитываем инверсированные полиномы;
				pim[i1 = (i+1)*(i+1)] = ((2.*i+1.)*Z*pim[ii = i*i]-i*rr*pim[(i-1)*(i-1)])/(i+1.); 
				pim[i1+1] = (FF_i = (i*Z*pim[ii+1]+2.*pim[ii]*imag(z))*(hh = 1./(i+2.)));
				pim[i1+2] = (FF = (i*Z*pim[ii+2]+2.*pim[ii]*real(z))*hh);
				for (m = 1; m <= i; m++) {
					pim[i1+1+m*2] = (FF_i = ((i-m)*Z*pim[ii+1+m*2]+2.*(m+1.)*(pim[ii+m*2]*imag(z)+pim[ii-1+m*2]*real(z)))*(hh = 1./(i+m+2.)));
					pim[i1+2+m*2] = (FF = ((i-m)*Z*pim[ii+2+m*2]+2.*(m+1.)*(pim[ii+m*2]*real(z)-pim[ii-1+m*2]*imag(z)))*hh);
				}
			}
			for (i = 0; i <= this->N+m_dop; i++) {//...extern system;
				sum = (f = 1.); k = 1; i1 = 2*i; ii = i*i;
				do {
					sum += (f *= q/(k*(-i-.5+k))); ++k;
				}
				while (k <= m || fabs(f) > EE);
				for (m = 0; m <= i1; m++) pim[ii+m] *= sum;
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CHeat2DStreep<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! this->px && ! this->py) {
		this->px = (T *)new_struct((this->NN_grad = this->NN_dop)*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		deriv_X(this->px);
		deriv_Y(this->py);
	}
}

template <typename T>
void CHeat2DCircle<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m_dop))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		this->N += m_dop;
		deriv_X(this->px);
		deriv_Y(this->py);
		this->N -= m_dop;
	}
}

template <typename T>
void CHeat3DSphere<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m_dop))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
		this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py && this->pz) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		this->N += m_dop;
		deriv_X(this->px);
		deriv_Y(this->py);
		deriv_Z(this->pz);
		this->N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CHeat2DStreep<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = (T *)new_struct((this->NN_hess = this->NN_grad)*sizeof(T));
		this->pxy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyy = (T *)new_struct( this->NN_hess*sizeof(T));
	}
	if (P && this->pxx && this->pxy && this->pyy) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
	}
}

template <typename T>
void CHeat2DCircle<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = (T *)new_struct((this->NN_hess = freedom(this->N+m_dop))*sizeof(T));
		this->pxy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyy = (T *)new_struct( this->NN_hess*sizeof(T));
	}
	if (P && this->pxx && this->pxy && this->pyy) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		this->N += m_dop;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		this->N -= m_dop;
	}
}

template <typename T>
void CHeat3DSphere<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = (T *)new_struct((this->NN_hess = freedom(this->N+m_dop))*sizeof(T));
		this->pxy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyy = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pxz = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pyz = (T *)new_struct( this->NN_hess*sizeof(T));
		this->pzz = (T *)new_struct( this->NN_hess*sizeof(T));
	}
	if (P && this->pxx && this->pxy && this->pyy && this->pxz && this->pyz && this->pzz) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		memset (this->pxz, 0, this->NN_hess*sizeof(T));
		memset (this->pyz, 0, this->NN_hess*sizeof(T));
		memset (this->pzz, 0, this->NN_hess*sizeof(T));
		this->N += m_dop;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		swap(this->p, this->pz);
			deriv_X(this->pxz);
			deriv_Y(this->pyz);
			deriv_Z(this->pzz);
		swap(this->p, this->pz);
		this->N -= m_dop;
	}
}

///////////////////////////////////////////////////////////////////
//...copy of the composite multipoles for the cell with inclusions;
template <typename T>
void CHeat2DCircle<T>::cpy(T * p_ext)
{
	if (p_ext && this->p && E1) {
		if (! this->id_inverse) //...internal system;
		for (int i = 0; i < this->NN; i++) p_ext[i] = this->p[i]*real(E1[(i+1)/2]);
		else if (pim) //...external system;
		for (int i = 0; i < this->NN; i++) p_ext[i] = this->p[i]+pim[i]*imag(E1[(i+1)/2]);
	}
}

template <typename T>
void CHeat3DSphere<T>::cpy(T * p_ext)
{
	if (p_ext && this->p && E1) {
		if (! this->id_inverse) //...internal system;
		for (int i = this->N;				 i >= 0; i--)
		for (int m = i*2, ii = i*i; m >= 0; m--) p_ext[ii+m] = this->p[ii+m]*real(E1[i]);
		else if (pim) //...external system;
		for (int i = this->N;				 i >= 0; i--)
		for (int m = i*2, ii = i*i; m >= 0; m--) p_ext[ii+m] = this->p[ii+m]+pim[ii+m]*imag(E1[i]);
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat2DStreep<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int i = this->N; i > 0; i--) {
		deriv[i*2+1] += i*this->p[i*2-1]*f;
		deriv[i*2]	 += i*this->p[i*2-2]*f;
	}
}

template <typename T>
void CHeat2DCircle<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	if (! this->id_inverse) { //...internal system;
		double kappa = sqr(this->param[0]*this->param[2]);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-3]-.25/(m+1)*kappa*this->p[m*2+1])*real(E1[m])*f;
			 deriv[m*2]   += (m*this->p[m*2-2]-.25/(m+1)*kappa*this->p[m*2+2])*real(E1[m])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] -= .125*kappa*this->p[m*2+1]*real(E1[1])*f;
			 deriv[m*2]   += (m*this->p[m*2-2]-.125*kappa*this->p[m*2+2])*real(E1[1])*f;
		}
		deriv[0] -= .5*kappa*this->p[2]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		double kappa = sqr(this->param[0]), g0 = sqr(this->param[1]);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-3]-.25/(m+1)*kappa*this->p[m*2+1])*f;
			 deriv[m*2]   += (m*this->p[m*2-2]-.25/(m+1)*kappa*this->p[m*2+2])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] -= .125*kappa*this->p[m*2+1]*f;
			 deriv[m*2]   += (m*this->p[m*2-2]-.125*kappa*this->p[m*2+2])*f;
		}
		deriv[0] -= .5*kappa*this->p[2]*f;

		f /= g0; kappa *= sqr(g0);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] -= (m*pim[m*2+1]-.25/(m-1)*kappa*pim[m*2-3])*imag(E1[m])*f;
			 deriv[m*2]   -= (m*pim[m*2+2]-.25/(m-1)*kappa*pim[m*2-2])*imag(E1[m])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] -=  pim[m*2+1]*imag(E1[1])*f;
			 deriv[m*2]   -= (pim[m*2+2]+.25*kappa*pim[m*2-2])*imag(E1[1])*f;
		}
		deriv[0] += 2.*pim[2]*imag(E1[0])*f;
	}
}

template <typename T>
void CHeat3DSphere<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int i, m, ii = sqr(this->N), jj, ll = sqr(this->N+1);
	if (! this->id_inverse) { //...internal system;
		double kappa = sqr(this->param[0]*this->param[2]), ff, dd;
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i+3.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-(ff =.25/(m+1.)*(i-m-1.)*(i-m))*this->p[jj+m*2+1]+
									  d*(m*this->p[ll+m*2-3]-(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+1]))*real(E1[i])*f;
				deriv[ii+m*2] += (m*this->p[jj+m*2-2]-ff*this->p[jj+m*2+2]+ d*(m*this->p[ll+m*2-2]-dd*this->p[ll+m*2+2]))*real(E1[i])*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] -= ((ff =.25/(m+1.)*(i-m)*(i-m-1))*this->p[jj+m*2+1]+d*(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+1])*real(E1[i])*f;
				deriv[ii+m*2] += (m*this->p[jj+m*2-2]-ff*this->p[jj+m*2+2]+d*(m*this->p[ll+m*2-2]+dd*this->p[ll+m*2+2]))*real(E1[i])*f;
			}
			deriv[ii] -= (i*(i-1)*this->p[jj+2]+d*(i+1.)*(i+2.)*this->p[ll+2])*.5*real(E1[i])*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll+2]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		double kappa = sqr(this->param[0]), g0 = sqr(this->param[1]), ff, dd;
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i+3.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-(ff =.25/(m+1.)*(i-m-1.)*(i-m))*this->p[jj+m*2+1]+
									  d*(m*this->p[ll+m*2-3]-(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+1]))*f;
				deriv[ii+m*2] += (m*this->p[jj+m*2-2]-ff*this->p[jj+m*2+2]+ d*(m*this->p[ll+m*2-2]-dd*this->p[ll+m*2+2]))*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] -= ((ff =.25/(m+1.)*(i-m)*(i-m-1))*this->p[jj+m*2+1]+d*(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+1])*f;
				deriv[ii+m*2] += (m*this->p[jj+m*2-2]-ff*this->p[jj+m*2+2]+d*(m*this->p[ll+m*2-2]+dd*this->p[ll+m*2+2]))*f;
			}
			deriv[ii] -= (i*(i-1.)*this->p[jj+2]+d*(i+1.)*(i+2.)*this->p[ll+2])*.5*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll+2]*f;

		f /= g0; kappa *= sqr(g0);
		for (ll = sqr(this->N+1), ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i-1.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*pim[ll+m*2-3]-(ff =.25/(m+1.)*(i+m+1.)*(i+m+2))*pim[ll+m*2+1]+
									  d*(m*pim[jj+m*2-3]-(dd =.25/(m+1.)*(i-m-1.)*(i-m))*pim[jj+m*2+1]))*imag(E1[i])*f;
				deriv[ii+m*2] += (m*pim[ll+m*2-2]-ff*pim[ll+m*2+2]+ d*(m*pim[jj+m*2-2]-dd*pim[jj+m*2+2]))*imag(E1[i])*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] -= ((ff =.25/(m+1.)*(i+m+1.)*(i+m+2))*pim[ll+m*2+1]+d*(dd =.25/(m+1.)*(i-m-1.)*(i-m))*pim[jj+m*2+1])*imag(E1[i])*f;
				deriv[ii+m*2] += (m*pim[ll+m*2-2]-ff*pim[ll+m*2+2]+d*(m*pim[jj+m*2-2]+dd*pim[jj+m*2+2]))*imag(E1[i])*f;
			}
			deriv[ii] -= ((i+1.)*(i+2.)*pim[ll+2]+d*i*(i-1.)*pim[jj+2])*.5*imag(E1[i])*f;
		}
		deriv[ii] -= pim[ll+2]*imag(E1[0])*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat2DStreep<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! pim) return;
	f *= this->R0_inv;
	int m, i, i2, m2, mmm = this->N/2;
	if (! this->id_inverse) { //...internal system;
		double q,  kap = sqr(this->param[0]*this->param[2]);
		for (m = 0; m <= mmm; m++, f *= (m2+1)*(m2+2))
		for (q = f, i2 = (i = (m2 = m*2))*2;  i <= this->N; i++, i2 = i*2) {
			deriv[i2] -= q*pim[i2-2*m2+1]*(real(E1[m2])*kap+(m2 ? real(E1[m2-2]) : 0.));
			deriv[i2+1] += q*pim[i2-2*m2]* real(E1[m2+1]);
			q *= (i+1.)/(i-m2+1.);
		}
	}
	else { //...external system;
		double q,  kap = sqr(this->param[0]);
		for (i = this->N; i >= 0; i--) {
			deriv[i*2+1] += pim[i*2]*f;
			deriv[i*2]	 -= pim[i*2+1]*kap*f;
			if (i > 1) deriv[i*2] -= i*(i-1)*pim[i*2-3]*f;
		}
		for (m = 0; m <= mmm; m++, f *= (m2+1)*(m2+2))
		for (q = f*this->id_inverse, i2 = (i = (m2 = m*2))*2;  i <= this->N; i++, i2 = i*2) {
			deriv[i2] += q*pim[i2-2*m2]* imag(E1[m2]);
			deriv[i2+1] -= q*pim[i2-2*m2+1]*(imag(E1[m2+1])*kap+(m2 ? imag(E1[m2-1]) : 0.));
			q *= (i+1.)/(i-m2+1.);
		}
	}
}

template <typename T>
void CHeat2DCircle<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	if (! this->id_inverse) { //...internal system;
		double kappa = sqr(this->param[0]*this->param[2]);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-2]+.25/(m+1)*kappa*this->p[m*2+2])*real(E1[m])*f;
			 deriv[m*2]   -= (m*this->p[m*2-3]+.25/(m+1)*kappa*this->p[m*2+1])*real(E1[m])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += (m*this->p[m*2-2]+.125*kappa*this->p[m*2+2])*real(E1[1])*f;
			 deriv[m*2]   -= .125*kappa*this->p[m*2+1]*real(E1[1])*f;
		}
		deriv[0] -= .5*kappa*this->p[1]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		double kappa = sqr(this->param[0]), g0 = sqr(this->param[1]);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-2]+.25/(m+1)*kappa*this->p[m*2+2])*f;
			 deriv[m*2]   -= (m*this->p[m*2-3]+.25/(m+1)*kappa*this->p[m*2+1])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += (m*this->p[m*2-2]+.125*kappa*this->p[m*2+2])*f;
			 deriv[m*2]   -= .125*kappa*this->p[m*2+1]*f;
		}
		deriv[0] -= .5*kappa*this->p[1]*f;

		f /= g0; kappa *= sqr(g0);
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*pim[m*2+2]+.25/(m-1)*kappa*pim[m*2-2])*imag(E1[m])*f;
			 deriv[m*2]   -= (m*pim[m*2+1]+.25/(m-1)*kappa*pim[m*2-3])*imag(E1[m])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += (pim[m*2+2]-.25*kappa*pim[m*2-2])*imag(E1[1])*f;
			 deriv[m*2]   -=  pim[m*2+1]*imag(E1[1])*f;
		}
		deriv[0] += 2.*pim[1]*imag(E1[0])*f;
	}
}

template <typename T>
void CHeat3DSphere<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int i, m, ii = sqr(this->N), jj, ll = sqr(this->N+1);
	if (! this->id_inverse) { //...internal system;
		double kappa = sqr(this->param[0]*this->param[2]), ff, dd;
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i+3.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff =.25/(m+1.)*(i-m-1.)*(i-m))*this->p[jj+m*2+2]+
									  d*(m*this->p[ll+m*2-2]+(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+2]))*real(E1[i])*f;
				deriv[ii+m*2] -= (m*this->p[jj+m*2-3]+ff*this->p[jj+m*2+1]+ d*(m*this->p[ll+m*2-3]+dd*this->p[ll+m*2+1]))*real(E1[i])*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff =.25/(m+1.)*(i-m-1.)*(i-m))*this->p[jj+m*2+2]+
									  d*(m*this->p[ll+m*2-2]+(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+2]))*real(E1[i])*f;
				deriv[ii+m*2] -=(ff*this->p[jj+m*2+1]+d*dd*this->p[ll+m*2+1])*real(E1[i])*f;
			}
			deriv[ii] -= (i*(i-1)*this->p[jj+1]+d*(i+1.)*(i+2)*this->p[ll+1])*.5*real(E1[i])*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll+1]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		double kappa = sqr(this->param[0]), g0 = sqr(this->param[1]), ff, dd;
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i+3.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff =.25/(m+1.)*(i-m-1.)*(i-m))*this->p[jj+m*2+2]+
									  d*(m*this->p[ll+m*2-2]+(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+2]))*f;
				deriv[ii+m*2] -= (m*this->p[jj+m*2-3]+ff*this->p[jj+m*2+1]+ d*(m*this->p[ll+m*2-3]+dd*this->p[ll+m*2+1]))*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff =.25/(m+1.)*(i-m)*(i-m-1))*this->p[jj+m*2+2]+
									  d*(m*this->p[ll+m*2-2]+(dd =.25/(m+1.)*(i+m+1.)*(i+m+2))*this->p[ll+m*2+2]))*f;
				deriv[ii+m*2] -=(ff*this->p[jj+m*2+1]+d*dd*this->p[ll+m*2+1])*f;
			}
			deriv[ii] -= (i*(i-1)*this->p[jj+1]+d*(i+1.)*(i+2)*this->p[ll+1])*.5*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll+1]*f;

		f /= g0; kappa *= sqr(g0);
		for (ll = sqr(this->N+1), ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/((2.*i+1.)*(2.*i-1.));
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*pim[ll+m*2-2]+(ff =.25/(m+1.)*(i+m+1.)*(i+m+2))*pim[ll+m*2+2]+
									  d*(m*pim[jj+m*2-2]+(dd =.25/(m+1.)*(i-m-1.)*(i-m))*pim[jj+m*2+2]))*imag(E1[i])*f;
				deriv[ii+m*2] -= (m*pim[ll+m*2-3]+ff*pim[ll+m*2+1]+ d*(m*pim[jj+m*2-3]+dd*pim[jj+m*2+1]))*imag(E1[i])*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] += (m*pim[ll+m*2-2]+(ff =.25/(m+1.)*(i+m+1.)*(i+m+2.))*pim[ll+m*2+2]+
									  d*(m*pim[jj+m*2-2]+(dd =.25/(m+1.)*(i-m-1.)*(i-m))*pim[jj+m*2+2]))*imag(E1[i])*f;
				deriv[ii+m*2] -=(ff*pim[ll+m*2+1]+d*dd*pim[jj+m*2+1])*imag(E1[i])*f;
			}
			deriv[ii] -= ((i+1.)*(i+2)*pim[ll+1]+d*i*(i-1)*pim[jj+1])*.5*imag(E1[i])*f;
		}
		deriv[ii] -= pim[ll+1]*imag(E1[0])*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat3DSphere<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int i, m, ii = sqr(this->N), jj, ll = sqr(this->N+1);
	if (! this->id_inverse) { //...internal system;
		double kappa = sqr(this->param[0]*this->param[2]);
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/( (2.*i+1.)*(2.*i+3.));
			deriv[ii+i*2-1] -= (2.*i+1.)*d*this->p[ll+i*2-1]*real(E1[i])*f;
			deriv[ii+i*2]   -= (2.*i+1.)*d*this->p[ll+i*2]*real(E1[i])*f;
			for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
				deriv[ii+m*2-1] += ((i-m)*this->p[jj+m*2-1]-(i+m+1.)*d*this->p[ll+m*2-1])*real(E1[i])*f;
				deriv[ii+m*2]   += ((i-m)*this->p[jj+m*2]	 -(i+m+1.)*d*this->p[ll+m*2])*real(E1[i])*f;
			}
			deriv[ii] += (i*this->p[jj]-(i+1.)*d*this->p[ll])*real(E1[i])*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		double kappa = sqr(this->param[0]), g0 = sqr(this->param[1]);
		for ( i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/( (2.*i+1.)*(2.*i+3.));
			deriv[ii+i*2-1] -= (2.*i+1.)*d*this->p[ll+i*2-1]*f;
			deriv[ii+i*2]   -= (2.*i+1.)*d*this->p[ll+i*2]*f;
			for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
				deriv[ii+m*2-1] += ((i-m)*this->p[jj+m*2-1]-(i+m+1.)*d*this->p[ll+m*2-1])*f;
				deriv[ii+m*2]   += ((i-m)*this->p[jj+m*2]	 -(i+m+1.)*d*this->p[ll+m*2])*f;
			}
			deriv[ii] += (i*this->p[jj]-(i+1.)*d*this->p[ll])*f;
		}
		deriv[ii] -= kappa/3.*this->p[ll]*f;

		f /= g0; kappa *= sqr(g0);
		for (ll = sqr(this->N+1), ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double d = kappa/( (2.*i+1.)*(2.*i-1.));
			deriv[ii+i*2-1] -= (2.*i+1.)*pim[ll+i*2-1]*imag(E1[i])*f;
			deriv[ii+i*2]   -= (2.*i+1.)*pim[ll+i*2]*imag(E1[i])*f;
			for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
				deriv[ii+m*2-1] -= ((i+m+1.)*pim[ll+m*2-1]-(i-m)*d*pim[jj+m*2-1])*imag(E1[i])*f;
				deriv[ii+m*2]   -= ((i+m+1.)*pim[ll+m*2]  -(i-m)*d*pim[jj+m*2])*imag(E1[i])*f;
			}
			deriv[ii] -= ((i+1.)*pim[ll]-i*d*pim[jj])*imag(E1[i])*f;
		}
		deriv[ii] -= pim[ll]*imag(E1[0])*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CHeat2DStreep<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= -this->param[4];
	for (int i = this->N; i >= 0; i--) {
		deriv[i*2+1] += this->p[i*2+1]*f;
		deriv[i*2]	 += this->p[i*2]*f;
	}
}

template <typename T>
void CHeat2DCircle<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= -this->param[4];

	if (! this->id_inverse) { //...internal system;
		for (int m = this->N; m > 0; m--) {
			 deriv[m*2-1] += this->p[m*2-1]*real(E1[m])*f;
			 deriv[m*2]   += this->p[m*2]*real(E1[m])*f;
		}
		deriv[0] += this->p[0]*real(E1[0])*f;
	}
	else if (pim) { //...external system;
		for (int m = this->N; m > 0; m--) {
			 deriv[m*2-1] += (this->p[m*2-1]+pim[m*2-1]*imag(E1[m]))*f;
			 deriv[m*2]   += (this->p[m*2]+pim[m*2]*imag(E1[m]))*f;
		}
		deriv[0] += (this->p[0]+pim[0]*imag(E1[0]))*f;
	}
}

template <typename T>
void CHeat3DSphere<T>::deriv_t(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= -this->param[4];

	if (! this->id_inverse) //...internal system;
	for (int i = this->N;				 i >= 0; i--)
	for (int m = i*2, ii = i*i; m >= 0; m--) deriv[ii+m] += this->p[ii+m]*real(E1[i])*f;
	else if (pim) //...external system;
	for (int i = this->N;				 i >= 0; i--)
	for (int m = i*2, ii = i*i; m >= 0; m--) deriv[ii+m] += (this->p[ii+m]+pim[ii+m]*imag(E1[i]))*f;
}
#endif
