/*========================================*/
/*                 CEXPP                  */
/*========================================*/
#ifndef ___CEXPP___
#define ___CEXPP___

#include "cshapes.h"

////////////////////////////////////////////////////
//														        //
//  EXPONENTIAL SYSTEM of FUNDAMENTAL MULTIPOLES  //
//														        //
////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//...class of skin multipoles for interior and exterior fundamental solution;
template <typename T>
class CSkin3DExpp : public CShape<T> {
protected:
		T * sk;
public:
		int type()          { return SK3D_EXPP_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(2);}
public:
//...initialization of multipoles;
		void release() {
			 delete_struct(sk); 
			 CShape<T>::release();
		}
		void set_shape(double R0, double kk = 0., double p1 = 0., double p2 = 0., double p3 = 0., double p4 = 0.) {
			 CShape<T>::R0 = R0;
			 this->R0_inv  = R0 > EE ? 1./R0 : 1.;
			 int   k    = -1;
			 if (++k < size_of_param()) this->param[k] = (Param)kk; else return;
			 if (++k < size_of_param()) this->param[k] = (Param)p1; else return;
			 if (++k < size_of_param()) this->param[k] = (Param)p2; else return;
			 if (++k < size_of_param()) this->param[k] = (Param)p3; else return;
			 if (++k < size_of_param()) this->param[k] = (Param)p4; else return;
		}
public:
//...calculation of multipoles;
		void parametrization(double * P = NULL, int m_dop = 0);
//...integration;
		void primitive(T * deriv);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CSkin3DExpp (int inverse = 0) {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			this->id_inverse = inverse;
			sk = NULL;
       }
};

////////////////////////////////////////////////////////////////
//...class of regular fundamental solutions for sphere interior;
template <typename T>		
class CSkin3DZoom : public CSkin3DExpp<T> {
public:
		int type()			 { return SK3D_ZOOM_SHAPE;}
		int freedom(int m) { return(sqr(m+1));}
public:
//...initialization and calculation of multipoles;
		void parametrization(double * P = NULL, int m = 0);
//...integration;
		void primitive(T * deriv);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(0);
#define SHAPE_normalization_rad         0. //...normalization radius;
#undef  SHAPE_normalization_rad            //...param(1);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CSkin3DExpp<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
///////////////////////
//...technical support;
			delete_struct(sk); sk = (T *)new_struct((this->N+m_dop+1)*sizeof(T));
	}
	if (P && this->p) {
		double  Z = P[2], f3 = Z*Z, rr, r2, r1, r0 = sqr(0.5*this->param[1]);
		T		  h1, h2, hh, F1, F2, F1_i, F2_i, zm, zm_i;
		complex z = comp(P[0], P[1]);
		int     i, m; 
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		this->cs[3] = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		rr	= sqrt(f3 += real(z)*real(z)+imag(z)*imag(z));

///////////////////////////////////
//...prepare of the skin functions;
		//if (this->id_inverse) {//...гиперболический синус;
		//	r0 = sqr(this->param[1]*this->param[0]);
		//	hh = exp( this->param[0]*rr)*.5; 
		//	h1 = exp(-this->param[0]*rr)*.5; 
		//	h2 = (hh+h1)/(r0*this->param[0]);
		//	h1 = (hh-h1)/rr; 
		//	this->cs[4] = (sqr(rr*this->param[0])+2.*(rr*this->param[0])+2.)/sqr(rr*this->param[0]);
		//	this->cs[5] = 2./(sqr(rr*this->param[0])*rr);
		//}
		//else {//...гиперболический косинус;
		//	r0 = sqr(this->param[1]/sqr(this->param[0]));
		//	hh = exp( this->param[0]*rr)*.5; 
		//	h1 = exp(-this->param[0]*rr)*.5; 
		//	h2 = (hh-h1)/(r0*param[0]);
		//	h1 = (hh+h1)/rr; 
		//	this->cs[4] = (sqr(rr*this->param[0])-2.*(rr*this->param[0])+2.)/sqr(rr*this->param[0]);
		//	this->cs[5] = 2./(sqr(rr*this->param[0])*rr);
		//}
		r1 = this->R0_inv/r0; //r1 = sqr(0.5*this->param[1]);
		Z *= r1;	f3 *= sqr(r1);
		z *= r1;
		r1	= r0/sqr(rr);
		r2	= r0*r1;
		if (this->id_inverse) {//...убывающая экспонента;
			hh = exp(-this->param[0]*(rr-this->param[1])); 
			h2 = -hh/(r0*this->param[0]);
			h1 = hh/rr; 
			this->cs[4] = (sqr(rr*this->param[0])+2.*(rr*this->param[0])+2.)/sqr(rr*this->param[0]);
			this->cs[5] = 2./(sqr(rr*this->param[0])*rr);
		}
		else {//...растущая экспонента;
			hh = exp(this->param[0]*(rr-this->param[1])); 
			h2 = hh/(r0*this->param[0]);
			h1 = hh/rr; 
			this->cs[4] = (sqr(rr*this->param[0])-2.*(rr*this->param[0])+2.)/sqr(rr*this->param[0]);
			this->cs[5] = 2./(sqr(rr*this->param[0])*rr);
		}
////////////////////////////////////////////////////////////
//...рабочая прогонка специальных функций в оба направления;
		//int M = 15;
		//T * sk1 = (T *)new_struct(M*sizeof(T)),
		//  * sk2 = (T *)new_struct(M*sizeof(T));

		//for (sk1[m = 0] = h1; m < M; sk1[++m] = h1) {
		//	  h2 = (sqr(this->param[0])*h2*r0-(2.*m+1.)*h1)*r0/sqr(rr); 
		//	  swap(h1, h2);
		//}
		//for ( sk2[m = M-1] = sk1[M-1], sk2[--m] = sk1[M-2], --m; m >= 0; --m)
		//		sk2[m] = (sk2[m+2]*sqr(rr)/r0+(2.*m+3.)*sk2[m+1])/(sqr(this->param[0])*r0); 
//...рабочая прогонка специальных функций в оба направления;
////////////////////////////////////////////////////////////
		for ( sk[m = 0] = h1; m < this->N+m_dop; sk[++m] = h1) {
			  h2 = sqr(this->param[0])*h2*r2-(2.*m+1.)*h1*r1; 
			  swap(h1, h2);
		}
///////////////////////////////////
//...calculation of the multipoles;
		for (F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N+m_dop; i++) {
			  F2 = ((2.*i+1.)*Z*F1-i*f3*F2)/(i+1.); this->p[i*i] = F1*sk[i]; 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1*sk[i],	m = 1; m <= this->N+m_dop; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N+m_dop; i++) {
					this->p[i*i+m*2-1] = F1_i*sk[i];
					this->p[i*i+m*2]   = F1*sk[i];
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i*sk[i];
			 this->p[i*i+m*2]   = F1*sk[i];
		}
	}
}
template <typename T>
void CSkin3DZoom<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
///////////////////////
//...technical support;
			delete_struct(this->sk); this->sk = (T *)new_struct((this->N+m_dop+2)*sizeof(T));
	}
	if (P && this->p) {
		double   Z = P[2]*this->R0_inv;
		T			h1, h2, hh, f3 = Z*Z, rr, r2, F1, F2, F1_i, F2_i, zm, zm_i;
		complex  z = comp(P[0], P[1])*this->R0_inv;
		int      i, m; 

		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		this->cs[3] = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		f3   += real(z)*real(z)+imag(z)*imag(z);
		rr		= sqrt(to_double(f3))*this->R0*this->param[0];
		r2		= rr*rr;

////////////////////////////////////////////////////////////////////////////////////////////////
//...prepare of the skin functions: sk[m] = Norm*(1/r*d/dr)**m [sh(kappa*r)/r]*(kappa)*(-2*m-1);
		hh = fabs(this->param[0]*this->param[1]);
		if ( to_double(hh) > 5e-4) hh = this->param[0]*this->param[1]/(1.-exp(to_double(-2.*this->param[0]*this->param[1]))); 
		else hh = .5+hh*hh*.25*(1.-hh*hh*.1);
		if (to_double(r2) > 16.) {
			h2 = exp(to_double( rr-this->param[0]*this->param[1]))*hh; 
			hh = exp(to_double(-rr-this->param[0]*this->param[1]))*hh; 
			h1 = (h2-hh)/rr; 
			h2 = (h2+hh);
			for ( this->sk[m = 0] = h1; m <= this->N+m_dop; this->sk[++m] = h1) {
				  hh = h2-(2.*m+1.)*h1; 
				  h2 = hh/r2; 
				  swap(h1, h2);
			}
		}
		else
		for ( h1 = ((h2 = hh) *= 2.), m = 0; m <= this->N+m_dop+1; m++) {
			this->sk[m] = h1;
			for (i = 1; (fabs(h1) > 1e-18) && i <= 14; i++) 
				this->sk[m] += (h1 *= r2/((2.*i)*(2.*(i+m)+1.)));
			h1 = (h2 /= (2.*m+3.));
		}

///////////////////////////////////
//...calculation of the multipoles;
		for (F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N+m_dop; i++) {
			  F2 = ((2.*i+1.)*Z*F1-i*f3*F2)/(i+1.); this->p[i*i] = F1*this->sk[i]; 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1*this->sk[i],	m = 1; m <= this->N+m_dop; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N+m_dop; i++) {
					this->p[i*i+m*2-1] = F1_i*this->sk[i];
					this->p[i*i+m*2]   = F1*this->sk[i];
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i*this->sk[i];
			 this->p[i*i+m*2]   = F1*this->sk[i];
		}
	}
}

////////////////////////////
//...integration multipoles;
template <typename T>
void CSkin3DExpp<T>::primitive(T * deriv)
{
	if (! deriv || ! this->p) return;
	int  jj, ii, i, m; 
	for (jj = 0, ii = i = 1; i <= this->N; i++) {
		deriv[ii++] = (this->p[jj]*this->cs[2]-(i-1.)*.5*(this->p[jj+2]*this->cs[0]+this->p[jj+1]*this->cs[1]))*this->R0_inv; jj++;
		for (m = 1; m < i; m++) { 
			deriv[ii++] = (this->p[jj]*this->cs[2]-(i-m-1.)/(m+1.)*.5*(this->p[jj+2]*this->cs[0]-this->p[jj+3]*this->cs[1]))*this->R0_inv; jj++;
			deriv[ii++] = (this->p[jj]*this->cs[2]-(i-m-1.)/(m+1.)*.5*(this->p[jj+2]*this->cs[0]+this->p[jj+1]*this->cs[1]))*this->R0_inv; jj++;
		}
		if (1 < i) {
			deriv[ii++] = (this->p[jj-2]*this->cs[0]+this->p[jj-1]*this->cs[1])*this->R0_inv;
			deriv[ii++] = (this->p[jj-1]*this->cs[0]-this->p[jj-2]*this->cs[1])*this->R0_inv;
		}
		else {
			deriv[ii++] = this->p[jj-1]*this->cs[1]*this->R0_inv;
			deriv[ii++] = this->p[jj-1]*this->cs[0]*this->R0_inv;
		}
	}
////////////////////////
//...нулевой мультиполь;
	deriv[0] = (this->p[0]*(1.-this->cs[4])+this->cs[5])*this->cs[3]*.5;
}
template <typename T>
void CSkin3DZoom<T>::primitive(T * deriv)
{
	if (! deriv || ! this->p) return;
	double f0 = this->R0_inv/sqr(this->param[0]);
	int  jj, ii, i, m; 
	for (jj = 0, ii = i = 1; i <= this->N; i++) {
		deriv[ii++] = (this->p[jj]*this->cs[2]-(i-1.)*.5*(this->p[jj+2]*this->cs[0]+this->p[jj+1]*this->cs[1]))*f0; jj++;
		for (m = 1; m < i; m++) { 
			deriv[ii++] = (this->p[jj]*this->cs[2]-(i-m-1.)/(m+1.)*.5*(this->p[jj+2]*this->cs[0]-this->p[jj+3]*this->cs[1]))*f0; jj++;
			deriv[ii++] = (this->p[jj]*this->cs[2]-(i-m-1.)/(m+1.)*.5*(this->p[jj+2]*this->cs[0]+this->p[jj+1]*this->cs[1]))*f0; jj++;
		}
		if (1 < i) {
			deriv[ii++] = (this->p[jj-2]*this->cs[0]+this->p[jj-1]*this->cs[1])*f0;
			deriv[ii++] = (this->p[jj-1]*this->cs[0]-this->p[jj-2]*this->cs[1])*f0;
		}
		else {
			deriv[ii++] = this->p[jj-1]*this->cs[1]*f0;
			deriv[ii++] = this->p[jj-1]*this->cs[0]*f0;
		}
	}
////////////////////////
//...нулевой мультиполь;
	deriv[0] = (this->p[0]-this->sk[0]+2.*this->sk[1])*this->cs[3]*.5;
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CSkin3DExpp<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0_inv);
	f *= this->R0;

////////////////////////////////////
//...recurrent derivation (about X);
	deriv[0] += this->p[3]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.);
		deriv[ii++] += this->p[ll+2]*f; ll++;
		if (i > 1) {
			T ff = .5*i*(i-1.)*f0;
			deriv[ii-1] += ff*(this->p[ll+1]-kk*this->p[jj+2])*f; jj++;
			if (i > 2) {
				T ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] += (.5*this->p[ll+2]-ff*(kk*this->p[jj+2]-this->p[ll+2])*f0)*f; ll++;
				deriv[ii++] += (.5*this->p[ll+2]+((kk*this->p[jj-1]-this->p[ll-2])-ff*(kk*this->p[jj+3]-this->p[ll+2]))*f0)*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] += (.5*this->p[ll+2]+(m*(kk*this->p[jj]-this->p[ll-2])-ff*(kk*this->p[jj+4]-this->p[ll+2]))*f0)*f; jj++; ll++;
					deriv[ii++] += (.5*this->p[ll+2]+(m*(kk*this->p[jj]-this->p[ll-2])-ff*(kk*this->p[jj+4]-this->p[ll+2]))*f0)*f; jj++; ll++;
				}
				ff = (i-1.)*f0;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
				ff = i*f0;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
			}
			else {
				ff = (i-1.)*f0;
				deriv[ii++] += 0.5*this->p[ll+2]*f; ll++;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj-1]-this->p[ll-2]))*f; ll++;
				ff = i*f0;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+2]+ff*(kk*this->p[jj++]-this->p[ll-2]))*f; ll++;
			}
		}
		else {
			deriv[ii++] += 0.5*this->p[ll+2]*f; ll++;
			deriv[ii++] += (.5*this->p[ll+2]+(kk*this->p[jj++]-this->p[ll-2])*f0)*f; ll++;
		}
	}
}
template <typename T>
void CSkin3DZoom<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0);
	f *= this->R0_inv;

////////////////////////////////////
//...recurrent derivation (about X);
	deriv[0] += kk*this->p[3]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.);
		deriv[ii++] += kk*this->p[ll+2]*f; ll++;
		if (i > 1) {
			T ff = .5*i*(i-1.)*f0;
			deriv[ii-1] += ff*(kk*this->p[ll+1]-this->p[jj+2])*f; jj++;
			if (i > 2) {
				T ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] += (.5*kk*this->p[ll+2]-ff*(this->p[jj+2]-kk*this->p[ll+2])*f0)*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+2]+((this->p[jj-1]-kk*this->p[ll-2])-ff*(this->p[jj+3]-kk*this->p[ll+2]))*f0)*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] += (.5*kk*this->p[ll+2]+(m*(this->p[jj]-kk*this->p[ll-2])-ff*(this->p[jj+4]-kk*this->p[ll+2]))*f0)*f; jj++; ll++;
					deriv[ii++] += (.5*kk*this->p[ll+2]+(m*(this->p[jj]-kk*this->p[ll-2])-ff*(this->p[jj+4]-kk*this->p[ll+2]))*f0)*f; jj++; ll++;
				}
				ff = (i-1.)*f0;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
				ff = i*f0;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
			}
			else {
				ff = (i-1.)*f0;
				deriv[ii++] += 0.5*kk*this->p[ll+2]*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj-1]-kk*this->p[ll-2]))*f; ll++;
				ff = i*f0;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+2]+ff*(this->p[jj++]-kk*this->p[ll-2]))*f; ll++;
			}
		}
		else {
			deriv[ii++] += 0.5*kk*this->p[ll+2]*f; ll++;
			deriv[ii++] += (.5*kk*this->p[ll+2]+(this->p[jj++]-kk*this->p[ll-2])*f0)*f; ll++;
		}
	}
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CSkin3DExpp<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0_inv);
	f *= this->R0;

////////////////////////////////////
//...recurrent derivation (about Y);
	deriv[0] += this->p[2]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.);
		deriv[ii++] += this->p[ll+1]*f; ll++;
		if (i > 1) {
			T ff = .5*i*(i-1.)*f0;
			deriv[ii-1] += ff*(this->p[ll]-kk*this->p[jj+1])*f; jj++;
			if (i > 2) {
				T ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] -= (.5*this->p[ll+3]-((kk*this->p[jj-1]-this->p[ll-1])+ff*(kk*this->p[jj+3]-this->p[ll+3]))*f0)*f; ll++;
				deriv[ii++] += (.5*this->p[ll+1]-ff*(kk*this->p[jj+2]-this->p[ll+1])*f0)*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] -= (.5*this->p[ll+3]-(m*(kk*this->p[jj+1]-this->p[ll-1])+ff*(kk*this->p[jj+5]-this->p[ll+3]))*f0)*f; jj++; ll++;
					deriv[ii++] += (.5*this->p[ll+1]-(m*(kk*this->p[jj-1]-this->p[ll-3])+ff*(kk*this->p[jj+3]-this->p[ll+1]))*f0)*f; jj++; ll++;
				}
				ff = (i-1.)*f0;
				deriv[ii++] -= (.5*this->p[ll+3]-ff*(kk*this->p[(++jj)  ]-this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+1]-ff*(kk*this->p[(++jj)-2]-this->p[ll-3]))*f; ll++;
				ff = i*f0;
				deriv[ii++] -= (.5*this->p[ll+3]-ff*(kk*this->p[(++jj)  ]-this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+1]-ff*(kk*this->p[(++jj)-2]-this->p[ll-3]))*f; ll++;
			}
			else {
				ff = (i-1.)*f0;
				deriv[ii++] -= (.5*this->p[ll+3]-ff*(kk*this->p[jj-1]-this->p[ll-1]))*f; ll++;
				deriv[ii++] += 0.5*this->p[ll+1]*f; ll++;
				ff = i*f0;
				deriv[ii++] -= (.5*this->p[ll+3]-ff*(kk*this->p[(++jj)  ]-this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*this->p[ll+1]-ff*(kk*this->p[(++jj)-2]-this->p[ll-3]))*f; ll++;
			}
		}
		else {
			deriv[ii++] -= (.5*this->p[ll+3]-(kk*this->p[jj++]-this->p[ll-1])*f0)*f; ll++;
			deriv[ii++] += 0.5*this->p[ll+1]*f; ll++;
		}
	}
}
template <typename T>
void CSkin3DZoom<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0);
	f *= this->R0_inv;

////////////////////////////////////
//...recurrent derivation (about Y);
	deriv[0] += kk*this->p[2]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.);
		deriv[ii++] += kk*this->p[ll+1]*f; ll++;
		if (i > 1) {
			T ff = .5*i*(i-1.)*f0;
			deriv[ii-1] += ff*(kk*this->p[ll]-this->p[jj+1])*f; jj++;
			if (i > 2) {
				T ff = (i-2.)*(i-1.)*.125;
				deriv[ii++] -= (.5*kk*this->p[ll+3]-((this->p[jj-1]-kk*this->p[ll-1])+ff*(this->p[jj+3]-kk*this->p[ll+3]))*f0)*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+1]-ff*(this->p[jj+2]-kk*this->p[ll+1])*f0)*f; ll++;
				for (int m = 2; m < i-1; m++) {
					ff = (i-m-1.)*(i-m)*.25/(m+1.);
					deriv[ii++] -= (.5*kk*this->p[ll+3]-(m*(this->p[jj+1]-kk*this->p[ll-1])+ff*(this->p[jj+5]-kk*this->p[ll+3]))*f0)*f; jj++; ll++;
					deriv[ii++] += (.5*kk*this->p[ll+1]-(m*(this->p[jj-1]-kk*this->p[ll-3])+ff*(this->p[jj+3]-kk*this->p[ll+1]))*f0)*f; jj++; ll++;
				}
				ff = (i-1.)*f0;
				deriv[ii++] -= (.5*kk*this->p[ll+3]-ff*(this->p[(++jj)  ]-kk*this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+1]-ff*(this->p[(++jj)-2]-kk*this->p[ll-3]))*f; ll++;
				ff = i*f0;
				deriv[ii++] -= (.5*kk*this->p[ll+3]-ff*(this->p[(++jj)  ]-kk*this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+1]-ff*(this->p[(++jj)-2]-kk*this->p[ll-3]))*f; ll++;
			}
			else {
				ff = (i-1.)*f0;
				deriv[ii++] -= (.5*kk*this->p[ll+3]-ff*(this->p[jj-1]-kk*this->p[ll-1]))*f; ll++;
				deriv[ii++] += 0.5*kk*this->p[ll+1]*f; ll++;
				ff = i*f0;
				deriv[ii++] -= (.5*kk*this->p[ll+3]-ff*(this->p[(++jj)  ]-kk*this->p[ll-1]))*f; ll++;
				deriv[ii++] += (.5*kk*this->p[ll+1]-ff*(this->p[(++jj)-2]-kk*this->p[ll-3]))*f; ll++;
			}
		}
		else {
			deriv[ii++] -= (.5*kk*this->p[ll+3]-(this->p[jj++]-kk*this->p[ll-1])*f0)*f; ll++;
			deriv[ii++] += 0.5*kk*this->p[ll+1]*f; ll++;
		}
	}
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CSkin3DExpp<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0_inv);
	f *= this->R0;

////////////////////////////////////
//...recurrent derivation (about Z);
	deriv[0] += this->p[1]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.), ff = i*f0;
		deriv[ii++] += (this->p[ll]+ff*(kk*this->p[jj++]-this->p[ll]))*f; ll++;
		for (int m = 1; m < i; m++) {
			ff = (i-m)*f0;
			deriv[ii++] += (this->p[ll]+ff*(kk*this->p[jj++]-this->p[ll]))*f; ll++;
			deriv[ii++] += (this->p[ll]+ff*(kk*this->p[jj++]-this->p[ll]))*f; ll++;
		}
		deriv[ii++] += this->p[ll++]*f;
		deriv[ii++] += this->p[ll++]*f;
	}
}
template <typename T>
void CSkin3DZoom<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	double kk = sqr(this->param[0]*this->R0);
	f *= this->R0_inv;

////////////////////////////////////
//...recurrent derivation (about Z);
	deriv[0] += kk*this->p[1]*f;
	for ( int ii = 1, jj = 0, ll = 4, i = 1; i <= this->N; ll += 2, i++) {
		T f0 = 1./(2.*i+1.), ff = i*f0;
		deriv[ii++] += (kk*this->p[ll]+ff*(this->p[jj++]-kk*this->p[ll]))*f; ll++;
		for (int m = 1; m < i; m++) {
			ff = (i-m)*f0;
			deriv[ii++] += (kk*this->p[ll]+ff*(this->p[jj++]-kk*this->p[ll]))*f; ll++;
			deriv[ii++] += (kk*this->p[ll]+ff*(this->p[jj++]-kk*this->p[ll]))*f; ll++;
		}
		deriv[ii++] += kk*this->p[ll++]*f;
		deriv[ii++] += kk*this->p[ll++]*f;
	}
}

//////////////////////////////////////////////////
//													         //
//    BASIC SYSTEM of ELLIPSOIDAL MULTIPOLES    //
//														      //
//////////////////////////////////////////////////
//////////////////////////////////////////////
//...class of harmonic ellipsoidal multipoles;
template <typename T>
class CMapi3DEll : public CShape<T> {
public:
		int type() {return MP3D_ELLI_SHAPE;}
		int freedom(int m) {return (m+1)*(m+2);}
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CMapi3DEll () {
			this->param = (Param *)new_struct(this->size_of_param()*sizeof(Param));
		}
};

/////////////////////////////////////////////////////
//...class of skin multipoles, based on one exponent;
template <typename T>
class CSkin2DEll : public CShape<T> {
public:
		int type()          { return SK2D_ELLI_SHAPE;}
		int freedom (int m) { return (m+1)*2;}
		int size_of_param() { return(5);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
//...constructor;
		CSkin2DEll() {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized exponential number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized exponential number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized exponential number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized exponential number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_kappa_shift               0. //...shift-normalized exponential number;
#undef  SHAPE_kappa_shift                  //...param(4);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CMapi3DEll<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
		this->p = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
	}
	if (P) { //...индексация -- [n*(n+1)+m*2], [n*(n+1)+m*2+1];
		int m0(0), m1(2), m2, i, i0, i1, i2, m;
		double  X  = P[0]*this->R0_inv,
				  Y  = P[1]*this->R0_inv,
				  Z  = P[2]*this->R0_inv;
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];

//////////////////////////////////////////////////
//...симметричное и антисимметричное по Z решение;
		this->p[0] = 1; 
		this->p[1] = Z; 
		for (m = 0; m < this->N; m++, m0 = m1, m1 = m2) {
			for (i0 = (i1 = m0)-m*2, i = m; i < this->N; i++, i0 = i1, i1 = i2) {
				i2 = i1+2*(i+1);
				this->p[i2+m*2+1] = ((Z*(this->p[i2+m*2] = X*this->p[i1+m*2]-(i-m)*Z*this->p[i0+m*2+1])+(i-m+1.)*X*this->p[i1+m*2+1])-
								(m-1.)*m*(this->p[i2+m*2-3]-X*this->p[i1+m*2-3])/(i-m+2.))/(i-m+2.);
			}
			m2 = m1+2*(m+2);
			this->p[m2-1] = (Z*(this->p[m2-2] = Y*this->p[m1-2]-m*Z*this->p[m0-1])+(m+1.)*Y*this->p[m1-1])/(m+2.);
		}
	}
}

template <typename T>
void CSkin2DEll<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
		this->p = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
	}
	if (P) {
		int  i, i2;
		double  X  = P[0]*this->R0_inv,
				  Y  = P[1]*this->R0_inv, ff;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

//////////////////////////////////////////////
//...calculation direction of the integration;
		if (this->param[0] > EE_dop) { //...индексация -- [n*2], [n*2+1];

//////////////////////////////////////////////////
//...симметричное и антисимметричное по X решение;
			this->p[0] = (exp(this->param[0]*(X-this->param[4]*this->R0_inv))+(ff = exp(-this->param[0]*(X+this->param[4]*this->R0_inv))))*.5; 
			this->p[1] = (this->p[0]-ff)*this->param[1]; 
			for (i2 = i = 0; i < this->N; i++, i2 = i*2) {
				this->p[i2+2] = Y*this->p[i2]-i*X*this->p[i2-1];
				this->p[i2+3] = Y*this->p[i2+1]-i*sqr(this->param[1])*(X*this->p[i2-2]-i*this->p[i2-1]+(i-1)*Y*this->p[i2-3]);
			}
		}
		else {
///////////////////////////////////////////////////////////////
//...гармоничекое симметричное и антисимметричное по X решение;
			for ( this->p[0] = 1, this->p[1] = X, i2 = i = 0; i < this->N; i++, i2 = i*2) {
				this->p[i2+2] =  Y*this->p[i2]-i*X*this->p[i2-1];
				this->p[i2+3] = (X*this->p[i2+2]+(i+1)*Y*this->p[i2+1])/(i+2.);
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CMapi3DEll<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! this->px && ! this->py && ! this->pz) {
		this->px = (T *)new_struct((this->NN_grad = this->NN_dop)*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
		this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py && this->pz) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		deriv_X(this->px);
		deriv_Y(this->py);
		deriv_Z(this->pz);
	}
}

template <typename T>
void CSkin2DEll<T>::parametrization_grad(double * P, int m_dop)
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

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CMapi3DEll<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! this->pxx && ! this->pxy && ! this->pyy && ! this->pxz && ! this->pyz && ! this->pzz) {
		this->pxx = (T *)new_struct((this->NN_hess = this->NN_grad)*sizeof(T));
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
	}
}

template <typename T>
void CSkin2DEll<T>::parametrization_hess(double * P, int m_dop)
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

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi3DEll<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = this->N*(this->N+1), i = this->N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*i, m = i-1; m >= 0; m--) {
			deriv[ii+m*2+1] += (i-m)*this->p[jj+m*2+1]*f;
			deriv[ii+m*2]   += (i-m)*this->p[jj+m*2]*f;
		}
	}
}

template <typename T>
void CSkin2DEll<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int i = this->N; i >= 0; i--) {
		deriv[i*2+1] += this->p[i*2]*f;
		deriv[i*2]	 += this->p[i*2+1]*this->param[2]*f;
		if (i > 1) deriv[i*2] -= i*(i-1)*this->p[i*2-3]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi3DEll<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = this->N*(this->N+1), i = this->N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*i, m = i; m > 0; m--) {
			deriv[ii+m*2+1] += m*this->p[jj+m*2-1]*f;
			deriv[ii+m*2]   += m*this->p[jj+m*2-2]*f;
		}
	}
}

template <typename T>
void CSkin2DEll<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int i = this->N; i > 0; i--) {
		deriv[i*2+1] += i*this->p[i*2-1]*f;
		deriv[i*2]	 += i*this->p[i*2-2]*f;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi3DEll<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = this->N*(this->N+1), i = this->N; i >= 0; i--, ii = i*(i+1)) {
		for (jj = (i-2)*(i-1), m = i; m >= 0; m--) {
			deriv[ii+m*2+1] += this->p[ii+m*2]*f;
			if (m < i-1) deriv[ii+m*2] -= (i-m)*(i-m-1)*this->p[jj+m*2+1]*f;
			if (m > 1)	 deriv[ii+m*2] -= m*(m-1)*this->p[jj+m*2-3]*f;
		}
	}
}

//////////////////////////////////////////////////////////////
//																				//
//    BASIC SYSTEM of SKIN MULTIPOLES - EXPONENTIAL BEAM    //
//																				//
//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...class of skin multipoles of 2D beams type with exponntial along X;
template <typename T>
class CSkin2DBeamZ : public CShape<T> {
public:
		int type()          { return SK2D_BEAMZSHAPE;}
		int freedom (int m) { return((m+1)*2);}
		int size_of_param() { return(2);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
//...constructor;
		CSkin2DBeamZ () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number (directional);
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);

/////////////////////////////////////////////////////////////////////////////////
//...class of skin multipoles of beams type with two exponential characteristics;
template <typename T>
class CSkin3DBeam : public CShape<T> {
protected:
		T * sk, * sz, * pim, * pxim, * pyim, * pzim;
		T * E1, * E1_i;
public:
		int type()          { return SK3D_BEAM_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(6);};
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0., double kk_dop = 0., double L1 = 0., double L2 = 0.) {
			CSkin3DBeam<T>::R0 = R0;
			this->R0_inv     = R0 > EE ? 1./R0 : 1.;
			int   k    = -1;
			if (++k < size_of_param()) this->param[k] = (Param)((kk = sqrt(fabs(kk*kk-kk_dop*kk_dop)))*this->R0);
			if (++k < size_of_param()) this->param[k] = (Param)(this->param[0] > EE ? 1./this->param[0] : 1.);
			if (++k < size_of_param()) this->param[k] = (Param)(this->param[0]*this->param[0]);
			if (++k < size_of_param()) this->param[k] = (Param)kk;
			if (++k < size_of_param()) this->param[k] = (Param)((kk_dop = fabs(kk_dop))*this->R0);
			if (++k < size_of_param()) this->param[k] = (Param)kk_dop;
		}
		void release() {
			delete_struct(pim);
			delete_struct(pxim);
			delete_struct(pyim);
			delete_struct(pzim);
			delete_struct(sk);
			delete_struct(sz);
			delete_struct(E1);
			CShape<T>::release();
		}
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CSkin3DBeam () {
			delete_struct(this->param);
			this->param   = (Param *)new_struct(size_of_param()*sizeof(Param));
			sk = sz = pim = pxim = pyim = pzim = NULL;
			E1 = E1_i = NULL;
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized cylindrical wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_norm_kappa_dop            0. //...normalized directional wave number;
#undef  SHAPE_norm_kappa_dop               //...param(4);
#define SHAPE_kappa_dop                 0. //...non-normalized wave number;

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CSkin2DBeamZ<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
	}
	if (P && this->p) {
		T   X = P[0]*this->R0_inv, Y = P[1]*this->R0_inv,
			F1 = exp(to_double(-this->param[0]*X)), F2 = F1*Y, F0 = 0., F3 = F1, F4 = F2;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		this->p[0] = F1; this->p[1] = F2;
		if (this->param[0] > EE_dop)
		for (int m = 1; m <= this->N; m++) {
			F1 =  F1*X+(this->param[0]*F2-(m-1.)*F0)*Y;
			F0 = (F1*Y+(this->param[1]*.5*(F4-F3*Y)+F2*X*2.)*m)/(m*2.+1.); swap(F0, F2);
			F3 =  F1+F3*this->param[1]*.5*m;
			F4 =  F2+F4*this->param[1]*.5*m;
			this->p[m*2] = F1; this->p[m*2+1] = F2;
		}
		else
		for (int m = 1; m <= this->N; m++) {
			F1 =  F1*X-(m-1.)*F0*Y;
			F0 = (F1*Y+F2*X*m)/(m+1.); swap(F0, F2);
			this->p[m*2] = F1; this->p[m*2+1] = F2;
		}
	}
}

template <typename T>
void CSkin3DBeam<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p   = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
			pim = (T *)new_struct( this->NN_dop*sizeof(T));
///////////////////////
//...technical support;
			delete_struct(sk); sk = (T *)new_struct((this->N+1+m_dop)*sizeof(T));
			delete_struct(sz); sz = (T *)new_struct((this->N+1+m_dop)*sizeof(T));
			delete_struct(E1); E1 = (T *)new_struct((this->N+1+m_dop)*sizeof(T));
			delete_struct(E1_i); E1_i = (T *)new_struct((this->N+1+m_dop)*sizeof(T));
	}
	if (P && this->p && pim) {
		int  i, j, l, ii, nm, m;
		complex w = comp(P[0], P[1])*this->R0_inv;
		T	zm = 1., zm_i = 0., F, A, F_i, A_i, d, f, ff, hh;
		double rr = abs(w),
				Z   = P[2]*this->R0_inv,
				Co  = exp( this->param[4]*Z),
				Si  = exp(-this->param[4]*Z);
		if (rr > EE)
			w *= (1./rr); else w  = comp(1.);
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		for (sz[i = 0] = 1.; i < this->N+m_dop; i++) sz[i+1] = sz[i]*Z;
		if (sk) {
			T  z = this->param[0]*rr, h0, rr0 = z*z*.25;
			int m;
			for ( m = 0, h0 = 1.; m <= this->N+m_dop; m++, h0 *= rr) {
				T   f = 1., P = f; 
				int k = 1;
				do {
					 P += (f *= rr0/((T)k*(m+k))); ++k;
				}
				while (fabs(f) > EE);
				sk[m] = h0*P;
			}
		}
		for (i = 0; i <= this->N+m_dop; i++) {
			memset(E1, 0, (this->N+m_dop+1)*sizeof(T));
			memset(E1_i, 0, (this->N+m_dop+1)*sizeof(T));
			E1[nm = i] = 1.; E1_i[nm] = 1.;
			F = (ff = .5)*E1[nm]*sz[nm]*sk[0]; F_i = ff*E1_i[nm]*sz[nm]*sk[0];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/sqr(j);
				A = 0.; A_i = 0.;
				for (d = sz[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(-2.*this->param[4]*E1[l+1]-(l+2.)*E1[l+2])*f;
					E1_i[l] = (l+1.)*(2.*this->param[4]*E1_i[l+1]-(l+2.)*E1_i[l+2])*f;
					A += E1[l]*d; A_i += E1_i[l]*d;
					d *= Z;
				}
				E1[l] = (l+1.)*(-2.*this->param[4]*E1[l+1])*f;
				E1_i[l] = (l+1.)*(2.*this->param[4]*E1_i[l+1])*f;
				A += E1[l]*d; A_i += E1_i[l]*d;
				F += (ff *= rr)*A*sk[j]; F_i += ff*A_i*sk[j];
			}
			this->p  [ii] = F*Co+F_i*Si;
			pim[ii] = F*Co-F_i*Si;
		}
		for (m = 1; m <= this->N+m_dop; m++)
		for (zm = (hh = zm)*real(w)-zm_i*imag(w), zm_i = hh*imag(w)+zm_i*real(w), i = m; i <= this->N+m_dop; i++) {
			memset(E1, 0, (this->N+m_dop+1)*sizeof(T)); memset(E1_i, 0, (this->N+m_dop+1)*sizeof(T)); 

			E1[nm = i-m] = 1.; E1_i[nm] = 1.;
			F = (ff = .5)*E1[nm]*sz[nm]*sk[m]; F_i = ff*E1_i[nm]*sz[nm]*sk[m];
			for (ii = i*i, j = 1; j <= nm; j++) {
				f = .25/(j*(m+j));
				A = 0.; A_i = 0.;
				for (d = sz[l = max(0, nm-2*j)]; l <= nm-j-1; l++) {
					E1[l] = (l+1.)*(-2.*this->param[4]*E1[l+1]-(l+2.)*E1[l+2])*f;
					E1_i[l] = (l+1.)*(2.*this->param[4]*E1_i[l+1]-(l+2.)*E1_i[l+2])*f;
					A    += E1[l]*d;	
					A_i  += E1_i[l]*d;
					d    *= Z;
				}
				E1[l] = (l+1.)*(-2.*this->param[4]*E1[l+1])*f;
				E1_i[l] = (l+1.)*(2.*this->param[4]*E1_i[l+1])*f;
				A += E1[l]*d; A_i += E1_i[l]*d;
				F += (ff *= rr)*A*sk[m+j];	F_i += ff*A_i*sk[m+j];
			}
			this->p  [ii+m*2-1] = (F*Co+F_i*Si)*zm_i;
			this->p  [ii+m*2]   = (F*Co+F_i*Si)*zm;
			pim[ii+m*2-1] = (F*Co-F_i*Si)*zm_i;
			pim[ii+m*2]   = (F*Co-F_i*Si)*zm;
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CSkin2DBeamZ<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! P) {
		delete_struct(this->px);
		delete_struct(this->py);
	}
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
void CSkin3DBeam<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py && ! this->pz && ! pxim && ! pyim && ! pzim) {
		this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m_dop))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
		this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
		pxim = (T *)new_struct( this->NN_grad*sizeof(T));
		pyim = (T *)new_struct( this->NN_grad*sizeof(T));
		pzim = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py && this->pz && pxim && pyim && pzim) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		memset (pxim, 0, this->NN_grad*sizeof(T));
		memset (pyim, 0, this->NN_grad*sizeof(T));
		memset (pzim, 0, this->NN_grad*sizeof(T));
		this->N += m_dop;
		swap(this->p, pim); this->param[4] = -this->param[4];
			deriv_X(pxim);
			deriv_Y(pyim);
			deriv_Z(pzim);
		swap(this->p, pim); this->param[4] = -this->param[4];
			deriv_X(this->px);
			deriv_Y(this->py);
			deriv_Z(this->pz);
		this->N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CSkin2DBeamZ<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! P) {
		delete_struct(this->pxx);
		delete_struct(this->pxy);
		delete_struct(this->pyy);
	}
	if (! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = (T *)new_struct((this->NN_hess = this->NN_grad)*sizeof(T));
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
void CSkin3DBeam<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! this->pxx && ! this->pxy && ! this->pyy && ! this->pxz && ! this->pyz && ! this->pzz) {
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
		swap(this->p, this->px); swap(pim, pxim);
			deriv_X(this->pxx);
		swap(this->p, this->px); swap(pim, pxim);
		swap(this->p, this->py); swap(pim, pyim);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py); swap(pim, pyim);
		swap(this->p, this->pz); swap(pim, pzim);
			deriv_X(this->pxz);
			deriv_Y(this->pyz);
			deriv_Z(this->pzz);
		swap(this->p, this->pz); swap(pim, pzim);
		this->N -= m_dop;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DBeamZ<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m = this->N; m > 0; m--) {
		deriv[m*2+1] += (m*this->p[m*2-1]-this->param[0]*this->p[m*2+1])*f;
		deriv[m*2]   += (m*this->p[m*2-2]-this->param[0]*this->p[m*2])*f;
	}
	deriv[1] -= this->param[0]*this->p[1]*f;
	deriv[0] -= this->param[0]*this->p[0]*f;
}

template <typename T>
void CSkin3DBeam<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p || ! pim) return;
	f *= this->R0_inv;
	int ll = sqr(this->N+1);
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-.25/(m+1.)*((i-m)*((i-m-1)*this->p[jj+m*2+1]+this->param[4]*pim[ii+m*2+1]*2.)-this->param[2]*this->p[ll+m*2+1]))*f;
			deriv[ii+m*2]   += (m*this->p[jj+m*2-2]-.25/(m+1.)*((i-m)*((i-m-1)*this->p[jj+m*2+2]+this->param[4]*pim[ii+m*2+2]*2.)-this->param[2]*this->p[ll+m*2+2]))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1.)*this->p[jj+i*2-5]-.25/i*(this->param[4]*pim[ii+i*2-1]*2.-this->param[2]*this->p[ll+i*2-1]))*f;
			deriv[ii+i*2-2] += ((i-1.)*this->p[jj+i*2-4]-.25/i*(this->param[4]*pim[ii+i*2]*2.-this->param[2]*this->p[ll+i*2]))*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += (i*this->p[jj+i*2-3]+.25/(i+1.)*this->param[2]*this->p[ll+i*2+1])*f;
			deriv[ii+i*2]   += (i*this->p[jj+i*2-2]+.25/(i+1.)*this->param[2]*this->p[ll+i*2+2])*f;

			deriv[ii+1] -= (.125*((i-1)*((i-2)*this->p[jj+3]+this->param[4]*pim[ii+3]*2.)-this->param[2]*this->p[ll+3]))*f;
			deriv[ii+2] -= (.125*((i-1)*((i-2)*this->p[jj+4]+this->param[4]*pim[ii+4]*2.)-this->param[2]*this->p[ll+4])-this->p[jj])*f;

			deriv[ii] -= .5*(i*((i-1)*this->p[jj+2]+this->param[4]*pim[ii+2]*2.)-this->param[2]*this->p[ll+2])*f;
		}
		if (i == 1) {
			deriv[ii+1] += (.125*this->param[2]*this->p[ll+3])*f;
			deriv[ii+2] += (.125*this->param[2]*this->p[ll+4]+this->p[jj])*f;

			deriv[ii] -= .5*(this->param[4]*pim[ii+2]*2.-this->param[2]*this->p[ll+2])*f;
		}
	}
	deriv[0] += .5*this->param[2]*this->p[ll+2]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DBeamZ<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m = this->N; m > 0; m--) {
		deriv[m*2+1] += this->p[m*2]*f;
		deriv[m*2]   += m*this->p[m*2-1]*this->param[0]*2.*f;
		if (m > 1) deriv[m*2] -= m*(m-1)*this->p[m*2-3]*f;
	}
	deriv[1] += this->p[0]*f;
}

template <typename T>
void CSkin3DBeam<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p || ! pim) return;
	f *= this->R0_inv;
	int ll = sqr(this->N+1);
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-2; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+.25/(m+1.)*((i-m)*((i-m-1.)*this->p[jj+m*2+2]+this->param[4]*pim[ii+m*2+2]*2.)-this->param[2]*this->p[ll+m*2+2]))*f;
			deriv[ii+m*2]   -= (m*this->p[jj+m*2-3]+.25/(m+1.)*((i-m)*((i-m-1.)*this->p[jj+m*2+1]+this->param[4]*pim[ii+m*2+1]*2.)-this->param[2]*this->p[ll+m*2+1]))*f;
		}
		if (i > 2) {
			deriv[ii+i*2-3] += ((i-1)*this->p[jj+i*2-4]+.25/i*(this->param[4]*pim[ii+i*2]*2.-this->param[2]*this->p[ll+i*2]))*f;
			deriv[ii+i*2-2] -= ((i-1)*this->p[jj+i*2-5]+.25/i*(this->param[4]*pim[ii+i*2-1]*2.-this->param[2]*this->p[ll+i*2-1]))*f;
		}
		if (i > 1) {
			deriv[ii+i*2-1] += (i*this->p[jj+i*2-2]-.25/(i+1.)*this->param[2]*this->p[ll+i*2+2])*f;
			deriv[ii+i*2]   -= (i*this->p[jj+i*2-3]-.25/(i+1.)*this->param[2]*this->p[ll+i*2+1])*f;

			deriv[ii+1] += (.125*((i-1.)*((i-2.)*this->p[jj+4]+this->param[4]*pim[ii+4]*2.)-this->param[2]*this->p[ll+4])+this->p[jj])*f;
			deriv[ii+2] -= (.125*((i-1.)*((i-2.)*this->p[jj+3]+this->param[4]*pim[ii+3]*2.)-this->param[2]*this->p[ll+3]))*f;

			deriv[ii] -= .5*(i*((i-1.)*this->p[jj+1]+this->param[4]*pim[ii+1]*2.)-this->param[2]*this->p[ll+1])*f;
		}
		if (i == 1) {
			deriv[ii+1] -= (.125*this->param[2]*this->p[ll+4]-this->p[jj])*f;
			deriv[ii+2] += (.125*this->param[2]*this->p[ll+3])*f;

			deriv[ii] -= .5*(this->param[4]*pim[ii+1]*2.-this->param[2]*this->p[ll+1])*f;
		}
	}
	deriv[0] += .5*this->param[2]*this->p[ll+1]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin3DBeam<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p || ! pim) return;
	f *= this->R0_inv;
	for (int m, jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m >  0; m--) {
			deriv[ii+m*2-1] += ((i-m)*this->p[jj+m*2-1]+this->param[4]*pim[ii+m*2-1])*f;
			deriv[ii+m*2]   += ((i-m)*this->p[jj+m*2]+this->param[4]*pim[ii+m*2])*f;
		}
		deriv[ii] += (i*this->p[jj]+this->param[4]*pim[ii])*f;
	}
}
#endif
