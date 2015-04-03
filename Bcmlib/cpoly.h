/*========================================*/
/*                 CPOLY                  */
/*========================================*/
#ifndef ___CPOLY___
#define ___CPOLY___

#include "cshapes.h"

////////////////////////////////////////////////////////
//          													   //
//   MAPI MULTIPOLES WITH POLYNOMIAL CHARACTERISTIC   //
//          													   //
////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//...class of harmonical multipoles with polynomial behaviour;
template <typename T>
class CMapi2DPoly : public CShape<T> {
public:
		int type() {return MP2D_POLY_SHAPE;}
		int freedom(int m) { return m*2+1;}
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
//...constructor;
		CMapi2DPoly() {
			this->param = (Param *)new_struct(this->size_of_param()*sizeof(Param));
		}
};

/////////////////////////////////////////////////////////////
//...class of 3D harmonical multipoles (spherical functions);
template <typename T>
class CMapi3DPoly : public CShape<T> {
public:
		int type() {return MP3D_POLY_SHAPE;}
		int freedom(int m) { return sqr(m+1);}
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...integration;
		void primitive(T * deriv);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CMapi3DPoly() {
			this->param = (Param *)new_struct(this->size_of_param()*sizeof(Param));
		}
};

///////////////////////////////////////////////////////////
//...class of 3D harmonical multipoles for external sphere;
template <typename T>
class CMapi3DZoom : public CShape<T> {
protected:
       T * ap;
public:
		int type() {return MP3D_ZOOM_SHAPE;}
		int freedom(int m) { return sqr(m+1);}
public:
//...initialization of multipoles;
		void release() {
			delete_struct(ap);
			CShape<T>::release();
		}
public:
//...calculation of multipoles;
		void parametrization(double * P = NULL, int m = 0);
//...integration multipoles;
		void primitive  (T * deriv);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CMapi3DZoom () {
			this->param = (Param *)new_struct(this->size_of_param()*sizeof(Param));
			ap = NULL;
		}
};

////////////////////////////////////////////////////////////////////////////////////
//...class of 2D harmonical multipoles for circle exteriour with boundary condition;
template <typename T>
class CMapi2DCircle : public CShape<T> {
public:
		int type()			 { return MP2D_CIRCLE_SHAPE;}
		int size_of_param(){ return(1);}
		int freedom(int m) { return(2*m+1);}
public:
//...calculation of multipoles;
	void parametrization(double * P = NULL, int m = 0){}
//...differentiation;
		void deriv_X(T * deriv, double f = 1.){}
		void deriv_Y(T * deriv, double f = 1.){}
//...constructor;
		CMapi2DCircle() {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

////////////////////////////////////////////////////////////////////////////////////
//...class of 3D harmonical multipoles for sphere exteriour with boundary condition;
template <typename T>
class CMapi3DSphere : public CShape<T> {
public:
		int type()			 { return MP3D_SPHERE_SHAPE;}
		int size_of_param(){ return(1);}
		int freedom(int m) { return sqr(m+1);}
public:
//...calculation of multipoles;
		void parametrization		 (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CMapi3DSphere() {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_sphere_rad					1. //...радиус окружности или сферы;
#undef  SHAPE_sphere_rad						//...param(0);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CMapi2DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
	}
	if (P && this->p) {
		complex z = comp(P[0], P[1])*this->R0_inv;
		T		  hh, zm = 1., zm_i = 0.;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		int m;
		for ( this->p[0] = zm, m = 1; m <= this->N; m++) {
				this->p[m*2-1] = (zm_i = zm*imag(z)+(hh = zm_i)*real(z));
				this->p[m*2]   = (zm = zm*real(z)-hh*imag(z));
		}
	}
}

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CMapi3DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
	}
	if (P && this->p) {
		double   Z = P[2]*this->R0_inv;
		T			h1, h2, hh, f3 = Z*Z, F1, F2, F1_i, F2_i, zm, zm_i;
		complex  z = comp(P[0], P[1])*this->R0_inv;
		int      i, m; 

		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		f3   += real(z)*real(z)+imag(z)*imag(z);

///////////////////////////////////
//...calculation of the multipoles;
		for (F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N; i++) {
			  F2 = ((2.*i+1.)*Z*(this->p[i*i] = F1)-i*f3*F2)/(i+1.); 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1,	m = 1; m <= this->N; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N; i++) {
					this->p[i*i+m*2-1] = F1_i;
					this->p[i*i+m*2]   = F1;
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i;
			 this->p[i*i+m*2]   = F1;
		}
/*		int ii, i1; //...вычисление функций по слоям однородной степени;
		for (this->p[0] = 1., i = 0; i < this->N; i++) {
			T f0 = 1./(i+1.);
			this->p[i1 = (i+1)*(i+1)] = ((2.*i+1.)*Z*this->p[ii = i*i]-i*f3*this->p[(i-1)*(i-1)])*f0; 
			this->p[i1+1] = (i*this->p[ii+1]*Z+2.*this->p[ii]*imag(z))*f0;
			this->p[i1+2] = (i*this->p[ii+2]*Z+2.*this->p[ii]*real(z))*f0;
			for (m = 1; m <= i; m++) {
				this->p[i1+1+m*2] = ((i-m)*this->p[ii+1+m*2]*Z+2*(m+1)*(this->p[ii+m*2]*imag(z)+this->p[ii-1+m*2]*real(z))*(hh = 1./(i+m+1.));
				this->p[i1+2+m*2] = ((i-m)*this->p[ii+2+m*2]*Z+2*(m+1)*(this->p[ii+m*2]*real(z)-this->p[ii-1+m*2]*imag(z))*hh;
			}
		}/**/
	}
}

template <typename T>
void CMapi3DZoom<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
/////////////////////////////////////////
//...technical support (for primitive_Z);
			delete_struct(ap); ap = (T *)new_struct(this->N*sizeof(T));        
	}
	if (P && this->p) {
		double   Z = P[2]*this->R0_inv, f3 = Z*Z, rad;
		T			h1, h2, hh, F1, F2, F1_i, F2_i, zm, zm_i;
		complex  z = comp(P[0], P[1])*this->R0_inv;
		int      i, m; 

		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		this->cs[3] = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		rad = 1./sqrt(f3 += real(z)*real(z)+imag(z)*imag(z));
		Z  *=(f3 = rad*rad);
		z  *= f3; 

///////////////////////////////////
//...calculation of the multipoles;
		for (F1 = (zm = rad), F2 = (zm_i = 0.), i = 0; i < this->N+m_dop; i++) {
			  F2 = ((2.*i+1.)*Z*(this->p[i*i] = F1)-i*f3*F2)/(i+1.); 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1,	m = 1; m <= this->N+m_dop; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N+m_dop; i++) {
					this->p[i*i+m*2-1] = F1_i;
					this->p[i*i+m*2]   = F1;
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i;
			 this->p[i*i+m*2]   = F1;
		}
	}
}

///////////////////////////////////////////////////////////////
//...parametrization of the multipoles with boundary condition;
template <typename T>
void CMapi3DSphere<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = this->NN+freedom(this->N+m_dop))*sizeof(T));
	}
	if (P && this->p) {
		double   Z = P[2]*this->R0_inv, f3 = Z*Z, rad;
		T			h1, h2, hh, F1, F2, F1_i, F2_i, zm, zm_i;
		complex  z = comp(P[0], P[1])*this->R0_inv;
		int      i, m, ii; 

		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
		this->cs[3] = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		f3   += real(z)*real(z)+imag(z)*imag(z);

////////////////////////////////////////////
//...calculation of the internal multipoles;
		for (F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N; i++) {
			  F2 = ((2.*i+1.)*Z*(this->p[i*i] = F1)-i*f3*F2)/(i+1.); 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1,	m = 1; m <= this->N; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N; i++) {
					this->p[i*i+m*2-1] = F1_i;
					this->p[i*i+m*2]   = F1;
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i;
			 this->p[i*i+m*2]   = F1;
		}

////////////////////////////////////////////
//...calculation of the external multipoles;
		rad = this->param[0]*this->R0_inv/sqrt(f3); Z *= (f3 = rad*rad); z *= f3; f3 *= sqr(this->param[0]*this->R0_inv);
		for (F1 = (zm = rad), F2 = (zm_i = 0.), i = 0; i < this->N+m_dop; i++) {
			  F2 = ((2.*i+1.)*Z*(this->p[this->NN+i*i] = F1)-i*f3*F2)/(i+1.); 
			  swap(F1, F2);
		}
		for (this->p[this->NN+i*i] = F1,	m = 1; m <= this->N+m_dop; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F2 = 0., 
					F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), F2_i = 0., i = m; i < this->N+m_dop; i++) {
					this->p[this->NN+i*i+m*2-1] = F1_i;
					this->p[this->NN+i*i+m*2]   = F1;
					F2 = ((h1 = (2.*i+1.)*Z)*F1-(h2 = (i-m)*f3)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h1*F1_i-h2*F2_i)*hh; 
					swap(F1, F2);
					swap(F1_i, F2_i);
			 }
			 this->p[this->NN+i*i+m*2-1] = F1_i;
			 this->p[this->NN+i*i+m*2]   = F1;
		}

///////////////////////////////////////////
//...correction heading part of multipoles;
	   for (i = 0,   ii = 0; i < this->N;  i++)
		for (m = 2*i; m >= 0; m--) this->p[ii++] += this->p[this->NN+ii]*i/(i+1.);
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CMapi2DPoly<T>::parametrization_grad(double * P, int m_dop)
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
void CMapi3DPoly<T>::parametrization_grad(double * P, int m_dop)
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
void CMapi3DSphere<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py && ! this->pz) {
		this->px = (T *)new_struct((this->NN_grad = this->NN+freedom(this->N+m_dop))*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
		this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py && this->pz) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		this->set_inverse(ERR_STATE);
		deriv_X(this->px);
		deriv_Y(this->py);
		deriv_Z(this->pz);
		this->N += m_dop; this->set_inverse(OK_STATE);
		deriv_X(this->px);
		deriv_Y(this->py);
		deriv_Z(this->pz);
		this->N -= m_dop; this->change();
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CMapi2DPoly<T>::parametrization_hess(double * P, int m_dop)
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
void CMapi3DPoly<T>::parametrization_hess(double * P, int m_dop)
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
void CMapi3DSphere<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! this->pxx && ! this->pxy && ! this->pyy && ! this->pxz && ! this->pyz && ! this->pzz) {
		this->pxx = (T *)new_struct((this->NN_hess = this->NN+freedom(this->N+m_dop))*sizeof(T));
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
		this->set_inverse(ERR_STATE);
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
		this->N += m_dop; this->set_inverse(OK_STATE);
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
		this->N -= m_dop; this->change();
	}
}

////////////////////////////
//...integration multipoles;
template <typename T>
void CMapi3DPoly<T>::primitive(T * deriv)
{
	if (! deriv || ! this->p || ! this->N) return;
	int i, m, ii, jj; 
	T   R2 = this->p[1]*this->p[1]+this->p[2]*this->p[2]+this->p[3]*this->p[3], f0 = this->R0*this->cs[2], f;

	for (jj = 1, ii = i = 0; i < this->N; i++, jj += 2) {
		deriv[ii++] = this->p[jj++]/(i+1.)*f0;
		for (m = 1; m <= i; m++) {
			deriv[ii++] = this->p[jj++]/(i-m+1.)*f0;
			deriv[ii++] = this->p[jj++]/(i-m+1.)*f0;
		}
	}
	deriv[ii++] = ((2.*i+1.)*this->p[1]*this->p[ii]-i*R2*this->p[(jj = ii-2*i+1)++])*f0/sqr(i+1.);
	for (m = 1; m < i; m++) {
		deriv[ii++] = ((2.*i+1.)*this->p[1]*this->p[ii]-(i-m)*R2*this->p[jj++])*(f = f0/(sqr(i+1.)-sqr(m)));
		deriv[ii++] = ((2.*i+1.)*this->p[1]*this->p[ii]-(i-m)*R2*this->p[jj++])*f;
	}
	deriv[ii++] = this->p[1]*this->p[ii]*f0;
	deriv[ii++] = this->p[1]*this->p[ii]*f0;
}

////////////////////////////
//...integration multipoles;
template <typename T>
void CMapi3DZoom<T>::primitive(T * deriv)
{
	if (! deriv || ! this->p) return;
	int  jj, ii, i, m; 
	for (jj = 0, ii = i = 1; i <= this->N; i++) {
		T f0 = this->R0/(1.-2.*i);
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
	deriv[0] = this->p[0]*this->cs[3]*.5;
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CMapi2DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	int m;
	for (f *= this->R0_inv, m = this->N; m > 1; m--) {
		deriv[m*2-1] += m*this->p[m*2-3]*f;
		deriv[m*2]   += m*this->p[m*2-2]*f;
	}
	if (m > 0)
		deriv[m*2]   += m*this->p[m*2-2]*f;
}
template <typename T>
void CMapi3DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
		T ff;
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-(ff = .25*(i-m)*(i-m-1)/(m+1.))*this->p[jj+m*2+1])*f;
			deriv[ii+m*2]   += (m*this->p[jj+m*2-2]- ff*this->p[jj+m*2+2])*f;
		}
		if (m > 0) {
			deriv[ii+m*2-1] -= (ff = .25*(i-m)*(i-m-1)/(m+1.))*this->p[jj+m*2+1]*f;
			deriv[ii+m*2]   += (m*this->p[jj+m*2-2]- ff*this->p[jj+m*2+2])*f;
		}
		deriv[ii] -= .5*i*(i-1.)*this->p[jj+2]*f;
	}
}
template <typename T>
void CMapi3DZoom<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
		T ff;
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[ll+m*2-3]-(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[ll+m*2+1])*f;
			deriv[ii+m*2]   += (m*this->p[ll+m*2-2]- ff*this->p[ll+m*2+2])*f;
		}
		if (m > 0) {
			deriv[ii+m*2-1] -= (ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[ll+m*2+1]*f;
			deriv[ii+m*2]   += (m*this->p[ll+m*2-2]- ff*this->p[ll+m*2+2])*f;
		}
		deriv[ii] -= .5*(i+1.)*(i+2.)*this->p[ll+2]*f;
	}
}
template <typename T>
void CMapi3DSphere<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	if (this->id_inverse == OK_STATE) {//...инверсированные функции;
		f *= this->R0/sqr(this->param[0]);
		for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
			T ff;
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[this->NN+ii+m*2-1] += (m*this->p[this->NN+ll+m*2-3]-(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+1])*f;
				deriv[this->NN+ii+m*2]   += (m*this->p[this->NN+ll+m*2-2]- ff*this->p[this->NN+ll+m*2+2])*f;
			}
			if (m > 0) {
				deriv[this->NN+ii+m*2-1] -= (ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+1]*f;
				deriv[this->NN+ii+m*2]   += (m*this->p[this->NN+ll+m*2-2]- ff*this->p[this->NN+ll+m*2+2])*f;
			}
			deriv[this->NN+ii] -= .5*(i+1.)*(i+2.)*this->p[this->NN+ll+2]*f;
		}
	   for (int i = 0,   ii = 0; i < this->N; i++)
		for (int m = 2*i; m >= 0; m--) deriv[ii++] += deriv[this->NN+ii]*i/(i+1.);
	}
	else {//...прямые функции;
		f *= this->R0_inv;
		for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
			T ff;
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*(this->p[jj+m*2-3]-this->p[this->NN+jj+m*2-3]*(i-1.)/i)-(ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+1]-this->p[this->NN+jj+m*2+1]*(i-1.)/i))*f;
				deriv[ii+m*2]   += (m*(this->p[jj+m*2-2]-this->p[this->NN+jj+m*2-2]*(i-1.)/i)- ff*(this->p[jj+m*2+2]-this->p[this->NN+jj+m*2+2]*(i-1.)/i))*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] -= (ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+1]-this->p[this->NN+jj+m*2+1]*(i-1.)/i)*f;
				deriv[ii+m*2]   += (m*(this->p[jj+m*2-2]-this->p[this->NN+jj+m*2-2]*(i-1.)/i)- ff*(this->p[jj+m*2+2]-this->p[this->NN+jj+m*2+2]*(i-1.)/i))*f;
			}
			deriv[ii] -= .5*i*(i-1.)*(this->p[jj+2]-this->p[this->NN+jj+2]*(i-1.)/i)*f;
		}
		if (! this->id_inverse) { //...корректируем головную часть инверсированными функциями;
			f /= sqr(this->param[0]*this->R0_inv);
			for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
				T ff;
				for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
					deriv[ii+m*2-1] += (m*this->p[this->NN+ll+m*2-3]-(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+1])*f*i/(i+1.);
					deriv[ii+m*2]   += (m*this->p[this->NN+ll+m*2-2]- ff*this->p[this->NN+ll+m*2+2])*f*i/(i+1.);
				}
				if (m > 0) {
					deriv[ii+m*2-1] -= (ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+1]*f*i/(i+1.);
					deriv[ii+m*2]   += (m*this->p[this->NN+ll+m*2-2]- ff*this->p[this->NN+ll+m*2+2])*f*i/(i+1.);
				}
				deriv[ii] -= .5*(i+1.)*(i+2.)*this->p[this->NN+ll+2]*f*i/(i+1.);
			}
		}
	}
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CMapi2DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	int m;
	for (f *= this->R0_inv, m = this->N; m > 1; m--) {
		deriv[m*2-1] += m*this->p[m*2-2]*f;
		deriv[m*2]   -= m*this->p[m*2-3]*f;
	}
	if (m > 0)
		deriv[m*2-1] += m*this->p[m*2-2]*f;
}
template <typename T>
void CMapi3DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
		T ff;
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff = .25*(i-m)*(i-m-1)/(m+1.))*this->p[jj+m*2+2])*f;
			deriv[ii+m*2]   -= (m*this->p[jj+m*2-3]+ ff*this->p[jj+m*2+1])*f;
		}
		if (m > 0) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+(ff = .25*(i-m)*(i-m-1)/(m+1.))*this->p[jj+m*2+2])*f;
			deriv[ii+m*2]   -= ff*this->p[jj+m*2+1]*f;
		}
		deriv[ii] -= .5*i*(i-1.)*this->p[jj+1]*f;
	}
}
template <typename T>
void CMapi3DZoom<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, ll = sqr((this->N+1)), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
		T ff;
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[ll+m*2+2])*f;
			deriv[ii+m*2]   -= (m*this->p[ll+m*2-3]+ ff*this->p[ll+m*2+1])*f;
		}
		if (m > 0) {
			deriv[ii+m*2-1] += (m*this->p[ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[ll+m*2+2])*f;
			deriv[ii+m*2]   -= ff*this->p[ll+m*2+1]*f;
		}
		deriv[ii] -= .5*(i+1.)*(i+2.)*this->p[ll+1]*f;
	}
}
template <typename T>
void CMapi3DSphere<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	if (this->id_inverse == OK_STATE) {//...инверсированные функции;
		f *= this->R0/sqr(this->param[0]);
		for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
			T ff;
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[this->NN+ii+m*2-1] += (m*this->p[this->NN+ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+2])*f;
				deriv[this->NN+ii+m*2]   -= (m*this->p[this->NN+ll+m*2-3]+ ff*this->p[this->NN+ll+m*2+1])*f;
			}
			if (m > 0) {
				deriv[this->NN+ii+m*2-1] += (m*this->p[this->NN+ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+2])*f;
				deriv[this->NN+ii+m*2]   -= ff*this->p[this->NN+ll+m*2+1]*f;
			}
			deriv[this->NN+ii] -= .5*(i+1.)*(i+2.)*this->p[this->NN+ll+1]*f;
		}
	   for (int i = 0,   ii = 0; i < this->N; i++)
		for (int m = 2*i; m >= 0; m--) deriv[ii++] += deriv[this->NN+ii]*i/(i+1.);
	}
	else {//...прямые функции;
		f *= this->R0_inv;
		for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
			T ff;
			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				deriv[ii+m*2-1] += (m*(this->p[jj+m*2-2]-this->p[this->NN+jj+m*2-2]*(i-1.)/i)+(ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+2]-this->p[this->NN+jj+m*2+2]*(i-1.)/i))*f;
				deriv[ii+m*2]   -= (m*(this->p[jj+m*2-3]-this->p[this->NN+jj+m*2-3]*(i-1.)/i)+ ff*(this->p[jj+m*2+1]-this->p[this->NN+jj+m*2+1]*(i-1.)/i))*f;
			}
			if (m > 0) {
				deriv[ii+m*2-1] += (m*(this->p[jj+m*2-2]-this->p[this->NN+jj+m*2-2]*(i-1.)/i)+(ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+2]-this->p[this->NN+jj+m*2+2]*(i-1.)/i))*f;
				deriv[ii+m*2]   -= ff*(this->p[jj+m*2+1]-this->p[this->NN+jj+m*2+1]*(i-1.)/i)*f;
			}
			deriv[ii] -= .5*i*(i-1.)*(this->p[jj+1]-this->p[this->NN+jj+1]*(i-1.)/i)*f;
		}
		if (! this->id_inverse) { //...корректируем головную часть инверсированными функциями;
			f /= sqr(this->param[0]*this->R0_inv);
			for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
				T ff;
				for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
					deriv[ii+m*2-1] += (m*this->p[this->NN+ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+2])*f*i/(i+1.);
					deriv[ii+m*2]   -= (m*this->p[this->NN+ll+m*2-3]+ ff*this->p[this->NN+ll+m*2+1])*f*i/(i+1.);
				}
				if (m > 0) {
					deriv[ii+m*2-1] += (m*this->p[this->NN+ll+m*2-2]+(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->p[this->NN+ll+m*2+2])*f*i/(i+1.);
					deriv[ii+m*2]   -= ff*this->p[this->NN+ll+m*2+1]*f*i/(i+1.);
				}
				deriv[ii] -= .5*(i+1.)*(i+2.)*this->p[this->NN+ll+1]*f*i/(i+1.);
			}
		}
	}
}

////////////////////////////////
//...differentiation multipoles;
template <typename T>
void CMapi3DPoly<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
			deriv[ii+m*2-1] += (i-m)*this->p[jj+m*2-1]*f;
			deriv[ii+m*2]   += (i-m)*this->p[jj+m*2]*f;
		}
		deriv[ii] += i*this->p[jj]*f;
	}
}
template <typename T>
void CMapi3DZoom<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 0; m--) {
			deriv[ii+m*2-1] -= (i+m+1)*this->p[ll+m*2-1]*f;
			deriv[ii+m*2]   -= (i+m+1)*this->p[ll+m*2]*f;
		}
		deriv[ii] -= (i+1)*this->p[ll]*f;
	}
}
template <typename T>
void CMapi3DSphere<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	if (this->id_inverse == OK_STATE) {//...инверсированные функции;
		f *= this->R0/sqr(this->param[0]);
		for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
			for (jj = (i-1)*(i-1), m = i; m > 0; m--) {
				deriv[this->NN+ii+m*2-1] -= (i+m+1)*this->p[this->NN+ll+m*2-1]*f;
				deriv[this->NN+ii+m*2]   -= (i+m+1)*this->p[this->NN+ll+m*2]*f;
			}
			deriv[this->NN+ii] -= (i+1)*this->p[this->NN+ll]*f;
		}
	   for (int i = 0,   ii = 0; i < this->N; i++)
		for (int m = 2*i; m >= 0; m--) deriv[ii++] += deriv[this->NN+ii]*i/(i+1.);
	}
	else {//...прямые функции;
		f *= this->R0_inv;
		for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ii = jj) {
			for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
				deriv[ii+m*2-1] += (i-m)*(this->p[jj+m*2-1]-this->p[this->NN+jj+m*2-1]*(i-1.)/i)*f;
				deriv[ii+m*2]   += (i-m)*(this->p[jj+m*2]-this->p[this->NN+jj+m*2]*(i-1.)/i)*f;
			}
			deriv[ii] += i*(this->p[jj]-this->p[this->NN+jj]*(i-1.)/i)*f;
		}
		if (! this->id_inverse) { //...корректируем головную часть инверсированными функциями;
			f /= sqr(this->param[0]*this->R0_inv);
			for (int m, ll = sqr(this->N+1), jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ll = ii, ii = jj) {
				for (jj = (i-1)*(i-1), m = i; m > 0; m--) {
					deriv[ii+m*2-1] -= (i+m+1)*this->p[this->NN+ll+m*2-1]*f*i/(i+1.);
					deriv[ii+m*2]   -= (i+m+1)*this->p[this->NN+ll+m*2]*f*i/(i+1.);
				}
				deriv[ii] -= (i+1)*this->p[this->NN+ll]*f*i/(i+1.);
			}
		}
	}
}

///////////////////////////////////////////////////////////////
//          																 //
//   EXPONENTIAL MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
//          																 //
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//...class of skin multipoles with polynomial behaviour;
template <typename T>
class CSkin2DPoly : public CShape<T> {
protected:
      T * sk;
		static double NORMALIZATION_EXPO_LIMIT;
public:
		int type()          { return SK2D_POLY_SHAPE;}
		int freedom (int m) { return m*2+1;}
		int size_of_param() { return(4);}
		void release() {
			 delete_struct(sk);
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
//...constructor;
		CSkin2DPoly () {
			delete_struct(this->param);
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			sk    = NULL;
		}
};
template <typename T> double CSkin2DPoly<T>::NORMALIZATION_EXPO_LIMIT = 5.;

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);

//////////////////////////////////////////////////////////////////////
//...class of skin multipoles with polynomial behaviour along the line;
template <typename T>
class CSkin3DPoly : public CShape<T> {
protected:
		T * sk, * E1;
public:
		int type()          { return SK3D_POLY_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);};
		void release() {
			 delete_struct(sk);
			 delete_struct(E1);
			 CShape<T>::release();
		}
public:
//...initialization and calculation of multipoles;
		void parametrization(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CSkin3DPoly () {
			delete_struct(this->param);
			this->param   = (Param *)new_struct(size_of_param()*sizeof(Param));
			sk = E1 = NULL;
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CSkin2DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
///////////////////////
//...technical support;
			delete_struct(sk); sk = (T *)new_struct((this->N+1+m_dop)*sizeof(T));
	}
	if (P && this->p) {
		complex w = comp(P[0], P[1])*this->R0_inv;
		T		 zm = 1., zm_i = 0., hh;
		double rr = abs(w);	int m;
		if (rr > EE)
			w *= (1./rr); else w = comp(1.);
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		if (sk) {
			T z = this->param[0]*rr, d = z*.5,  rr0 = d*d, z_inv, h0, h1, h2;

			if (this->param[0] <= NORMALIZATION_EXPO_LIMIT)
			for ( m = 0, h0 = 1.; m <= this->N+m_dop; m++, h0 *= rr) {
				T   f = 1., P = f; 
				int k = 1;
				do {
					 P += (f *= rr0/((T)k*(m+k))); ++k;
				}
				while (fabs(f) > EE);
				sk[m] = h0*P;
			}
			else {
				if ( to_double(rr0) >= 100. && to_double(z) >= this->N+m_dop+5.0) {
					z_inv = .5/z;
					h1 = sqrt(to_double(2.*z_inv*M_2_PI))*exp(to_double( z-this->param[0]));
					h2 = sqrt(to_double(2.*z_inv*M_2_PI))*exp(to_double(-z-this->param[0]));
					for ( m = 0; m <= this->N+m_dop; m++) { 
						T   f = 1., P = f, Q = 0., h;
						int k = 1;
						do {
							  h  =  f;
							  Q += (f *= (m-k+.5)*(m+k-.5)*z_inv/(T)k); ++k;
							  P += (f *= (m-k+.5)*(m+k-.5)*z_inv/(T)k); ++k;
						}
						while (k <= m-.5 || fabs(f) > EE && fabs(f) < fabs(h));
						sk[m] = .5*(h1+h2)*P-.5*(h1-h2)*Q;
					}
				}
				else
				for ( m = 0, h0 = exp(to_double(-this->param[0])); m <= this->N+m_dop; h0 *= d/(m+1.), m++) {
					T   f = 1., P = f; 
					int k = 1;
					do {
						 P += (f *= rr0/((T)k*(m+k))); ++k;
					}
					while (fabs(f) > EE);
					sk[m] = h0*P;
				}
			}
		}
		for ( this->p[0] = sk[0], m = 1; m <= this->N+m_dop; m++) {
				this->p[m*2-1] = sk[m]*(zm_i = zm*imag(w)+(hh = zm_i)*real(w));
				this->p[m*2]   = sk[m]*(zm = zm*real(w)-hh*imag(w));
		}
	}
}

template <typename T>
void CSkin3DPoly<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
///////////////////////
//...technical support;
			delete_struct(sk); sk = (T *)new_struct(( this->N+1+m_dop)*sizeof(T));        
			delete_struct(E1); E1 = (T *)new_struct(((this->N+m_dop)/2+1)*sizeof(T));
	}
	if (P && this->p) {
		int  i, l, nm, m;
		complex w = comp(P[0], P[1])*this->R0_inv;
		T		 zm = 1., zm_i = 0., d, F, hh;
		double rr = abs(w),
				r0  = rr*.25,
				Z   = P[2]*this->R0_inv,
				Z0  = Z*Z, f1, f2;
		if (rr > EE)
			w *= (1./rr); else w = comp(1.);
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];

///////////////////////////////////
//...calculation of the multipoles;
		if (sk) {
			T  z = this->param[0]*rr, rr0 = z*z*.25, h0;
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
		for (f1 = 1., f2 = Z, i = 0; i <= this->N+m_dop; i++) {
			for ( E1[l = 0] = 1., nm = i/2; l < nm; l++)
					E1[l+1] = -r0*(i-2*l-1.)*(i-2*l)/sqr(l+1.)*E1[l];
			for ( F  = E1[nm]*sk[nm]*(d  = f1), l = nm-1; l >= 0; l--)
					F += E1[l ]*sk[l ]*(d *= Z0);
			this->p[i*i] = F; swap(f1, f2);
		}
		for (m = 1; m <= this->N+m_dop; m++)
		for (zm_i = zm*imag(w)+(hh = zm_i)*real(w), zm = zm*real(w)-hh*imag(w), 
							 f1 = 1., f2 = Z, i = m; i <= this->N+m_dop; i++) {
			for ( E1[l = 0] = 1., nm = (i-m)/2;  l <  nm;      l++)
					E1[l+1] = -r0*(i-m-2*l-1.)*(i-m-2*l)/((l+1.)*(m+l+1.))*E1[l];
			for ( F  = E1[nm]*sk[m+nm]*(d  = f1), l = nm-1; l >= 0; l--)
					F += E1[l ]*sk[m+l ]*(d *= Z0);
			this->p[i*i+m*2-1] = F*zm_i;
			this->p[i*i+m*2]   = F*zm; swap(f1, f2);
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CSkin2DPoly<T>::parametrization_grad(double * P, int m_dop)
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

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CSkin2DPoly<T>::parametrization_hess(double * P, int m_dop)
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

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	if (this->param[0] > NORMALIZATION_EXPO_LIMIT) {
		f *= .5*this->param[0];
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (this->p[m*2-3]+this->p[m*2+1])*f;
			 deriv[m*2]   += (this->p[m*2-2]+this->p[m*2+2])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] +=  this->p[m*2+1]*f;
			 deriv[m*2]   += (this->p[m*2-2]+this->p[m*2+2])*f;
		}
		deriv[0] += this->p[2]*2.*f;
	}
	else {
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-3]+.25*this->param[2]/(m+1.)*this->p[m*2+1])*f;
			 deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]/(m+1.)*this->p[m*2+2])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += .25*this->param[2]/(m+1.)*this->p[m*2+1]*f;
			 deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]/(m+1.)*this->p[m*2+2])*f;
		}
		deriv[0] += .5*this->param[2]*this->p[2]*f;
	}
}

template <typename T>
void CSkin3DPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int ll = sqr(this->N+1);
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-3]-.25/(m+1.)*((i-m)*(i-m-1.)*this->p[jj+m*2+1]-this->param[2]*this->p[ll+m*2+1]))*f;
			deriv[ii+m*2]   += (m*this->p[jj+m*2-2]-.25/(m+1.)*((i-m)*(i-m-1.)*this->p[jj+m*2+2]-this->param[2]*this->p[ll+m*2+2]))*f;
		}
		if (m > 0) {
			deriv[ii+1] -= (.125*((i-1.)*(i-2.)*this->p[jj+3]-this->param[2]*this->p[ll+3]))*f;
			deriv[ii+2] -= (.125*((i-1.)*(i-2.)*this->p[jj+4]-this->param[2]*this->p[ll+4])-this->p[jj])*f;
		}
		deriv[ii] -= .5*(i*(i-1.)*this->p[jj+2]-this->param[2]*this->p[ll+2])*f;
	}
	deriv[0] += .5*this->param[2]*this->p[ll+2]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	if (this->param[0] > NORMALIZATION_EXPO_LIMIT) {
		f *= .5*this->param[0];
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (this->p[m*2-2]-this->p[m*2+2])*f;
			 deriv[m*2]   -= (this->p[m*2-3]-this->p[m*2+1])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += (this->p[m*2-2]-this->p[m*2+2])*f;
			 deriv[m*2]   +=  this->p[m*2+1]*f;
		}
		deriv[0] += this->p[1]*2.*f;
	}
	else {
		for (m = this->N; m > 1; m--) {
			 deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]/(m+1.)*this->p[m*2+2])*f;
			 deriv[m*2]   -= (m*this->p[m*2-3]-.25*this->param[2]/(m+1.)*this->p[m*2+1])*f;
		}
		if (m > 0) {
			 deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]/(m+1.)*this->p[m*2+2])*f;
			 deriv[m*2]   += .25*this->param[2]/(m+1.)*this->p[m*2+1]*f;
		}
		deriv[0] += .5*this->param[2]*this->p[1]*f;
	}
}

template <typename T>
void CSkin3DPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int ll = sqr(this->N+1);
	for (int m, jj, ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
		for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
			deriv[ii+m*2-1] += (m*this->p[jj+m*2-2]+.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+2]-this->param[2]*this->p[ll+m*2+2]))*f;
			deriv[ii+m*2]   -= (m*this->p[jj+m*2-3]+.25/(m+1.)*((i-m)*(i-m-1)*this->p[jj+m*2+1]-this->param[2]*this->p[ll+m*2+1]))*f;
		}
		if (m > 0) {
			deriv[ii+1] += (.125*((i-1)*(i-2)*this->p[jj+4]-this->param[2]*this->p[ll+4])+this->p[jj])*f;
			deriv[ii+2] -= (.125*((i-1)*(i-2)*this->p[jj+3]-this->param[2]*this->p[ll+3]))*f;
		}
		deriv[ii] -= .5*(i*(i-1)*this->p[jj+1]-this->param[2]*this->p[ll+1])*f;
	}
	deriv[0] += .5*this->param[2]*this->p[ll+1]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin3DPoly<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int m, jj, ii = sqr(this->N), i = this->N; i >= 0; i--, ii = jj) {
		for (jj = (i-1)*(i-1), m = i-1; m > 0; m--) {
			deriv[ii+m*2-1] += (i-m)*this->p[jj+m*2-1]*f;
			deriv[ii+m*2]   += (i-m)*this->p[jj+m*2]*f;
		}
		deriv[ii] += i*this->p[jj]*f;
	}
}

//////////////////////////////////////////////////////////////
//																				//
//        MAPI MULTIPOLES FOR CELL WITH INCLUSIONS          //
//																				//
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//...class of mapi multipoles for circle layer inclusion;
template <typename T>
class CMapi2DClayer : public CShape<T> {
protected:
		T * pim;
		complex * E1;
public:
		int type() {return MP2D_CLAYER_SHAPE;}
		int freedom(int m) { return m*2+1;}
		int size_of_param() { return(5);}
public:
//...initialization of multipoles;
		void set_shape(double R0 = 1., double K1 = 1., double K2 = 1., double K3 = 1., double R1 = 1., double R2 = 1) {
			CMapi2DClayer<T>:: R0 = R0;
			this->R0_inv     = R0 > EE ? 1./R0 : 1.;
			int   k    = -1;
			if (++k < size_of_param()) this->param[k] = (Param)(K1);
			if (++k < size_of_param()) this->param[k] = (Param)(K2);
			if (++k < size_of_param()) this->param[k] = (Param)(K3);
			if (++k < size_of_param()) this->param[k] = (Param)(R1);
			if (++k < size_of_param()) this->param[k] = (Param)(R2);
		}		
		void release() {
			delete_struct(pim);
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
//...constructor;
		CMapi2DClayer() {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			pim = NULL;
			E1  = NULL;
		};
};

////////////////////////////////////////////////////////////
//...class of mapi multipoles for spherical layer inclusion;
template <typename T>
class CMapi3DClayer : public CMapi2DClayer<T> {
public:
		int type() {return MP3D_CLAYER_SHAPE;}
		int freedom(int m) { return sqr(m+1);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
		void deriv_Z(T * deriv, double f = 1.);
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_mapi_K1					0. //...heat transfer coefficient in inclusion;
#undef  SHAPE_mapi_K1						//...param(0);
#define SHAPE_mapi_K2               0. //...heat transfer coefficient in layer;
#undef  SHAPE_mapi_K2                  //...param(1);
#define SHAPE_mapi_K3               0. //...heat transfer coefficient in matrix;
#undef  SHAPE_mapi_K3                  //...param(2);
#define SHAPE_mapi_R1               0. //...radius of inclusion;
#undef  SHAPE_mapi_R1                  //...param(3);
#define SHAPE_mapi_R2               0. //...radius of layer zone;
#undef  SHAPE_mapi_R2                  //...param(4);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CMapi2DClayer<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
			pim 	= (T *)new_struct( this->NN_dop*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(E1); E1 = (complex *)new_struct((this->N+m_dop+1)*sizeof(complex));
	}
	if (P && this->p && pim) {
		complex z = this->R0_inv*comp(P[0], P[1])/this->param[4], w;
		T		zm = 1., zm_i = 0., wm = 1., wm_i = 0., h1, h2;
		double  q = this->param[3]/this->param[4], p1 = this->param[1]/this->param[0], p2 = this->param[1]/this->param[2], p3 = 1.-p1, p4 = 1.-p2, f = 1., Q, A, B;
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		int m;
		if ((this->id_inverse = ((A = abs(z)*this->R0) > 1.))) { //...external system;
			w = this->R0_inv*this->param[4]/comp(P[0], -P[1]);
			Q = (1.+p1)*(1.+p2)-p3*p4;
			B = (p4*(1.+p1)-p3*(1.+p2))/Q; E1[0] = comp(1., B);

			for ( this->p[0] = 1.+(pim[0] = B), m = 1; m <= this->N+m_dop; m++) {
					f *= q*q;
					Q = (1.+p1)*(1.+p2)- p3*p4*f;
					B = (p4*(1.+p1)-p3*(1.+p2)*f)/Q; E1[m] = comp(1., B);

					this->p[m*2-1] = (zm_i = zm*imag(z)+(h1 = zm_i)*real(z))+(pim[m*2-1] = (wm_i = wm*imag(w)+(h2 = wm_i)*real(w))*B);
					this->p[m*2]   = (zm = zm*real(z)-h1*imag(z))+(pim[m*2] = (wm = wm*real(w)-h2*imag(w))*B);
			}
		}
		else
		if (! (this->id_inverse = (A > q))) { //...internal system;
			Q = (1.+p1)*(1.+p2)-p3*p4;
			A = 4.*p1/Q; E1[0] = comp(A);

			for ( this->p[0] = A+(pim[0] = 0.), m = 1; m <= this->N+m_dop; m++) {
					f *= q*q;
					Q = (1.+p1)*(1.+p2)-p3*p4*f;
					A =  4.*p1/Q; E1[m] = comp(A);

					this->p[m*2-1] = (zm_i = zm*imag(z)+(h1 = zm_i)*real(z))*A+(pim[m*2-1] = 0.);
					this->p[m*2]   = (zm = zm*real(z)-h1*imag(z))*A+(pim[m*2] = 0.);
			}
		}
		else { //...intermediate system;
			this->id_inverse = -1;
			w = this->R0_inv*this->param[3]/comp(P[0], -P[1]);
			Q = (1.+p1)*(1.+p2)-p3*p4;
			A = (1.+p1)*2./Q;
			B = -p3*2./Q; E1[0] = comp(A, B);

			for ( this->p[0] = A+(pim[0] = B), m = 1; m <= this->N+m_dop; m++) {
					f *= q;
					Q = (1.+p1)*(1.+p2)-p3*p4*f*f;
					A = (1.+p1)*2./Q;
					B = -f*p3*2./Q; E1[m] = comp(A, B);

					this->p[m*2-1] = (zm_i = zm*imag(z)+(h1 = zm_i)*real(z))*A+(pim[m*2-1] = (wm_i = wm*imag(w)+(h2 = wm_i)*real(w))*B);
					this->p[m*2]   = (zm = zm*real(z)-h1*imag(z))*A+(pim[m*2] = (wm = wm*real(w)-h2*imag(w))*B);
			}
		}
	}
}

template <typename T>
void CMapi3DClayer<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p   = (T *)new_struct((this->NN_dop = freedom(this->N+m_dop))*sizeof(T));
			this->pim = (T *)new_struct( this->NN_dop*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(this->E1); this->E1 = (complex *)new_struct((this->N+m_dop+1)*sizeof(complex));
	}
	if (P && this->p && this->pim) {
		int i, m, k, m2;
		double  R = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])), RR = 1., Z = R < EE ? 1.: P[2]*(RR = 1./R), f1, f2, 
				  q = this->param[3]/this->param[4], p1 = this->param[1]/this->param[0], p2 = this->param[1]/this->param[2], p3 = 1.-p1, p4 = 1.-p2, f, Q, A, B, pp, d;
		T		  zm, zm_i, F1, F2, F1_i, F2_i, hh, h0;
		complex z = R < EE ? comp(1.) : comp(P[0], P[1])*RR;
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = P[5];
			
/////////////////////////////////////////////////
//...calculation of the angle part of multipoles;
		for (F1 = (zm = 1.),  F2 = (zm_i = 0.), i = 0; i < this->N+m_dop; i++) {
			  F2 = ((2.*i+1.)*Z*(this->p[i*i] = F1)-i*F2)/(i+1.); 
			  swap(F1, F2);
		}
		for (this->p[i*i] = F1, m = 1; m <= this->N+m_dop; m++) {
			 for (F1 = (zm = (hh = zm)*real(z)-zm_i*imag(z)), F1_i = (zm_i = hh*imag(z)+zm_i*real(z)), 
					F2 = 0., F2_i = 0., i = m; i < this->N+m_dop; i++) {
					this->p[i*i+m*2-1] = F1_i;
					this->p[i*i+m*2]   = F1;
					F2 = ((h0 = (2.*i+1.)*Z)*F1-(i-m)*F2)*(hh = 1./(i+m+1.)); 
					F2_i = (h0*F1-(i-m)*F2)*hh; 
					swap(F1, F2); swap(F1_i, F2_i);
			 }
			 this->p[i*i+m*2-1] = F1_i;
			 this->p[i*i+m*2]   = F1;
		}
		
///////////////////////////////////
//...calculation of the multipoles;
		this->E1[0] = comp(1.);	 R *= this->R0_inv/this->param[4];
		if ((this->id_inverse = ((A = R*this->R0) > 1.))) { //...external system;
			f = q; f1 = 1.; f2 = (RR *= this->R0_inv*this->param[4])*this->R0;
			for ( this->p[0] = 1.+(this->pim[0] = 0.), k = m = 1; m <= this->N+m_dop; m++) {
					f *= q*q;
					Q = (pp = (1.+p1*(m+1.)/m))*(1.+p2*m/(m+1.))-p3*p4*f;
					B = (pp-p3*f)/Q*(2.*m+1.)/(m+1.)-1.; this->E1[m] = comp(1., B);

					for (f1 *= R, f2 *= RR, m2 = m*2, i = 0; i <= m2; i++, k++)
						this->p[k] = f1*this->p[k]+( this->pim[k] = f2*this->p[k]*B);
			}
		}
		else
		if (! (this->id_inverse = (A > q))) { //...internal system;
			f = q; f1 = 1.;
			for ( this->p[0] = 1.+(this->pim[0] = 0.), k = m = 1; m <= this->N+m_dop; m++) {
					f *= q*q;
					Q = (pp = (1.+p1*(m+1.)/m))*(1.+p2*m/(m+1.))-p3*p4*f;
					A = (pp-p3)/Q*(2.*m+1.)/(m+1.); this->E1[m] = comp(A);

					for (f1 *= R, m2 = m*2,   i = 0; i <= m2; i++, k++)
						this->p[k] = f1*this->p[k]*A+(this->pim[k] = 0.);
			}
		}
		else { //...intermediate system;
			this->id_inverse = -1;
			f = 1.; f1 = 1.; f2 = (RR *= this->R0_inv*this->param[3])*this->R0;
			for (this->p[0] = 1.+(this->pim[0] = 0.), k = m = 1; m <= this->N+m_dop; m++) {
					f *= q;
					Q = (pp = (1.+p1*(m+1.)/m))*(1.+p2*m/(m+1.))-p3*p4*f*f*q;
					A =  pp/Q*( d = (2.*m+1.)/(m+1.));
					B = -f*p3/Q*d; this->E1[m] = comp(A, B);

					for (f1 *= R, f2 *= RR,  m2 = m*2, i = 0; i <= m2; i++, k++)
						this->p[k] = f1*this->p[k]*A+(this->pim[k] = f2*this->p[k]*B);
			}
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CMapi2DClayer<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py) {
			this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m_dop))*sizeof(T));
			this->py = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));

////////////////////////////////////////////
//...calculation gradient of the multipoles;
		this->N += m_dop;
		if (this->id_inverse > 0) { //...external system;
			for ( int m = 1; m <= this->N; m++) {
				double A = real(this->E1[m])/(real(this->E1[m-1])*this->param[4])*m*this->R0_inv, 
					B = fabs(imag(this->E1[m+1])) > EE ? imag(this->E1[m])/(imag(this->E1[m+1])*this->param[4])*m*this->R0 : 0.;
				if (m != 1)
				this->px[m*2-1]  = A*(this->p[m*2-3]-this->pim[m*2-3]);
				this->px[m*2-1] -= B*this->pim[m*2+1];
				this->py[m*2] = -this->px[m*2-1]-B*2.*this->pim[m*2+1];
				this->px[m*2] =  this->py[m*2-1] = A*(this->p[m*2-2]-this->pim[m*2-2])-B*this->pim[m*2+2];
				this->py[m*2-1] += B*2.*this->pim[m*2+2];
			}
		}
		else
		if (! this->id_inverse) { //...internal system;
			for ( int m = 1; m <= this->N; m++) {
				double A = real(this->E1[m])/(real(this->E1[m-1])*this->param[4])*m*this->R0_inv;
				if (m != 1)
				this->px[m*2-1] = A*(this->p[m*2-3]-this->pim[m*2-3]);
				this->py[m*2] = -this->px[m*2-1];
				this->px[m*2] =  this->py[m*2-1] = A*(this->p[m*2-2]-this->pim[m*2-2]);
			}
		}
		else { //...intermediate system;
			for ( int m = 1; m <= this->N; m++) {
				double A = real(this->E1[m])/(real(this->E1[m-1])*this->param[4])*m*this->R0_inv, 
					B = fabs(imag(this->E1[m+1])) > EE ? imag(this->E1[m])/(imag(this->E1[m+1])*this->param[3])*m*this->R0 : 0.;
				if (m != 1)
				this->px[m*2-1]  = A*(this->p[m*2-3]-this->pim[m*2-3]);
				this->px[m*2-1] -= B*this->pim[m*2+1];
				this->py[m*2] = -this->px[m*2-1]-B*2.*this->pim[m*2+1];
				this->px[m*2] =  this->py[m*2-1] = A*(this->p[m*2-2]-this->pim[m*2-2])-B*this->pim[m*2+2];
				this->py[m*2-1] += B*2.*pim[m*2+2];
			}
		}
		this->N -= m_dop;
	}
}

template <typename T>
void CMapi3DClayer<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! this->px && ! this->py && ! this->pz) {
			this->px = (T *)new_struct((this->NN_grad = freedom(this->N+m_dop))*sizeof(T));
			this->py = (T *)new_struct( this->NN_grad*sizeof(T));
			this->pz = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		memset (this->pz, 0, this->NN_grad*sizeof(T));

//////////////////////////////////////////////////
//...calculation gradient of the space multipoles;
		double B; int ii, jj, ll, i, m; this->N += m_dop;	
		for (ll = sqr(this->N+1), ii = sqr(this->N), i = this->N; i > 0; i--, ll = ii, ii = jj) {
			double A = real(this->E1[i])/(real(this->E1[i-1])*this->param[4])*this->R0_inv, ff;
			if (this->id_inverse > 0 && fabs(imag(this->E1[i+1])) > EE) //...external system;
				B = imag(this->E1[i])/(imag(this->E1[i+1])*this->param[4])*this->R0;	else
			if (this->id_inverse < 0 && fabs(imag(this->E1[i+1])) > EE) //...intermediate system;
				B = imag(this->E1[i])/(imag(this->E1[i+1])*this->param[3])*this->R0; else B = 0.;

			for (jj = (i-1)*(i-1), m = i; m > 1; m--) {
				this->px[ii+m*2-1] =  A*(m*(this->p[jj+m*2-3]-this->pim[jj+m*2-3])-(ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+1]-this->pim[jj+m*2+1]));
				this->px[ii+m*2]   =  A*(m*(this->p[jj+m*2-2]-this->pim[jj+m*2-2])- ff*(this->p[jj+m*2+2]-this->pim[jj+m*2+2]));
				this->py[ii+m*2-1] =  A*(m*(this->p[jj+m*2-2]-this->pim[jj+m*2-2])+ ff*(this->p[jj+m*2+2]-this->pim[jj+m*2+2]));
				this->py[ii+m*2]   = -A*(m*(this->p[jj+m*2-3]-this->pim[jj+m*2-3])+ ff*(this->p[jj+m*2+1]-this->pim[jj+m*2+1]));
				this->pz[ii+m*2-1] =  A*(i-m)*(this->p[jj+m*2-1]-this->pim[jj+m*2-1]);
				this->pz[ii+m*2]   =  A*(i-m)*(this->p[jj+m*2]-this->pim[jj+m*2]);

				this->px[ii+m*2-1] += B*(m*this->pim[ll+m*2-3]-(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->pim[ll+m*2+1]);
				this->px[ii+m*2]   += B*(m*this->pim[ll+m*2-2]- ff*this->pim[ll+m*2+2]);
				this->py[ii+m*2-1] += B*(m*this->pim[ll+m*2-2]+ ff*this->pim[ll+m*2+2]);
				this->py[ii+m*2]   -= B*(m*this->pim[ll+m*2-3]+ ff*this->pim[ll+m*2+1]);
				this->pz[ii+m*2-1] -= B*(i+m+1)*this->pim[ll+m*2-1];
				this->pz[ii+m*2]   -= B*(i+m+1)*this->pim[ll+m*2];

			}
			if (m > 0) {
				this->px[ii+m*2-1] =										  -A*(ff = .25*(i-m)*(i-m-1)/(m+1.))*(this->p[jj+m*2+1]-this->pim[jj+m*2+1]);
				this->px[ii+m*2]   = A*(m*(this->p[jj+m*2-2]-this->pim[jj+m*2-2])-ff*(this->p[jj+m*2+2]-this->pim[jj+m*2+2]));
				this->py[ii+m*2-1] = A*(m*(this->p[jj+m*2-2]-this->pim[jj+m*2-2])+ff*(this->p[jj+m*2+2]-this->pim[jj+m*2+2]));
				this->py[ii+m*2]   =											-A*ff*(this->p[jj+m*2+1]-this->pim[jj+m*2+1]);
				this->pz[ii+m*2-1] =  A*(i-m)*(this->p[jj+m*2-1]-this->pim[jj+m*2-1]);
				this->pz[ii+m*2]   =  A*(i-m)*(this->p[jj+m*2]-this->pim[jj+m*2]);

				this->px[ii+m*2-1] -=						B*(ff = .25*(i+m+1)*(i+m+2)/(m+1.))*this->pim[ll+m*2+1];
				this->px[ii+m*2]   += B*(m*this->pim[ll+m*2-2]- ff*this->pim[ll+m*2+2]);
				this->py[ii+m*2-1] += B*(m*this->pim[ll+m*2-2]+ ff*this->pim[ll+m*2+2]);
				this->py[ii+m*2]   -=						 B*ff*this->pim[ll+m*2+1];
				this->pz[ii+m*2-1] -= B*(i+m+1)*this->pim[ll+m*2-1];
				this->pz[ii+m*2]   -= B*(i+m+1)*this->pim[ll+m*2];
			}
			this->px[ii] = -A*.5*i*(i-1.)*(this->p[jj+2]-this->pim[jj+2])-B*.5*(i+1.)*(i+2.)*this->pim[ll+2];
			this->py[ii] = -A*.5*i*(i-1.)*(this->p[jj+1]-this->pim[jj+1])-B*.5*(i+1.)*(i+2.)*this->pim[ll+1];
			this->pz[ii] =  A*i*(this->p[jj]-this->pim[jj])-B*(i+1)*this->pim[ll];
		}
		if (this->id_inverse > 0 && fabs(imag(this->E1[1])) > EE) B = imag(this->E1[0])/(imag(this->E1[1])*this->param[4])*this->R0; else
		if (this->id_inverse < 0 && fabs(imag(this->E1[1])) > EE) B = imag(this->E1[0])/(imag(this->E1[1])*this->param[3])*this->R0; else B = 0.;
		this->px[0] = -B*this->pim[3];
		this->py[0] = -B*this->pim[2];
		this->pz[0] = -B*this->pim[1];
		this->N -= m_dop;
	}
}

//////////////////////////////////////////////////////////////////////////
//...parametrization hessian of the multipoles (гессиан не реализован!!!);
template <typename T>
void CMapi2DClayer<T>::parametrization_hess(double * P, int m_dop)
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

///////////////////////////////////////////
//...calculation hessian of the multipoles;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi2DClayer<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->px) return;
	for (int m = this->N; m > 0; m--) {
		deriv[m*2-1] += this->px[m*2-1]*f;
		deriv[m*2]   += this->px[m*2]*f;
	}
}

template <typename T>
void CMapi3DClayer<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->px) return;
	for (int m = freedom(this->N)-1; m >= 0; m--)
		deriv[m] += this->px[m]*f;
}

template <typename T>
void CMapi2DClayer<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->py) return;
	for (int m = this->N; m > 0; m--) {
		deriv[m*2-1] += this->py[m*2-1]*f;
		deriv[m*2]   += this->py[m*2]*f;
	}
}

template <typename T>
void CMapi3DClayer<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->py) return;
	for (int m = freedom(this->N)-1; m >= 0; m--)
		deriv[m] += this->py[m]*f;
}

template <typename T>
void CMapi3DClayer<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->pz) return;
	for (int m = freedom(this->N)-1; m >= 0; m--)
		deriv[m] += this->pz[m]*f;
}
#endif
