/*========================================*/
/*                 CBEAM                  */
/*========================================*/
#ifndef ___CBEAM___
#define ___CBEAM___

#include "cshapes.h"

/////////////////////////////////////////////////////
//																	//
//       BEAM AND ROUNDED WEDGE MULTIPOLES         //
//																	//
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//...class of beams polynomial multipoles (sopromat);
template <typename T>
class CBeamShape : public CShape<T> {
public:
		Num_Shape	 type() { return BEAM_POLY_SHAPE;}
		int size_of_param() { return(11);};
		int freedom (int m) { return(4); };
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_Z(T * deriv, double f = 1.);
//...constructor;
		CBeamShape() {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			this->param[0] = this->param[1] = this->param[3] = this->param[4] = 
			this->param[5] = this->param[6] = 1.;
			this->param[2] = .3;
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE__E            2.00109e11 //...Yung modulus of the section;
#undef  SHAPE__E								//...param(0);
#define SHAPE__G							1. //...shear modulus of the section;
#undef  SHAPE__G								//...param(1);
#define SHAPE__nju					  0.3 //...Poisson coefficient of the section;
#undef  SHAPE__nju							//...param(2);
#define SHAPE__GF_inv					1. //...inverse tension modulus of the section;
#undef  SHAPE__GF_inv						//...param(3);
#define SHAPE__EIy_inv					1. //...inverse X-bending modulus of the section;
#undef  SHAPE__EIy_inv						//...param(4);
#define SHAPE__EIx_inv					1. //...inverse Y-bending modulus of the section;
#undef  SHAPE__EIx_inv						//...param(5);
#define SHAPE__GJ_inv					1. //...inverse torsional modulus of the section;
#undef  SHAPE__GJ_inv						//...param(6);
#define SHAPE__mx_J						0. //...part of torsional modulus by the X-bending;
#undef  SHAPE__mx_J							//...param(7);
#define SHAPE__my_J						0. //...part of torsional modulus by the Y-bending;
#undef  SHAPE__my_J							//...param(8);
#define SHAPE__cx_I						0. //...part of X-bending by the torsion;
#undef  SHAPE__cx_I							//...param(9);
#define SHAPE__cy_I						0. //...part of Y-bending by the torsion;
#undef  SHAPE__cy_I							//...param(10);

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CBeamShape<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N = 4))*sizeof(T));
	}
	if (P && this->p) {
		double Z = P[2]*this->R0_inv;
		this->cs[0] = 0.;
		this->cs[1] = 0.;
		this->cs[2] = 1.;

///////////////////////////////////
//...calculation of the multipoles;
		int m;
		for ( this->p[0] = 1., m = 1; m < this->N; m++) {
				this->p[m] = this->p[m-1]*Z;
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CBeamShape<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! this->pz) {
		this->pz = (T *)new_struct((this->NN_grad = this->NN_dop)*sizeof(T));
	}
	if (P && this->pz) {
		memset (this->pz, 0, this->NN_grad*sizeof(T));
		deriv_Z(this->pz);
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CBeamShape<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P);
	if (! this->pzz) {
		this->pzz = (T *)new_struct((this->NN_hess = this->NN_grad)*sizeof(T));
	}
	if (P && this->pzz) {
		memset (this->pzz, 0, this->NN_hess*sizeof(T));
		swap(this->p, this->pz);
			deriv_Z(this->pzz);
		swap(this->p, this->pz);
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CBeamShape<T>::deriv_Z(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	for (int i = this->NN-1; i >= 0; i--) {
		deriv[i] += i*this->p[i-1]*f;
	}
}

////////////////////////////////////////////////////////
//															         //
//       CLASSICAL MULTIPOLES of ROUNDED WEDGE        //
//															         //
////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//...class of classical multipoles by conformal mapping of rounded wedge;
template <typename T>
class CMapi2DCorner : public CShape<T> {
protected:
		double * b, * B;
		double alpha, sigma, l, L, test1, test2, test3;//...internal paramters of conformal maqpping;
		complex Conform, DConform, DDConform;			  //...conformal maqpping and derivatives;
		void   mapping();
		double gamma(double z);
		static const double c0[15];
public:
		Num_Shape    type() { return MP2D_CORNER_SHAPE;}
		int size_of_param() { return(4);};
		int freedom (int m) { return m*2+1; }
public:
//...initialization of multipoles;
		void init1(int N, int dim) { init3(N, 0, 0, dim);}
		void init3(int N, int N1, int N2, int dim) {
			CShape<T>::N  = N;
			CShape<T>::N1 = N1;
			CShape<T>::N2 = N2;
			CShape<T>::id_dim = dim;
			this->NN = freedom(N);
			release();
			//this->NN = (CMapi2DCorner::N1*2+1)*(CMapi2DCorner::N2*2+1);
		}
		void set_shape(double R0 = 0., double beta = 0., double radc = 0., double L1 = 0., double L2 = 0., double L3 = 0.) {
			CMapi2DCorner::R0 = R0;
			this->R0_inv  = R0 > EE ? 1./R0 : 1.;
			int   k = -1;
			if (++k < size_of_param()) this->param[k] = (Param)(1./beta);
			if (++k < size_of_param()) this->param[k] = (Param)radc;
			if (++k < size_of_param()) this->param[k] = (Param)(radc/cos(beta*M_PI_2));
			if (++k < size_of_param()) this->param[k] = (Param)(fabs(radc*tan(beta*M_PI_2)));
		}
		void release() {
			delete_struct(b);
			delete_struct(B);
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
		CMapi2DCorner() {
			this->param = (Param *)new_struct(size_of_param()*sizeof(Param));
			b = B = NULL;
			alpha = sigma = l = L = test1 = test2 = test3 = 0.;
			Conform = DConform = DDConform = comp(0.);
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_corner_inv_beta           1. //...inverse of opening angle;
#undef  SHAPE_corner_inv_beta              //...param(0);
#define SHAPE_corner_radc               0. //...rounding radius;
#undef  SHAPE_corner_radc                  //...param(1);
#define SHAPE_corner_O1                 0. //...origine of rounding circle;
#undef  SHAPE_corner_O1                    //...param(2);
#define SHAPE_corner_R0                 0. //...radius of rounding zone;
#undef  SHAPE_corner_R0                    //...param(3);

///////////////////////////////////////////////////////////////////////////////////
//...коэффициенты в разложениях по полиномам Чебышева для вычисления гамма-функции;
template <typename T> const double CMapi2DCorner<T>::c0[15] = {
3.65738772508338243850,  1.95754345666126826928,  3.3829711382616038916E-1,
  4.208951276557549199E-2,  4.28765048212908770E-3,  3.6521216929461767E-4,
     2.740064222642200E-5,     1.81240233365124E-6,     1.0965775865997E-7,
	 5.98718404552E-9,         3.0769080535E-10,        1.431793030E-11,
	     6.5108773E-13,            2.595850E-14,            1.10789E-15};

////////////////////////////////////////
//...процедура вычисления гамма-функции;
template <typename T>
double CMapi2DCorner<T>::gamma(double z)
{
  double sum = c0[0], f = 1., t0, t1 = 1., t2;
  while (z < 3.) {
         f /= z; z += 1.;
  }
  while (z > 4.) {
         z -= 1.; f *= z;
  }
  t2 = 2.*(z-3.)-1.; z = 2.*t2;

////////////////////////////////////////////////////////////////
//...суммиpуем pяд из полиномов Чебышева и возвpащаем pезультат;
  for (int i = 1; i <= 14; i++) {
       sum += t2*c0[i]; t0 = t2; t2 = z*t2-t1; t1 = t0;
  }
  return(f*sum);
}

//////////////////////////////////////////////////////////////////////////////////////
//...процедура вычисления коэффициентов для конфоpмного отобpажения уголковой области;
template <typename T>
void CMapi2DCorner<T>::mapping()
{
	double beta = 1./this->param[0];
	if (beta <= 0. || beta >= 2.0) return;

	b[0] = B[0] = 1.;
	if (beta == 1.) return;

///////////////////////////////////////
//...инициализируем рабочие переменные;
	double nju = .5*(beta-1.), g1 = .5*(gamma(.5-.5*nju)/gamma(1.-.5*nju))*
												  (gamma(1.+.5*nju)/gamma(1.5+.5*nju)),
		g2 = 2.*pow((gamma(1.-nju)/gamma(.5-nju))*(gamma(1.5+nju)/gamma(2.+nju)), (1./beta)),
		r1 = g1*g1, r2 = g2*g2, q = 4./beta*(beta+2.)/(beta-1.)/(beta+1.)*.8, fm = 1.-2.*nju,
		f1, f2, fr = 1., fk;
	int  k, j, m, kk;

//////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем скалярные параметры структуры (параметр деформации l может быть отриц. числом);
	L = pow(fabs(tan(nju*M_PI)), (1./beta));
	l = 1./tan(nju*M_PI_2);
	alpha = L*g2; sigma = alpha*g1;

/////////////////////////////////////////////////////
//...распределяем в общей памяти два рабочих массива;
	double * S = (double *)new_struct((2+max(this->N1, this->N2*(this->N2+1)/2))*sizeof(double)),
			 * C = (double *)new_struct((2+max(this->N1, this->N2))*sizeof(double));

/////////////////////////////////////////////////////////////////////
//...делим два гипергеометрических ряда (для внутреннего разложения);
	f1 = f2 = S[0] = 1.;
	for (k = 1; k <= this->N1; k++) {
		f1 *= ((2.*k-2.+nju)/k)*.5*((2.*k-3.-nju)/(2.*k-1.)); S[k] = f1;
		f2 *= ((2.*k-2.-nju)/k)*.5*((2.*k-1.+nju)/(2.*k+1.)); C[k] = f2;
		for (j = 1; j <= k; j++) S[k] -= C[j]*S[k-j];
	}

/////////////////////////////////////////////////////
//...возводим полученный результат во вторую степень;
	for (k = 1; k <= this->N1; k++) {
		C[k] = S[k];
		for (m = 1; m <= k; m++) C[k] += S[m]*S[k-m];
	}

///////////////////////////////////
//...обращаем ряд Бюрмана-Лагранжа;
	for (j = 1; j <= this->N1; j++) {
		for (k = this->N1; k >= 0; k--) {
			for (m = 1; m <= k; m++) S[k] += C[m]*S[k-m];
			S[k] *= r1;
		}
		b[j] = S[j]/(2.*j+1.); test1 += b[j];
	}
	test1 += (1.-1./g1);

/*=================================================================================*/
/*      Повторяем процедуру для вычисления коэффициентов внешнего разложения.      */
/*=================================================================================*/
////////////////////////////////////////
//...делим два гипергеометрических ряда;
	f1 = .5*nju*(nju+1.); f2 = f1/(3.+2.*nju); S[1] = (f1-fm*f2)*q; C[1] = f2;
	for (k = 2; k <= this->N2; k++) {
		f1 *= ((2.*k-2.-nju)/k)*.5*((2.*k-3.-nju)/(2.*k-1.-2.*nju));
		f2 *= ((2.*k-2.+nju)/k)*.5*((2.*k-1.+nju)/(2.*k+1.+2.*nju)); C[k] = f2;
		kk = k*(k-1)/2; S[kk+1] = (f1-fm*f2)*q;
		for (j = 1; j <= k-1; j++) S[kk+1] -= S[j*(j-1)/2+1]*C[k-j];
	}

///////////////////////////////////
//...обращаем ряд Бюрмана-Лагранжа;
	for (k = 1; k <= this->N2; k++) {
		fk = (2.*k-1.)/beta; fr *= r2; kk = k*(k-1)/2;
		f1 = (fr/fm)*(fk/q);
		f2 = S[kk+1]*f1;
		for (m = 2; m <= k; m++) {
			f1 *= (fk-m+1.)/fm/(m*q); S[kk+m] = 0.;
			for (j = m-1; j <= k-1; j++) S[kk+m] += S[j*(j-1)/2+m-1]*S[(k-j-1)*(k-j)/2+1];
			f2 += S[kk+m]*f1;
		}
		B[k] = -f2/(2.*k-1.); test2 += B[k];
	}
	test2 += 1.-g2; test3 = test2*L-test1*sigma;

///////////////////////////////////////////////////////////////////
//...освобождаем память от рабочих массивов и выходим из пpогpаммы;
	delete_struct(C); 
	delete_struct(S);
}
/*============================================================================================*/
/*  Для внутренней области G1 представление для конформного отображения имеет следующий вид:  */
/*                                                                       		                */
/*           z=F(w) = sigma*<сумма k=0,EE>[ b[k]*(i*l*(W-N)/(W-M1))**(2k+1) ],    	          */
/*                                                                       		                */
/*  где W = w*exp(-i*pi*beta/2), N и M1 - точки на окружности.           		                */
/*                                                                       		                */
/*  Для внешней области G2 представление для конформного отображения имеет следующий вид:     */
/*                                                                       		                */
/*           z=F(w) = w**(1/beta)*<сумма k=0,EE>*[ B[k]*(L/w**(1/beta))**(2k) ].  	          */
/*============================================================================================*/

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CMapi2DCorner<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
			this->p  = (T *)new_struct((this->NN_dop = freedom(this->N))*sizeof(T));
/////////////////////////////////////////////////////////////////////
//...technical support (calculation of the conjunction coeffiсients);
			delete_struct(b); b = (double  *)new_struct((this->N1+1)*sizeof(double));
			delete_struct(B); B = (double  *)new_struct((this->N2+1)*sizeof(double));
			mapping();
	}
	if (P && this->p) {
		complex z = comp(P[0], P[1]), w0;
		T      zm = 1., zm_i = 0., hh;
		double ff = sigma*pow(fabs(this->param[1])*this->R0_inv, this->param[0]),
				 LL = pow(this->param[3]*this->R0_inv, this->param[0]), R = abs(z), fi = arg2(z);
		this->cs[0] = P[3];
		this->cs[1] = P[4];
		this->cs[2] = 0.;

//////////////////////////////////////////////////
//...конфоpмное отобpажение во внутpенней области;
		Conform = comp(0.);

		if (R < this->param[3]) {
			z -= this->param[2]; Conform = (z-this->param[1])/(z+this->param[1])*comp(0., l); w0 = Conform*Conform;
			z = comp(0.);

			for (int j = this->N1; j >= 1; j--) z = b[j]+w0*z;
			Conform *= (1.+w0*z)*ff;
		}
		else

///////////////////////////////////////////////
//...конфоpмное отобpажение во внешней области;
		if (R >= EE) {
			Conform = polar(pow(R*this->R0_inv, this->param[0]), fi*this->param[0])*comp(0., 1.); 

			if (this->param[3] > 0.) {
				w0 = LL/Conform; w0 *= w0; z = comp(0.);

				for (int j = this->N2; j >= 1; j--) z = B[j]+w0*z;
				Conform *= (1.+w0*z);
			}
		}

///////////////////////////////////
//...calculation of the multipoles;
		int m;
		for ( this->p[0] = zm, m = 1; m <= this->N; m++) {
				this->p[m*2-1] = (zm = (hh = zm)*imag(Conform)+zm_i*real(Conform));
				this->p[m*2]   = (zm_i = hh*real(Conform)-zm_i*imag(Conform));
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CMapi2DCorner<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P);
	if (! this->px && ! this->py) {
		this->px = (T *)new_struct((this->NN_grad = this->NN_dop)*sizeof(T));
		this->py = (T *)new_struct( this->NN_grad*sizeof(T));
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));

		complex z = comp(P[0], P[1]), w0, w;
		double ff = sigma*pow(fabs(this->param[1])*this->R0_inv, this->param[0]),
				 LL = pow(this->param[3]*this->R0_inv, this->param[0]), R = abs(z), fi = arg2(z);

///////////////////////////////////////////////////////////////
//...пpоизводная конфоpмного отобpажения во внутpенней области;
		DConform = comp(0.);

		if (R < this->param[3]) {
			z -= this->param[2]; w0 = (z-this->param[1])/(z+this->param[1])*comp(0., l); w0 *= w0;
			DConform = comp(0., 2.*this->param[1]*l)/((z+this->param[1])*(z+this->param[1])); z = comp(0.);

			for (int j = this->N1; j >= 1; j--) z = (2.*j+1.)*b[j]+w0*z;
			DConform *= (1.+w0*z)*ff;
		}
		else

////////////////////////////////////////////////////////////
//...пpоизводная конфоpмного отобpажения во внешней области;
		if (this->param[0] == 1. && R < EE) DConform = comp(0., this->R0_inv); 	else 
		if (R >= EE) {
			DConform = (w = polar(pow(R*this->R0_inv, this->param[0]), fi*this->param[0])*comp(0., 1.))*this->param[0]/z; 

			if (this->param[3] > 0.) {
				w0 = LL/w; w0 *= w0; z = comp(0.);

				for (int j = this->N2; j >= 1; j--) z = (2.*j-1.)*B[j]+w0*z;
				DConform *= (1.-w0*z);
			}
		}
		DConform *= comp(this->CZ, -this->SZ); //...переводим в общую систему координат!!!

////////////////////////////////////////
//...calculation grad of the multipoles;
		int m;
		for ( this->px[0] = this->py[0] = 0., m = this->N; m > 1; m--) {
				this->px[m*2-1] = m*(this->p[m*2-2]*imag(DConform)+this->p[m*2-3]*real(DConform));
				this->px[m*2]	 = m*(this->p[m*2-2]*real(DConform)-this->p[m*2-3]*imag(DConform));
				this->py[m*2-1] =  this->px[m*2];
				this->py[m*2]   = -this->px[m*2-1];
		}
		if (m) {
				this->px[m*2-1] = imag(DConform);
				this->px[m*2]   = real(DConform);
				this->py[m*2-1] =  this->px[m*2];
				this->py[m*2]   = -this->px[m*2-1];
		}
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CMapi2DCorner<T>::parametrization_hess(double * P, int m_dop)
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

		complex z = comp(P[0], P[1]), w0, nm, w;
		double ff = sigma*pow(fabs(this->param[1])*this->R0_inv, this->param[0]),
				 LL = pow(this->param[3]*this->R0_inv, this->param[0]), R = abs(z), fi = arg2(z);

//////////////////////////////////////////////////////////////////////
//...вторая пpоизводная конфоpмного отобpажения во внутpенней области;
		DDConform = comp(0.);

		int j, m;
		if (R < this->param[3]) {
			z -= this->param[2]; w = (z-this->param[1])/(z+this->param[1])*comp(0., l); w0 = w*w;
			DDConform = (nm = comp(0., 2.*this->param[1]*l)/((z+this->param[1])*(z+this->param[1])))*(-2./(z+this->param[1]));
			nm *= nm; z = comp(0.);

			for (j = this->N1; j >= 1; j--) z = (2.*j+1.)*(2.*j)*b[j]+w0*z;
			nm *= w*z*ff; z = comp(0.);

			for (j = this->N1; j >= 1; j--) z = (2.*j+1.)*b[j]+w0*z;
			DDConform *= (1.+w0*z)*ff;
			DDConform += nm;
		}
		else

///////////////////////////////////////////////////////////////////
//...вторая пpоизводная конфоpмного отобpажения во внешней области;
		if (R >= EE) {
			DDConform  = (w = polar(pow(R*this->R0_inv, this->param[0]), fi*this->param[0])*comp(0., 1.))*(nm = 1./z);
			DDConform *= (this->param[0]-1.)*(nm *= this->param[0]);

			if (this->param[3] > 0.) {
				nm *= nm*w; w0 = LL/w; w0 *= w0; z = comp(0.);

				for (j = this->N2; j >= 1; j--) z = (2.*j-1.)*(2.*j)*B[j]+w0*z;
				nm *= w0*z; z = comp(0.);

				for (j = this->N2; j >= 1; j--) z = (2.*j-1.)*B[j]+w0*z;
				DDConform *= (1.-w0*z);
				DDConform += nm;
			}
		}
		DDConform *= comp(sqr(this->CZ)-sqr(this->SZ), -2.*this->CZ*this->SZ); //...переводим в общую систему координат!!!

///////////////////////////////////////////
//...calculation hessian of the multipoles;
		for ( this->pxx[0] = this->pxy[0] = this->pyy[0] = 0., m = this->N; m > 1; m--) {
				this->pxx[m*2-1] = m*((this->px[m*2-2]*imag(DConform)+this->px[m*2-3]*real(DConform))+(this->p[m*2-2]*imag(DDConform)+this->p[m*2-3]*real(DDConform)));
				this->pxx[m*2]   = m*((this->px[m*2-2]*real(DConform)-this->px[m*2-3]*imag(DConform))+(this->p[m*2-2]*real(DDConform)-this->p[m*2-3]*imag(DDConform)));
				this->pxy[m*2-1] =  this->pxx[m*2];
				this->pxy[m*2]   = -this->pxx[m*2-1];
				this->pyy[m*2-1] = -this->pxx[m*2-1];
				this->pyy[m*2]   = -this->pxx[m*2];
		}
		if (m) {
				this->pxx[m*2-1] = imag(DDConform);
				this->pxx[m*2]   = real(DDConform);
				this->pxy[m*2-1] =  this->pxx[m*2];
				this->pxy[m*2]   = -this->pxx[m*2-1];
				this->pyy[m*2-1] = -this->pxx[m*2-1];
				this->pyy[m*2]   = -this->pxx[m*2];
		}
/////////////////////
//	cpy(this->pxx, p_cpy);//
/////////////////////
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi2DCorner<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->px) return;

	for (int m = this->N; m > 0; m--) {
		deriv[m*2-1] += this->px[m*2-1]*f;
		deriv[m*2]   += this->px[m*2]*f;
	}
	deriv[0] += this->px[0]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CMapi2DCorner<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->py) return;

	for (int m = this->N; m > 0; m--) {
		deriv[m*2-1] += this->py[m*2-1]*f;
		deriv[m*2]   += this->py[m*2]*f;
	}
	deriv[0] += this->py[0]*f;
}
#endif
