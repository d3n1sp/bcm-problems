/*==========================================*/
/*                 csmixer.h                  */
/*==========================================*/
#ifndef ___CSMIXER___
#define ___CSMIXER___

#include <typeinfo>

#include "cgrid.h"
#include "cshapes.h"

//////////////////////////////////////////////////////
//...class of the combination various shape functions;
template <typename T>
class CShapeMixer {
protected:
		int			N_shape, nm;							//...number of any type and usual shapes;
		CShape<T> ** shape;									//...list of used shape functions;
		double CX, CY, CZ, SX, SY, SZ, P0[3], cs[6];	//...local coordinate system;
public:
		T ** A;							//...saved potentials;
		int NN_local, NN, n;			//...number of all freedom, saved potentials, classical shapes;
		int num_shape() {return N_shape;}
		int num_usual() {return nm;}
public:
		static T * deriv, * p_cpy;	//...derivative and function (technical arrays);
		static T * d_cmp, * p_cmp;	//...another part (real or imag) of multipokes (technical arrays);
		static int N_mx;				//...max local freedom of technical arrays;
public:
		Num_Shape type	(int k = 0) { return (shape && k < N_shape ? shape[k]->type() : NULL_SHAPE);}
		int freedom	   (Num_State id_full = NULL_STATE);
public:
//...initialization of multipoles;
		void add_shape(		 CShape<T> * sp, int usual = 1);
		void add_shape(int k, CShape<T> * sp, int usual = 1);
		void init1(			int N, int m, int dim);
		void init1(int k, int N, int m, int dim);
		void init2(			int N, int M, int m, int dim);
		void init2(int k, int N, int M, int m, int dim);
		void init3(			int N, int N1, int N2, int m, int dim);
		void init3(int k, int N, int N1, int N2, int m, int dim);
		void release();
		void release(int k) { if (shape && k < N_shape) shape[k]->release();}
		void set_shape(		 double R0, double kk = 0., double p1 = 0., double p2 = 0., double p3 = 0., double p4 = 0.);
		void set_shape(int k, double R0, double kk = 0., double p1 = 0., double p2 = 0., double p3 = 0., double p4 = 0.);
//...setting of coordinate system;
		void set_local (double * P = NULL) {
			if (P) { 
				CZ = cos(P[3]); SZ = sin(P[3]);
				CY = cos(P[4]); SY = sin(P[4]);
				CX = cos(P[5]); SX = sin(P[5]); memcpy(P0, P, 3*sizeof(double));
			}
			else {
				CZ = CY = CX = 1.;
				SZ = SY = SX = 0.; memset(P0, 0, 3*sizeof(double));
			}
		}
		void set_local	  (int k, double * P = NULL) { if (shape && k < N_shape) shape[k]->set_local(P);}
		void set_local_P0(		 double * P) { if (P) memcpy(P0, P, 3*sizeof(double));}
		void set_local_P0(int k, double * P) {	if (shape && k < N_shape) shape[k]->set_local_P0(P);}
		void make_local (double * P) {
			point_shift<double>(P, P0, OK_STATE);
			point_iso<double>  (P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
		}
		void make_common(double  * P) { point_iso<double>(P, P0, CZ, SZ, CY, SY, CX, SX);}
		void norm_local (double  * P) { point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);}
		void norm_common(double  * P) { point_iso<double>(P, NULL, CZ,  SZ, CY,  SY, CX,  SX);}
		void set_norm_cs(double  * P)	{ cs[0] = P[3]; cs[1] = P[4]; cs[2] = P[5];}
		void norm_local_T (T * P) { point_iso<T>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);}
		void norm_common_T(T * P) { point_iso<T>(P, NULL, CZ,  SZ, CY,  SY, CX,  SX);}
//...parameters of shape functions;
		CShape<T> * get_shape(int k = 0) { return (shape ? shape[min(N_shape-1, k)] : NULL);}
		int	 get_NN		(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->NN : 0);}
		int	 get_N		(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->N : 0);}
		double get_R		(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->R0 : 1.);}
		double get_R_inv	(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->R0_inv : 1.);}
		int	 cmpl			(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->cmpl() : 0);}
		int	 inverse		(int k = 0) {	return(shape ? shape[min(N_shape-1, k)]->inverse() : 0);}
		void	 change		(int k = 0) { if (shape && k < N_shape) shape[k]->change();}
		void	 change_cmpl(double *& pp, double *& dd) { if (cmpl() == 1 && dd)  swap(pp, dd);};
		void	 set_cmpl	(int k, int m_cmpl) { if (shape && k < N_shape) shape[k]->set_cmpl(m_cmpl);}
		void	 set_cmpl	(int m_cmpl = 0) { set_cmpl(0, m_cmpl);}
		void   set_R		(			double value); //...устанавливаем на все;
		void   set_R		(int k,	double value);
		void   set_param	(		   int m, double value) { set_param(0, m, value);}
		void   set_param	(int k, int m, double value) { if (shape && k < N_shape) shape[k]->set_param(m, value);}
		double get_param	(		   int m) { return get_param(0, m);}
		double get_param	(int k, int m) { return(shape ? shape[min(N_shape-1, k)]->get_param(m) : 0.);}
//...calculation shape functions;
		void parametrization     (			double * P = NULL, int m = 0); //...вычисляем все;
		void parametrization     (int k, double * P = NULL, int m = 0);
		void parametrization_grad(			double * P = NULL, int m = 0);
		void parametrization_grad(int k, double * P = NULL, int m = 0);
		void parametrization_hess(			double * P = NULL, int m = 0);
		void parametrization_hess(int k, double * P = NULL, int m = 0);
//...differentiation of shape functions;
		void deriv_X (			T * deriv, double f = 1.); //...дифференцируем все;
		void deriv_X (int k, T * deriv, double f = 1.);
		void deriv_Y (			T * deriv, double f = 1.);
		void deriv_Y (int k, T * deriv, double f = 1.);
		void deriv_Z (			T * deriv, double f = 1.);
		void deriv_Z (int k, T * deriv, double f = 1.);
		void deriv_t (			T * deriv, double f = 1.);
		void deriv_t (int k, T * deriv, double f = 1.);
		void deriv_N ();
		void deriv_N (int k);
//...integration of multipoles;
		void primitive(		 T * deriv); //...интегрируем все;
		void primitive(int k, T * deriv);
//...setting potentials;
		void   init_potential();
		void delete_potential();
		void get_potential(T * AA, int m) {	if (A && AA && m < n) memcpy(AA, A[m], NN*sizeof(T));}
		void set_potential(T * AA, int m) {	if (A && AA && m < n) memcpy(A[m], AA, NN*sizeof(T));}
//...calculation potentials;
		T potential(		 T * p_cpy, int m);
		T potential(int k, T * p_cpy, int m);
//...constructor and destructor;
		CShapeMixer() {
			N_shape = nm = N_mx = NN_local = NN = n = 0;
			A		  = NULL;
			shape   = NULL;
			set_local();
		}
		virtual ~CShapeMixer(void);
//...auxilliary functions;
		T * FULL(T * p_ext, int k = 0, int l = 0, int shift = 1); 
//...copy;
		void cpy	  (		 T * p_ext);
		void cpy	  (int k, T * p_ext);
		void cpy	  (		 T * pp, T * p_ext) { if (p_ext && pp) memcpy(p_ext, pp, NN_local*sizeof(T));}
		void cpy	  (int k, T * pp, T * p_ext);
		void cpy_x (		 T * p_ext = NULL);
		void cpy_x (int k, T * p_ext = NULL);
		void cpy_y (		 T * p_ext = NULL);
		void cpy_y (int k, T * p_ext = NULL);
		void cpy_z (		 T * p_ext = NULL);
		void cpy_z (int k, T * p_ext = NULL);
		void cpy_xx(		 T * p_ext = NULL);
		void cpy_xx(int k, T * p_ext = NULL);
		void cpy_xy(		 T * p_ext = NULL);
		void cpy_xy(int k, T * p_ext = NULL);
		void cpy_yy(		 T * p_ext = NULL);
		void cpy_yy(int k, T * p_ext = NULL);
		void cpy_xz(		 T * p_ext = NULL);
		void cpy_xz(int k, T * p_ext = NULL);
		void cpy_yz(		 T * p_ext = NULL);
		void cpy_yz(int k, T * p_ext = NULL);
		void cpy_zz(		 T * p_ext = NULL);
		void cpy_zz(int k, T * p_ext = NULL);
//...linear combination;
		void adm	  (		 T * p_ext, T adm);
		void adm	  (int k, T * p_ext, T adm);
		void adm_x (		 T * p_ext, T adm);
		void adm_x (int k, T * p_ext, T adm);
		void adm_y (		 T * p_ext, T adm);
		void adm_y (int k, T * p_ext, T adm);
		void adm_z (		 T * p_ext, T adm);
		void adm_z (int k, T * p_ext, T adm);
		void adm_xx(		 T * p_ext, T adm);
		void adm_xx(int k, T * p_ext, T adm);
		void adm_xy(		 T * p_ext, T adm);
		void adm_xy(int k, T * p_ext, T adm);
		void adm_yy(		 T * p_ext, T adm);
		void adm_yy(int k, T * p_ext, T adm);
		void adm_xz(		 T * p_ext, T adm);
		void adm_xz(int k, T * p_ext, T adm);
		void adm_yz(		 T * p_ext, T adm);
		void adm_yz(int k, T * p_ext, T adm);
		void adm_zz(		 T * p_ext, T adm);
		void adm_zz(int k, T * p_ext, T adm);
//...admittance combination;
		void admittance(		  T adm);
		void admittance(int k, T adm);
		void admittance(		  T * dd, T * pp, T adm_dd, T adm_pp);
		void admittance(int k, T * dd, T * pp, T adm_dd, T adm_pp);
//...testing functions;
		void TestShape (FILE * TST, int i, int j, double eps, int i_dop);
		void TestShape (char * ch_POTEN, int i = -1, int j = -1, double eps = EE_ker, int i_dop = -1);
		void TestShape (const char * ch_POTEN, int i = -1, int j = -1, double eps = EE_ker, int i_dop = -1);
		void GetSurferFormat(FILE * SURF1, FILE * SURF2, double * par, int id_axis);
		void GetSurferFormat(char * SURF_FILE, double * par, int id_axis);
		void GetSurferFormat(const char * SURF_FILE, double * par, int id_axis);
};
//////////////////////////////////////////
//...complete deleting of abstract shapes;
template <typename T>
inline void delete_shapes(CShapeMixer<T> *& sp) {
	delete sp; sp = NULL;
}

////////////////////////////////////////////////////
//          TEMPLATE VARIANT OF CMIXER            //
////////////////////////////////////////////////////
template <typename T> int CShapeMixer<T>::N_mx  = 0;
template <typename T> T * CShapeMixer<T>::deriv = NULL;
template <typename T> T * CShapeMixer<T>::p_cpy = NULL;

//////////////////////////////////////////
//...destructor of the class CShapeMixer;
template <typename T>
CShapeMixer<T>::~CShapeMixer(void)
{
	delete_potential();
   if (shape)
		for (int k = 0; k < N_shape; k++) 
			delete shape[k];
	delete_struct(shape);
}

/////////////////////////////////
//...mixer of shapes composition;
template <typename T>
void CShapeMixer<T>::add_shape(CShape<T> * shape, int usual)
{
	CShape<T> ** new_shape = (CShape<T> **)new_struct(++N_shape*sizeof(CShape<T> *));
	memcpy (new_shape, this->shape, (N_shape-1)*sizeof(CShape<T> *));
	delete [] this->shape;
	this->shape = new_shape;
	this->shape[N_shape-1] = shape;
	nm += usual;
}
template <typename T>
void CShapeMixer<T>::add_shape(int k, CShape<T> * shape, int usual)
{
	if (shape && k < N_shape) {
		delete_potential();
		delete shape[k]; shape[k] = shape;
		nm += usual;
	}
}

///////////////////////////////
//...initialization potentials;
template <typename T>
void CShapeMixer<T>::init_potential()
{
	if (N_mx < NN_local) {
		delete_struct(deriv); deriv = (T *)new_struct(NN_local*sizeof(T));
		delete_struct(p_cpy); p_cpy = (T *)new_struct(NN_local*sizeof(T)); N_mx = NN_local;
	}
	set_matrix(A, n, NN);
}
template <typename T>
void CShapeMixer<T>::delete_potential()
{
	release();
	delete_struct(deriv);
	delete_struct(p_cpy);
	delete_struct(A);	N_mx = NN_local = NN = n = 0;
}

//////////////////////////
//...auxilliary functions;
template <typename T>
int CShapeMixer<T>::freedom(Num_State id_full)
{
	int k, NN_local = 0;
	if (shape) {
		if (id_full == NULL_STATE)
		for (k = 0; k < N_shape; k++)	NN_local += shape[k]->NN;	else
		for (k = 0; k < N_shape; k++)	NN_local += shape[k]->NN*shape[k]->dim();
	}
	return NN_local;
}

/////////////////////////////////
//...initialization of the shape;
template <typename T>
void CShapeMixer<T>::init1(int N, int m, int dim)
{
	this->n = m;
   if (shape)
		for (int k = 0; k < N_shape; k++)
			shape[k]->init1(N, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::init1(int k, int N, int m, int dim)
{
	this->n = m;
	if (shape && k < N_shape) shape[k]->init1(N, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::init2(int N, int M, int m, int dim)
{
	this->n = m;
   if (shape)
		for (int k = 0; k < N_shape; k++)
			shape[k]->init2(N, M, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::init2(int k, int N, int M, int m, int dim)
{
	this->n = m;
   if (shape && k < N_shape) shape[k]->init2(N, M, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::init3(int N, int N1, int N2, int m, int dim)
{
	this->n = m;
   if (shape)
		for (int k = 0; k < N_shape; k++)
			shape[k]->init3(N, N1, N2, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::init3(int k, int N, int N1, int N2, int m, int dim)
{
	this->n = m;
   if (shape && k < N_shape) shape[k]->init3(N, N1, N2, dim);
	NN			= freedom(OK_STATE); 
	NN_local = freedom(NULL_STATE); 
}
template <typename T>
void CShapeMixer<T>::release()
{
   if (shape)
		for (int k = 0; k < N_shape; k++) 
			shape[k]->release();
}

/////////////////////////////////////
//...setting of the shape parameters;
template <typename T>
void CShapeMixer<T>::set_shape(double R0, double kk, double p1, double p2, double p3, double p4)
{
   if (shape)
		for (int k = 0; k < N_shape; k++) 
			shape[k]->set_shape(R0, kk, p1, p2, p3, p4);
}
template <typename T>
void CShapeMixer<T>::set_shape(int k, double R0, double kk, double p1, double p2, double p3, double p4)
{
   if (shape && k < N_shape) shape[k]->set_shape(R0, kk, p1, p2, p3, p4);
}

///////////////////////////////////
//...parameters of shape functions;
template <typename T>
void CShapeMixer<T>::set_R(double value)
{
	if (shape)
		for (int k = 0;  k < N_shape; k++) {
			shape[k]->R0		= value;
			shape[k]->R0_inv = value > EE ? 1./value : 1.;
		}
}
template <typename T>
void CShapeMixer<T>::set_R(int k, double value)
{
	if (shape && k < N_shape) {
		 shape[k]->R0		= value;
		 shape[k]->R0_inv = value > EE ? 1./value : 1.;
	}
}

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CShapeMixer<T>::parametrization(double * P, int m)
{
	if (P && shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		for (int k = 0; k < N_shape; k++) {
			shape[k]->make_local(P);
			shape[k]->parametrization(P, m);
			shape[k]->make_common(P);
		}
		cpy(p_cpy);
	}
}
template <typename T>
void CShapeMixer<T>::parametrization(int k, double * P, int m)
{
	if (P && shape && k < N_shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		shape[k]->make_local(P);
		shape[k]->parametrization(P, m);
		shape[k]->make_common(P);
		cpy(k, p_cpy);
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CShapeMixer<T>::parametrization_grad(double * P, int m)
{
	if (P && shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		for (int k = 0; k < N_shape; k++) {
			shape[k]->make_local(P);
			shape[k]->parametrization_grad(P, m);
			shape[k]->make_common(P);
		}
		cpy(p_cpy);
	}
}
template <typename T>
void CShapeMixer<T>::parametrization_grad(int k, double * P, int m)
{
	if (P && shape && k < N_shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		shape[k]->make_local(P);
		shape[k]->parametrization_grad(P, m);
		shape[k]->make_common(P);
		cpy(k, p_cpy);
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CShapeMixer<T>::parametrization_hess(double * P, int m)
{
	if (P && shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		for (int k = 0; k < N_shape; k++) {
			shape[k]->make_local(P);
			shape[k]->parametrization_hess(P, m);
			shape[k]->make_common(P);
		}
		cpy(p_cpy);
	}
}
template <typename T>
void CShapeMixer<T>::parametrization_hess(int k, double * P, int m)
{
	if (P && shape && k < N_shape) {
		memcpy(cs, P+3, 3*sizeof(double));
		shape[k]->make_local(P);
		shape[k]->parametrization_hess(P, m);
		shape[k]->make_common(P);
		cpy(k, p_cpy);
	}
}

///////////////////////////////////////////
//...differentiation of shapes composition;
template <typename T>
void CShapeMixer<T>::deriv_X(T * deriv, double f)
{
	if (shape && deriv)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->deriv_X(deriv, f);
			deriv += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::deriv_X(int k, T * deriv, double f)
{
   if (shape && k < N_shape) {
		if (deriv) 
				for (int i = 0; i < k; i++)
					deriv += shape[i]->NN;
		shape[k]->deriv_X(deriv, f);
	}
}
template <typename T>
void CShapeMixer<T>::deriv_Y(T * deriv, double f)
{
	if (shape && deriv)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->deriv_Y(deriv, f);
			deriv += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::deriv_Y(int k, T * deriv, double f)
{
   if (shape && k < N_shape) {
		if (deriv) 
			for (int i = 0; i < k; i++)
				deriv += shape[i]->NN;
		shape[k]->deriv_Y(deriv, f);
	}
}
template <typename T>
void CShapeMixer<T>::deriv_Z(T * deriv, double f)
{
	if (shape && deriv)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->deriv_Z(deriv, f);
			deriv += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::deriv_Z(int k, T * deriv, double f)
{
   if (shape && k < N_shape) {
		if (deriv) 
			for (int i = 0; i < k; i++)
				deriv += shape[i]->NN;
		shape[k]->deriv_Z(deriv, f);
	}
}
template <typename T>
void CShapeMixer<T>::deriv_t(T * deriv, double f)
{
	if (shape && deriv)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->deriv_t(deriv, f);
			deriv += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::deriv_t(int k, T * deriv, double f)
{
   if (shape && k < N_shape) {
		if (deriv) 
			for (int i = 0; i < k; i++)
				deriv += shape[i]->NN;
		shape[k]->deriv_t(deriv, f);
	}
}

///////////////////////
//...normal derivation;
template <typename T>
void CShapeMixer<T>::deriv_N()
{
	if (deriv) memset(deriv, 0, NN_local*sizeof(T));
	deriv_X(deriv, cs[0]);
	deriv_Y(deriv, cs[1]);
	deriv_Z(deriv, cs[2]);
}
template <typename T>
void CShapeMixer<T>::deriv_N(int k)
{
	admittance(k, deriv, NULL, 0., 0.);	
	deriv_X(k, deriv, cs[0]);
	deriv_Y(k, deriv, cs[1]);
	deriv_Z(k, deriv, cs[2]);
}

///////////////////////////////////////////
//...integration of multipoles composition;
template <typename T>
void CShapeMixer<T>::primitive(T * deriv)
{
	if (shape && deriv)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->primitive(deriv);
			deriv += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::primitive(int k, T * deriv)
{
   if (shape && k < N_shape) {
		if (deriv) 
			for (int i = 0; i < k; i++)
				deriv += shape[i]->NN;
		shape[k]->primitive(deriv);
	}
}

//////////////////////////////////////////////////
//...calculation potentials on collocation vector;
template <typename T>
T CShapeMixer<T>::potential(T * p_cpy, int m)
{
	T f(0);
	if (A && p_cpy && m < n)
	for (int k = 0; k < NN; k++) f += A[m][k]*p_cpy[k];
	return f;
}
template <typename T> //...потенциал k-ой компоненты;
T CShapeMixer<T>::potential(int k, T * p_cpy, int m)
{
	T f(0);
	if (A && p_cpy && m < n && shape) {
		int NN_ini = 0, NN_fin = 0, i, j, l;
	
		for (j = i = 0; i < N_shape; j++) 
		if ((j += shape[i]->dim()) < k) NN_ini += shape[i]->NN*shape[i]->dim();	
		else break;

		if (i < N_shape) {
			NN_ini += shape[i]->NN*(k-j+shape[i]->dim());
			NN_fin += shape[i]->NN+NN_ini;
		}
		for (l = NN_ini; l < NN_fin; l++) f += A[m][l]*p_cpy[l];
	}
	return f;
}

/////////////////////////////////////////////
//...auxilliary pointer on colocation vector;
template <typename T>
T * CShapeMixer<T>::FULL(T * p_ext, int k, int l, int shift) 
{
	int i;
   if (p_ext && shape && k < N_shape) {
		if (shift < 0 && k+shift >= 0) for (i = -1; i >= shift; i--) 
			p_ext -= shape[k+i]->NN*shape[k+i]->dim();
		else {
			for (i = 0; i < k; i++)
				p_ext += shape[i]->NN*(shape[i]->dim()-shift);
				p_ext += shape[k]->NN*l;	
		}
	}
	return p_ext;
};

///////////////////////////////////////
//...auxilliary copy colocation vector;
template <typename T>
void CShapeMixer<T>::cpy(T * p_ext) 
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy(p_ext);
			p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy(int k, T * p_ext) 
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++)
			p_ext += shape[i]->NN;
		shape[k]->cpy(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy(int k, T * pp, T * p_ext) 
{
   if (shape && k < N_shape && pp && p_ext) {
		for (int i = 0; i < k; i++) {
			pp		+= shape[i]->NN;
			p_ext += shape[i]->NN;
		}
		shape[k]->cpy(pp, p_ext);
	}
}

//////////////////////////////////////////////////////
//...auxilliary copy derivatives of colocation vector;
template <typename T>
void CShapeMixer<T>::cpy_x(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_x(p_ext);
			if (p_ext) 
				p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_x(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_x(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_y(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_y(p_ext);
			if (p_ext) 
				p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_y(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_y(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_z(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_z(p_ext);
			if (p_ext) 
				p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_z(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_z(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_xx(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_xx(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_xx(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_xx(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_xy(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_xy(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_xy(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_xy(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_yy(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_yy(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_yy(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_yy(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_xz(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_xz(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_xz(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_xz(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_yz(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_yz(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_yz(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_yz(p_ext);
	}
}
template <typename T>
void CShapeMixer<T>::cpy_zz(T * p_ext) 
{
   if (shape)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->cpy_zz(p_ext);
			if (p_ext) 
				 p_ext += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::cpy_zz(int k, T * p_ext) 
{
   if (shape && k < N_shape) {
		if (p_ext) 
			for (int i = 0; i < k; i++)
				p_ext += shape[i]->NN;
		shape[k]->cpy_zz(p_ext);
	}
}

////////////////////////////////////////////////
//...linear combinations with colocation vector;
template <typename T>
void CShapeMixer<T>::adm(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_x(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_x(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_x(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_x(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_y(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_y(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_y(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_y(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_z(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_z(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_z(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_z(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_xx(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_xx(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_xx(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_xx(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_xy(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_xy(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_xy(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_xy(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_yy(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_yy(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_yy(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_yy(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_xz(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_xz(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_xz(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_xz(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_yz(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_yz(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_yz(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_yz(p_ext, adm);
	}
}
template <typename T>
void CShapeMixer<T>::adm_zz(T * p_ext, T adm)
{
   if (shape && p_ext)
		for (int k = 0; k < N_shape; k++) {
			shape[k]->adm_zz(p_ext, adm);
			p_ext  += shape[k]->NN;
		}
}
template <typename T>
void CShapeMixer<T>::adm_zz(int k, T * p_ext, T adm)
{
   if (shape && k < N_shape && p_ext) {
		for (int i = 0; i < k; i++) 
			p_ext += shape[i]->NN;
		shape[k]->adm_zz(p_ext, adm);
	}
}

/////////////////////////////
//...attmidance combinations;
template <typename T>
void CShapeMixer<T>::admittance (T adm)
{
   if (deriv && p_cpy)
	for (int k = 0; k < NN_local; k++) deriv[k] += adm*p_cpy[k];
}
template <typename T>
void CShapeMixer<T>::admittance (int k, T adm)
{
   if (shape && k < N_shape && deriv && p_cpy) {
		for (int i = 0; i < N_shape; i++) {
			 deriv += shape[i]->NN;
			 p_cpy += shape[i]->NN;
		}
		shape[k]->admittance(deriv, p_cpy, 1., adm);
	}
}
template <typename T>
void CShapeMixer<T>::admittance(T * dd, T * pp, T adm_dd, T adm_pp)
{
	int k;
	if (! dd) return;
	for (k = 0;       k < NN_local; k++) dd[k] *= adm_dd;
	for (k = 0; pp && k < NN_local; k++) dd[k] += adm_pp*pp[k];
}
template <typename T>
void CShapeMixer<T>::admittance(int k, T * dd, T * pp, T adm_dd, T adm_pp)
{
   if (shape && k < N_shape && dd) {
		for (int i = 0; i < k; i++) {
			dd += shape[i]->NN;
			if (pp) 
				 pp += shape[i]->NN;
		}
		shape[k]->admittance(dd, pp, adm_dd, adm_pp);
	}
}

/////////////////////////////////////////////////
//...тестовая печать потенциалов для мультиполей;
template <typename T>
void CShapeMixer<T>::TestShape(FILE * TST, int i, int j, double eps, int i_dop)
{
  int k, m;
  if (j >= 0) fprintf(TST, "\nPotentials for %d block:",   j);
  else        fprintf(TST, "\nPotentials for unknown block:");

  if (i >= n)
  for (k = 0; k < NN; k++) {
		 fprintf(TST, "\nA[%2d] =", k);
		 for (m = 0; m < n;  m++) 
			if (typeid(T) == typeid(complex))
				fprintf(TST, "(%12.5lg,%12.5lg)", filtr_Re(A[m][k], eps), filtr_Im(A[m][k], eps)); else
				fprintf(TST, "%12.5lg ", filtr_Re(A[m][k], eps));
  }    
  else if (i >= 0 && i_dop < i)
  for (    k  = 0; k < NN; k++)
		if (typeid(T) == typeid(complex))
			fprintf(TST, "\nA[%d][%d] = \t(%12.5lg,%12.5lg)", i, k, filtr_Re(A[i][k], eps), filtr_Im(A[i][k], eps)); else
			fprintf(TST, "\nA[%d][%d] = \t%12.5lg ", i, k, filtr_Re(A[i][k], eps));
  else if (i >= 0 && i_dop >= i && i_dop < n)
  for (k = 0; k < NN; k++) {
       fprintf(TST, "\nA[%2d][%d-%d] =", k, i, i_dop);
       for (m = i; m <= i_dop;  m++) 
		 if (typeid(T) == typeid(complex))
			 fprintf(TST, "(%12.5lg,%12.5lg)", filtr_Re(A[m][k], eps), filtr_Im(A[m][k], eps)); else
			 fprintf(TST, "%12.5lg ", filtr_Re(A[m][k], eps));
  }    
  fprintf(TST, "\n");
}

////////////////////////////////
//...строковые варианты функции;
template <typename T>
void CShapeMixer<T>::TestShape(char * ch_POTEN, int i, int j, double eps, int i_dop)
{
	char buff[1000], str[40]; ::strcpy (buff, ch_POTEN);
	if (j >= 0) {
		sprintf(str, "%d", j); strcat(buff, str);
	}
	FILE *	 TST = fopen(buff, "w");
	TestShape(TST, i, j, eps, i_dop);
	fclose	(TST);
}
template <typename T>
void CShapeMixer<T>::TestShape(const char * ch_POTEN, int i, int j, double eps, int i_dop)
{
	char buff[1000], str[40]; ::strcpy (buff, ch_POTEN);
	if (j >= 0) {
		sprintf(str, "%d", j); strcat(buff, str);
	}
	FILE *	 TST = fopen(buff, "w");
	TestShape(TST, i, j, eps, i_dop);
	fclose	(TST);
}

////////////////////////////////////////////////////////
//...тестовая визуализация мультиполей в формате Surfer;
template <typename T>
void CShapeMixer<T>::GetSurferFormat(FILE * SURF, FILE * SURF1, double * par, int id_axis)
{
	size_t res;
	int NX = 100, NY = 100;
	if (SURF && SURF1 && par) {
		short int i0 = (short int)(NX+1), 
					 j0 = (short int)(NY+1);
		double pp[7] = {0., 0., 0., 0., 1., 0., 0.}, dd,
				  min1F = 0., max1F = 1., 
				  min2F = 0., max2F = 1., f = 1./NX, g = 1./NY;
		T *	  out_F = (T *)new_struct(3*sizeof(T));

		res = fwrite("DSBB", sizeof(char)*4,  1, SURF);
		res = fwrite(& i0,   sizeof(short int), 1, SURF); res = fwrite(& j0, sizeof(short int), 1, SURF);
		res = fwrite(par,    sizeof(double),  1, SURF); res = fwrite(par+1,  sizeof(double),   1, SURF);
		res = fwrite(par+2,  sizeof(double),  1, SURF); res = fwrite(par+3,  sizeof(double),   1, SURF);
		res = fwrite(& min1F, sizeof(double), 1, SURF); res = fwrite(& max1F, sizeof(double),  1, SURF);

      res = fwrite("DSBB", sizeof(char)*4,  1, SURF1);
      res = fwrite(& i0,   sizeof(short int), 1, SURF1); res = fwrite(& j0, sizeof(short int), 1, SURF1);
      res = fwrite(par,    sizeof(double),  1, SURF1); res = fwrite(par+1,  sizeof(double),   1, SURF1);
      res = fwrite(par+2,  sizeof(double),  1, SURF1); res = fwrite(par+3,  sizeof(double),   1, SURF1);
      res = fwrite(& min2F, sizeof(double), 1, SURF1); res = fwrite(& max2F, sizeof(double),  1, SURF1);

      min1F = MAX_HIT; 
		max1F = MIN_HIT;
      for (int j = 0; j <= NY; j++)
      for (int i = 0; i <= NX; i++) {
			if (id_axis == AXIS_X) {
				pp[1] = par[0]+(par[1]-par[0])*f*i;
				pp[2] = par[2]+(par[3]-par[2])*g*j;
				pp[0] = par[4];
			}
			else
			if (id_axis == AXIS_Y) {
				pp[2] = par[0]+(par[1]-par[0])*f*i;
				pp[0] = par[2]+(par[3]-par[2])*g*j;
				pp[1] = par[4];
			}
			else
			if (id_axis == AXIS_Z) {
				pp[0] = par[0]+(par[1]-par[0])*f*i;
				pp[1] = par[2]+(par[3]-par[2])*g*j;
				pp[2] = par[4];
			}
			else
			if (id_axis == AXIS_SPH) {
				pp[2] = cos(pp[1] = par[2]+(par[3]-par[2])*g*j)*par[4];
				pp[0] = cos(par[0]+(par[1]-par[0])*f*i)*sin(pp[1])*par[4];
				pp[1] = sin(par[0]+(par[1]-par[0])*f*i)*sin(pp[1])*par[4];
			}
			int   lll = sqr(pp[0])+sqr(pp[1])+sqr(pp[2]) < sqr(get_param(0)); if (id_axis == AXIS_SPH) lll = 0;
			if (! lll){
				parametrization_hess(pp);
				out_F[0] = potential(p_cpy, 0);
			}

/////////////////////////////////////////////////////////////////////////
//deriv_N();//...тестирование сшивки производных;
//out_F[1] = -potential(deriv, 0); if (! inverse()) out_F[1] *= get_param(3);
/////////////////////////////////////////////////////////////////////////
//			parametrization_hess(pp);
//			out_F[1] = potential(p_cpy, 0);
//			change();
/////////////////////////////////////////////////////////////////////////
			float  ff = NOT_HIT;
			if (lll) {
				dd  = real(out_F[0]);
				res = fwrite(& (ff = static_cast<float>(dd)), sizeof(float), 1, SURF);
				if (min1F > dd) min1F = dd;
				if (max1F < dd) max1F = dd;

				if (typeid(T) == typeid(complex))
					dd = imag(out_F[0]); else
					dd = real(out_F[1]);
				res = fwrite(& (ff = static_cast<float>(dd)), sizeof(float), 1, SURF1);
				if (min2F > dd) min2F = dd;
				if (max2F < dd) max2F = dd;
			}
      }
      delete_struct(out_F);

//////////////////////////////////////////////////////////////
//...перезапись максимального и минимального значения функции;
      res = fseek(SURF, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
      res = fwrite(& min1F, sizeof(double), 1, SURF);
      res = fwrite(& max1F, sizeof(double), 1, SURF);
      res = fseek(SURF, 0L, SEEK_END);

      res = fseek(SURF1, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
      res = fwrite(& min2F, sizeof(double), 1, SURF1);
      res = fwrite(& max2F, sizeof(double), 1, SURF1);
      res = fseek(SURF1, 0L, SEEK_END);
	}
}

////////////////////////////////
//...строковые варианты функции;
template <typename T>
void CShapeMixer<T>::GetSurferFormat(char * SURF_FILE, double * par, int id_axis)
{
	char buff[1000]; ::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	FILE * SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	FILE * SURF1 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, par, id_axis);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
}
template <typename T>
void CShapeMixer<T>::GetSurferFormat(const char * SURF_FILE, double * par, int id_axis)
{
	char buff[1000]; ::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	FILE * SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	FILE * SURF1 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, par, id_axis);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
}
#endif
