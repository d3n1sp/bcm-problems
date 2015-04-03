/*=========================================*/
/*                 CSHAPE                  */
/*=========================================*/
#ifndef ___CSHAPE___
#define ___CSHAPE___

#include "cbase.h"
#include "kernel.h"
#include "dd_real.h"
#include "qd_real.h"

/////////////////////////////////////////////////
//...description of existing types of multipoles;
enum Num_Shape {
						NULL_SHAPE = 0,
///////////////////////////
//...гармонические функции;
				 MP2D_POLY_SHAPE,
				 MP3D_POLY_SHAPE,
				 MP3D_ZOOM_SHAPE,
				 MP3D_ELLI_SHAPE,
				 BEAM_POLY_SHAPE,
			  MP2D_CORNER_SHAPE,
			  MP2D_CLAYER_SHAPE,
			  MP2D_CIRCLE_SHAPE,
			  MP3D_CLAYER_SHAPE,
			  MP3D_SPHERE_SHAPE,
			MP3D_SPHEROID_SHAPE,
	 MP3D_SPHEROID_FULL_SHAPE,
//////////////////////
//...тепловые функции;
				 HT1D_POLY_SHAPE,
				 HT2D_POLY_SHAPE,
				 HT3D_POLY_SHAPE,
			  HT2D_STREEP_SHAPE,
			  HT2D_CIRCLE_SHAPE,
			  HT3D_SPHERE_SHAPE,
//////////////////////////////
//...экспоненциальные функции;
				 SK2D_BEAMZSHAPE,
				 SK3D_BEAM_SHAPE,
				 SK2D_ELLI_SHAPE,
				 SK2D_POLY_SHAPE,
				 SK3D_POLY_SHAPE,
				 SK3D_ZOOM_SHAPE,
				 SK3D_EXPP_SHAPE,
//////////////////////
//...волновые функции;
				 WV2D_POLY_SHAPE,
				 WV3D_POLY_SHAPE,
///////////////////////////
//...осцилл€цонные функции;
				 AU3D_ELLI_SHAPE,
				 AU3D_POLY_SHAPE,
				 AU3D_ZOOM_SHAPE,
				 AU3D_WAVE_SHAPE,
				 AU3D_BEAM_SHAPE,
				 AU2D_BEAMZSHAPE,
				 AU3D_BEAMZSHAPE,
				  AU_YJUMP_SHAPE,
			  AU3D_CORNER_SHAPE,
						NUMS_SHAPE
};

/////////////////////////////////////////////////
//...базовый класс дл€ специальных функций формы;
template <typename T>
class CShape : public CBase {
public:
		int	 NN;				//...number of local freedom;
		int    N, M, N1, N2;	//...degree of multipoles;
		double R0, R0_inv;	//...normalizing radius;
protected:
		int id_inverse, id_dim;								//...internal or external shapes and dimension;
		int NN_dop, NN_grad, NN_hess;						//...internal degrees of freedom;
		T * p, * px, * py,  * pz, * pxx, * pxy,		//...parametrization of multipoles;                
						 * pyy, * pxz, * pyz, * pzz;
		double CX, CY, CZ, SX, SY, SZ, P0[3], cs[6];	//...local coordinate system;
public:
		virtual int type() {return NULL_SHAPE;}
		virtual int freedom(int m) {return 0;}
//...auxilliary functions;
		virtual void release() { 
			delete_struct(p);   delete_struct(px);  
			delete_struct(py);  delete_struct(pz);
			delete_struct(pxx); delete_struct(pxy); delete_struct(pyy);  
			delete_struct(pxz); delete_struct(pyz); delete_struct(pzz);
		}
		virtual void cpy(T * p_ext) {
			if (p_ext && p)
				memcpy(p_ext, p, NN*sizeof(T));
		}
		void cpy_x(T * p_ext = NULL) {
			if (! p || ! px) return;
			if (! p_ext) swap(p, px); else 
				memcpy(p_ext, px, NN*sizeof(T));
		}
		void cpy_y(T * p_ext = NULL) {
			if (! p || ! py) return;
			if (! p_ext) swap(p, py);	else 
				memcpy(p_ext, py, NN*sizeof(T));
		}
		void cpy_z(T * p_ext = NULL) {
			if (! p || ! pz) return;
			if (! p_ext) swap(p, pz); else 
				memcpy(p_ext, pz, NN*sizeof(T));
		}
		void cpy_xx(T * p_ext = NULL) {
			if (! p || ! pxx) return;
			if (! p_ext) swap(p, pxx); else 
				memcpy(p_ext, pxx, NN*sizeof(T));
		}
		void cpy_xy(T * p_ext = NULL) {
			if (! p || ! pxy) return;
			if (! p_ext) swap(p, pxy); else 
				memcpy(p_ext, pxy, NN*sizeof(T));
		}
		void cpy_yy(T * p_ext = NULL) {
			if (! p || ! pyy) return;
			if (! p_ext) swap(p, pyy); else 
				memcpy(p_ext, pyy, NN*sizeof(T));
		}
		void cpy_xz(T * p_ext = NULL) {
			if (! p || ! pxz) return;
			if (! p_ext) swap(p, pxz); else 
				memcpy(p_ext, pxz, NN*sizeof(T));
		};
		void cpy_yz(T * p_ext = NULL) {
			if (! p || ! pyz) return;
			if (! p_ext) swap(p, pyz); else 
				memcpy(p_ext, pyz, NN*sizeof(T));
		}
		void cpy_zz(T * p_ext = NULL) {
			if (! p || ! pzz) return;
			if (! p_ext) swap(p, pzz); else 
				memcpy(p_ext, pzz, NN*sizeof(T));
		}
		void cpy(T * pp, T * p_ext) {
			if (p_ext && pp) 
				memcpy(p_ext, pp, NN*sizeof(T));
		}
		void adm	  (T * p_ext, T adm);
		void adm_x (T * p_ext, T adm);
		void adm_y (T * p_ext, T adm);
		void adm_z (T * p_ext, T adm);
		void adm_xx(T * p_ext, T adm);
		void adm_xy(T * p_ext, T adm);
		void adm_yy(T * p_ext, T adm);
		void adm_xz(T * p_ext, T adm);
		void adm_yz(T * p_ext, T adm);
		void adm_zz(T * p_ext, T adm);
		void admittance(T * dd, T * pp, T adm_dd, T adm_pp);
//...dimension and inverse;
		int  dim			 (){ return id_dim;}
		int  inverse	 (){ return id_inverse;}
		void change		 (){ id_inverse = id_inverse ? NULL_STATE : OK_STATE;}
		void set_inverse(int m_inv = 0){ id_inverse = m_inv;}
//...setting of coordinate system;
		void set_local(double * P = NULL) {
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
		void set_local_P0(double * P) { if (P) memcpy(P0, P, 3*sizeof(double));}
		void make_local  (double * P) {
			point_shift<double>(P, P0, OK_STATE);
			point_iso<double>  (P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
		}
		void make_common (double * P) { point_iso<double>(P, P0, CZ, SZ, CY, SY, CX, SX);}
		void norm_local  (double * P) { point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);}
		void norm_common (double * P) { point_iso<double>(P, NULL, CZ,  SZ, CY,  SY, CX,  SX);}
//...constructor and destructor;
		CShape() {
			N = M  = N1 = N2 = NN = 0;
			NN_dop = NN_grad = NN_hess = 0;
			R0 = R0_inv  = 1.;
			p  = px = py = pz = pxx = pxy = pyy = pxz = pyz = pzz = NULL;
			id_inverse = id_dim = 0; set_local();
		}
		virtual ~CShape(void) { release();}
//...initialization multipoles;
		virtual void init1(int N, int dim) { this->N = N; this->id_dim = dim; NN = freedom(N); release();}
		virtual void init2(int N, int M, int dim) { init1(M, dim);}
		virtual void init3(int N, int N1, int N2, int dim) { init1(N, dim);}
		virtual void set_shape(double R0, double kk = 0., double p1 = 0., double p2 = 0., double p3 = 0., double p4 = 0.){
			int   k  = -1;
			this->R0 = R0;
			R0_inv   = R0 > EE ? 1./R0 : 1.;
			if (++k < this->size_of_param()) this->param[k] = (Param)(kk*R0); else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)(this->param[0] > EE ? 1./this->param[0] : 1.); else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)(this->param[0]*this->param[0]); else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)kk; else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)p1; else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)p2; else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)p3; else return;
			if (++k < this->size_of_param()) this->param[k] = (Param)p4; else return;
		}
//...calculation multipoles;
		virtual void parametrization     (double * P = NULL, int m_dop = 0){}
		virtual void parametrization_grad(double * P = NULL, int m_dop = 0);
		virtual void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation multipoles;
		virtual void deriv_X(T * deriv, double f = 1.){}
		virtual void deriv_Y(T * deriv, double f = 1.){}
		virtual void deriv_Z(T * deriv, double f = 1.){}
		virtual void deriv_t(T * deriv, double f = 1.){}
//...integration multipoles;
		virtual void primitive(T * deriv){}
};
///////////////////////////////////////////////////////////
//...construction of all existing types of shape functions;
CShape<double>  * CreateShapeS(Num_Shape id_SHAPE, int id_dop = 0);
CShape<double>  * CreateShapeD(Num_Shape id_SHAPE, int id_dop = 0);
CShape<complex> * CreateShapeC(Num_Shape id_SHAPE, int id_dop = 0);
//CShape<dd_real> * CreateShapeR(Num_Shape id_SHAPE, int id_dop = 0);
//CShape<qd_real> * CreateShapeQ(Num_Shape id_SHAPE, int id_dop = 0);

////////////////////////////////////////////////////
//           TEMPLATE VARIANT OF SHAPES           //
////////////////////////////////////////////////////
//////////////////////////////////////////////
//...combination of functions with attmidance;
template <typename T>
void CShape<T>::adm(T * p_ext, T adm)
{
	if (p_ext && p)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*p[k];
}
template <typename T>
void CShape<T>::adm_x(T * p_ext, T adm)
{
	if (p_ext && px)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*px[k];
}
template <typename T>
void CShape<T>::adm_y(T * p_ext, T adm)
{
	if (p_ext && py)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*py[k];
}
template <typename T>
void CShape<T>::adm_z(T * p_ext, T adm)
{
	if (p_ext && pz)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pz[k];
}
template <typename T>
void CShape<T>::adm_xx(T * p_ext, T adm)
{
	if (p_ext && pxx)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pxx[k];
}
template <typename T>
void CShape<T>::adm_xy(T * p_ext, T adm)
{
	if (p_ext && pxy)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pxy[k];
}
template <typename T>
void CShape<T>::adm_yy(T * p_ext, T adm)
{
	if (p_ext && pyy)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pyy[k];
}
template <typename T>
void CShape<T>::adm_xz(T * p_ext, T adm)
{
	if (p_ext && pxz)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pxz[k];
}
template <typename T>
void CShape<T>::adm_yz(T * p_ext, T adm)
{
	if (p_ext && pyz)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pyz[k];
}
template <typename T>
void CShape<T>::adm_zz(T * p_ext, T adm)
{
	if (p_ext && pzz)
		for (int k = 0; k < NN; k++) p_ext[k] += adm*pzz[k];
}
template <typename T>
void CShape<T>::admittance(T * dd, T * pp, T adm_dd, T adm_pp)
{
	int k;
	if (! dd) return;
	for (k = 0;       k < NN; k++) dd[k] *= adm_dd;
	for (k = 0; pp && k < NN; k++) dd[k] += adm_pp*pp[k];
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CShape<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! px && ! py && ! pz) {
		px = (T *)new_struct((NN_grad = freedom(N+m_dop))*sizeof(T));
		py = (T *)new_struct( NN_grad*sizeof(T));
		pz = (T *)new_struct( NN_grad*sizeof(T));
	}
	if (P && px && py && pz) {
		memset (px, 0, NN_grad*sizeof(T));
		memset (py, 0, NN_grad*sizeof(T));
		memset (pz, 0, NN_grad*sizeof(T));
		N += m_dop;
		deriv_X(px);
		deriv_Y(py);
		deriv_Z(pz);
		N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CShape<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! pxx && ! pxy && ! pyy && ! pxz && ! pyz && ! pzz) {
		pxx = (T *)new_struct((NN_hess = freedom(N+m_dop))*sizeof(T));
		pxy = (T *)new_struct( NN_hess*sizeof(T));
		pyy = (T *)new_struct( NN_hess*sizeof(T));
		pxz = (T *)new_struct( NN_hess*sizeof(T));
		pyz = (T *)new_struct( NN_hess*sizeof(T));
		pzz = (T *)new_struct( NN_hess*sizeof(T));
	}
	if (P && pxx && pxy && pyy && pxz && pyz && pzz) {
		memset (pxx, 0, NN_hess*sizeof(T));
		memset (pxy, 0, NN_hess*sizeof(T));
		memset (pyy, 0, NN_hess*sizeof(T));
		memset (pxz, 0, NN_hess*sizeof(T));
		memset (pyz, 0, NN_hess*sizeof(T));
		memset (pzz, 0, NN_hess*sizeof(T));
		N += m_dop;
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
		N -= m_dop;
	}
}
#endif
