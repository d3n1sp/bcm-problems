/*=========================================*/
/*                CGRID_QG                 */
/*=========================================*/
#ifndef ___CGRID_QG___
#define ___CGRID_QG___

#include "cgrid.h"

//////////////////////////////////////////////////////////////////////
//...class of gaussian quadratures for all existing types of elements;
class CGrid_QG : public CGrid {
		static const double cc[][2];
public:
		virtual int type() { return GRID_QG_NODES;}
public:
//...normalization and weihgts/nodes of quadratures;
		double QG_normal(double & nX, double & nY, double & nZ) {
			double t, F;
			if ((F = sqrt(sqr(nX)+sqr(nY)+sqr(nZ))) > EE) {
				 nX *= (t = 1./F); nY *= t; nZ *= t;
			}
			else {
				 nX = nY = 0.; nZ = 1.;
			}
			return(F);
		}
//...функция, генерирующая квадратуры Гаусса;
		void QG_generate(int N, long double cc[][2]);
protected:
//...gaussian quadratures for polygonal elements;
		void QG2M_line (double * P, int M, const double c2M [][2], int normal);
		void QG2M_tria (double * P, int M, const double c2M [][2], int normal);
		void QG2M_quad (double * P, int M, const double c2M [][2], int normal);
		void QG2M1_line(double * P, int M, const double c2M1[][2], int normal);
		void QG2M1_tria(double * P, int M, const double c2M1[][2], int normal);
		void QG2M1_quad(double * P, int M, const double c2M1[][2], int normal);
//...gaussian quadratures for dimensional elements;
		void QG_sheet(CMap * mp, int N_elem, int N_max);
		void QG_sheet(CMap * mp, int N_elem);
		void QG_ring_segment	(CMap * mp, int N_elem, int N_max);
		void QG_ring_segment	(CMap * mp, int N_elem);
		void QG_cyl_segment	(CMap * mp, int N_elem, int N_max);
		void QG_cyl_segment	(CMap * mp, int N_elem);
		void QG_sph_segment	(CMap * mp, int N_elem, int N_max);
		void QG_sph_segment	(CMap * mp, int N_elem);
		void QG_spr_segment	(CMap * mp, int N_elem, int N_max);
		void QG_spr_segment	(CMap * mp, int N_elem);
		void QG_sht_intrusion(CMap * mp, int N_elem, int N_max);
		void QG_sph_intrusion(CMap * mp, int N_elem, int N_max);
		void QG_sph_intrusion(CMap * mp, int N_elem);
		void QG_cone_segment (CMap * mp, int N_elem, int N_max);
		void QG_cone_segment	(CMap * mp, int N_elem);
		void QG_torus_segment(CMap * mp, int N_elem, int N_max);
		void QG_torus_segment(CMap * mp, int N_elem);
		void QG_ugolok_cell	(CMap * mp, int N_elem);
//...correction weigts and nodes for curvilinear surface;
		void QG_circle  (CMap * mp);
		void QG_ellipse (CMap * mp);
		void QG_tria_circle (CMap * mp, double * Po, int N_ini = 0);
		void QG_tria_ellipse(CMap * mp, double * Po, int N_ini = 0);
		void QG_quad_circle (CMap * mp1, CMap * mp2, double * Po, int N_ini = 0);
		void QG_quad_ellipse(CMap * mp1, CMap * mp2, double * Po, int N_ini = 0);
		void QG_sphere  (CMap * mp);
		void QG_cylinder(CMap * mp);
		void QG_spheroid(CMap * mp);
		void QG_tria_sphere (CMap * mp, double * Po, int N_ini = 0);
		void QG_quad_sphere (CMap * mp1, CMap * mp2, double * Po, int N_ini = 0);
public:
//...gaussian quadratures;
		void facet_QG(double * P, int N_elem, int mode = NULL_STATE, int normal = OK_STATE);
		void segms_QG(CMap * mp, int N_elem, int N_max);
//...special gaussian quadratures;
		void sphere_intrusion_QG(CMap * mp, int N_elem, double L);
//...gaussian quadratures corrections;
		void QG_curve	 (CMap * mp);
		void QG_surface (CMap * mp);
		void QG_tria_curve  (CMap * mp, double * Po, int N_ini = 0);
		void QG_quad_curve  (CMap * mp1, CMap * mp2, double * Po, int N_ini = 0);
		void QG_tria_surface(CMap * mp, double * Po, int N_ini = 0);
		void QG_quad_surface(CMap * mp1, CMap * mp2, double * Po, int N_ini = 0);
		void QG_quad_bi(double * Po, int N_ini = 0);
};
#endif
