/*=========================================*/
/*                CGRID_EL                 */
/*=========================================*/
#ifndef ___CGRID_EL___
#define ___CGRID_EL___

#include "cgrid.h"

////////////////////////////////////////////////////////////
//...class of grid nodes for all existing types of elements;
class CGrid_el : public CGrid {
public:
		virtual int type() { return GRID_EL_NODES;}
public:
//...uniform grid for parametric boxes and cells;
		void grid_box1D		 (int N1, int shift = 0, int id_link = 0, double fy = 0., double fz = 0.);
		void grid_box2D		 (int N1, int N2, int shift = 0, int id_link = 0, double fz = 0.);
		void grid_box3D		 (int N1, int N2, int N3, int shift = 0, int id_link = 0);
		void grid_tria_3D		 (int & N_elem, double * P, int shift);
		void grid_hexa_refine (int N_elem, double * P);
		void grid_tetra_refine(int N_elem, double * P);
//...uniform grid for meshless block method;
		int  grid_sph_2D    (int N_elem, CMap * mp);
		int  grid_sph_3D    (int N_elem, CMap * mp);
		int  grid_cyl_3D    (int N_elem, CMap * mp, int id_flag);
		int  grid_box_3D    (int N_elem, CMap * mp, int dir);
		int  grid_box_box_2D(int N_elem, CMap * mp, CMap * mp_ext);
		void grid_box_box_3D(int N_elem, CMap * mp, CMap * mp_ext, int id_flag, int & dir);
		void grid_sph_sph_2D(int N_elem, CMap * mp, CMap * mp_ext);
		void grid_sph_sph_3D(int N_elem, CMap * mp, CMap * mp_ext);
		int  grid_sph_sph_3D_diameter(int N_elem, CMap * mp, CMap * mp_ext, double theta0);
//...uniform grid of BARs elements;
		complex make_common(complex z0, complex t0, complex z) { return z*t0+z0;}
		complex make_local (complex z0, complex t0, complex z) { return (z-z0)*conj(t0);}
		int grid_quad		 (complex * w, double *& X, double *& Y, double h1, double h2,	int *& geom, int & N, int shift, int Max_N1, int Max_N2);
		int grid_knee		 (complex * w, double *& X, double *& Y, double h1, double h2, int *& geom, int & N, int shift, int Max_N1, int Max_N2);
		void set_geom		 (int *& geom, int &N1, int &N2, int &m1, int &N, int &shift, double &f1, double &f2);
protected:
		int in_circle	(CMap * mp, double X = 0., double Y = 0., double eps = EE);
		int in_ellipse (CMap * mp, double X = 0., double Y = 0., double eps = EE);
		int in_sphere	(CMap * mp, double X = 0., double Y = 0., double Z = 0., double eps = EE);
		int in_spheroid(CMap * mp, double X = 0., double Y = 0., double Z = 0., double eps = EE);
public:
//...geometrical map interior;
		int maps_in(CMap * mp, double * P, double eps = EE) {
		if ( in_circle  (mp, P[0], P[1], eps) ||
			  in_ellipse (mp, P[0], P[1], eps) ||
			  in_sphere  (mp, P[0], P[1], P[2], eps) ||
			  in_spheroid(mp, P[0], P[1], P[2], eps)
			 ) return(1);
		else  return(0);
		}
protected:
		int in_sheet(CMap * mp, double X = 0., double Y = 0., double eps = EE);
		int in_box  (CMap * mp, double X = 0., double Y = 0., double Z = 0., double eps = EE);
public:
//...dimension cells interior;
		int segms_in(CMap * mp, double * P, double eps = EE) {
		if ( in_sheet(mp, P[0], P[1], eps) ||
			  in_box  (mp, P[0], P[1], P[2], eps)/* ||
			  in_ring (mp, P[0], P[1], P[2], eps))*/
			 ) return(1);
		else  return(0);
		}
};
#endif
