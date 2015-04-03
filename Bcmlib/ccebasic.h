/*============================================*/
/*                  CCEBASIC                  */
/*============================================*/
#ifndef ___CCEBASIC___
#define ___CCEBASIC___

#include "ccells.h"

//////////////////////////////////////
//...class of basic geometrical cells;
class CCeBasic : public CCells {
		virtual int type() { return BASIC_CELLS;}
public:
//...plane facet compositions;
      void get_facet_directly(double * P, int N, int id_plane = NULL_STATE, int id_fast = NULL_STATE, int id_dbl_point = NULL_STATE, double eps = EE_ker);
      void get_tria_facet(double * P1, double * P2, double * P3, int id_fast = OK_STATE);
      void get_quad_facet(double * P1, double * P2, double * P3, double * P4, int id_fast = OK_STATE);
 //...parameters of plane facet compositions;
      int  SetFacetParam (int k, double square, double * pp = NULL, int id_property = DEFAULT_BND);
      void SetFacetXParam(int k, double X0, double * pp, int id_property) {
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[4])-1.) < EE_ker && fabs(ce[k]->mp[1]-X0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
      void SetFacetYParam(int k, double Y0, double * pp, int id_property) {
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[5])-1.) < EE_ker && fabs(ce[k]->mp[2]-Y0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
      void SetFacetZParam(int k, double Z0, double * pp, int id_property) {
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[6])-1.) < EE_ker && fabs(ce[k]->mp[3]-Z0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
 //...parameters of node compositions;
		int  SetNodeParam	  (double * pp = NULL, int id_property = DEFAULT_BND);
		void SetNodeXParam  (double X0, double * pp, int id_property);
		void SetNodeYParam  (double Y0, double * pp, int id_property);
		void SetNodeZParam  (double Z0, double * pp, int id_property);
		void SetNodeXYParam (double X0, double Y0, double * pp, int id_property);
		void SetNodeXZParam (double X0, double Z0, double * pp, int id_property);
		void SetNodeYZParam (double Y0, double Z0, double * pp, int id_property);
		void SetNodeXYZParam(double X0, double Y0, double Z0, double * pp, int id_property);
//...polygonal compositions and boundary conditions of polygonal compositions;
      void get_nd_bar_directly(CGrid * nd, int k);
      int  set_nd_bar_condit(int k, CGrid * nd, int k1, int k2, int id_property, double * pp, int * id_pp, double & S, int id_facet = -1);
      void set_nd_bar_condit(int k, int id_property, double * pp, int * id_pp, double S);
  //...library of basic elements;
      void get_point				(double * P);
      void get_point				(double X, double Y, double Z);
      void get_circle			(double R);
      void get_disk				(double R);
      void get_ring				(double rad, double R);
      void get_sphere			(double R);
		void get_spheroid			(double A, double B);
		void get_ellipsoid		(double A, double B, double C);
      void get_ball				(double R);
      void get_cylinder			(double R);
      void get_roll				(double R);
      void get_space				(int N_dim);
      void get_line				(double * P1, double * P2, int id_cell = OK_STATE, int id_fast = NULL_STATE);
      void get_arc				(double R, double f0, double f1, int id_cell = OK_STATE);
      void get_arc				(double R, complex z1, complex z2, int topo);
      void get_ellipt_arc		(double A, double B, double f0, double f1, int id_cell = OK_STATE);
      void get_sheet				(double A, double B);
		void get_sheet_intrusion(double A, double B, double rad);
      void get_box				(double A, double B, double C);
		void get_sph_intrusion	(double R, double L);
//...plane and space composition of cells;
      void get_polygon_directly    (double * X, double * Y, int N, int id_dbl = 0, double eps = EE_ker);
      void get_arc_polygon_directly(double * X, double * Y, double * R, int * topo, int N, int id_dbl = 0, double eps = EE_ker);
      void get_line_strip_directly (double * X, double * Y, double * Z, int * geom, int * mask = NULL, int shift = 2, int id_fast = OK_STATE, int id_dbl = NULL_STATE, double eps = EE_ker);
      void get_beam                (CCells * f0, double length);
      void get_cyl_beam            (CCells * f0, double fi);
//...geometrical segments;
      void get_ring_segment (              double f0, double f1, double t0, double t1);
      void get_cyl_segment  (double R,     double f0, double f1, double t0, double t1);
      void get_sph_segment  (double R,     double f0, double f1, double t0, double t1);
      void get_cone_segment (double theta, double f0, double f1, double t0, double t1);
      void get_torus_segment(double r0,    double r1, double f0, double f1, double t0, double t1);
      void get_ugolok_cell  (double A1, double A2, double B1, double B2, double radc,
                             double rad, double beta, double x0 = 0., double y0 = 0., double f0 = 0.);
//...special compositions;
		void get_circle_profile	  (double r);
      void get_ugolok           (double beta, double radc, double rad, double A1, double B1, double A2, double B2);
      void get_blend_cyl_segment(double beta, double R, double L, int * mm = NULL);
};
////////////////////////////////////////
//...complete deleting of abstract cell;
inline void delete_cells(CCeBasic *& ce) {
	ce->zero_cells(); delete ce; ce = NULL;
}
#endif
