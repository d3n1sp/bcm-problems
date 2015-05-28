/*============================================*/
/*                  CCELLS                    */
/*============================================*/
#ifndef ___CCELLS___
#define ___CCELLS___

#include "cgrid.h"

////////////////////////////////////////////////////////
//...description of existing types of geometrical cells;
enum Num_Cells {
        NULL_CELLS = 0,
       //BASIC_CELLS,
        NUMS_CELLS
};

///////////////////////////////////////////////////////////////////////////
//...the main object of geometrical kernel -- basic class of spatial cells;
class CCells {
public:
		CCells ** ce;    //...list of all elements of the cell (including itself);
		Topo   *  graph; //...list of links (topology of boundary);
		CMap   *  mp;    //...geometrical map of the cell;
		CMap   ** pm;    //...geometrical map of the cell in 2D parametric space;
public:
		virtual int type() { return NULL_CELLS;}
//...constructors/destructor;
       CCells(int N_ce = 0) {
              ce    = NULL;
              graph = NULL;
              mp    = NULL;
              pm    = NULL;
              if (N_ce >= 0) cells_new(N_ce+1, 2, 0);
       }
      virtual ~CCells(void) { //...здесь нельз€ уничтожать все €чейки "автоматом" функцией zero_cells()!!!
		}
public:
		inline void add_pm(int k, CMap *& mp) {
			if (graph && 0 <= k && k < graph[1] && pm) { pm[k] = mp; mp = NULL;}
		}
		inline void delete_cells(CCells *& ce, Num_State id_zero = OK_STATE) {
			if (id_zero == OK_STATE) ce->zero_cells(); 
			delete ce; ce = NULL;
		}
		void cells_new (int N_ce = 0, int N_graph = 0, int N_mp = 0);
		void zero_cells();
//...cell printing and reading in ASCII;
		void cells_out(char * ch_CELLS, int id_long = 0);
		void cells_out(FILE * id_CELLS, int id_long = 0);
		int  cells_in (char * id_CELLS, unsigned long & count, unsigned long upper_limit);
		void cells_in (char * ch_CELLS);
//...dimension characteristics of the cell;
		inline int cells_genus() {
			return map_genus(mp);
		}
		inline int cells_dim() {
			return map_dim(mp);
		}
		inline void cells_to(int genus) {
			map_to(mp, map_dim(mp), genus);
		}
		int bar_dim();
		int arcs_number();
		int arcs_number(int id_arc);
		int circ_number(int id_arc);
//...isometric transformations;
		void cells_iso(double * P, double &CZ, double &SZ, double &CY, double &SY, double &CX, double &SX);
		void cells_iso(double * P, double fi = 0., double theta = 0., double fX = 0.);
		void	 bar_iso(double * P, double fi = 0., double theta = 0., double fX = 0.);
		inline void cells_to_origine() {
			if (mp) {
				double P[3] = {-mp[1], -mp[2], -mp[3]};
				cells_iso(P);
				if (ID_MAP(0, NULL_GENUS) != mp[0])
					cells_iso (NULL, -mp[6], -mp[5], -mp[4]);
			}
		}
		inline void cells_local(CMap * ext_mp) {
			if (ext_mp) {
				double P[3] = {-ext_mp[1], -ext_mp[2], -ext_mp[3]};
				cells_iso (P);
				if (ID_MAP(0, NULL_GENUS) != ext_mp[0])
					cells_iso (NULL, -ext_mp[6], -ext_mp[5], -ext_mp[4]);
			}
		}
		inline void cells_common(CMap * ext_mp) {
			if (ext_mp) {
				double P[3] = {ext_mp[1], ext_mp[2], ext_mp[3]};
				if (ID_MAP(0, NULL_GENUS) != ext_mp[0])
				cells_iso (NULL, ext_mp[4], ext_mp[5], ext_mp[6]);
				cells_iso (P);
			}
		}
		inline void bar_local(CMap * ext_mp) {
			if (ext_mp) {
				double P[3] = {-ext_mp[1], -ext_mp[2], -ext_mp[3]};
				bar_iso (P);
				if (ID_MAP(0, NULL_GENUS) != ext_mp[0])
					bar_iso (NULL, -ext_mp[6], -ext_mp[5], -ext_mp[4]);
			}
		}
		inline void bar_common(CMap * ext_mp) {
			if (ext_mp) {
				double P[3] = {ext_mp[1], ext_mp[2], ext_mp[3]};
				if (ID_MAP(0, NULL_GENUS) != ext_mp[0])
				bar_iso (NULL, ext_mp[4], ext_mp[5], ext_mp[6]);
				bar_iso (P);
			}
		}
//...spatial arc of the cell;
		inline double cells_length() {
			double length = 0.;
			if (mp) return(length);
			if (mp[0] == ID_MAP(1, NULL_GENUS)) {
				length = abs_point(ce[get_num(graph, 0)]->mp+1,
										 ce[get_num(graph, 1)]->mp+1);
			}
			else  
			if (mp[0] == ID_MAP(1, SPHERE_GENUS)) {
				length = fabs(mp[7]*2.*mp[8]);
			}
			else  
			if (mp[0] == ID_MAP(1, CYL_GENUS)) {
				length = 0.;
			}
			else  
			if (mp[0] == ID_MAP(1, CONE_GENUS)) {
				length = 0.;
			}
			else  
			if (mp[0] == ID_MAP(1, ELL_CONE_GENUS)) {
				length = 0.;
			}
			return(length);
		}
		inline void edge_correct() {
			if (mp && graph[1] > 0 && mp[0] == ID_MAP(1, NULL_GENUS)) {
				int m0 = get_num(graph, 0), 
					 m1 = get_num(graph, 1);
				mp[1] = (CMap)(ce[m0]->mp[1]+ce[m1]->mp[1])*.5;
				mp[2] = (CMap)(ce[m0]->mp[2]+ce[m1]->mp[2])*.5;
				mp[3] = (CMap)(ce[m0]->mp[3]+ce[m1]->mp[3])*.5;
			}
			else  
			if (mp[0] == ID_MAP(1, SPHERE_GENUS)) {
				int m1 = get_num(graph, 0),
					 m2 = get_num(graph, 1);
				double  CZ = cos(mp[6]), SZ = -sin(mp[6]),
						  CY = cos(mp[5]), SY = -sin(mp[5]),
						  CX = cos(mp[4]), SX = -sin(mp[4]), f1, f2,
						  p1[3] = { ce[m1]->mp[1]-mp[1],
										ce[m1]->mp[2]-mp[2],
										ce[m1]->mp[3]-mp[3] },
						  p2[3] = { ce[m2]->mp[1]-mp[1],
										ce[m2]->mp[2]-mp[2],
										ce[m2]->mp[3]-mp[3] };
				complex z1, z2;

				point_iso<double>(p1, NULL, CZ, SZ, CY, SY, CX, SX); f1 = arg2(z1 = comp(p1[0], p1[1]));
				point_iso<double>(p2, NULL, CZ, SZ, CY, SY, CX, SX); f2 = arg2(z2 = comp(p2[0], p2[1]));
				if (z1 == z2 || f2+EE_ker < f1) f2 += 2.*M_PI;

				mp[6] += .5*(f1+f2);
				mp[8]  = .5*fabs(f2-f1);
			}
			else  
			if (mp[0] == ID_MAP(1, CYL_GENUS)) {
			}
			else  
			if (mp[0] == ID_MAP(1, CONE_GENUS)) {
			}
			else  
			if (mp[0] == ID_MAP(1, ELL_CONE_GENUS)) {
			}
		}
//...cells identification;
		int   topo_id(Topo some_element, Topo exc_element = -1);
		int common_id(CCells * ce);
//...cells ordering;
		inline void bar_mutat(int k, int m) {
			if (! mp && k < graph[0] && m < graph[0]) {
				graph[0]--;
				topo_correct(k, -m-1); 
				topo_correct(m, -k-1); topo_correct(); 
				graph[0]++;
				swap(ce[k], ce[m]);//(CCells * temp_ce = ce[k]; ce[k] = ce[m]; ce[m] = temp_ce;)
			}
		}
		inline void cells_mutat(int k, int m) {
			if (mp && k < graph[0] && m < graph[0]) {
				topo_correct(k, -m-1); 
				topo_correct(m, -k-1); topo_correct(); 
				swap(ce[k], ce[m]);//(CCells * temp_ce = ce[k]; ce[k] = ce[m]; ce[m] = temp_ce;)
			}
		}
		void   bar_ord();
		void cells_ord();
		void cells_ord0(int id_ord = 0);
		void bar_invers();
		void cells_invers1();
		void cells_invers2();
		void cells_invers(int id_sub_element = NULL_STATE);
		void cells_cyclic_shift(int m = 1);
//...cells including in the surface structure;
		void topo_correct	  (int N = -1,  int N_new = -1, int id_list = OK_STATE);
		int  cells_id		  (CCells * ce, int id_num_correct = NULL_STATE);
		int  cells_in_struct(CCells * ce, int id_num_correct = NULL_STATE);
		void delete_cells_in_struct();
		void search_cells_in_struct(CCells * add_ce, int & i, int id_cell, CCells **& dop_ce, int & N_buf, int & buf_size, int buf_delta = 20);
		void search_cells_element	(int *& id, int & N_buf, int & buf_size, int buf_delta = 20);
		int   bar_add(CCells * ext_ce, int id_cell = OK_STATE, int buf_delta = 20);
		void trim_add(CCells *& trim,  int id_mp   = OK_STATE);
//...cells extraction from surface geometry;
		inline CMap * map_cpy(int N) {
			if (0 <= N && N < graph[0]) return(::map_cpy(ce[N]->mp));
			else                        return(NULL);
		}
		inline CMap ** pm_cpy(int N) {
			if (0 <= N && N < graph[0] && ce[N] && ce[N]->graph && ce[N]->pm)
				return(::pm_cpy(ce[N]->pm, ce[N]->graph[1])); else
				return(NULL);
		}
		inline Topo  * graph_cpy(int N) {
			if (N >= graph[0] || N < 0) return(NULL);
			Topo * new_graph = (Topo *)new_struct((ce[N]->graph[1]+2)*sizeof(Topo));
			if ( ! new_graph) return(new_graph);
			memcpy(new_graph, ce[N]->graph, (ce[N]->graph[1]+2)*sizeof(Topo));
			return(new_graph);
		}
		int		bar_span(CMap *& ext_mp);
		CCells * bar_cpy(int N, int id_origine = NULL_STATE, int buf_delta = 20);
		CCells * bar_sub(int N, int id_origine = NULL_STATE, int buf_delta = 20);
		CCells * bar_cpy();
		void		bar_exc(int N);
		void install_struct();
		void  bar_generate				(double eps = EE_ker);
		void  bar_correct_double_point(double eps = EE_ker);
		void cell_correct_double_point(double eps = EE_ker);
protected:
//...identification of geometrical segments;
		int id_sheet();
		int id_ring_segment();
		int id_cyl_segment();
		int id_sph_segment();
		int id_spr_segment();
		int id_cone_segment();
		int id_torus_segment();
		int id_ugolok_cell();
		int id_curve1_tria();
		int id_curve2_tria();
		int id_curve1_quad();
		int id_curve2_quad();
		int id_curve11quad();
		int id_curve3_quad();
		int id_curve4_quad();
public:
		int segms_id ();
		int segms_bar();
protected:
//...uniform grid for one dimension cells;
		int  grid_line     (CGrid * nd, double h, int Max_N, int hit = 1);
		int  grid_circ     (CGrid * nd, double h, int Max_N, int hit = 0);
		int  grid_ellipt   (CGrid * nd, double h, int Max_N, int hit = 1);
public:
		int  grid_cells1   (CGrid * nd, double h = 0., int Max_knee = 20);
		int  grid_skeleton (CGrid * nd, double h = 0., int Max_knee = 10, int id_full = 0, Topo * link = NULL);
protected:
//...cells interior;
		int  in_plane_facet(double * P, double *  pm, int   N_arc, int k1, int k2);
		int  in_plane_facet(double * P, double *& pm, int & N_ini, int * mm, int id_fast);
		int  in_half_plane (double X, double Y, double Z, int id_local = OK_STATE, double eps = EE_ker);
public:
		int  in_poly_body	(double X, double Y, double Z, int id_local = NULL_STATE, int * mm = NULL, double eps = EE_ker);
		int  in_poly_facet(double X, double Y, double Z, int id_fast = OK_STATE, double eps = EE_ker);
		int  in_bar_facet	(double X, double Y, double Z, int id_fast = ZERO_STATE);
		int  in_bar			(double X, double Y);
		int  in_bar_cell	(double X, double Y);
//...auxilliary functions for lines and circular arcs;
		inline complex get_arc_center() {
			if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0]) return comp(mp[1], mp[2]);
			else return comp(0.);
		}
		inline complex get_arc_point(int m) {
			if (1 == cells_dim()) {
				if (m >= graph[1]) m = graph[1]-1;
				if (m <  0)        m = 0;
				return comp(ce[graph[2+m]]->mp[1], ce[graph[2+m]]->mp[2]);
			}
			return comp(0.);
		}
		inline complex get_arc_point(int m, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX) {
			if (1 == cells_dim()) {
				if (m >= graph[1]) m = graph[1]-1;
				if (m <  0)        m = 0;
				double p[3] = {ce[graph[2+m]]->mp[1], ce[graph[2+m]]->mp[2], ce[graph[2+m]]->mp[3]};

				point_iso<double>(p, NULL, CZ, SZ, CY, SY, CX, SX);
				return comp(p[0], p[1]);
			}
			return comp(0.);
		}
		inline void set_arc_point (int m, complex z) {
					if (1 == cells_dim()) {
						if (m >= graph[1]) m = graph[1]-1;
						if (m <  0)        m = 0;
						ce[graph[2+m]]->mp[1] = real(z);
						ce[graph[2+m]]->mp[2] = imag(z);
						ce[graph[2+m]]->mp[3] = 0.;
					}
				}
		inline void set_arc_center(complex z) {
			if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0]) {
				mp[1] = real(z);
				mp[2] = imag(z);
				mp[3] = 0.;
			}
		}
		inline void line_correct(int id_fast = NULL_STATE) {
			if (mp && ID_MAP(1, NULL_GENUS) == mp[0] && graph[1] > 0) {
				double X, Y, Z, fi,
					P1[3] = { X = ce[graph[2]]->mp[1], Y = ce[graph[2]]->mp[2], Z = ce[graph[2]]->mp[3]},
					P2[3] = { X, Y, Z };

				if (graph[1] != 1) {
					P2[0] = ce[graph[3]]->mp[1];
					P2[1] = ce[graph[3]]->mp[2];
					P2[2] = ce[graph[3]]->mp[3];
				}
				X -= P2[0];
				Y -= P2[1];
				Z -= P2[2];
		////////////////////////////////////
		//...коррекци€ геометрической карты;
				mp[1] = (CMap)(P1[0]+P2[0])*.5;
				mp[2] = (CMap)(P1[1]+P2[1])*.5;
				mp[3] = (CMap)(P1[2]+P2[2])*.5;
				mp[4] = (CMap)0.;
				mp[5] = (CMap)0.;
				mp[6] = (CMap)0.;
				if (id_fast == NULL_STATE) {
					mp[4] = (CMap)(arg0(comp(-X, -Y)));
					mp[5] = (CMap)(arg0(comp(-Z, sqrt(X*X+Y*Y))));
					mp[6] = (CMap)(0.);
				}
				else {
					mp[4] = 0.;
					mp[5] = 0.;
					mp[6] = 1.;
					if ((fi = sqrt(sqr(X)+sqr(Y)+sqr(Z))) > EE) {
						mp[4] = (CMap)(-X*(fi = 1./fi));
						mp[5] = (CMap)(-Y* fi);
						mp[6] = (CMap)(-Z* fi);
					}
				}
			}
		}
		inline void circ_correct(double R, int topo, int points_inverse = NULL_STATE) {
			if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0] && graph[1] > 0) {
				if (! topo) {
					mp[0] = ID_MAP(1, NULL_GENUS); line_correct();
				}
				else {
					double X, Y, Z, fi, theta, beta, sg, f0, f,
							 P1[3] = { X = ce[graph[2]]->mp[1], Y = ce[graph[2]]->mp[2], Z = 0.},
							 P2[3] = { X, Y, Z };

					if (graph[1] != 1) {
						 P2[0] = ce[graph[3]]->mp[1];
						 P2[1] = ce[graph[3]]->mp[2];
						 P2[2] = Z;
					}
					if (points_inverse) {
						swap(P1[0], P2[0]); X = P1[0];
						swap(P1[1], P2[1]); Y = P1[1];
						swap(P1[2], P2[2]); Z = P1[2];
					}
					X -= P2[0];
					Y -= P2[1];
					Z -= P2[2];
					R  = (f = sqrt(X*X+Y*Y+Z*Z)) < EE_ker ? 0. : fabs(R);
					f0 = R == 0. ? 0. : R*R-f*f*.25;

					if ((f0 = filtr_(f0, EE_ker)) < 0.) {
						 mp[0] = ID_MAP(1, NULL_GENUS); line_correct();
						 return;
					}
					beta  = arg0(comp(sqrt(f0), f*.5)),
					fi    = arg0(comp(-X, -Y))-M_PI_2,
					sg    = topo > 0 ? -1. : 1.,
					theta = topo > 0 ?  0. : M_PI; if (2 == abs(topo)) beta = M_PI-beta;

					mp[1] = (CMap)(P1[0]+sg*R*cos(fi+sg*beta));
					mp[2] = (CMap)(P1[1]+sg*R*sin(fi+sg*beta));
					mp[5] = (CMap)theta;
					mp[6] = (CMap)(cos(theta)*fi); //...устанавливаем "доворот" в mp[6];
					mp[7] = (CMap)R;
					mp[8] = (CMap)beta;
					mp[3] =
					mp[4] = (CMap)(fi = 0.);
				}
			}
		}
protected:
//...quadratures for dimensional cells;
		void curve1_tria_QG(CGrid * nd, int N_elem);
		void curve2_tria_QG(CGrid * nd, int N_elem);
		void curve1_quad_QG(CGrid * nd, int N_elem);
		void curve2_quad_QG(CGrid * nd, int N_elem);
		void curve11quad_QG(CGrid * nd, int N_elem);
		void curve3_quad_QG(CGrid * nd, int N_elem);
		void curve4_quad_QG(CGrid * nd, int N_elem);
public:
		void facet_QG(CGrid * nd, int N_elem, int N_max = 1, double * pp_cond = NULL, int id_fast = OK_STATE);
		void segms_QG(CGrid * nd, int N_elem, int N_max = 1, double * pp_cond = NULL, int id_fast = OK_STATE);
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
//...special compositions;
		void get_circle_profile	(double r);
      void get_ugolok			(double beta, double radc, double rad, double A1, double B1, double A2, double B2);
      void get_ugolok_cell		(double A1, double A2, double B1, double B2, double radc,
										 double rad, double beta, double x0 = 0., double y0 = 0., double f0 = 0.);
//...geometrical segments;
      void get_blend_cyl_segment	(double beta, double R, double L, int * mm = NULL);
      void get_torus_segment		(double r0,    double r1, double f0, double f1, double t0, double t1);
      void get_ring_segment		(              double f0, double f1, double t0, double t1);
      void get_cone_segment		(double theta, double f0, double f1, double t0, double t1);
      void get_cyl_segment			(double R,     double f0, double f1, double t0, double t1);
      void get_sph_segment			(double R,     double f0, double f1, double t0, double t1);
//...plane and space composition of cells;
      void get_polygon_directly    (double * X, double * Y, int N, int id_dbl = 0, double eps = EE_ker);
      void get_arc_polygon_directly(double * X, double * Y, double * R, int * topo, int N, int id_dbl = 0, double eps = EE_ker);
      void get_line_strip_directly (double * X, double * Y, double * Z, int * geom, int * mask = NULL, int shift = 2, int id_fast = OK_STATE, int id_dbl = NULL_STATE, double eps = EE_ker);
      void get_beam                (CCells * f0, double length);
      void get_cyl_beam            (CCells * f0, double fi);
//...plane facet compositions;
      void get_facet_directly(double * P, int N, int id_plane = NULL_STATE, int id_fast = NULL_STATE, int id_dbl_point = NULL_STATE, double eps = EE_ker);
      void get_tria_facet(double * P1, double * P2, double * P3, int id_fast = OK_STATE);
      void get_quad_facet(double * P1, double * P2, double * P3, double * P4, int id_fast = OK_STATE);
//...plane facet and node compositions;
      int  SetFacetParam (int k, double square, double * pp = NULL, int id_property = DEFAULT_BND);
      void SetFacetXParam(int k, double X0, double * pp, int id_property) {
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[4])-1.) < EE_ker && fabs(ce[k]->mp[1]-X0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
      void SetFacetYParam(int k, double Y0, double * pp, int id_property){
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[5])-1.) < EE_ker && fabs(ce[k]->mp[2]-Y0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
      void SetFacetZParam(int k, double Z0, double * pp, int id_property) {
			if (ce[k]->mp[8] > 0. && fabs(fabs(ce[k]->mp[6])-1.) < EE_ker && fabs(ce[k]->mp[3]-Z0) < EE_ker)
				SetFacetParam(k, ce[k]->mp[8], pp, id_property);
		}
 		int  SetNodeParam	  (double * pp = NULL, int id_property = DEFAULT_BND);
		void SetNodeXParam  (double X0, double * pp, int id_property);
		void SetNodeYParam  (double Y0, double * pp, int id_property);
		void SetNodeZParam  (double Z0, double * pp, int id_property);
		void SetNodeXYParam (double X0, double Y0, double * pp, int id_property);
		void SetNodeXZParam (double X0, double Z0, double * pp, int id_property);
		void SetNodeYZParam (double Y0, double Z0, double * pp, int id_property);
		void SetNodeXYZParam(double X0, double Y0, double Z0, double * pp, int id_property);
//...polygonal compositions and boundary conditions;
      void get_nd_bar_directly(CGrid * nd, int k);
      int  set_nd_bar_condit(int k, CGrid * nd, int k1, int k2, int id_property, double * pp, int * id_pp, double & S, int id_facet = -1);
      void set_nd_bar_condit(int k, int id_property, double * pp, int * id_pp, double S);
};
////////////////////////////////////////
//...complete deleting of abstract cell;
inline void delete_cells(CCells *& ce) {
	ce->zero_cells(); delete ce; ce = NULL;
}
#endif
