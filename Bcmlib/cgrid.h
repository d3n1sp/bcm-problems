/*=========================================*/
/*                 CGrid                   */
/*=========================================*/
#ifndef ___CGrid___
#define ___CGrid___

#include "kernel.h"

////////////////////////////////////////////////////////
//...description of existing types of blocks partitions;
enum Num_Nodes {
        NULL_NODES = 0,
     GRID_EL_NODES,
     GRID_QG_NODES,
        NUMS_NODES
};

//////////////////////////////////////////////////////////
//...description of existing types of boundary conditions;
enum Num_Bnd {
		  DEFAULT_BND = 0,
			 RIGID_BND,
//...for CAcou3D;
				VEL_BND,
				 PN_BND,
				 VN_BND,
		  VNTABLE_BND,
			ABSORB_BND,
	 ABSORBTABLE_BND,
			SOURCE_BND,
		  BSOURCE_BND,
		  SOURCE2_BND,
		 BSOURCE2_BND,
			EXTERN_BND,
//...for CLame3D;
		  DISPL3D_BND,
		  FORCE3D_BND,
			DISPLN_BND,
			FORCEN_BND,
//...for filtration;
		  ADHESION_BND,
		  PRESSURE_BND,
		  VELOCITY_BND,
		VELOCITY_N_BND,
//...for homogenization;
			FNORMS_BND,
			FSKEWS_BND,
		  UPTAKES_BND,
//...for CBeam3D;
  RIGID_DISPL3D_BND,
		 MOMENT3D_BND,
//...free boundary condition (signed for special aids);
		  SPECIAL_BND,
		  MARKER1_BND,
		  MARKER2_BND,
		  MARKER3_BND,
		  MARKER4_BND,
		  MARKER5_BND,
		  MARKER6_BND,
		  MARKER7_BND,
		  MARKER8_BND,
			  NUMS_BND,
	//...signs for special conditions (must be less than NUMS_BND+SPECIAL_BND); 
		  MOMENTS_BND,
			 SKEWS_BND,
			 NORMS_BND,
			UPTAKE_BND,
};

////////////////////////////////////////////
//...type of element for OpenGL (and other);
#define GL_POINTS                     0x0000
#define GL_LINES                      0x0001
#define GL_LINE_LOOP                  0x0002
#define GL_LINE_STRIP                 0x0003
#define GL_TRIANGLES                  0x0004
#define GL_TRIANGLE_STRIP             0x0005
#define GL_TRIANGLE_FAN               0x0006
#define GL_QUADS                      0x0007
#define GL_QUAD_STRIP                 0x0008
#define GL_POLYGON                    0x0009
#define GL_BOXS                       0x000A
#define GL_PENTA                      0x000B
#define GL_TETRA                      0x000C
#define GL_PYRAMID                    0x000D
#define GL_SPDC                       0x000E
#define GL_PLOAD4                     0x000F
#define GL_NODES                      0x0010
#define GL_NODES_GENERATE             0x0011
#define GL_ELEMENTS                   0x0012
#define GL_ELEMENTS_GENERATE          0x0013
#define FV_BOXS                       0x0014


///////////////////////////////////////////
//...axis for visualization grid in Surfer;
#define AXIS_X                              1
#define AXIS_Y                              2
#define AXIS_Z                              3
#define AXIS_SPH                            4
#define AXIS_TIME                           5

/////////////////////////////////////////
//...dividing of circle in Bars elements;
#define MAX_KNEE                            10

//////////////////////////////////////////////////////////////////////////////
//...макросы для поиска ключевого слова и для отделения одной строки в потоке;
#define PPOS_CUR(id_STREAM, pchar, ppos_cur, upper_limit, id_STRSTR)\
if ((pchar = strstr(id_STREAM+ppos_cur, id_STRSTR)) != NULL) \
ppos_cur = (ppos_cur = (unsigned long)(pchar-id_STREAM)+strlen(id_STRSTR)) < upper_limit ? ppos_cur : upper_limit; else \
ppos_cur = upper_limit;

#define ONE_LINE(id_STREAM, pchar, ppos_cur, upper_limit, one_line, STR_SIZE)\
if (ppos_cur < upper_limit && (pchar = strstr(id_STREAM+ppos_cur, "\xA")) != NULL) {\
	unsigned long count = (unsigned long)(pchar-id_STREAM)+1;\
	memcpy(one_line, id_STREAM+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));\
	one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0'; ppos_cur = count; \
} else one_line[0] = '\x0';

///////////////////////////////
//...basic class of grid nodes;
class CGrid {
public:
		int N, N1, N2, N3, N_par;									//...number of grid nodes and addiional parameters;
		int * hit, * geom, * geom_ptr, * cond, * cond_ptr; //...appearence of node, topology and boundary conditions;
		double  * X, * Y, * Z, * nX, * nY, * nZ;				//...coordinates and normals of grid nodes;
static int SPTR_BUF, STRU_BUF;
protected:
		double ** pp;	  											  //...table of addiional parameters;
		int buf_count, buf_X, buf_Y, buf_Z, N_buf, N1buf; //...bufferization;
		int bufferization (int buf_size, int mm = OK_STATE);
		int bufferizat_X  (int buf_size, int mm = OK_STATE);
		int bufferizat_Y  (int buf_size, int mm = OK_STATE);
		int bufferizat_Z  (int buf_size, int mm = OK_STATE);
public:
		virtual int type() { return NULL_NODES;}
		void reset_hit	 () { if (hit) memset(hit, -1, N*sizeof(int));}  
//...constructor/destructor;
		CGrid() {
			N = N1 = N2 = N3 = N_par = buf_count = buf_X = buf_Y = buf_Z = 0;
			X = Y = Z =	nX = nY = nZ = NULL;
			geom  = geom_ptr = NULL;
			cond  = cond_ptr = NULL;
			hit   = NULL;
			pp    = NULL; 
			N_buf = N1buf = 2000;
		};
      virtual ~CGrid (void) {
			zero_grid();
		}
//...quadratures;
		virtual void facet_QG(double * P, int N_elem, int mode = NULL_STATE, int normal = OK_STATE) {}
		virtual void segms_QG(CMap * mp, int N_elem, int N_max) {}
		virtual void sphere_intrusion_QG(CMap * mp, int N_elem, double L) {}
//...quadratures corrections;
		virtual void QG_curve  (CMap * mp) {}
		virtual void QG_surface(CMap * mp) {}
		virtual void QG_tria_curve(CMap * mp, double * P, int N_ini = 0) {}
		virtual void QG_quad_curve(CMap * mp1, CMap * mp2, double * Po, int N_ini = 0) {}
		virtual void QG_tria_surface(CMap * mp, double * Po, int N_ini = 0) {}
		virtual void QG_quad_surface(CMap * mp1, CMap * mp2, double * Po, int N_ini = 0) {}
		virtual void QG_quad_bi(double * Po, int N_ini = 0) {}
//...elements subdivision;
		virtual void grid_box1D(int N1, int shift = 0, int id_link = 0, double fy = 0., double fz = 0.) {}
		virtual void grid_box2D(int N1, int N2, int shift = 0, int id_link = 0, double fz = 0.) {}
		virtual void grid_box3D(int N1, int N2, int N3, int shift = 0, int id_link = 0) {}
//...geometrical map interior;
		virtual int maps_in (CMap * mp, double * P, double eps = EE) { return(0);}
		virtual int segms_in(CMap * mp, double * P, double eps = EE) { return(0);}
public:
		void zero_grid();
		int  add_new_point  (double X, double Y, double Z, double nX = 0., double nY = 0., double nZ = 1., double * pp = NULL);
		int  add_new_point  (double * P, double * pp = NULL);
		int  add_new_point_X(double val, double eps = EE_ker);
		int  add_new_point_Y(double val, double eps = EE_ker);
		int  add_new_point_Z(double val, double eps = EE_ker);
		void set_buf_size	  (int buf_size) { N_buf = buf_size;}
		void add_buffer     (int m) {
			N -= (m = min(max(0, m), N));
			buf_count += m;
		}
		void add_params (int m = 1, int id_bufferizat = OK_STATE);
		void set_param  (double * pp);
		void set_param  (int m, int k, double value) {
			if (0 <= m && m < N_par && pp[m] && 0 <= k && k < N)
				 pp[m][k] = value;
		}
		double get_param(int m, int k) {
			if (0 <= m && m < N_par && pp[m] && 0 <= k && k < N) return pp[m][k];
			else return(0.);
		}
		void grid_add(CGrid * ext_nd);
		int  grid_add(double *& X0, double *&Y0,double *& Z0, double *& nX0, double *& nY0, 
											 double *& nZ0, int N0, int *& gm, int id_block);
		void grid_iso(int k0, int k1, double * P0,
												double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
		void grid_tria_3D_refine(int N_elem, double * P, double * pp = NULL, int id_block = -1);
		void grid_quad_3D_refine(int N_elem1, int N_elem2, double * P, double * pp = NULL, int id_block = -1);
		void grid_quad_3D_refine_spec(int N_elem1, int N_elem2, double * P, double * pp = NULL, int id_block = -1);
//...operations with geometry;
		void init_geom   (int m = 1) { delete_struct(geom); geom = (int *)new_struct((m+1)*sizeof(int)); };
		int  stru_install(int k, int ** mm_facet);
		int  stru_shift  (int shift_num = 0);
		int  link_number (Topo * link, int id_arc);
		int  geom_element(int k);
		int  geom_add_ptr(int N, int type = GL_POINTS) {
			if (! geom_ptr) {//...инициализация;
					geom_ptr = (int *)new_struct(sizeof(int)*2); geom_ptr[0] = 1;
					delete_struct(geom); geom = (int *)new_struct(sizeof(int)*2);
			}
			if (! geom_ptr[geom[0]+1]) {//...буфферизация указателя;
			int *	geom_ptr_new = (int *)new_struct((geom[0]+2)*sizeof(int));
					memcpy(geom_ptr_new, geom_ptr, (geom[0]+1)*sizeof(int)); geom_ptr_new[geom[0]+1] = SPTR_BUF;
					delete_struct(geom_ptr); geom_ptr = geom_ptr_new;
			}
			if (  geom[geom_ptr[geom[0]]] < 2+N) {//...буфферизация геометрии;
					int buff = max (STRU_BUF, 2+N), * geom_new;
					geom_new = (int *)new_struct((geom_ptr[geom[0]]+buff+1)*sizeof(int));
					memcpy(geom_new, geom, geom_ptr[geom[0]]*sizeof(int)); geom_new[geom_ptr[geom[0]]] = buff;
					delete_struct(geom); geom = geom_new;
			}
			geom_ptr[geom[0]+2] = geom_ptr[geom[0]+1]-1;
			geom_ptr[geom[0]+1] = geom_ptr[geom[0]]+N+2;
			geom[geom_ptr[geom[0]+1]] = geom[geom_ptr[geom[0]]]-N-2;
			geom[geom_ptr[geom[0]]]	  = type;
			geom[geom_ptr[geom[0]]+1] = N;
			geom[0] += 1;
			return(geom_ptr[geom[0]-1]);
		};
		void geom_ptr_add(int N, int & N_geom) {
			if (! geom_ptr) {//...инициализация;
					geom_ptr = (int *)new_struct(sizeof(int)*2); 
					geom_ptr[N_geom = 0] = 1;
			}
			if (N >= 0) {
				if (! geom_ptr[N_geom+1]) {//...буфферизация указателя;
				int *	geom_ptr_new = (int *)new_struct((N_geom+SPTR_BUF+2)*sizeof(int));
						memcpy(geom_ptr_new, geom_ptr, (N_geom+1)*sizeof(int)); geom_ptr_new[N_geom+1] = SPTR_BUF;
						delete_struct(geom_ptr); geom_ptr = geom_ptr_new;
				}
				geom_ptr[N_geom+2] = geom_ptr[N_geom+1]-1;
				geom_ptr[N_geom+1] = geom_ptr[N_geom]+N+2; N_geom++;
			}
		};
		int  cond_add_ptr(int N, int type = GL_POINTS) {
			if (! cond_ptr) {//...инициализация;
					cond_ptr = (int *)new_struct(sizeof(int)*2); cond_ptr[0] = 1;
					delete_struct(cond); cond = (int *)new_struct(sizeof(int)*2);
			}
			if (! cond_ptr[cond[0]+1]) {//...буфферизация указателя;
			int *	cond_ptr_new = (int *)new_struct((cond[0]+SPTR_BUF+2)*sizeof(int));
					memcpy(cond_ptr_new, cond_ptr, (cond[0]+1)*sizeof(int)); cond_ptr_new[cond[0]+1] = SPTR_BUF;
					delete_struct(cond_ptr); cond_ptr = cond_ptr_new;
			}
			if (  cond[cond_ptr[cond[0]]] < 2+N) {//...буфферизация геометрии;
					int buff = max (STRU_BUF, 2+N), * cond_new;
					cond_new = (int *)new_struct((cond_ptr[cond[0]]+buff+1)*sizeof(int));
					memcpy(cond_new, cond, cond_ptr[cond[0]]*sizeof(int)); cond_new[cond_ptr[cond[0]]] = buff;
					delete_struct(cond); cond = cond_new;
			}
			cond_ptr[cond[0]+2] = cond_ptr[cond[0]+1]-1;
			cond_ptr[cond[0]+1] = cond_ptr[cond[0]]+N+2;
			cond[cond_ptr[cond[0]+1]] = cond[cond_ptr[cond[0]]]-N-2;
			cond[cond_ptr[cond[0]]]	  = type;
			cond[cond_ptr[cond[0]]+1] = N;
			cond[0] += 1;
			return(cond_ptr[cond[0]-1]);
		};
		void cond_ptr_add(int N, int & N_cond) {
			if (! cond_ptr) {//...инициализация;
					cond_ptr = (int *)new_struct(sizeof(int)*2); 
					cond_ptr[N_cond = 0] = 1;
			}
			if (N >= 0) {
				if (! cond_ptr[N_cond+1]) {//...буфферизация указателя;
				int *	cond_ptr_new = (int *)new_struct((N_cond+SPTR_BUF+2)*sizeof(int));
						memcpy(cond_ptr_new, cond_ptr, (N_cond+1)*sizeof(int)); cond_ptr_new[N_cond+1] = SPTR_BUF;
						delete_struct(cond_ptr); cond_ptr = cond_ptr_new;
				}
				cond_ptr[N_cond+2] = cond_ptr[N_cond+1]-1;
				cond_ptr[N_cond+1] = cond_ptr[N_cond]+N+2; N_cond++;
			}
		};
		void quick_sort			  (int first, int last, int * index);
		void quick_sort_inverse	  (int first, int last, int * index);
		int  binary_search		  (int first, int last, int elem, int * index);
		int  binary_search_inverse(int first, int last, int elem, int * index);
//...nodes reading from NASTRAN and ABAQUS input file;
		int element_type(int & elem, char * id_NODES = NULL) {
			int quad = 0;
			if (id_NODES) { //...,определяем тип элемента и его квадратичную характеристику; 
				if (! ::strncmp(id_NODES, "B31",   3) ||
					 ! ::strncmp(id_NODES, "B21",   3) ||
					 ! ::strncmp(id_NODES, "T3D2",  4)) elem = GL_LINE_STRIP; else
				if (! ::strncmp(id_NODES, "CPE3",  4) ||
					 ! ::strncmp(id_NODES, "CPS3",  4) || 
					 ! ::strncmp(id_NODES, "S3",	  2)) elem = GL_TRIANGLES; else
				if (! ::strncmp(id_NODES, "CPE4",  4) || 
					 ! ::strncmp(id_NODES, "CPS4",  4) || 
					 ! ::strncmp(id_NODES, "S4",	  2)) elem = GL_QUADS; else 
				if (! ::strncmp(id_NODES, "C3D4",  4)) elem = GL_TETRA; else 
				if (! ::strncmp(id_NODES, "C3D8",  4)) elem = GL_BOXS; else 
				if (! ::strncmp(id_NODES, "C3D6",  4)) elem = GL_PENTA; 
				else {
					quad = 1;
					if (! ::strncmp(id_NODES, "B32",	  3) ||
						 ! ::strncmp(id_NODES, "B22",   3)) elem = GL_LINE_STRIP; else
					if (! ::strncmp(id_NODES, "CPE6",  4) || 
						 ! ::strncmp(id_NODES, "CPS6",  4) || 
						 ! ::strncmp(id_NODES, "S6",	  2)) elem = GL_TRIANGLES; else
					if (! ::strncmp(id_NODES, "CPE8",  4) || 
						 ! ::strncmp(id_NODES, "CPS8",  4) || 
						 ! ::strncmp(id_NODES, "S8",	  2)) elem = GL_QUADS; else 
					if (! ::strncmp(id_NODES, "C3D10", 5)) elem = GL_TETRA; else 
					if (! ::strncmp(id_NODES, "C3D20", 5)) elem = GL_BOXS; else 
					if (! ::strncmp(id_NODES, "C3D15", 5)) elem = GL_PENTA; else elem = ERR_STATE;
				}
			}
			int elem_params = 0;
			if (elem != ERR_STATE) { //...,определяем число параметров для зачитывания элемента;
				if (elem == GL_LINE_STRIP) elem_params = quad ?  7 : 6; else
				if (elem == GL_TRIANGLES)  elem_params = quad ? 10 : 7; else
				if (elem == GL_QUADS) elem_params = quad ? 12 : 8; else 
				if (elem == GL_TETRA) elem_params = quad ? 14 : 8; else 
				if (elem == GL_BOXS)  elem_params = quad ? 24 : 12; else 
				if (elem == GL_PENTA) elem_params = quad ? 19 : 10;
			}
			return elem_params;
		}
		int  converts_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long);
		void converts_nas(char * ch_NODES, int id_long = 0);
		void converts_nas(const char * ch_NODES, int id_long = 0);

		int  condit_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long);
		void condit_nas(char * ch_NODES, int id_long = 0);
		void condit_nas(const char * ch_NODES, int id_long = 0);

		int  nodes_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long, int id_status = OK_STATE);
		void nodes_nas(char * ch_NODES, int id_long = 0);
		void nodes_nas(const char * ch_NODES, int id_long = 0);

		int  converts_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit, int numeration_shift = 0);
		void converts_inp(char * ch_NODES);
		void converts_inp(const char * ch_NODES);

		int  condit_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit, int numeration_shift = 0, int max_phase = 2, int id_status = OK_STATE);
		void condit_inp(char * ch_NODES);
		void condit_inp(const char * ch_NODES);

		int  nodes_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit);
		void nodes_inp(char * ch_NODES, int numeration_shift = 0, int max_phase = 2);
		void nodes_inp(const char * ch_NODES, int numeration_shift = 0, int max_phase = 2);

		void nodes_in (char * ch_NODES, int max_phase = 2) {
			char * pchar = ::strrchr(ch_NODES, '.');
			if (pchar && (! ::strncmp(pchar, ".nas", 4) || ! ::strncmp(pchar, ".NAS", 4))) nodes_nas(ch_NODES); else
			if (pchar && (! ::strncmp(pchar, ".inp", 4) || ! ::strncmp(pchar, ".INP", 4))) nodes_inp(ch_NODES, 1, max_phase);
		};
		void nodes_in (const char * ch_NODES, int max_phase = 2) {
			const char * pchar = ::strrchr(ch_NODES, '.');
			if (pchar && (! ::strncmp(pchar, ".nas", 4) || ! ::strncmp(pchar, ".NAS", 4))) nodes_nas(ch_NODES); else
			if (pchar && (! ::strncmp(pchar, ".inp", 4) || ! ::strncmp(pchar, ".INP", 4))) nodes_inp(ch_NODES, 1, max_phase);
		};
//...testing functions;
		void set_point(double * P, int k,
							double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX) {
			P[0] = X ? X[k] : 0.;
			P[1] = Y ? Y[k] : 0.;
			P[2] = Z ? Z[k] : 0.;
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
		}
		void set_point(double * P, int i, int j, int l,
							double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX) {
			P[0] = X ? X[i] : 0.;
			P[1] = Y ? Y[j] : 0.;
			P[2] = Z ? Z[l] : 0.;
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
		}
		void set_point(double * P, int k, int dir,
							double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX) {
			double A = nX ? nX[k] : 0.,
					 B = nY ? nY[k] : 0.,
					 C = nZ ? nZ[k] : 0.;
			P[0] = X ? X[k] : 0.;
			P[1] = Y ? Y[k] : 0.;
			P[2] = Z ? Z[k] : 0.;
			switch (dir) {
				case 1: P[0] -= A; P[1] -= B; P[2] -= C; break;
				case 2: P[0] += A; P[1] -= B; P[2] -= C; break;
				case 3: P[0] += A; P[1] += B; P[2] -= C; break;
				case 4: P[0] -= A; P[1] += B; P[2] -= C; break;
				case 5: P[0] -= A; P[1] -= B; P[2] += C; break;
				case 6: P[0] += A; P[1] -= B; P[2] += C; break;
				case 7: P[0] += A; P[1] += B; P[2] += C; break;
				case 8: P[0] -= A; P[1] += B; P[2] += C; break;
			}
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
		}
		void grid_out (FILE * GR, double fi, double theta, double fX, int m_axis, int shift, int m_contour, int element);
		void grid_out (char * GR_FILE, double fi = 0., double theta = 0.,
												 double fX = 0., int m_axis = AXIS_Z, int shift = 0, int m_contour = 0, int element = 0);
		void grid_out (const char * GR_FILE, double fi = 0., double theta = 0.,
												 double fX = 0., int m_axis = AXIS_Z, int shift = 0, int m_contour = 0, int element = 0);
		void TestGrid (FILE * GR_FILE, double mtb, double fi, double theta, double fX, int m_axis, int id_normal);
		void TestGrid (char * GR_FILE, double mtb, double fi = 0., double theta = 0.,
																 double fX = 0., int m_axis = AXIS_Z, int id_normal = 0);
		void TestGrid (const char * GR_FILE, double mtb, double fi = 0., double theta = 0.,
																 double fX = 0., int m_axis = AXIS_Z, int id_normal = 0);
};
/////////////////////////////////////////
//...complete deleting of abstract nodes;
inline void delete_nodes(CGrid *& nd) {
	delete nd; nd = NULL;
}
/////////////////////////////////////////////////////
//...global function for construction all grid nodes;
CGrid * CreateNodes(Num_Nodes id_NODES = NULL_NODES);
#endif
