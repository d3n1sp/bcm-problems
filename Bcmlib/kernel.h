/*============================================*/
/*                  KERNEL                    */
/*============================================*/
#ifndef ___KERNEL___
#define ___KERNEL___

#include "utils.h"

#define  MAX_DIM               4
#define  ERR_DIM              -1
#define  ERR_GRAPH            -1
#define  MAX_SPLINE            9
#define  ID_MAP(map_dim, map_genus)  ((int)map_dim+(int)map_genus*MAX_DIM)

///////////////////////////////////////////////////////////////////////////
//...список кодов элементов геометрического ядра, используемых в программе;
//enum Num_Genus {
const int
     ERR_GENUS          = -1,
     NULL_GENUS         =  0,
     CYL_GENUS          =  1,
     ELL_CYL_GENUS      =  2,
     SPHERE_GENUS       =  3,
     SPHEROID_GENUS     =  4,
     ELLIPSOID_GENUS    =  5,
     CONE_GENUS         =  6,
     ELL_CONE_GENUS     =  7,
     TORUS_GENUS        =  8,
     PARAB_GENUS        =  9,
     ELL_PARAB_GENUS    = 10,
     BI_HYP_GENUS       = 11,
     ELL_BI_HYP_GENUS   = 12,
     HYPERB_GENUS       = 13,
     ELL_HYPERB_GENUS   = 14,
     NURBS_GENUS        = 15,
     NURBS_NULL_GENUS   = 16,
     NURBS_CYL_GENUS    = 17,
     BEZIER_GENUS       = 18, 
     WEDGE_GENUS        = 19;
//};
typedef int    Topo;
typedef double CMap;

////////////////////////////////////////////////////////////////////
//...основной объект геометрического ядра - пространственная ячейка;
struct Cells {
			Cells ** ce;    //...общий список всех элементов ячейки (включая себя);
			Topo  *  graph; //...список связей ячейки;
			CMap  *  mp;    //...каpта геометpии ячейки;
			CMap  ** pm;    //...каpта геометpии ячейки в пространстве параметров;
};
typedef  Cells Bar;		 //...поверхность тела (mp == NULL, число элементов поверхности = graph[0]);

///////////////////////////////////
//...идентификация размерных ячеек;
enum Num_Cell {
     NULL_CELL				= -1,
///////////////////////////////////
//...некоторые общепринятые ячейки;
     SHEET_CELL			= -2,
     BOX_CELL				= -3,
     UGOLOK_CELL			= -4,
     FACET_CELL			= -5,
     SPH_INTRUSION_CELL = -6,
     SHT_INTRUSION_CELL	= -7,
     CURVE1_TRIA_CELL	= -8,
     CURVE2_TRIA_CELL	= -9,
     CURVE1_QUAD_CELL	= -10,
     CURVE2_QUAD_CELL	= -11,
     CURVE11QUAD_CELL	= -12,
     CURVE3_QUAD_CELL	= -13,
     CURVE4_QUAD_CELL	= -14,
     B1SPLINE_CELL		= -15,
     B2SPLINE_CELL		= -16,
     B1POLY_CELL			= -17,
//////////////////////////////////////////////
//...узловые элементы с граничными значениями;
     BOUNDARY_NODE		= -50,
//////////////////////////
//...специфические кривые;
     ELLIPT_ARC			= -101,
////////////////////////////////////////
//...ячейки в прямоугольных координатах;
     RING_SEGMENT			= -102,
     CYL_SEGMENT			= -103,
     SPH_SEGMENT			= -104,
     SPR_SEGMENT			= -105,
     CONE_SEGMENT			= -106,
     TORUS_SEGMENT		= -107,
//////////////////////
//...ячейки - шапочки;
     SPH_BLEND				= -108,
     SPH_BLEND_DOP		= -109,
     SHEET_BLEND			= -110,
     SHEET_BLEND_DOP		= -111,
/////////////////////////////////////////////////////
//...специальные ячейки для задания геометрии блоков;
     SPECIAL_SHEET		= -112,
     SPECIAL_BOX			= -113,
     BASIS_SHEET			= -114,
     BASIS_BOX				= -115,
};

//////////////////////////////////////////
//...извлечение рода геометрической карты;
inline int map_genus(CMap * mp)
{
	return mp ? (int)mp[0] / MAX_DIM : ERR_GENUS;
}

/////////////////////////////////////////////////
//...извлечение размерности геометрической карты;
inline int map_dim(CMap * mp)
{
	return mp ? (int)mp[0] % MAX_DIM : ERR_DIM;
}

/////////////////////////////////////////////////////////////////////////////////////////
//...таблица размеров геометрической карты пространственной ячейки в зависимости от рода;
inline int size_of_map(int dim, int genus)
{
	switch (genus) {
		case         NULL_GENUS: return 0 == dim ? 4  : 7;  //...прямая;
		case          CYL_GENUS: return 1 == dim ? 9  : 10; //...эллипс;
		case       SPHERE_GENUS:                            //...окpужность;
		case         CONE_GENUS: return 1 == dim ? 9  : 8;  //...гипеpбола;
		case     ELL_CONE_GENUS: return 1 == dim ? 8  : 10; //...паpабола;
		case        PARAB_GENUS: return 8;
		case        TORUS_GENUS:
		case       BI_HYP_GENUS:
		case       HYPERB_GENUS:
		case      ELL_CYL_GENUS:
		case     SPHEROID_GENUS:
		case    ELL_PARAB_GENUS: return 9;
		case    ELLIPSOID_GENUS:
		case   ELL_BI_HYP_GENUS:
		case   ELL_HYPERB_GENUS:
		case   NURBS_NULL_GENUS:                            //...цилиндрическая поверхность;
		case    NURBS_CYL_GENUS: return 10;                 //...осесимметричная поверхность;
		case        NURBS_GENUS: return 1 == dim ? 10 : 12; //...B-кривая или B-поверхность;
		case       BEZIER_GENUS: return 8;
		case        WEDGE_GENUS: return 9;
		default                : return 0;                  //...неизвестная геометрическая карта;
	}
}

inline int size_of_map(CMap * mp)
{
	return size_of_map(map_dim(mp), map_genus(mp));
}

///////////////////////////////////////////////////////////////////////////
//...число дополнительных параметров геометрической карты размерной ячейки;
inline int size_of_dop(int id_CELL, CMap * mp = NULL)
{
	switch (id_CELL) {
		case			 NULL_CELL : return 0;
		case			SHEET_CELL : return 2;
		case SHT_INTRUSION_CELL :
		case			  BOX_CELL : return 3;
		case		  UGOLOK_CELL : return 7;
		case			FACET_CELL : 
		case	CURVE1_TRIA_CELL : 
		case	CURVE2_TRIA_CELL : return 5; //...square and 4 parameters of boundary condition;
		case SPH_INTRUSION_CELL :
		case		BOUNDARY_NODE : return 4; //...4 parameters of boundary condition;
		case		B1SPLINE_CELL : return mp ? ((int)mp[7]/2)*(int)mp[9] +(int)mp[8] +2 : 0;
		case		B2SPLINE_CELL : return mp ? ((int)mp[7]/4)*(int)mp[10]*(int)mp[11]+(int)mp[8]+(int)mp[9]+4 : 0;
		case		  B1POLY_CELL : return mp ? ((int)mp[7]/2)*(int)mp[9] : 0;
		case			ELLIPT_ARC : return 2;
		case		  CYL_SEGMENT : return 2;
		case		 RING_SEGMENT :
		case		  SPH_SEGMENT :
     case		  SPR_SEGMENT :
		case		 CONE_SEGMENT :
		case		TORUS_SEGMENT :
		case		  SHEET_BLEND : return 3;
		case			 SPH_BLEND : return 4;
		case	 SHEET_BLEND_DOP : return 2;
		case		SPH_BLEND_DOP : return 3;
		case		  BASIS_SHEET :
		case		SPECIAL_SHEET : return 4;
		case			 BASIS_BOX :
		case		  SPECIAL_BOX : return 6;
		default					  : return 0;  //...неизвестная дополнительная часть геометрической карты;
	}
}

inline int size_of_dop(CMap * mp)
{
	return size_of_dop((int)mp[size_of_map(mp)], mp);
}

////////////////////////////////////////////////
//..."генус" (или род) пространственной ячейки;
inline int cells_genus(Cells * ce)
{
	return ce ? map_genus(ce->mp) : ERR_GENUS;
}

/////////////////////////////////////////
//...размерность пространственной ячейки;
inline int cells_dim(Cells * ce)
{
	return ce ? map_dim(ce->mp) : ERR_DIM;
}

//////////////////////////////////////////////////////////
//...распределение пустой ячейки (с нулевыми параметрами);
inline Cells * new_cells()
{
  return(Cells *)new_struct(sizeof(Cells));
}

/////////////////////////////////////////////////////////
//...уничтожение внутренней структуры абстрактной ячейки;
void zero_cells(Cells * ce);

////////////////////////////////////////////////////////////////////
//...динамическое распределение абстрактной (геометрической) ячейки;
inline void cells_new (Cells * ce, int N_ce, int N_graph, int N_mp) 
{
	if (ce) {
		zero_cells(ce);

		ce->ce    = (Cells **)new_struct(N_ce*sizeof(Cells *));
		ce->graph = (Topo   *)new_struct(N_graph*sizeof(Topo));
		ce->mp    = (CMap   *)new_struct(N_mp*sizeof(CMap));
		ce->pm    = NULL;
		if (ce->mp) ce->mp[N_mp-1] = (CMap)NULL_CELL;
		if (! ce->ce    && N_ce    > 0 ||
			 ! ce->graph && N_graph > 0 ||
			 ! ce->mp    && N_mp    > 0) zero_cells(ce);
		if (ce->ce) 
			 ce->ce[N_ce-1] = ce;
	}
	return;
}

/////////////////////////////////////////////////////////////
//...динамическое распределение абстрактной поверхности тела;
inline Bar * get_bar(int N_ce)
{
	Cells *  bar; cells_new(bar = new_cells(), N_ce+1, 2, 0);
	if (bar) bar->graph[0] = N_ce;
	return(Bar *)bar;
}

//////////////////////////////////////////////////////////
//...уничтожение абстрактной ячейки и поверхности целиком;
inline void delete_cells(Cells *& ce, int id_zero = OK_STATE)
{
	if (id_zero == OK_STATE) zero_cells(ce); 
	delete_struct(ce);
	return;
}

inline void delete_bar (Bar *& bar)
{
	Cells * ce; delete_cells(ce = (Cells *)bar); 
	bar = NULL;	return;
}

///////////////////////////////////////////////
//...уничтожение составной поверхности целиком;
void delete_bar(int N, Bar ** & bar);

//////////////////////////////////
//...размерность поверхности тела;
int bar_dim(Bar * bar);

///////////////////////////////////////////////////////////////////////////////////////////
//...зачитывание геометрической карты и топологии ячейки в кодах ASCII (рабочие программы);
CMap *  map_in(char * id_MAP,  unsigned long & count, unsigned long upper_limit);
Topo * topo_in(char * id_GRAF, unsigned long & count, unsigned long upper_limit, int N = 0);

void map_out (FILE * id_MAP,  CMap * mp, int id_long = NULL_STATE);
void topo_out(FILE * id_GRAF, Topo * graph);
void bar_out (char * ch_BAR,  Bar  * bar);

///////////////////////////////////////////////////////////////////
//...добавление описания заданной ячейки в пространстве параметров;
inline void add_pm(Cells * ce, int k, CMap *& mp)
{
	if (ce && ce->graph && 0 <= k && k < ce->graph[1] && ce->pm) {
		ce->pm[k] = mp; mp = NULL;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////
//...изометрические преобразования над точкой (два поворота вокруг Z и один - вокруг Y);
template <typename T>
inline void point_rotY(T * P, double C, double S)
{
	T		 X = P[0]*C+P[2]*S, Z = P[2]*C-P[0]*S;
	P[0] = X;  P[2] = Z;
}
template <typename T>
inline void point_rotZ(T * P, double C, double S)
{
   T		 X = P[0]*C-P[1]*S, Y = P[1]*C+P[0]*S;
   P[0] = X;  P[1] = Y;
}
template <typename T>
inline void point_shift(T * P, T * P_shift, int id_reverse = NULL_STATE)
{
	if (id_reverse) {
		P[0] -= P_shift[0]; 
		P[1] -= P_shift[1]; 
		P[2] -= P_shift[2];
	}
	else {
		P[0] += P_shift[0]; 
		P[1] += P_shift[1]; 
		P[2] += P_shift[2];
	}
}
template <typename T>
inline void point_iso(T * P, T * P_shift, double CZ, double SZ,
							 double CY, double SY,  double CX, double SX)
{
   point_rotZ(P, CX, SX); //...Euler angles in inverse order !!!
   point_rotY(P, CY, SY);
   point_rotZ(P, CZ, SZ);
   if (P && P_shift) point_shift(P, P_shift);
}

////////////////////////////////////////////////////////////////////////////
//...относительное перемещение и поворот точки на углы, заданные в pадианах;
inline void point_iso(double * P, double * P0, double fi = 0., double theta = 0., double fX = 0.)
{
	if (P) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		point_iso(P, P0, CZ, SZ, CY, SY, CX, SX);
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////
//...изометрические преобразования над геометрической картой, ячейкой и поверхностью;
void  map_normalizat      (CMap * mp, double * P, double * ort);
void  map_normalizat_facet(CMap * mp, double * P, double * ort);
void  map_iso  (CMap  * mp, double * P, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
void  map_iso  (CMap  * mp, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  cells_iso(Cells * ce, double * P, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
void  cells_iso(Cells * ce, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  bar_iso  (Bar  * bar, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  cells_to_origine(Cells * ce);

///////////////////////////////////////////////////////////////////////////////////////////////////
//...перевод геометpических элементов в локальную систему кооpдинат и наобоpот
void  make_local  (double * P, CMap * mp, int id_shift = OK_STATE);
void  make_common (double * P, CMap * mp, int id_shift = OK_STATE);
void  map_local   (CMap   * mp,  CMap * ext_mp);
void  map_common  (CMap   * mp,  CMap * ext_mp);
void  cells_local (Cells * ce,  CMap * ext_mp);
void  cells_common(Cells * ce,  CMap * ext_mp);
void  bar_local   (Bar   * bar, CMap * ext_mp);
void  bar_common  (Bar   * bar, CMap * ext_mp);

///////////////////////////////////////////////////////////////////////////////////////////
//...длина дуги одномерной пространственной ячейки и коррекция локальной системы координат;
inline int get_num(Topo * graph, int m)
{
	if (graph) {
		if (m >= graph[1]) m = graph[1]-1; else
		if (m <  0)        m = 0;
		if (m >= 0) return graph[2+m];
	}
	return(-1/*0 -- до 11.03.2010*/);
}
double cells_length(Cells * ce);
void   edge_correct(Cells * ce);

////////////////////////////////////////////////////////////
//...порождение геометрической карты с нулевыми параметрами;
inline CMap * get_map(int dim, int genus, int id_cell = NULL_CELL)
{
	int     k = size_of_map(dim, genus), id_dop = size_of_dop(id_cell);
	CMap * mp = k && id_dop >= 0 ? (CMap *)new_struct((k+1+id_dop)*sizeof(CMap)) : NULL;
	if (mp) {
		mp[0] = ID_MAP(dim, genus);
		mp[k] = (CMap)id_cell;
	}
	return mp;
}

inline int add_new_maps(CMap *& mp, int N_map, int N_new_maps)
{
	CMap * new_map = N_new_maps > 0 ? (CMap *)new_struct((N_map+N_new_maps)*sizeof(CMap)) : NULL;
	if (new_map) {
		if (mp) memcpy(new_map, mp, N_map*sizeof(CMap));
		delete_struct(mp); mp = new_map; return(1);
	}
	return(0);
}

inline int add_new_cell_maps(CMap *& mp, int id_CELL)
{
	int m = size_of_map(mp);
	CMap * new_map = (CMap *)new_struct((m+1+size_of_dop(id_CELL, mp))*sizeof(CMap));
	if (new_map) {
		memcpy(new_map, mp, m*sizeof(CMap)); new_map[m] = id_CELL;
		delete_struct(mp); mp = new_map; 
		return(1);
	}
	return(0);
}

//////////////////////////////////////////////////////////////
//...идентификация геометрической карты, ячейки и поверхности;
void map_to(CMap *& mp, int dim, int genus);
inline void cells_to(Cells * ce, int genus) { if (ce) map_to(ce->mp, map_dim(ce->mp), genus);}
int  map_id    (CMap  * mp, CMap   * ext_mp);
int  map_into  (CMap  * mp, double X, double Y, double Z, int id_special = OK_STATE);
int  map1into  (CMap  * mp, double X, double Y, double Z);
int  cells_id  (Cells * ce, Cells * ext_ce, int id_num_correct = NULL_STATE);
int  common_id (Cells * ce, Cells * ext_ce);
int  cells_in_struct(Cells * ce, Cells * ext_ce, int id_num_correct);
int  cells_in_bar	  (Bar   * bar, Cells * ce);
void search_cells_element(Cells * ce, int *& id, int & N_buf, int & buf_size, int buf_delta = 20);

//////////////////////////////////////////////
//...добавление нового элемента в граф связей;
int   topo_insert	(Topo *& graph, Topo new_element);
void  topo_exc		(Topo  * graph, Topo some_element);
void	topo_correct(Cells * ce, int N = -1, int N_new = -1, int id_list = OK_STATE);

////////////////////////////////////////////////////////
//...идентификация топологии и упоpядочивание элементов;
void  topo_ord (Topo  * graph);
int   topo_pad (Topo *& graph, int N_pad);
int   topo_id  (Topo  * graph, Topo some_element);
int   topo_id  (Cells * ce,    Topo some_element, Topo exc_element = -1);
void  bar_mutat  (Bar * bar,  int k, int m);
void  cells_mutat(Cells * ce, int k, int m);
void  bar_ord    (Bar   * bar);
void  cells_ord  (Cells * ce);
void  cells_ord0 (Cells * ce, int id_ord = 0);
void  bar_invers (Bar   * bar);

//////////////////////////////////////////////////////////////////////////////////////////////
//...изменение ориентации и циклическая сдвижка элементов границы общей геометрической ячейки;
void cells_invers      (Cells * ce, int id_all = NULL_STATE);
void cells_cyclic_shift(Cells * ce, int      m = 1);

/////////////////////////////////////////////////////////////////////////////
//...извлечение геометрической карты заданной ячейки из описания поверхности;
CMap * map_cpy	 (CMap  * mp, int id_dop = OK_STATE);
CMap ** pm_cpy	 (CMap ** pm, int l);
CMap * map_cpy	 (Bar * bar, int N);
CMap ** pm_cpy	 (Bar * bar, int N);
Topo * graph_cpy(Bar * bar, int N);

////////////////////////////////////////////////////////
//...композиция поверхности тела из отдельных элементов;
void delete_cells_in_struct(Cells * ce);
void search_cells_in_struct(Cells * ce, Cells * ext_ce, int & i, int id_cell, Cells **& dop_ce, int & N_buf, int & buf_size, int buf_delta = 20);
int     bar_add (Bar * bar, Cells *& ext_ce, int id_cell, int buf_delta);
int     bar_add (Bar * bar, Cells *& ext_ce, int id_cell = OK_STATE);
int     bar_span(Bar * bar, CMap  *& mp);
void    trim_add(Cells * ce, Cells *& trim, int id_mp = OK_STATE);
void    bar_exc (Bar * bar, int N);
Cells * bar_cpy (Bar * bar, int N, int id_origine = NULL_STATE);
Cells * bar_sub (Bar * bar, int N, int id_origine = NULL_STATE);

///////////////////////////////////////////////////////////////////////////////////////////////////
//...удаление вырожденных элементов единичной размерности и исключение двойной нумерации точки
void  install_struct				(Bar * bar);
void  bar_generate            (Bar * bar,  double eps = EE_ker);
void  bar_correct_double_point(Bar * bar,  double eps = EE_ker);
void cell_correct_double_point(Cells * ce, double eps = EE_ker);

////////////////////////////////////////////
//...операции над точками плоских элементов;
complex get_arc_center(Cells * ce);
complex get_arc_point (Cells * ce, int m);
complex get_arc_point (Cells * ce, int m, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
double  arc_delta		 (complex z0, complex z1, complex z2, int topo);
double  arc_length	 (complex z0, complex z1, complex z2, int topo);
void set_arc_point    (Cells * ce, int m, complex z);
void set_arc_center   (Cells * ce, complex z);

//////////////////////////////////////////////////////////////////////////////////////
//...коррекция описания прямой или окружности (установка локальной системы координат);
void line_correct(Cells * ce, int id_fast= NULL_STATE);
void circ_correct(Cells * ce, double R, int topo, int points_inverse = NULL_STATE);

//////////////////////////////////////////////////////
//...описание общей части геометрической карты ячейки;
//   mp[0]               - идентификация карты;
//   mp[1], mp[2], mp[3] - координаты опорной точки;
//   mp[4], mp[5], mp[6] - Эйлеровы углы локальной системы координат;

////////////////////////////////////////////////////////////////////////////
//...описание дополнительной части (параметров) геометрической карты ячеек;
//   ID_MAP(1,       NULL_GENUS): отсутствуют;
//
//   ID_MAP(2,       NULL_GENUS): отсутствуют;
//
//   ID_MAP(3,       NULL_GENUS): отсутствуют;
//
//   ID_MAP(1,      NURBS_GENUS): mp[7] - размерность контрольных точек + периодичность кривой;
//                                mp[8] - число узлов;
//                                mp[9] - число контрольных точек;
//
//   ID_MAP(2,      NURBS_GENUS): mp[7]  - размерность контр. точек + периодичность в U- и V-направлениях;
//                                mp[8]  - число узлов в U-направлении;
//                                mp[9]  - число узлов в V-направлении;
//                                mp[10] - число контрольных точек в U-направлении;
//                                mp[11] - число контрольных точек в V-направлении;
//
//   ID_MAP(2, NURBS_NULL_GENUS): mp[7]  - размерность контрольных точек + периодичность кривой;
//                                mp[8]  - число узлов;
//                                mp[9]  - число контрольных точек;
//
//   ID_MAP(2,  NURBS_CYL_GENUS): mp[7]  - размерность контрольных точек + периодичность кривой;
//                                mp[8]  - число узлов;
//                                mp[9]  - число контрольных точек;
//
//   ID_MAP(1,     SPHERE_GENUS): mp[7] - радиус дуги окружности R >= 0;
//                                mp[8] - полудлина дуги окружности в радианах;
//
//   ID_MAP(2,     SPHERE_GENUS): mp[7] - радиус сферы(R < 0 - внутренняя нормаль);
//
//   ID_MAP(2,   SPHEROID_GENUS): mp[7] - большая полуось сфероида;
//                                mp[8] - малая полуось сфероида (м.б. > большой полуоси);
//
//   ID_MAP(2,  ELLIPSOID_GENUS): mp[7] -
//                                mp[8] -
//                                mp[9] - полуоси сфероида;
//
//   ID_MAP(1,        CYL_GENUS): mp[7] - большая полуось эллипса;
//                                mp[8] - малая полуось эллипса;
//
//   ID_MAP(2,        CYL_GENUS): mp[7] - радиус цилиндра (R < 0 - внутренняя нормаль);
//                                mp[8] -
//                                mp[9] - начало и конец образующей цилиндра (если есть);
//
//   ID_MAP(2,       CONE_GENUS): mp[7] - полураствор угла конуса в радианах  (mp[7] < 0 - внутренняя нормаль);
//
//   ID_MAP(2,      TORUS_GENUS): mp[7] - большой радиус R тороидальных координат (R < 0 - внутренняя нормаль);
//                                mp[8] - малый радиус r
//                                (параметр alpha > 0 тороидальных координат z+ir = i th((alpha+ibeta)/2 -- отсутствует);

///////////////////////////////////////
//...описание геометрических элементов;
//       NULL_GENUS : 0 - точка; 1 - прямая;     2 - плоскость; 3 - пространство;
//     SPHERE_GENUS :            1 - окружность; 2 - сфера;
//        CYL_GENUS :            1 - эллипс;     2 - цилиндр;
//    ELL_CYL_GENUS :                            2 - эллиптический цилиндр;
//       CONE_GENUS :            1 - гипербола;  2 - конус;
//   ELL_CONE_GENUS :            1 - парабола;   2 - эллиптический конус;
//      TORUS_GENUS :                            2 - тор;
//   SPHEROID_GENUS :                            2 - сфероид;
//  ELLIPSOID_GENUS :                            2 - эллипсоид;
//      PARAB_GENUS :                            2 - параболоид;
//  ELL_PARAB_GENUS :                            2 - эллиптический параболоид;
//     BI_HYP_GENUS :                            2 - двуполостной гиперболоид;
// ELL_BI_HYP_GENUS :                            2 - эллиптический двуполостной гиперболоид;
//     HYPERB_GENUS :                            2 - гиперболоид;
// ELL_HYPERB_GENUS :                            2 - эллиптический гиперболоид;
//      NURBS_GENUS :            1 - B-кривая;   2 - B-поверхность;
//  NURBS_CYL_GENUS :                            2 - B-осесимметричная поверхность;
// NURBS_NULL_GENUS :                            2 - B-цилиндрическая поверхность;
#endif
