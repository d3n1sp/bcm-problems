/*=========================================*/
/*                  UTILS                  */
/*=========================================*/
#ifndef ___UTILS___
#define ___UTILS___

///////////////////////////////////////////////////
//...standart definitions and including for SYSTEM;
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <iostream>
#include <fstream>
#include "complex.h"
#include "dd_real.h"
#include "qd_real.h"

///////////////////////////////////////////
//...identification of the program version;
const unsigned long ID_VERSION  = 20121108;
inline int Is_Consistent(unsigned long version) { return version <= ID_VERSION; }

////////////////////////////
//...creating any structure;
inline char * new_struct(int size_of_char)
{
  if (size_of_char <= 0) return(NULL);
  char * m = new char[size_of_char];
  memset(m, 0, size_of_char);
  return m;
}
///////////////////////////////////////////
//...template definitions (min/max, swap and delete_struct);
template <typename T>
inline T min(T A, T B) { return A < B ? A : B; }
template <typename T>
inline T max(T A, T B) { return A > B ? A : B; }
template <typename T>
inline void swap(T & A, T & B) { 
	T swap_temp = A; 
	A = B; 
	B = swap_temp; 
}
template <typename T>
inline void delete_struct(T *& m) { delete[] m; m = NULL; }

///////////////////////////////////////////////////////////
//...инициализация и обнуление матрицы (рабочие программы);
template <typename T> 
void set_matrix(T **& m, int dim_N, int dim_M, int dim_buf = 0) {
	delete_struct(m);
	if (dim_N > 0 && dim_M > 0 &&
		(m = (T **)new char [dim_N*(sizeof(T *)+dim_M*sizeof(T))+sizeof(T *)+dim_buf*sizeof(T)]) != NULL) {
		for (int k = 0; k < dim_N; k++) 
		m[k] = (T *)(m+dim_N+1)+k*dim_M; 
		m[dim_N] = (T *)(m+dim_N+1)+dim_N*dim_M+dim_buf;
		memset(m[0], 0, (dim_N*dim_M+dim_buf)*sizeof(T));
	}
}

//////////////////////////////////////
//...обнуление всех элементов матрицы;
template <typename T>
void clean_matrix(T ** m, int dim_N, int dim_M)
{
	if (m) for (int k = 0; k < dim_N; k++)
		memset(m[k], 0,  dim_M*sizeof(T));
}

//////////////////////////
//...standart definitions;
#define  MAX_CMBL        10 //...common block size (for data transfer);
#define  MAX_TIME        20 //...size of array for process timing;
#define  MAX_ITER       100
#define  BUF_SIZE      2000
#undef	MAX_HIT
#define	MAX_HIT        1.111111e38
#undef	MIN_HIT
#define	MIN_HIT       -1.111111e38
#undef	NOT_HIT
#define	NOT_HIT        1.70141e38f

//////////////////////////////////////////
//...константы, определяющие маску ошибки;
#define ERR_NULL         0x00000000
#define ERR_MASK         0xFFFF0000
#define ERR_CODE         0x0000FFFF

////////////////////////////////
//...несколько уровней точности;
extern double EE_fine;
extern double EE;
extern double EE_ker;
extern double EE_dop;
extern double EE_usl;
extern double EE_fuz;
extern double EE_res;

/////////////////////////////////////////////////////////////////
//...общие блоки для таймирования процесса и для передачи данных;
extern clock_t inter_time[MAX_TIME];
extern int         i_parm[MAX_CMBL];
extern double      d_parm[MAX_CMBL];
extern complex     c_parm[MAX_CMBL];

///////////////////////////////////////////
//...description the states of the process;
enum Num_State {
         END_STATE = -777,
         SRF_STATE = -20,
		BOUND1_STATE = -12,
		BOUND2_STATE = -11,
		BOUND3_STATE = -10,
		BOUND4_STATE = -9,
		BOUND5_STATE = -8,
		BOUND6_STATE = -7,
		BOUND7_STATE = -6,
		BOUND8_STATE = -5,
		BOUND9_STATE = -4,
	   BOUND10STATE = -3,
	   BOUND11STATE = -2,
         ERR_STATE = -1,
        NULL_STATE = 0,
		    OK_STATE = 1,
		    TR_STATE = 2,
			 TL_STATE = 3,
	  SCTTYPE_STATE = 12,
		  ZERO_STATE = 1000,
     INITIAL_STATE,
  ADDITIONAL_STATE,
     SPECIAL_STATE,
     SPECIAL2STATE,
     SPECIAL3STATE,
     SPECIAL4STATE,
      EXTERN_STATE,
		ANALYT_STATE,
       CLEAN_STATE,
       BREAK_STATE,
          NO_STATE,
		   GV1_STATE,
		   GV2_STATE,
		   GV3_STATE,
		   GV4_STATE,
		 LOOP1_STATE,
		 LOOP3_STATE,
	  PZS_PZS_STATE = 2000,
	   PZS_NZ_STATE,
	    NZ_NZ_STATE,
	   NZ_PZS_STATE,
};

///////////////////////////////////////////////////////////////////////
//...фиксирование времени работы процесса и его перевод в 60-ю систему;
inline void timing_process(clock_t start, clock_t inter, int * hund, int * sec, int * min, int * hour)
{
	if ((hund || sec || min || hour) && inter >= start) {
		double        f     = (double)(inter-start)/CLOCKS_PER_SEC+.005;
		unsigned long sec0  = (unsigned long)f,
						  hund0 = (unsigned long)(100.*(f-sec0));
		if (hour) {
		  * hour = (int)sec0/3600; sec0 -= 3600*(* hour);
		}
		if (min ) {
		  * min  = (int)sec0/60;   sec0 -= 60*(* min);
		}
		if (sec ) * sec  = (int)sec0; else hund0 += 100*sec0;
		if (hund) * hund = (int)hund0;
	}
	return;
}
inline clock_t timing_process(clock_t start = 0, int * hund = NULL, int * sec  = NULL,
														int * min  = NULL, int * hour = NULL)
{
	clock_t inter = clock();
	timing_process(start, inter, hund, sec, min, hour);
	return(inter);
}

//////////////////////////////
//...pаспутывание комментаpия;
char user_Filtr(char * FILE, unsigned long & k, unsigned long upper_limit, Num_State id_toupper = OK_STATE);
int  user_Read (char * buf,  char * FILE, unsigned long & k, unsigned long upper_limit, Num_State id_toupper = OK_STATE);
int  user_Read (char * buf,  char * FILE, unsigned long & k, unsigned long upper_limit, int m, Num_State id_toupper = OK_STATE);
int  user_Count(char * FILE, unsigned long k, unsigned long & upper_limit, char punct = '\xFF');

////////////////////////////////////////
//...чтение с учетом десятичной запятой;
inline double user_strtod(char * buf)
{
	char * pchar = buf, * pnext;
	while ( (pnext = strstr(pchar, ".")) != NULL) {
			pnext[0] = localeconv()->decimal_point[0];
			pchar = pnext+1;
	}
	return(strtod(buf, NULL));
}

/////////////////////////////////////////////
//...длина строки с учетом нулеого указателя;
inline size_t user_strlen(		  char * buf) { return(buf ? strlen(buf) : 0); }
inline size_t user_strlen(const char * buf) { return strlen(buf); }

////////////////////////////////////////////////////
//...проверки и считывание данных из входного файла;
unsigned long length_ascii(		char * cfg_name);
unsigned long length_ascii(const char * cfg_name);
char   * read_struct_ascii(		char * cfg_name);
char   * read_struct_ascii(const char * cfg_name);
int				exist_ascii(		char * cfg_name);
int				exist_ascii(const char * cfg_name);

///////////////////////////////////////
//...упаковка и распаковка целых чисел;
inline void pack_ints(int m1, int m2, double & A)
{
	union { double A; int m[2]; } temp;
	temp.m[0] = m1;
	temp.m[1] = m2;
	A = temp.A;
}

inline void unpack_ints(double A, int & m1, int & m2)
{
	union { double A; int m[2]; } temp;
	temp.A = A;
	m1 = temp.m[0];
	m2 = temp.m[1];
}

inline double PackInts(int m1, int m2)
{
	double A; pack_ints(m1, m2, A);
	return A;
}

inline int UnPackInts(double A, int m = 0)
{
	int m1, m2; unpack_ints(A, m1, m2);
	return (m ? m2 : m1);
}

///////////////////////////////////////
//...преобразования вещественного типа;
inline double to_double(const double &a) { return a;}
inline int to_int (const double &a) { return static_cast<int>(a);}
inline double real(const double &a) { return a;}
inline double imag(const double &a) { return 0.0;}

/////////////////////////////
//...фильтры для малых чисел;
template <typename T>
inline T filtr_(T A, double eps = 0.0) { return (fabs(A) < eps ? T(0) : A);}
template <typename T>
inline double filtr_Re(T A, double eps = 0.0) { return real(filtr_(A, eps));}
template <typename T>
inline double filtr_Im(T A, double eps = 0.0) { return imag(filtr_(A, eps));}

////////////////////////////////////////////////////////////////////////////////////
//...заполнение кооpдинат и определение расстояния между точками (pабочие пpогpаммы);
inline void   set_point(double * P, double * ext_P, double f = 1.) { P[0] = ext_P[0]*f; P[1] = ext_P[1]*f; P[2] = ext_P[2]*f;}
inline double abs_point(double * P, double * ext_P) { return sqrt(sqr(P[0]-ext_P[0])+sqr(P[1]-ext_P[1])+sqr(P[2]-ext_P[2]));}
inline double max_point(double * P, double * ext_P) { return max(fabs(P[0]-ext_P[0]), max(fabs(P[1]-ext_P[1]), fabs(P[2]-ext_P[2])));}
inline void   sup_point(complex z, double & X1, double & X2, double & Y1, double & Y2) {
  X1 = min(X1, real(z)); X2 = max(X2, real(z));
  Y1 = min(Y1, imag(z)); Y2 = max(Y2, imag(z));
}

/////////////////////////////////////////////////////////////////
//...линейная аппроксимация табличных данных (рабочая программа);
double table_approx(double arg, double table[][2], int N_table, int decrease_flag = NULL_STATE);

//////////////////////////////////////////////////////////////////////////////////////
//...пpоцедуpы целочисленного возведения в степень вещественного и комплексного числа;
inline double powI(double x, int m)
{
  double f = 1.;
  for (int l = 1; l <= abs(m); l++) f *= x;
  return f;
}

inline complex powC(complex x, int m)
{
  complex f = comp(1);
  for (int l = 1; l <= abs(m); l++) f *= x;
  return f;
}

///////////////////////////////////////////////////////////////
//...корректное вычисление значения аpгумента: arg0([0, 2*pi));
inline double arg0(double y, double x)
{
	double fi = 0.;
	if (fabs(x) > fabs(y)) fi = atan(fabs(y/x)); else if (y != 0.) fi = M_PI_2-atan(fabs(x/y));
	if (x < 0.) if (y <= 0.) fi = M_PI+fi; else fi = M_PI-fi; else if (y < 0.) fi = 2.*M_PI-fi;
	return(fi);
}
inline double arg0(complex z)
{
	return(arg0(imag(z), real(z)));
}

///////////////////////////////////////////////////////////////
//...корректное вычисление значения аpгумента: arg2([-pi, pi));
inline double arg2(double y, double x)
{
	double fi = 0.;
	if (fabs(x) > fabs(y)) fi = atan(fabs(y/x)); else if (y != 0.) fi = M_PI_2-atan(fabs(x/y));
	if (x < 0.) if (y <= 0.) fi = -M_PI+fi; else fi = M_PI-fi; else if (y <= 0.) fi = -fi;
	return fi;
}
inline double arg2(complex z)
{
	return (arg2(imag(z), real(z)));
}

/////////////////////////////////////////////////////////////////////////////
//...вспомогательные функции для работы с файлом геометрии массива элементов;
  int  geom_pad    (int *& geom, int k, int N_pad);
  void geom_pad1   (int *& geom,        int N_pad, int i, int & m);
  void geom_pad2   (int *& geom,        int N_pad, int i, int & m, int Ng_buf);
  int  geom_unpad  (int *& geom, int k, int N_pad);
  int  geom_insert (int *& geom, int k, int new_element, int start = 1);
  void geom_insert1(int *& geom,        int new_element, int start, int i, int & m);
  void geom_insert2(int *& geom,        int new_element, int start, int i, int & m, int N_geom);
  int  geom_indent (int *& geom, int k, int new_element, int start = 1);
  int  geom_add    (int *& geom, int k, int type, int N);
  void geom_add1   (int *& geom,        int type, int N, int i, int & m);
  int  geom_exl    (int *& geom, int k);
  int  graph_insert(int *& graph, int start, int element);

//////////////////////////////////////////////////////
//...идентификация элементов в топологическом формате;
int geom_link_id4 (int m1, int m2, int m3, int m4, int i, int * geom, int id_pos = 2);
int geom_link_id3 (int m1, int m2, int m3,         int i, int * geom, int id_pos = 2);
int geom_link_id1 (int m1,                         int i, int * geom, int id_pos = 2);
int geom_link_id4s(int m1, int m2, int m3, int m4, int i, int * geom, int id_pos = 2);
int geom_link_id3s(int m1, int m2, int m3,         int i, int * geom, int id_pos = 2);
int geom_link_id2s(int m1, int m2,                 int i, int * geom, int id_pos = 2);
int geom_link_id1s(int m1,									int i, int * geom, int id_pos = 2);
int geom_link_id4 (int m1, int m2, int m3, int m4, int k1, int k2);
int geom_link_id3 (int m1, int m2, int m3,         int k1, int k2);
int geom_link_id3 (int m1, int m2, int m3,         int k1);
int geom_link_id2 (int m1, int m2,						int k1);

///////////////////////////////////
//...алгоритмы определения соседей;
void ListBlocks3D_utils (int _nelem, double *_xyzrelem, int _nballs, double *_xyzrballs, // Compute the list of blocks to be checked
								 int *&_ilinks, int *&_jlinks);

////////////////////////////////////////////////////////////////////////////////////////////////
//...алгоритмы перенумерации симметричной, разреженной матрицы (reverse Cuthill-Mckee ordering);
void genrcm_utils(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask);
void subrcm_utils(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n);
/////////////////////////
//...описание параметров;
//   n          - the dimension of the matrix;
//   xadj, iadj - the matrix structure: xadj[n+1], iadj[*]; information about row i is stored 
//                in xadj[i-1] -- xadj[i]-1 of the the adjacency structure iadj[*];
//                for each row, it contains the column indices of the nonzero entries;
//   perm[n]    - contains the rcm ordering;
//   mask[n]    - marks variables that have been numbered (working array);
//   xls[n+1]   - the index vector for a level structure; the level structure is stored 
//                in the currently unused spaces in the permutation vector perm;
//   nsubg      - the size of the subgraph;
//   subg[n]    - contains the nodes in subgraph (which may be disconnected);
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//...алгоритм вычисления с.з. и с.ф. общей матрицы QR-методом;
int HQRfunction(double  ** A, double  * H_re, double  * H_im, int dim_N, int Max_iter = 30);
int HQRfunction(dd_real ** A, dd_real * H_re, dd_real * H_im, int dim_N, int Max_iter = 30);
int HQRfunction(qd_real ** A, qd_real * H_re, qd_real * H_im, int dim_N, int Max_iter = 30);
int HQRfunction(complex ** A, complex * H_re, complex * H_im, int dim_N, int Max_iter = 30);
#endif
