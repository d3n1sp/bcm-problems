#include "stdafx.h"
#include "cgrid.h"

////////////////////////////////
//...несколько уровней точности;
double EE_fine= 1e-21; //...fine accuracy;
double EE     = 1e-15; //...refine accuracy;
double EE_ker = 1e-10; //...accuracy of geometrical kernel operations;
double EE_dop = 1e-8;  //...middle accuracy for geometrical checkings;
double EE_usl = 1e-5;  //...middle accuracy for geometrical checkings;
double EE_fuz = 1e-3;  //...accuracy for chain closing;
double EE_res = 2e-3;  //...for logarithmic brand under resonance;

////////////////////////////////////////////////////////////////
//...общие блоки (для таймирования процесса и  передачи данных);
clock_t inter_time[MAX_TIME];
int         i_parm[MAX_CMBL];
double      d_parm[MAX_CMBL];
complex     c_parm[MAX_CMBL];

//////////////////////////////
//...pаспутывание комментаpия;
char user_Filtr(char * FILE, unsigned long & k, unsigned long upper_limit, Num_State id_toupper)
{
	if (k >= upper_limit || ! FILE) return('\xFF');
	char buf = id_toupper == OK_STATE ? (char)toupper(FILE[k++]) : (char)(FILE[k++]);
	if (buf == ';') {
		while (buf != '\n' && k < upper_limit) buf = FILE[k++];
		return('\xFF');
	}
	if (isalpha(buf) || isdigit(buf) || buf == '+'  || buf == ':' && (k == 1 || isalpha(FILE[k-2])) || 
			buf == '('  || buf == ')'  || buf == '-'  || buf == '.' || buf == '#' || buf == '_' || 
													buf == '\x1C' || buf == '!' || buf == '\x2F') return(buf);
	else return('\xFF');
}

int user_Read(char * buf, char * FILE, unsigned long & k, unsigned long upper_limit, Num_State id_toupper)
{
	int l = 0;
	while (l < BUF_SIZE && k < upper_limit) {
		char buf0 = user_Filtr(FILE, k, upper_limit, id_toupper);
		if (buf0 != '\xFF') buf[l++] = buf0; else
		if (l != 0) {
			buf[l] = '\x0'; return(1);
		}
	}
	buf[l] = '\x0'; return (l != 0);
}

int user_Read(char * buf, char * FILE, unsigned long & k, unsigned long upper_limit, int m, Num_State id_toupper)
{
	int l = 0;
	for (m--; m >= 0; m--) if (l < BUF_SIZE && k < upper_limit)
	buf [l++]  = id_toupper == OK_STATE ? (char)toupper(FILE[k++]) : (char)(FILE[k++]);
	buf [l]    = '\x0';
	return (l != 0);
}

int user_Count(char * FILE, unsigned long k, unsigned long & upper_limit, char punct)
{
	char buf = FILE[k];
	while (punct == '\xFF' && (isalpha(buf) || buf == '#') ||
			punct != '\xFF' && buf == punct  && buf != '\x0') buf = FILE[++k];

	upper_limit = k;
	if (buf == '\x0') return(0); else
	while (! (punct == '\xFF' && (isalpha(buf) || buf == '#') &&
			(k == 0 || toupper(buf) != 'E' || ! isdigit(FILE[k-1]) && FILE[k-1] != '.') ||
			punct != '\xFF' && buf == punct  || buf == '\x0')) {
		if (buf == ';') while(buf != '\n' && buf != '\x0') buf = FILE[++k];
		if (buf != '\x0') buf = FILE[++k];
	}
	upper_limit = k; return(1);
}

////////////////////////////////////////////
//...определение длины файла входных данных;
unsigned long length_ascii(FILE * TXT)
{
	unsigned long length = 0;
	if (TXT) {
		fseek (TXT, 0L, SEEK_END);
		length = (unsigned long) ftell(TXT);
		fclose(TXT); 
	}  
	return(length);
}
unsigned long length_ascii(char * cfg_name)
{
	FILE * TXT = fopen(cfg_name, "rt");
	return(length_ascii(TXT));
}
unsigned long length_ascii(const char * cfg_name)
{
	FILE * TXT = fopen(cfg_name, "rt");
	return(length_ascii(TXT));
}

/////////////////////////////////////
//...считывание файла входных данных;
void read_struct_ascii(FILE * device, char * ascii, unsigned long length)
{
	if (device) {
		size_t res = fread(ascii, sizeof(char), length, device);
		fclose(device);
	}
}
char * read_struct_ascii(char * cfg_name)
{
	unsigned long length = length_ascii(cfg_name);
	char *         ascii = NULL;
	if (length && (ascii = (char *)new_struct((length+1)*sizeof(char))) != NULL) {
		FILE * device = fopen(cfg_name, "rb");
		read_struct_ascii(device, ascii, length);
	}
	return(ascii);
}
char * read_struct_ascii(const char * cfg_name)
{
	unsigned long length = length_ascii(cfg_name);
	char *         ascii = NULL;
	if (length && (ascii = (char *)new_struct((length+1)*sizeof(char))) != NULL) {
		FILE * device = fopen(cfg_name, "rb");
		read_struct_ascii(device, ascii, length);
	}
	return(ascii);
}

////////////////////////////////////////////////////
//...проверка наличия файла входных данных на диске;
int exist_ascii(FILE * device)
{
	int m = device != 0;
	if (m) fclose (device);
	return(m);
}
int exist_ascii(char * cfg_name)
{
	FILE * device = fopen(cfg_name, "rb");
	return(exist_ascii(device));
}
int exist_ascii(const char * cfg_name)
{
	FILE * device = fopen(cfg_name, "rb");
	return(exist_ascii(device));
}

/////////////////////////////////////////////
//...линейная аппроксимация табличных данных;
double table_approx(double arg, double table[][2], int N_table, int decrease_flag)
{
	double value = 0.; int i;
	if (N_table > 0 && table) 
	if (! decrease_flag) {
		for (value = table[i = 0][1]; i < N_table; value = table[i++][1])
			if (arg < table[i][0]) break;
		if (0 < i && i < N_table) value += (arg-table[i-1][0])*(table[i][1]-table[i-1][1])/(table[i][0]-table[i-1][0]);
	}
	else {
		for (value = table[i = N_table-1][1]; i >= 0; value = table[i--][1])
			if (arg > table[i][0]) break;
		if (0 <= i && i < N_table-1) value -= (table[i+1][0]-arg)*(table[i+1][1]-table[i][1])/(table[i+1][0]-table[i][0]);
	}
	return value;
}

/////////////////////////////////////////////////////////////////////////////////
//...увеличение размера описании геометрии массива элементов (рабочая программа);
int geom_pad(int *& geom, int k, int N_pad)
{
	int i = 0, m, l;
	if (geom && geom[0] > 0 && N_pad > 0 &&  k <= geom[0] && k > 0) {
		for (i = 1, l = 0; l < k-1;     l++) i += geom[++i]+1;
		for (m = i;        l < geom[0]; l++) m += geom[++m]+1;

      int * new_geom = (int *)new_struct((m+N_pad)*sizeof(int));
      if (! new_geom) return(0);

      memcpy(new_geom,           geom,       (i+2)*sizeof(int));
      memcpy(new_geom+i+2+N_pad, geom+i+2, (m-i-2)*sizeof(int));

      delete_struct(geom); geom = new_geom;
      geom[i+1] += N_pad;
	}
	return(i);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//...увеличение размера описании геометрии массива элементов по имеющимся данным о расположении раздела;
void geom_pad1(int *& geom, int N_pad, int i, int & m)
{
	int * new_geom = (int *)new_struct((m+N_pad)*sizeof(int));
	if (! new_geom) return;

	memmove(new_geom,           geom,       (i+2)*sizeof(int));
	memmove(new_geom+i+2+N_pad, geom+i+2, (m-i-2)*sizeof(int));

	delete_struct(geom); geom = new_geom;
	geom[i+1] += N_pad;
	m         += N_pad;
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...увеличение размера геометрии массива элементов без дополнительного распределения памяти;
void geom_pad2(int *& geom, int N_pad, int i, int & m, int Ng_buf)
{
	memmove(geom+i+2+N_pad, geom+i+2, (Ng_buf-i-2)*sizeof(int));
	geom[i+1] += N_pad;
	m         += N_pad;
	return;
}

/////////////////////////////////////////////////////////////////////////////////
//...сокращение размера описании геометрии массива элементов (рабочая программа);
int geom_unpad(int *& geom, int k, int N_pad)
{
	int i = 0, m, l;
	if (geom && geom[0] > 0 && N_pad > 0 &&  k <= geom[0] && k > 0) {
      for (i = 1, l = 0; l < k-1;     l++) i += geom[++i]+1;
      for (m = i;        l < geom[0]; l++) m += geom[++m]+1; if (N_pad > geom[i+1]) return(0);

      int * new_geom = (int *)new_struct((m -= N_pad)*sizeof(int));
      if (! new_geom) return(0);

      memcpy(new_geom,     geom,             (i+2)*sizeof(int));
      memcpy(new_geom+i+2, geom+i+2+N_pad, (m-i-2)*sizeof(int));

      delete_struct(geom); geom = new_geom;
      geom[i+1] -= N_pad;
	}
	return(i);
}

//////////////////////////////////////////////////////////////////////////////////
//...добавление нового элемента в геометрию массива элементов (рабочая программа);
int geom_insert(int *& geom, int k, int new_element, int start)
{
	int i = 0, m, l, j;
	if (geom && geom[0] > 0 && k < geom[0] && k >= 0) {
      for (i = 1, l = 0; l < k; l++) i += geom[++i]+1;

////////////////////////////////////////////
//...проверяем наличие элемента в геометрии;
      j = geom[i+1]+1;
      while (start < j && new_element != geom[i+j]) j--;

//////////////////////////////////////////
//...вставляем элемент в начало геометрии;
		if (j == start) {
         for (m = i; l < geom[0]; l++) m += geom[++m]+1;

         int * new_geom = (int *)new_struct((m+1)*sizeof(int));
         if (! new_geom) return(0);

         memcpy(new_geom,           geom,           (i+1+start)*sizeof(int));
         memcpy(new_geom+i+2+start, geom+i+1+start, (m-i-1-start)*sizeof(int));

         delete_struct(geom); geom = new_geom;
         geom[i+1]       += 1;
         geom[i+1+start]  = new_element;
		}
	}
	return(i);
}

///////////////////////////////////////////////////////////////////////////
//...добавление нового элемента по имеющимся данным о расположении раздела;
void geom_insert1(int *& geom, int new_element, int start, int i, int & m)
{
////////////////////////////////////////////
//...проверяем наличие элемента в геометрии;
	int   j = geom[i+1]+1;
	while (start < j && new_element != geom[i+j]) j--;

/////////////////////////////////////////
//...вставляем элемент в начало геометрии;
	if (j == start) {
      int * new_geom = (int *)new_struct((m+1)*sizeof(int));
      if (! new_geom) return;

      memmove(new_geom,           geom,           (i+1+start)*sizeof(int));
      memmove(new_geom+i+2+start, geom+i+1+start, (m-i-1-start)*sizeof(int));

      delete_struct(geom); geom = new_geom;
      geom[i+1]       += 1;
      geom[i+1+start]  = new_element;
      m               += 1;
	}
	return;
}

/////////////////////////////////////////////////////////
//...добавление нового элемента без распределения памяти;
void geom_insert2(int *& geom, int new_element, int start, int i, int & m, int N_geom)
{
////////////////////////////////////////////
//...проверяем наличие элемента в геометрии;
	int   j = geom[i+1]+1;
	while (start < j && new_element != geom[i+j]) j--;

/////////////////////////////////////////
//...вставляем элемент в начало геометрии;
	if (j == start) {
      memmove(geom+i+2+start, geom+i+1+start, (N_geom-i-1-start)*sizeof(int));

      geom[i+1]       += 1;
      geom[i+1+start]  = new_element;
      m               += 1;
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////
//...идентификация элемента в геометрии массива элементов (рабочая программа);
int geom_indent(int *& geom, int k, int new_element, int start)
{
	int i, l, j;
	if (geom && geom[0] > 0 && k < geom[0] && k >= 0) {
		for (i = 1, l = 0; l < k; l++) i += geom[++i]+1;

////////////////////////////////////////////
//...проверяем наличие элемента в геометрии;
      j = geom[i+1]+1;
      while (start < j && new_element != geom[i+j]) j--;

      if (j != start) return(i);
	}
	return(0);
}

///////////////////////////////////////////////////////////////
//...добавление нового раздела в геометрию (рабочая программа);
int geom_add(int *& geom, int k, int type, int N)
{
	int i = 0, m, l;
	if (geom && N >= 0 && 0 <= k && k <= geom[0]) {
		for (i = 1, l = 0; l < k;       l++)
           i += geom[++i]+1;
      for (m = i;        l < geom[0]; l++) 
           m += geom[++m]+1;

      int * new_geom = (int *)new_struct((m+2+N)*sizeof(int));
      if (! new_geom) return(0);

      memcpy(new_geom,       geom,       i*sizeof(int));
      memcpy(new_geom+i+2+N, geom+i, (m-i)*sizeof(int));

      delete_struct(geom); geom = new_geom;
      geom[i]   = type;
      geom[i+1] = N;
      geom[0]  += 1;
	}
	return(i);
}

//////////////////////////////////////////////////////////////////////////////////////
//...добавление нового раздела в геометрию по имеющимся данным о расположении раздела;
void geom_add1(int *& geom, int type, int N, int i, int & m)
{
	int * new_geom = (int *)new_struct((m+2+N)*sizeof(int));
	if (! new_geom) return;

	memmove(new_geom,       geom,       i*sizeof(int));
	memmove(new_geom+i+2+N, geom+i, (m-i)*sizeof(int));

	delete_struct(geom); geom = new_geom;
	geom[i]   = type;
	geom[i+1] = N;
	geom[0]  += 1;
	m        += 2+N;
	return;
}

////////////////////////////////////////////////////////////////////////////////////
//...добавление нового раздела в геометрию без дополнительного распределения памяти;
void geom_add2(int *& geom, int type, int N, int i, int & m)
{
	memmove(geom+i+2+N, geom+i, (m-i)*sizeof(int));

	geom[i]   = type;
	geom[i+1] = N;
	geom[0]  += 1;
	m        += 2+N;
	return;
}

/////////////////////////////////////////////////////////
//...исключение раздела из геометрии (рабочая программа);
int geom_exl(int *& geom, int k)
{
	int i = 0, m, l;
	if (geom && geom[0] > 0 && k > 0 && k <= geom[0]) {
      for (i = 1, l = 0; l < k-1;     l++) i += geom[++i]+1;
      for (m = i;        l < geom[0]; l++) m += geom[++m]+1;

      int * new_geom = (int *)new_struct((m -= 2+geom[i+1])*sizeof(int));
      if (! new_geom) return(0);

      memcpy(new_geom,   geom,                   i*sizeof(int));
      memcpy(new_geom+i, geom+i+2+geom[i+1], (m-i)*sizeof(int));

      delete_struct(geom); geom = new_geom;
      geom[0] -= 1;
	}
	return(i);
}

//////////////////////////////////////////////////////////////////
//...добавление элемента в описание топологии (рабочая программа);
int graph_insert(int *& graph, int start, int element)
{
	if (graph &&  graph[1] > 0 && start <= graph[1]+1 && start > 0) {
      int * new_graph =   (int *)new_struct((graph[1]+3)*sizeof(int));
      if (! new_graph) return(0);

      memcpy (new_graph,         graph,         (start+1)*sizeof(int));
      memcpy (new_graph+start+2, graph+start+1, (graph[1]-start+1)*sizeof(int));
      delete_struct(graph);      graph = new_graph;

      graph[start+1] = element;
      graph[0]      += 1;
      graph[1]      += 1;
	}
	return(0);
}
/*============================================================*/
/*      ИДЕНТИФИКАЦИЯ ЭЛЕМЕНТОВ В ТОПОЛОГИЧЕСКОМ ФОРМАТЕ      */
/*============================================================*/
///////////////////////////////////////////////////////////////////////////////////////
//...идентификация четырехугольной грани в пространственных элементах геометрии ячейки;
int geom_link_id4(int m1, int m2, int m3, int m4, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_BOXS:  
		case GL_PENTA:
		case GL_PYRAMID: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m2) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m3) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m4) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

///////////////////////////////////////////////////////////////////////////////////
//...идентификация треугольной грани в пространственных элементах геометрии ячейки;
int geom_link_id3(int m1, int m2, int m3, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_PENTA:  
		case GL_TETRA: 
		case GL_PYRAMID: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m2) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m3) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

///////////////////////////////////////////////////////////////////////
//...идентификация точки в пространственных элементах геометрии ячейки;
int geom_link_id1(int m1, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_BOXS:  
		case GL_PENTA:  
		case GL_TETRA:
		case GL_PYRAMID: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

////////////////////////////////////////////////////////////////////////////////////
//...идентификация четырехугольной грани в поверхностных элементах геометрии ячейки;
int geom_link_id4s(int m1, int m2, int m3, int m4, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_QUADS: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m2) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m3) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m4) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

////////////////////////////////////////////////////////////////////////////////
//...идентификация треугольной грани в поверхностных элементах геометрии ячейки;
int geom_link_id3s(int m1, int m2, int m3, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_TRIANGLES: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m2) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m3) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

////////////////////////////////////////////////////////////////////
//...идентификация линии в поверхностных элементах геометрии ячейки;
int geom_link_id2s(int m1, int m2, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_TRIANGLES:
		case GL_QUADS: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m2) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

///////////////////////////////////////////////////////////////
//...идентификация точки в линейных элементах геометрии ячейки;
int geom_link_id1s(int m1, int i, int * geom, int id_pos)
{
	int l;
	switch (geom[i]) {
		case GL_LINE_STRIP: {
			for (l = geom[i+1]; l > id_pos; l--) if (geom[i+l+1] == m1) break; if (l == id_pos) return(0);
			return(1);
		}
	}
	return(0);
}

/////////////////////////////////////////////////
//...идентификация точек в четырехугольной грани;
int geom_link_id4(int m1, int m2, int m3, int m4, int k1, int k2)
{
	int  l, mm[] = {m1, m2, m3, m4};
	for (l = 4; l > 0; l--) if (mm[l-1] == k1) break; if (! l) return(0);
	for (l = 4; l > 0; l--) if (mm[l-1] == k2) break; if (! l) return(0);
	return(1);
}

/////////////////////////////////////////////
//...идентификация точек в треугольной грани;
int geom_link_id3(int m1, int m2, int m3, int k1, int k2)
{
	int  l, mm[] = {m1, m2, m3};
	for (l = 3; l > 0; l--) if (mm[l-1] == k1) break; if (! l) return(0);
	for (l = 3; l > 0; l--) if (mm[l-1] == k2) break; if (! l) return(0);
	return(1);
}

/////////////////////////////////////////////
//...идентификация точки в треугольной грани;
int geom_link_id3(int m1, int m2, int m3, int k1)
{
	int  l, mm[] = {m1, m2, m3};
	for (l = 3; l > 0; l--) if (mm[l-1] == k1) break; if (! l) return(0);
	return(1);
}

/////////////////////////////////
//...идентификация точки в линии;
int geom_link_id2(int m1, int m2, int k1)
{
	int  l, mm[] = {m1, m2};
	for (l = 2; l > 0; l--) if (mm[l-1] == k1) break; if (! l) return(0);
	return(1);
}

/*=====================================================================*/
/*                   АЛГОРИТМЫ ОПРЕДЕЛЕНИЯ СОСЕДЕЙ                     */
/*=====================================================================*/
//========================================================================
// Compare two integer values
int compint_utils (const void *arg1, const void *arg2) { // Compare two integer values
	int *iarg1= (int *) arg1;
	int *iarg2= (int *) arg2;
	if (*iarg1 <  *iarg2) {
		return -1;
	} else if (*iarg1 == *iarg2) {
		return 0;
	} else {
		return 1;
	};
};

//========================================================================
// Compute the list of blocks to be checked
void ListBlocks3D_utils (int _nelem, double *_xyzrelem, int _nballs, double *_xyzrballs, // Compute the list of blocks to be checked
				         int *&_ilinks, int *&_jlinks) {

// Compute the size of the integer 3D mesh

	int ncube;

	double daux, daux1;

	daux = pow ((double) _nelem,1.0e0/3.0e0);

	ncube = (int) daux;

	if (ncube < 1) ncube = 1;

// Compute the x, y and z min/max values

	int ix=0, iy=1, iz=2, ir=3;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	double xminl, xmaxl, yminl, ymaxl, zminl, zmaxl;

	int ielem, iball, i, j, k, ii, jj;

	ielem = 0;

	xmin = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
	xmax = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];

	ymin = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
	ymax = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];

	zmin = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
	zmax = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];

	for (ielem=1;ielem<_nelem;ielem++) {
		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		if (daux < xmin) xmin = daux;
		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		if (daux > xmax) xmax = daux;
		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		if (daux < ymin) ymin = daux;
		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		if (daux > ymax) ymax = daux;
		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		if (daux < zmin) zmin = daux;
		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		if (daux > zmax) zmax = daux;
	};

// Set x, y and z linear interpolation data

	double deltax = (xmax-xmin) / ((double) ncube);
	double deltay = (ymax-ymin) / ((double) ncube);
	double deltaz = (zmax-zmin) / ((double) ncube);

// For each subcube count the number of nodes

	int ncube3 = ncube * ncube * ncube;

	int ncube2 = ncube * ncube;

	int *icube, *jcube, *iptr;

	icube = new int [ncube3+1];
	iptr = new int [ncube3];

	for (i=0;i<ncube3+1;i++) icube[i] = 0;

	int ixmin, ixmax, iymin, iymax, izmin, izmax;

	for (ielem=0;ielem<_nelem;ielem++) {

// For current element determine the set of 3D indices that contains the element

		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Register that all detected subcubes contain current element

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					icube[i*ncube2+j*ncube+k+1]++;
				};
			};
		};

	};

	for (i=0;i<ncube3;i++) icube[i+1] += icube[i];

// Complete computation of the lists

	int nlist = icube[ncube3];

	jcube = new int [nlist];

	for (i=0;i<ncube3;i++) iptr[i] = icube[i];

	for (ielem=0;ielem<_nelem;ielem++) {

// For current element again determine the set of 3D indices that contains the element

		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Register that all detected subcubes contain current element

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					int ip = iptr[i*ncube2+j*ncube+k];
					jcube[ip] = ielem;
					iptr[i*ncube2+j*ncube+k]++;
				};
			};
		};

	};

// Prepare mask arrays

	int *imask, *list;

	imask = new int [_nelem];
	list = new int [_nelem];

	int icycle = -1;

	for (i=0;i<_nelem;i++) imask[i] = icycle;

// For each ball compute the set of elements

	_ilinks = new int [_nballs+1];

	int **jlinksarr;

	jlinksarr = new int * [_nballs];

	_ilinks[0] = 0;

	for (iball=0;iball<_nballs;iball++) {

// For current ball determine the set of 3D indices that contains the ball

		xminl = _xyzrballs[iball*4+ix] - _xyzrballs[iball*4+ir];
		daux = xminl;
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		xmaxl = _xyzrballs[iball*4+ix] + _xyzrballs[iball*4+ir];
		daux = xmaxl;
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		yminl = _xyzrballs[iball*4+iy] - _xyzrballs[iball*4+ir];
		daux = yminl;
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		ymaxl = _xyzrballs[iball*4+iy] + _xyzrballs[iball*4+ir];
		daux = ymaxl;
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		zminl = _xyzrballs[iball*4+iz] - _xyzrballs[iball*4+ir];
		daux = zminl;
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		zmaxl = _xyzrballs[iball*4+iz] + _xyzrballs[iball*4+ir];
		daux = zmaxl;
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Scan the list of elements in each subcube

		icycle++;

		nlist = 0;

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					int ind = i*ncube2+j*ncube+k;
					for (ii=icube[ind];ii<icube[ind+1];ii++) {
						jj = jcube[ii];
						if (imask[jj] != icycle) {

							imask[jj] = icycle;

							daux = _xyzrelem[jj*4+ix] - _xyzrelem[jj*4+ir];
							if (daux > xmaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+ix] + _xyzrelem[jj*4+ir];
							if (daux < xminl) goto ExitElem;
							daux = _xyzrelem[jj*4+iy] - _xyzrelem[jj*4+ir];
							if (daux > ymaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+iy] + _xyzrelem[jj*4+ir];
							if (daux < yminl) goto ExitElem;
							daux = _xyzrelem[jj*4+iz] - _xyzrelem[jj*4+ir];
							if (daux > zmaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+iz] + _xyzrelem[jj*4+ir];
							if (daux < zminl) goto ExitElem;

							list[nlist] = jj;
							nlist++;

ExitElem:;

						};
					};
				};
			};
		};

		qsort (list, nlist, sizeof(int), compint_utils);

		int *jlinkloc;

		jlinkloc = new int [nlist];

		for (i=0;i<nlist;i++) jlinkloc[i] = list[i];

		jlinksarr[iball] = jlinkloc;

		_ilinks[iball+1] = _ilinks[iball] + nlist;

	};

// Transform data

	int nz = _ilinks[_nballs];

	_jlinks = new int [nz];

	for (iball=0;iball<_nballs;iball++) {
		int nzloc = _ilinks[iball+1]-_ilinks[iball];
		int ibs = _ilinks[iball];
		int *jlinkloc = jlinksarr[iball];
		for (i=0;i<nzloc;i++) _jlinks[ibs+i] = jlinkloc[i];
	};

// Free work arrays

	delete [] icube;
	delete [] iptr;
	delete [] jcube;
	delete [] imask;
	delete [] list;
	for (i=0;i<_nballs;i++) {
		delete [] jlinksarr[i];
	};
	delete [] jlinksarr;

};

/*=====================================================================*/
/*                     АЛГОРИТМЫ ПЕРЕНУМЕРАЦИИ                         */
/*=====================================================================*/
////////////////////////////////////////////////////////////////////
//...generates the connected level structure rooted at a given node;
void rootls_utils(int root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int n)
{
   int i, iccsze, j, jstop, jstrt, lbegin, lvlend, nbr, node;

   mask[root-1] = 0;
   ls[0]        = root;
   nlvl         = 0;
   lvlend       = 0;
   iccsze       = 1;

   do {
       lbegin      = lvlend + 1;
       lvlend      = iccsze;
       xls[nlvl++] = lbegin;

/////////////////////////////////////////////////////////////////////////////////////////////
//...Generate the next level by finding all the masked neighbors of nodes in the current level;
       for (i = lbegin; i <= lvlend;  i++) {
            jstrt = xadj[(node = ls[i-1])-1];
            jstop = xadj[ node]-1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ls[iccsze++] = nbr;
                 mask[nbr-1]  = 0;
            }
	   }
   }
   while (iccsze-lvlend > 0);

////////////////////////////////////////////////////////////
//...Reset MASK to one for the nodes in the level structure;
   for (xls[nlvl] = lvlend+1, i = 1; i <= iccsze; i++)
       mask[(node = ls[i-1])-1] = 1;
}

///////////////////////////////////
//...finds pseudo-peripheral nodes;
void fnroot_utils(int & root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int n)
{
   int iccsze, j, jstrt, k, kstop, kstrt, mindeg, ndeg, node, nunlvl;

   rootls_utils(root, xadj, iadj, mask, nlvl, xls, ls, n);
   if (nlvl == 1 || nlvl == (iccsze = xls[nlvl]-1)) return;

   do {
       mindeg = iccsze;
       root   = ls[(jstrt = xls[nlvl-1])-1];

       if (iccsze > jstrt) {
           for (j = jstrt; j <= iccsze; j++) {
                ndeg  = 0;
                kstrt = xadj[(node = ls[j-1])-1];
                kstop = xadj[ node]-1;

                for (k = kstrt; k <= kstop; k++) 
				if (mask[iadj[k-1]-1] > 0 ) ndeg++;

                if (ndeg < mindeg) {
					root   = node;
                    mindeg = ndeg;
                }
           }
       }

       rootls_utils (root, xadj, iadj, mask, nunlvl, xls, ls, n);
       if (nunlvl <= nlvl) return;
   }
   while ((nlvl = nunlvl) < iccsze);
}

//////////////////////////////////////////////////////////////////
//...computes the degrees of the nodes in the connected component;
void degree_utils (int root, int * xadj, int * iadj, int * mask, int * deg, int & iccsze, int * ls, int n)
{
	int i, ideg, j, jstop, jstrt, lbegin, lvlend, nbr, node;

   ls[0]        = root;
   lvlend       = 0;
   iccsze       = 1;
   xadj[root-1] = -xadj[root-1];

   do {
       lbegin = lvlend+1;
       lvlend = iccsze;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = -xadj[(node = ls[i-1])-1];
            jstop =  abs(xadj[node])-1;
            ideg  = 0;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ideg = ideg+1;
                 if (xadj[nbr-1] >= 0) {
                     xadj[nbr-1]  = -xadj[nbr-1];
                     ls[iccsze++] = nbr;
                 }
            }
            deg[node-1] = ideg;
       }
   }
   while (iccsze - lvlend > 0);

///////////////////////////////////////////////
//...Reset XADJ to its correct sign and return;
   for (i = 1; i <= iccsze; i++) {
        node         = ls[i-1];
        xadj[node-1] = -xadj[node-1];
   }
}

////////////////////////////////////////////////
//...reverses the elements of an integer vector;
void ivec_reverse_utils (int n, int * a)
{
  int  m, i;
  for (m = n/2, i = 1; i <= m; i++)
       swap(a[i-1], a[n-i]);
}

/////////////////////////////////////////////////////////////////////////////
//...numbers a connected component using the reverse Cuthill McKee algorithm;
void rcm_utils(int root, int * xadj, int * iadj, int * mask, int * perm, int & iccsze, int * deg, int n)
{
   int fnbr, i, j, jstop, jstrt, k, l, lbegin, lnbr, lperm, lvlend, nbr, node;

   degree_utils (root, xadj, iadj, mask, deg, iccsze, perm, n);
   mask[root-1] = 0;
   if ( iccsze <= 1) return;

   lvlend = 0;
   lnbr   = 1;

   do {
       lbegin = lvlend+1;
       lvlend = lnbr;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = xadj[(node = perm[i-1])-1];
            jstop = xadj[ node]-1;
            fnbr  = lnbr+1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 mask[ nbr-1] = 0;
                 perm[lnbr++] = nbr;
            }
/////////////////////////////////////////////////////////////
//...Sort the neighbors of node in increasing order by degree;
            if (fnbr < lnbr) {
                k = fnbr;
                do {
                    l   = k;
                    nbr = perm[k++];
label40:
                    if (l > fnbr && deg[(lperm = perm[l-1])-1] > deg[nbr-1]) {
                        perm[l--] = lperm;
                        goto label40;
                    }
                    perm[l] = nbr;
                }
                while (k < lnbr);
            }
       }
   }
   while (lnbr > lvlend);

////////////////////////////////////////
//...Reverse the Cuthill-McKee ordering;
   ivec_reverse_utils(iccsze, perm);
}

//////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a general graph;
void genrcm_utils(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask)
{
   int i, iccsze, nlvl, num, root;

   for (num = 1, i = 1; i <= n; i++)
   if  (mask[i-1] != 0) {
       fnroot_utils(root = i, xadj, iadj, mask, nlvl, xls,  perm+num-1,  n);
       rcm_utils   (root,     xadj, iadj, mask, perm+num-1, iccsze, xls, n);

       if ((num += iccsze) > n) return;
   }
}

///////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a given subgraph;
void subrcm_utils(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n)
{
   int  i, iccsze, nlvl, node, num;
   for (i = 1; i <= nsubg; i++)
        mask[(node = subg[i-1])-1] = 1;

   for (num = 0, i = 1; i <= nsubg; i++)
   if  (mask[(node = subg[i-1])-1] > 0 ) {
       fnroot_utils(node, xadj, iadj, mask, nlvl, xls,  perm+num,  n);
       rcm_utils   (node, xadj, iadj, mask, perm+num, iccsze, xls, n);

       if ((num += iccsze) >= nsubg ) return;
   }
}

//////////////////////////////////////////////////////////////
//...алгоритм вычисления с.з. и с.ф. общей матрицы QR-методом;
int HQRfunction(double ** A, double * H_re, double * H_im, int dim_N, int Max_iter)
{//...индексы в циклах были смещены в оригинале на единицу (фортрановский стандрат)!!!
  if (! A || ! H_re || ! H_im) return(0);

  int    i, j, ip, iq, it = 0, l, k, m, last = 0, nn, mmin, sorting = 1;
  double f, s, g, c, r, q, p, radix = 2.,
         sqrdx = radix*radix, t = 1./radix, d = 1./sqrdx, y, x, z, w, v, u;

/////////////
//...balance;
  while (last == 0) for (last = i = 1; i <= dim_N; i++) {
       for (r = c = 0., j = 1; j <= dim_N; j++) if (j != i) {
           c += fabs(A[j-1][i-1]);
           r += fabs(A[i-1][j-1]);
       }
       if (c && r) {
           g = r*t;
           f = 1.;
           s = c+r;
           while (c < g) {
                  f *= radix;
                  c *= sqrdx;
           }
           g = r*radix;
           while (c > g) {
                  f *= t;
                  c *= d;
           }
           if ((c+r)*(g = 1./f) < 0.95*s) {
              for (last = 0, j = 1; j <= dim_N; j++) A[i-1][j-1] *= g;
              for (          j = 1; j <= dim_N; j++) A[j-1][i-1] *= f;
           }
       }
  }

////////////
//...elmhes;
  for (m = 2; m < dim_N; m++) {
      for (x = 0., i = j = m; j <= dim_N; j++)
      if  (fabs(A[j-1][m-2]) > fabs(x)) {
           x = A[j-1][m-2];
           i = j;
      }
      if (i != m) {
         for (j = m-1;j <= dim_N; j++) swap(A[i-1][j-1], A[m-1][j-1]);
         for (j = 1;  j <= dim_N; j++) swap(A[j-1][i-1], A[j-1][m-1]);
      }
      if (x) for (i = m+1; i <= dim_N; i++)
      if ((y = A[i-1][m-2]) != 0.) {
          A[i-1][m-2] = (y /= x);
          for (j = m; j <= dim_N; j++) A[i-1][j-1] -= y*A[m-1][j-1];
          for (j = 1; j <= dim_N; j++) A[j-1][m-1] += y*A[j-1][i-1];
      }
  }

/////////
//...hqr;
  for (f = fabs(A[0][0]), i = 2;   i <= dim_N; i++)
  for (                   j = i-1; j <= dim_N; j++) f += fabs(A[i-1][j-1]);
  nn = dim_N; t = 0.;
  while (nn >= 1) {
         it  = 0;
         do {
            for (l = nn; l >= 2; l--) {
                 s = fabs(A[l-2][l-2])+fabs(A[l-1][l-1]);
                 if (s == 0.) s = f;
                 if (fabs(A[l-1][l-2]) + s == s) break;
            }
            x = A[nn-1][nn-1];
            if (l == nn) {
                H_re[--nn] = x+t;
                H_im[  nn] = 0.;
            }
            else {
                y = A[nn-2][nn-2];
                w = A[nn-1][nn-2]*A[nn-2][nn-1];
                if (l == (nn-1)) {
                    p = .5*(y-x);
                    q = p*p+w;
                    z = sqrt(fabs(q));
                    x += t;
                    if (q >= 0.) {
                        z  = p+(p > 0. ? fabs(z) : -fabs(z));
                        H_re[nn-2] = H_re[nn-1] = x+z; if (z) H_re[nn-1] = x-w/z;
                        H_im[nn-2] = H_im[nn-1] = 0.;
                    }
                    else {
                        H_re[nn-2] =  H_re[nn-1] = x+p;
                        H_im[nn-2] =-(H_im[nn-1] = z);
                    }
                    nn -= 2;
                }
                else {
                    if (it == Max_iter) printf("Too many iterations in HQR");
                    if (it == 10 || it == 20) {
                        t += x;
                        for (i = 1; i <= nn; i++) A[i-1][i-1] -= x;
                        s = fabs(A[nn-1][nn-2])+fabs(A[nn-2][nn-3]);
                        y = x = .75*s;
                        w =-.4375*s*s;
                    }
                    for (it++, m = nn-2; m >= l; m--) {
                         z  = A[m-1][m-1];
                         r  = x-z;
                         s  = y-z;
                         p  = (r*s-w)/A[m][m-1]+A[m-1][m];
                         q  = A[m][m]-z-r-s;
                         r  = A[m+1][m];
                         p *= (s = 1./(fabs(p)+fabs(q)+fabs(r)));
                         q *=  s;
                         r *=  s; if (m == l) break;
                         u  = fabs(A[m-1][m-2])*(fabs(q)+fabs(r));
                         v  = fabs(p)*(fabs(A[m-2][m-2])+fabs(z)+fabs(A[m][m]));
                         if (u+v == v) break;
                    }
                    for (i = m+2; i <= nn; i++) {
                         A[i-1][i-3] = 0.;
                         if (i != (m+2)) A[i-1][i-4] = 0.;
                    }
                    for (k = m; k <= nn-1; k++) {
                         if (k != m) {
                             p  = A[k-1][k-2];
                             q  = A[k][k-2];
                             if (k != nn-1) r = A[k+1][k-2]; else r = 0.;
                             if ((x = fabs(p)+fabs(q)+fabs(r)) != 0.) {
                                p *= (d = 1./x);
                                q *=  d;
                                r *=  d;
                             }
                         }
                         s = sqrt(p*p+q*q+r*r);
                         if ((s = (p > 0. ? s : -s)) != 0.) {
                            if (k != m) A[k-1][k-2] = -s*x; else
                            if (l != m) A[k-1][k-2] = -A[k-1][k-2];
                                p += s;
                                x  = p*(s = 1./s);
                                y  = q*s;
                                z  = r*s;
                                q *= (p = 1./p);
                                r *=  p;
                                for (j = k; j <= nn; j++) {
                                     p = A[k-1][j-1]+q*A[k][j-1];
                                     if (k != (nn-1)) {
                                         p += r*A[k+1][j-1];
                                         A[k+1][j-1] -= p*z;
                                     }
                                     A[k][j-1]   -= p*y;
                                     A[k-1][j-1] -= p*x;
                                }
                                for (mmin = nn < k+3 ? nn : k+3, i = l; i <= mmin; i++) {
                                     p = x*A[i-1][k-1]+y*A[i-1][k];
                                     if (k != nn-1) {
                                         p += z*A[i-1][k+1];
                                         A[i-1][k+1] -= p*r;
                                     }
                                     A[i-1][k]   -= p*q;
                                     A[i-1][k-1] -= p;
                                }
                         }
                    }
                }
            }
         }
         while (l < nn-1);
  }

///////////////////////////////////////////////////////////////////////////
//...упорядочиваем найденные собственные значения в порядке возрастания Re;
  for (iq = 1; sorting && iq < dim_N; iq++) {
       for (d = abs(comp(H_re[(ip = iq)-1], H_im[(ip = iq)-1])), j = iq+1; j <= dim_N; j++) if (abs(comp(H_re[j-1], H_im[j-1])) <= d) d = abs(comp(H_re[(ip = j)-1], H_im[(ip = j)-1]));
       if (ip != iq) {
           swap(H_im[ip-1], H_im[iq-1]);
           swap(H_re[ip-1], H_re[iq-1]);
       }
  }
  return(1);
}
int HQRfunction(dd_real ** A, dd_real * H_re, dd_real * H_im, int dim_N, int Max_iter)
{
  if (! A || ! H_re || ! H_im) return(0);
  return(1);
}
int HQRfunction(qd_real ** A, qd_real * H_re, qd_real * H_im, int dim_N, int Max_iter)
{
  if (! A || ! H_re || ! H_im) return(0);
  return(1);
}
int HQRfunction(complex ** A, complex * H_re, complex * H_im, int dim_N, int Max_iter)
{
  if (! A || ! H_re || ! H_im) return(0);
  return(1);
}

