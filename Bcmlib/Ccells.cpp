#include "stdafx.h"
#include "ccells.h"

/////////////////////////////////////////////////////////
//...test printing surface cells structure in ASCII code;
void CCells::cells_out(FILE * id_CELLS, int id_long)
{
	if (! id_CELLS || ! graph) return;
	if (! mp) {
		fprintf (id_CELLS, "0:| # Bar %d\n", graph[0]);
		for (int k = 0; k < graph[0]; k++) {
			fprintf (id_CELLS, "%d:| ", k);  map_out(id_CELLS, ce[k]->mp, id_long);
			fprintf (id_CELLS, "%d:| ", k); topo_out(id_CELLS, ce[k]->graph);
//////////////////////////////
//...cell in parametric space;
			for (int m = ce[k]->pm && ce[k]->graph ? ce[k]->graph[1] : 0, l = 0; l < m; l++) {
				fprintf(id_CELLS, "%d:| ", k);  map_out(id_CELLS, ce[k]->pm[l], id_long);
			}
      }
	}
	else {
		fprintf (id_CELLS, "0:| # Cells %d\n", graph[0]);
		for (int k = 0; k <= graph[0]; k++) {
			fprintf (id_CELLS, "%d:| ", k);  map_out(id_CELLS, ce[k]->mp, id_long);
			fprintf (id_CELLS, "%d:| ", k); topo_out(id_CELLS, ce[k]->graph);
//////////////////////////////
//...cell in parametric space;
			for (int m = ce[k]->pm && ce[k]->graph ? ce[k]->graph[1] : 0, l = 0; l < m; l++) {
				fprintf(id_CELLS, "%d:| ", k);  map_out(id_CELLS, ce[k]->pm[l], id_long);
			}
      }
	}
}

void CCells::cells_out(char * ch_CELLS, int id_long)
{
	FILE  *   id_CELLS = fopen(ch_CELLS, "w");
	cells_out(id_CELLS, id_long);
	fclose   (id_CELLS);
}

////////////////////////////////////////////////////////
//...test reading surface cells structure in ASCII code;
int  CCells::cells_in(char * id_CELLS, unsigned long & count, unsigned long upper_limit)
{
	if (id_CELLS && user_Count(id_CELLS, count, count, '#') && count < upper_limit) {
		char buf[BUF_SIZE+1];
		int  N_ce = 0, k, l, m, j, id_bar;

		if (user_Read (buf, id_CELLS, count += 1, upper_limit)) {
			if (! ::strcmp(buf, "BAR"  )) id_bar = 1; else
			if (! ::strcmp(buf, "CELLS")) id_bar = 0; else return(0);
		}
		if (user_Read (buf, id_CELLS, count, upper_limit)) N_ce = atoi(buf);

      for (cells_new(N_ce+1), N_ce += (id_bar ? 0 : 1), k = 0; k < N_ce; k++) {
			if (! user_Count(id_CELLS, count, count, '|')) return(0); 
			if (  id_bar || k != N_ce-1) 
			ce[k] = new CCells(-1);
			ce[k]->ce = ce; 
			if ((ce[k]->mp = map_in(id_CELLS, ++count, upper_limit)) == NULL) return(0);

			if (! user_Count(id_CELLS, count, count, '|')) return(0);
			if ((ce[k]->graph = topo_in(id_CELLS, ++count, upper_limit)) == NULL) return(0);
//////////////////////////////
//...cell in parametric space;
			if (2 == map_dim(ce[k]->mp)) {
				if (user_Read (buf, id_CELLS, count, upper_limit)) j = atoi(buf); else j = -1;
				count--;
				if (k == j) {
					if (! user_Count(id_CELLS, count, count, '|')) return(0); 
					ce[k]->pm = (CMap **)new_struct(ce[k]->graph[1]*sizeof(CMap *));
				}
				for (m = ce[k]->pm && ce[k]->graph ? ce[k]->graph[1] : 0, l = 0; l < m; l++) {
					if (user_Read (buf, id_CELLS, count, upper_limit)) j = atoi(buf); else j = -1;
					if (user_Count(id_CELLS, count, count, '|')) {
						if (k !=  j) break;
						ce[k]->pm[l] = map_in(id_CELLS, ++count, upper_limit);
					}
				}
			}
		}
		if (id_bar && N_ce && (graph = (int *)new_struct(2*sizeof(int))) != NULL) graph[0] = N_ce;
	}
	return(1);
}

void CCells::cells_in(char * ch_CELLS)
{
	unsigned long count    = 0, upper_limit;
	char        * id_CELLS = read_struct_ascii(ch_CELLS);
	if         (! id_CELLS) return;
	user_Count   (id_CELLS, 0, upper_limit, '\x0');
	cells_in     (id_CELLS, count, upper_limit);
	delete_struct(id_CELLS);
}

////////////////////////////////////////////////////////////////
//...распределение внутренней структуры пространственной ячейки;
void CCells::cells_new (int N_ce, int N_graph, int N_mp) 
{
	zero_cells();

	ce    = (CCells **)new_struct(N_ce*sizeof(CCells *));
	graph = (Topo    *)new_struct(N_graph*sizeof(Topo));
	mp    = (CMap    *)new_struct(N_mp*sizeof(CMap));
	if (! ce    && N_ce    > 0 ||
		 ! graph && N_graph > 0 ||
		 ! mp    && N_mp    > 0) 
	zero_cells();

	if (mp) mp[N_mp-1] = (CMap)NULL_CELL;
	if (ce) ce[N_ce-1] = this;
}

//////////////////////////////////////////////////////////////////////////
//...удаление из памяти всей внутренней структуры пространственной ячейки;
void CCells::zero_cells()
{
	for (int l, m = graph ? graph[0] : -1, k = 0; k <= m; k++)
	if (ce && ce[k]) {
		delete_struct(ce[k]->mp); l = ce[k]->pm && ce[k]->graph ? ce[k]->graph[1] : 0;
		delete_struct(ce[k]->graph);
		
		for (l--; l >= 0; l--)
			delete_struct(ce[k]->pm[l]);
			delete_struct(ce[k]->pm);

		if (k != m) delete(ce[k]);
	}
	delete_struct(ce);
}

//////////////////////////////////
//...размерность поверхности тела;
int CCells::bar_dim()
{
	int dim = ERR_DIM, k;
	if (! mp && graph)
	for (k = 0; k < graph[0]; k++)
		dim = max(dim, ce[k]->cells_dim());
	return dim;
}

/////////////////////////////////////////////////////////////////////////
//...вспомогательная функция, определяющая число элементов границы блока;
int CCells::arcs_number()
{
  int      id_arc = -1, dim = bar_dim();
  while (++id_arc  < graph[0] && dim == ce[id_arc]->cells_dim());
  return(  id_arc);
}

/////////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция, определяющая ближайший граничный элемент;
int CCells::arcs_number(int id_arc)
{
	int dim = bar_dim();
	while (++id_arc  < graph[0] && dim != ce[id_arc]->cells_dim());
	return(  id_arc == graph[0] ? -1 : id_arc);
}

/////////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция, определяющая ближайшую окружность внутри геометрии;
int CCells::circ_number(int id_arc)
{
	while (++id_arc  < graph[0] && ID_MAP(1, SPHERE_GENUS) != ce[id_arc]->mp[0]);
	return(  id_arc == graph[0] ? -1 : id_arc);
}

////////////////////////////////////////////////////////////////////////
//...относительное перемещение и поворот ячейки (рекурсивная процедура);
void CCells::cells_iso(double * P, double &CZ, double &SZ, double &CY, double &SY, double &CX, double &SX)
{
	map_iso(mp, P, CZ, SZ, CY, SY, CX, SX);
	for (int k = graph[1]+1; k > 1; k--)
		ce[graph[k]]->cells_iso(P, CZ, SZ, CY, SY, CX, SX);
}

////////////////////////////////////////////////
//...относительное перемещение и поворот ячейки;
void CCells::cells_iso(double * P, double fi, double theta, double fX)
{
	if (mp) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		for (int k = 0; k <= graph[0]; k++) 
			if (ce[k]) map_iso(ce[k]->mp, P, CZ, SZ, CY, SY, CX, SX);
	}
}

/////////////////////////////////////////////////////
//...относительное перемещение и поворот поверхности;
void CCells::bar_iso(double * P, double fi, double theta, double fX)
{
	double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
			 CX = cos(fX), SX = sin(fX);
	for (int k = 0; k < graph[0]; k++) 
		if (ce[k]) map_iso(ce[k]->mp, P, CZ, SZ, CY, SY, CX, SX);
}

////////////////////////////////////////////////////
//...идентификация элемента в списках связей ячейки;
int CCells::topo_id(Topo some_element, Topo exc_element)
{
	if (! mp) return(0);
	int id = ::topo_id(graph, some_element), k;
	for (k = graph[1]+1; ! id && k > 1; k--) //...сравниваем списки подэлементов;
	if (graph[k] != exc_element) id = ce[graph[k]]->topo_id(some_element, exc_element);
	return id;
}

///////////////////////////////////////////////
//...cells identification in the cell boundary;
int CCells::common_id(CCells * ext_ce)
{
  if (! ext_ce)									return(1);
  if (ce != ext_ce->ce || ext_ce == this) return(0);
  int id, k;

/////////////////////////////////////////////////////
//...comparison of the elements in the cell boundary;
  for (   id = 0,  k = graph[1]+1; ! id && k > 1; k--)
  for (int l = ext_ce->graph[1]+1; ! id && l > 1; l--) id = graph[k] == ext_ce->graph[l];
  return id*k;
}

//////////////////////////////////////////////////////////////////////////////
//...упорядочивание элементов общего списка поверхности и всех списков связей;
void CCells::bar_ord()
{
	if (! mp)
		for (int  j = 0, k, dim = bar_dim(); dim >= 0; dim--)
			for (k = j; k < graph[0]; k++)
				if (dim == ce[k]->cells_dim()) bar_mutat(k, j++);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...упорядочивание элементов общего списка ячейки (не включенной в список) и всех списков связей;
void CCells::cells_ord()
{
	if (mp) {
		int j, k, m, dim;
		for (     j = graph[0], dim = cells_dim()-1; dim >= 0; dim--)
			for (k = j-1; k >= 0; k--)
				if (dim == ce[k]->cells_dim()) cells_mutat(k, --j);

////////////////////////////////////////////////
//...зеркально отражаем элементы границы ячейки;
		m = 0; while (ce[m]->cells_dim() < cells_dim()-1) m++;
		for (k = 0; k < (graph[0]-m)/2; k++) cells_mutat(m+k, graph[0]-1-k);
	}
}

////////////////////////////////////////////////////////////////////////////////
//...упорядочивание элементов размерности 0 в общем списке ячейки размерности 2;
void CCells::cells_ord0(int id_ord)
{
	int  k, l, j = 0, arc, prev;
	if (id_ord) cells_ord();
	if (cells_dim() == 2) {
		prev = graph[graph[1]+1];
		for (l = 0; l < graph[1]; l++, prev = arc) //...цикл по всем элементам границы ячейки;
			if (ce[arc = graph[l+2]]->graph[1] == 2) {
			if (! ce[prev]->topo_id(k = ce[arc]->graph[2])) k = ce[arc]->graph[3];
			cells_mutat(k, j++);
		}
	}
}

///////////////////////////////////////////////////
//...зеркальное отражение элементов границы ячейки;
void CCells::bar_invers()
{
	int dim, k, m;
	if (! mp) {
		bar_ord();

		dim = bar_dim(); m = graph[0]-1;
		while (m > 0 && ce[m]->cells_dim() < dim) m--;

		for (k = 0; k < (m+1)/2; k++) bar_mutat(k, m-k);
	}
}

////////////////////////////////////////////
//...изменение ориентации одномерной ячейки;
void CCells::cells_invers1()
{
	if (1 == cells_dim()) {
		int N = graph[1], k;

////////////////////////////////////////////////
//...зеркально отражаем элементы границы ячейки;
		for (k = 0; k < N/2; k++) swap(graph[k+2], graph[graph[1]+1-k]);

/////////////////////////////////////////////
//...преобразуем геометрическую карту ячейки;
		if (mp[0] == ID_MAP(1, SPHERE_GENUS) ||
			 mp[0] == ID_MAP(1,   NULL_GENUS)) { //...пpеобpазование оставляет ось X на месте;
			 mp[6]  = M_PI - mp[6];
			 mp[5] += M_PI;
		}
	}
}

///////////////////////////////////////////
//...изменение ориентации двумерной ячейки;
void CCells::cells_invers2()
{
	if (2 == cells_dim()) {
		int N = graph[1], k;

//////////////////////////////////////////////////////////////////////
//...зеркально отражаем и циклически сдвигаем элементы границы ячейки;
		for (k = 0; k < N/2; k++) swap(graph[k+2], graph[N+1-k]);
		//for (k = 0; k < N-1; k++) swap(graph[(k*(N-1))%N+2], graph[((k+1)*(N-1))%N+2]);//...зачем сдвигать элементы на одну позицию???

/////////////////////////////////////////////
//...преобразуем геометрическую карту ячейки;
		if (mp[0] == ID_MAP(2, NULL_GENUS)) { //...пpеобpазование оставляет ось X на месте;
			 mp[6]  = M_PI - mp[6];
			 mp[5] += M_PI;
		}
		else  
		if (mp[0] == ID_MAP(2, SPHERE_GENUS) ||
			 mp[0] == ID_MAP(2,    CYL_GENUS) ||
			 mp[0] == ID_MAP(2,   CONE_GENUS) ||
			 mp[0] == ID_MAP(2,  TORUS_GENUS)) { //...изменяем знак параметра, отвечающего за нормаль;
			 mp[7]  = -mp[7];
		}
	}
}

//////////////////////////////////////////////////////
//...изменение ориентации общей геометрической ячейки;
void CCells::cells_invers(int id_sub_element)
{
	if (mp) {
		switch(cells_dim()) {
			case 1: cells_invers1(); break;
			case 2: cells_invers2(); break;
		}
		if (id_sub_element)
		for (int i = 0; i < graph[1]; i++) ce[graph[i+2]]->cells_invers();
	}
}

////////////////////////////////////////////////////////
//...циклическая сдвижка элементов границы общей ячейки;
void CCells::cells_cyclic_shift(int m)
{
	int k, l, p, N = graph[1];

//////////////////////////////////////////
//...нормируем модуль циклической сдвижки;
	while ((m %= N) <= -1) m += N;
	while ( m       >=  N) m -= N;

/////////////////////////////////////////////////
//...циклически сдвигаем элементы границы ячейки;
	for (l = 0;                    m && l < abs(m); l++)
	for (k = 0, p = max(N/abs(m)-1, 1); k < p;      k++) 
		swap(graph[(k*(N-m)+l)%N+2], graph[((k+1)*(N-m)+l)%N+2]);
}

//////////////////////////////////////////
//...коррекция связей в списках топологии;
void CCells::topo_correct(int N, int N_new, int id_list)
{
	int k, i;
	if (id_list || ! mp) { //...корректируем или устанавливаем номер по общему списку;
		for (k = graph[0]; k >= 0; k--) if (ce[k])
		for (i = ce[k]->graph[1]+1; i > 1; i--)
			if (N == -1 && N_new  == -1 && ce[k]->graph[i] < 0) ++(ce[k]->graph[i]) *= -1; else
			if (N ==  ce[k]->graph[i]) ce[k]->graph[i] = N_new;
	}
	else 
	if (N == -1 && N_new == -1) { //...корректируем номер (для ячейки) по структуре
		for (i = 0; i < graph[1]; i++)
			if (graph[i+2] < 0) ++(graph[i+2]) *= -1; 
			else ce[graph[i+2]]->topo_correct(N, N_new);
	}
	else //...устанавливаем номер (для ячейки) по структуре;
	for (i = 0; i < graph[1]; i++)
	if (graph[i+2] >= 0) {
		if (graph[i+2] == N) graph[i+2] = N_new; 
		else ce[graph[i+2]]->topo_correct(N, N_new);
	}
}

/////////////////////////////////////////////////////////////////
//...идентификация ячейки (для реверсированных ячеек -- id = -1);
int CCells::cells_id(CCells * ext_ce, int id_num_correct)
{//...правило: ячейка с отрицательной границей в пределах структуры -- идентифицирована,
 //...с отрицательной границей вне структуры -- не идентифицирована;
	if (! ext_ce)				  return(1);
	if (! ext_ce->mp || ! mp) return(0);

/////////////////////////////
//...проверка границы ячейки;
	int id = graph[1] == ext_ce->graph[1], i, k, l, m_ce;
	for (k = graph[1]+1; id && k > 1; k--) if (ext_ce->graph[k] <= -graph[0]-1) id = 0;
  
//////////////////////////////////////////
//...проверка геометрической карты ячейки;  
	if (id)
      id *= map_id(mp, ext_ce->mp);

/////////////////////////////////////////////////////////////////////////
//...сравнение элементов границы ячейки и маркировка совпавших элементов;
	if (id) {
		for (l = ext_ce->graph[1]+1; id && l > 1; l--) {
			for (k = graph[1]+1; k > 1; k--) if (graph[k] >= 0)
				if ((m_ce  = ext_ce->graph[l]) < 0 && (i = (m_ce+graph[k]+1) == 0) != 0 ||
					  m_ce >= 0 && (i = ce[graph[k]]->cells_id(ext_ce->ce[m_ce])) != 0) break;

			id *= abs(i);
			if (k > 1) { //...маркируем границу и отмечаем в ячейке совпавшие элементы;
				graph[k] = -graph[k]-1;
				if (id_num_correct && m_ce >= 0) { 
					ext_ce->ce[ext_ce->graph[0]]->topo_correct(m_ce, graph[k], NULL_STATE);
					if (m_ce >= 0) ext_ce->graph[l] = m_ce;
				}
			}
		}

//////////////////////////////////////////
//... восстановление отмеченных элементов;
		for (k = graph[1]+1; k > 1; k--) 
			if (graph[k] < 0)	graph[k] = -graph[k]-1;
	}
	return id;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...определение номера ячейки внутри структуры (с отрицательным знаком для противоположно ориентированных);
int CCells::cells_in_struct(CCells * ext_ce, int id_num_correct)
{//...правило: нулевая ячейка -- идентифицируется автоматически;
	if (! ext_ce)		return(1);
	if (! ext_ce->mp) return(0);

	int N, i;
	if (! mp)
	for ( N = i = 0; i < graph[0]; i++) {
			N = i; if ((++N *= ce[i]->cells_id(ext_ce, id_num_correct)) != 0) break;
	}
	else {
		for ( N = i = 0; i < graph[1]; i++) {
				N = graph[i+2]; if ((++N *= ce[i]->cells_id(ext_ce, id_num_correct)) != 0) break;
		}
		if (! N)
			for (i = 0; i < graph[1]; i++)
				if ((N = ce[graph[i+2]]->cells_in_struct(ext_ce, id_num_correct)) != 0) break;
	}
	return N;
}

////////////////////////////////////////////////////////////////
//...вспомогательная функция -- удаление элементов подструктуры;
void CCells::delete_cells_in_struct()
{
	for ( int k = 0; k < graph[1]; k++) if (graph[k+2] >= 0 && ce[graph[k+2]]) {
			int l = pm ? graph[1] : 0;
			delete_struct(mp);

			for (l--; l >= 0; l--)
				delete_struct(pm[l]);
				delete_struct(pm);

			ce[graph[k+2]]->delete_cells_in_struct();

			delete_struct(ce[graph[k+2]]->graph);
			delete_cells (ce[graph[k+2]], NULL_STATE);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//...идентификация элементов подструктуры, внесение новых элементов в список и удаление повторяющихся;
void CCells::search_cells_in_struct(CCells * ext_ce, int & i, int id_cell, CCells **& dop_ce, int & N_buf, int & buf_size, int buf_delta)
{
	int k, m, l, m_ce;
	for ( k = 0; k < graph[1]; k++) if ((m_ce = graph[k+2]) >= 0 && ce[m_ce]) {
			m = id_cell != NULL_STATE  ? ext_ce->cells_in_struct(ce[m_ce], OK_STATE) : 0;
			if (! m) { //...запоминаем новый номер в структуре поверхности;
				if (buf_size == N_buf) {
					CCells ** new_ce = (CCells **)new_struct((buf_size += buf_delta)*sizeof(CCells *));
					memcpy	(new_ce, dop_ce, N_buf*sizeof(CCells *));	delete_struct(dop_ce); dop_ce = new_ce;
				}
				ce[m_ce]->search_cells_in_struct(ext_ce, i, id_cell, dop_ce, N_buf, buf_size, buf_delta);
				ce[graph[0]]->topo_correct(m_ce, -(++i), NULL_STATE); dop_ce[N_buf++] = ce[m_ce];
			}
			else if (ce[m_ce]) { //...удаляем повторяющиеся структуры;
				l = ce[m_ce]->pm ? ce[m_ce]->graph[1] : 0;
				delete_struct(ce[m_ce]->mp);

				for (l--; l >= 0; l--)
					delete_struct(ce[m_ce]->pm[l]);
					delete_struct(ce[m_ce]->pm);

				ce[m_ce]->delete_cells_in_struct();
				ce[graph[0]]->topo_correct(m_ce, -abs(m), NULL_STATE);
			
				delete_struct(ce[m_ce]->graph);
				delete_cells (ce[m_ce], NULL_STATE);
			}
	}
}

/////////////////////////////////////////////////////////////////////////
//...вспомогательная функция -- построение списка элементов подструктуры;
void CCells::search_cells_element(int *& id, int & N_buf, int & buf_size, int buf_delta)
{
	if (graph)
	for ( int k = 0; k < graph[1]; k++) {
			if (buf_size == N_buf) {
				int * new_id = (int *)new_struct((buf_size += buf_delta)*sizeof(int));
				memcpy(new_id, id, N_buf*sizeof(int));
				delete_struct(id); id = new_id;
			}
			id[N_buf++] = graph[k+2];
			ce[graph[k+2]]->search_cells_element(id, N_buf, buf_size, buf_delta);
	}
}

///////////////////////////////////////////////////////////////////////////////
//...cells including in the surface structure (return number of included cell);
int CCells::bar_add(CCells * ext_ce, int id_cell, int buf_delta)
{
///////////////////////////////////////////////////
//...начальные проверки и инициализация параметров;
	if (! ext_ce) return(1);
	if (! ext_ce->mp || mp) return(0);

	int  i = graph[0], k, l, m, N, N_buf = 0, buf_size = 0;
	CCells ** dop_ce = (CCells **)new_struct((buf_size += buf_delta)*sizeof(CCells *));

///////////////////////////////////////////////////////////////////
//...идентификация геометрической ячейки в целом и коррекция связей;
	m = id_cell != NULL_STATE ? cells_in_struct(ext_ce, OK_STATE) : 0;

////////////////////////////////////////////////////
//...запоминаем новый номер в структуре поверхности;
	if (! m) {
		dop_ce[N_buf++] = ext_ce; N = ++i;
		ext_ce->search_cells_in_struct(this, i, id_cell, dop_ce, N_buf, buf_size, buf_delta);
		ext_ce->topo_correct();
	}

////////////////////////////////////
//...удаляем повторяющуюся структур;
	else {
		N = m; l = ext_ce->pm ? ext_ce->graph[1] : 0;
		delete_struct(ext_ce->mp);

		for (l--; l >= 0; l--)
			delete_struct(ext_ce->pm[l]);
			delete_struct(ext_ce->pm);

		ext_ce->delete_cells_in_struct();

		delete_struct(ext_ce->graph);
		delete_cells (ext_ce, NULL_STATE);
	}

/////////////////////////////////////////////
//...формирование нового общего списка ячеек;
	if (i > graph[0]) {
		CCells ** new_ce = (CCells **)new_struct((i+1)*sizeof(CCells *));
		for (k = 0; k < graph[0]; k++) {
			new_ce[k] = ce[k]; new_ce[k]->ce = new_ce;
		}
		for (m = 0; m < N_buf; m++) {
			new_ce[k] = dop_ce[m]; new_ce[k++]->ce = new_ce;
		}
		delete_struct(ce); ce = new_ce;
		ce[graph[0] = i] = this;
	}
	install_struct();	

//////////////////////////
//...выходим из программы;
	delete_struct(dop_ce);
	return abs(N);
}

/////////////////////////////////////////////////////////////////////
//...parametric cells description including in the surface structure;
void CCells::trim_add(CCells *& trim, int id_mp)
{
	if (trim && mp && trim->mp && ! pm && graph && trim->graph && (graph[1] <= trim->graph[1])) {
		int k, N_pad;

////////////////////
//...align topology;
		if  (topo_pad(graph, N_pad = trim->graph[1]-graph[1]))
		for (k = 0; k < N_pad;  k++) graph[k+2] = graph[N_pad+2]/*ERR_GRAPH*/;

/////////////////////////////////////////////
//...allocation and transfer parametric maps;
		for (pm = (CMap **)new_struct(graph[1]*sizeof(CMap *)), k = 0; pm && k < graph[1]; k++)
			swap(pm[k], trim->ce[trim->graph[k+2]]->mp);

///////////////////////////////////////////////////////
//...transfer parametrization and delete rest elements;
		if (id_mp) {
			delete_struct(mp); swap(mp, trim->mp);
		}
		delete_cells(trim);
	}
}

/////////////////////////////////////////////////////
//...натягивание на повеpхность геометpической каpты;
int CCells::bar_span(CMap *& ext_mp)
{
	if (! ext_mp || mp) return(1);

	int  k, j = 1, m = map_dim(ext_mp);
	if ((k = bar_dim()) != m-1 && k != ERR_DIM) j = 0;

///////////////////////////////////////////////////////////////////////////////////////////
//...просматриваем все элементы размерности на единицу меньше и заносим их в список связей;
	for (k = 0; j && k < graph[0]; k++)
		if (ce[k]->cells_dim() == m-1 && ! topo_insert(graph, k)) j = 0;

////////////////////////////////////////////////////////////////////////////////////
//...подключаем геометpическую каpту, упорядочиваем элементы и выходим из пpогpаммы;
	mp = ext_mp; ext_mp = NULL; cells_ord();
	return(j);
}

//////////////////////////////////////////////////////////////
//...копирование ячейки по номеру (в том числе и всей ячейки);
CCells * CCells::bar_cpy(int N, int id_origine, int buf_delta)
{
	if (--N >= 0 && N <= graph[0] && ce[N]->mp) {
		int * id, i, k, N_buf = 0, buf_size = 0;

////////////////////////////////////////////////////////////////////////
//...формируем головной элемент копируемой ячейки и список подэлементов;
		CCells * cpy_ce = new CCells(-1); if (! cpy_ce) return(NULL);
		cpy_ce->graph	 = graph_cpy(N);
		cpy_ce->mp		 = map_cpy(N);
		cpy_ce->pm		 = pm_cpy (N);

		id = (int *)new_struct((buf_size += buf_delta)*sizeof(int));
		ce[N]->search_cells_element(id, N_buf, buf_size, buf_delta);

/////////////////////////////////////////////////////////////////////
//...отмечаем и вычеркиваем повторяющиеся элементы (с помощью маски);
		int min_id = id[0], max_id = id[0];
		for (k = 1; k < N_buf; k++)
		if (id[k] < min_id) min_id = id[k]; else
		if (id[k] > max_id) max_id = id[k];

		int * mask = (int *)new_struct((max_id-min_id+1)*sizeof(int));
		for (k = 0; k < N_buf; k++)
		if ( mask[id[k]-min_id]) id[k] = -id[k]-1; 
		else mask[id[k]-min_id]++;

		for (k = 0; k < N_buf; k++) 
		if (id[k] < 0) {
			for (i = N_buf-1; i > k; i--) 
			if (id[i] >= 0) {
				swap(id[i], id[k]); break;
			}
			N_buf = i;; 
		}
		delete_struct(mask);

///////////////////////////////////////////////
//...заполняем данные общего списока элементов;
		cpy_ce->ce = (CCells **)new_struct((N_buf+1)*sizeof(CCells *));
		for (i = 0; i < N_buf; i++) {
			cpy_ce->ce[i] = new CCells(-1);
			cpy_ce->ce[i]->graph = graph_cpy(id[i]);
			cpy_ce->ce[i]->mp    = map_cpy  (id[i]);
			cpy_ce->ce[i]->pm    = pm_cpy   (id[i]);
			cpy_ce->ce[i]->ce    = cpy_ce->ce;
			cpy_ce->ce[i]->graph[0] = N_buf;
		}
		cpy_ce->ce[N_buf] = cpy_ce; 
		cpy_ce->graph[0]  = N_buf;

////////////////////////////////////
//...перенумеровываем списки связей;
		for (i = 0; i < N_buf; i++) cpy_ce->topo_correct(id[i], -(i+1));
		cpy_ce->topo_correct();
		cpy_ce->cells_ord();
		
		if (id_origine) cpy_ce->cells_to_origine();

//////////////////////////
//...выходим из программы;
		delete_struct(id); return(cpy_ce);
	}
	return(NULL);
}

//////////////////////////////////////////////////////////
//...extraction cell from surface description (by number);
CCells * CCells::bar_sub(int N, int id_origine, int buf_delta)
{
  CCells * ext_ce = bar_cpy(N, id_origine, buf_delta);
  if      (ext_ce)  bar_exc(N);
  return  (ext_ce);
}

//////////////////////////////////////
//...копирование всей структуры ячеек;
CCells * CCells::bar_cpy()
{
	CCells * bar_N = new CCells(-1);
	bar_N->cells_new(1, 2, 0);

	for (int i = 0; i < graph[0]; i++)
		bar_N->bar_add(bar_cpy(i+1), OK_STATE);

	bar_N->bar_ord();	return bar_N;
}

/////////////////////////////////////////////
//...удаление ячейки из описания поверхности;
void CCells::bar_exc(int N)
{
	if (N-- > graph[0] || N < 0 || mp) return;
	int i, j, k, l; 

///////////////////////////////////////////////////////////////////
//...просматриваем общий список ячеек и удаляем не нужные элементы;
	for (k = 0; k < graph[0]; k++) if (ce[N]->topo_id(k)) {
		for (j = i = 0; ! j && i < graph[0]; i++)
		if ( i != N) j = ce[i]->topo_id(k, N);

////////////////////////////////////////////////////////////////////////
//...удаляем элемент, который не входит во все ячейки, кроме bar->ce[N];
      if (! j && k != N) {
			l = ce[k]->pm && ce[k]->graph ? ce[k]->graph[1] : 0;
			delete_struct(ce[k]->mp); 

			for (l--; l >= 0; l--)
				delete_struct(ce[k]->pm[l]);
				delete_struct(ce[k]->pm);

         delete_struct(ce[k]->graph);
         delete_cells (ce[k], NULL_STATE);

///////////////////////////////////////////////////////////
//...изменяем нумерацию и удаляем элемент из общего списка;
         for (i = k; i < graph[0];	 i++) ce[i] = ce[i+1]; graph[0] -= 1;
         for (i = 0; i < graph[0];	 i++) topo_exc (ce[i]->graph, k);
			for (i = k; i < graph[0]+1; i++) topo_correct(i+1, i);
		}
	}

///////////////////////////
//...удаляем элемент ce[N];
	l = ce[N]->pm && ce[N]->graph ? ce[N]->graph[1] : 0;
	delete_struct(ce[N]->mp); 

	for (l--; l >= 0; l--)
		delete_struct(ce[N]->pm[l]);
		delete_struct(ce[N]->pm);

	delete_struct(ce[N]->graph);
	delete_cells (ce[N], NULL_STATE);

/////////////////////////////////////////////////////////////////
//...изменяем нумерацию и удаляем элемент ce[N] из общего списка;
	for (i = N; i < graph[0];	 i++) ce[i] = ce[i+1]; graph[0] -= 1;
	for (i = 0; i < graph[0];	 i++) topo_exc (ce[i]->graph, N);
	for (i = N; i < graph[0]+1; i++) topo_correct(i+1, i);
 
///////////////////////////////////////////////////
//...устанавливаем структуру, выходим из программы;
	install_struct();
}

//////////////////////////////////////////////////////////////////
//...устанавливаем размерность общей структуры для всех элементов;
void CCells::install_struct()
{
	if (! mp)
	for ( int i = 0; i < graph[0]; i++) 
		if (ce[i]) ce[i]->graph[0] = graph[0];
}

//////////////////////////////////////////////////////////////////
//...excluding 1D or 2D degenerate elements from boundary surface;
void CCells::bar_generate(double eps)
{
	for (int m, l, k = graph[0]-1; k >= 0; k--)
	if (1 == ce[k]->cells_dim() && 
		(1 == ce[k]->graph[1] ||
		 2 == ce[k]->graph[1] && ce[k]->graph[2] == ce[k]->graph[3])) {
		if (ce[k]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
			if (fabs(ce[k]->mp[7]) < eps ||
				 fabs(ce[k]->mp[8]) < eps) bar_exc(k+1);
		}
		else
		if (ce[k]->mp[0] == ID_MAP(1, NULL_GENUS)) {
			bar_exc(k+1);
		}
	}
//////////////////////////////////////////////////////////
//...excluding 2D degenerate elements (remove repetition);
	else
	if (2 == ce[k]->cells_dim())
	for (l = ce[k]->graph[1]-1; l > 0; l--) if (ce[k]->graph[2+l] == ce[k]->graph[1+l]) {
		for (m = l; m < ce[k]->graph[1]; m++) {
			if (ce[k]->pm)
			ce[k]->pm[m-1]		= ce[k]->pm[m];
			ce[k]->graph[1+m] = ce[k]->graph[2+m];
		}
		ce[k]->graph[1]--;
	}
}

///////////////////////////////////////////////////////////////////
//...исключение точек с двойной нумерацией из описания поверхности;
void CCells::bar_correct_double_point(double eps)
{
	for (int k = graph[0]-1, l, m; k >= 0; k--)
		if (1 == ce[k]->cells_dim() && 
			(2 == ce[k]->graph[1] && 
			(l  = ce[k]->graph[2]) != (m = ce[k]->graph[3]) && max_point(ce[l]->mp+1, ce[m]->mp+1) < eps)) {
			topo_correct(m, l);
			bar_exc(m+1);
		}
}

/////////////////////////////////////////////////////////////////
//...исключение точек с двойной нумерацией из описания 1D ячейки;
void CCells::cell_correct_double_point(double eps)
{
	int l, m;
	if (1 == cells_dim() && 
		(2 == graph[1] && 
		(l  = graph[2]) != (m = graph[3]) && max_point(ce[l]->mp+1, ce[m]->mp+1) < eps)) {

			for (topo_correct(m, l); m < graph[0]; m++) {
				CCells * temp_ce = ce[m]; ce[m] = ce[m+1]; ce[m+1] = temp_ce;
			}
			graph[0]--;

			delete_struct(ce[m]->mp);
			delete_struct(ce[m]->graph);
			delete_cells (ce[m]);

			install_struct();
	}
}

/*=================================================================*/
/*                  ИДЕНТИФИКАЦИЯ ЭЛЕМЕНТОВ                        */
/*=================================================================*/
/////////////////////////////////////////////////////////////////
//...идентификация геометрической ячейки "плоский прямоугольник";
int CCells::id_sheet()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
		graph && graph[1] == 4                     &&
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1, NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1, NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, NULL_GENUS)) {
      double A = ce[m1]->cells_length(),
             B = ce[m2]->cells_length(), ort1[9], ort2[9], ort3[9];

///////////////////////////////////////////////////////
//...проверяем параллельность и ортогональность сторон;
		ce[m1]->line_correct();	map_normalizat(ce[m1]->mp, NULL, ort1);
		ce[m3]->line_correct();	map_normalizat(ce[m3]->mp, NULL, ort2);
		if (fabs(fabs(ort1[0]*ort2[0]+ort1[1]*ort2[1]+ort1[2]*ort2[2])-1.) >= EE_ker) return(0);

		ce[m2]->line_correct();	map_normalizat(ce[m2]->mp, NULL, ort3);
		ce[m4]->line_correct();	map_normalizat(ce[m4]->mp, NULL, ort2);
		if (fabs(ort1[0]*ort2[0]+ort1[1]*ort2[1]+ort1[2]*ort2[2]) >= EE_ker ||
	       fabs(fabs(ort3[0]*ort2[0]+ort3[1]*ort2[1]+ort3[2]*ort2[2])-1.) >= EE_ker) return(0);

///////////////////////////////////////////////////////////
//...добавляем параметры (и копируем геометрическую карту);
		CMap * map = (CMap *)new_struct(((m = size_of_map(2, NULL_GENUS))+1)*sizeof(CMap));
		if ( ! map || ! add_new_maps(mp, m+1, 2)) {
			delete[] map; return(0);
		}
		memcpy(map, mp, (m+1)*sizeof(CMap));
		mp[  m] = (CMap)SHEET_CELL;
		mp[++m] = (CMap)A;
		mp[++m] = (CMap)B;

//////////////////////////////////////////////
//...корректируем локальную систему координат;
		if (ce[m4]->topo_id(ce[m1]->graph[3])) {
			ort1[0] = -ort1[0];
			ort1[1] = -ort1[1];
			set_point(ort2, ce[ce[m1]->graph[2]]->mp+1);
		}
		else set_point(ort2, ce[ce[m1]->graph[3]]->mp+1);

		if ( ce[m2]->topo_id(ce[m3]->graph[3]))
			  set_point(ort3, ce[ce[m3]->graph[2]]->mp+1);
		else set_point(ort3, ce[ce[m3]->graph[3]]->mp+1);

		mp[1] = mp[2] = mp[3] =
		mp[4] = mp[5] = mp[6] = map[1] = map[2] = map[3] = 0.; make_local(ort1, map);

		ort2[0] = (ort2[0]+ort3[0])*.5;
		ort2[1] = (ort2[1]+ort3[1])*.5;
		ort2[2] = (ort2[2]+ort3[2])*.5; set_point(map+1, ort2);

		map_iso   (mp, NULL, arg2(comp(ort1[0], ort1[1])));
		map_common(mp, map);

		delete[] map;
	}
	return(1);
}

///////////////////////////////////////
//...идентификация кольцевого сегмента;
int CCells::id_ring_segment()
{
	int m1, m2, m3, m4, k1, k2, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && ( //...четыре ортогональных дуги;
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS))) {

//////////////////////////////////////////////
//...проверяем параллельность дуг окружностей;
		if (fabs(ce[m2]->mp[1]-ce[m4]->mp[1]) >= EE_ker ||
	       fabs(ce[m2]->mp[2]-ce[m4]->mp[2]) >= EE_ker ||
	       fabs(ce[m2]->mp[3]-ce[m4]->mp[3]) >= EE_ker) return(0);

////////////////////////////
//...проверяем длины прямых;
		k1 = get_num(ce[m1]->graph, 0);
		k2 = get_num(ce[m1]->graph, 1);
		if (fabs(abs_point(ce[k1]->mp+1, ce[k2]->mp+1)-
			 fabs(fabs(ce[m2]->mp[7])-fabs(ce[m4]->mp[7]))) >= EE_ker) return(0);

		k1 = get_num(ce[m3]->graph, 0);
		k2 = get_num(ce[m3]->graph, 1);
		if (fabs(abs_point(ce[k1]->mp+1, ce[k2]->mp+1)-
			 fabs(fabs(ce[m2]->mp[7])-fabs(ce[m4]->mp[7]))) >= EE_ker) return(0);

////////////////////////////////////////////////////////////////////
//...добавляем параметры и корректируем локальную систему координат;
		if (fabs(ce[m2]->mp[7]) < fabs(ce[m4]->mp[7])) swap(m2, m4);
		if (add_new_maps(mp, (m = size_of_map(2, NULL_GENUS))+1, size_of_dop(RING_SEGMENT))) {
			mp[  m] = (CMap)RING_SEGMENT;
			mp[++m] = (CMap)fabs(ce[m2]->mp[8]*2.);
			mp[++m] = (CMap)fabs(ce[m4]->mp[7]);
			mp[++m] = (CMap)fabs(ce[m2]->mp[7]);
		}
		set_point(mp+1, ce[m2]->mp+1);
		set_point(mp+4, ce[m2]->mp+4);
	}
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 3 && ( //...две прямые и ортогональная дуга;
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS))) {
		if ((m4 = ce[m1]->common_id(ce[m3])) != 0) m4 = ce[m1]->graph[m4]; else return(1);

///////////////////////////////////////////////
//...проверяем ортогональность дуги окружности;
		if (fabs(ce[m2]->mp[1]-ce[m4]->mp[1]) >= EE_ker ||
	       fabs(ce[m2]->mp[2]-ce[m4]->mp[2]) >= EE_ker ||
	       fabs(ce[m2]->mp[3]-ce[m4]->mp[3]) >= EE_ker) return(0);

////////////////////////////////////////////////////////////////////
//...добавляем параметры и корректируем локальную систему координат;
		if (add_new_maps(mp, (m = size_of_map(2, NULL_GENUS))+1, size_of_dop(RING_SEGMENT))) {
			mp[  m] = (CMap)RING_SEGMENT;
			mp[++m] = (CMap)fabs(ce[m2]->mp[8]*2.);
			mp[++m] = (CMap)0.;
			mp[++m] = (CMap)fabs(ce[m2]->mp[7]);
		}
		set_point(mp+1, ce[m2]->mp+1);
		set_point(mp+4, ce[m2]->mp+4);
	}
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 3 && ( //...две совпадающих дуги и прямая;
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1, SPHERE_GENUS))) {

//////////////////////////////////////////
//...проверяем совпадение дуг окружностей;
		if (! map_id(ce[m2]->mp, ce[m3]->mp)) return(0);

////////////////////////////
//...проверяем длины прямых;
		k1 = get_num(ce[m1]->graph, 0);
		k2 = get_num(ce[m1]->graph, 1);
		if (fabs(abs_point(ce[k1]->mp+1, ce[k2]->mp+1)-fabs(ce[m2]->mp[7]*2.)) >= EE_ker) return(0);

////////////////////////////////////////////////////////////////////
//...добавляем параметры и корректируем локальную систему координат;
		if (fabs(ce[m2]->mp[7]) < fabs(ce[m4]->mp[7])) swap(m2, m4);
		if (add_new_maps(mp, (m = size_of_map(2, NULL_GENUS))+1, size_of_dop(RING_SEGMENT))) {
			mp[  m] = (CMap)RING_SEGMENT;
			mp[++m] = (CMap)M_PI;
			mp[++m] = (CMap)0.;
			mp[++m] = (CMap)fabs(ce[m2]->mp[7]);
		}
		set_point(mp+1, ce[m2]->mp+1);
		set_point(mp+4, ce[m2]->mp+4);

		map_local (mp, ce[m2]->mp);
		mp[4] -=  (M_PI_2-ce[m2]->mp[8])*cos(ce[m2]->mp[5]);
		map_common(mp, ce[m2]->mp);
	}
	return(1);
}

///////////////////////////////////////////////////////////
//...идентификация прямоугольного цилиндрического сегмента;
int CCells::id_cyl_segment()
{
	int m1, m2, m3, m4, k, m;
	if (mp && mp[0] == ID_MAP(2, CYL_GENUS) && graph)
	if (graph[1] == 0 && add_new_maps(mp, (m = size_of_map(2, CYL_GENUS))+1, size_of_dop(CYL_SEGMENT))) { //...добавляем параметры;
		mp[  m] = (CMap)CYL_SEGMENT;
		mp[++m] = (CMap)M_PI*2.;
		mp[++m] = (CMap)fabs(mp[7])*2.;//...задаем длину, равную радиусу;
	}
	else
	if (graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS))) {
      double A = ce[m1]->cells_length(), ort1[9], ort2[9], P[3], f;

///////////////////////////////////////////////////////////
//...добавляем параметры (и копируем геометрическую карту);
		CMap * map = (CMap *)new_struct(((m = size_of_map(2, CYL_GENUS))+1)*sizeof(CMap));
		if (! map || ! add_new_maps(mp, m+1, size_of_dop(CYL_SEGMENT))) {
			delete[] map; return(0);
		}
		for ( k = 0; k <= m; k++) map[k] = mp[k];
		mp[  m] = (CMap)CYL_SEGMENT;
		mp[++m] = (CMap)ce[m2]->mp[8]*2.;
		mp[++m] = (CMap)A;

//////////////////////////////////////////////
//...корректируем локальную систему координат;
		map_normalizat(ce[m1]->mp, NULL, ort1);
		map_normalizat(ce[m2]->mp, NULL, ort2); map[1] = map[2] = map[3] = 0.;

		make_local(ort1,   map);
		make_local(ort2,   map);
		make_local(ort2+6, map); set_point(map+1, mp+1); f = arg2(comp(ort2[6], ort2[7]));
		mp[1] = mp[2] = mp[3] =
		mp[4] = mp[5] = mp[6] = 0.;

		if (ce[m4]->topo_id(ce[m1]->graph[3])) {
			set_point(P, ce[ ce[m1]->graph[2]]->mp+1);
			ort1[2] = -ort1[2];
		}
		else set_point(P, ce[ce[m1]->graph[3]]->mp+1);
		make_local    (P, mp);

		k     = fabs(1.+ort1[2]*ort2[2]) < EE_dop; m = mp[7] < 0.; if (m) k = 1-k;
		P[0]  =
		P[1]  = 0.;
		P[2] -= A*(k ? .5 : -.5); make_common(P, map); set_point(map+1, P);

		map_iso   (mp, NULL, 0., k ? 0. : M_PI, k ? f : M_PI-f);
		map_common(mp, map);

		delete[] map;
  }
  return(1);
}

////////////////////////////////////////////////////////
//...идентификация прямоугольного сферического сегмента;
int CCells::id_sph_segment()
{
	int /*m1, m2, m3, m4, k, */m;
	if (mp && mp[0] == ID_MAP(2, SPHERE_GENUS) && graph)
	if (graph[1] == 0 && add_new_maps(mp, (m = size_of_map(2, SPHERE_GENUS))+1, size_of_dop(SPH_SEGMENT))) { //...добавляем параметры;
		mp[  m] = (CMap)SPH_SEGMENT;
		mp[++m] = (CMap)M_PI*2.;
		mp[++m] = (CMap)0.;
		mp[++m] = (CMap)M_PI;
	}
#ifdef ___CONSTRUCTION___
	else
	if (graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS))) {
      double A = ce[m1]->cells_length(), ort1[9], ort2[9], P[3], f;

///////////////////////////////////////////////////////////
//...добавляем параметры (и копируем геометрическую карту);
		CMap * map = (CMap *)new_struct(((m = size_of_map(2, CYL_GENUS))+1)*sizeof(CMap));
		if (! map || ! add_new_maps(mp, m+1, size_of_dop(SPH_SEGMENT))) {
			delete[] map; return(0);
		}
		for ( k = 0; k <= m; k++) map[k] = mp[k];
		mp[  m] = (CMap)CYL_SEGMENT;
		mp[++m] = (CMap)ce[m2]->mp[8]*2.;
		mp[++m] = (CMap)A;

//////////////////////////////////////////////
//...корректируем локальную систему координат;
		map_normalizat(ce[m1]->mp, NULL, ort1);
		map_normalizat(ce[m2]->mp, NULL, ort2); map[1] = map[2] = map[3] = 0.;

		make_local(ort1,   map);
		make_local(ort2,   map);
		make_local(ort2+6, map); set_point(map+1, mp+1); f = arg2(comp(ort2[6], ort2[7]));
		mp[1] = mp[2] = mp[3] =
		mp[4] = mp[5] = mp[6] = 0.;

		if (ce[m4]->topo_id(ce[m1]->graph[3])) {
			set_point(P, ce[ ce[m1]->graph[2]]->mp+1);
			ort1[2] = -ort1[2];
		}
		else set_point(P, ce[ce[m1]->graph[3]]->mp+1);
		make_local    (P, mp);

		k     = fabs(1.+ort1[2]*ort2[2]) < EE_dop; m = mp[7] < 0.; if (m) k = 1-k;
		P[0]  =
		P[1]  = 0.;
		P[2] -= A*(k ? .5 : -.5); make_common(P, map); set_point(map+1, P);

		map_iso   (mp, NULL, 0., k ? 0. : M_PI, k ? f : M_PI-f);
		map_common(mp, map);

		delete[] map;
	}
#endif
	return(1);
}

/////////////////////////////////////////////////////
//...частичная идентификация сфероидального сегмента;
int CCells::id_spr_segment()
{
	int /*m1, m2, m3, m4, k,*/ m;
	if (mp && mp[0] == ID_MAP(2, SPHEROID_GENUS) && graph)
	if (graph[1] == 0 && add_new_maps(mp, (m = size_of_map(2, SPHEROID_GENUS))+1, size_of_dop(SPR_SEGMENT))) { //...добавляем параметры;

		mp[  m] = (CMap)SPR_SEGMENT;
		mp[++m] = (CMap)M_PI*2.;
		mp[++m] = (CMap)0.;
		mp[++m] = (CMap)M_PI;
	}
	return(1);
}

//////////////////////////////////////////////////
//...идентификация прямоугольного сегмента конуса;
int CCells::id_cone_segment()
{
	int /*m1, m2, m3, m4, k, */m;
	if (mp && mp[0] == ID_MAP(2, CONE_GENUS) && graph)
	if (graph[1] == 0 && add_new_maps(mp, (m = size_of_map(2, CONE_GENUS))+1, size_of_dop(CONE_SEGMENT))) { //...добавляем параметры;
		mp[  m] = (CMap)CONE_SEGMENT;
		mp[++m] = (CMap)M_PI*2.;
		mp[++m] = (CMap)0.;
		mp[++m] = (CMap)1.;//...задаем длину, равную единице;
	}
#ifdef ___CONSTRUCTION___
	else
	if (graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1, SPHERE_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] == ID_MAP(1, SPHERE_GENUS))) {
      double A = ce[m1]->cells_length(), ort1[9], ort2[9], P[3], f;

///////////////////////////////////////////////////////////
//...добавляем параметры (и копируем геометрическую карту);
		CMap * map = (CMap *)new_struct(((m = size_of_map(2, CYL_GENUS))+1)*sizeof(CMap));
		if (! map || ! add_new_maps(mp, m+1, size_of_dop(CONE_SEGMENT))) {
			delete[] map; return(0);
		}
		for ( k = 0; k <= m; k++) map[k] = mp[k];
		mp[  m] = (CMap)CYL_SEGMENT;
		mp[++m] = (CMap)ce[m2]->mp[8]*2.;
		mp[++m] = (CMap)A;

//////////////////////////////////////////////
//...корректируем локальную систему координат;
		map_normalizat(ce[m1]->mp, NULL, ort1);
		map_normalizat(ce[m2]->mp, NULL, ort2); map[1] = map[2] = map[3] = 0.;

		make_local(ort1,   map);
		make_local(ort2,   map);
		make_local(ort2+6, map); set_point(map+1, mp+1); f = arg2(comp(ort2[6], ort2[7]));
		mp[1] = mp[2] = mp[3] =
		mp[4] = mp[5] = mp[6] = 0.;

		if (ce[m4]->topo_id(ce[m1]->graph[3])) {
			set_point(P, ce[ ce[m1]->graph[2]]->mp+1);
			ort1[2] = -ort1[2];
		}
		else set_point(P, ce[ce[m1]->graph[3]]->mp+1);
		make_local    (P, mp);

		k     = fabs(1.+ort1[2]*ort2[2]) < EE_dop; m = mp[7] < 0.; if (m) k = 1-k;
		P[0]  =
		P[1]  = 0.;
		P[2] -= A*(k ? .5 : -.5); make_common(P, map); set_point(map+1, P);

		map_iso   (mp, NULL, 0., k ? 0. : M_PI, k ? f : M_PI-f);
		map_common(mp, map);

		delete[] map;
	}
#endif
	return(1);
}

////////////////////////////////////////////////
//...идентификация прямоугольного сегмента тора;
int CCells::id_torus_segment()
{
	int /*m1, m2, m3, m4, k, */m;
	if (mp && mp[0] == ID_MAP(2, TORUS_GENUS) && graph)
	if (graph[1] == 0 && add_new_maps(mp, (m = size_of_map(2, TORUS_GENUS))+1, size_of_dop(TORUS_SEGMENT))) { //...добавляем параметры;
		mp[  m] = (CMap)TORUS_SEGMENT;
		mp[++m] = (CMap)M_PI*2.;
		mp[++m] = (CMap)0.;
		mp[++m] = (CMap)M_PI*2.;
	}
  return(1);
}

/////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной уголком;
int CCells::id_ugolok_cell()
{
  return(1);
}

/////////////////////////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной треугольником с одной дугой;
int CCells::id_curve1_tria()
{
	int m1, m2, m3, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 3 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE1_TRIA_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE1_TRIA_CELL;
	}
	return(1);
}

//////////////////////////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной треугольником с двумя дугами;
int CCells::id_curve2_tria()
{
	int m1, m2, m3, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 3 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE2_TRIA_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE2_TRIA_CELL)))
			 mp[m] = (CMap)CURVE2_TRIA_CELL;
	}
	return(1);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной четырехугольником с одной дугой;
int CCells::id_curve1_quad()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[5]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[2]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[5]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[3]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE1_QUAD_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE1_QUAD_CELL;
	}
	return(1);
}

///////////////////////////////////////////////////////////////////
//...идентификация четырехугольника с двумя противолежащими дугами;
int CCells::id_curve2_quad()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[5]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[2]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[5]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[3]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE2_QUAD_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE2_QUAD_CELL;
	}
	return(1);
}

/////////////////////////////////////////////////////////////
//...идентификация четырехугольника с двумя соседними дугами;
int CCells::id_curve11quad()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[5]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[2]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[5]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[3]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE11QUAD_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE11QUAD_CELL;
	}
	return(1);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной четырехугольником с тремя дугами;
int CCells::id_curve3_quad()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[3]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[4]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[5]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[2]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[4]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[5]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[2]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[3]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS) ||
      ce[m1  = graph[5]]->graph[1] == 2 && ce[m1]->mp[0] == ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[2]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[3]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[4]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE3_QUAD_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE3_QUAD_CELL;
	}
	return(1);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...идентификация пространственной ячейки, ограниченной полностью кривым четырехугольником;
int CCells::id_curve4_quad()
{
	int m1, m2, m3, m4, m;
	if (mp   && mp[0]    == ID_MAP(2, NULL_GENUS) &&
      graph && graph[1] == 4 && (
      ce[m1  = graph[2]]->graph[1] == 2 && ce[m1]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m2  = graph[3]]->graph[1] == 2 && ce[m2]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m3  = graph[4]]->graph[1] == 2 && ce[m3]->mp[0] != ID_MAP(1,   NULL_GENUS) &&
      ce[m4  = graph[5]]->graph[1] == 2 && ce[m4]->mp[0] != ID_MAP(1,   NULL_GENUS))) {

/////////////////////////
//...добавляем параметры;
		if (mp[m  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) 
			 mp[m] = (CMap)CURVE4_QUAD_CELL; else
		if (add_new_maps(mp, m+1, size_of_dop(CURVE1_TRIA_CELL)))
			 mp[m] = (CMap)CURVE4_QUAD_CELL;
	}
	return(1);
}

///////////////////////////////////
//...идентификация размерных ячеек;
int CCells::segms_id()
{
  int k, m;
  if (2 == (k = cells_dim())) {
     if (mp[m = size_of_map(k, map_genus(mp))] == (CMap)NULL_CELL/* || 
		   mp[m] == (CMap)FACET_CELL || mp[m] == (CMap)SPH_SEGMENT*/) {
         if (! id_sheet         () || //...определяется;
             ! id_ring_segment  () || //...определяется;
             ! id_cyl_segment   () || //...определяется;
             ! id_sph_segment   () || //...определяется частично;
             ! id_spr_segment   () || //...определяется частично;
             ! id_cone_segment  () || //...определяется частично;
             ! id_torus_segment () || //...не определяется;
             ! id_ugolok_cell   () || //...не определяется;
             ! id_curve1_tria   () || //...определяется;
             ! id_curve2_tria   () || //...определяется;
             ! id_curve1_quad   () || //...определяется;
             ! id_curve2_quad   () || //...определяется;
             ! id_curve11quad   () || //...определяется;
             ! id_curve3_quad   () || //...определяется;
             ! id_curve4_quad   ()    //...определяется;
				) return(0);
     }
  }
  return(1);
}

///////////////////////////////////////////////////
//...идентификация размерных ячеек для поверхности;
int CCells::segms_bar()
{
  if (! mp) {
     for (int k = 0; k < graph[0]; k++)
         if (! ce[k]->segms_id()) return(0);
  }
  return(1);
}

/////////////////////////////////////////////////////
//...auxilliary function for grid nodes of line cell;
int CCells::grid_line(CGrid * nd, double h, int Max_N, int hit)
{
	if (nd && mp && ID_MAP(1, NULL_GENUS) == mp[0]) {
		int * gm, N0, j, l = nd->N,
				m1 = get_num(graph, 0),
				m2 = get_num(graph, 1);
		double * X0, * Y0, * Z0, * nX0, * nY0, * nZ0, fx = ce[m2]->mp[1]-ce[m1]->mp[1],
				CY = cos(mp[5]), SY = sin(mp[5]), fy = ce[m2]->mp[2]-ce[m1]->mp[2],
				CX = cos(mp[4]), SX = sin(mp[4]), fz = ce[m2]->mp[3]-ce[m1]->mp[3], f;

		if (! nd->geom) nd->init_geom();
		if (Max_N <= 0) Max_N = 1;

///////////////////////////////////////
//...определяем количество узлов сетки;
		if ((f = sqrt(fx*fx+fy*fy+fz*fz)) < EE_ker) f = 0.;
		N0 = h > 0. ? 2+(int)ceil(f/h) : 1+abs(Max_N);
		if (f == 0.) return(1);

///////////////////////////////////////////////////////////////////////
//...распределяем массивы точек для всего контура и описание геометрии;
		X0  = (double *)new_struct (N0*sizeof(double));
		Y0  = (double *)new_struct (N0*sizeof(double));
		Z0  = (double *)new_struct (N0*sizeof(double));
		nX0 = (double *)new_struct (N0*sizeof(double));
		nY0 = (double *)new_struct (N0*sizeof(double));
		nZ0 = (double *)new_struct (N0*sizeof(double));
		gm  = (int    *)new_struct((N0+3)*sizeof(int));
		if (! X0 || ! Y0 || ! Z0 || ! nX0 || ! nY0 || ! nZ0 || ! gm || m1 < 0 || m2 < 0) {
			delete_struct(X0); delete_struct(nX0);
			delete_struct(Y0); delete_struct(nY0);
			delete_struct(Z0); delete_struct(nZ0);
			delete_struct(gm); return(0);
		}

///////////////////////////////////////////////////////
//...инициализиpуем узлы сетки и касательную на прямой;

/*		double ort[3] = {0., 1., 0.};
		point_iso(ort, NULL, mp[4], mp[5], mp[6]);*/

		X0[0] = ce[m1]->mp[1]; fx /= N0-1; nX0[0] = CX*SY; /*nX0[0] = ort[0];*/ X0[N0-1] = ce[m2]->mp[1]; nX0[N0-1] = nX0[0];
		Y0[0] = ce[m1]->mp[2]; fy /= N0-1; nY0[0] = SX*SY; /*nY0[0] = ort[1];*/ Y0[N0-1] = ce[m2]->mp[2]; nY0[N0-1] = nY0[0];
		Z0[0] = ce[m1]->mp[3]; fz /= N0-1; nZ0[0] = CY;    /*nZ0[0] = ort[2];*/ Z0[N0-1] = ce[m2]->mp[3]; nZ0[N0-1] = nZ0[0];
		gm[0]  = 1;  gm[1] = GL_LINE_STRIP;
		gm[2]  = N0; gm[3] = l; gm[2+N0] = N0-1+l;
		for (j = 1; j < N0-1; j++) {
			X0[j] = X0[j-1]+fx; nX0[j] = nX0[0];
			Y0[j] = Y0[j-1]+fy; nY0[j] = nY0[0];
			Z0[j] = Z0[j-1]+fz; nZ0[j] = nZ0[0]; gm[3+j] = j+l;
		}

///////////////////////////////////////////////////////////
//...добавляем построенные элементы сетки к общему формату;
		if (nd->grid_add (X0, Y0, Z0, nX0, nY0, nZ0, N0, gm, hit)) return(1);
		else {
			delete_struct(X0); delete_struct(nX0);
			delete_struct(Y0); delete_struct(nY0);
			delete_struct(Z0); delete_struct(nZ0);
			delete_struct(gm); nd->zero_grid(); return(0);
		}
	}
	return(0);
}

///////////////////////////////////////////////////////
//...auxilliary function for grid nodes of circle cell;
int CCells::grid_circ(CGrid * nd, double h, int Max_N, int hit)
{
	if (nd && mp && ID_MAP(1, SPHERE_GENUS) == mp[0]) {
		int * gm, N0, j, l = nd->N,
			  m1 = get_num(graph, 0),
			  m2 = get_num(graph, 1);
		double * X0, * Y0, * Z0, * nX0, * nY0, * nZ0, R = fabs(mp[7]), fi = 2.*fabs(mp[8]), f,
			  CZ = cos(mp[4]), SZ = sin(mp[4]),
			  CY = cos(mp[5]), SY = sin(mp[5]),
			  CX = cos(mp[6]), SX = sin(mp[6]), P[3];

		if (! nd->geom) nd->init_geom();
		if (Max_N <= 0) Max_N = 1;

///////////////////////////////////////
//...определяем количество узлов сетки;
		if ((f = fabs(R*fi)) < EE_ker) f = 0.;
		N0 = h > 0. ? 2+(int)ceil(f/h) : 1+abs(Max_N);
		if (f == 0.) return(1);
		f = fabs(mp[8]);

///////////////////////////////////////////////////////////////////////
//...распределяем массивы точек для всего контура и описание геометрии;
		X0  = (double *)new_struct (N0*sizeof(double));
		Y0  = (double *)new_struct (N0*sizeof(double));
		Z0  = (double *)new_struct (N0*sizeof(double));
		nX0 = (double *)new_struct (N0*sizeof(double));
		nY0 = (double *)new_struct (N0*sizeof(double));
		nZ0 = (double *)new_struct (N0*sizeof(double));
		gm  = (int    *)new_struct((N0+3)*sizeof(int));
		if (! X0 || ! Y0 || ! Z0 || ! nX0 || ! nY0 || ! nZ0 || ! gm) {
			delete_struct(X0); delete_struct(nX0);
			delete_struct(Y0); delete_struct(nY0);
			delete_struct(Z0); delete_struct(nZ0);
			delete_struct(gm); return(0);
		}

////////////////////////////////////////////////////////////////
//...инициализиpуем узлы сетки и касательную на дуге окружности;
		X0[0] =  R*(nY0[0] = cos(f)); nY0[N0-1] =  nY0[0];
		Y0[0] = -R*(nX0[0] = sin(f)); nX0[N0-1] = -nX0[0]; fi /= N0-1;

		gm[0]  = 1;  gm[1] = GL_LINE_STRIP;
		gm[2]  = N0; gm[3] = l; gm[2+N0] = N0-1+l;
		for (f = fi-f, j = 1; j < N0; j++, f += fi) {
			X0[j] =  R*(nY0[j] =  cos(f));
			Y0[j] = -R*(nX0[j] = -sin(f)); gm[3+j] = j+l;
		}

///////////////////////////////////////////////////////////
//...добавляем построенные элементы сетки к общему формату;
		if (! nd->grid_add (X0, Y0, Z0, nX0, nY0, nZ0, N0, gm, hit)) {
			delete_struct(X0); delete_struct(nX0);
			delete_struct(Y0); delete_struct(nY0);
			delete_struct(Z0); delete_struct(nZ0);
			delete_struct(gm); nd->zero_grid(); return(0);
		}

/////////////////////////////////////////////////////////////////////////
//...коррекция начальной и конечной точки, чтобы не терять значащих цифр;
		int id_correct = 0;//...работает  с ошибками!!!
		if (id_correct && 0 <= m1 && 0 <= m2) { 
			P[0] = nd->X[l]; P[1] = nd->Y[l];
			P[2] = nd->Z[l]; if (abs_point(P, ce[m1]->mp+1) >= EE_dop) swap(m1, m2);

			nd->X[l] = ce[m1]->mp[1]; nd->X[nd->N-1] = ce[m2]->mp[1];
			nd->Y[l] = ce[m1]->mp[2]; nd->Y[nd->N-1] = ce[m2]->mp[2];
			nd->Z[l] = ce[m1]->mp[3]; nd->Z[nd->N-1] = ce[m2]->mp[3];
		}

///////////////////////////////////////////////////////////////
//...pазвоpачиваем сетку в пpостpанстве и выходим из программы;
		nd->grid_iso(l, nd->N, mp+1, CZ, SZ, CY, SY, CX, SX);
		return(1);
	}
	return(0);
}

//////////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция, стpоящая геометpическую сетку на эллиптической дуге;
int CCells::grid_ellipt(CGrid * nd, double h, int Max_N, int hit)
{
  int k;
  if (nd && mp && ID_MAP(1, CYL_GENUS) == mp[0] && mp[k = size_of_map(1, CYL_GENUS)] == (CMap)ELLIPT_ARC) {
     int    * gm, N0, j, l = nd->N,
              m1 = get_num(graph, 0),
              m2 = get_num(graph, 1);
     double * X0, * Y0, * Z0, * nX0, * nY0, * nZ0, A = fabs(mp[7]), B = fabs(mp[8]),
              C  = sqrt(fabs(A*A-B*B)), f, f0, f1, fi, R, ee = C/A, ff = B*B,
              CZ = cos(mp[4]), SZ = sin(mp[4]),
              CY = cos(mp[5]), SY = sin(mp[5]),
              CX = cos(mp[6]), SX = sin(mp[6]), P[3];

     if (! nd->geom) nd->init_geom();
     if (Max_N <= 0) Max_N = 1;

////////////////////////////////////////
//...извлекаем дополнительные параметры;
     f0 = mp[++k];
     f1 = mp[++k];

///////////////////////////////////////
//...определяем количество узлов сетки;
     if ((f = fabs(.5*(A+B)*(fi = f1-f0))) < EE_ker) f = 0.;
     N0 = h > 0. ? 2+(int)ceil(f/h) : 1+abs(Max_N);
     if (f == 0.) return(1);
     f = f0;

///////////////////////////////////////////////////////////////////////
//...распределяем массивы точек для всего контура и описание геометрии;
     X0  = (double *)new_struct (N0*sizeof(double));
     Y0  = (double *)new_struct (N0*sizeof(double));
     Z0  = (double *)new_struct (N0*sizeof(double));
     nX0 = (double *)new_struct (N0*sizeof(double));
     nY0 = (double *)new_struct (N0*sizeof(double));
     nZ0 = (double *)new_struct (N0*sizeof(double));
     gm  = (int    *)new_struct((N0+3)*sizeof(int));
     if (! X0 || ! Y0 || ! Z0 || ! nX0 || ! nY0 || ! nZ0 || ! gm) {
         delete_struct(X0); delete_struct(nX0);
         delete_struct(Y0); delete_struct(nY0);
         delete_struct(Z0); delete_struct(nZ0);
         delete_struct(gm); return(0);
     }

////////////////////////////////////////////////////////////////
//...инициализиpуем узлы сетки и касательную на дуге окружности;
     X0[0] = (R = ff/(A+C*(nY0[0] = cos(f))))*nY0[0]+C; nY0[0] += ee; fi /= N0-1;
     Y0[0] = -R*(nX0[0] = -sin(f)); nX0[0] *= (R = 1./sqrt(nX0[0]*nX0[0]+nY0[0]*nY0[0])); nY0[0] *= R;

     gm[0]   = 1;  gm[1] = GL_LINE_STRIP;
     gm[2]   = N0; gm[3] = l; gm[2+N0] = N0-1+l;
     for (f += fi, j = 1; j < N0; j++, f += fi) {
          X0[j] = (R = ff/(A+C*(nY0[j] = cos(f))))*nY0[j]+C; nY0[j] += ee; gm[3+j] = j+l;
          Y0[j] = -R*(nX0[j] = -sin(f)); nX0[j] *= (R = 1./sqrt(nX0[j]*nX0[j]+nY0[j]*nY0[j])); nY0[j] *= R;
     }

///////////////////////////////////////////////////////////
//...добавляем построенные элементы сетки к общему формату;
     if (! nd->grid_add (X0, Y0, Z0, nX0, nY0, nZ0, N0, gm, hit)) {
         delete_struct(X0); delete_struct(nX0);
         delete_struct(Y0); delete_struct(nY0);
         delete_struct(Z0); delete_struct(nZ0);
         delete_struct(gm); nd->zero_grid(); return(0);
     }

/////////////////////////////////////////////////////////////////////////
//...коррекция начальной и конечной точки, чтобы не терять значащих цифр;
		int id_correct = 0;
		if (id_correct && 0 <= m1 && 0 <= m2) {  
			P[0] = nd->X[l]; P[1] = nd->Y[l];
			P[2] = nd->Z[l]; if (abs_point(P, ce[m1]->mp+1) >= EE_dop) swap(m1, m2);

			nd->X[l] = ce[m1]->mp[1]; nd->X[nd->N-1] = ce[m2]->mp[1];
			nd->Y[l] = ce[m1]->mp[2]; nd->Y[nd->N-1] = ce[m2]->mp[2];
			nd->Z[l] = ce[m1]->mp[3]; nd->Z[nd->N-1] = ce[m2]->mp[3];
		}

///////////////////////////////////////////////////////////////
//...pазвоpачиваем сетку в пpостpанстве и выходим из программы;
		nd->grid_iso(l, nd->N, mp+1, CZ, SZ, CY, SY, CX, SX);
		return(1);
  }
  return(0);
}

/////////////////////////////////////////
//...grid nodes for one dimensional cell;
int CCells::grid_cells1(CGrid * nd, double h, int Max_knee)
{
  if (1 == cells_dim())
  if (grid_line  (nd, h, Max_knee) ||
      grid_circ  (nd, h, Max_knee) ||
      grid_ellipt(nd, h, Max_knee)) 
	return(1);
	return(0);
}

/////////////////////////////////////////////////////////
//...grid nodes for all elements of skeleton of the cell;
int CCells::grid_skeleton(CGrid * nd, double h, int Max_knee, int id_full, Topo * link)
{
  if (nd) {
     int  k, j, l;
     for (k = 0; k < graph[0]; k++) {
          for (j = l = 0; link && l < k; l++)
				 if ((j = ::topo_id(ce[l]->graph, k)) != 0) break;

          if (! link || j && l < link[0] && link[l+1] >= 0)
          if (1 == ce[k]->cells_dim()) {
              ce[k]->grid_line  (nd, h, id_full ? Max_knee : 1, -k-1);
              ce[k]->grid_circ  (nd, h, Max_knee, -k-1);
              ce[k]->grid_ellipt(nd, h, Max_knee, -k-1);
          }
     }
  }
  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////
//...проверка попадания одной точки (лежащей в плоскости фасеты) внутрь плоской фасеты;
int CCells::in_plane_facet(double * P, double * pm, int N_arc, int k1, int k2)
{
  double f, X = P[k1], Y = P[k2];
  int    k, i = 0, l = 0, m1 = 0, m3 = X < pm[k1],
            j = 0, m = 0, m2 = 0, m4 = Y < pm[k2]; swap(m1, m3); swap(m2, m4);

//////////////////////////////////////////////////////////////////////////////////
//...определяем число пересечений контура с лучом по вещественной и по мнимой оси;
  for (k = 1; k < N_arc; k++) {
      m3 = X < pm[3*k+k1]; swap(m1, m3);
      m4 = Y < pm[3*k+k2]; swap(m2, m4);

      if (m1 != m3 && m2 != m4) {
          f  = (pm[3*(k-1)+k1]-X)*(pm[3*k+k2]-pm[3*(k-1)+k2])+
               (Y-pm[3*(k-1)+k2])*(pm[3*k+k1]-pm[3*(k-1)+k1]);
          i += (f > 0. && m2 == 1) || (f < 0. && m2 == 0);
          l += (f > 0. && m2 == 0) || (f < 0. && m2 == 1);
          j += (f > 0. && m1 == 0) || (f < 0. && m1 == 1);
          m += (f > 0. && m1 == 1) || (f < 0. && m1 == 0);
      }
      else {
          if (m2 != m4) {
              i += m1 == 1;
              l += m1 == 0;
          }
          if (m1 != m3) {
              j += m2 == 1;
              m += m2 == 0;
          }
      }
  }
  return((i % 2) || (j % 2) || (l % 2) || (m % 2));
}

/////////////////////////////////////////////////////////////
//...проверка попадания проекции точки внутрь плоской фасеты;
int CCells::in_plane_facet(double * P, double *& pm, int & N_ini, int * mm, int id_fast)
{
  if (mp &&   ID_MAP(2, NULL_GENUS)  == mp[0] && (
      mp[size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL)) {
     int N_arc, arc, prev = graph[N_arc = graph[1]+1], k, i, m1, m2, k1, k2, k3;
     double * pm_new, nX, nY, nZ, f, d[3], p[3];

/////////////////////
//...проверяем буфер;
     if (N_ini < N_arc) {
         if (NULL != (pm_new = (double *)new_struct(3*N_arc*sizeof(double)))) {
             delete_struct(pm); pm = pm_new; N_ini = N_arc;
         } else return(0);
     }

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
     for (k = 0, i = 2; i <= N_arc; i++, prev = arc) {
			arc = graph[i];
			if (id_fast == OK_STATE) m1 = arc; else {
				m1  = get_num(ce[arc]->graph, 0),
				m2  = get_num(ce[arc]->graph, 1);
				if (! ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			pm[k++] = ce[m1]->mp[1];
			pm[k++] = ce[m1]->mp[2];
			pm[k++] = ce[m1]->mp[3];
     }
     pm[k++] = pm[0];
     pm[k++] = pm[1];
     pm[k++] = pm[2];

///////////////////////////////////////////////////////
//...определяем наиболее подходящую компоненту нормали;
	  if (id_fast == NULL_STATE) {
		  nX = cos(mp[4])*(f = sin(mp[5]));
		  nY = sin(mp[4])* f;
		  nZ = cos(mp[5]);
	  }
	  else {
		  nX = mp[4];
		  nY = mp[5];
		  nZ = mp[6];
	  }
     if (fabs(nY) > fabs(nZ)) {
         k = 1; f = fabs(nY);
     }
     else {
         k = 2; f = fabs(nZ);
     }
     if (fabs(nX) > f) k = 0; k1 = (k+1)%3; k2 = (k+2)%3; k3 = k;

/////////////////////////////////////////////////////////////////////////////////////////
//...определяем попадание проекции точки внутрь проекции плоского многоугольного контура;
     if (fabs(f = (mp[1]-P[0])*nX+(mp[2]-P[1])*nY+(mp[3]-P[2])*nZ) < EE_ker) {
         d[0] =
         d[1] =
         d[2] = 0.;
     }
     else {
         d[0] = fabs(nX) < EE_ker ? MAX_HIT : f/nX;
         d[1] = fabs(nY) < EE_ker ? MAX_HIT : f/nY;
         d[2] = fabs(nZ) < EE_ker ? MAX_HIT : f/nZ;
     }

////////////////////
//...направление k3;
     if (d[k3] !=  MAX_HIT) {
         ::set_point(p, P); i = in_plane_facet(p, pm, N_arc, k1, k2);
         if (d[k3] <= 0.) mm[2*k3]   += i;
         if (d[k3] >= 0.) mm[2*k3+1] += i;
     }

////////////////////
//...направление k1;
     if (d[k1] == 0.) {
         mm[2*k1]   += i;
         mm[2*k1+1] += i;
     }
     else if (d[k1] != MAX_HIT) {
         ::set_point(p, P); p[k1] += d[k1]; k = in_plane_facet(p, pm, N_arc, k1, k2);
          if (d[k1] < 0.)  mm[2*k1  ] += k;
          else             mm[2*k1+1] += k;
     }

////////////////////
//...направление k2;
     if (d[k2] == 0.) {
         mm[2*k2]   += i;
         mm[2*k2+1] += i;
     }
     else if (d[k2] != MAX_HIT) {
         ::set_point(p, P); p[k2] += d[k2]; k = in_plane_facet(p, pm, N_arc, k1, k2);
          if (d[k2] < 0.)  mm[2*k2  ] += k;
          else             mm[2*k2+1] += k;
     }
  }
  return(1);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...проверка попадания точки в пространственный размерный обpазец, заданный форматом фасет;
int CCells::in_bar_facet(double X, double Y, double Z, int id_fast)
{
  int		mm[6] = {0, 0, 0, 0, 0, 0}, k, N_ini = 6;
  double *	pm = (double *)new_struct(3*N_ini*sizeof(double)), P[3] = {X, Y, Z};
  if (pm && graph && ! mp) {
      for (k = 0; k < graph[0]; k++)
          ce[k]->in_plane_facet(P, pm, N_ini, mm, id_fast);
  }
  else
  if (pm && graph && mp) {
      for (k = 0; k < graph[1]; k++)
          ce[graph[k+2]]->in_plane_facet(P, pm, N_ini, mm, id_fast);
  }
  delete_struct(pm);
  return((mm[0] % 2) || (mm[1] % 2) || (mm[2] % 2) || (mm[3] % 2) || (mm[4] % 2) || (mm[5] % 2));
}

////////////////////////////////////////////////////////
//...проверка принадлежности точки нижней полуплоскости;
int CCells::in_half_plane(double X, double Y, double Z, int id_local, double eps)
{
  if (mp && ID_MAP(2, NULL_GENUS) == mp[0]) {
      double p[3] = { X, Y, Z };
      if (id_local) make_local(p, mp); 
      else          p[2] = (p[0]-mp[1])*mp[4]+(p[1]-mp[2])*mp[5]+(p[2]-mp[3])*mp[6];
      return       (p[2] < eps);
  }
  return(1);
}

///////////////////////////////////////////////////////////////
//...проверка попадания точки внутрь (выпуклого) многогранника;
int CCells::in_poly_body(double X, double Y, double Z, int id_local, int * mm, double eps)
{
  int k, m = 1;
  if (graph && ! mp) {
      for (k = 0;  m && k < graph[0]; k++)
      if  (! mm ||  k != mm[0] && k != mm[1]) m &= ce[k]->in_half_plane(X, Y, Z, id_local, eps);
      if  (m && mm) m &= ce[mm[0]]->in_half_plane(X, Y, Z, id_local, eps) ||
                         ce[mm[1]]->in_half_plane(X, Y, Z, id_local, eps);
  }
  else 
  if (graph && mp) {
      for (k = 0;  m && k < graph[1]; k++)
      if  (! mm ||  k != mm[0] && k != mm[1]) m &= ce[graph[k+2]]->in_half_plane(X, Y, Z, id_local, eps);
      if  (m && mm) m &= ce[graph[mm[0]+2]]->in_half_plane(X, Y, Z, id_local, eps) ||
                         ce[graph[mm[1]+2]]->in_half_plane(X, Y, Z, id_local, eps);
  }
  return(m);
}

////////////////////////////////////////////////////
//...проверка попадания точки внутрь плоской фасеты;
int CCells::in_poly_facet(double X, double Y, double Z, int id_fast, double eps)
{
	int l = 0;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) {
		double nX, nY, nZ, px, py, pz, P[3], f; 

///////////////////////////////////////////////////////////////////////
//...определяем плоскость нормали и проверяем принадлежность плоскости;
		if (id_fast == NULL_STATE) {
			nX = cos(mp[4])*(nZ = sin(mp[5]));
			nY = sin(mp[4])* nZ;
			nZ = cos(mp[5]);
		}
		else {
			nX = mp[4]; 
			nY = mp[5]; 
			nZ = mp[6]; 
		}
		if (fabs(f = (mp[1]-X)*nX+(mp[2]-Y)*nY+(mp[3]-Z)*nZ) < eps) l = 1;

/////////////////////////////////////////
//...пробегаем все звенья многоугольника;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], m1, m2;
      for (int k = 1; l && k <= N_arc; k++, prev = arc) {
			arc = graph[k+1];
			if (id_fast == OK_STATE) {
				m1  = prev;
				m2  = arc;
			}
			else {
				m1  = get_num(ce[arc]->graph, 0),
				m2  = get_num(ce[arc]->graph, 1);
				if (! ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			px = ce[m2]->mp[1]-ce[m1]->mp[1];
			py = ce[m2]->mp[2]-ce[m1]->mp[2];
			pz = ce[m2]->mp[3]-ce[m1]->mp[3];
			P[0] = X-ce[m1]->mp[1];
			P[1] = Y-ce[m1]->mp[2];
			P[2] = Z-ce[m1]->mp[3];

			l &= ((f = (P[2]*py-P[1]*pz)*nX+(P[0]*pz-P[2]*px)*nY+(P[1]*px-P[0]*py)*nZ+eps) > 0);
		}
	}
	return(l);
}

//////////////////////////////////////////////////////////////
//...проверка попадания точки внутрь кругового многоугольника;
int  CCells::in_bar(double X, double Y)
{
	if (! mp) {
		complex C = comp(X, Y), z0, z1, z2;
		double  f = 0.;

////////////////////////////////////////////////////////////////
//...определяем число полных обходов при движении вдоль границы;
		for (int k = 0; k < graph[0]; k++)
		if (ce[k]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
			z0 = ce[k]->get_arc_center();
			z1 = ce[k]->get_arc_point(0);
			z2 = ce[k]->get_arc_point(1);

			f += arc_delta(z0-C, z1-C, z2-C, ce[k]->mp[5] == M_PI);
		}
		else
		if (ce[k]->mp[0] == ID_MAP(1, NULL_GENUS)) {
			z1 = ce[k]->get_arc_point(0);
			z2 = ce[k]->get_arc_point(1);

			f += arc_delta(z1-C, z1-C, z2-C, 0);
		}
		return(fabs(f) > M_PI_2);
  }
  return(0);
}

////////////////////////////////////////////////////////////////////////
//...проверка попадания точки внутрь плоской, почти произвольной ячейки;
int  CCells::in_bar_cell(double X, double Y)
{
	if (mp && mp[0] == ID_MAP(2, NULL_GENUS)) {
		complex C = comp(X, Y), z0, z1, z2;
		double  f = 0.;

////////////////////////////////////////////////////////////////
//...определяем число полных обходов при движении вдоль границы;
		for (int k, i = 0; i < graph[1]; i++)
		if (ce[k = graph[i+2]]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
			z0 = ce[k]->get_arc_center();
			z1 = ce[k]->get_arc_point(0);
			z2 = ce[k]->get_arc_point(1);

			f += arc_delta(z0-C, z1-C, z2-C, ce[k]->mp[5] == M_PI);
		}
		else
		if (ce[k]->mp[0] == ID_MAP(1, NULL_GENUS)) {
			z1 = ce[k]->get_arc_point(0);
			z2 = ce[k]->get_arc_point(1);

			f += arc_delta(z1-C, z1-C, z2-C, 0);
		}
		return(fabs(f) > M_PI_2);
  }
  return(0);
}

/////////////////////////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур для плоской фасеты с одной кривой дугой;
void CCells::curve1_tria_QG(CGrid * nd, int N_elem)
{
	int l;
	if (nd && ID_MAP(2, NULL_GENUS) == mp[0] && mp[size_of_map(2, NULL_GENUS)] == (CMap)CURVE1_TRIA_CELL &&
		 graph && graph[1] == 3) {
		double pm[9];

/////////////////////////////////////////
//...определяем положение кривой границы;
		for (l = 1; l <= graph[1]; l++)
			if ( ce[graph[l+1]]->mp[0] != ID_MAP(1, NULL_GENUS)) break;
		if (l > graph[1]) l = 2;

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], m1, m2, j;
		for (j = 1; j <= N_arc; j++, prev = arc) {
			arc = graph[j+1];
			m1  = get_num(ce[arc]->graph, 0),
			m2  = get_num(ce[arc]->graph, 1);
			if (! ce[prev]->topo_id(m1)) swap(m1, m2);

			pm[3*((j-l+N_arc+1)%N_arc)]   = ce[m1]->mp[1];
			pm[3*((j-l+N_arc+1)%N_arc)+1] = ce[m1]->mp[2];
			pm[3*((j-l+N_arc+1)%N_arc)+2] = ce[m1]->mp[3];
		}
		nd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
		nd->QG_tria_curve(ce[graph[l+1]]->mp, pm);
	}
}

///////////////////////////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур для плоской фасеты с двумя кривыми дугами;
void CCells::curve2_tria_QG(CGrid * nd, int N_elem)
{
	int l;
	if (nd && ID_MAP(2, NULL_GENUS) == mp[0] && mp[size_of_map(2, NULL_GENUS)] == (CMap)CURVE2_TRIA_CELL &&
		 graph && graph[1] == 3) {
		double pm[12];

/////////////////////////////////////////
//...определяем положение кривой границы;
		for (l = 1; l <= graph[1]; l++)
			if ( ce[graph[l+1]]->mp[0] == ID_MAP(1, NULL_GENUS)) break;
		if (l > graph[1]) l = 2;

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], m1, m2, j, N_ini;
		for (j = 1; j <= N_arc; j++, prev = arc) {
			arc = graph[j+1];
			m1  = get_num(ce[arc]->graph, 0),
			m2  = get_num(ce[arc]->graph, 1);
			if (! ce[prev]->topo_id(m1)) swap(m1, m2);

			pm[3*((j-l+N_arc+1)%N_arc)+3] = ce[m1]->mp[1];
			pm[3*((j-l+N_arc+1)%N_arc)+4] = ce[m1]->mp[2];
			pm[3*((j-l+N_arc+1)%N_arc)+5] = ce[m1]->mp[3];
		}
		pm[0] = (pm[6]+pm[9] )*.5;
		pm[1] = (pm[7]+pm[10])*.5;
		pm[2] = (pm[8]+pm[11])*.5;

		nd->facet_QG(pm, (int)((N_elem+1)/1.4), NULL_STATE, NULL_STATE); N_ini = nd->N;
		nd->QG_tria_curve(ce[graph[(l+1)%N_arc+2]]->mp, pm);

		swap(pm[6], pm[9]);  swap(pm[3], pm[6]);
		swap(pm[7], pm[10]); swap(pm[4], pm[7]);
		swap(pm[8], pm[11]); swap(pm[5], pm[8]);

		nd->facet_QG(pm, (int)((N_elem+1)/1.4), NULL_STATE, NULL_STATE);
		nd->QG_tria_curve(ce[graph[l%N_arc+2]]->mp, pm, N_ini);
	}
}

/////////////////////////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур для плоской фасеты с одной кривой дугой;
void CCells::curve1_quad_QG(CGrid * nd, int N_elem)
{
	int l;
	if (nd && ID_MAP(2, NULL_GENUS) == mp[0] && mp[size_of_map(2, NULL_GENUS)] == (CMap)CURVE1_QUAD_CELL &&
		 graph && graph[1] == 4) {
		double pm[12];

/////////////////////////////////////////
//...определяем положение кривой границы;
		for (l = 1; l <= graph[1]; l++)
			if ( ce[graph[l+1]]->mp[0] != ID_MAP(1, NULL_GENUS)) break;
		if (l > graph[1]) l = 2;

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], m1, m2, j;
		for (j = 1; j <= N_arc; j++, prev = arc) {
			arc = graph[j+1];
			m1  = get_num(ce[arc]->graph, 0),
			m2  = get_num(ce[arc]->graph, 1);
			if (! ce[prev]->topo_id(m1)) swap(m1, m2);

			pm[3*((j-l+N_arc-1)%N_arc)]   = ce[m1]->mp[1];
			pm[3*((j-l+N_arc-1)%N_arc)+1] = ce[m1]->mp[2];
			pm[3*((j-l+N_arc-1)%N_arc)+2] = ce[m1]->mp[3];
		}
		swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
		swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);

		nd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
		nd->QG_quad_curve(ce[graph[l+1]]->mp, NULL, pm);
	}
}

///////////////////////////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур для плоской фасеты с двумя кривыми дугами;
void CCells::curve11quad_QG(CGrid * nd, int N_elem)
{
	int l;
	if (nd && ID_MAP(2, NULL_GENUS) == mp[0] && mp[size_of_map(2, NULL_GENUS)] == (CMap)CURVE11QUAD_CELL &&
		 graph && graph[1] == 4) {
		double pm[12];

/////////////////////////////////////////
//...определяем положение кривой границы;
		for (l = 1; l <= graph[1]; l++)
			if ( ce[graph[l+1]]->mp[0] != ID_MAP(1, NULL_GENUS)) break;
		if (l > graph[1]) l = 2;

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], m1, m2, j, N_ini;
		for (j = 1; j <= N_arc; j++, prev = arc) {
			arc = graph[j+1];
			m1  = get_num(ce[arc]->graph, 0),
			m2  = get_num(ce[arc]->graph, 1);
			if (! ce[prev]->topo_id(m1)) swap(m1, m2);

			pm[3*((j-l+N_arc+1)%N_arc)]   = ce[m1]->mp[1];
			pm[3*((j-l+N_arc+1)%N_arc)+1] = ce[m1]->mp[2];
			pm[3*((j-l+N_arc+1)%N_arc)+2] = ce[m1]->mp[3];
		}
		nd->facet_QG(pm, (int)((N_elem+1)/1.4), NULL_STATE, NULL_STATE); N_ini = nd->N;
		nd->QG_tria_curve(ce[graph[l+1]]->mp, pm);

		swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]); 

		nd->facet_QG(pm+3, (int)((N_elem+1)/1.4), NULL_STATE, NULL_STATE);
		nd->QG_tria_curve(ce[graph[l%N_arc+2]]->mp, pm+3, N_ini);
	}
}

////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур для плоской фасеты;
void CCells::facet_QG(CGrid * nd, int N_elem, int N_max, double * pp_cond, int id_fast)
{
	int k, l;
	if (nd && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) {
		double pm[12], P[12];
		int	 N_ini = 4;	if (N_max <= 0) N_max = 1;

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = graph[(N_arc = graph[1])+1], 
			 N_min = min(N_arc, N_ini), m1, m2, j;
		for (j = 1; j <= N_min; j++, prev = arc) {
			arc = graph[j+1];
			if (id_fast != OK_STATE) {
				m1  = get_num(ce[arc]->graph, 0),
				m2  = get_num(ce[arc]->graph, 1);
				if (! ce[prev]->topo_id(m1)) swap(m1, m2);

				pm[3*(j-1)]   = ce[m1]->mp[1];
				pm[3*(j-1)+1] = ce[m1]->mp[2];
				pm[3*(j-1)+2] = ce[m1]->mp[3];
			}
			else {
				pm[3*(j-1)]   = ce[arc]->mp[1];
				pm[3*(j-1)+1] = ce[arc]->mp[2];
				pm[3*(j-1)+2] = ce[arc]->mp[3];
			}
		}

/////////////////////////
//...tria and quad facet;
		if (N_arc == 3) {
			for (j = 1; j <= N_arc; j++)
			memcpy(P+3*(j-1), pm+3*((j+2)%N_arc), 3*sizeof(double));
			nd->grid_tria_3D_refine(N_max, P);
		}
		else
		if (N_arc == 4) {
			for (j = 1; j <= N_arc; j++)
			memcpy(P+3*(j-1), pm+3*((j+1)%N_arc), 3*sizeof(double));
			nd->grid_quad_3D_refine(N_max, N_max, P);
		}

/////////////////////////////////////
//...dividing any facet on triangles;
		else {
			memcpy(P, pm, 3*sizeof(double));
			for (j = 3; j <= N_min; j++) {
				memcpy(P+3, pm+3*((j+1)%3), 3*sizeof(double));
				memcpy(P+6, pm+3*((j+2)%3), 3*sizeof(double));

				nd->grid_tria_3D_refine(N_max, P);
			}
			memcpy(P, pm, 3*sizeof(double));
			memcpy(P+3, pm+3*(N_min-2), 6*sizeof(double));
			for (; j <= N_arc; j++, prev = arc) {
				arc = graph[j+1];
				m1  = get_num(ce[arc]->graph, 0),
				m2  = get_num(ce[arc]->graph, 1);
				if (! ce[prev]->topo_id(m1)) swap(m1, m2);

				P[3] = P[6]; P[6] = ce[m1]->mp[1];
				P[4] = P[7]; P[7] = ce[m1]->mp[2];
				P[5] = P[8]; P[8] = ce[m1]->mp[3];

				nd->grid_tria_3D_refine(N_max, P);
			}
		}
		if (pp_cond) {
			pp_cond[0] = mp[k+2];
			pp_cond[1] = mp[k+3];
			pp_cond[2] = mp[k+4];
			pp_cond[3] = mp[k+5];
		}

////////////////////////////////////////////////////
//...запоминаем данные и строим квадратуры на сетке;
		int * nd_geom = NULL;
		double * nd_X = NULL, * nd_Y = NULL, * nd_Z = NULL;

		swap(nd_X, nd->X);
		swap(nd_Y, nd->Y);
		swap(nd_Z, nd->Z);
		swap(nd_geom, nd->geom);
		nd->zero_grid();
		nd->add_params(1);

/////////////////////////////////
//...накапливаем граничные точки;
		if (nd_geom)
		for (k = 1, l = 0; l <  nd_geom[0]; l++, k += nd_geom[++k]+1) {
			int num = k, num_n = nd_geom[num+1], num_f = num_n+num, cnt = 0;

			if (nd_geom[num] == GL_TRIANGLES) {
				for (; num < num_f; num++) {
					P[cnt++] = nd_X[nd_geom[num+2]];
					P[cnt++] = nd_Y[nd_geom[num+2]];
					P[cnt++] = nd_Z[nd_geom[num+2]];
				}
				nd->facet_QG(P, N_elem, NULL_STATE, NULL_STATE);
			}
			else 
			if (nd_geom[num] == GL_QUAD_STRIP) {
				for (; num < num_f-3; num += 2) {
					P[0] = nd_X[nd_geom[num+2]];
					P[1] = nd_Y[nd_geom[num+2]];
					P[2] = nd_Z[nd_geom[num+2]];
					P[3] = nd_X[nd_geom[num+3]];
					P[4] = nd_Y[nd_geom[num+3]];
					P[5] = nd_Z[nd_geom[num+3]];
					P[6] = nd_X[nd_geom[num+4]];
					P[7] = nd_Y[nd_geom[num+4]];
					P[8] = nd_Z[nd_geom[num+4]];
					P[9] = nd_X[nd_geom[num+5]];
					P[10] = nd_Y[nd_geom[num+5]];
					P[11] = nd_Z[nd_geom[num+5]];

					nd->facet_QG(P, N_elem, OK_STATE, NULL_STATE);
				}
			}
		}
		for (l = 0; l < nd->N; l++) { //...правим нормали;
			nd->nX[l] = mp[4];
			nd->nY[l] = mp[5];
			nd->nZ[l] = mp[6];
		}
		delete_struct(nd_X);
		delete_struct(nd_Y);
		delete_struct(nd_Z);
		delete_struct(nd_geom);
	}
}

/////////////////////////////////////////////////////////
//...веса и узлы гауссовых квадратур размерных элементов;
void CCells::segms_QG(CGrid * nd, int N_elem, int N_max, double * pp_cond, int id_fast) 
{
	if (pp_cond) pp_cond[0] = pp_cond[1] = pp_cond[2] = pp_cond[3] = 0;
	if (2 == map_dim(mp)) {	//...некоторые усложненные квадратуры;
		curve1_tria_QG(nd, N_elem); //...сделано;
		curve2_tria_QG(nd, N_elem); //...сделано, но путем разбиения;
		curve1_quad_QG(nd, N_elem); //...сделано;
		curve11quad_QG(nd, N_elem); //...сделано путем разбиения;
		facet_QG		  (nd, N_elem, N_max, pp_cond, id_fast); //...сделано;
	}
	nd->segms_QG(mp, N_elem, N_max); //...набор квадратур для размерных сегментов;
}
