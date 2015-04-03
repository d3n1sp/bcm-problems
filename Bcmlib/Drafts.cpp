#include "stdafx.h"
#include "cdraft.h"
#include "clame2d.h"
#include "cheat2d.h"
#include "cheat3d.h"
#include "ccohes2d.h"
#include "ccohes3d.h"
#include "cmindl2d.h"
#include "cmindl3d.h"
#include "chydro3d.h"
#include "cvisco2d_grad.h"

//////////////////////////////////////////////////////////
//           CONSTRUCTION of the ALL DRAFTS             //
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////
//...global function for construction all drafts;
CDraft<double> * CreateDraft(Num_Draft id_DRAFT, int id_dop)
{
	switch (id_DRAFT) {
      case	  NULL_DRAFT: return new CDraft<double>;
      case	LAME2D_DRAFT: return new CLame2D(id_dop);
      case	LAME3D_DRAFT: return new CLame3D(id_dop);
      case	HEAT2D_DRAFT: return new CHeat2D(id_dop);
      case	HEAT3D_DRAFT: return new CHeat3D(id_dop);
      case COHES2D_DRAFT: return new CCohes2D(id_dop);
      case COHES3D_DRAFT: return new CCohes3D(id_dop);
      case MINDL2D_DRAFT: return new CMindl2D(id_dop);
      case MINDL3D_DRAFT: return new CMindl3D(id_dop);
      case HYDRO3D_DRAFT: return new CHydro3D(id_dop);
	}
   return NULL;
}

/////////////////////////////////
//...construction complex drafts;
CDraft<complex> * CreateDraftC(Num_Draft id_DRAFT, int id_dop)
{
	switch (id_DRAFT) {
      case			 NULL_DRAFT: return new CDraft<complex>;
//      case VISCO2D_GRAD_DRAFT: return new CVisco2D_grad(id_dop);
	}
   return NULL;
}

///////////////////////////////////////
//...construction double-double drafts;
CDraft<dd_real> * CreateDraftDD(Num_Draft id_DRAFT, int id_dop)
{
	switch (id_DRAFT) {
      case	  NULL_DRAFT: return new CDraft<dd_real>;
      //case HYDRO3D_DRAFT: return new CHydro3D(id_dop);
	}
   return NULL;
}

/////////////////////////////////////
//...construction quad-double drafts;
CDraft<qd_real> * CreateDraftQD(Num_Draft id_DRAFT, int id_dop)
{
	switch (id_DRAFT) {
      case	  NULL_DRAFT: return new CDraft<qd_real>;
      //case HYDRO3D_DRAFT: return new CHydro3D(id_dop);
	}
   return NULL;
}

/////////////////////////////////////////////////
//          КОРРЕКЦИЯ СЕТКИ ДЛЯ ABAQUS         //
/////////////////////////////////////////////////
/////////////////////////////////////
//...добавление новых узлов в модель;
void Inp_nodes_add(char * ch_NODES, CGrid * nd, int ID_node_set, int ID_element_set)
{
	const int STR_SIZE = 250;
	char * id_NODES = read_struct_ascii(ch_NODES);
	if (nd && id_NODES) {
		char  buff[2000], one_line[STR_SIZE+1], * pchar, temp = 0;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, end_element;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, "_add_nodes.inp");
		FILE * INP = fopen(buff, "w");
		if (INP) {
			int i, j, k, l, m = 0, NN = nd->N, elem;

/////////////////////////////////////////
//...отмечаем узлы с заданными условиями;
			if (nd->hit && nd->cond && nd->cond_ptr && ID_node_set < 0) {
				for (k = 0; k < nd->cond[0]; k++)
				if (nd->cond[l = nd->cond_ptr[k]] == (int)GL_NODES && nd->cond[l+2] == ID_node_set) {
					for (j = 2; j <= nd->cond[l+1]; j++) if (nd->cond[l+1+j] < nd->N) 
						if (nd->add_new_point(nd->X[nd->cond[l+1+j]], nd->Y[nd->cond[l+1+j]], nd->Z[nd->cond[l+1+j]])) 
							 nd->hit[nd->N-1] = nd->cond[l+1+j]+1; //...временная нумерация;
				}
				else
				if (nd->cond[l] == (int)GL_NODES_GENERATE && nd->cond[l+2] == ID_node_set) {
					for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) if (j < nd->N) 
						if (nd->add_new_point(nd->X[nd->cond[l+1+j]], nd->Y[nd->cond[l+1+j]], nd->Z[nd->cond[l+1+j]])) 
							 nd->hit[nd->N-1] = nd->cond[l+1+j]+1;
				}
			}	
		
////////////////////////////////////////
//...запись дополнительных узлов модели;
			FILE * INP_ADD = fopen("INP_ADD", "w");
			if (INP_ADD) {
				for (k = NN; k < nd->N; k++) //...дополнительные узлы;
				fprintf(INP_ADD, "%7i, %12g, %12g, %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);
				fclose (INP_ADD);
			}

///////////////////////////////
//...коррекция элеметов модели;
			if (nd->cond && nd->cond_ptr && ID_element_set < 0) {

////////////////////////////////////////////
//...ищем коррекируемое множество элементов;
				for (k = 0; ! m && k < nd->cond[0]; k++)
				if (nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS && nd->cond[l+2] == ID_element_set) m = 1;
				else
				if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE && nd->cond[l+2] == ID_element_set) m = 1;

////////////////////////////////////////////////
//...ищем корректируемые узлы в элеметах модели;
				if (nd->geom && nd->geom_ptr && m)
				for (k = NN; k < nd->N; k++) { //...дополнительные узлы;
					for (j = nd->cond[l+1]; j > 1; j--) if (nd->cond[l+1+j] < nd->geom[0]) //...корректируемые элементы; 
					for (m = nd->geom[(i = nd->geom_ptr[ nd->cond[l+1+j]])+1]; m > 2; m--) 
					if (nd->hit[k] == nd->geom[i+1+m]+1) nd->geom[i+1+m] = k; //...пишем во внутренней нумерации;
				}
				for (k = 0; k < nd->N; k++) nd->hit[k] = k+1; //...восстанавливаем нумерацию;
			}

//////////////////////////////////////////////////////////////////////////
//...запись скорректированных данных в первый *Part входного файла модели;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");

////////////////////////////////////////////////////////////////////////////////
//...ищем вхождение узлов в разделе *Part и печатаем все до начала узлов и узлы;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Node");
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE); swap(pchar[0], temp);
	
				fprintf(INP, "%s", id_NODES); swap(pchar[0], temp);
				fprintf(INP, "%s", one_line);
				for (k = 0; k < nd->N; k++)
				fprintf(INP, "%7i, %12g, %12g, %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);

////////////////////////////////////////////////////////////////////
//...ищем все вхождения элементов и печатаем их с внешними номерами;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_element = (upper_element -= (upper_element < upper ? 1 : 0));
					if (! ::strncmp(id_NODES+ppos_cur, "B31",   3) || //...,определяем тип элемента; 
						 ! ::strncmp(id_NODES+ppos_cur, "B32",	  3)) elem = GL_LINE_STRIP; else
					if (! ::strncmp(id_NODES+ppos_cur, "CPE3",  4) ||
						 ! ::strncmp(id_NODES+ppos_cur, "CPE6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS3",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "S3",	  2)) elem = GL_TRIANGLES; else
					if (! ::strncmp(id_NODES+ppos_cur, "CPE4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPE8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "S4R",	  3)) elem = GL_QUADS; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D10", 5)) elem = GL_TETRA; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D20", 5)) elem = GL_BOXS; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D15", 5)) elem = GL_PENTA; else elem = ERR_STATE;
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					fprintf(INP, "*Element, type=%s", one_line);

////////////////////////////////////////////
//...печатаем все элементы данной топологии;
					for (l = k = 0; k < nd->geom[0]; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+1]);					
							for (j = 2; j < nd->geom[l]; j++)
							fprintf(INP, ",%7i", nd->geom[l+1+j]+1);					
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				}

////////////////////////////////////////////////
//...печатаем оставшийся хвостик входного файла;
				fprintf(INP, "%s", id_NODES+end_element);
			}
		}
		fclose(INP);
	}
	delete_struct(id_NODES);
}

//////////////////////////////////////////////////////////////////////////////////
//...добавление новых элементов в модель (нумерация элементов и узлов -- внешняя);
void Inp_elements_add(char * ch_NODES, CGrid * nd, int * ID_elements, int ID_element_part)
{
	const int STR_SIZE = 250;
	char * id_NODES = read_struct_ascii(ch_NODES);
	if (nd && id_NODES) {
		char  buff[2000], one_line[STR_SIZE+1], * pchar, temp = 0;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, end_element, count_set, end_set, begin_set;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, "_add_elements.inp");
		FILE * INP = fopen(buff, "w");
		if (INP && ID_element_part > 0) {
			int i, j, j_max = 17, k, l, l0, m, set, elem, N_geom, NN_geom, N_set = 0, N_part = ID_element_part, 
				N_part_end, numeration_shift, max_element, max_element_shift;

/////////////////////////////////////////////////////////////////////////////////////
//...ищем конец элементов заданного *Part (и определяем максимальный номер элемента);
			for ( j = 0; j < nd->geom[0]; j++) 
			if (nd->geom[nd->geom_ptr[j]+2]+1 == 0) {
				if (N_part--) {
					numeration_shift = j;				
					max_element = 1;
				}
				else break;
			}
			else if (max_element < -nd->geom[nd->geom_ptr[j]+2]) max_element = -nd->geom[nd->geom_ptr[j]+2];

			max_element_shift = nd->geom[0]-max_element+numeration_shift; //...нумерация дополнительных элементов - после максимального!!!
			N_part  = nd->geom[0]-(N_part_end = j)+numeration_shift; //...смещение на фактически последний элемент (без учета пропусков)!!!
			NN_geom = nd->geom[0];

//////////////////////////////////////////////////////////////////
//...создаем табличку перенумерации элементов для заданного *Part;
			int *	geom_pernum = (int *)new_struct(max_element*sizeof(int));
			for ( j = numeration_shift; j < N_part_end; j++)
				geom_pernum[-nd->geom[nd->geom_ptr[j]+2]-1] = j;

//////////////////////////////////
//...добавление элеметов в модель;
			if (nd->cond && nd->cond_ptr && ID_elements) 
			for (set = 1; set <= ID_elements[0]; set++) {

//////////////////////////////////////////
//...ищем добавляемое множество элементов;
				for (m = k = 0; ! m && k < nd->cond[0]; k++)
				if ((nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS || nd->cond[l] == (int)GL_ELEMENTS_GENERATE) 
				  && nd->cond[l+2] == ID_elements[set]) m = 1;

//////////////////////////////////////////////////////////////////////////
//...добавляем элеметы из указанного множества в заданный *Part элементов;
				if (nd->geom && nd->geom_ptr && m) {
					N_geom = nd->geom[0];

//////////////////////////////////////
//...просчитываем геометрию элементов;
					if (nd->cond[l] == (int)GL_ELEMENTS) { //...добавляемое множество;
						for ( j = 2; j <= nd->cond[l+1]; j++) if (nd->cond[l+1+j] <= max_element)
							nd->geom_ptr_add(nd->geom[nd->geom_ptr[geom_pernum[nd->cond[l+1+j]-1]]+1], N_geom);
					}
					else
					if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE) { //...добавляемое множество;
						for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) if (j <= max_element)
							nd->geom_ptr_add(nd->geom[nd->geom_ptr[geom_pernum[j-1]]+1], N_geom);
					}

//////////////////////////////////////////
//...перераспределяем геометрию элементов;
					int * geom_new = NULL;
					if ( (geom_new = (int *)new_struct((nd->geom_ptr[N_geom]+1)*sizeof(int))) != NULL) {
						memcpy(geom_new, nd->geom, (nd->geom_ptr[nd->geom[0]]+1)*sizeof(int)); delete_struct(nd->geom);
						nd->geom = geom_new;

//////////////////////
//...переносим данные;
						N_geom = nd->geom[0];
						if (nd->cond[l] == (int)GL_ELEMENTS) {
							for ( j = 2; j <= nd->cond[l+1]; j++) if (nd->cond[l+1+j] <= max_element) { //...;добавляемое множество 
								elem = geom_pernum[nd->cond[l+1+j]-1];
								l0 = nd->geom[nd->geom_ptr[elem]+1]+2;
								memcpy(nd->geom+(i = nd->geom_ptr[nd->geom[0]]), nd->geom+nd->geom_ptr[elem], l0*sizeof(int));
								nd->geom[i+2] = -(++nd->geom[0])+max_element_shift;
							}
						}
						else
						if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE) {
							for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) if (j <= max_element) { 
								elem = geom_pernum[j-1];
								l0 = nd->geom[nd->geom_ptr[elem]+1]+2;
								memcpy(nd->geom+(i = nd->geom_ptr[nd->geom[0]]), nd->geom+nd->geom_ptr[elem], l0*sizeof(int));
								nd->geom[i+2] = -(++nd->geom[0])+max_element_shift;
							}
						}
					}
				}

/////////////////////////////////////////////////////////////////
//...просчитываем и добавляем дополнительное множество элементов;
				if (m && N_geom < nd->geom[0]) {
					for (l = j = k = 0; k < nd->cond[0]; k++) {
						if (nd->cond[(++l)++] == GL_ELEMENTS || nd->cond[l-1] == GL_ELEMENTS_GENERATE) --j;
						l += nd->cond[l];
					}
					nd->cond_ptr_add(4, nd->cond[0]);

//////////////////////////////////////////////////////////////////////////
//...перераспределяем условия и заносим данные о дополнительном множестве;
					int * cond_new = NULL;
					if ( (cond_new = (int *)new_struct((nd->cond_ptr[nd->cond[0]]+1)*sizeof(int))) != NULL) {
						memcpy(cond_new, nd->cond, (nd->cond_ptr[nd->cond[0]-1]+1)*sizeof(int)); delete_struct(nd->cond);
						nd->cond = cond_new;

//////////////////////
//...переносим данные;
						nd->cond[i = nd->cond_ptr[nd->cond[0]-1]] = GL_ELEMENTS_GENERATE;
						nd->cond[i+1] = 4;
						nd->cond[i+2] = --j;
						nd->cond[i+3] = N_geom+1-max_element_shift;
						nd->cond[i+4] = nd->geom[0]-max_element_shift;
						nd->cond[i+5] = 1;
						N_set++;
					}
				}
			}

////////////////////////////////////////////////////////////////////////////
//...запись скорректированных данных в заданный *Part входного файла модели;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = begin_set = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
				PPOS_CUR(id_NODES, pchar, begin_set, upper, "*Elset, elset=");
			}
			for ( j = 1; j < ID_element_part; j++) 
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = count_set = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
				PPOS_CUR(id_NODES, pchar, count_set, upper, "*Elset, elset=");

/////////////////////////////////////////////////////////////////////////////////
//...ищем вхождение элементов в разделе *Part и печатаем все до начала элементов;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type="); swap(pchar[0], temp);
				fprintf(INP, "%s", id_NODES); swap(pchar[0], temp);

////////////////////////////////////////////////////////////////////
//...ищем все вхождения элементов и печатаем их с внешними номерами;
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_element = (upper_element -= (upper_element < upper ? 1 : 0));
					
					nd->element_type(elem, id_NODES+ppos_cur);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					
					fprintf(INP, "*Element, type=%s", one_line);

///////////////////////////////////////////////////////////////
//...печатаем все элементы данной топологии из заданного *Part;
					for (l = nd->geom_ptr[k = numeration_shift]-1; k < N_part_end; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+1]);					
							for (j = 2; j < nd->geom[l]; j++) {
								if (j != j_max) fprintf(INP, ",%7i", nd->geom[l+1+j]);
								else			 fprintf(INP, ",\n%15i", nd->geom[l+1+j]);					
							}
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					for (l = nd->geom_ptr[k = NN_geom]-1; k < nd->geom[0]; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+1]);					
							for (j = 2; j < nd->geom[l]; j++) {
								if (j != j_max) fprintf(INP, ",%7i", nd->geom[l+1+j]);					
								else			 fprintf(INP, ",\n%15i", nd->geom[l+1+j]);					
							}
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				}

/////////////////////////////////////////////////////////////////////////////////////////////
//...просматриваем все подмножества элементов и ищем место, куда дописать новые подмножества;
				ppos_cur = count_set;
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_set = (upper_element -= (upper_element < upper ? 1 : 0));
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
				}
				swap(id_NODES[end_set], temp); fprintf(INP, "%s", id_NODES+end_element);
				swap(id_NODES[end_set], temp);

///////////////////////////////////////////////////////////////////////
//...пишем добавляемые множества элементов с модифицированными именами;
				if (nd->cond && nd->cond_ptr && ID_elements) 
				for (set = 1; set <= ID_elements[0]; set++) {

//////////////////////////////////////////
//...ищем добавляемое множество элементов;
				for (m = j = k = 0; ! m && k < nd->cond[0]-N_set; k++)
				if ((nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS || nd->cond[l] == (int)GL_ELEMENTS_GENERATE) && --j &&
					  nd->cond[l+2] == ID_elements[set]) m = 1;

///////////////////////////////////////////////
//...ищем имя добавляемого множества элементов и пишем множества элементов;
					if (m && N_geom < nd->geom[0]) {
						ppos_cur = begin_set;
						while ((upper_element = ppos_cur) < upper && ++j) {
							PPOS_CUR(id_NODES, pchar, upper_element, upper, "*");
							PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
						}
						if (! j) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
							if (pchar = strstr(one_line, ",")) {
								swap(pchar[0], temp); fprintf(INP, "*Elset, elset=%s_add, generate", one_line);
								swap(pchar[0], temp); fprintf(INP, "\xA");
							}
							else
							if (pchar = strstr(one_line, "\xA")) {
								swap(pchar[0], temp); fprintf(INP, "*Elset, elset=%s_add, generate", one_line);
								swap(pchar[0], temp); fprintf(INP, "%s", pchar);
							}
							fprintf(INP, "", pchar);  i = nd->cond_ptr[nd->cond[0]-N_set];
							fprintf(INP, "%7i,%7i,%7i\n", nd->cond[i+3], nd->cond[i+4], nd->cond[i+5]);					
							N_set--;
						}
					}
				}

////////////////////////////////////////////////
//...печатаем оставшийся хвостик входного файла;
				fprintf(INP, "%s", id_NODES+end_set);
			}
			delete_struct(geom_pernum);
		}
		fclose(INP);
	}
	delete_struct(id_NODES);
}

/////////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR CONVERTING INP-FORMAT          //
/////////////////////////////////////////////////////////////////////
//...подготовка списка зачитанных условий (множества узлов);
unsigned long * Condit_list_nodes(char * id_NODES, unsigned long count, unsigned long upper_limit, int id_status = OK_STATE)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, pposs, upper, upper_element, N_buf = 50, buf_incr = 20,
						* nodes_list = (unsigned long *)new_struct((N_buf*2+1)*sizeof(unsigned long));
		const int STR_SIZE = 250;
		int   elem, j = -1, j_elem = -1;
		char one_line[STR_SIZE+1], * pchar;

//////////////////////////////////////////////////////////////
//...зачитываем специальные подмножеcтва из раздела *Assembly;
		if (id_status) {
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Assembly, name="); upper = ppos_cur;
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Assembly"); upper -= sizeof("*End Assembly");
			
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Instance, name="); upper_element = ppos_cur;
			while (upper_element < upper) {
				PPOS_CUR(id_NODES, pchar, ppos_cur,	upper, "*End Instance");
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*Instance, name=");
			}

//////////////////////////////////
//...заполняем подмножества узлов;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			while ((upper_element = pposs = ppos_cur) < upper) {
				if (nodes_list[0] == N_buf) {
					unsigned long * new_nodes_list = nodes_list; nodes_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(nodes_list, new_nodes_list, (new_nodes_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_nodes_list);
				}
				nodes_list[nodes_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

////////////////////////////////////////////////////////
//...,определяем тип узла (можно расширить этот список);
				if (! ::strncmp(id_NODES+pposs, "CONTACT",  7) || 
						! ::strncmp(id_NODES+pposs, "LOADING",  7) || 
						! ::strncmp(id_NODES+pposs, "LATERAL_BOUND", 13) ||
						! ::strncmp(id_NODES+pposs, "BOTTOM_BOUND",  12)) elem = GL_NODES; else elem = ERR_STATE;
				if (elem == GL_NODES && strstr(one_line, "internal")) elem = ERR_STATE; 
				if (elem == GL_NODES && strstr(one_line, "generate")) elem = GL_NODES_GENERATE; 
	
				if (elem == GL_NODES || elem == GL_NODES_GENERATE)
				if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					nodes_list[nodes_list[0]*2+2] = strlen(one_line);
					nodes_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			}
		}

/////////////////////////////////////
//...просматриваем все разделы *Part;
		ppos_cur = count;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");

//////////////////////////////////
//...заполняем подмножества узлов;
			while ((upper_element = ppos_cur) < upper) {
				if (nodes_list[0] == N_buf) {
					unsigned long * new_nodes_list = nodes_list; nodes_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(nodes_list, new_nodes_list, (new_nodes_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_nodes_list);
				}
				nodes_list[nodes_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////
//...,определяем тип узла;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_NODES_GENERATE; else elem = GL_NODES;

				if (elem == GL_NODES || elem == GL_NODES_GENERATE) 
				if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					nodes_list[nodes_list[0]*2+2] = strlen(one_line);
					nodes_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
		return(nodes_list);
	}
	return(NULL); 
}

////////////////////////////////////////////////////////////////
//...подготовка списка зачитанных условий (множества елементов);
unsigned long * Condit_list_elements(char * id_NODES, unsigned long count, unsigned long upper_limit)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, upper, upper_element, N_buf = 50, buf_incr = 20,
						* elements_list = (unsigned long *)new_struct((N_buf*2+1)*sizeof(unsigned long));
		const int STR_SIZE = 250;
		int  elem, j = -1, j_elem = -1;
		char one_line[STR_SIZE+1], * pchar;

/////////////////////////////////////
//...просматриваем все разделы *Part;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");

//////////////////////////////////////
//...заполняем подмножества элементов;
			while ((upper_element = ppos_cur) < upper) {
				if (elements_list[0] == N_buf) {
					unsigned long * new_elements_list = elements_list; elements_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(elements_list, new_elements_list, (new_elements_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_elements_list);
				}
				elements_list[elements_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////////
//...,определяем тип элемента;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_ELEMENTS_GENERATE; else elem = GL_ELEMENTS; 

				if (elem == GL_ELEMENTS || elem == GL_ELEMENTS_GENERATE) 
				if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					elements_list[elements_list[0]*2+2] = strlen(one_line);
					elements_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
		return(elements_list);
	}
	return(NULL); 
}

////////////////////////////////////////////////////////////
//...конвертация геометрии из 3D модели Abaqus (файл *.gmt);
void Convert3D_gmt(char * ch_NODES, CGrid * nd, int ID_part)
{
	char * id_NODES = read_struct_ascii(ch_NODES), buff[2000];
	if (nd && id_NODES) {
		const int STR_SIZE = 250, ID_hexa_elem = 300;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, count, count_part, count_length, * elements_list = NULL;
		char one_line[STR_SIZE+1], * pchar, temp = '\0';
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, ".gmt");
		FILE * GMT = fopen(buff, "w");
		if (GMT) {
			int j, k, l, m, N_part = ID_part, N_part_nodes_beg, N_part_nodes_end, N_part_element_beg, N_part_element_end, dim = 3;

//////////////////////////
//...запись имени проекта;
			char job_name[STR_SIZE+1], model_name[STR_SIZE+1];
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "** Job name: ");
			ONE_LINE(id_NODES, pchar, ppos_cur, upper_limit, job_name, STR_SIZE); count = 0; upper = (unsigned long)strlen(job_name);
			PPOS_CUR(job_name, pchar, count, upper, " Model name: "); job_name[upper-2] = 0;

			if (pchar) pchar[0] = 0; strcpy(model_name, job_name+count); 	
			fprintf(GMT, "[MESH NAME]\n%s\n%s\n<END>\n%i   // размерность пространства\n", job_name, model_name, dim);
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

//////////////////////////////////////
//...ищем конец узлов заданного *Part;
			for ( j = 0; j < nd->N; j++) 
			if (nd->hit[j] == 1) {//...первый узел part'а (с номером 1);
				if (N_part--) N_part_nodes_beg = j;				
				else break;
			}
			N_part_nodes_end = j;

//////////////////////////////////////////
//...ищем конец элементов заданного *Part;
			for (N_part = ID_part, j = 0; j < nd->geom[0]; j++) 
			if (nd->geom[nd->geom_ptr[j]+2]+1 == 0) {//...первый элемент part'а (с номером 1);
				if (N_part--) N_part_element_beg = j;				
				else break;
			}
			N_part_element_end = j;

/////////////////////////////////
//...поиск начала нужного part'а;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_part = ppos_cur;
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); 
			}
			for ( j = 1; j < ID_part; j++) {
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
				if ((upper = ppos_cur) < upper_limit) {
					PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); count_part = ppos_cur; 
				}
			}

////////////////////////////////////////////////
//...запись аттрибутов материала (визуализация);
			unsigned long Color[4] = {0x32C8C8C8, 0x00808080, 0x0080FFFF, 0xC80000FF}, index = 0, index_max = 4, count_light = 0, count_max = 7;
			elements_list = Condit_list_elements(id_NODES, 0, upper_limit); ppos_cur = count_part;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");

			fprintf(GMT, "[MATERIAL COLORS]\n");
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
			
				if ( (pchar = strstr(one_line, ", material")) != NULL) {
					pchar[0] = temp; count_length = strlen(one_line);

					for (count = 0; count < elements_list[0]; count++)  
						if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;
					
					if (count < elements_list[0]) {
						count_light++;
						float blue  = (float)(((count_light/1)%2) ? ((Color[index] & 0x000000ff) >>  0) : (Color[index] & 0xff000000) >> 24)/255.f,
								green = (float)(((count_light/2)%2) ? ((Color[index] & 0x0000ff00) >>  8) : (Color[index] & 0xff000000) >> 24)/255.f,
								red   = (float)(((count_light/4)%2) ? ((Color[index] & 0x00ff0000) >> 16) : (Color[index] & 0xff000000) >> 24)/255.f;
						fprintf(GMT, "%7i   %3.2g   %3.2g   %3.2g\n", (int)count+1, red, green, blue);
						fprintf(GMT, "Материал ID=%i   // ассоциируется с областью: %s\n<END>\n", (int)count+1, one_line);
						if (count_light == count_max) {
							 count_light = 0; index++;
							 if (index == index_max-1) {
								  count_max = 4; count_light = 2;
							 }
							 if (index == index_max) {
								  index = 0; count_max = 7;
								  Color[0] -= 0x10202020;
								  Color[1] -= 0x00101010;
								  Color[2] -= 0x00102020;
								  Color[3] -= 0x20000020;
							 }
						}
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			}
			fprintf(GMT, "0   // признак конца набора\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

//////////////////////////////////
//...запись фигур (пустой раздел);
			fprintf(GMT, "[FIGURES]\n0   // признак конца набора\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

/////////////////////////
//...запись узлов модели;
			fprintf(GMT, "[NODES]\n// масштабирующие множители:\n//__________X_______Y_______Z:\n           1.0     1.0     1.0\n");
			fprintf(GMT, "// информация об узлах:\n//_ID_______X_______Y_______Z:\n");
			for (k = N_part_nodes_beg; k < N_part_nodes_end; k++) //...базовые узлы;
			fprintf(GMT, "%7i   %12g    %12g    %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);
			fprintf(GMT, "0   // признак конца набора\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

////////////////////////////
//...запись элеметов модели;
			fprintf(GMT, "[ELEMENTS]\n// информация об элементах:\n");
			fprintf(GMT, "//_ID__mat__type_____1_____2_____3_____4_____5_____6_____7_____8    базовые узлы\n");
			fprintf(GMT, "// дополнительные узлы (по 10)\n//              _____9____10____11____12____13____14____15____16____17____18____19\n//              ____20____21 ... номера узлов\n");
			fprintf(GMT, "//кол-во добавочных узлов на рёбрах:   __1__2__3__4__5:\n");
			if (nd->geom && nd->geom_ptr) {
				for (l = nd->geom_ptr[k = N_part_element_beg];	k < N_part_element_end; l = nd->geom_ptr[++k]) { //...базовые элементы;
					if (nd->geom[l] == (int)GL_BOXS)	 m = ID_hexa_elem; else m = ERR_STATE;
					if (m != ERR_STATE) {
						fprintf(GMT, "%7i   %7i   %7i", -nd->geom[l+2], -nd->geom[l+3], m);
						for (j = 2; j < nd->geom[l+1]; j++)
						fprintf(GMT, "   %7i", nd->geom[l+j+2]);
						fprintf(GMT, "\n");
					}
				}
			}
			fprintf(GMT, "0   // признак конца набора\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");
			fclose (GMT);
		}
	}
	delete_struct(id_NODES);
}

/////////////////////////////////////////////////////////////////////////
//...конвертация нагрузок и закреплений из 3D модели Abaqus (файл *.prb);
void Convert3D_prb(char * ch_NODES, CGrid * nd, int ID_part)
{
	char * id_NODES = read_struct_ascii(ch_NODES), buff[2000];
	if (nd && id_NODES) {
		const int STR_SIZE = 250; char one_line[STR_SIZE+1], material_line[STR_SIZE+1], * pchar, * ppos, * material_name, temp = '\0';
		unsigned long ppos_cur = 0, ppos_material, ppos_save, upper_limit, upper, upper_element, upper_material, 
						  count, count_part, count_beg, count_length, * nodes_list = NULL, * elements_list = NULL;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, ".prb");
		FILE * PRB = fopen(buff, "w");
		if (PRB) {
			int k, l, j, j_max = 10, n_steps = 0, buf_steps = 20, buf_incr = 10, m[3];
			double * steps = (double *)new_struct(buf_steps*sizeof(double)), sum_steps = 0.;

//////////////////////////
//...запись имени проекта;
			char job_name[STR_SIZE+1], model_name[STR_SIZE+1];
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "** Job name: "); count_beg = ppos_cur;
			ONE_LINE(id_NODES, pchar, ppos_cur, upper_limit, job_name, STR_SIZE); count = 0; upper = (unsigned long)strlen(job_name);
			PPOS_CUR(job_name, pchar, count, upper, " Model name: "); job_name[upper-2] = 0;

			if (pchar) pchar[0] = 0; strcpy(model_name, job_name+count); 	
			fprintf(PRB, "[PROBLEM NAME]\n%s\n%s\n<END>\n", job_name, model_name);
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

//////////////////////////////////////////////////
//...запись аттрибутов задачи (старт, тип задачи);
			fprintf(PRB, "[START]\n0   // 0 - счёт с \"0\", иначе - описание ниже\n// файл с информацией для рестарта:\nD:\\Progs\\UWAY\\DATA\\Mike_Test\\Name_RESTART.dat\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

			fprintf(PRB, "[PROBLEM]\n1 3   // тип задачи и вид напряжённого состояния\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

////////////////////////////////////////////////
//...подсчет и запись аттрибутов времени (шаги);
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			while ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Step"); upper -= sizeof("*End Step");
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				if (n_steps == buf_steps) {
					double * prev_steps = steps; steps = (double *)new_struct((buf_steps += buf_incr)*sizeof(double));
					memcpy(steps, prev_steps, n_steps*sizeof(double));
					delete_struct(prev_steps);
				}
				if ((ppos = strstr(one_line, ",")) != NULL) steps[n_steps] = user_strtod(ppos+1);
				else													  steps[n_steps] = 1.;
				sum_steps += steps[n_steps++];
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			}
			fprintf(PRB, "[TIMES]\n0.0	// время начала    (T0)\n");
			fprintf(PRB, "%2.1f	// время окончания (Tn)\n",   sum_steps);
			fprintf(PRB, "%2.1f	// шаг по времени  (Step)\n", sum_steps/n_steps);
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

///////////////////////
//...вывод результатов;
			fprintf(PRB, "[RESULT OUTS]\n0		// \"0\" - вывод результатов в текущую директорию\n");
			fprintf(PRB, "2		// тип вывода\n// 0 - узлы + точки интегрирования (элементы)\n// 1 - узлы + узлы (элементы)\n// 2 - узлы + центры тяжести (элементы)\n// 3 - все типы\n");
			fprintf(PRB, "0		// тип формата вывода результатов\n// 0 - формат UWay  (параллельный)\n// 1 - формат FEMAP (последовательный)\n");
			fprintf(PRB, "// Моменты времени, в которые будет производиться вывод результатов расчётов:\n");
			fprintf(PRB, "//____1_______2_______3_______4_______5_______6_______7_______8_______9______10:");
			for (j = 0; j < n_steps; j++) {
				if (j%j_max) fprintf(PRB, "%8.1f", steps[j]);
				else			 fprintf(PRB, "\n%8.1f", steps[j]);
			}
			if (j%j_max) fprintf(PRB, "      0      //<<< \"0\" - the end of data\n");
			else			 fprintf(PRB, "\n      0      //<<< \"0\" - the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			
///////////////////////
//...начальные условия;
			fprintf(PRB, "[STATE 0]\n0   // 0 - все начальные условия = 0.0;\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

/////////////////////////////////////////////
//...запись условий закрепления узлов модели;
			fprintf(PRB, "[FIXED NODES]\n// момент изменения закреплений узлов и комментарий:\n0.0   Комментарий к t=0.0\n");
			fprintf(PRB, "//_N/F__ID__X__Y__Z:   закрепления (N/F: \"0\"-узлы; \"1\"-фигуры; ID-узел/фигура)\n");

//////////////////////////////////////////////////////////////////////////////
//...поиск и запись условий закрепления до первого step'а (начальные условия);
			nodes_list = Condit_list_nodes(id_NODES, 0, upper_limit); ppos_cur = upper = count_beg;
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*Step, name=");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Boundary"); m[0] = m[1] = m[2] = 0;

			while (ppos_cur < upper) {
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);

				if ((ppos = strstr(one_line, ".")) != NULL)
				if ((pchar = strstr(ppos, ", ")) != NULL) {
					swap(temp, pchar[0]); count_length = strlen(ppos+1);
					for (count = 0; count < nodes_list[0]; count++)  
						if (nodes_list[count*2+2] == count_length && ! strncmp(ppos+1, id_NODES+nodes_list[count*2+1], count_length)) break;
					swap(temp, pchar[0]); k = (int)count; 
					if ((ppos = strstr(pchar, ", ")) != NULL) {
						ppos[0] = temp; j = atoi(pchar+1);
						if (0 < j && j < 4) m[j-1] = 1;
					}
					ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
					if (one_line[0] != '*' && (pchar = strstr(one_line, ",")) != NULL) {
						if ((ppos = strstr(pchar, ", ")) != NULL) {
							ppos[0] = temp; j = atoi(pchar+1);
							if (0 < j && j < 4) m[j-1] = 1;
						}
						ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
						if (one_line[0] != '*' && (pchar = strstr(one_line, ",")) != NULL) {
							if ((ppos = strstr(pchar, ", ")) != NULL) {
								ppos[0] = temp; j = atoi(pchar+1);
								if (0 < j && j < 4) m[j-1] = 1;
							}
						}
					}
				}
///////////////////////////////////
//...запись конкретных закреплений;
				if (nd->cond_ptr && nd->cond) {
					l = nd->cond_ptr[k];
					if (nd->cond[l] == (int)GL_NODES) {
						for (j = 2; j <= nd->cond[l+1]; j++) 
							fprintf(PRB, "    0   %7i   %7i   %7i   %7i\n", nd->cond[l+1+j], m[0], m[1], m[2]);
					}
					else
					if (nd->cond[l] == (int)GL_NODES_GENERATE) {
						for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) 
							fprintf(PRB, "    0   %7i   %7i   %7i   %7i\n", j, m[0], m[1], m[2]);
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Boundary"); m[0] = m[1] = m[2] = 0;
			}
			fprintf(PRB, "<END>\n");
			fprintf(PRB, "<END>   //<<< the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			delete_struct(steps);

/////////////////////////////////
//...поиск начала нужного part'а;
			ppos_cur = 0;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_part = ppos_cur;
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); 
			}
			for ( j = 1; j < ID_part; j++) {
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
				if ((upper = ppos_cur) < upper_limit) {
					PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); count_part = ppos_cur; 
				}
			}

/////////////////////////////////////////////////////////////////////////////
//...запись аттрибутов материала (физические параметры материалов и грунтов);
			int material_first = 1;
			elements_list = Condit_list_elements(id_NODES, 0, upper_limit); ppos_cur = count_part;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			
			fprintf(PRB, "[MATERIAL DATA]\n");
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
			
				if ( (pchar = strstr(one_line, ", material=")) != NULL) {
					swap(pchar[0], temp); count_length = strlen(one_line);

					for (count = 0; count < elements_list[0]; count++)  
						if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;
					
					swap(pchar[0], temp); material_name = pchar+strlen(", material=");
					if ( (ppos = strstr(material_name, "\xD")) != NULL || (ppos = strstr(material_name, "\xA")) != NULL) ppos[0] = temp;
					count_length = strlen(material_name);

					if (count < elements_list[0]) {
						if (! material_first)
						fprintf(PRB, "//................................................\n");

//////////////////////////////////////////////
//...зачитываем материал и определяем его тип;
						ppos_material = ppos_cur;
						PPOS_CUR(id_NODES, pchar, ppos_material, upper_limit, "*Material, name=");

						while ((upper_material = ppos_material) < upper_limit) {
							PPOS_CUR(id_NODES, pchar, upper_material, upper_limit, "*Material, name=");  upper_material -= strlen("*Material, name=");
							ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE); ppos_save = ppos_material;

							if ((ppos = strstr(material_line, "\xD")) != NULL || (ppos = strstr(material_line, "\xA")) != NULL) ppos[0] = temp;
							if (strlen(material_line) == count_length && ! strncmp(material_name, material_line, count_length)) {
								double Gamma = 0.;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Density");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										Gamma = user_strtod(material_line); 
									}
								}
								double E = 0., nju = 0.; int  model = 0, aniso = 0;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Elastic");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										E = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju = user_strtod(pchar+1); 
									}
									model = 1; 
								}
								double E1 = 0., E2 = 0., E3 = 0., nju1 = 0., nju2 = 0., nju3 = 0., G1 = 0., G2 = 0., G3 = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Elastic, type=ENGINEERING CONSTANTS");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										E1 = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										E2 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										E3 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju1 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										nju2 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju3 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										G1 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										G2 = user_strtod(pchar+1); 
									}
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										G3 = user_strtod(material_line); 
									}
									aniso = 2; 
								}
								double angle_frict = 0., flow_ratio = 0., angle_dilat = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Drucker Prager");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										angle_frict = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										flow_ratio = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										angle_dilat = user_strtod(ppos+1); 
									}
								   model = 9;
								}
								double yield_stress = 0., plastic_strain = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Drucker Prager Hardening");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										yield_stress = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										plastic_strain = user_strtod(pchar+1); 
									}
								}

//////////////////////////////
//...раздел печати материалов;
								fprintf(PRB, "%i	%i  %i     // номера: материала, модели, анизотропии\n", (int)count+1, model, aniso);
								switch (model) {
								case 1: switch (aniso) {
									case 0: fprintf(PRB, "// E      v       Gamma:\n");
									fprintf(PRB, "%g     %g     %g\n", E, nju, Gamma);
									break;
									case 2: fprintf(PRB, "// E1      E2      E3      v1      v2      v3      G1      G2      G3       Gamma:\n");
									fprintf(PRB, "%g     %g     %g     %g     %g     %g     %g     %g     %g     %g\n", E1, E2, E3, nju1, nju2, nju3, 
										G1, G2, G3, Gamma);
									}
								break;
								case 9:  switch (aniso) {
									case 0: fprintf(PRB, "// E      v      C     fi    Gamma\n");
									fprintf(PRB, "%g     %g     %g     %g     %g\n", E, nju, yield_stress, angle_frict, Gamma);
									break;
									case 2: fprintf(PRB, "// E1      E2      E3      v1      v2      v3      G1      G2      G3      C     fi       Gamma:\n");
									fprintf(PRB, "%g     %g     %g     %g     %g     %g     %g     %g     %g     %g     %g     %g\n", E1, E2, E3, nju1, nju2, nju3, 
										G1, G2, G3, yield_stress, angle_frict, Gamma);
									}
								}
							}
							ppos_material = upper_material+strlen("*Material, name=");
						}
						fprintf(PRB, "%s		// название материала\n<END>\n", material_name);
						material_first = 0;
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			}
			fprintf(PRB, "0   //<<< признак конца набора\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

//////////////////////////////////////////
//...смена типа материала (пустой раздел);
			fprintf(PRB, "[VARIABLE MATERIALS]\n");
			//fprintf(PRB, "//0   // есть/нет изменяемые материалы (\"0\" - нет)\n");
			fprintf(PRB, "0   //<<< признак конца набора\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

////////////////////////////////////////////////////////////
////...запись условий нагружения элементов или узлов модели;
			fprintf(PRB, "[LOADS]\n");
//			fprintf(PRB, "// \"нагрузка\" № 1    // f(x,t)=f(t)*f(x)\n103	               // тип нагрузки - объёмная\n");
//			fprintf(PRB, "<Собственный вес грунта>    // название нагрузки\n0	// признак подкачки с диска: 0-RAM; 1-ASCII; 2-Binary\n");
//			fprintf(PRB, "0							// ф-ция нагружения от пространственных координат - f(x)\n3							// кол-во точек по времени\n");
//			fprintf(PRB, "// t(0.0 -:- +1.0)  p(-1.0 -:- +1.0)\n			0.0					0.0\n			1.0					1.0\n			3.0 					1.0\n");
//			fprintf(PRB, "0.0 -1.0					// направления\n9.81						// множитель\n1.0						// кол-во повторов за расчётный период\n");
//			fprintf(PRB, "0		// кол-во загруженных объектов (0 или >= объектов - для всех объектов)\n");
//			fprintf(PRB, "//..............................................................................\n");
//			fprintf(PRB, "// \"нагрузка\" № 2:	// f(x,t)=f(t)*f(x)\n106	               // тип нагрузки - равномерно распределённая\n");
//			fprintf(PRB, "<нагрузка от К6/ж>  // название нагрузки\n0	// признак подкачки с диска: 0-RAM; 1-ASCII; 2-Binary\n");
//			fprintf(PRB, "0							// ф-ция нагружения от пространственных координат - f(x)\n2							// кол-во точек по времени\n");
//			fprintf(PRB, "// t(0.0 -:- +1.0)  p(-1.0 -:- +1.0)\n			2.0   	  			0.0\n			3.0 					1.0\n");
//			fprintf(PRB, "0.0 -1.0					// направления\n0.24						// множитель\n1.0						// кол-во повторов за расчётный период\n");
//
//////////////////////////////////////////////////////////////
////...вычисление количества и идетификация нагруженных узлов;
//			for (j_node = 0, k = 0; k < nd->N; k++)
//			if (nd->hit[k] == ID_loading_node) j_node++;
//			fprintf(PRB, "%i		// кол-во загруженных объектов (0 или >= объектов - для всех объектов)\n", j_node);
//			fprintf(PRB, "1		// \"+1\" - для объектов, \"-1\" - исключая объекты\n");
//			fprintf(PRB, "// №№ загруженных объектов:\n");
//			fprintf(PRB, "//_____1_____2_____3_____4_____5_____6_____7_____8_____9____10:");
//
//			for (j_node = 0, k = 0; k < nd->N; k++) //...индексы нагуженных узлов;
//			if (nd->hit[k] == ID_loading_node) {
//				if (! (j_node % 10))	fprintf(PRB, "\n ");
//				fprintf(PRB, " %7i", k+1);	j_node++;
//			}
//			fprintf(PRB, "\n");
			fprintf(PRB, "0		//<<<<<<<<<<<<<<<<<<<<<<< the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			fclose (PRB);
		}
		delete_struct(nodes_list);
		delete_struct(elements_list);
	}
	delete_struct(id_NODES);
}

////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR ESHELBY PROBLEMS          //
////////////////////////////////////////////////////////////////
double TakeLayer_EH(int N, double * ff, double * kk, double * ll)
{
	CCohes3D draft;
	CLame3D ldraft; double EH;
	if (ll == NULL) EH = ldraft.TakeLayer_kk(N, ff, kk);
	else				 EH =  draft.TakeLayer_kk(N, ff, kk, ll);
	return(EH);
}

double TakeCylinder_EH(int N, double * ff, double * kk)
{
	CLame2D ldraft; double EH;
	EH = ldraft.TakeLayer_kk(N, ff, kk);
	return(EH);
}

double TakeCylinder_KH(int N, double * ff, double * kp, double * mu)
{
	CLame2D ldraft; double KH;
	KH = ldraft.TakeLayer_kk(N, ff, kp, mu);
	return(KH);
}

double TakeCylinder_GH(int N, double * ff, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, mu, nj);
	return(GH);
}

double TakeCylinder_GH(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, kp, mu, nj);
	return(GH);
}

double TakeCylinder_SH(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_sh(N, ff, kp, mu, nj);
	return(GH);
}

double TakeCylinder_LH(int N, double * ff, double * lm, double * mu, double * nj)
{
	CLame2D ldraft; double LH;
	LH = ldraft.TakeLayer_kk(N, ff, lm, mu, nj);
	return(LH);
}

double TakeSphere_KH(int N, double * ff, double * kv, double * mu)
{
	CLame3D ldraft; double KH;
	KH = ldraft.TakeLayer_kk(N, ff, kv, mu);
	return(KH);
}

double TakeSphere_GH(int N, double * ff, double * kv, double * mu, double * nj)
{
	CLame3D ldraft; double GH;
	GH = ldraft.TakeLayer_GG(N, ff, kv, mu, nj);
	return(GH);
}

double TakeSphere_GH_det(double c0, double nju1, double nju2, double E1, double E2)
{
	CLame3D ldraft; double GH;
	ldraft.set_fasa_hmg(nju2, nju1, E2/(1.+nju2)*.5, E1/(1.+nju1)*.5);
	GH = ldraft.TakeEshelby_shear_det(c0);
	return(GH);
}

double TakeSphere_volm_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		KH = draft.TakeEshelby_volm_two(c0);
	}
	return(KH);
}

double TakeSphere_shear_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		GH = ldraft.TakeEshelby_shear_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		GH = draft.TakeEshelby_shear_two(c0);
	}
	return(GH);
}

double TakeSphere_volm_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		KH = draft.TakeEshelby_volm_sym(c0);
	}
	return(KH);
}

double TakeSphere_shear_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 1e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) <= 0.01 && fabs(l2) <= 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear_old(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}

double TakeSphere_shear(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 2e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) < 0.01 && fabs(l2) < 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}
