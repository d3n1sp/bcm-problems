#include "stdafx.h"

#include "unit_mes.h"
#include "cprofile.h"
#include "cpro_bars.h"
#include "cpro_wire.h"
#include "cpro_cable.h"
#include "cpro_tower.h"
#include "cpro_fitting.h"

extern const char *     COLONSTR = "";    //...символьная строка, обозначающая начало записи;
extern const char * SEMICOLONSTR = ";";   //...символьная строка, обозначающая конец записи;
extern const char * STRSEPARATOR = "\xA"; //...символьная строка, обозначающая конец строки;
extern const char * PARSEPARATOR = " ";   //...символьная строка, обозначающая разделитель параметров;
extern const char *     TABLESTR = ">>";  //...символьная строка, обозначающая начало таблицы;
extern const char *   VERSIONSTR = "##";  //...символьная строка, обозначающая вхождение кода версии таблицы;

///////////////////////////////////////////////////////////////////////////////////////////////////
//...проверка попадания точки в пространственный размерный обpазец
int Sample3DIn(void * context, int & k, double * P) 
{
//  Context * cont = (Context *)context;
/* if (cont) return sample3D_in(cont->sm, k, P);
  else*/      return (-1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...вспомогательные функции для абсолютной идентификации образца
int is_Lame3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_LM3D_SAMPLES &&
                    cont->N <  SHIFT_LM3D_SAMPLES+NUM_LM3D_SAMPLES);
  else      return (0);
}

int is_Maxw3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_MX3D_SAMPLES &&
                    cont->N <  SHIFT_MX3D_SAMPLES+NUM_MX3D_SAMPLES);
  else      return (0);
}

int is_Acou3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_AU3D_SAMPLES &&
                    cont->N <  SHIFT_AU3D_SAMPLES+NUM_AU3D_SAMPLES);
  else      return (0);
}

int is_Mapi3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_MP3D_SAMPLES &&
                    cont->N <  SHIFT_MP3D_SAMPLES+NUM_MP3D_SAMPLES);
  else      return (0);
}

int is_Acou2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_AU2D_SAMPLES &&
                    cont->N <  SHIFT_AU2D_SAMPLES+NUM_AU2D_SAMPLES);
  else      return (0);
}

int is_Skin2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_SK2D_SAMPLES &&
                    cont->N <  SHIFT_SK2D_SAMPLES+NUM_SK2D_SAMPLES);
  else      return (0);
}

int is_Mapi2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_MP2D_SAMPLES &&
                    cont->N <  SHIFT_MP2D_SAMPLES+NUM_MP2D_SAMPLES);
  else      return (0);
}

int is_BarsGOST2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_BARS_GOST_SAMPLES &&
                    cont->N <  SHIFT_BARS_GOST_SAMPLES+NUM_BARS_GOST_SAMPLES);
  else      return (0);
}

int is_Bars2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_BARS_SAMPLES &&
                    cont->N <  SHIFT_BARS_SAMPLES+NUM_BARS_SAMPLES);
  else      return (0);
}

int is_Wire2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_WIRE_SAMPLES &&
                    cont->N <  SHIFT_WIRE_SAMPLES+NUM_WIRE_SAMPLES);
  else      return (0);
}

int is_Cable2D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_CABLE_SAMPLES &&
                    cont->N <  SHIFT_CABLE_SAMPLES+NUM_CABLE_SAMPLES);
  else      return (0);
}

int is_Tower3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_TOWER_SAMPLES &&
                    cont->N <  SHIFT_TOWER_SAMPLES+NUM_TOWER_SAMPLES);
  else      return (0);
}

int is_Fitting(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N >= SHIFT_FITTING_SAMPLES &&
                    cont->N <  SHIFT_FITTING_SAMPLES+NUM_FITTING_SAMPLES);
  else      return (0);
}

int is_AuGp3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_AUGPF3D);
  else      return (0);
}

int is_AuIb3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_AUGPF3D+1);
  else      return (0);
}

int is_LmGp3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_LMGPF3D);
  else      return (0);
}

int is_LmIb3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_LMGPF3D+1);
  else      return (0);
}

int is_MpGp3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_MPGPF3D);
  else      return (0);
}

int is_MpIb3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_MPGPF3D+1);
  else      return (0);
}

int is_MxGp3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_MXGPF3D);
  else      return (0);
}

int is_MxIb3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_MXGPF3D+1);
  else      return (0);
}

int is_SkGp3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_SKGPF3D);
  else      return (0);
}

int is_SkIb3D(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_SKGPF3D+1);
  else      return (0);
}

int is_Chain(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return (cont->N == NUM_SAMPLE_CHAIN);
  else      return (0);
}

///////////////////////////////////
//...уничтожение контекста образца;
void DeleteContext(void *& context)
{
  Context * cont = (Context *)context;
  if (cont) {
     delete_table (cont->table, cont->static_char);
     delete_struct(cont->param);
     delete_struct(cont->sample_comment);
     delete_struct(cont->sample_description);
	  if (cont->static_char == NULL_STATE) {
		  delete_struct(cont->sample_name);
		  delete_struct(cont->GOST_name);
	  }
     if (cont->sm) delete(cont->sm);
     delete[](cont); context = NULL;
  }
  return;
}

////////////////////////////////////////////////////////////////////
//...вспомогательная функция освобождения памяти от таблицы образца;
void delete_table(Table *& tab, int id_static_char)
{
	if (tab) {
		int i, * m;
		Record * rec;
		delete_struct(tab->table_units);
		for ( i = 0; i < tab->N; i++) {
			rec = &tab->table[i];
			if (ERR_TYPE_RECORD < rec->type && rec->type < INT_TYPE_RECORD) {
				delete_struct(((char **)rec->parm_values)[0]);
				if (id_static_char < ADDITIONAL_STATE)
					for (int j = 1; j <= tab->N_group; j++)
						delete_struct(((char **)rec->parm_values)[j]);
			}
			if (id_static_char == NULL_STATE)
				delete_struct(m = (int *)rec->parm_name);
			if (id_static_char < SPECIAL_STATE)
				delete_struct(m = (int *)rec->parm_values);
		}
		if (tab->table_names && id_static_char == NULL_STATE)
			for (i = 0; i < tab->N_group; i++)
				delete_struct(tab->table_names[i]);
		if (id_static_char < ADDITIONAL_STATE)
		delete_struct(m = (int *)tab->table_names);
		delete_struct(m = (int *)tab); tab = NULL;
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция установки системы единиц во всех сечениях образца;
void set_table_units(Table * tab, int id_units)
{
  for (int k = 0;
       tab && tab->table_units && k < tab->N_group; k++) tab->table_units[k] = id_units;
  return;
}

//////////////////////////////////////////////////////////////////////////
//...вспомогательная функция установки одного параметра в таблице образца;
void set_table_param(Table * tab, int k, Num_type_Record type, char * parm_name, void * parm_values)
{
  if (tab && --k >= 0) {
      tab->table[k].type        = type;
      tab->table[k].parm_name   = parm_name;
      tab->table[k].parm_values = parm_values;
  }
  return;
}

/////////////////////////////////////////////////
//...добавление и удаление параметров из таблицы;
void add_table_param(Table *& tab, int N, Shablon * records, int N_ini)
{
	if (tab && N > 0) {
		if (N_ini  < 0 || N_ini > tab->N) N_ini = tab->N;
		Table * table = (Table *)new_struct(sizeof(Table)+((N += tab->N)-1)*sizeof(Record));

		memcpy(table, tab, sizeof(Table)+(N_ini-1)*sizeof(Record));
		memcpy((char *)table+sizeof(Table)+(N-tab->N+N_ini-1)*sizeof(Record), 
				 (char *)tab+sizeof(Table)+(N_ini-1)*sizeof(Record), (tab->N-N_ini)*sizeof(Record));
		delete [] tab; tab = table;

		for (int i_last = 0, i = tab->N; i < N; i++) {
			if (records[i_last].parm_type < CHAR_TYPE_RECORD); else
			if (records[i_last].parm_type <  INT_TYPE_RECORD) {
					char ** parm_values = (char **)new_struct((table->N_group+1)*sizeof(char *));
					set_table_param(tab, N_ini+i-tab->N+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
			}
			else
			if (records[i_last].parm_type < FLOAT_TYPE_RECORD) {
					int * parm_values = (int *)new_struct((table->N_group+1)*sizeof(int));
					set_table_param(tab, N_ini+i-tab->N+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
			}
			else
			if (records[i_last].parm_type < DOUBLE_TYPE_RECORD) {
					float * parm_values = (float *)new_struct((table->N_group+1)*sizeof(float));
					set_table_param(tab, N_ini+i-tab->N+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
			}
			else
			if (records[i_last].parm_type < NUM_OF_TYPE_RECORD) {
					double * parm_values = (double *)new_struct((table->N_group+1)*sizeof(double));
					set_table_param(tab, N_ini+i-tab->N+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
			}
			if (i+1 < N && records[i_last+1].parm_type != ERR_TYPE_RECORD) i_last = i+1-tab->N;
		}
		tab->N = N;
	}
}

void del_table_param(Table *& tab, int N, int N_fin, int id_static_char)
{
	if (tab && N > 0) { //...удаляем "с хвоста";
		if (N_fin  < 0 || N_fin > tab->N) N_fin = tab->N;
		int m = max(N_fin-N, 0);

		Record * rec;
		for (int i = m; i < N_fin; i++) {
			rec = &tab->table[i];
			if (ERR_TYPE_RECORD < rec->type && rec->type < INT_TYPE_RECORD) {
				delete_struct(((char **)rec->parm_values)[0]);
				if (id_static_char < ADDITIONAL_STATE)
					for (int j = 1; j <= tab->N_group; j++)
						delete_struct(((char **)rec->parm_values)[j]);
			}
			if (id_static_char == NULL_STATE) 
				delete [] rec->parm_name;
			if (id_static_char < SPECIAL_STATE) 
				delete [] rec->parm_values;
		}
		Table * table = (Table *)new_struct(sizeof(Table)+((N = tab->N-N_fin+m)-1)*sizeof(Record));

		memcpy(table, tab, sizeof(Table)+(m-1)*sizeof(Record));
		memcpy((char *)table+sizeof(Table)+(m-1)*sizeof(Record), 
				 (char *)tab+sizeof(Table)+(N_fin-1)*sizeof(Record), (tab->N-N_fin)*sizeof(Record));
		delete [] tab; tab = table;
		tab->N = N;
	}
}

///////////////////////////////////////////////////////
//...добавление и удаление групп параметров из таблицы;
void add_table_group(Table * tab, int N_group, char ** table_names, int N_ini, int id_static_char)
{
	if (tab && N_group > 0 && id_static_char != SPECIAL_STATE) {
		if (N_ini  < 0 || N_ini > tab->N_group) N_ini = tab->N_group;
		N_group += tab->N_group;

		for (int i = 0; i < tab->N; i++) {
			if (tab->table[i].type < CHAR_TYPE_RECORD); else
			if (tab->table[i].type <  INT_TYPE_RECORD) {
				char ** new_parm_values = (char **)new_struct((N_group+1)*sizeof(char *));
				memcpy(new_parm_values, tab->table[i].parm_values, (N_ini+1)*sizeof(char *));
				memcpy(new_parm_values+N_group-tab->N_group+N_ini+1, (char **)tab->table[i].parm_values+N_ini+1, (tab->N_group-N_ini)*sizeof(char *));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < FLOAT_TYPE_RECORD) {
				int * new_parm_values = (int *)new_struct((N_group+1)*sizeof(int));
				memcpy(new_parm_values, tab->table[i].parm_values, (N_ini+1)*sizeof(int));
				memcpy(new_parm_values+N_group-tab->N_group+N_ini+1, (int *)tab->table[i].parm_values+N_ini+1, (tab->N_group-N_ini)*sizeof(int));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < DOUBLE_TYPE_RECORD) {
				float * new_parm_values = (float *)new_struct((N_group+1)*sizeof(float));
				memcpy(new_parm_values, tab->table[i].parm_values, (N_ini+1)*sizeof(float));
				memcpy(new_parm_values+N_group-tab->N_group+N_ini+1, (float *)tab->table[i].parm_values+N_ini+1, (tab->N_group-N_ini)*sizeof(float));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < NUM_OF_TYPE_RECORD) {
				double * new_parm_values = (double *)new_struct((N_group+1)*sizeof(double));
				memcpy(new_parm_values, tab->table[i].parm_values, (N_ini+1)*sizeof(double));
				memcpy(new_parm_values+N_group-tab->N_group+N_ini+1, (double *)tab->table[i].parm_values+N_ini+1, (tab->N_group-N_ini)*sizeof(double));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
		}
		if (tab->table_units) {
			int * new_table_units = (int *)new_struct(N_group*sizeof(int));
			memcpy(new_table_units, tab->table_units, N_ini*sizeof(int));
			memcpy(new_table_units+N_group-tab->N_group+N_ini, tab->table_units+N_ini, (tab->N_group-N_ini)*sizeof(int));
			delete [] tab->table_units; tab->table_units = new_table_units;
		}
		if (tab->table_names && id_static_char != ADDITIONAL_STATE) {
			char ** new_table_names = (char **)new_struct(N_group*sizeof(char *));
			memcpy(new_table_names, tab->table_names, N_ini*sizeof(char *));
			memcpy(new_table_names+N_group-tab->N_group+N_ini, tab->table_names+N_ini, (tab->N_group-N_ini)*sizeof(char *));
			if (table_names)
				memcpy(new_table_names+N_ini, table_names, (N_group-tab->N_group)*sizeof(char *));
			delete [] tab->table_names; tab->table_names = new_table_names;
		}
		else 	
		tab->table_names = NULL;
		tab->index = N_ini+1;
		tab->N_group = N_group;
	}
}

void del_table_group(Table * tab, int N_group, int N_ini, int id_static_char)
{
	if (tab && N_group > 0 && id_static_char != SPECIAL_STATE) {
		if (N_ini  < 0 || N_ini > tab->N_group) N_ini = tab->N_group;
		int m = max(N_ini-N_group, 0), i;

		N_group = tab->N_group-N_ini+m;
		for ( i = 0; i < tab->N; i++) {
			if (tab->table[i].type < CHAR_TYPE_RECORD); else
			if (tab->table[i].type <  INT_TYPE_RECORD) {
				char ** new_parm_values = (char **)new_struct((N_group+1)*sizeof(char *));
				memcpy(new_parm_values, tab->table[i].parm_values, (m+1)*sizeof(char *));
				memcpy(new_parm_values+m+1, ((char **)tab->table[i].parm_values)+N_ini+1, (tab->N_group-N_ini)*sizeof(char *));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < FLOAT_TYPE_RECORD) {
				int * new_parm_values = (int *)new_struct((N_group+1)*sizeof(int));
				memcpy(new_parm_values, tab->table[i].parm_values, (m+1)*sizeof(int));
				memcpy(new_parm_values+m+1, ((int *)tab->table[i].parm_values)+N_ini+1, (tab->N_group-N_ini)*sizeof(int));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < DOUBLE_TYPE_RECORD) {
				float * new_parm_values = (float *)new_struct((N_group+1)*sizeof(float));
				memcpy(new_parm_values, tab->table[i].parm_values, (m+1)*sizeof(float));
				memcpy(new_parm_values+m+1, ((float *)tab->table[i].parm_values)+N_ini+1, (tab->N_group-N_ini)*sizeof(float));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
			else
			if (tab->table[i].type < NUM_OF_TYPE_RECORD) {
				double * new_parm_values = (double *)new_struct((N_group+1)*sizeof(double));
				memcpy(new_parm_values, tab->table[i].parm_values, (m+1)*sizeof(double));
				memcpy(new_parm_values+m+1, ((double *)tab->table[i].parm_values)+N_ini+1, (tab->N_group-N_ini)*sizeof(double));
				delete [] tab->table[i].parm_values; tab->table[i].parm_values = new_parm_values;
			}
		}
		if (tab->table_units) {
			int * new_table_units = (int *)new_struct(N_group*sizeof(int));
			memcpy(new_table_units, tab->table_units, m*sizeof(int));
			memcpy(new_table_units+m, ((int *)tab->table_units)+N_ini, (tab->N_group-N_ini)*sizeof(int));
			delete [] tab->table_units; tab->table_units = new_table_units;
		}
		if (tab->table_names) {
			if (id_static_char != ADDITIONAL_STATE) {
				char ** new_table_names = (char **)new_struct(N_group*sizeof(char *));
				memcpy(new_table_names, tab->table_names, m*sizeof(char *));
				memcpy(new_table_names+m, ((char **)tab->table_names)+N_ini, (tab->N_group-N_ini)*sizeof(char *));
				delete [] tab->table_names; tab->table_names = new_table_names;
			}
			else 
				for (i = m; i < N_ini; i++)
					tab->table_names[i] = tab->table_names[tab->N_group-N_ini+i];
		}
		if (tab->index > m) tab->index = m+max(tab->index-N_ini, 0);
		tab->N_group = N_group;
	}
}

////////////////////////////////////////////////////
//...перемена порядка группы параметров на обратный;
void inverse_table_group(Table * tab)
{
	if (tab) {
		int N_group2 = tab->N_group/2, j;
		for (int i = 0; i < tab->N; i++) {
			if (tab->table[i].type < CHAR_TYPE_RECORD); else
			if (tab->table[i].type <  INT_TYPE_RECORD) {
				for (j = 1; j <= N_group2; j++)
				swap(((char **)tab->table[i].parm_values)[j],
					   ((char **)tab->table[i].parm_values)[tab->N_group-j+1]);
			}
			else
			if (tab->table[i].type < FLOAT_TYPE_RECORD) {
				for (j = 1; j <= N_group2; j++)
				swap(((int *)tab->table[i].parm_values)[j],
					  ((int *)tab->table[i].parm_values)[tab->N_group-j+1]);
			}
			else
			if (tab->table[i].type < DOUBLE_TYPE_RECORD) {
				for (j = 1; j <= N_group2; j++)
				swap(((float *)tab->table[i].parm_values)[j],
					  ((float *)tab->table[i].parm_values)[tab->N_group-j+1]);
			}
			else
			if (tab->table[i].type < NUM_OF_TYPE_RECORD) {
				for (j = 1; j <= N_group2; j++)
				swap(((double *)tab->table[i].parm_values)[j],
					  ((double *)tab->table[i].parm_values)[tab->N_group-j+1]);
			}
		}
		for (j = 1; tab->table_units && j <= N_group2; j++)
		swap(((int *)tab->table_units)[j-1],
			  ((int *)tab->table_units)[tab->N_group-j]);

		for (j = 1; tab->table_names && j <= N_group2; j++)
		swap(((char **)tab->table_names)[j-1],
			   ((char **)tab->table_names)[tab->N_group-j]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...функция копирования значений из таблицы в другую таблицу (статическую или динамическую);
void table_cpy(Table * tab, Table * tab_source, int id_static_char)
{
	if (tab && tab_source && tab->N <= tab_source->N && tab->N_group <= tab_source->N_group) {
		for (int i = 0; i < tab->N; i++)
		if ( tab->table[i].type == tab_source->table[i].type) {

			if (tab->table[i].type < CHAR_TYPE_RECORD); else
			if (tab->table[i].type <  INT_TYPE_RECORD)
			for (int j = 0; j <= tab->N_group; j ++) {
				if (id_static_char < ADDITIONAL_STATE) {
					char * sss;
					delete_struct(sss = ((char **)tab->table[i].parm_values)[j]);
					if ((sss = ((char **)tab_source->table[i].parm_values)[j]) != NULL) {
						( (char **)tab->table[i].parm_values)[j] = new_struct(((int)strlen(((char **)tab_source->table[i].parm_values)[j])+1)*sizeof(char));

						::strcpy(((char **)tab->table[i].parm_values)[j], ((char **)tab_source->table[i].parm_values)[j]);
					}
				}
//				else
//					((char **)tab->table[i].parm_values)[j] = ((char **)tab_source->table[i].parm_values)[j];
			}
			else
			if (tab->table[i].type < FLOAT_TYPE_RECORD)
				memcpy(tab->table[i].parm_values, tab_source->table[i].parm_values, (tab->N_group+1)*sizeof(int));
			else
			if (tab->table[i].type < DOUBLE_TYPE_RECORD)
				memcpy(tab->table[i].parm_values, tab_source->table[i].parm_values, (tab->N_group+1)*sizeof(float));
			else
			if (tab->table[i].type < NUM_OF_TYPE_RECORD)
				memcpy(tab->table[i].parm_values, tab_source->table[i].parm_values, (tab->N_group+1)*sizeof(double));
		}
		if (tab->table_units && tab_source->table_units)
			memcpy(tab->table_units, tab_source->table_units, tab->N_group*sizeof(int));

      tab->index = (tab_source->index <= tab->N_group ? (tab_source->index < 0 ? 0 : tab_source->index) : tab->N_group);
	}
}

/////////////////////////////////////////////////////////
//...вспомогательная функция образования таблицы образца;
Table * get_table(int N, int N_group, char * table_names[])
{
	if (N_group < 1) return(NULL);
	Table * table = (Table *)new_struct(sizeof(Table)+(N-1)*sizeof(Record));
	if (table) {
		table->N           = N;
		table->N_group     = N_group;
		table->table_names = table_names;
		if ((table->table_units = (int *)new_struct(N_group*sizeof(int))) == NULL)
			delete_table(table); else
			set_table_units(table, UNIT_SI);
	}
	return table;
}

//////////////////////////////////////////////////////////
//...функция порождения статической таблицы по ее шаблону;
Table * get_shablon_table(int N, int N_group, Shablon * records, char * table_names[], int id_static_char)
{
	Table * table = get_table(N, N_group, table_names);
	if (table)
	for (int i_last = 0, i = 0; i < table->N; i++) {
		if (records[i_last].parm_type < CHAR_TYPE_RECORD); else
		if (records[i_last].parm_type <  INT_TYPE_RECORD) {
				char ** parm_values = id_static_char == SPECIAL_STATE ? NULL : (char **)new_struct((table->N_group+1)*sizeof(char *));
				set_table_param(table, i+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
		}
		else
		if (records[i_last].parm_type < FLOAT_TYPE_RECORD) {
				int * parm_values = id_static_char == SPECIAL_STATE ? NULL : (int *)new_struct((table->N_group+1)*sizeof(int));
				set_table_param(table, i+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
		}
		else
		if (records[i_last].parm_type < DOUBLE_TYPE_RECORD) {
				float * parm_values = id_static_char == SPECIAL_STATE ? NULL : (float *)new_struct((table->N_group+1)*sizeof(float));
				set_table_param(table, i+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
		}
		else
		if (records[i_last].parm_type < NUM_OF_TYPE_RECORD) {
				double * parm_values = id_static_char == SPECIAL_STATE ? NULL : (double *)new_struct((table->N_group+1)*sizeof(double));
				set_table_param(table, i+1, records[i_last].parm_type, records[i_last].parm_name, parm_values);
		}
		if (i+1 < table->N && records[i_last+1].parm_type != ERR_TYPE_RECORD) i_last = i+1;
	}
	return table;
}

/////////////////////////////////////////////////////////////////////////////////
//...функция порождения таблицы (статической или динамической) по другой таблице;
Table * get_shablon_table(Table * table, int id_static_char)
{
	Table * tab_shablon = NULL;
	if (table) {
		Shablon * records = (Shablon *)new_struct(table->N*sizeof(Shablon));
		if (records) {
			for (int i = 0; i < table->N; i++) {
				char * parm_name = table->table[i].parm_name;
				if (id_static_char == NULL_STATE && table->table[i].parm_name) {
					parm_name = (char *)new_struct(((int)strlen(table->table[i].parm_name)+1)*sizeof(char));
					::strcpy(parm_name, table->table[i].parm_name);
				}
				records[i].parm_name = parm_name;
				records[i].parm_type = (Num_type_Record)table->table[i].type;
			}
			char ** table_names = table->table_names;
			if (id_static_char < ADDITIONAL_STATE && table->table_names) {
				table_names = (char **)new_struct(table->N_group*sizeof(char *));
				for (int j = 0; table_names && j < table->N_group; j++) {
					table_names[j] = table->table_names[j];
					if (id_static_char == NULL_STATE && table->table_names[j]) {
							table_names[j] = (char *)new_struct(((int)strlen(table->table_names[j])+1)*sizeof(char));
							::strcpy(table_names[j], table->table_names[j]);
					}
				}
			}
			tab_shablon = get_shablon_table(table->N, table->N_group, records, table_names, id_static_char);
		}
		delete[] records;
	}
	return(tab_shablon);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка индекса таблицы в контексте
void SetTableIndex(void * context, int index)
{
  Context * cont = (Context *)context;
  if (cont && cont->table) {
	  index = min(max(0, index), cont->table->N_group);

      cont->table->index = index;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка значения параметра в таблице контекста
void SetTableParam(Table * table, int index, int N, double f, int id_static_char)
{
  if (table && index >= 0 && index <= table->N_group && N >= 0 && N < table->N) {
      Record * rec = &table->table[N];
		if (rec->type < CHAR_TYPE_RECORD); else
		if (rec->type <  INT_TYPE_RECORD) {
			if (id_static_char < ADDITIONAL_STATE) {
				delete_struct(((char **)rec->parm_values)[index]);
				char * sss = (char *)(long long)f;
				if (sss) {
					((char **)rec->parm_values)[index] = new_struct(((int)strlen(sss)+1)*sizeof(char));
					::strcpy (((char **)rec->parm_values)[index], sss);
				}
			}
//			else
//				((char **)rec->parm_values)[index] = (char *)((unsigned long)f);
		}
		else
		if (rec->type < FLOAT_TYPE_RECORD)
			((int *)rec->parm_values)[index] = (int)f; 
		else
		if (rec->type < DOUBLE_TYPE_RECORD)
			((float *)rec->parm_values)[index] = (float)f; 
		else
		if (rec->type < NUM_OF_TYPE_RECORD)
			((double *)rec->parm_values)[index] = (double)f; 
  }
}

void SetTableParam(void * context, int index, int N, double f)
{
  Context * cont = (Context *)context;
  if (cont) SetTableParam(cont->table, index, N, f, cont->static_char);
}

void SetTableParam(Table * table, int N, double f, int id_static_char)
{
  if (table) SetTableParam(table, table->index, N, f, id_static_char);
}

void SetTableParam(void * context, int N, double f)
{
  Context * cont = (Context *)context;
  if (cont && cont->table)
		SetTableParam(context, cont->table->index, N, f);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...извлечение значения параметра из таблицы контекста (произвольный и текущий индекс)
double GetTableParam(Table * table, int index, int N)
{
  if (table) {
      index = (index <= table->N_group ? (index < 0 ? 0 : index) : table->N_group);
      N = (N < table->N ? (N < 0 ? 0 : N) : table->N-1);

      Record * rec = &table->table[N];
		if (rec->type < CHAR_TYPE_RECORD); else
		if (rec->type <  INT_TYPE_RECORD)
			return (double)((long long)((char **)rec->parm_values)[index]);
		else
		if (rec->type < FLOAT_TYPE_RECORD)
			return (double)((int *)rec->parm_values)[index];
		else
		if (rec->type < DOUBLE_TYPE_RECORD)
			return (double)((float *)rec->parm_values)[index];
		else
		if (rec->type < NUM_OF_TYPE_RECORD) {
			double fff = (double)((double *)rec->parm_values)[index];
			return (double)((double *)rec->parm_values)[index];
		}
  }
  return(0.);
}

double GetTableParam(void * context, int index, int N)
{
  Context * cont = (Context *)context;
  if (cont) return(GetTableParam(cont->table, index, N));
  else      return(0.);
}

double GetTableParam(Table * table, int N)
{
  if (table) return(GetTableParam(table, table->index, N));
  else		 return(0.);
}

double GetTableParam(void * context, int N)
{
  Context * cont = (Context *)context;
  if (cont  && cont->table) 
		return(GetTableParam(context, cont->table->index, N));

  return(0.);
}

int GetTableIndex(void * context)
{
  Context * cont = (Context *)context;
  if (cont  && cont->table) 
		return(cont->table->index);

  return(0);
}

///////////////////////////////////////////////////////////////////
//...извлечение числа параметров и числа груп параметров в таблице;
int GetNParam(void * context)
{
  Context * cont = (Context *)context;
  if (cont && cont->table) {
      return  cont->table->N;
  }
  return(0);
}

int GetNGroup(void * context)
{
  Context * cont = (Context *)context;
  if (cont && cont->table) {
      return  cont->table->N_group;
  }
  return(0);
}

//////////////////////////////////////////////////////////
//...загрузка параметров таблицы в позицию "по-умолчанию";
void LoadingParameters(Table * table, int index)
{
	Record * rec;
	char   * str;
	if (table) {
      index = min(max(0, index), table->N_group);
   	for (int len, i = 0; i < table->N; i++) {
   		rec = &table->table[i];
   		if (rec->type < CHAR_TYPE_RECORD); else
   		if (rec->type <  INT_TYPE_RECORD) {
   			delete_struct(((char **)rec->parm_values)[0]);
			if ((str = ((char **)rec->parm_values)[index]) != NULL) {
   				len = (int)strlen(str);

   				((char **)rec->parm_values)[0] = (char *)new_struct((len+1)*sizeof(char));
   				memcpy(((char **)rec->parm_values)[0], str, len*sizeof(char));
   				str =  ((char **)rec->parm_values)[0];
   			}
   		}
   		else
			if (rec->type < FLOAT_TYPE_RECORD)
				((int *)rec->parm_values)[0] = ((int *)rec->parm_values)[index];
   		else
   		if (rec->type < DOUBLE_TYPE_RECORD); else
   		if (rec->type < NUM_OF_TYPE_RECORD)
				((double *)rec->parm_values)[0] = ((double *)rec->parm_values)[index];
   	}
	}
}

void LoadingParameters(void * context, int index)
{   
	Table * tab = (Table *)GetTable(context);
	LoadingParameters(tab, index);

	SetTableIndex(context, 0);
	SetUnits     (context, tab ? tab->table_units[index-1] : UNIT_SI);
}

/////////////////////////////////////////////////
//...обнуление параметров таблицы "по-умолчанию";
void ZeroCustomParameters(Table * table)
{
	Record * rec;
	for (int i = 0; table && i < table->N; i++) {
   	rec = &table->table[i];
   	if (rec->type < CHAR_TYPE_RECORD); else
   	if (rec->type <  INT_TYPE_RECORD)
   		delete_struct(((char **)rec->parm_values)[0]);
   	else
		if (rec->type < FLOAT_TYPE_RECORD)
			((int *)rec->parm_values)[0] = 0;
   	else
   	if (rec->type < DOUBLE_TYPE_RECORD); else
   	if (rec->type < NUM_OF_TYPE_RECORD)
			((double *)rec->parm_values)[0] = 0.;
   }
}
void ZeroCustomParameters(void * context){   
	ZeroCustomParameters((Table *)GetTable(context));
}

///////////////////////////////////////////////////////////////
//...копирование параметров "по-умолчанию" в группу параметров;
void CopyTableParameters(Table * table, int index, int index_ini, int N_group)
{   
	if (! table) return;

	if (index < index_ini) {
		for (int j = 0; j < N_group; j++)
		for (int i = table->N-1; i >= 0; i--)
			SetTableParam(table, index+j, i, GetTableParam(table, index_ini+j, i));
	}
	else {
		for (int j = N_group-1;  j >= 0; j--)
		for (int i = table->N-1; i >= 0; i--)
			SetTableParam(table, index+j, i, GetTableParam(table, index_ini+j, i));
	}
}

void CopyTableParameters(void * context, int index, int index_ini, int N_group)
{   
	CopyTableParameters((Table *)GetTable(context), index, index_ini, N_group);
}

///////////////////////////////////////////////////////////////////////////////////
//...добавление параметров "по-умолчанию" в группу параметров (следом за индексом);
void AddCustomParameters(Table * table, char * table_name, int index, int id_static_char)
{  
	add_table_group(table, 1, table_name ? & table_name : NULL, index, id_static_char);
	CopyTableParameters(table, index+1);
}

void AddCustomParameters(void * context, char * table_name, int index)
{  
	AddCustomParameters((Table *)GetTable(context), table_name, index, GetIdStaticChar(context));
	CopyTableParameters(context, index+1);
}

void AddTableParameters(void * context, int N_group, char ** table_names, int index)
{   
	add_table_group((Table *)GetTable(context), N_group, table_names, index, GetIdStaticChar(context));
}

///////////////////////////////////////////////////////
//...преобразование статической таблицы в динамическую;
void ConvertStaticTable(void * context, int id_static_char)
{   
	Table * tab = (Table *)GetTable(context);
	int static_char = GetIdStaticChar(context);
	Record * rec;											

	if (id_static_char == SPECIAL_STATE || static_char <= id_static_char || ! tab) return;
	if (id_static_char <  SPECIAL_STATE) { //,,,распределяем оболочку строковых переменных;
		for (int i = 0; i < tab->N; i++) {
			rec = &tab->table[i];
			if (ERR_TYPE_RECORD < rec->type && rec->type < INT_TYPE_RECORD && rec->parm_values) {
				char ** new_parm_values = (char **)new_struct((tab->N+1)*sizeof(char *));
				memcpy (new_parm_values, rec->parm_values, (tab->N+1)*sizeof(char *));
				rec->parm_values = new_parm_values;
			}
		}
		((Context *)context)->static_char = ADDITIONAL_STATE;
	}
	if (id_static_char < ADDITIONAL_STATE) { //,,,распределяем оболочку строковых имен групп параметров;
		if (tab->table_names) {
			char ** new_table_names = (char **)new_struct(tab->N_group*sizeof(char *));
			memcpy (new_table_names, tab->table_names, tab->N_group*sizeof(char *));
			tab->table_names = new_table_names;
		}
		((Context *)context)->static_char = OK_STATE;
	}
	if (id_static_char < OK_STATE) { //,,,распределяем строковые имена и переменные;
		for (int i = 0; i < tab->N; i++) {
			rec = &tab->table[i];
			if (ERR_TYPE_RECORD < rec->type && rec->type < INT_TYPE_RECORD && rec->parm_values) {
				for (int j = 1; j <= tab->N_group; j++)
					if (((char **)rec->parm_values)[j]) {
						int m = (int)strlen(((char **)rec->parm_values)[j]);
						char * new_parm_value = (char *)new_struct((m+1)*sizeof(char));
						memcpy(new_parm_value, ((char **)rec->parm_values)[j], m*sizeof(char));
						((char **)rec->parm_values)[j] = new_parm_value;
					}
			}
			if (rec->parm_name) { //,,,строковые имена параметров;
				int m = (int)strlen(rec->parm_name);
				char * new_parm_name = (char *)new_struct((m+1)*sizeof(char));
				memcpy(new_parm_name, rec->parm_name, m*sizeof(char));
				rec->parm_name = new_parm_name;
			}
		}
		if (tab->table_names) { //,,,строковые имена групп параметров;
			for (int i = 0; i < tab->N_group; i++)
				if (tab->table_names[i]) {
					int m = (int)strlen(tab->table_names[i]);
					char * new_table_name = (char *)new_struct((m+1)*sizeof(char));
					memcpy(new_table_name, tab->table_names[i], m*sizeof(char));
					tab->table_names[i] = new_table_name;
				}
		}
		if (((Context *)context)->sample_name) { //,,,строковое именя образца;
			int m = (int)strlen(((Context *)context)->sample_name);
			char * new_sample_name = (char *)new_struct((m+1)*sizeof(char));
			memcpy(new_sample_name, ((Context *)context)->sample_name, m*sizeof(char));
			((Context *)context)->sample_name = new_sample_name;
		}
		if (((Context *)context)->GOST_name) { //,,,строковое именя ГОСТа образца;
			int m = (int)strlen(((Context *)context)->GOST_name);
			char * new_GOST_name = (char *)new_struct((m+1)*sizeof(char));
			memcpy(new_GOST_name, ((Context *)context)->GOST_name, m*sizeof(char));
			((Context *)context)->GOST_name = new_GOST_name;
		}
		((Context *)context)->static_char = NULL_STATE;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...concatenate string to the input string of the table in a sample (working procedure)
void str_reinst(char *& str, char * str_concat, int & buf_count, int id_reinst, int buf_length)
{
	int  k = (int)user_strlen(str), m = (int)user_strlen(str_concat);
	while (k+m+1 > buf_count) {
		char * new_str = new_struct((buf_count += buf_length)*sizeof(char));

		if (str) ::strcpy(new_str, str);
		delete_struct(str); str = new_str;	
	}
	if (id_reinst) str[0] = '\x0';

	if (str_concat) strcat(str, str_concat);
}

void str_reinst(char *& str, const char * str_concat, int & buf_count, int id_reinst, int buf_length)
{
	int  k = (int)user_strlen(str), m = (int)user_strlen(str_concat);
	while (k+m+1 > buf_count) {
		char * new_str = new_struct((buf_count += buf_length)*sizeof(char));

		if (str) ::strcpy(new_str, str);
		delete_struct(str); str = new_str;	
	}
	if (id_reinst) str[0] = '\x0';

	if (str_concat) strcat(str, str_concat);
}

//////////////////////////////////////////////////////////////////////////////////////
//...запись выделенной строки параметров в таблицу (метод для записи плавающих строк);
char * SetTableParamsAsString(Table * table, int index, char * pchar, char *& p_end, const char * ColonStr, const char * SemiColonStr)
{
	char * pnext, * pname, * plast;
	int lenColon = (int)user_strlen(ColonStr), lenSemiColon = (int)user_strlen(SemiColonStr), lenStrSeparator = (int)user_strlen(STRSEPARATOR), lenName, i = -1;

	if (pchar)
	while ( pchar < p_end && ++i < table->N &&
		  ((pnext = strstr(pchar, SemiColonStr)) != NULL && pnext+lenSemiColon <= p_end ||
		   (pnext = p_end) != NULL)) {
		char Temp0 = pnext[0];
		pnext[0]   = '\x0';

//////////////////////////////////////////////////////////////////////////////
//...выделяем имя параметра и заносим выделенное значение параметра в таблицу;
		if ((pname = strstr(pchar, ColonStr)) != NULL && pname+lenColon <= pnext ||
			(pname = pnext) != NULL) {
			char Temp2 = pname[0];
			pname[0] = '\x0';
			if (! table->table[i].parm_name && (lenName = (int)strlen(pchar)) > 0) { // pchar -- имя параметра;
				table->table[i].parm_name = (char *)new_struct((lenName+1)*sizeof(char));
				memcpy(table->table[i].parm_name, pchar, lenName*sizeof(char));
			}
			pname[0] = Temp2;
			if ((pchar = pname) < pnext) pchar += lenColon;
			if (table->table[i].type < CHAR_TYPE_RECORD); else
			if (table->table[i].type <  INT_TYPE_RECORD) {
				delete_struct(((char **)(table->table[i].parm_values))[index]);
				if (table->table[i].type != COMMENT_CHAR_TYPE_RECORD) {
					if ( pchar < pnext && strncmp(pchar, STRSEPARATOR, lenStrSeparator) && pchar[0] != '\xD') {
						((char **)table->table[i].parm_values)[index] = (char *)new_struct(((int)strlen(pchar)+1)*sizeof(char));
						strcat(((char **)(table->table[i].parm_values))[index], pchar);
					}
				}
				else {
					int m = 0; pname = pchar;
					while ((plast = strchr(pname, '\"')) != NULL) {
						pname = plast+1; m++;
					}
					while (pname && m%2) {//...ищем конец комментария;
						pnext[0] = Temp0;
						if ((plast = strchr(pname, '\"')) == NULL) pname = plast;
						else {
							pname = plast+1; m++;
							if (plast[1] == '\"') {
								pname = plast+2; m++;
							}
							pnext = strstr(pname, SemiColonStr);
						}
						if ((p_end = strstr(pnext, STRSEPARATOR)) == NULL)
							 p_end = pnext+strlen(pnext);
						else p_end += lenStrSeparator;

						Temp0 = pnext[0];	pnext[0] = '\x0';
					}
					if ( pchar < pnext-1) {//...записываем комментарий (в том числе и нулевой длины);
						((char **)table->table[i].parm_values)[index] = (char *)new_struct((m = (int)strlen(pchar)-1)*sizeof(char));
						memcpy(pname = ((char **)(table->table[i].parm_values))[index], pchar+1, m-1);

						while ((plast = strchr(pname, '\"')) != NULL && plast < pname+m-2) {
							memmove(plast, plast+1, m -= (int)(plast-pname)+1); 
							pname = plast+1;
						}
					}
				}
			}
			else
			if (table->table[i].type < FLOAT_TYPE_RECORD) {
				int i_value;
				if (sscanf(pchar, "%i", &i_value) == 1)
					((int *)(table->table[i].parm_values))[index] = i_value; else
					((int *)(table->table[i].parm_values))[index] = 0;
			}
			else
			if (table->table[i].type < DOUBLE_TYPE_RECORD) {
				char * ppp = pchar;
				while ((ppp = strstr(ppp, ".")) != NULL) { ppp[0] = localeconv()->decimal_point[0]; ppp++;}

				float f_value;
				if (sscanf(pchar, "%g", &f_value) == 1)
					((float *)(table->table[i].parm_values))[index] = f_value; else
					((float *)(table->table[i].parm_values))[index] = 0.f;

				while ((pchar = strstr(pchar, localeconv()->decimal_point)) != NULL) { pchar[0] = '.'; pchar++;}
			}
			else
			if (table->table[i].type < NUM_OF_TYPE_RECORD) {
				char * ppp = pchar;
				while ((ppp = strstr(ppp, ".")) != NULL) { ppp[0] = localeconv()->decimal_point[0]; ppp++;}
	
				double d_value;
				if (sscanf(pchar, "%lg", &d_value) == 1)
					((double *)(table->table[i].parm_values))[index] = d_value; else
					((double *)(table->table[i].parm_values))[index] = 0.;

				while ((pchar = strstr(pchar, localeconv()->decimal_point)) != NULL) { pchar[0] = '.'; pchar++;}
			}
		}
		pnext[0] = Temp0;
		if ((pchar = pnext) < p_end) pchar += lenSemiColon;
	}
	return(pchar);
}

////////////////////////////////////////////////////////////////////////////////////////
//...запись выделенной строки параметров в таблицу (метод для записи константных строк);
const char * SetTableParamsAsString(Table * table, int index, const char * strParams, char *& strWork, int & buf_count, const char * ColonStr, const char * SemiColonStr)
{
	int posColon, posSemiColon, lenColon = (int)user_strlen(ColonStr), lenSemiColon = (int)user_strlen(SemiColonStr), shift = 0;
	char * str, * strTemp;

	for (int i = 0; i < table->N; i++) {
		str_reinst(strWork, table->table[i].parm_name, buf_count, 1);

		strTemp = (char *)(unsigned long *)strstr(strParams+shift, strWork);
		str_reinst(strWork, SemiColonStr, buf_count, 1);

		str = strstr(strTemp, strWork); posSemiColon  = (int)(str-strParams);
		str_reinst(strWork, ColonStr, buf_count, 1);

		str = strstr(strTemp, strWork); posColon = (int)(str-strParams);
		str_reinst(strWork, strParams+posColon+lenColon, buf_count, 1);

		strWork[posSemiColon-posColon-lenColon] = '\x0';
		shift = posSemiColon+lenSemiColon; 

		if (table->table[i].type < CHAR_TYPE_RECORD); else
		if (table->table[i].type <  INT_TYPE_RECORD) {
			delete_struct(((char **)(table->table[i].parm_values))[index]);
			if (posSemiColon-posColon-lenColon > 0) {
				((char **)table->table[i].parm_values)[index] = (char *)new_struct((posSemiColon-posColon-lenColon+1)*sizeof(char));
				strcat(((char **)(table->table[i].parm_values))[index], strWork);
			}
		}
		else
		if (table->table[i].type < FLOAT_TYPE_RECORD) {
			int i_value;
			if (sscanf(strWork, "%i", &i_value) == 1)
				((int *)(table->table[i].parm_values))[index] = i_value; else
				((int *)(table->table[i].parm_values))[index] = 0;
		}
		else
		if (table->table[i].type < DOUBLE_TYPE_RECORD) {
			float f_value;
			if (sscanf(strWork, "%g", &f_value) == 1)
				((float *)(table->table[i].parm_values))[index] = f_value; else
				((float *)(table->table[i].parm_values))[index] = 0.f;
		}
		else
		if (table->table[i].type < NUM_OF_TYPE_RECORD) {
			double d_value;
			if (sscanf(strWork, "%lg", &d_value) == 1)
				((double *)(table->table[i].parm_values))[index] = d_value; else
				((double *)(table->table[i].parm_values))[index] = 0.;
		}
	}
	return(strParams+shift);
}

////////////////////////////////////////////////////////////////////////
//...преобразование строки параметров таблицы в текстовое представление;
void GetTableParamsAsString(Table * table, int index, char *& str, int & buf_count, const char * ColonStr, const char * SemiColonStr)
{
	char tempStr[201];
	str_reinst(str, (char *)NULL, buf_count, 1);

	for (int i = 0; i < table->N; i++) {
		str_reinst(str, table->table[i].parm_name, buf_count);
		str_reinst(str, ColonStr, buf_count);

		if (table->table[i].type < CHAR_TYPE_RECORD); else
		if (table->table[i].type <  INT_TYPE_RECORD) {
			if (table->table[i].type != COMMENT_CHAR_TYPE_RECORD)
				str_reinst(str, ((char **)(table->table[i].parm_values))[index], buf_count);
			else {
				char * pchar = ((char **)(table->table[i].parm_values))[index], * pnext;
				if (pchar) {
					str_reinst(str, "\"", buf_count);
					while ((pnext = strchr(pchar, '\"')) != NULL) {
						pnext[0] = '\x0';
						str_reinst(str, pchar,  buf_count);
						str_reinst(str, "\"\"", buf_count);
						pnext[0] = '\"';
						pchar = pnext+1;
					}
					str_reinst(str, pchar, buf_count);
					str_reinst(str, "\"",  buf_count);
				}
			}
		}
		else
		if (table->table[i].type < FLOAT_TYPE_RECORD)
			str_reinst(str, _itoa(((int *)(table->table[i].parm_values))[index], tempStr, 10), buf_count);
		else
		if (table->table[i].type < DOUBLE_TYPE_RECORD) {
			sprintf(tempStr, "%g", ((float *)(table->table[i].parm_values))[index]);

			char * ppp = tempStr;
			while ((ppp = strstr(ppp, localeconv()->decimal_point)) != NULL) { ppp[0] = '.'; ppp++;}

			str_reinst(str, tempStr, buf_count);
		}
		else
		if (table->table[i].type < NUM_OF_TYPE_RECORD) {
			sprintf(tempStr, "%0.15lg", ((double *)(table->table[i].parm_values))[index]);

			char * ppp = tempStr;
			while ((ppp = strstr(ppp, localeconv()->decimal_point)) != NULL) { ppp[0] = '.'; ppp++;}

			str_reinst(str, tempStr, buf_count);
		}
		str_reinst(str, SemiColonStr, buf_count);
  }
}

///////////////////////////////
//...сохранение версии формата;
void version_out(FILE * device, int wrongVersion)
{
	if (device) {
		char buff[201]; ::strcpy(buff,  VERSIONSTR); //...запись версии формата проектируемой линии;
		_ultoa(ID_VERSION+wrongVersion, buff+strlen(VERSIONSTR), 10);
		
		fwrite(buff, sizeof(char), strlen(buff), device);
		fwrite(STRSEPARATOR, sizeof(char), strlen(STRSEPARATOR), device);
	}
}

/////////////////////////////////////////////
//...сохранение таблицы из контекста в файле;
void table_out(FILE * device, void * context, int head_writing, int id_csv_state)
{
  Context * cont = (Context *)context;
  if (cont && device) {
		char * str = NULL, * strWork = NULL, buff[201];
		int buf_count = 0, lenStrSeparator = (int)strlen(STRSEPARATOR), i;

		if (head_writing) {
			strWork = GetSampleComment(context);
			if (strWork) 
				fwrite (strWork, sizeof(char), strlen(strWork), device);
			fwrite (STRSEPARATOR, sizeof(char), lenStrSeparator, device);

			strWork = GetSampleDescription(context);
			if (strWork) 
				fwrite (strWork, sizeof(char), strlen(strWork), device);
			fwrite (STRSEPARATOR, sizeof(char), lenStrSeparator, device);
		}

		if (cont->table) {
			sprintf(buff, "%s%i", TABLESTR, cont->table->N_group);
			fwrite (buff, sizeof(char), strlen(buff), device);
			fwrite ("\n", sizeof(char), 1, device);

			for (i = 1; i <= cont->table->N_group; i++)
			{
				switch (id_csv_state) {
				case	  NULL_STATE: GetTableParamsAsString(cont->table, i, str, buf_count); break;
				case		 OK_STATE: GetTableParamsAsString(cont->table, i, str, buf_count, " "); break;
				case SPECIAL_STATE: int N_temp = cont->table->N, N_points = ((int *)cont->table->table[1].parm_values)[i];
					cont->table->N = N_points < 0 ? 2-N_points*2 : 2+N_points;
					GetTableParamsAsString(cont->table, i, str, buf_count); 
					cont->table->N = N_temp;
				break;
				}
				fwrite (str, sizeof(char), strlen(str), device);

				strWork = GetProfileName(context, i);
				if (strWork) 
					fwrite (strWork, sizeof(char), strlen(strWork), device);
				fwrite (STRSEPARATOR, sizeof(char), lenStrSeparator, device);
			}
			delete_struct(str);
		}
		else {
			sprintf(buff, "%s%i", TABLESTR, 0);
			fwrite (buff, sizeof(char), strlen(buff), device);
			fwrite (STRSEPARATOR, sizeof(char), lenStrSeparator, device);
		}
	}
}

int TableOut(char * ch_TABLE, void * context, int head_writing)
{
	FILE * device = fopen(ch_TABLE, "w+b");
	table_out(device, context, head_writing);

	fclose(device);
	return(device != 0);
}

int TableOut(const char * ch_TABLE, void * context, int head_writing)
{
	FILE * device = fopen(ch_TABLE, "w+b");
	table_out(device, context, head_writing);

	fclose(device);
	return(device != 0);
}

///////////////////////////////////////////////////////////////
//...шаблон CSV таблицы для зачитывания и записи Excel формата;
Table * GetCSVTable(int N_group, int N)
{
	Shablon records[] = {{CHAR_TYPE_RECORD, NULL},//(0) -- First position и т.д.
								{ ERR_TYPE_RECORD, NULL},//(1) -- признак продолжения...
	};
	Table * table = get_shablon_table(N, N_group, records, NULL);
	return table;
}

////////////////////////////////////////////////////////
//...шаблон таблиц для описания образцов из модуля BARs;
Table * GetBARsTable(int N_group, int N)
{
	Shablon records[] = {{	 INT_TYPE_RECORD, NULL},//(0) -- First position
								{DOUBLE_TYPE_RECORD, NULL},//(1) -- вещественный параметр и т.д.
								{	 ERR_TYPE_RECORD, NULL},//(2) -- признак продолжения...
	};
	Table * table = get_shablon_table(N, N_group, records, NULL);
	return table;
}

// ======================================================================================================
// Version reading  ========================================================================================
// ======================================================================================================
unsigned long VersionReading(char * pchar)
{
	char * p_end, * pnext, 
		  * p_eof = pchar+strlen(pchar);
	unsigned long version = 0;

///////////////////////////////////
//...зачитываем код версии формата;
	if ((p_end = strstr(pchar, VERSIONSTR)) != NULL &&
	   ((pnext = strstr(p_end, STRSEPARATOR)) != NULL || (pnext = p_eof) != NULL)) {
			char Temp0 = pnext[0];
			pnext[0]   = '\x0';

			version = strtoul(p_end+strlen(VERSIONSTR), NULL, 10);

			pnext[0] = Temp0;
	}
	return(version);
}

// ======================================================================================================
// CSV table reading  ===================================================================================
// ======================================================================================================
void TableCSVReading(void * context, char *& pchar, char * p_end)
{
	char * pnext = pchar;
	int i, N_group;

///////////////////////////////////////
//...вычисляем число строчек в таблице;
	N_group = 0;
	while ((pnext = strstr(pnext, "\xD\xA")) != NULL && (pnext += 2) <= p_end) N_group++;

//////////////////////////////////////////////
//...порождение новой таблицы для csv-формата;
	Table * table = GetCSVTable(N_group, 12);
//	table->table[9].type = COMMENT_CHAR_TYPE_RECORD;
	if (table) {
		SetTable(context, table);

//////////////////////////////////////////////
//...заполнение новой таблицы данными с диска;
		i = 0;
		while ((pnext = strstr(pchar, "\xD\xA")) != NULL && (pnext += 2) <= p_end) {
				SetTableParamsAsString(table, ++i, pchar, pnext);
				pchar = pnext;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...уничтожение только образца внутри контекста
void DeleteSample(void * context)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm) {
      delete (cont->sm); cont->sm = NULL;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...извлечение образца изнутри контекста
void * GetSample(void * context)
{
  Context * cont   = (Context *)context;
  void    * sample = NULL; 
  if (cont) {
    sample   = (void *)cont->sm;
    cont->sm = NULL;
  }
  return sample;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка существующего образца внутрь контекста
#ifndef ___ABRIDGE_PROFILE_MODE___
void SetSample(void * context, void * sample)
{
  Context * cont = (Context *)context;
  if (cont) {
  if (cont->sm) delete(cont->sm);
      cont->sm = (CSample *)sample;
  }
  return;
}

/////////////////////////////////////////
//...запись образца из контекста на диск;
void OutSample(void * context, FILE * device)
{
  Context * cont = (Context *)context;
  if (cont) {
	  cont->sm->sample_out(device);
  }
}

/////////////////////////////////////////
//...извлечение образца из его контекста;
CSample * GetContextSample(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->sm);
  else      return(NULL);
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция, устанавливающая параметры по умолчанию в контекст
void set_default(Context * cont)
{
  if (cont) {
//		cont->units    = UNIT_SI;
		cont->units    = UNIT_CGS;
		cont->E        = 1.;
		cont->nju      = 0.3;
		cont->Ro       = 1.225;
		cont->Hz       = 0.;
		cont->velocity = 340.;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка массива параметров в контекст
void SetContextParam(void * context, Param * param, int id_free)
{
  Context * cont = (Context *)context;
  if (cont) {
      if (id_free) delete_struct(cont->param); cont->param = param;
  }
  return;
}

Param * GetContextParam(void * context)
{
  Context * cont = (Context *)context;
  return cont ? cont->param : NULL;
}

void SetSampleNumber(void * context, int Num)
{
  Context * cont = (Context *)context;
  if (cont) cont->N = Num;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка параметров материала в контекст
int SetMaterial(void * context, double E, double nju, double Ro)
{
  Context * cont = (Context *)context;
  if (cont && E > EE && 0. <= nju && nju < 0.5-EE) {
      cont->E   = E;
      cont->nju = nju;
      cont->Ro  = Ro;

      return(1);
  }
  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка скорости звука в акустической среде в контекст
int SetVelocity(void * context, double velocity)
{
  Context * cont = (Context *)context;
  if (cont && velocity >= 0.) {
      cont->velocity    = velocity;
      return(1);
  }
  return(0);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка системы единиц в контекст
void SetUnits(void * context, int units)
{
  Context * cont = (Context *)context;
  if (cont && 0 <= units && units < NUM_UNITS) {
      cont->units = units;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка признака замыкания образца в контекст
void SetConcat(void * context, int concat)
{
  Context * cont = (Context *)context;
  if (cont) {
      cont->concat = concat;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка признака типа используемой блочной структуры
void SetIdStruc(void * context, int id_struc)
{
  Context * cont = (Context *)context;
  if (cont) {
      cont->id_struc = id_struc;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...определение метода решения блочной системы уравнений в контексте
void SetMethod(void * context, int method, double * param)
{
  Context * cont = (Context *)context;
  if (cont) {
      cont->method = method;
      if (param) switch (cont->method) {
          case 1: cont->Om  =      param[0];
                  cont->Eps =      param[1];
                  cont->Do  = (int)param[2];
          break;
      }
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установление параметров метода покрытия блоками параллелепипедного представления
void SetGridRect(void * context, double h, double h_min, double h_med, int id_struc)
{
  Context * cont = (Context *)context;
  if (cont) {
      cont->h        = h;
      cont->h_min    = h_min;
      cont->h_med    = h_med;
      cont->id_struc = id_struc;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...установка окаймляющего параллелепипеда в контекст
void SetFrameContext(void * context, double X1, double X2,
                                     double Y1, double Y2, double Z1, double Z2)
{
  Context * cont = (Context *)context;
  if (cont) {
      cont->left   = X1;
      cont->right  = X2;
      cont->bottom = Y1;
      cont->top    = Y2;
      cont->back   = Z1;
      cont->front  = Z2;
  }
  return;
}

/////////////////////////////////////////////////////
//...извлечение характеристик материала из контекста;
void GetMaterial(void * context, double & E, double & nju, double & Ro)
{
  Context * cont = (Context *)context;
  if (cont) {
      E   = cont->E;
      nju = cont->nju;
      Ro  = cont->Ro;
  }
  return;
}

//////////////////////////////////////////////////////////
//...извлечение признака статической таблицы из контекста;
int GetIdStaticChar(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->static_char;
  }
  return(ERR_STATE);
}

////////////////////////////////////////////
//...извлечение системы единиц из контекста;
int GetUnits(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->units;
  }
  return(0);
}

///////////////////////////////////////////
//...извлечение признака замыкания образца;
int GetConcat(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->concat;
  }
  return(0);
}

/////////////////////////////////////
//...извлечение границ рамки образца;
double Get2DLeft(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->left;
  }
  return(0.);
}

double Get2DRight(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->right;
  }
  return(0.);
}

double Get2DBottom(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->bottom;
  }
  return(0.);
}

double Get2DTop(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->top;
  }
  return(0.);
}

double Get3DNear(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->front;
  }
  return(0.);
}

double Get3DFar(void * context)
{
  Context * cont = (Context *)context;
  if (cont) {
      return cont->back;
  }
  return(0.);
}

#ifdef ___PRO_SHEA_cpp___
///////////////////////////////////////////////////////////////////////////////////////////////////
//...количество связных компонент контура образца, извлекаемого из контекста
  int Sample2DContourLinks(void * context)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm) return max(cont->sm->N_link, 1);
  else                  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...количество точек, описывающих одну компоненту контуpа всего образца, извлекаемого из контекста
  int Sample2DContourPoints(void * context, int id_link, int Max_knee)
{
  Context * cont = (Context *)context;
  if (cont) {
      return (sample_contour_points(cont->sm, id_link, Max_knee));
  }
  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...описание одной компоненты контуpа всего образца (массив типа double), извлекаемое из контекста
  void Sample2DContour(void * context, int id_link, int Max_knee, double * X, double * Y)
{
  Context * cont = (Context *)context;
  if (cont) {
      int count = 0;
      sample_contour(cont->sm, id_link, Max_knee, count, X, Y);
  }
  return;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//...извлечение имени обpазца, его номера, таблицы, описания и степени решенности из контекста
char * GetSampleName(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->sample_name);
  else      return(NULL);
}

char * GetGOSTName(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->GOST_name);
  else      return(NULL);
}

char * GetProfileName(void * context, int index)
{
  Context * cont = (Context *)context;
  if (cont && cont->table && cont->table->table_names) {
	  index = min(max(1, index), cont->table->N_group);

	  return(cont->table->table_names[index-1]);
  }
  else return(NULL);
}


int GetSampleNumber(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->N);
  else      return(NUM_SAMPLE_ERR);
}

int GetSampleGpfNumber(int id_problem)
{
  switch (id_problem) {
  case 0: return(NUM_SAMPLE_AUGPF3D);
  case 1: return(NUM_SAMPLE_LMGPF3D);
  case 2: return(NUM_SAMPLE_MPGPF3D);
  case 3: return(NUM_SAMPLE_MXGPF3D);
  }
  return(NUM_SAMPLE_ERR);
}

int GetSampleIbrNumber(int id_problem)
{
  switch (id_problem) {
  case 0: return(NUM_SAMPLE_AUGPF3D+1);
  case 1: return(NUM_SAMPLE_LMGPF3D+1);
  case 2: return(NUM_SAMPLE_MPGPF3D+1);
  case 3: return(NUM_SAMPLE_MXGPF3D+1);
  }
  return(NUM_SAMPLE_ERR);
}

int GetSampleChainNumber()
{
  return(NUM_SAMPLE_CHAIN);
}

void * GetTable(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return((void *)cont->table);
  else      return(NULL);
}

char * GetSampleDescription(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->sample_description);
  else      return(NULL);
}

char * GetSampleComment(void * context)
{
  Context * cont = (Context *)context;
  if (cont) return(cont->sample_comment);
  else      return(NULL);
}

#ifndef ___ABRIDGE_PROFILE_MODE___
int GetSampleBlockNum(void * context)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm) return(cont->sm->N);
  else                  return(0);
}

//////////////////////////////////////////////////////////////////////////////
//...извлечение геометрии обpазца (или одного из его блоков) по его контексту;
CCells * get_bar(void * context, int id_block)
{
  Context * cont = (Context *)context;
  CCells  * bar  = NULL;
  if (cont && cont->sm) {
     if (id_block < 0 || id_block >= cont->sm->N) bar = cont->sm->bar;
     else                                         bar = cont->sm->B[id_block].bar;
  }
  return bar;
}

void * GetContextBar(void * context, int id_block)
{
  return (void *)get_bar(context, id_block);
}

//////////////////////////////////////////////////////////////////////////////////////
//...извлечение локальной системы координат одного из блоков образца по его контексту;
CMap * GetContextMap(void * context, int id_block)
{
  Context * cont = (Context *)context;
  CMap    * map  = NULL;
  if (cont && cont->sm) {
     if (0 <= id_block && id_block < cont->sm->N) map = cont->sm->B[id_block].mp;
  }
  return map;
}

/////////////////////////////////////////////////////
//...установка существующей геометрии внутрь образца;
void SetContextBar(void * context, void * bar, int id_block)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm && bar) {
      if (id_block < 0 || id_block >= cont->sm->N) {
          if (cont->sm->bar) {
              cont->sm->bar->zero_cells(); 
              delete cont->sm->bar; 
          }
          cont->sm->bar = (CCells *)bar;
      }
      else {
          if (cont->sm->B[id_block].bar) {
              cont->sm->B[id_block].bar->zero_cells(); 
              delete cont->sm->B[id_block].bar; 
          }
          cont->sm->B[id_block].bar = (CCells *)bar;
      }
  }
  return;
}
#endif

////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для любого образца по его абсолютному номеру;
void * CreateContext(int N_sm)
{
  if (N_sm >= SHIFT_LM3D_SAMPLES && N_sm < SHIFT_LM3D_SAMPLES+NUM_LM3D_SAMPLES)
#ifdef ___PRO_LAME3D_cpp___
      return CreateLM3DContext(N_sm-SHIFT_LM3D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_MX3D_SAMPLES && N_sm < SHIFT_MX3D_SAMPLES+NUM_MX3D_SAMPLES)
#ifdef ___PRO_MAXW3D_cpp___
      return CreateMX3DContext(N_sm-SHIFT_MX3D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_AU3D_SAMPLES && N_sm < SHIFT_AU3D_SAMPLES+NUM_AU3D_SAMPLES)
#ifdef ___PRO_ACOU3D_cpp___
      return CreateAU3DContext(N_sm-SHIFT_AU3D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_MP3D_SAMPLES && N_sm < SHIFT_MP3D_SAMPLES+NUM_MP3D_SAMPLES)
#ifdef ___PRO_MAPI3D_cpp___
      return CreateMP3DContext(N_sm-SHIFT_MP3D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_AU2D_SAMPLES && N_sm < SHIFT_AU2D_SAMPLES+NUM_AU2D_SAMPLES)
#ifdef ___PRO_ACOU2D_cpp___
      return CreateAU2DContext(N_sm-SHIFT_AU2D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_SK2D_SAMPLES && N_sm < SHIFT_SK2D_SAMPLES+NUM_SK2D_SAMPLES)
#ifdef ___PRO_SKIN2D_cpp___
      return CreateSK2DContext(N_sm-SHIFT_SK2D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_MP2D_SAMPLES && N_sm < SHIFT_MP2D_SAMPLES+NUM_MP2D_SAMPLES)
#ifdef ___PRO_MAPI2D_cpp___
      return CreateMP2DContext(N_sm-SHIFT_MP2D_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_BARS_SAMPLES && N_sm < SHIFT_BARS_SAMPLES+NUM_BARS_SAMPLES)
#ifdef ___PRO_BARS2D_cpp___
      return CreateBARSContext(N_sm-SHIFT_BARS_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_WIRE_SAMPLES && N_sm < SHIFT_WIRE_SAMPLES+NUM_WIRE_SAMPLES)
#ifdef ___PRO_WIRE2D_cpp___
      return CreateWIREContext(N_sm-SHIFT_WIRE_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_REDUCEDWIRE_SAMPLES && N_sm < SHIFT_REDUCEDWIRE_SAMPLES+NUM_REDUCEDWIRE_SAMPLES)
#ifdef ___PRO_WIRE2D_cpp___
      return CreateREDUCEDWIREContext(N_sm-SHIFT_REDUCEDWIRE_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_TOWER_SAMPLES && N_sm < SHIFT_TOWER_SAMPLES+NUM_TOWER_SAMPLES)
#ifdef ___PRO_TOWER3D_cpp___
      return CreateTOWERContext(N_sm-SHIFT_TOWER_SAMPLES);
#else
      return(NULL);
#endif
  if (N_sm >= SHIFT_FITTING_SAMPLES && N_sm < SHIFT_FITTING_SAMPLES+NUM_FITTING_SAMPLES)
#ifdef ___PRO_FITTING_cpp___
      return CreateFITTINGContext(N_sm-SHIFT_FITTING_SAMPLES);
#else
      return(NULL);
#endif
  return(NULL);
}
////////////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для любого образца по контексту (копироание таблицы);
void * CreateContext(void * context)
{
	Context * cont = (Context *)new_struct(sizeof(Context));
	if (! cont) return(NULL);

	cont->N           = GetSampleNumber(context);
	cont->static_char = INITIAL_STATE;
	cont->sample_name = GetSampleName(context);
	cont->units       = GetUnits(context);
	cont->table			= get_shablon_table((Table *)GetTable(context), cont->static_char);
	table_cpy(cont->table, (Table *)GetTable(context), cont->static_char);

	set_default   (cont);
	return((void *)cont);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция образования описания образца (рабочая функция)
void * CreateDescript(char * description)
{
  if (! description) return(NULL);
  int k = 0;

///////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление длины файла входных данных, образование контекста и копирование входных данных;
  while (description[k] != '\x0') k++;

  Context * cont  = (Context *)new_struct(sizeof(Context));
  char    * ascii = (char    *)new_struct((k+1)*sizeof(char));
  if (! cont || ! ascii) {
     delete_struct(ascii);
     delete_struct(ascii = (char *)cont); 
     return(NULL);
  }
  for (; k >= 0; k--) ascii[k] = description[k];

//////////////////////////
//...заполнение контекста;
  cont->N                  = NUM_SAMPLE_ERR;
  cont->units              = UNIT_SI;
  cont->sample_description = ascii;

  set_default(cont);
  return (void *)cont;
}

void SetSampleDescription(void * context, char * description)
{
  Context * cont = (Context *)context;
  if (! cont || ! description) return;
  int k = (int)strlen (description);

  delete_struct(cont->sample_description);
  if ((cont->sample_description = (char *)new_struct((k+1)*sizeof(char))) != NULL)
	  memcpy(cont->sample_description, description, k*sizeof(char));
}

void SetSampleDescription(void * context, const char * description)
{
  Context * cont = (Context *)context;
  if (! cont || ! description) return;
  int k = (int)strlen (description);

  delete_struct(cont->sample_description);
  if ((cont->sample_description = (char *)new_struct((k+1)*sizeof(char))) != NULL)
	  memcpy(cont->sample_description, description, k*sizeof(char));
}

void SetSampleComment(void * context, char * comment)
{
  Context * cont = (Context *)context;
  if (! cont || ! comment) return;
  int k = (int)strlen (comment);

  delete_struct(cont->sample_comment);
  if ((cont->sample_comment = (char *)new_struct((k+1)*sizeof(char))) != NULL)
	  memcpy(cont->sample_comment, comment, k*sizeof(char));
}

void SetSampleComment(void * context, const char * comment)
{
  Context * cont = (Context *)context;
  if (! cont || ! comment) return;
  int k = (int)strlen (comment);

  delete_struct(cont->sample_comment);
  if ((cont->sample_comment = (char *)new_struct((k+1)*sizeof(char))) != NULL)
	  memcpy(cont->sample_comment, comment, k*sizeof(char));
}

void SetTable(void * context, Table * table, int id_static_char)
{
  Context * cont = (Context *)context;
  if (! cont) return;

  if (cont->static_char == NULL_STATE && id_static_char != NULL_STATE) {
	  delete_struct(cont->sample_name);
	  delete_struct(cont->GOST_name);
  }
  if (id_static_char == NULL_STATE && cont->static_char != NULL_STATE) {
	  cont->sample_name = NULL;
	  cont->GOST_name = NULL;
  }
  delete_table(cont->table, cont->static_char);
  cont->static_char = id_static_char;
  cont->table = table;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для описания Gpf-образца
void * CreateGPContext(char * description, int id_problem)
{
  Context * cont        = (Context *)CreateDescript(description);
  if (cont) {
      switch (id_problem) {
      case 0: cont->N   = NUM_SAMPLE_LMGPF3D; break;
      case 1: cont->N   = NUM_SAMPLE_MXGPF3D; break;
      case 2: cont->N   = NUM_SAMPLE_AUGPF3D; break;
      case 3: cont->N   = NUM_SAMPLE_MPGPF3D; break;
      }
      cont->sample_name = "Gpf-sample";
  }
  return (void *)cont;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для описания Ibr-образца
void * CreateIBContext(char * description, int id_problem)
{
  Context * cont        = (Context *)CreateDescript(description);
  if (cont) {
      switch (id_problem) {
      case 0: cont->N   = NUM_SAMPLE_LMGPF3D+1; break;
      case 1: cont->N   = NUM_SAMPLE_MXGPF3D+1; break;
      case 2: cont->N   = NUM_SAMPLE_AUGPF3D+1; break;
      case 3: cont->N   = NUM_SAMPLE_MPGPF3D+1; break;
      }
      if (id_problem == 2) {
          cont->param = (Param *)new_struct(2*sizeof(Param));

      }
      cont->sample_name = "Inter_Bar Sample";
  }
  return (void *)cont;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция образования контекста для динамического образца (в теории кручения и изгиба стержней)
  void * CreateDSContext(char * description)
{
  Context * cont = (Context *)CreateDescript(description);
  if (cont) {
      cont->N    = NUM_SAMPLE_CHAIN;
  }
  return (void *)cont;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...предварительная инициализация для статического 3D образца
#ifndef ___ABRIDGE_PROFILE_MODE___
int sample3D_init(void * context, CGrid * block_nd, int type)
{
	Context * cont = (Context *)context;
	if (! cont) return(0);

#ifdef ___PRO_LAME3D_cpp___
	if (is_Lame3D(cont)) return lame3D_init(cont, block_nd); else
	if (is_LmGp3D(cont)) return lmgp3D_init(cont);           else
	if (is_LmIb3D(cont)) return lmib3D_init(cont);           else
#endif
#ifdef ___PRO_MAXW3D_cpp___
	if (is_Maxw3D(cont)) return maxw3D_init(cont, block_nd); else
	if (is_MxGp3D(cont)) return mxgp3D_init(cont);           else
	if (is_MxIb3D(cont)) return mxib3D_init(cont);           else
#endif
#ifdef ___PRO_ACOU3D_cpp___
	if (is_Acou3D(cont)) return acou3D_init(cont, block_nd, type); else
#endif
#ifdef ___PRO_ACOU2D_cpp___
	if (is_Acou2D(cont)) return acou2D_init(cont, block_nd); else
#endif
#if defined ___PRO_ACOU3D_cpp___ || defined ___PRO_ACOU2D_cpp___
	if (is_AuGp3D(cont)) return augp3D_init(cont);           else
	if (is_AuIb3D(cont)) return auib3D_init(cont, type);     else
#endif
#ifdef ___PRO_SKIN2D_cpp___
	if (is_Skin2D(cont)) return skin2D_init(cont, block_nd); else
#endif
#ifdef ___PRO_MAPI3D_cpp___
	if (is_Mapi3D(cont)) return mapi3D_init(cont, block_nd); else
#endif
#ifdef ___PRO_MAPI2D_cpp___
	if (is_Mapi2D(cont)) return mapi2D_init(cont, block_nd); else
#endif
#if defined ___PRO_MAPI3D_cpp___ || defined ___PRO_MAPI2D_cpp___
	if (is_MpGp3D(cont)) return mpgp3D_init(cont);           else
	if (is_MpIb3D(cont)) return mpib3D_init(cont);           else
#endif
#if defined ___PRO_BARS2D_cpp___
	if (is_BarsGOST2D(cont)) return bars2D_init(cont);       else
	if (is_Bars2D(cont)) return bars2D_init(cont);           else
#endif
#if defined ___PRO_WIRE2D_cpp___
	if (is_Wire2D(cont))  return wire2D_init(cont);          else
#endif
#if defined ___PRO_CABLE2D_cpp___
	if (is_Cable2D(cont))  return cable2D_init(cont);        else
#endif
#if defined ___PRO_TOWER3D_cpp___
	if (is_Tower3D(cont)) return tower3D_init(cont);         else
#endif
#if defined ___PRO_FITTING_cpp___
	if (is_Fitting(cont)) return fitting_init(cont);         else
#endif
	return(1);
}
#endif

///////////////////////////////////////
//...list of several problems creating;
void ** CreateNewList(int iSamplesCount)
{
  return(void **)new_struct((iSamplesCount+1)*sizeof(void *));
}

/////////////////////////////////////////
//...list of several problems destroying;
void DestroyList(void **& list)
{
  if (list) {
	  int k;
	  for ( k = 0;  NULL != list[k]; k++);
	  for ( k--;    k    >= 0;       k--) DeleteContext(list[k]);
      delete[](list); list = NULL;
  }
  return;
}

////////////////////////////////////////
//...creating the list of skin problems;
void ** CreateSKList()
{
#ifdef ___PRO_SKIN2D_cpp___
  int     k, N = GetSK2DSampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateSK2DContext(k);
       Table * table = (Table *)GetTable(list[k]);
  }
  return(list);
#else
  return(NULL);
#endif
}

////////////////////////////////////////
//...creating the list for GOST sections;
void ** CreateBARSGOSTList()
{
#ifdef ___PRO_BARS2D_cpp___
  int     k, N = GetBARSGOSTSampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateBARSGOSTContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

////////////////////////////////////////
//...creating the list of Bars25 module;
void ** CreateBARSList()
{
#ifdef ___PRO_BARS2D_cpp___
  int     k, N = GetBARSSampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateBARSContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

//////////////////////////////////////////////
//...creating the list of the wires data base;
void ** CreateWIREList()
{
#ifdef ___PRO_WIRE2D_cpp___
  int     k, N = GetWIRESampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateWIREContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

/////////////////////////////////////////////////////
//...creating the list of the reduced wire data base;
void ** CreateReducedWIREList()
{
#ifdef ___PRO_WIRE2D_cpp___
  int     k, N = GetREDUCEDWIRESampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateREDUCEDWIREContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

//////////////////////////////////////////////
//...creating the list of the cable data base;
void ** CreateCABLEList()
{
#ifdef ___PRO_CABLE2D_cpp___
  int     k, N = GetCABLESampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateCABLEContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

//////////////////////////////////////////////
//...creating the list of the tower data base;
void ** CreateTOWERList()
{
#ifdef ___PRO_TOWER3D_cpp___
  int     k, N = GetTOWERSampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateTOWERContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}

/////////////////////////////////////////////////
//...creating the list of the fittings data base;
void ** CreateFITTINGList()
{
#ifdef ___PRO_FITTING_cpp___
  int     k, N = GetFITTINGSampleCount();
  void ** list = CreateNewList(N);

  for (k = 0; list && k < N; k++) {
       list[k] = CreateFITTINGContext(k);
  }
  return(list);
#else
  return(NULL);
#endif
}
