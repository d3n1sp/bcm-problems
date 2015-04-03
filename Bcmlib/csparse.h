/*===========================================*/
/*                  CSPARSE                  */
/*===========================================*/
#ifndef ___CSPARSE___
#define ___CSPARSE___

#include "utils.h"

///////////////////////////////////////
//...базовый класс разреженной матрицы;
template <typename T>
class CSparse {
public:
		int    N;  //...dimension of structure;
		int  * p;  //...tables of renumbering;
		int ** JR; //...right part index of sparse matrix;
		int ** IL; //...left part index of sparse matrix (optional);
		T ** TR, ** TL;		  //...right and left part of sparse matrix;
		T ** TT, ** TD, ** hh; //...diagonal part of sparse matrix as partional elements and right-hand part;
static int JR_BUFF, JR_DIAG, JR_SHIFT, hh_numbs, band_buf;
protected:
		T * new_T(int N_stru){ T * stru = new T[N_stru]; memset(stru, 0, N_stru*sizeof(T)); return(stru);}
		int  N_band, N_buf, buf_count; //...bufferization structure (band and elements);
		void release_struct(int k);
		void release_struct();
		int band_bufferizat(int k);
public:
		int  inverse_index();
		void reset_indexes();
		void release_IL	();
		void release_matrix (T **& matr);
		void release_matrix (int k);
		void restore_matrix (int k);
		int  struct_permutat(int id_action);
//...constructor and destructor;
		CSparse() {
			JR = IL = NULL; p = NULL;
			TR = TL = TT = TD = hh = NULL;
			N = buf_count = 0;
			N_band = band_buf; 
			N_buf = 50;
		};
		~CSparse(void) {
			release_struct();
		};
//...setting of sparse structure;
		void add_element ();
		void set_elements(int N_sm = 0);
		void add_link (int k, int id_link);
		void set_links(int k, int * links = NULL, int shift = 1);
};

//////////////////////////////////////////////////////
//        TEMPLATE VARIANT OF SPARSE MATRIX         //
//////////////////////////////////////////////////////
template <typename T> int CSparse<T>::JR_BUFF  = 1;
template <typename T> int CSparse<T>::JR_DIAG  = 2;
template <typename T> int CSparse<T>::JR_SHIFT = 3;
template <typename T> int CSparse<T>::hh_numbs = 4;
template <typename T> int CSparse<T>::band_buf = 20;

////////////////////////////////
//...deleting of sparse structure;
template <typename T>
void CSparse<T>::release_struct()
{
	for (int k = 0; k < N; k++) release_struct(k);
	delete_struct(p);
	delete_struct(JR);
	delete_struct(IL);

	delete_struct(TR);
	delete_struct(TL);
	delete_struct(TT);
	delete_struct(TD);
	delete_struct(hh);
	N = buf_count = 0;
	N_band = band_buf; 
}

//////////////////////////////////////////////
//...deleting of k-th row of sparse structure;
template <typename T>
void CSparse<T>::release_struct(int k)
{
	if (k < N && 0 <= k) {
		if (JR) delete_struct(JR[k]); 
		if (IL) delete_struct(IL[k]); 
		release_matrix(k);
	}
}

/////////////////////////////////////////////////////
//...deleting element of sparse structure as a whole;
template <typename T>
void CSparse<T>::release_matrix(T **& matr) 
{
	if (matr)
	for (int k = 0; k < N; k++)	
	delete_struct(matr[k]);
	delete_struct(matr);
}

//////////////////////////////////////////////////////
//...deleting elements in one row of sparse structure;
template <typename T>
void CSparse<T>::release_matrix(int k)
{
	if (k < N && 0 <= k) {
		if (TR) delete_struct(TR[k]);
		if (TL) delete_struct(TL[k]); 
		if (TT) delete_struct(TT[k]); 
		if (TD) delete_struct(TD[k]); 
		if (hh) delete_struct(hh[k]);
	}
}

/////////////////////////////////////////////////////
//...restore elements in one row of sparse structure;
template <typename T>
void CSparse<T>::restore_matrix(int k)
{
	release_matrix(k);
	if (k < N && 0 <= k && JR && JR[k]) {
		if (TR) TR[k] = new_T(JR[k][0]);
		if (TL) TL[k] = new_T(JR[k][0]);
		if (TT) TT[k] = new_T(JR[k][0]);
		if (TD) TD[k] = new_T(JR[k][0]);
		if (hh) hh[k] = new_T(hh_numbs);
	}
}

//////////////////////////////////////////////////////////////////
//...reset of indexes with the account of table of permutation;
template <typename T>
void CSparse<T>::reset_indexes()
{
	if (p && JR && TR && TL && TT && TD) {
		int ** new_JR = (int **)new_struct(N*sizeof(int *)), j, k, l;
		T ** new_TR = (T **)new_struct(N*sizeof(T *)),
		  ** new_TL = (T **)new_struct(N*sizeof(T *)),
		  ** new_TT = (T **)new_struct(N*sizeof(T *)),
		  ** new_TD = (T **)new_struct(N*sizeof(T *));

/////////////////////////////////////////////////////////
//...count sparse structure and allocation new arrays;
      for (k = 0; k < N; k++)
		for (j = 0; j < JR[k][0]; j++)
			if ( p[N+JR[k][j+JR_SHIFT]] < p[N+k]) new_JR[JR[k][j+JR_SHIFT]] = (int *)((long long)new_JR[JR[k][j+JR_SHIFT]]+1); 
			else new_JR[k] = (int *)((long long)new_JR[k]+1);
		for (k = 0; k < N; k++) {
			new_JR[k] = (int *)new_struct( ((l = (int)(long long)new_JR[k])+JR_SHIFT)*sizeof(int));
			new_JR[k][JR_BUFF] = l+JR_SHIFT;
			new_TR[k] = new_T(l);
			new_TL[k] = new_T(l);
			new_TT[k] = new_T(l);
			new_TD[k] = new_T(l);
		}

///////////////////////////////////////////////////////
//...packing new values of indexes and sparse elements;
     for (k = 0; k < N; k++)
		for (j = 0; j < JR[k][0]; j++)
			if (p[N+JR[k][j+JR_SHIFT]] < p[N+k]) {
				new_JR[JR[k][j+JR_SHIFT]][(l = new_JR[JR[k][j+JR_SHIFT]][0]++)+JR_SHIFT] = k; 
				new_TR[JR[k][j+JR_SHIFT]][l] = TR[k][j];
				new_TL[JR[k][j+JR_SHIFT]][l] = TL[k][j];
				new_TT[JR[k][j+JR_SHIFT]][l] = TT[k][j];
				new_TD[JR[k][j+JR_SHIFT]][l] = TD[k][j];
			} 
			else {
				new_JR[k][(l = new_JR[k][0]++)+JR_SHIFT] = JR[k][j+JR_SHIFT];
				if (JR[k][j+JR_SHIFT] == k) new_JR[k][JR_DIAG] = l+JR_SHIFT;
				new_TR[k][l] = TR[k][j]; 
				new_TL[k][l] = TL[k][j]; 
				new_TT[k][l] = TT[k][j]; 
				new_TD[k][l] = TD[k][j]; 
			}

///////////////////////
//...exchange pointers;
		for (N_band = band_buf, k = 0; k < N; k++) {
			memset(TR[k], 0, JR[k][0]*sizeof(T)); delete_struct(TR[k]); TR[k] = new_TR[k];
			memset(TL[k], 0, JR[k][0]*sizeof(T)); delete_struct(TL[k]); TL[k] = new_TL[k];
			memset(TT[k], 0, JR[k][0]*sizeof(T)); delete_struct(TT[k]); TT[k] = new_TT[k];
			memset(TD[k], 0, JR[k][0]*sizeof(T)); delete_struct(TD[k]); TD[k] = new_TD[k];
			delete_struct(JR[k]); JR[k] = new_JR[k]; N_band = max(N_band, JR[k][0]);
		}
		delete_struct(new_JR);
		delete_struct(new_TR);
		delete_struct(new_TL);
		delete_struct(new_TT);
		delete_struct(new_TD);
	}
}

////////////////////////////////////////////////////////////////
//...construction of indexes of block strucure (vertical incoming);
template <typename T>
int CSparse<T>::inverse_index()
{
	release_IL();
	if (! p || ! JR  || ! (IL = (int **)new_struct(N*sizeof(int *)))) return(0);
////////////////////////////////////////////////////////////////////
//...count structure of vertical index and allocation new arrays;
   for (int k = 0; k < N; k++)
	for (int j = 0; j < JR[k][0]; j++) 
		IL[JR[k][j+JR_SHIFT]] = (int *)((long long)IL[JR[k][j+JR_SHIFT]]+1);
	for (int k = 0; k < N; k++)
	if (! (IL[k] = (int *)new_struct(((int)(long long)IL[k]+JR_SHIFT)*sizeof(int)))) {
		release_IL(); 
		return(0);
	}
////////////////////////////
//...filling vertical index;
	for (int k = 0; k < N; k++)
	for (int j = 0; j < JR[p[k]][0]; j++) 
		IL[JR[p[k]][j+JR_SHIFT]][JR_SHIFT+IL[JR[p[k]][j+JR_SHIFT]][0]++] = p[k];
	return(1);
}

////////////////////////////////////////////////
//...deleting incoming index of sparse strucure;
template <typename T>
void CSparse<T>::release_IL() 
{
	if (IL)
	for (int k = 0; k < N; k++)
	delete_struct(IL[k]);
	delete_struct(IL);
}

////////////////////////////////////
//...permutation of block strucutre;
template <typename T>
int CSparse<T>::struct_permutat(int id_action)
{
  if (N > 0 && (p || (p = (int *)new_struct((2*N+2)*sizeof(int))) != NULL)) {
      int id_err  = 0, k, i, j;
		for (p[N*2] = 0, k = 0; k < N; k++)	//...initial band width;
		for (j = 0; j < JR[k][0]; j++) p[N*2] = max(p[N*2], max(JR[k][j+JR_SHIFT]-k, k-JR[k][j+JR_SHIFT]));

//...renumbering algorithm RCM;
		if (id_action != NULL_STATE) {
			int * xadj = (int *)new_struct((N+1)*sizeof(int)),
				 * xls  = (int *)new_struct((N+1)*sizeof(int)), * iadj, m; 
			for (id_err = 1, m = k = 0; k < N; k++) m += JR[k][0]*2-1;

			if  (xadj && xls && (iadj = (int *)new_struct(m*sizeof(int))) != NULL) { //...reflect sample structure onto (xadj,iadj) - arrays;
			 	for (xadj[0] = (m = 0)+1, k = 0; k < N; k++) { //...fill arrays, mapped rows;
					for (i = max(0, k-p[N*2]); i < k; i++) {
						for (j = JR[i][0]; j > 0; j--) 
						if (JR[i][j+JR_SHIFT-1] == k) break;
						if (j) iadj[m++] = i+1;
					}
					for (j = 0; j < JR[k][0]; j++) iadj[m++] = JR[k][j+JR_SHIFT]+1;
					p[N+k] = (xadj[k+1] = m+1)-xadj[k];
				}
				 genrcm_utils(N, xadj, iadj, p, xls, (N+p)); id_err = 0;//...reverse Cuthill-Mckee ordering;
				 for (k = 0; k < N; k++) p[--p[k]+N] = k;

				 reset_indexes();
				 for (p[N*2+1] = 0, k = 0; k < N; k++)	//...renumbered band width;
				 for (j = 0; j < JR[k][0]; j++) p[N*2+1] = max(p[N*2+1], max(p[N+JR[k][j+JR_SHIFT]]-p[N+k], p[N+k]-p[N+JR[k][j+JR_SHIFT]]));
			}
			else
			for (k = 0; k < N; k++) p[k] = p[N+k] = k;
			delete_struct(xadj);
			delete_struct(iadj);
			delete_struct(xls);
		}
		else
      for (k = 0; k < N; k++) p[k] = p[N+k] = k;
      return(! id_err);
  }
  return(0);
}

////////////////////////////////////////
//...adding element to sparse structure;
template <typename T>
void CSparse<T>::add_element()
{
	if (buf_count == 0) {
		int ** new_JR = (int **)new_struct((N+N_buf)*sizeof(int *));
		if (JR) {
			memcpy(new_JR, JR, N*sizeof(int *)); delete_struct(JR);
		}
		JR  = new_JR;
		T ** new_TR = (T **)new_struct((N+N_buf)*sizeof(T *));
		if (TR) {
			memcpy (new_TR, TR, N*sizeof(T *)); delete_struct(TR);
		} 
		TR  = new_TR;
		T ** new_TL = (T **)new_struct((N+N_buf)*sizeof(T *));
		if (TL) {
			memcpy (new_TL, TL, N*sizeof(T *)); delete_struct(TL);
		}
		TL  = new_TL;
		T ** new_TT = (T **)new_struct((N+N_buf)*sizeof(T *));
		if (TT) {
			memcpy (new_TT, TT, N*sizeof(T *)); delete_struct(TT);
		}
		TT  = new_TT;
		T ** new_TD = (T **)new_struct((N+N_buf)*sizeof(T *));
		if (TD) {
			memcpy (new_TD, TD, N*sizeof(T *)); delete_struct(TD);
		}
		TD  = new_TD; 
		T ** new_hh = (T **)new_struct((N+N_buf)*sizeof(T *));
		if (hh) {
			memcpy(new_hh, hh, N*sizeof(T *)); delete_struct(hh);
		} 
		hh = new_hh; buf_count = N_buf;
	}
	N++; buf_count--;
	int j = 0; set_links(N-1, &j); //...initial settings;
}

//////////////////////////////////
//...distribution set of elements;
template <typename T>
void CSparse<T>::set_elements(int N_sm)
{
	release_struct();
	if ((N = N_sm) > 0) {
		JR = (int **)new_struct(N*sizeof(int *));
		TR = (T **)new_struct(N*sizeof(T *));
		TL = (T **)new_struct(N*sizeof(T *));
		TT = (T **)new_struct(N*sizeof(T *));
		TD = (T **)new_struct(N*sizeof(T *));
		hh = (T **)new_struct(N*sizeof(T *));
//////////////////////
//...initial settings;
		p = (int *)new_struct((2*N+2)*sizeof(int));
		for (int k = 0, j = 0; k < N; k++) {
			set_links(k, &j); p[k] = p[N+k] = k; 
		}
	}
}

/////////////////////////////////////////
//...bufferization row in band structure;
template <typename T>
int  CSparse<T>::band_bufferizat(int k)
{
	int m = 1;
	if (k < N && 0 <= k && JR) {
		int * new_JR_k = (int *)new_struct((N_band+JR[k][JR_BUFF])*sizeof(int));
		if (      JR[k]) memcpy(new_JR_k, JR[k], JR[k][JR_BUFF]*sizeof(int));
		if (! new_JR_k) m = 0; else {
			delete_struct(JR[k]); JR[k] = new_JR_k; JR[k][JR_BUFF] += N_band;
		}
		if (m && TR) {
			T * new_TR_k = new_T(N_band+JR[k][0]);
			if (		 TR[k]) memcpy(new_TR_k, TR[k], JR[k][0]*sizeof(T)); 
			if (! new_TR_k) m = 0; else {
				memset(TR[k], 0, JR[k][0]*sizeof(T));
				delete_struct(TR[k]); TR[k] = new_TR_k;
			}
		}
		if (m && TL) {
			T * new_TL_k = new_T(N_band+JR[k][0]);
			if (		 TL[k]) memcpy(new_TL_k, TL[k], JR[k][0]*sizeof(T)); 
			if (! new_TL_k) m = 0; else {
				memset(TL[k], 0, JR[k][0]*sizeof(T));
				delete_struct(TL[k]); TL[k] = new_TL_k;
			}
		}
		if (m && TT) {
			T * new_TT_k = new_T(N_band+JR[k][0]);
			if (		 TT[k]) memcpy(new_TT_k, TT[k], JR[k][0]*sizeof(T)); 
			if (! new_TT_k) m = 0; else {
				memset(TT[k], 0, JR[k][0]*sizeof(T));
				delete_struct(TT[k]); TT[k] = new_TT_k;
			}
		}
		if (m && TD) {
			T * new_TD_k = new_T(N_band+JR[k][0]);
			if (		 TD[k]) memcpy(new_TD_k, TD[k], JR[k][0]*sizeof(T)); 
			if (! new_TD_k) m = 0; else {
				memset(TD[k], 0, JR[k][0]*sizeof(T));
				delete_struct(TD[k]); TD[k] = new_TD_k;
			}
		}
	}
	return(m);
}

////////////////////////////////////////////////
//...setting one link into the sparse structure;
template <typename T>
void CSparse<T>::add_link(int k, int id_link)
{
	if (JR && k < N && 0 <= k) {
		int j;
		if (id_link < k && 0 <= id_link && (JR[id_link][0]+JR_SHIFT < JR[id_link][JR_BUFF] || band_bufferizat(id_link))) {
			for (j = JR[id_link][0]; j > 0; j--) 
			if (JR[id_link][j+JR_SHIFT-1] == k) break;
			if (! j) {
				JR[id_link][JR_SHIFT+JR[id_link][0]++] = k; 
				if (JR[id_link][JR_DIAG]) { //...меняем местами последний элемент и диагональ;
					int temp = JR[id_link][JR_SHIFT+JR[id_link][0]-1]; 
					JR[id_link][JR_SHIFT+JR[id_link][0]-1] = JR[id_link][JR[id_link][JR_DIAG]]; 
					JR[id_link][JR[id_link][JR_DIAG]] = temp;
					JR[id_link][JR_DIAG] = JR[id_link][0]+JR_SHIFT-1;
				}
			}
		}
		if (k <= id_link && 0 <= id_link && (JR[k][0]+JR_SHIFT < JR[k][JR_BUFF] || band_bufferizat(k))) {
			for (j = JR[k][0]; j > 0; j--) 
			if (JR[k][j+JR_SHIFT-1] == id_link) break;
			if (! j) {
				JR[k][JR_SHIFT+JR[k][0]++] = id_link;
				if (JR[k][JR_DIAG]) { //...меняем местами последний элемент и диагональ;
					int temp = JR[k][JR_SHIFT+JR[k][0]-1]; 
					JR[k][JR_SHIFT+JR[k][0]-1] = JR[k][JR[k][JR_DIAG]]; 
					JR[k][JR[k][JR_DIAG]] = temp;
					JR[k][JR_DIAG] = JR[k][0]+JR_SHIFT-1;
				}
			}
		}
	}
}

/////////////////////////////////////////////
//...setting links into the sparse structure;
template <typename T>
void CSparse<T>::set_links (int k, Topo * links, int shift)
{
	release_struct(k);
	if (k < N && 0 <= k && links && shift > 0 && p && JR) {
		int  kR, i;
      for (kR = 1, i = 0; i < links[0]; i++)
      if  (0 <= links[i+shift] && p[N+links[i+shift]] > p[N+k]) kR++;
/////////////////////////////////////////////////////////
//...распределяем и заполняем индексный массив элементов;
		if (( JR[k] = (int *)new_struct((kR+JR_SHIFT)*sizeof(int))) != NULL) {
			JR[k][0] = kR;	  N_band = max(N_band, JR[k][0]);
			JR[k][JR_BUFF] = kR+JR_SHIFT;
			JR[k][JR_DIAG] = kR+JR_SHIFT-1;
			JR[k][JR[k][JR_DIAG]] = k;
			for (i = kR = 0; i < links[0]; i++)
			if (0 <= links[i+shift] && p[N+links[i+shift]] > p[N+k])	JR[k][JR_SHIFT+kR++] = links[i+shift];
		}
	}
	restore_matrix(k);
}
#endif
