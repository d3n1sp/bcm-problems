/*============================================*/
/*                  CCOMPUT3D                 */
/*============================================*/
#ifndef ___CCOMPUT3D___
#define ___CCOMPUT3D___

#include "ccomput2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

/////////////////////////////////////////////////////
//...,базовый класс для счетных схем блочного метода;
template <typename T>
class CComput3D : public CComput2D<T> {
public:
		int  block_comput(CGrid * bnd, int k, double square, int i, int & j_surf, int Max_N, int ANY = NULL_STATE, int id_fast = OK_STATE);
		int  block_comput(CGrid * bnd, int k, int m, double square, int & j_surf, int Max_N, int id_fast = OK_STATE);
		void block_comput(CGrid * bnd, int k, double square, int Max_N, int * descr = NULL,  int id_fast = OK_STATE);
		void block_comput(CGrid * bnd, int k, int m, double square, int Max_N, int * descr = NULL, int id_fast = OK_STATE);
protected:
		void comput1(int opt);//...счетная схема формирования блочной матрицы для конечно-элементной модели;
		void comput2(int opt);//...счетная схема формирования блочной матрицы для аналитической модели;
		void comput3(int opt);//...счетная схема (аналитическая) для периодической задачи в одном блоке;
		void comput4(int opt);//...счетная схема (аналитическая) для задачи Эшелби в одном блоке;
		void comput5(int opt, T * K, Num_Value _FMF, int id_variant);//...счетная схема для интегрирования эффективных характеристик по границе блока;
		void comput6(int opt, T * K, Num_Value _FMF, int id_variant);//...счетная схема для интегрирования эффективных характеристик по объему  блока;
};

/////////////////////////////////////////////////////////////
//          TEMPLATE VARIANT OF COUNTING SCHEME            //
/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//...counting of the boundary points for individual facet;
template <typename T>
int CComput3D<T>::block_comput(CGrid * bnd, int k, double square, int i, int & j_surf, int Max_N, int ANY, int id_fast)
{
	if (! bnd || k < 0 || this->N <= k || ! this->B[k].link || ! this->B[k].bar || ! this->B[k].bar->graph ||
					 i < 0 || this->B[k].bar->graph[0] <= i) return(0);
	double P[15], pp_bnd[4], px, py, pz, pm[12];
	int	i_surf, j, l, N_ini = 4, m = 0;
	if (Max_N <= 0) Max_N = 1;

//////////////////////////////////////////
//...ищем первую криволинейную поверхность;
	for (i_surf = ERR_STATE, l = 0; l < this->B[k].link[0]; l++)
	if  (this->B[k].link[l+1] <= SRF_STATE && l+1 != this->NUM_PHASE) {
		j_surf = this->B[k].link[(i_surf = l)+1];
		for (l = this->NUM_PHASE; l < this->B[k].link[0]; l++)
		if  (this->B[k].link[l+1] == -this->B[k].link[i_surf+1]+SRF_STATE && this->B[k].link[l+1] >= 0) j_surf = this->B[k].link[l+2];
		break;
	}
	if (0 <= i_surf && this->B[k].bar->graph[0] <= i_surf) return(0);

///////////////////////////////
//..generation boundary points;
	if (bnd->N_par < 4) bnd->add_params(4-bnd->N_par);
	IF_FACET(this->B[k].bar, i, l, ANY) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1],
			 N_min = min(N_arc, N_ini), m1, m2;
		for (j = 1; j <= N_min; j++, prev = arc) {
			arc = this->B[k].bar->ce[i]->graph[j+1];
			if (id_fast == OK_STATE) m1 = arc; else {
				m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
				m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
				if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			pm[3*(j-1)]   = this->B[k].bar->ce[i]->ce[m1]->mp[1];
			pm[3*(j-1)+1] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
			pm[3*(j-1)+2] = this->B[k].bar->ce[i]->ce[m1]->mp[3];
		}

///////////////////////////////////////////////////
//...снимаем данные о граничных условиях и нормали;
		pp_bnd[0] = this->B[k].bar->ce[i]->mp[l+2];
		pp_bnd[1] = this->B[k].bar->ce[i]->mp[l+3];
		pp_bnd[2] = this->B[k].bar->ce[i]->mp[l+4];
		pp_bnd[3] = this->B[k].bar->ce[i]->mp[l+5];
		if (square == 0. && N_arc == 4) {
			P[12] = this->B[k].bar->ce[i]->mp[4];
			P[13] = this->B[k].bar->ce[i]->mp[5];
			P[14] = this->B[k].bar->ce[i]->mp[6];
		}
		else {
			P[9]  = this->B[k].bar->ce[i]->mp[4];
			P[10] = this->B[k].bar->ce[i]->mp[5];
			P[11] = this->B[k].bar->ce[i]->mp[6];
		}

/////////////////////////////////////////
//...определяем положение кривой границы;
		if (0 <= i_surf) m = this->B[k].bar->ce[i]->common_id(this->B[k].bar->ce[i_surf]);
		l = m ? m : 2;

/////////////////////////
//...tria and quad facet;
		if (square == 0. && N_arc == 3) {
			for (j = 1; j <= N_arc; j++)
			memcpy(P+3*(j-1), pm+3*((j+l)%N_arc), 3*sizeof(double));
			bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
		}
		else
		if (square == 0. && N_arc == 4) {
			for (j = 1; j <= N_arc; j++)
			memcpy(P+3*(j-1), pm+3*((j+l-1)%N_arc), 3*sizeof(double));
			bnd->grid_quad_3D_refine(Max_N, Max_N, P, pp_bnd);
		}

/////////////////////////////////////
//...dividing any facet on triangles;
		else {
			memcpy(P, pm+3*((1+l)%3), 3*sizeof(double));
			for (j = 3; j <= N_min; j++) {
				memcpy(P+3, pm+3*((j+l-1)%3), 3*sizeof(double));
				memcpy(P+6, pm+3*((j+l)%3), 3*sizeof(double));

				px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
				py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
				pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
				if (square > 0.)
					bnd->grid_tria_3D_refine(max(Max_N, (int)(sqrt(max(pz, max(px, py))/square))), P, pp_bnd);
				else
					bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
			}
			memcpy(P, pm, 3*sizeof(double));
			memcpy(P+3, pm+3*(N_min-2), 6*sizeof(double));
			for (; j <= N_arc; j++, prev = arc) {
				arc = this->B[k].bar->ce[i]->graph[j+1];
				if (id_fast == OK_STATE) m1 = arc; else {
					m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
					m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
					if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
				}
				P[3] = P[6]; P[6] = this->B[k].bar->ce[i]->ce[m1]->mp[1];
				P[4] = P[7]; P[7] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
				P[5] = P[8]; P[8] = this->B[k].bar->ce[i]->ce[m1]->mp[3];

				px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
				py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
				pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
				if (square > 0.)
					bnd->grid_tria_3D_refine(max(Max_N, (int)(sqrt(max(pz, max(px, py))/square))), P, pp_bnd);
				else
					bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
			}
		}
	}
	return(m);
}

/////////////////////////////////////////////////////////////////////////////////
//...counting of the boundary points on lateral links sides (with curve surface);
template <typename T>
int CComput3D<T>::block_comput(CGrid * bnd, int k, int k_ext, double square, int & j_surf, int Max_N, int id_fast)
{
	if (! bnd || k < 0 || this->N <= k || ! this->B[k].link || ! this->B[k].bar || ! this->B[k].bar->graph || k_ext < 0 || this->N <= k_ext) return(0);
	double P[15], px, py, pz, pm[12];
	int	i_surf, i, j, l, N_ini = 4, m = 0;
	if (Max_N <= 0) Max_N = 1;

//////////////////////////////////////////
//...ищем первую криволинейную поверхность; 
	for (i_surf = ERR_STATE, l = 0; l < this->B[k].link[0]; l++)
	if  (this->B[k].link[l+1] <= SRF_STATE && l+1 != this->NUM_PHASE) { 
		j_surf = this->B[k].link[(i_surf = l)+1]; 
		for (l = this->NUM_PHASE; l < this->B[k].link[0]; l++)
		if  (this->B[k].link[l+1] == -this->B[k].link[i_surf+1]+SRF_STATE && this->B[k].link[l+1] >= 0) j_surf = this->B[k].link[l+2]; 
		break;
	}
	if (0 <= i_surf && this->B[k].bar->graph[0] <= i_surf) return(0);

///////////////////////////////
//..generation boundary points;.
	LOOP_LINK_FACET(this->B[k].bar, i, l, k_ext) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
		int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1], 
			 N_min = min(N_arc, N_ini), m1, m2;
		for (j = 1; j <= N_min; j++, prev = arc) {
			arc = this->B[k].bar->ce[i]->graph[j+1];
			if (id_fast == OK_STATE) m1 = arc; else {
				m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
				m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
				if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
			}
			pm[3*(j-1)]   = this->B[k].bar->ce[i]->ce[m1]->mp[1];
			pm[3*(j-1)+1] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
			pm[3*(j-1)+2] = this->B[k].bar->ce[i]->ce[m1]->mp[3];
		}

//////////////////////////////
//...снимаем данные о нормали;
		if (square == 0. && N_arc == 4) {
			P[12] = this->B[k].bar->ce[i]->mp[4];
			P[13] = this->B[k].bar->ce[i]->mp[5];
			P[14] = this->B[k].bar->ce[i]->mp[6];
		}
		else {
			P[9]  = this->B[k].bar->ce[i]->mp[4];
			P[10] = this->B[k].bar->ce[i]->mp[5];
			P[11] = this->B[k].bar->ce[i]->mp[6];
		}

/////////////////////////////////////////
//...определяем положение кривой границы;
		if (0 <= i_surf) m = this->B[k].bar->ce[i]->common_id(this->B[k].bar->ce[i_surf]);
		l = m ? m : 2;

/////////////////////////
//...tria and quad facet;
		if (square == 0. && N_arc == 3) {
			for (j = 1; j <= N_arc; j++) 
			memcpy(P+3*(j-1), pm+3*((j+l)%N_arc), 3*sizeof(double));
			bnd->grid_tria_3D_refine(Max_N, P);
		}
		else
		if (square == 0. && N_arc == 4) {
			for (j = 1; j <= N_arc; j++) 
			memcpy(P+3*(j-1), pm+3*((j+l-1)%N_arc), 3*sizeof(double));
			bnd->grid_quad_3D_refine(Max_N, Max_N, P);
		}

/////////////////////////////////
//...dividing any facet on triangles;
		else {
			memcpy(P, pm+3*((1+l)%3), 3*sizeof(double));
			for (j = 3; j <= N_min; j++) {
				memcpy(P+3, pm+3*((j+l-1)%3), 3*sizeof(double));
				memcpy(P+6, pm+3*((j+l)%3), 3*sizeof(double));

				px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
				py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
				pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
				if (square > 0.) 
					bnd->grid_tria_3D_refine(max(Max_N, (int)(sqrt(max(pz, max(px, py))/square))), P);
				else
					bnd->grid_tria_3D_refine(Max_N, P);

			}
			memcpy(P, pm, 3*sizeof(double));
			memcpy(P+3, pm+3*(N_min-2), 6*sizeof(double));
			for (; j <= N_arc; j++, prev = arc) {
				arc = this->B[k].bar->ce[i]->graph[j+1];
				if (id_fast == OK_STATE) m1 = arc; else {
					m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
					m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
					if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
				}
				P[3] = P[6]; P[6] = this->B[k].bar->ce[i]->ce[m1]->mp[1];
				P[4] = P[7]; P[7] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
				P[5] = P[8]; P[8] = this->B[k].bar->ce[i]->ce[m1]->mp[3];

				px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
				py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
				pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
				if (square > 0.)
					bnd->grid_tria_3D_refine(max(Max_N, (int)(sqrt(max(pz, max(px, py))/square))), P);
				else
					bnd->grid_tria_3D_refine(Max_N, P);
			}
		}
	}
	return(m);
}

/////////////////////////////////////////////////////
//...counting of the boundary points of block matrix;
template <typename T>
void CComput3D<T>::block_comput(CGrid * bnd, int k, double square, int Max_N, int * descr, int id_fast)
{
    if (! bnd || k < 0 || this->N <= k || ! this->B[k].link || ! this->B[k].bar) return;

    int    i, j, l, N_elem, N_ini = 4;
    double P[15], pp_bnd[4], px, py, pz, pm[12];
    if (Max_N <= 0) Max_N = 1;

///////////////////////////////
//..generation boundary points;
    if (bnd->N_par < 4) bnd->add_params(4-bnd->N_par); 
	 LOOP_FACET(this->B[k].bar, i, l) 
		 if (! descr || i < descr[0] && ! descr[i+1]) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
         int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1], 
             N_min = min(N_arc, N_ini), m1, m2;
         for (j = 1; j <= N_min; j++, prev = arc) {
				arc = this->B[k].bar->ce[i]->graph[j+1];
				if (id_fast == OK_STATE) m1 = arc; else {
					  m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
					  m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
					  if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
				}
				pm[3*(j-1)]   = this->B[k].bar->ce[i]->ce[m1]->mp[1];
				pm[3*(j-1)+1] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
				pm[3*(j-1)+2] = this->B[k].bar->ce[i]->ce[m1]->mp[3];
			}

///////////////////////////////////////////////////
//...снимаем данные о граничных условиях и нормали;
         pp_bnd[0] = this->B[k].bar->ce[i]->mp[l+2];
         pp_bnd[1] = this->B[k].bar->ce[i]->mp[l+3];
         pp_bnd[2] = this->B[k].bar->ce[i]->mp[l+4];
         pp_bnd[3] = this->B[k].bar->ce[i]->mp[l+5];
         if (square == 0. && N_arc == 4) {
             P[12] = this->B[k].bar->ce[i]->mp[4];
             P[13] = this->B[k].bar->ce[i]->mp[5];
             P[14] = this->B[k].bar->ce[i]->mp[6];
         }
         else {
             P[9]  = this->B[k].bar->ce[i]->mp[4];
             P[10] = this->B[k].bar->ce[i]->mp[5];
             P[11] = this->B[k].bar->ce[i]->mp[6];
         }

/////////////////////////
//...tria and quad facet;
         if (square == 0. && N_arc == 3) {
             memcpy(P, pm, 9*sizeof(double));
            bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
         }
         else
         if (square == 0. && N_arc == 4) {
             memcpy(P, pm, 12*sizeof(double));
            bnd->grid_quad_3D_refine(Max_N, Max_N, P, pp_bnd);
         }

/////////////////////////////////////
//...dividing any facet on triangles;
         else {
             memcpy(P, pm, 3*sizeof(double));
             for (j = 3; j <= N_min; j++) {
                 memcpy(P+3, pm+3*(j-2), 6*sizeof(double));
/////////////////////////////////////////////////////////////////////////
//...defining max diagonal and number of elements (when positive square);
                 if (square > 0.) {
                     px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                     py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                     pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                     if (px > py && px > pz) {
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                     }
                     if (py > px && py > pz) {
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                     }
                     N_elem = max(Max_N, (int)(sqrt(pz/square)));

                     bnd->grid_tria_3D_refine(N_elem, P, pp_bnd);
                 }
                 else {
                     px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                     py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                     pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                     if (px > py && px > pz) {
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                     }
                     if (py > px && py > pz) {
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                     }

                     bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
                 }
             }
             for (; j <= N_arc; j++, prev = arc) {
                  arc = this->B[k].bar->ce[i]->graph[j+1];
						if (id_fast == OK_STATE) m1 = arc; else {
							m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
							m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
							if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
						}
                  P[3] = P[6]; P[6] = this->B[k].bar->ce[i]->ce[m1]->mp[1];
                  P[4] = P[7]; P[7] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
                  P[5] = P[8]; P[8] = this->B[k].bar->ce[i]->ce[m1]->mp[3];

                  if (square > 0.) {
                      px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                      py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                      pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                      if (px > py && px > pz) {
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                      }
                      if (py > px && py > pz) {
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                      }
                      N_elem = max(Max_N, (int)(sqrt(pz/square)));

                      bnd->grid_tria_3D_refine(N_elem, P, pp_bnd);
                  }
                  else {
                      px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                      py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                      pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                      if (px > py && px > pz) {
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                      }
                      if (py > px && py > pz) {
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                      }
 
                      bnd->grid_tria_3D_refine(Max_N, P, pp_bnd);
                  }
             }
         }
    }

////////////////////////////////////
//...special case of infinity block;
    if (descr)
	 LOOP_FACET(this->B[k].bar, i, l) 
		if (i < descr[0] && descr[i+1] > 0) {

         pp_bnd[0] = this->B[k].bar->ce[i]->mp[l+2];
         pp_bnd[1] = this->B[k].bar->ce[i]->mp[l+3];
         pp_bnd[2] = this->B[k].bar->ce[i]->mp[l+4];
         pp_bnd[3] = this->B[k].bar->ce[i]->mp[l+5];

         int  num = -1; 
         for (num = 0; num < this->B[k].bar->graph[0]; num++) 
         if  (this->B[k].bar->ce[num]->cells_dim() == 2 && this->B[k].bar->ce[num]->mp && this->B[k].bar->ce[num]->graph && 
              this->B[k].bar->ce[num]->mp[l = size_of_map (this->B[k].bar->ce[num]->mp)] == FACET_CELL && 
              this->B[k].bar->ce[num]->mp[l+1] > 0. && 
              this->B[k].bar->ce[num]->mp[l+2] == MIN_HIT) break;

         double P[17]; P[15] = this->Acou3D_homotat; P[16] = this->Acou3D_L_inf;
         int descript[4]; 
         if (num < this->B[k].bar->graph[0] && this->B[k].bar->ce[num]->cells_dim() == 2 && num != i) {
             int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1], m1, m2;
             if (N_arc == 4) {
						for (l = 1; l <= N_arc; l++, prev = arc) {
							arc = this->B[k].bar->ce[i]->graph[l+1];
							if (id_fast == OK_STATE) m1 = arc; else {
								 m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
								 m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
								 if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
							}
							P[3*(l-1)]    = this->B[k].bar->ce[i]->ce[m1]->mp[1];
							P[3*(l-1)+1]  = this->B[k].bar->ce[i]->ce[m1]->mp[2];
							P[3*(l-1)+2]  = this->B[k].bar->ce[i]->ce[m1]->mp[3];
							if (! this->B[k].bar->ce[num]->topo_id(m1)) descript[l-1] = m1;
							else											  descript[l-1] = -m1-1;
						}
                  num = N_arc-1; while (num >= 0 && descript[num] >= 0) num--;
                  if (num > 0) {
                      double pp[12];
                      int    dd[4];
							 if  (descript[num] < 0) num--;
                      for (l = 0; l < N_arc; l++) {
                           dd[ (l+N_arc-num-1)%N_arc] = descript[l];
                           memcpy(pp+((l+N_arc-num-1)%N_arc)*3, P+l*3, 3*sizeof(double));
                      }
                      memcpy(descript, dd, 4*sizeof(int));
                      memcpy(P,    pp, 12*sizeof(double));
                  }
                  P[12] = this->B[k].bar->ce[i]->mp[4];
                  P[13] = this->B[k].bar->ce[i]->mp[5];
                  P[14] = this->B[k].bar->ce[i]->mp[6];

                  bnd->grid_quad_3D_refine_spec(this->Acou3D_N_elem2, Max_N, P, pp_bnd);
             }
         }
    }
}

////////////////////////////////////////////////////////////
//...counting of the boundary points on lateral links sides;
template <typename T>
void CComput3D<T>::block_comput(CGrid * bnd, int k, int m, double square, int Max_N, int * descr, int id_fast)
{
    if (! bnd || k < 0 || this->N <= k || ! this->B[k].link || ! this->B[k].bar || m < 0 || this->N <= m) return;

    int    i, j, l, N_elem, N_ini = 4;
    double P[15], px, py, pz, pm[12];
    if (Max_N <= 0) Max_N = 1;

///////////////////////////////
//..generation boundary points;.
	 LOOP_LINK_FACET(this->B[k].bar, i, l, m) 
		 if (! descr || i < descr[0] && ! descr[i+1]) {
         int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1], 
             N_min = min(N_arc, N_ini), m1, m2;
         for (j = 1; j <= N_min; j++, prev = arc) {
				arc = this->B[k].bar->ce[i]->graph[j+1];
				if (id_fast == OK_STATE) m1 = arc; else {
					m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
					m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
					if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
				}
				pm[3*(j-1)]   = this->B[k].bar->ce[i]->ce[m1]->mp[1];
				pm[3*(j-1)+1] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
				pm[3*(j-1)+2] = this->B[k].bar->ce[i]->ce[m1]->mp[3];
         }

///////////////////////////////////////////////////
//...снимаем данные о граничных условиях и нормали;
         if (square == 0. && N_arc == 4) {
             P[12] = this->B[k].bar->ce[i]->mp[4];
             P[13] = this->B[k].bar->ce[i]->mp[5];
             P[14] = this->B[k].bar->ce[i]->mp[6];
         }
         else {
             P[9]  = this->B[k].bar->ce[i]->mp[4];
             P[10] = this->B[k].bar->ce[i]->mp[5];
             P[11] = this->B[k].bar->ce[i]->mp[6];
         }

/////////////////////////
//...tria and quad facet;
         if (square == 0. && N_arc == 3) {
             memcpy(P, pm, 9*sizeof(double));
             bnd->grid_tria_3D_refine(Max_N, P);
         }
         else
         if (square == 0. && N_arc == 4) {
             memcpy(P, pm, 12*sizeof(double));
             bnd->grid_quad_3D_refine(Max_N, Max_N, P);
         }

/////////////////////////////////
//...dividing any facet on triangles;
         else {
             memcpy(P, pm, 3*sizeof(double));
             for (j = 3; j <= N_min; j++) {
                 memcpy(P+3, pm+3*(j-2), 6*sizeof(double));
                 if (square > 0.) {
                     px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                     py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                     pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                     if (px > py && px > pz) {
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                     }
                     if (py > px && py > pz) {
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                     }
                     N_elem = max(Max_N, (int)(sqrt(pz/square)));

                     bnd->grid_tria_3D_refine(N_elem, P);
                 }
                 else {
                     px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                     py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                     pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                     if (px > py && px > pz) {
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                     }
                     if (py > px && py > pz) {
                         swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                         swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                     }

                     bnd->grid_tria_3D_refine(Max_N, P);
                 }
             }
             for (; j <= N_arc; j++, prev = arc) {
                  arc = this->B[k].bar->ce[i]->graph[j+1];
						if (id_fast == OK_STATE) m1 = arc; else {
							m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
							m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
							if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
						}
                  P[3] = P[6]; P[6] = this->B[k].bar->ce[i]->ce[m1]->mp[1];
                  P[4] = P[7]; P[7] = this->B[k].bar->ce[i]->ce[m1]->mp[2];
                  P[5] = P[8]; P[8] = this->B[k].bar->ce[i]->ce[m1]->mp[3];

                  if (square > 0.) {
                      px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                      py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                      pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                      if (px > py && px > pz) {
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                      }
                      if (py > px && py > pz) {
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                      }
                      N_elem = max(Max_N, (int)(sqrt(pz/square)));

                      bnd->grid_tria_3D_refine(N_elem, P);
                  }
                  else {
                      px = sqr(P[0]-P[3])+sqr(P[1]-P[4])+sqr(P[2]-P[5]);
                      py = sqr(P[3]-P[6])+sqr(P[4]-P[7])+sqr(P[5]-P[8]);
                      pz = sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]);
                      if (px > py && px > pz) {
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                      }
                      if (py > px && py > pz) {
                          swap(P[3], P[6]); swap(P[4], P[7]); swap(P[5], P[8]); swap(py, pz);
                          swap(P[0], P[3]); swap(P[1], P[4]); swap(P[2], P[5]); swap(px, py);
                      }
 
                      bnd->grid_tria_3D_refine(Max_N, P);
                  }
             }
         }
    }

////////////////////////////////////
//...special case of infinity block;
    if (descr)
	 LOOP_LINK_FACET(this->B[k].bar, i, l, m) 
		if (i < descr[0] && descr[i+1] > 0) {
         int  num = -1; 
         for (num = 0; num < this->B[k].bar->graph[0]; num++) 
         if  (this->B[k].bar->ce[num]->cells_dim() == 2 && this->B[k].bar->ce[num]->mp && this->B[k].bar->ce[num]->graph && 
              this->B[k].bar->ce[num]->mp[l = size_of_map (this->B[k].bar->ce[num]->mp)] == FACET_CELL && 
              this->B[k].bar->ce[num]->mp[l+1] > 0. && 
              this->B[k].bar->ce[num]->mp[l+2] == MIN_HIT) break;

         double P[17]; P[15] = this->Acou3D_homotat; P[16] = this->Acou3D_L_inf;
         int descript[4]; 
         if (num < this->B[k].bar->graph[0] && this->B[k].bar->ce[num]->cells_dim() == 2 && num != i) {
             int N_arc, arc, prev = this->B[k].bar->ce[i]->graph[(N_arc = this->B[k].bar->ce[i]->graph[1])+1], m1, m2;
             if (N_arc == 4) {
					for (l = 1; l <= N_arc; l++, prev = arc) {
						arc = this->B[k].bar->ce[i]->graph[l+1];
						if (id_fast == OK_STATE) m1 = arc; else {
							 m1  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 0),
							 m2  = get_num(this->B[k].bar->ce[i]->ce[arc]->graph, 1);
							 if (! this->B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
						}
						P[3*(l-1)]    = this->B[k].bar->ce[i]->ce[m1]->mp[1];
						P[3*(l-1)+1]  = this->B[k].bar->ce[i]->ce[m1]->mp[2];
						P[3*(l-1)+2]  = this->B[k].bar->ce[i]->ce[m1]->mp[3];
						if (! this->B[k].bar->ce[num]->topo_id(m1)) descript[l-1] = m1;
						else											  descript[l-1] = -m1-1;
               }
               num = N_arc-1; while (num >= 0 && descript[num] >= 0) num--;
               if (num > 0) {
                   double pp[12];
                   int    dd[4];
						 if  (descript[num] < 0) num--;
                   for (l = 0; l < N_arc; l++) {
                        dd[ (l+N_arc-num-1)%N_arc] = descript[l];
                        memcpy(pp+((l+N_arc-num-1)%N_arc)*3, P+l*3, 3*sizeof(double));
                   }
                   memcpy(descript, dd, 4*sizeof(int));
                   memcpy(P,    pp, 12*sizeof(double));
               }
               P[12] = this->B[k].bar->ce[i]->mp[4];
               P[13] = this->B[k].bar->ce[i]->mp[5];
               P[14] = this->B[k].bar->ce[i]->mp[6];

               bnd->grid_quad_3D_refine_spec(this->Acou3D_N_elem2, Max_N, P);
				}
			}
		}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...основная счетная схема на основе анализа размерных элементов с учетом включения и периодического условия;
template <typename T>
void CComput3D<T>::comput1(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), id_dir, i, j, k, l, elem;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(6);

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(1);

	CGrid * block_phs = CreateNodes();
			  block_phs->add_params(1/*2*/);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int m_ilist = max(this->N/5, 1), nums_facet, m_cell, loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE);

	this->SetBlockBounding(par);

////////////////////////////////
//...discrete norm on inclusion;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete inclusion...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (l = this->NUM_PHASE; l < this->B[k].link[0]; l++) if (this->B[k].link[l+1] >= 0)
      for (i = 0; i < nums_facet; i++) 
		if (this->B[k].link[i+1] == -this->B[k].link[l+1]+SRF_STATE &&	this->B[k].bar->ce[i]->segms_id() && 
			(! this->solver.mode(NO_TR) || this->B[k].link[l+1] == this->solver.p[opt])) {
			this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);

			for (int lp = 0; lp < gauss_bnd->N; lp++) {
				pp[0] = gauss_bnd->get_param(0, lp);
				block_phs->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
												gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
			}
			gauss_bnd->add_buffer(gauss_bnd->N);

			if ((k % m_ilist) == 0) {
				sprintf(msg, "block %4i: block_phs->N = %i", k, block_phs->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING)
			this->TransferBB (block_phs, k, this->B[k].link[l+1], 0, E_PERIOD_COMPUT); else
			this->TransferBB (block_phs, k, this->B[k].link[l+1], 0,   PERIOD_COMPUT);

			if (! this->solver.mode(ACCUMULATION)) block_phs->add_buffer(block_phs->N);
		}
	}

////////////////////////////////////////////////
//...discrete norm on stitching sides of blocks;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete stitching sides...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++)
		if ((j = this->B[k].link[i+1]) >= 0 && this->B[k].bar->ce[i]->segms_id() &&
			 (! this->solver.mode(NO_TR) || j == this->solver.p[opt])) {
			this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);

			for (int lp = 0; lp < gauss_bnd->N; lp++) {
				pp[0] = gauss_bnd->get_param(0, lp);
				block_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
												gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
			}
			gauss_bnd->add_buffer(gauss_bnd->N);

			if ((k % m_ilist) == 0) {
				sprintf(msg, "block %4i: block_bnd->N = %i", k, block_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING)
			this->TransferBB (block_bnd, k, j, 0, E_BASIC_COMPUT); else
			this->TransferBB (block_bnd, k, j, 0,   BASIC_COMPUT);

			if (! this->solver.mode(ACCUMULATION)) block_bnd->add_buffer(block_bnd->N);
		}
   }

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++) 
		if ((j = this->B[k].link[i+1]) <= ERR_STATE && this->B[k].bar->ce[i]->segms_id()) {

/////////////////////////////////////////////////////////////////////
//...определяем есть ли связь с включением или периодической ячейкой;
			for (m_cell = BOUNDR_LINK, l = this->NUM_PHASE; l < this->B[k].link[0] && this->NUM_PHASE; l++)
			if ( this->B[k].link[l+1] == -this->B[k].link[i+1]+SRF_STATE && this->B[k].link[l+1] >= 0) m_cell = INCLUD_LINK;
			if (m_cell == BOUNDR_LINK && (elem = block_plink_3D(this->B[k], l = i, id_dir, par)) >= 0) m_cell = PERIOD_LINK;
		
/////////////////////////////////
//...накапливаем граничные точки;
			if  (m_cell != INCLUD_LINK && (! this->solver.mode(NO_TR) || elem == this->solver.p[opt])) {
				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max, pp);
				this->bnd_marker3D(j, pp);

				if (m_cell == PERIOD_LINK) { //...периодические граничные условия;
					pp[0] = par[1]-par[0];
					pp[1] = par[3]-par[2];
					pp[2] = par[5]-par[4];
					pp[3] = id_dir;
					pp[5] = elem;
				}
				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[4] = gauss_bnd->get_param(0, lp);
					bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);

				if ((k % m_ilist) == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				if (m_cell == PERIOD_LINK) { //...периодические граничные условия;
					if (this->solv >= ENERGY_SOLVING)
					this->GramAll(bound_bnd, k, 0, E_PERIOD_COMPUT); else
					this->GramAll(bound_bnd, k, 0,   PERIOD_COMPUT);
				}
				else { //...обычные граничные условия;
					if (this->solv >= ENERGY_SOLVING)
					this->GramAll(bound_bnd, k, 0, E_BASIC_COMPUT); else
					this->GramAll(bound_bnd, k, 0,   BASIC_COMPUT);
				}
				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
		}
	}
   delete bound_bnd;
   delete block_bnd;
   delete block_phs;
	delete gauss_bnd;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...счетная схема для конечно-элементной модели с учетом криволинейного включения и периодического условия;
template <typename T>
void CComput3D<T>::comput2(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), id_dir, i, j, k, l, m, elem;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(5);

	CGrid * bound_phs = CreateNodes();
			  bound_phs->add_params(6);

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(3);

	CGrid * block_phs = CreateNodes();
			  block_phs->add_params(2);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[6], P0[6], Po[12];
	int m_ilist = max(this->N/5, 1), j_surf, m_cell, nums_facet, loop;

	double par[6];	this->SetGeomBounding(par);

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]]);

///////////////////////////////////////////////////////
//...discrete norm on curvilinear or periodic boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++) if ((j = this->B[k].link[i+1]) <= ERR_STATE) {

/////////////////////////////////////////////////////////////////////
//...определяем есть ли связь с включением или периодической ячейкой;
			for (m_cell = BOUNDR_LINK, l = this->NUM_PHASE; l < this->B[k].link[0] && this->NUM_PHASE; l++)
			if  ( this->B[k].link[l+1] == -this->B[k].link[i+1]+SRF_STATE && this->B[k].link[l+1] >= 0) m_cell = INCLUD_LINK;
			if  (m_cell == BOUNDR_LINK && (elem = block_plink_3D(this->B[k], l = i, id_dir, par)) >= 0) m_cell = PERIOD_LINK;
		
			if  (m_cell != INCLUD_LINK && (! this->solver.mode(NO_TR) || elem == this->solver.p[opt])) {
				bnd->zero_grid();
				m = block_comput(bnd, k, sqr(this->get_param(this->NUM_QUAD+1)), i, j_surf, N_max);

/////////////////////////////////
//...накапливаем граничные точки;
				if (bnd->geom)
				for (l = 0; l <  bnd->geom[0]; l++) {
					int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
						num_f = num_n+num, cnt = 0;

					if (bnd->geom[num] == GL_TRIANGLES) {
						P0[3] = bnd->nX[bnd->geom[num+2]];
						P0[4] = bnd->nY[bnd->geom[num+2]];
						P0[5] = bnd->nZ[bnd->geom[num+2]];
						pp[0] = bnd->get_param(0, bnd->geom[num+2]);
						pp[1] = bnd->get_param(1, bnd->geom[num+2]);
						pp[2] = bnd->get_param(2, bnd->geom[num+2]);
						pp[3] = bnd->get_param(3, bnd->geom[num+2]);
						for (; num < num_f; num++) {
							Po[cnt++] = bnd->X[bnd->geom[num+2]];
							Po[cnt++] = bnd->Y[bnd->geom[num+2]];
							Po[cnt++] = bnd->Z[bnd->geom[num+2]];
						}
						gauss_bnd->facet_QG(Po, N_elem, NULL_STATE, j <= SRF_STATE ? OK_STATE : NULL_STATE);
						if (m_cell == PERIOD_LINK) { //...коррекция квадратур для периодической ячейки;
							pp[0] = par[1]-par[0];
							pp[1] = par[3]-par[2];
							pp[2] = par[5]-par[4];
							pp[3] = id_dir;
							pp[5] = elem;
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[4] = gauss_bnd->get_param(0, lp);
								bound_phs->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
							}
						}
						else {
							if (j <= SRF_STATE && this->bar && this->bar->graph && -j+SRF_STATE < this->bar->graph[0]) {
								gauss_bnd->QG_surface(this->bar->ce[-j+SRF_STATE]->mp); 
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[4] = gauss_bnd->get_param(0, lp);
									bound_bnd->add_new_point(gauss_bnd->X [lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																	 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
								}
							}
							else {
								if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_tria_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, Po);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[4] = gauss_bnd->get_param(0, lp);
									bound_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
								}
							}
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
					else 
					if (bnd->geom[num] == GL_QUAD_STRIP) {
						P0[3] = bnd->nX[bnd->geom[num+2]];
						P0[4] = bnd->nY[bnd->geom[num+2]];
						P0[5] = bnd->nZ[bnd->geom[num+2]];
						pp[0] = bnd->get_param(0, bnd->geom[num+2]);
						pp[1] = bnd->get_param(1, bnd->geom[num+2]);
						pp[2] = bnd->get_param(2, bnd->geom[num+2]);
						pp[3] = bnd->get_param(3, bnd->geom[num+2]);
						for (; num < num_f-3; num += 2) {
							Po[0] = bnd->X[bnd->geom[num+2]];
							Po[1] = bnd->Y[bnd->geom[num+2]];
							Po[2] = bnd->Z[bnd->geom[num+2]];
							Po[3] = bnd->X[bnd->geom[num+3]];
							Po[4] = bnd->Y[bnd->geom[num+3]];
							Po[5] = bnd->Z[bnd->geom[num+3]];
							Po[6] = bnd->X[bnd->geom[num+4]];
							Po[7] = bnd->Y[bnd->geom[num+4]];
							Po[8] = bnd->Z[bnd->geom[num+4]];
							Po[9] = bnd->X[bnd->geom[num+5]];
							Po[10] = bnd->Y[bnd->geom[num+5]];
							Po[11] = bnd->Z[bnd->geom[num+5]];

							gauss_bnd->facet_QG(Po, N_elem, OK_STATE, j <= SRF_STATE ? OK_STATE : NULL_STATE);
							if (m_cell == PERIOD_LINK) { //...коррекция квадратур для периодической ячеки;
								pp[0] = par[1]-par[0];
								pp[1] = par[3]-par[2];
								pp[2] = par[5]-par[4];
								pp[3] = id_dir;
								pp[5] = elem;
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[4] = gauss_bnd->get_param(0, lp);
									bound_phs->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
								}
							}
							else {
								if (j <= SRF_STATE && this->bar && this->bar->graph && -j+SRF_STATE < this->bar->graph[0]) {
									gauss_bnd->QG_surface(this->bar->ce[-j+SRF_STATE]->mp);
									for (int lp = 0; lp < gauss_bnd->N; lp++) {
										pp[4] = gauss_bnd->get_param(0, lp);
										bound_bnd->add_new_point(gauss_bnd->X [lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																		 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
									}
								}
								else {
									if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_quad_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, NULL, Po);
									for (int lp = 0; lp < gauss_bnd->N; lp++) {
										pp[4] = gauss_bnd->get_param(0, lp);
										bound_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
									}
								}
							}
							gauss_bnd->add_buffer(gauss_bnd->N);
						}
					}
				}
			}
			if (k%m_ilist == 0) {
				if (bound_phs->N) sprintf(msg, "block %4i: bound_phs->N = %i", k, bound_phs->N);
				else					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING) //...периодические условия скачка;
			this->GramAll(bound_phs, k, 0, E_PERIOD_COMPUT); else
			this->GramAll(bound_phs, k, 0,   PERIOD_COMPUT);
			if (! this->solver.mode(ACCUMULATION)) bound_phs->add_buffer(bound_phs->N);
			
			if (this->solv >= ENERGY_SOLVING) //...обычное граничное условие;
			this->GramAll(bound_bnd, k, 0, E_BASIC_COMPUT); else
			this->GramAll(bound_bnd, k, 0,   BASIC_COMPUT);
			if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
		}
	}

//////////////////////////////////////////////////////////////////
//...discrete norm on stitching sides of blocks (or on inclusion);
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete stitching sides or inclusion...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		for (i = 0; i < this->B[k].link[0]; i++) if ((j = this->B[k].link[i+1]) >= 0 && (! this->solver.mode(NO_TR) || j == this->solver.p[opt])) {
			bnd->zero_grid(); 
			m = block_comput(bnd, k, j, sqr(this->get_param(this->NUM_QUAD+1)), j_surf, N_max);

/////////////////////////////////
//...накапливаем граничные точки;
			if (bnd->geom)
			for (l = 0; l <  bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num, cnt = 0;

				if (bnd->geom[num] == GL_TRIANGLES) {
					P0[3] = bnd->nX[bnd->geom[num+2]];
					P0[4] = bnd->nY[bnd->geom[num+2]];
					P0[5] = bnd->nZ[bnd->geom[num+2]];
					for (; num < num_f; num++) {
						Po[cnt++] = bnd->X[bnd->geom[num+2]];
						Po[cnt++] = bnd->Y[bnd->geom[num+2]];
						Po[cnt++] = bnd->Z[bnd->geom[num+2]];
					}
					gauss_bnd->facet_QG(Po, N_elem, NULL_STATE, this->NUM_PHASE < i+1 ? OK_STATE : NULL_STATE);
					if (this->NUM_PHASE < i+1) { //...коррекция квадратур на границе фаз;
						if ((m = this->B[k].link[i+2]) <= SRF_STATE && this->bar && this->bar->graph && -m+SRF_STATE < this->bar->graph[0]) {
							gauss_bnd->QG_surface(this->bar->ce[-m+SRF_STATE]->mp);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[1] = j;
								pp[0] = gauss_bnd->get_param(0, lp);
								block_phs->add_new_point(gauss_bnd->X [lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
							}
						}
						else
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[1] = j;
							pp[0] = gauss_bnd->get_param(0, lp);
							block_phs->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp);
						}
					}
					else {
						if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_tria_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, Po);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[0] = gauss_bnd->get_param(0, lp);
							block_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
						}
					}
					gauss_bnd->add_buffer(gauss_bnd->N);
				}
				else 
				if (bnd->geom[num] == GL_QUAD_STRIP) {
					P0[3] = bnd->nX[bnd->geom[num+2]];
					P0[4] = bnd->nY[bnd->geom[num+2]];
					P0[5] = bnd->nZ[bnd->geom[num+2]];
					for (; num < num_f-3; num += 2) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = bnd->Z[bnd->geom[num+2]];
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = bnd->Z[bnd->geom[num+3]];
						Po[6] = bnd->X[bnd->geom[num+4]];
						Po[7] = bnd->Y[bnd->geom[num+4]];
						Po[8] = bnd->Z[bnd->geom[num+4]];
						Po[9] = bnd->X[bnd->geom[num+5]];
						Po[10] = bnd->Y[bnd->geom[num+5]];
						Po[11] = bnd->Z[bnd->geom[num+5]];

						gauss_bnd->facet_QG(Po, N_elem, OK_STATE, this->NUM_PHASE < i+1 ? OK_STATE : NULL_STATE);
						if (this->NUM_PHASE < i+1) { //...коррекция квадратур на границе фаз;
							if ((m = this->B[k].link[i+2]) <= SRF_STATE && this->bar && this->bar->graph && -m+SRF_STATE < this->bar->graph[0]) {
								gauss_bnd->QG_surface(this->bar->ce[-m+SRF_STATE]->mp);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[1] = j;
									pp[0] = gauss_bnd->get_param(0, lp);
									block_phs->add_new_point(gauss_bnd->X[ lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																	 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
								}
							}
							else
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[1] = j;
								pp[0] = gauss_bnd->get_param(0, lp);
								block_phs->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
							}
						}
						else {
							if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_quad_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, NULL, Po);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								block_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
							}
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
			}
			if (k%m_ilist == 0) {
				if (this->NUM_PHASE < i) sprintf(msg, "block %4i: block_phs->N = %i", k, block_phs->N);
				else					 sprintf(msg, "block %4i: block_bnd->N = %i", k, block_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING) //...граница фаз;;
			this->TransferBB (block_phs, k, j, 0, E_PERIOD_COMPUT); else
			this->TransferBB (block_phs, k, j, 0,   PERIOD_COMPUT);
			if (! this->solver.mode(ACCUMULATION)) block_phs->add_buffer(block_phs->N);

			if (this->solv >= ENERGY_SOLVING) //...обычная граница внутри материала;
			this->TransferBB (block_bnd, k, j, 0, E_BASIC_COMPUT); else
			this->TransferBB (block_bnd, k, j, 0,   BASIC_COMPUT);
			if (! this->solver.mode(ACCUMULATION)) block_bnd->add_buffer(block_bnd->N);
		}
	}
   delete bound_bnd;
	delete bound_phs;
   delete block_bnd;
   delete block_phs;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////////////////////////////
//...счетная схема для одного блока и периодического условия;
template <typename T>
void CComput3D<T>::comput3(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), id_dir, i, j, k, l;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(5);

	CGrid * block_phs = CreateNodes();
			  block_phs->add_params(1);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[5], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int m_ilist = max(this->N/5, 1), nums_facet, m_cell, loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE);

	this->SetBlockBounding(par);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	k = 0;
	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++) 
		if ((j = this->B[k].link[i+1]) <= ERR_STATE && this->B[k].bar->ce[i]->segms_id()) {

/////////////////////////////////////////////////////////////////////
//...определяем есть ли связь с включением или периодической ячейкой;
			for (m_cell = BOUNDR_LINK, l = this->NUM_PHASE; l < this->B[k].link[0] && this->NUM_PHASE; l++)
			if (this->B[k].link[l+1] == -this->B[k].link[i+1]+SRF_STATE && this->B[k].link[l+1] >= 0) m_cell = INCLUD_LINK;
			if (m_cell == BOUNDR_LINK && (id_dir = block_iddir_3D(this->B[k], i, par))) m_cell = PERIOD_LINK;

////////////////////////////////////////////////////////////////////
//...искусственная нумерация фасет для параллелепипеда с включением;
			if (this->BOX_LINK_PERIOD == OK_STATE && this->B[k].link[this->NUM_PHASE] == -1 && i > 0) { 
				m_cell = PERIOD_LINK;
				id_dir = i;
			}

/////////////////////////////////
//...накапливаем граничные точки;
			if (m_cell == PERIOD_LINK) { //...периодические граничные условия;
 				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);
				pp[0] = par[1]-par[0];
				pp[1] = par[3]-par[2];
				pp[2] = par[5]-par[4];
				pp[3] = id_dir;
				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[4] = gauss_bnd->get_param(0, lp);
					bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);

				if ((k % m_ilist) == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				this->gram2peri(bound_bnd, k, 0);

				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
			else
			if (m_cell == INCLUD_LINK) { //...связь с фазами материала;
 				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);

				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[0] = gauss_bnd->get_param(0, lp);
					block_phs->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);

				if ((k % m_ilist) == 0) {
					sprintf(msg, "block %4i: block_phs->N = %i", k, block_phs->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				this->TransferBB(block_phs, k, -this->B[k].link[i+1]+SRF_STATE, 0, PERIOD_COMPUT);

				if (! this->solver.mode(ACCUMULATION)) block_phs->add_buffer(block_phs->N);
			}
			else
			if (m_cell == BOUNDR_LINK) { //...связь с внутренней границей материала;
 				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max, pp);
				this->bnd_marker3D(j, pp);

				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[4] = gauss_bnd->get_param(0, lp);
					bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);

				if ((k % m_ilist) == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				this->GramAll(bound_bnd, k, 0, BASIC_COMPUT);

				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
		}
	}
   delete bound_bnd;
   delete block_phs;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////
//...счетная схема для задачи Эшелби;
template <typename T>
void CComput3D<T>::comput4(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), i, j, k, l;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * block_phs = CreateNodes();
			  block_phs->add_params(1);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[5], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int m_ilist = max(this->N/5, 1), nums_facet, m_cell, loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE);

	this->SetBlockBounding(par);

///////////////////////////////////////////
//...discrete norm on inclusion's boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete inclusion's boundary...");

	k = 0;
	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++) 
		if ((j = this->B[k].link[i+1]) <= ERR_STATE && this->B[k].bar->ce[i]->segms_id()) {

///////////////////////////////////////////
//...определяем есть ли связь с включением;
			for (m_cell = BOUNDR_LINK, l = this->NUM_PHASE; l < this->B[k].link[0] && this->NUM_PHASE; l++)
			if (this->B[k].link[l+1] == -this->B[k].link[i+1]+SRF_STATE && this->B[k].link[l+1] >= 0) m_cell = INCLUD_LINK;

/////////////////////////////////
//...накапливаем граничные точки;
			if (m_cell == INCLUD_LINK) { //...связь с фазами материала;
 				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);

				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[0] = gauss_bnd->get_param(0, lp);
					block_phs->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);

				if ((k % m_ilist) == 0) {
					sprintf(msg, "block %4i: block_phs->N = %i", k, block_phs->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				trans_esh(block_phs, k, -this->B[k].link[i+1]+SRF_STATE, 0);
				if (! this->solver.mode(ACCUMULATION)) block_phs->add_buffer(block_phs->N);
			}
		}
	}
   delete block_phs;
	delete gauss_bnd;
}

//////////////////////////////////////////////////////////////////////////////////
//...счетная схема для интегрирования эффективных характеристик по границе блоков;
template <typename T>
void CComput3D<T>::comput5(int opt, T * K, Num_Value _FMF, int id_variant)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), id_dir, i, j, k, l, m, elem;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(1);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[1], Pn[3], Po[12], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int m_ilist = max(this->N/5, 1), j_surf, m_cell, nums_facet, loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE);

	for (k = 0; k < this->N; k++) SkeletonBounding(this->B[k], par);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) if (this->B[k].bar && this->B[k].link) {
		LOOP_MASK  (opt, k);
		nums_facet = this->B[k].bar->arcs_number();

		for (i = 0; i < nums_facet; i++) {
			bnd->zero_grid();
			m = block_comput(bnd, k, sqr(this->get_param(this->NUM_QUAD+1)), i, j_surf, N_max, OK_STATE);

//////////////////////////////////////////////////////////////////////////
//...определяем внутреннюю границу, внешнюю границу, прямоугольную ячейку;
			for (m_cell = BOUNDR_LINK, l = this->NUM_PHASE; l < this->B[k].link[0] && this->NUM_PHASE; l++)
			if  ( this->B[k].link[l+1] == -this->B[k].link[i+1]+SRF_STATE && this->B[k].link[l+1] >= 0) m_cell = INCLUD_LINK;
			if  (m_cell == BOUNDR_LINK &&  this->B[k].link[i+1] >= 0)  m_cell = BLOCKS_LINK;
			if  (m_cell == BOUNDR_LINK && (elem = block_plink_3D(this->B[k], l = i, id_dir, par)) >= 0) m_cell = PERIOD_LINK;

/////////////////////////////////
//...накапливаем граничные точки;
			IF_ANY_FACET(this->B[k].bar, i) {
				if (bnd->geom)
				for (l = 0; l <  bnd->geom[0]; l++) {
					int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
						num_f = num_n+num, cnt = 0;

					if (bnd->geom[num] == GL_TRIANGLES) {
						Pn[0] = bnd->nX[bnd->geom[num+2]];
						Pn[1] = bnd->nY[bnd->geom[num+2]];
						Pn[2] = bnd->nZ[bnd->geom[num+2]];
						for (; num < num_f; num++) {
							Po[cnt++] = bnd->X[bnd->geom[num+2]];
							Po[cnt++] = bnd->Y[bnd->geom[num+2]];
							Po[cnt++] = bnd->Z[bnd->geom[num+2]];
						}
						gauss_bnd->facet_QG(Po, N_elem, NULL_STATE, (j = this->B[k].link[i+1]) <= SRF_STATE ? OK_STATE : NULL_STATE);
						
						if (m_cell != PERIOD_LINK && j <= SRF_STATE && this->bar && this->bar->graph && -j+SRF_STATE < this->bar->graph[0]) {
							gauss_bnd->QG_surface(this->bar->ce[-j+SRF_STATE]->mp); 
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								bound_bnd->add_new_point(gauss_bnd->X [lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
							}
						}
						else {
							if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_tria_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, Po);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								bound_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], Pn[0], Pn[1], Pn[2], pp); 
							}
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
					else 
					if (bnd->geom[num] == GL_QUAD_STRIP) {
						Pn[0] = bnd->nX[bnd->geom[num+2]];
						Pn[1] = bnd->nY[bnd->geom[num+2]];
						Pn[2] = bnd->nZ[bnd->geom[num+2]];
						for (; num < num_f-3; num += 2) {
							Po[0] = bnd->X[bnd->geom[num+2]];
							Po[1] = bnd->Y[bnd->geom[num+2]];
							Po[2] = bnd->Z[bnd->geom[num+2]];
							Po[3] = bnd->X[bnd->geom[num+3]];
							Po[4] = bnd->Y[bnd->geom[num+3]];
							Po[5] = bnd->Z[bnd->geom[num+3]];
							Po[6] = bnd->X[bnd->geom[num+4]];
							Po[7] = bnd->Y[bnd->geom[num+4]];
							Po[8] = bnd->Z[bnd->geom[num+4]];
							Po[9] = bnd->X[bnd->geom[num+5]];
							Po[10] = bnd->Y[bnd->geom[num+5]];
							Po[11] = bnd->Z[bnd->geom[num+5]];

							gauss_bnd->facet_QG(Po, N_elem, OK_STATE, (j = this->B[k].link[i+1]) <= SRF_STATE ? OK_STATE : NULL_STATE);
							
							
							if (m_cell != PERIOD_LINK && j <= SRF_STATE && this->bar && this->bar->graph && -j+SRF_STATE < this->bar->graph[0]) {
								gauss_bnd->QG_surface(this->bar->ce[-j+SRF_STATE]->mp);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[0] = gauss_bnd->get_param(0, lp);
									bound_bnd->add_new_point(gauss_bnd->X [lp], gauss_bnd->Y [lp], gauss_bnd->Z [lp], 
																	 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
								}
							}
							else {
								if (m && this->bar && this->bar->graph && -j_surf+SRF_STATE < this->bar->graph[0]) gauss_bnd->QG_quad_surface(this->bar->ce[-j_surf+SRF_STATE]->mp, NULL, Po);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[0] = gauss_bnd->get_param(0, lp);
									bound_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], Pn[0], Pn[1], Pn[2], pp); 
								}
							}
							gauss_bnd->add_buffer(gauss_bnd->N);
						}
					}
				}
			}
			else
///////////////////////////////////////////////////////////////////////
//...определяем связи с фазами материала и внутренней границей области;
			if ((m_cell == INCLUD_LINK || m_cell == BOUNDR_LINK) && this->B[k].bar->ce[i]->segms_id()) {
				this->B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem, N_max);

				for (int lp = 0; lp < gauss_bnd->N; lp++) {
					pp[0] = gauss_bnd->get_param(0, lp);
					bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
													gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp);
				}
				gauss_bnd->add_buffer(gauss_bnd->N);
			}

//////////////////////////////////
//...интегрирование границы блока;
 			if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
				sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (_FMF == ERR_VALUE) RigidyAll( bound_bnd, k, K, BASIC_COMPUT); else
			for (l = 0; l <  bound_bnd->N; l++)
			GetFuncAllValues(bound_bnd->X[l], bound_bnd->Y[l], bound_bnd->Z[l], K, k, _FMF, id_variant);

			if (! this->solver.mode(ACCUMULATION))  bound_bnd->add_buffer(bound_bnd->N);
		}
	}
   delete bound_bnd;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////////////////////////////////////////////////
//...счетная схема для интегрирования эффективных характеристик по объему блоков;
template <typename T>
void CComput3D<T>::comput6(int opt, T * K, Num_Value _FMF, int id_variant)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), i, j, k, l, m;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * gauss_bnd = CreateNodes(GRID_EL_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
	int m_ilist = max(this->N/5, 1), loop, N_ini = 4;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE);

/////////////////////////////////////
//...discrete norm on blocks volume;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete blocks volume...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {
////////////////////////////////////////////////
//...вставляем временную затычку для квадратуры;
		gauss_bnd->grid_box3D(N_elem, N_elem, N_elem); 
		l = 0;
		for (m = 0; m <= N_elem; m++) {
			double q_m = 1./N_elem; if (m == 0 || m == N_elem) q_m *= .5;
			for (j = 0; j <= N_elem; j++) {
				double q_j = 1./N_elem; if (j == 0 || j == N_elem) q_j *= .5;
				for (i = 0; i <= N_elem; i++) {
					double q_i = 1./N_elem; if (i == 0 || i == N_elem) q_i *= .5;
					//gauss_bnd->X[l]  -= .5;
					//gauss_bnd->Y[l]  -= .5;
					//gauss_bnd->Z[l]  -= .5;  
					gauss_bnd->hit[l] = -1; 
					gauss_bnd->set_param(0, l++, q_i*q_j*q_m);
				}
			}
		}

////////////////////////////////////////
//...корректируем для сферического слоя;
		double R1 = this->param[4], R2 = this->B[0].mp[7], fff, ddd, rrr;
		l = 0;
		for (m = 0; m <= N_elem; m++)
		for (j = 0; j <= N_elem; j++)
		for (i = 0; i <= N_elem; i++) {
			gauss_bnd->X[l] *= 2.*M_PI;
			gauss_bnd->Y[l] *=  M_PI;
			gauss_bnd->Z[l] *= R2-R1; gauss_bnd->Z[l] += R1;  

			fff = gauss_bnd->X[l]; ddd = gauss_bnd->Y[l]; rrr = gauss_bnd->Z[l];
			gauss_bnd->X[l] = rrr*cos(fff)*sin(ddd);
			gauss_bnd->Y[l] = rrr*sin(fff)*sin(ddd);
			gauss_bnd->Z[l] *= cos(ddd);

			gauss_bnd->set_param(0, l, gauss_bnd->get_param(0, l)*2.*sqr(M_PI)*(R2-R1)*rrr*sin(ddd)); l++;
		}

///////////////////////////////////////////
//...вносим данные в интегрируемую функцию;
		if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
			sprintf(msg, "block %4i: gauss_bnd->N = %i", k, gauss_bnd->N);
			if (! this->solver.mode(NO_MESSAGE)) Message(msg);
		}
		if (_FMF == ERR_VALUE) RigidyAll( gauss_bnd, k, K, VOLUME_COMPUT); else
		for (l = 0; l <  gauss_bnd->N; l++)
		GetFuncAllValues(gauss_bnd->X[l], gauss_bnd->Y[l], gauss_bnd->Z[l], K, k, _FMF, id_variant);

		if (! this->solver.mode(ACCUMULATION)) gauss_bnd->add_buffer(gauss_bnd->N);
	}
	delete gauss_bnd;
}
#undef  Message
#endif
