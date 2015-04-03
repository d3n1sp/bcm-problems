/*===========================================*/
/*                  CBLOCKED                  */
/*===========================================*/
#ifndef ___CBLOCKED___
#define ___CBLOCKED___

#include "csolver.h"

////////////////////////////////////////////////////////////////
//     TEMPLATE BLOCKED SOLVER FOR BLOCK SYSTEM EQUATIONS     //
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//...class of solver for block system of equations on the base of GauusJ elimination;
template <typename T>
class CBlockedSolver : public CSolver<T> {
protected:
		CMatrix<T> pivot_C, pivot_S;
		void trf_step(T ** t0, T ** t1, T ** t2, int dim_N0, int dim_N1, int dim_N2, int id_zero = NULL_STATE);
		void hrf_step(T ** t0, T ** h1, T ** h2, int dim_N0, int dim_N1, int m, int id_zero = NULL_STATE);
		void hrf_step(T ** t0, T ** h0, int dim_N, int m1, int m2, int id_zero = NULL_STATE);
		int solver1(int m);
		int solver2(int m);
};
//
//////////////////////////////////////
////...transfer matrix transformation;
template <typename T>
void CBlockedSolver<T>::trf_step(T ** t0, T ** t1, T ** t2, int dim_N0, int dim_N1, int dim_N2, int id_zero)
{
	if (! t0 || ! t1 || ! t2) return;
	if (id_zero == OK_STATE) { //...прямое суммирование;
		for (int k = 0; k < dim_N1; k++) memset(t2[k], 0,  dim_N2*sizeof(T));
		for (int i = 0; i < dim_N1; i++)
		for (int l = 0; l < dim_N0; l++)
		for (int j = 0; j < dim_N2; j++) t2[i][j] += t0[l][i]*t1[l][j];
		return;
	}
	if (id_zero == ADDITIONAL_STATE) { //...инверсированное суммирование;
		for (int k = 0; k < dim_N1; k++) memset(t2[k], 0,  dim_N2*sizeof(T));
		for (int i = 0; i < dim_N1; i++)
		for (int l = 0; l < dim_N0; l++)
		for (int j = 0; j < dim_N2; j++) t2[i][j] += t0[i][l]*t1[l][j];
		return;
	}
	for (int i = 0; i < dim_N1; i++)
	for (int l = 0; l < dim_N0; l++)
	for (int j = 0; j < dim_N2; j++) t2[i][j] -= t0[l][i]*t1[l][j];
}

////////////////////////////////////
//...right-hand side transformation;
template <typename T>
void CBlockedSolver<T>::hrf_step(T ** t0, T ** h1, T ** h2, int dim_N0, int dim_N1, int m, int id_zero)
{
	if (! t0 || ! h1 || ! h2) return;
	if (id_zero == OK_STATE) { //...прямое суммирование;
		memset(h2[m], 0,  dim_N1*sizeof(T));
		for (int i = 0; i < dim_N1; i++)
		for (int l = 0; l < dim_N0; l++) h2[m][i] += t0[l][i]*h1[m][l];	  
		return;
  }
  if (id_zero == ADDITIONAL_STATE) { //...инверсированное суммирование;
		for (int i = 0; i < dim_N0; i++)
		for (int l = 0; l < dim_N1; l++) h2[m][i] -= t0[i][l]*h1[m][l];	  
		return;
	}
	for (int i = 0; i < dim_N1; i++)
	for (int l = 0; l < dim_N0; l++) h2[m][i] -= t0[l][i]*h1[m][l];
}

template <typename T>
void CBlockedSolver<T>::hrf_step(T ** t0, T ** h0, int dim_N, int m1, int m2, int id_zero)
{
	if (! t0 || ! h0) return;
	memset(h0[m2], 0,  dim_N*sizeof(T));
	if (id_zero == SPECIAL_STATE) { //...инверсированное суммирование;
		for (int i = 0; i < dim_N; i++)
		for (int l = 0; l < dim_N; l++) h0[m2][i] += t0[i][l]*h0[m1][l];	  
		return;
	}
	for (int i = 0; i < dim_N; i++)
	for (int l = 0; l < dim_N; l++) h0[m2][i] += t0[l][i]*h0[m1][l];	  
}

///////////////////////////////////////////////////////////////////////////////////
//...direct method of solving block matrix by blocked variant of Gauss elimination;
template <typename T>
int CBlockedSolver<T>::solver1(int k)
{
	int i, j, l, j_diag, k_elim, l_elim, N_dim, M_dim, L_dim, m;
///////////////////////////////////////////
//...устанавливаем вспомогательные матрицы;
	if (k == 0) {
		for (m = 0, i = 0; i < this->N; i++) 
			if (m < this->dim[i]) m = this->dim[i];
		pivot_C.set_matrix(m, m);
		if (! this->mode(RAPID_MODE)) 
		pivot_S.set_matrix(m, m);
		this->pivot_init(m);
////////////////////////////////////////////////////
//...если нет левой части, то копируем ее из правой;
		this->to_presenceTL();
	}

////////////////////////////////////////////////////
//...forward elimination in blocked Gauss-procedure;
	if (k < this->N && 0 <= k && this->p && 
		this->JR[this->p[k]] && this->TR[this->p[k]] && this->TL[this->p[k]]) {
		j_diag = this->JR[this->p[k]][this->JR_DIAG]-this->JR_SHIFT;
		N_dim = this->dim[this->p[k]];

////////////////////////////////
//...estimation of eigen values;
		if (this->mode(REDUCED_PRINT) | this->mode(FULLY_MODE)) {
			memcpy(pivot_C.GetMatrix()[0],   this->TR[this->p[k]][j_diag].GetMatrix()[0], N_dim*sizeof(T));
			//HQRfunction(pivot_C.GetMatrix(), this->hh[this->p[k]][0].GetMatrix()[n-1], this->hh[this->p[k]][0].GetMatrix()[n-2], N_dim);
			
			FILE * TST = fopen("eigen_all", k ? "a" : "w");
			test_gram_matrix(TST, this->p[k], END_STATE, EE);
			fclose(TST);
		}

/////////////////////////////////////
//...inversing of the diagonal block;
		if (! this->mode(RAPID_MODE))
		if (! this->JR[this->p[k]][this->JR_DIAG] || ! this->GaussJ(this->TR[this->p[k]][j_diag].GetMatrix(), N_dim)) {
			pivot_C.set_matrix(0, 0);
			pivot_S.set_matrix(0, 0);
			return(0);
		}

///////////////////////////////////////
//...vertical cycle on column elements;
		for (j = 0; j < this->JR[this->p[k]][0]; j++) 
			if (j != j_diag) {

/////////////////////////////
//...preparing pivot element;
				k_elim = this->JR[this->p[k]][j+this->JR_SHIFT];
				trf_step(this->TR[this->p[k]][j_diag].GetMatrix(), this->TL[this->p[k]][j].GetMatrix(), pivot_C.GetMatrix(), N_dim, N_dim, M_dim = this->dim[k_elim], OK_STATE);
				trf_step(this->TR[this->p[k]][j_diag].GetMatrix(), this->TR[this->p[k]][j].GetMatrix(), pivot_S.GetMatrix(), N_dim, N_dim, M_dim, ADDITIONAL_STATE);

//////////////////////////////
//...elimination block column;
				if (! this->mode(RAPID_MODE))
				for (l = 0; l < this->JR[this->p[k]][0]; l++) 
					if (this->p[this->N+(l_elim = this->JR[this->p[k]][l+this->JR_SHIFT])] >= this->p[this->N+k_elim])	{
						i = -1; while (++i < this->JR[k_elim][0] && this->JR[k_elim][i+this->JR_SHIFT] != l_elim);
///////////////////////////////////////////////////////
//...distribution of a new element and its elimination;
						if (i >= this->JR[k_elim][0])	{
							if (i+this->JR_SHIFT >= this->JR[k_elim][this->JR_BUFF]) {
								int * new_JR = (int *)new_struct((this->JR[k_elim][this->JR_BUFF] += this->N_band)*sizeof(int));
								memcpy(new_JR, this->JR[k_elim], (this->JR[k_elim][this->JR_BUFF]-this->N_band)*sizeof(int));
								delete_struct (this->JR[k_elim]); this->JR[k_elim] = new_JR;

								CMatrix<T> * new_TR = new CMatrix<T>[this->JR[k_elim][this->JR_BUFF]-this->JR_SHIFT];
								memcpy(new_TR, this->TR[k_elim], this->JR[k_elim][0]*sizeof(CMatrix<T>)); 
								memset(this->TR[k_elim], 0, this->JR[k_elim][0]*sizeof(CMatrix<T>));
								delete_struct(this->TR[k_elim]); this->TR[k_elim] = new_TR;

								CMatrix<T> * new_TL = new CMatrix<T>[this->JR[k_elim][this->JR_BUFF]-this->JR_SHIFT];
								memcpy(new_TL, this->TL[k_elim], this->JR[k_elim][0]*sizeof(CMatrix<T>)); 
								memset(this->TL[k_elim], 0, this->JR[k_elim][0]*sizeof(CMatrix<T>));
								delete_struct(this->TL[k_elim]); this->TL[k_elim] = new_TL;
							}
							this->TR[k_elim][i].set_matrix(M_dim, L_dim = this->dim[l_elim]);
							this->TL[k_elim][i].set_matrix(M_dim, L_dim);
							this->JR[k_elim][i+this->JR_SHIFT] = l_elim;
							this->JR[k_elim][0]++;
						}
						trf_step(pivot_C.GetMatrix(), this->TR[this->p[k]][l].GetMatrix(), this->TR[k_elim][i].GetMatrix(), N_dim, M_dim, L_dim = this->dim[l_elim]);

						if (this->p[this->N+l_elim] > this->p[this->N+k_elim])
						trf_step(pivot_S.GetMatrix(), this->TL[this->p[k]][l].GetMatrix(), this->TL[k_elim][i].GetMatrix(), N_dim, M_dim, L_dim);
				}
///////////////////////////////////
//...elimination of right-hand part;
				for (m = 0; m < this->id_norm; m++)
					hrf_step(pivot_C.GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), this->hh[k_elim][0].GetMatrix(), N_dim, M_dim, m);
		}
	}

////////////////////////////////////////
//...уничтожаем вспомогательные матрицы;
	if (k == -this->N) {
		pivot_C.set_matrix(0, 0);
		pivot_S.set_matrix(0, 0);
		delete[] this->pivot_ii; this->pivot_ii = NULL;
		delete[] this->pivot_kk; this->pivot_kk = NULL;
		delete[] this->pivot_ll; this->pivot_ll = NULL;
	}

//////////////////////////////////////////////////////////////////////////////////////
//...backward elimination in blocked Gauss-procedure (right-hand side transformation);
	if (-this->N <= k && k < 0 && this->p) {
		if (this->JR[this->p[k = -k-1]] && this->TR[this->p[k]] && this->TL[this->p[k]]) {

/////////////////////////////////////
//...elimination all previous element;
			for (N_dim = this->dim[this->p[k]], j = 0; j < this->JR[this->p[k]][0]; j++) 
				if (this->p[this->N+(k_elim = this->JR[this->p[k]][j+this->JR_SHIFT])] > k)
				for (m = 0; m < this->id_norm; m++)
					hrf_step(this->TR[this->p[k]][j].GetMatrix(), this->hh[k_elim][0].GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), N_dim, M_dim = this->dim[k_elim], m, ADDITIONAL_STATE);

//////////////////////////////////
//...elimination diagonal element;
			for (j_diag = this->JR[this->p[k]][this->JR_DIAG]-this->JR_SHIFT, m = 0; m < this->id_norm; m++) {
				hrf_step(this->TR[this->p[k]][j_diag].GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), N_dim, m, this->id_norm, SPECIAL_STATE);
				T * temp_h = this->hh[this->p[k]][0][m]; //...переставляем строчки матрицы правых частей; 
				this->hh[this->p[k]][0][m] = this->hh[this->p[k]][0][this->id_norm]; 
				this->hh[this->p[k]][0][this->id_norm] = temp_h;
			}

/////////////////////////////////////////////
//...resetting all new elements in block row;
			if (this->mode(CLEAN_MODE))
			for ( j = j_diag+1; j < this->JR[this->p[k]][0]; j++) {
				this->TR[this->p[k]][j].clean_matrix(N_dim, M_dim = this->dim[this->JR[this->p[k]][j+this->JR_SHIFT]]);
				this->TL[this->p[k]][j].clean_matrix(N_dim, M_dim);
			}
		}
	}
	return(1);
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...direct method of solving symmetric block matrix by blocked variant of Gauss elimination;
template <typename T>
int CBlockedSolver<T>::solver2(int k)
{
	int i, j, l, m, j_diag, k_elim, l_elim, N_dim, M_dim;
///////////////////////////////////////////
//...устанавливаем вспомогательные матрицы;
	if (k == 0) {
		for (m = i = 0; i < this->N; i++) 
			if (m < (int)this->dim[i]) m = this->dim[i];
		pivot_C.set_matrix(m, m);
		this->pivot_init(m);
	}

////////////////////////////////////////////////////
//...forward elimination in blocked Gauss-procedure;
	if (0 <= k && k < this->N && this->p && 
		this->JR[this->p[k]] && this->TR[this->p[k]]) {
		j_diag = this->JR[this->p[k]][this->JR_DIAG]-this->JR_SHIFT;
		N_dim = this->dim[this->p[k]];

////////////////////////////////
//...estimation of eigen values;
		if (this->mode(REDUCED_PRINT) | this->mode(FULLY_MODE)) {
			memcpy(pivot_C.GetMatrix()[0],   this->TR[this->p[k]][j_diag].GetMatrix()[0], N_dim*sizeof(T));
			//HQRfunction(pivot_C.GetMatrix(), this->hh[this->p[k]][0].GetMatrix()[n-1], this->hh[this->p[k]][0].GetMatrix()[n-2], N_dim);
			
			FILE * TST = fopen("eigen_all", k ? "a" : "w");
			test_gram_matrix(TST, this->p[k], END_STATE, EE);
			fclose(TST);
		}

/////////////////////////////////////
//...inversing of the diagonal block;
		if (! this->mode(RAPID_MODE))
		if (! this->JR[this->p[k]][this->JR_DIAG] || ! this->GaussJ(this->TR[this->p[k]][j_diag].GetMatrix(), N_dim)) {
			pivot_C.set_matrix(0, 0);
			return(0);
		}

///////////////////////////////////////
//...vertical cycle on column elements;
		for (j = 0; j < this->JR[this->p[k]][0]; j++) 
			if (j != j_diag) {

/////////////////////////////
//...preparing pivot element;
				k_elim = this->JR[this->p[k]][j+this->JR_SHIFT];
				trf_step(this->TR[this->p[k]][j_diag].GetMatrix(), this->TR[this->p[k]][j].GetMatrix(), pivot_C.GetMatrix(), N_dim, N_dim, M_dim = this->dim[k_elim], OK_STATE);

//////////////////////////////
//...elimination block column;
				if (! this->mode(RAPID_MODE))
				for (l = 0; l < this->JR[this->p[k]][0]; l++) 
					if (this->p[this->N+(l_elim = this->JR[this->p[k]][l+this->JR_SHIFT])] >= this->p[this->N+k_elim])	{
						i = -1; while (++i < this->JR[k_elim][0] && this->JR[k_elim][i+this->JR_SHIFT] != l_elim);
///////////////////////////////////////////////////////
//...distribution of a new element and its elimination;
						if (i >= this->JR[k_elim][0])	{
							if (i+this->JR_SHIFT >= this->JR[k_elim][this->JR_BUFF]) {
								int * new_JR = (int *)new_struct((this->JR[k_elim][this->JR_BUFF] += this->N_band)*sizeof(int));
								memcpy(new_JR, this->JR[k_elim], (this->JR[k_elim][this->JR_BUFF]-this->N_band)*sizeof(int));
								delete_struct (this->JR[k_elim]); this->JR[k_elim] = new_JR;

								CMatrix<T> * new_TR = new CMatrix<T>[this->JR[k_elim][this->JR_BUFF]-this->JR_SHIFT];
								memcpy(new_TR, this->TR[k_elim], this->JR[k_elim][0]*sizeof(CMatrix<T>)); 
								memset(this->TR[k_elim], 0, this->JR[k_elim][0]*sizeof(CMatrix<T>));
								delete_struct(this->TR[k_elim]); this->TR[k_elim] = new_TR;
							}
							this->TR[k_elim][i].set_matrix(M_dim, this->dim[l_elim]);
							this->JR[k_elim][i+this->JR_SHIFT] = l_elim;
							this->JR[k_elim][0]++;
						}
						trf_step(pivot_C.GetMatrix(), this->TR[this->p[k]][l].GetMatrix(), this->TR[k_elim][i].GetMatrix(), N_dim, M_dim, this->dim[l_elim]);
				}
///////////////////////////////////
//...elimination of right-hand part;
				for (m = 0; m < (int)this->id_norm; m++)
					hrf_step(pivot_C.GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), this->hh[k_elim][0].GetMatrix(), N_dim, M_dim, m);
		}
	}

////////////////////////////////////////
//...уничтожаем вспомогательные матрицы;
	if (k == -this->N) {
		pivot_C.set_matrix(0, 0);
		delete[] this->pivot_ii; this->pivot_ii = NULL;
		delete[] this->pivot_kk; this->pivot_kk = NULL;
		delete[] this->pivot_ll; this->pivot_ll = NULL;
	}

//////////////////////////////////////////////////////////////////////////////////////
//...backward elimination in blocked Gauss-procedure (right-hand side transformation);
	if (-this->N <= k && k < 0 && this->p) {
		if (this->JR[this->p[k = -k-1]] && this->TR[this->p[k]]) {

/////////////////////////////////////
//...elimination all previous element;
			for (N_dim = this->dim[this->p[k]], j = 0; j < this->JR[this->p[k]][0]; j++) 
				if (this->p[this->N+(k_elim = this->JR[this->p[k]][j+this->JR_SHIFT])] > (int)k)
				for (m = 0; m < (int)this->id_norm; m++)
					hrf_step(this->TR[this->p[k]][j].GetMatrix(), this->hh[k_elim][0].GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), N_dim, M_dim = this->dim[k_elim], m, ADDITIONAL_STATE);

//////////////////////////////////
//...elimination diagonal element;
			for (j_diag = this->JR[this->p[k]][this->JR_DIAG], m = 0; m < this->id_norm; m++) {
				hrf_step(this->TR[this->p[k]][j_diag].GetMatrix(), this->hh[this->p[k]][0].GetMatrix(), N_dim, m, this->id_norm, SPECIAL_STATE);
				T * temp_h = this->hh[this->p[k]][0][m]; //...переставляем строчки матрицы правых частей; 
				this->hh[this->p[k]][0][m] = this->hh[this->p[k]][0][this->id_norm]; 
				this->hh[this->p[k]][0][this->id_norm] = temp_h;
			}

////////////////////////////////////////////
//...resetting all new element in block row;
			if (this->mode(CLEAN_MODE))
			for ( j = j_diag+1; j < this->JR[this->p[k]][0]; j++)
				this->TR[this->p[k]][j].clean_matrix(N_dim, M_dim = this->dim[this->JR[this->p[k]][j+this->JR_SHIFT]]);
		}
	}
	return(1);
}
#endif
