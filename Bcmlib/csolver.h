/*===========================================*/
/*                  CSOLVER                  */
/*===========================================*/
#ifndef ___CSOLVER___
#define ___CSOLVER___

#include "csparse.h"

///////////////////////////////////////////////////////
//...описание существующих типов связей между ячейками;
enum Num_Links {
			BOUNDR_LINK = 0,
			INCLUD_LINK,
			PERIOD_LINK,
			BLOCKS_LINK,
			  NUMS_LINK
};

//////////////////////////////////////////////////////////
//...набор двоичных масок для параллельной версии солвера;
enum Num_Parallel_Masks {
		NULL_MODE		 = 0x00000000L, //...нулевая мода - ничего не делаем;
		NO_TR				 = 0x00000001L, //...не заполняем массивы TR, TT и диагональ;
		NO_TL				 = 0x00000002L, //...не заполняем массивы TL, TD;
		NO_PHASE			 = 0x00000004L, //...отключение коррекции связей на границе фаз;
		ACCUMULATION	 = 0x00000008L, //...включение режима накапливания узлов квадратуры на звеньях;
		PROCESSOR_ID	 = 0x00000010L, //...включение режима записи номера процессора в h[1];
		REGULARIZATION  = 0x00000020L, //...включение регуляризации матрицы;
		REGUL_BOUNDARY  = 0x00000040L, //...включение регуляризации матрицы через граничное условие;
		REDUCED_MESSAGE = 0x00000080L, //...отключение подробной печати сообщений;
///////////////////////////////////////////////////////////////////////////////////
//...набор двоичных масок для организации управления солвером на внутреннем уровне;
		RAPID_MODE	  = 0x00000100L, //...включение быстрого режима (формируется только правая часть, матрица не обращается);
		CLEAN_MODE	  = 0x00000200L, //...включение режима обнуления дополнительно порожденных матриц;
		CONTI_MODE	  = 0x00000400L, //...включение режима продолжения;
		TESTI_MODE	  = 0x00000800L, //...включение режима тестирования;
		PRINT_MODE	  = 0x00001000L, //...включение режима тестовой печати промежуточной информации;
		FULLY_MODE	  = 0x00002000L, //...включение режима расширенного вывода промежуточной информации;
		MASKS_MODE	  = 0x00004000L, //...включение режима вывода маски подматриц-блоков;
		NO_MESSAGE	  = 0x00008000L, //...отключение печати сообщений;
		REDUCED_PRINT = 0x00010000L, //...включение режима сокращенного вывода промежуточной информации;
//////////////////////////////////////////////////////////////////
//...тестовое формирование матрицы Грама и энергетической матрицы;
		TESTI_GRAM	  = 0x00020000L, //...включение режима формирования только блочной матрицы Грама;
		TESTI_ENERGY  = 0x00040000L, //...включение режима формирования только энергетической матрицы;
};

/////////////////////////////////////////
//...базовый класс прямоугольной матрицы;
template <typename T>
class CMatrix {
protected:
		T ** m; //...объединенный массив элементов и указателей строчек (матрица);
public:
      inline T**& GetMatrix(){ return m;}
//...constructor and destructor;
		CMatrix()	 { m = NULL;}
	  ~CMatrix(void){ delete[] m;}
//...setting elements;
		void   set_matrix(int dim_N, int dim_M, int dim_buffer = 0) { ::set_matrix(m, dim_N, dim_M, dim_buffer);}
		void clean_matrix(int dim_N, int dim_M) { ::clean_matrix(m, dim_N, dim_M);}
		void clean_row_matrix(int k, int dim_M) {
			if (m) memset(m[k], 0, dim_M*sizeof(T));
		}
//...loaded operator [];
		inline T *& operator[](int _i) { 
			return m[_i];
		}
};

//////////////////////////////////////////////////////////////
//...базовый класс солвера блочной системы линейных уравнений;
template <typename T>
class CSolver : public CSparse<CMatrix<T> > {
public:
		int n, id_norm, id_change; //...number of potentials and type of solver;
		int * dim;					   //...blocks dimensions;
protected:
		unsigned long long mask; //...state mask of the solver;
		int * pivot_ii, * pivot_kk, * pivot_ll; //...for GaussJ;
		virtual int solver1(int m){ return(OK_STATE);}
		virtual int solver2(int m){ return(OK_STATE);}
public:
//...solver controls;
		void change_state(int set_change = OK_STATE) { id_change = set_change;}
		int  changed(int change) {return(change == id_change);}
		void set_mode	 (unsigned long long MODE = NULL_MODE) { mask |=  MODE;}
		void clean_mode (unsigned long long MODE = NULL_MODE) { mask &= ~MODE;}
		void change_mode(unsigned long long MODE1, unsigned long long MODE2) { mask = mask & ~MODE1 | MODE2;}
		int mode(unsigned long long MODE){ return((mask & MODE) == MODE);} 
public:
//...structure initialization;
		void reset_struct();
		void struct_init ();
		void struct_init (int k, int id_strong);
		void struct_clean(int k, int id_all);
		void struct_clean(int id_all = NULL_STATE);
		void set_blocks  (int N_sm, int m = 0) {
			this->set_elements(N_sm); n = (id_norm = m)+1;
			dim = (int *)new_struct(this->N*sizeof(int));
		}
		void set_dimension(int k, int dim_N) {
			if (k < this->N && 0 <= k && dim) dim[k] = dim_N;
		}
//...operations with solvers;
		int solver(int m) {
			if (id_change == OK_STATE)	return(solver1(m)); else
			if (id_change == TR_STATE)	return(solver2(m)); else 
												return(OK_STATE);
		}
//...constructor and destructor;
		CSolver() {
			mask = n = id_norm = 0;
			pivot_ii = pivot_kk = pivot_ll = dim = NULL;
			id_change = OK_STATE;
		}
		~CSolver(void) {
			delete[] pivot_ii;
			delete[] pivot_kk;
			delete[] pivot_ll;
		}
//...filling block matrix;
		void to_equationHH(int i, int j, T* F, T h);
		void to_equationHH(int i, int j, T* F, int NN, T h);
		void to_equationHL(int i, int j, T* F, T h);
		void to_equationEH(int i, int j, T* F, T h);
		void to_equationHD(int i, int j, T* F);
		void to_equationHD(int i, T* F1, T* F2, double f = 1.);
		void to_equationDD(int i, T* F1, T* F2, double f = 1.);
		void to_equationEE(int i, T* F1, T* F2, double f = 1.);
		void to_equationER(int i, T* F1, T* F2, double f = 1.);
		void to_equationEL(int i, T* F1, T* F2, double f = 1.);
		void to_transferTR(int i, int j, T* F1, T* F2, double f = 1.);
		void to_transferTL(int i, int j, T* F1, T* F2, double f = 1.);
		void to_transferTT(int i, int j, T* F1, T* F2, double f = 1.);
		void to_transferTD(int i, int j, T* F1, T* F2, double f = 1.);
		void to_transferDD(int i, int j, T* F1, T* F2, double f = 1.);
//...assembling block matrix;
		void to_presenceTL();
		void to_presenceTL(int i);
		void diagonal(int i, double f, int id_all = OK_STATE);
		void diagonal(int id_all = OK_STATE);
		void lagrangian(int i, double f = 1., double fEE = 1., int id_all = OK_STATE);
		void lagrangian(		  double f = 1., double fEE = 1., int id_all = OK_STATE);
		void admittance(int i, int m, T adm_dd);
		void admittance(int i, int m, T adm_dd, int k, T adm_pp);
		void admittance(T* d, T* p, int M, T adm_dd, T adm_pp);
public:
		void pivot_init(int dim_N) {
			delete[] pivot_ii; pivot_ii = new int[dim_N];
			delete[] pivot_kk; pivot_kk = new int[dim_N];
			delete[] pivot_ll; pivot_ll = new int[dim_N];
		}
//...inversing square matrix;
		int GaussJ   (T ** A, int dim_N);
		void trf_step(T ** A, int dim_N, T * h, T * X);
//...inversing diagonal block on the row;
		int GaussJ(int k);
		void hrf_step_solver(int k, int m);
//...converting block matrix to the extern solver;
		void Blocks_Row	 (int _BlkRow, T * _defc, double alpha);
		void Right_Handside(int _BlkRow, T * _refc, double alpha);
		void Initial_Guess (int _BlkRow, T * _refc);
//...testing elements of block matrix;
		void test_struct(FILE * id_STRU, int id_pernum);
		void test_struct(char * ch_STRU, int id_pernum = 1);
		void test_struct(const char * ch_STRU, int id_pernum = 1);
		void test_gram_matrix      (FILE * TST, int k, int variant = ERR_STATE, double eps = EE_ker);
		void test_energy_matrix    (FILE * TST, int k, int variant = ERR_STATE, double eps = EE_ker);
		void test_transfer_matrix  (FILE * TST, int k, double eps = EE_ker);
		void test_gram_symmetry    (FILE * TST, int k, double eps = EE_ker);
		void test_energy_symmetry  (FILE * TST, int k, double eps = EE_ker);
		void test_transfer_symmetry(FILE * TST, int k, double eps = EE_ker);
};

///////////////////////////////////////////////
//...restore initial state of sparse structure;
template <typename T>
void CSolver<T>::reset_struct()
{
	if (this->JR)
	for (int k = 0; k < this->N; k++) if (this->JR[k]) {
		if (this->TR && this->TR[k]) //...уничтожаем дополнительно порожденные блоки;
		for (int j = this->JR[k][this->JR_DIAG]-this->JR_SHIFT+1; j < this->JR[k][0]; j++) this->TR[k][j].set_matrix(0, 0);

		if (this->TL && this->TL[k] && id_change != TR_STATE) //...в симметричном солвере массив TL не порождаем!!!;
		for (int j = this->JR[k][this->JR_DIAG]-this->JR_SHIFT+1; j < this->JR[k][0]; j++) this->TL[k][j].set_matrix(0, 0);
	}
}

//////////////////////////////////////////////////////////
//...blocks structure initialization of auxilliary arrays;
template <typename T>
void CSolver<T>::struct_init()
{
	for (int k = 0; k < this->N; k++) 
	struct_init(k, OK_STATE);
}

/////////////////////////////////////////////////
//...one row initialization of auxilliary arrays;
template <typename T>
void CSolver<T>::struct_init(int k, int id_strong)
{//...последовательность: h, hL, hEE, HD;
	if (! dim || ! this->hh) return;
	if (id_strong || ! this->hh[k][0].GetMatrix()) this->hh[k][0].set_matrix(n, dim[k]);
	if (id_strong || ! this->hh[k][2].GetMatrix()) this->hh[k][2].set_matrix(id_norm, dim[k]);
}

/////////////////////////////////////////////////
//...cleaning auxilliary arrays in one block row;
template <typename T>
void CSolver<T>::struct_clean(int k, int id_all)
{
	if (! dim || ! this->hh) return;
	if (id_all == ADDITIONAL_STATE) this->hh[k][3].clean_matrix(dim[k], dim[k]);
	else {
		this->hh[k][0].clean_matrix(n, dim[k]);
		this->hh[k][1].clean_matrix(id_norm, dim[k]);
		this->hh[k][2].clean_matrix(id_norm, dim[k]);

		if ((id_all != NULL_STATE) && this->JR && this->JR[k]) {
			if (this->TR && this->TR[k]) for (int i = 0; i < this->JR[k][0]; i++)
				this->TR[k][i].clean_matrix(dim[k], dim[this->JR[k][i+this->JR_SHIFT]]);
			if (this->TL && this->TL[k]) for (int i = 0; i < this->JR[k][0]; i++)
				this->TL[k][i].clean_matrix(dim[k], dim[this->JR[k][i+this->JR_SHIFT]]);
			if (this->TT && this->TT[k]) for (int i = 0; i <= this->JR[k][this->JR_DIAG]-this->JR_SHIFT; i++) //...определяем длину массива по положению диагонали;
				this->TT[k][i].clean_matrix(dim[k], dim[k]);
			if (this->TD && this->TD[k]) for (int i = 0; i <= this->JR[k][this->JR_DIAG]-this->JR_SHIFT; i++)
				this->TD[k][i].clean_matrix(dim[this->JR[k][i+this->JR_SHIFT]], dim[this->JR[k][i+this->JR_SHIFT]]);
			if (id_all != ZERO_STATE) 
				this->hh[k][3].clean_matrix(dim[k], dim[k]);
		}
	}
}

//////////////////////////////////////////////////
//...cleaning auxilliary arrays in block strucure;
template <typename T>
void CSolver<T>::struct_clean(int id_all)
{
	for (int k = 0; k < this->N; k++) 
	struct_clean(k, id_all);
}

/////////////////////////////////////////////////////////
//              FILLING BLOCK MATRIXES                 //
/////////////////////////////////////////////////////////
template <typename T>
void CSolver<T>::to_equationHH(int i, int j, T* F, T h)
{
	if (mode(NO_TR) || ! this->hh[i][0].GetMatrix()) return;	int NN = dim[i]; T* m = this->hh[i][0][j];
	for (int k = 0; k < NN; k++) m[k] += F[k]*h;
}

template <typename T>
void CSolver<T>::to_equationHH(int i, int j, T* F, int NN, T h)
{
	if (mode(NO_TR) || ! this->hh[i][0].GetMatrix()) return;	T* m = this->hh[i][0][j];
	for (int k = 0; k < NN; k++) m[k] += F[k]*h;
}

template <typename T>
void CSolver<T>::to_equationHL(int i, int j, T* F, T h)
{
	if (mode(NO_TL) || ! this->hh[i][0].GetMatrix()) return;	int NN = dim[i];
	if (! this->hh[i][1].GetMatrix()) this->hh[i][1].set_matrix(id_norm, NN); T* m = this->hh[i][1][j]; /*было: T* m = this->hh[i][0][j]; - исправлено 09.06.2015*/
	for (int k = 0; k < NN; k++) m[k] += F[k]*h;
}

template <typename T>
void CSolver<T>::to_equationEH(int i, int j, T* F, T h)
{
	if (mode(NO_TR) || ! this->hh[i][2].GetMatrix() || j >= id_norm) return; int NN = dim[i]; T* m = this->hh[i][2][j];
	for (int k = 0; k < NN; k++) m[k] += F[k]*h;
}

template <typename T>
void CSolver<T>::to_equationHD(int i, int j, T* F)
{
	if (mode(NO_TR) || ! this->hh[i][0].GetMatrix() || ! this->hh[i][3].GetMatrix()) return; int NN = dim[i]; T** m = this->hh[i][3].GetMatrix(), * m2 = this->hh[i][2][j]; 
	for (int k = 0; k < NN; k++)
	for (int l = 0; l < NN; l++) m2[k] += m[k][l]*F[l];
}

template <typename T>
void CSolver<T>::to_equationHD(int i, T* F1, T* F2, double f)
{
	if (mode(NO_TR) || ! this->hh[i][0].GetMatrix()) return; int NN = dim[i];
	if (! this->hh[i][3].GetMatrix()) this->hh[i][3].set_matrix(NN, NN); T** m = this->hh[i][3].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int	l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_equationDD(int i, T* F1, T* F2, double f)
{
	if (mode(NO_TR) || ! this->JR[i][this->JR_DIAG]) return;	int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN = dim[i];
	if (! this->TR[i][j_diag].GetMatrix()) this->TR[i][j_diag].set_matrix(NN, NN); T** m = this->TR[i][j_diag].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int	l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_equationEE(int i, T* F1, T* F2, double f)
{
	if (mode(NO_TR) || ! this->JR[i][this->JR_DIAG]) return;	int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN = dim[i];
	if (! this->TL[i][j_diag].GetMatrix()) this->TL[i][j_diag].set_matrix(NN, NN); T** m = this->TL[i][j_diag].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int	l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_equationER(int i, T* F1, T* F2, double f)
{
	if (mode(NO_TR) || ! this->JR[i][this->JR_DIAG]) return;	int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN = dim[i];
	if (! this->TT[i][j_diag].GetMatrix()) this->TT[i][j_diag].set_matrix(NN, NN); T** m = this->TT[i][j_diag].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int	l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_equationEL(int i, T* F1, T* F2, double f)
{
	if (mode(NO_TL) || ! this->JR[i][this->JR_DIAG]) return;	int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN = dim[i];
	if (! this->TD[i][j_diag].GetMatrix()) this->TD[i][j_diag].set_matrix(NN, NN); T** m = this->TD[i][j_diag].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int	l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_transferTR(int i, int j, T* F1, T* F2, double f)
{
	if (mode(NO_TR)) return; int NN1 = dim[i], NN2 = dim[this->JR[i][j+this->JR_SHIFT]];
	if (! this->TR[i][j].GetMatrix()) this->TR[i][j].set_matrix(NN1, NN2); T** m = this->TR[i][j].GetMatrix();
	for (int k = 0; k < NN1; k++)
	for (int	l = 0; l < NN2; l++) m[k][l] -= F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_transferTL(int i, int j, T* F1, T* F2, double f)
{
	if (mode(NO_TL)) return; int NN1 = dim[i], NN2 = dim[this->JR[i][j+this->JR_SHIFT]];
	if (! this->TL[i][j].GetMatrix()) this->TL[i][j].set_matrix(NN1, NN2); T** m = this->TL[i][j].GetMatrix();
	for (int k = 0; k < NN1; k++)
	for (int	l = 0; l < NN2; l++) m[k][l] -= F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_transferTT(int i, int j, T* F1, T* F2, double f)
{
	if (mode(NO_TR)) return; int NN = dim[i];
	if (! this->TT[i][j].GetMatrix()) this->TT[i][j].set_matrix(NN, NN); T** m = this->TT[i][j].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_transferTD(int i, int j, T* F1, T* F2, double f)
{
	if (mode(NO_TL)) return; int NN = dim[this->JR[i][j+this->JR_SHIFT]];
	if (! this->TD[i][j].GetMatrix()) this->TD[i][j].set_matrix(NN, NN); T** m = this->TD[i][j].GetMatrix();
	for (int k = 0; k < NN; k++)
	for (int l = 0; l < NN; l++) m[k][l] += F1[k]*F2[l]*f;
}

template <typename T>
void CSolver<T>::to_transferDD(int i, int j, T* F1, T* F2, double f)
{
	to_transferTT(i, j, F1, F1, f);
	to_transferTD(i, j, F2, F2, f);
}

template <typename T>
void CSolver<T>::to_presenceTL()
{
	for (int i = 0; i < this->N; i++) to_presenceTL(i);
}

template <typename T>
void CSolver<T>::to_presenceTL(int i)
{
	int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN1 = dim[i], NN2;
	for (int j = 0; j < j_diag; j++) //...определяем длину массива по положению диагонали;
	if (this->TR[i][j].GetMatrix() && ! this->TL[i][j].GetMatrix() && (NN2 = dim[this->JR[i][j+this->JR_SHIFT]])) {
		this->TL[i][j].set_matrix(NN1, NN2);  T** m1 = this->TL[i][j].GetMatrix(), ** m2 = this->TR[i][j].GetMatrix(); 
		for (int k = 0; k < NN1; k++)
			memcpy(m1[k], m2[k], NN2*sizeof(T));
	}
}

template <typename T>
void CSolver<T>::diagonal(int i, double f, int id_all)
{
	if (this->JR[i][this->JR_DIAG]) {
		int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, j, k, m, l, NN = dim[i];
		if (id_all) {
			if (! this->TR[i][j_diag].GetMatrix()) this->TR[i][j_diag].set_matrix(NN, NN);
			if (! this->TL[i][j_diag].GetMatrix()) this->TL[i][j_diag].set_matrix(NN, NN);
			for (k = 0; k < this->p[this->N+i]; k++) {
				for (j = this->JR[this->p[k]][this->JR_DIAG]-this->JR_SHIFT+1; j > 0; j--) //...определяем длину массива по положению диагонали;
				if (this->JR[this->p[k]][j+this->JR_SHIFT-1] == i) break;

				if (j && this->TD[this->p[k]][j-1].GetMatrix())
				for (m = 0; m < NN; m++)
				for (l = 0; l < NN; l++) this->TR[i][j_diag][m][l] += this->TD[this->p[k]][j-1][m][l]*f;
			}
			for (j = 0; j <= this->JR[i][this->JR_DIAG]-this->JR_SHIFT; j++) if (j != j_diag) {
				if (this->TT[i][j].GetMatrix())
				for (m = 0; m < NN; m++)
				for (l = 0; l < NN; l++)	this->TR[i][j_diag][m][l] += this->TT[i][j][m][l]*f;
			}
			else {
				if (this->TT[i][j].GetMatrix())
				for (m = 0; m < NN; m++)
				for (l = 0; l < NN; l++)	this->TL[i][j_diag][m][l] += this->TT[i][j][m][l]*f;

				if (this->TD[i][j].GetMatrix())
				for (m = 0; m < NN; m++)
				for (l = 0; l < NN; l++)	this->TL[i][j_diag][m][l] += this->TD[i][j][m][l]*f;
			}
		}
		if (this->hh[i][1].GetMatrix())
		for (m = 0; m < id_norm; m++)
		for (l = 0; l < NN; l++) this->hh[i][0][m][l] += this->hh[i][1][m][l]*f;
	}
}

template <typename T>
void CSolver<T>::diagonal(int id_all)
{
	for (int i = 0; i < this->N; i++) if (this->JR[i][this->JR_DIAG]) {
		int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, k_diag, k, j, NN1 = dim[i], NN2;
		if (id_all) {
			if (! this->TR[i][j_diag].GetMatrix()) this->TR[i][j_diag].set_matrix(NN1, NN1);
			if (! this->TL[i][j_diag].GetMatrix()) this->TL[i][j_diag].set_matrix(NN1, NN1);
			for (  j = 0; j <= this->JR[i][this->JR_DIAG]-this->JR_SHIFT; j++) //...определяем длину массива по положению диагонали;
			if (this->JR[k = this->JR[i][j+this->JR_SHIFT]][this->JR_DIAG] && (NN2 = dim[k])) {
				if (! this->TR[k][k_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT].GetMatrix()) this->TR[k][k_diag].set_matrix(NN2, NN2);
				if (j != j_diag) { 
					if (this->TT[i][j].GetMatrix())
					for (int m = 0; m < NN1; m++)
					for (int l = 0; l < NN1; l++)	this->TR[i][j_diag][m][l] += this->TT[i][j][m][l];

					if (this->TD[i][j].GetMatrix())
					for (int m = 0; m < NN2; m++)
					for (int l = 0; l < NN2; l++)	this->TR[k][k_diag][m][l] += this->TD[i][j][m][l];
				}
				else {
					if (this->TT[i][j].GetMatrix())
					for (int m = 0; m < NN1; m++)
					for (int l = 0; l < NN1; l++)	this->TL[i][j_diag][m][l] += this->TT[i][j][m][l];

					if (this->TD[i][j].GetMatrix())
					for (int m = 0; m < NN1; m++)
					for (int l = 0; l < NN1; l++)	this->TL[i][j_diag][m][l] += this->TD[i][j][m][l];
				}
			}
		}
		if (this->hh[i][1].GetMatrix())
		for (int m = 0; m < id_norm; m++)
		for (int l = 0; l < NN1; l++)	this->hh[i][0][m][l] += this->hh[i][1][m][l];
	}
}

template <typename T>
void CSolver<T>::lagrangian(int i, double f, double fEE, int id_all)
{
	if (! this->JR[i][this->JR_DIAG]) return; int j_diag = this->JR[i][this->JR_DIAG]-this->JR_SHIFT, NN1 = dim[i];
	for (int m = 0; m < NN1; m++) {
		if (this->hh[i][0].GetMatrix() && this->hh[i][2].GetMatrix())
		for (int l = 0; l < id_norm; l++) { 
			this->hh[i][0][l][m] *= f;
			this->hh[i][0][l][m] += fEE*this->hh[i][2][l][m];
		}
		if (id_all) {
			for (int j = 0; j < j_diag; j++) { //...определяем длину массива по положению диагонали;
				if (this->TR[i][j].GetMatrix())
				for (int NN2 = dim[this->JR[i][j+this->JR_SHIFT]], l = 0; l < NN2; l++) this->TR[i][j][m][l] *= f;
				if (this->TL[i][j].GetMatrix())
				for (int NN2 = dim[this->JR[i][j+this->JR_SHIFT]], l = 0; l < NN2; l++) this->TL[i][j][m][l] *= f;
			}
			if (this->TR[i][j_diag].GetMatrix() && this->TL[i][j_diag].GetMatrix())
			for (int l = 0; l < NN1; l++) { 
				this->TR[i][j_diag][m][l] *= f;
				this->TR[i][j_diag][m][l] += fEE*this->TL[i][j_diag][m][l];
			}
		}
	}
}

template <typename T>
void CSolver<T>::lagrangian(double f, double fEE, int id_all)
{
	for (int i = 0; i < this->N; i++) lagrangian(i, f, fEE, id_all);
}

template <typename T>
void CSolver<T>::admittance(int i, int m, T adm_dd)
{
	m += id_norm;
	for (int NN = dim[i], l = 0; l < NN; l++) this->hh[i][0][m][l] *= adm_dd;
}

template <typename T>
void CSolver<T>::admittance(int i, int m, T adm_dd, int k, T adm_pp)
{
	m += id_norm;
	k += id_norm;
	for (int NN = dim[i], l = 0; l < NN; l++) this->hh[i][0][m][l] = adm_dd*this->hh[i][0][m][l]+adm_pp*this->hh[i][0][k][l];
}

template <typename T>
void CSolver<T>::admittance(T* d, T* p, int M, T adm_dd, T adm_pp)
{
	if (! d) return;
	for (int k = 0; k < M; k++) d[k] *= adm_dd;
	if (! p) return;
	for (int k = 0; k < M; k++) d[k] += adm_pp*p[k];
}

/////////////////////////////////////////////
//...inversing matrix by Gauss-Jordan method;
template <typename T>
int CSolver<T>::GaussJ(T ** A, int dim_N)
{
	if (! A || ! pivot_ii || ! pivot_kk || ! pivot_ll) return(0);
	int i, k, l, k0, l0, id_err = 0;
	memset(pivot_ii, 0, dim_N*sizeof(int));
 
///////////////////////////////////////
//...main loop for inversing of matrix;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;

//////////////////////////////////////////
//...look for position of maximal element;
		for (k = 0; k < dim_N; k++)
			if (pivot_ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! pivot_ii[l]) {
						double f_match = to_double(fabs(A[k][l]));
						if (f_match >= f) { 
							f = f_match; k0 = k; l0 = l; 
						}
					}
					else 
						if (pivot_ii[l] > 1) {
							id_err = -2; return(0);
						}

		++(pivot_ii[l0]);

////////////////////////////////////////////////////////////
//...row swapping for diagonal position for maximal element;
		if (k0 != l0) {T * temp_A = A[k0]; A[k0] = A[l0]; A[l0] = temp_A; }

///////////////////////////////////////////////////////////////////////////
//...storing position of maximal element and normalization by them all row;
		pivot_kk[i] = k0; 
		pivot_ll[i] = l0;
		if (A[l0][l0] == T(0.)) {
		  id_err = -3; return(0);
		}

////////////////////////////////
//...diagonal row normalization;
		T finv = 1./A[l0][l0]; A[l0][l0] = T(1.);
		for (l = 0; l < dim_N; l++) A[l0][l] *= finv;

/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = A[k][l0]; A[k][l0] = T(0.);
				for (l = 0; l < dim_N; l++) A[k][l] -= A[l0][l]*finv;
			}
	}
	l = id_err;

////////////////////////////////////////////////////////////////////////////////
//...reverse sorting of columns of inverse matrix and memory release and return;
	for (l = dim_N-1; l >= 0; l--)
		if (pivot_kk[l] != pivot_ll[l])
			for (k = 0; k < dim_N; k++) {
				T temp_A = A[k][pivot_kk[l]]; A[k][pivot_kk[l]] = A[k][pivot_ll[l]]; A[k][pivot_ll[l]] = temp_A;
			}
	return(1);
}

/////////////////////////////////
//...умножение матрицы на вектор;
template <typename T>
void CSolver<T>::trf_step(T ** A, int dim_N, T * h, T * X)
{
	int i, l;
	if  (A && h && X)
	for (					 i = 0; i < dim_N; i++)
	for (X[i] = T(0.), l = 0; l < dim_N; l++) X[i] += A[i][l]*h[l];
}

template <typename T>
int CSolver<T>::GaussJ(int k)
{
	if (k < this->N && 0 <= k && this->JR && this->JR[k] && this->JR[k][this->JR_DIAG] && this->TR && this->TR[k] && dim) {
		pivot_init(dim[k]);
		return GaussJ(this->TR[k][this->JR[k][this->JR_DIAG]-this->JR_SHIFT].GetMatrix(), dim[k]);
	}
	return(0);
}

template <typename T>
void CSolver<T>::hrf_step_solver(int k, int m)
{
	if (k < this->N && 0 <= k  && m < id_norm && 0 <= m && this->JR && this->JR[k] && this->JR[k][this->JR_DIAG] &&
		this->TR && this->TR[k] && this->hh && this->hh[k] && this->hh[k][0].GetMatrix() && dim) {
/////////////////////
//...находим решение;
		trf_step(this->TR[k][this->JR[k][this->JR_DIAG]-this->JR_SHIFT].GetMatrix(), dim[k], this->hh[k][0].GetMatrix()[m], this->hh[k][0].GetMatrix()[id_norm]);
///////////////////////////////////////////////////////////////
//...переставляем строчки матрицы правых частей;
		T * temp_h = this->hh[k][0].GetMatrix()[m]; 
		this->hh[k][0].GetMatrix()[m] = this->hh[k][0].GetMatrix()[id_norm]; 
		this->hh[k][0].GetMatrix()[id_norm] = temp_h;
	}
}

////////////////////////////////////////////////////////////
//...заполнение строки блочной матрицы для внешнего солвера;
template <typename T>
void CSolver<T>::Blocks_Row(int _BlkRow, T * _defc, double alpha)
{
	int i, j, k, l, m, ii;
	double beta = 1.-1./alpha = max(1., alpha);
	diagonal(this->p[_BlkRow], 1.);
	for (m = 0, i = 0; i < this->IL[this->p[_BlkRow]][0]; i++) {
		for (j = this->JR[ii = this->IL[this->p[_BlkRow]][i+this->JR_SHIFT]][0]; j > 0; j--)
		if ( ii != this->p[_BlkRow] && this->JR[ii][j+this->JR_SHIFT-1] == this->p[_BlkRow]) break;
		if (j)
		for (k = 0; k < dim[ii]; k++)
		for (l = 0; l < dim[this->p[_BlkRow]]; l++) 
			_defc[m++] = this->TL[ii][j-1][k][l]*alpha;
	}
	for (j = 0; j < this->JR[this->p[_BlkRow]][0]; j++) {
		if (j != this->JR[this->p[_BlkRow]][this->JR_DIAG]) {
			for (l = 0; l < dim[this->JR[this->p[_BlkRow]][j+this->JR_SHIFT]]; l++)
			for (k = 0; k < dim[this->p[_BlkRow]]; k++) 
				_defc[m++] =  this->TR[this->p[_BlkRow]][j][k][l]*alpha;
		}
		else {
			for (l = 0; l < dim[this->p[_BlkRow]]; l++)
			for (k = 0; k < dim[this->p[_BlkRow]]; k++) 
				_defc[m++] =  this->TR[this->p[_BlkRow]][j][k][l]*alpha+this->TL[this->p[_BlkRow]][j][k][l]*beta;
		}
	}
}

/////////////////////////////////////////////
//...заполнение правой части внешней матрицы;
template <typename T>
void CSolver<T>::Right_Handside(int _BlkRow, T * _refc, double alpha)
{
	double beta = 1.-1./(alpha = max(1., alpha));
	for (int k = 0; k < dim[this->p[_BlkRow]]; k++) 
		_refc[k] = this->hh[this->p[_BlkRow]][0][0][k]*alpha+this->hh[this->p[_BlkRow]][2][0][k]*beta;
}

/////////////////////////////////////
//...начальное приближение к решению;
template <typename T>
void CSolver<T>::Initial_Guess(int _BlkRow, T * _refc)
{
	for (int k = 0; k < dim[this->p[_BlkRow]]; k++) 
		_refc[k] = this->hh[this->p[_BlkRow]][0][id_norm][k];
}

////////////////////////////////////////////////////////////////////////
//...отображение структуры блочной матрицы (т.е. связей блоков) на диск;
template <typename T>
void CSolver<T>::test_struct(FILE * id_STRU, int id_pernum)
{
  if (! id_STRU) return;

  fprintf(id_STRU, "\nStructure of the block matrix (N = %i):", this->N);
  if (! this->JR) return; 

////////////////////////////////////////////////////////////////
//...заполняем матрицу на диске, отображающую структуру образца;
  int k_par, i_par, k, i, j;
  for (                        k_par = 0; k_par < this->N; k_par++)
  for (fprintf(id_STRU, "\n"), i_par = 0; i_par < this->N; i_par++) 
  if  (this->JR[k = id_pernum ? this->p[k_par] : k_par]) {
       i    = id_pernum ? this->p[i_par] : i_par;
		 if (this->p[this->N+i] < this->p[this->N+k]) {
			 for (j = this->JR[i][0]; j > 0; j--) 
				 if (this->JR[i][j+this->JR_SHIFT-1] == k) break;
		 }
		 else {
			 for (j = this->JR[k][0]; j > 0; j--) 
				 if (this->JR[k][j+this->JR_SHIFT-1] == i) break;
		 }
       if  (j) fprintf(id_STRU, " X");
       else    fprintf(id_STRU, /*" _"*/"  ");
  }
  fprintf(id_STRU, "\n");
}

////////////////////////////////
//...строковые варианты функции;
template <typename T>
void CSolver<T>::test_struct(char * ch_STRU, int id_pernum)
{
  if (! ch_STRU) return;
  FILE * id_STRU = fopen(ch_STRU, "w");
	test_struct(id_STRU, id_pernum);
  fclose(id_STRU);
}
template <typename T>
void CSolver<T>::test_struct(const char * ch_STRU, int id_pernum)
{
  if (! ch_STRU) return;
  FILE * id_STRU = fopen(ch_STRU, "w");
	test_struct(id_STRU, id_pernum);
  fclose(id_STRU);
}

////////////////////////////////////////////////////////////
//...testing output of diagonal matrix with right-hand part;
template <typename T>
void CSolver<T>::test_gram_matrix(FILE * TST, int k, int variant, double eps)
{
	if (! dim || k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! this->JR[k][this->JR_DIAG] || ! this->hh || ! TST) return;
	int NN = dim[k], i, j, j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT;

	fprintf(TST, "\nGram matrix for %d block (loc_%d):", k, this->p[this->N+k]);

	if (variant != END_STATE)
	for (i = 0; this->TR[k][j_diag].GetMatrix() && i < NN; i++) {
		fprintf(TST, "\n");
		if (mode(MASKS_MODE)) { //...печать маски подматриц-блоков;
			for (j = 0; j < NN; j++)
				if (fabs(this->TR[k][j_diag][i][j]) < eps) fprintf(TST, " _"); 
				else fprintf(TST, " X");
		}
		else { //...печать значений подматриц-блоков;
			for (j = 0; j < NN; j++) 
				if (typeid(T) == typeid(complex))
					fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->TR[k][j_diag][i][j], eps), filtr_Im(this->TR[k][j_diag][i][j], eps)); else
					fprintf(TST, "%0.15lg ", filtr_Re(this->TR[k][j_diag][i][j], eps));
		}
	}	fprintf(TST, "\n");  

	if (variant == END_STATE) {
		fprintf(TST, "\nEigenvalues vectors (Re, Im):");
		if (mode(FULLY_MODE)) { //...полная печать собственных значений;
			for (i = 0; i < NN; i++) {
				fprintf(TST, "\n");
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-1][i], eps));
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-2][i], eps));
			}  fprintf(TST, "\n");
		}
		else {
			for (i = 0; i < 3; i++) {
				fprintf(TST, "\n");
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-1][i], eps));
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-2][i], eps));
			}
			for (i = NN-3; i < NN; i++) {
				fprintf(TST, "\n");
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-1][i], eps));
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][n-2][i], eps));
			}  fprintf(TST, "\n");
		}
	}
	else
	if (0 <= variant && variant < n) { //...печать одного варианта правых частей;
		fprintf(TST, "\nRow vectors:");
		for (i = 0; this->hh[k][0].GetMatrix() && i < NN; i++) {
			if (typeid(T) == typeid(complex))
				fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->hh[k][0][variant][i], eps), filtr_Im(this->hh[k][0][variant][i], eps)); else
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][variant][i], eps));
		}  fprintf(TST, "\n");
	}
	else
	if (variant >= n) { //...печать всех вариантов правых частей;
		fprintf(TST, "\nRow vectors:");
		for (i = 0; i < NN; i++) {
			fprintf(TST, "\n");
			if (mode(MASKS_MODE)) { //...печать маски правых частей;
				for (j = 0; j < id_norm; j++)
					if (fabs(this->hh[k][0][j][i]) < eps) fprintf(TST, " _"); 
					else fprintf(TST, " X");
			}
			else { //...печать значений правых частей;
				for (j = 0; j < id_norm; j++)
					if (typeid(T) == typeid(complex))
						fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->hh[k][0][j][i], eps), filtr_Im(this->hh[k][0][j][i], eps)); else
						fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][0][j][i], eps));
			}
		}  fprintf(TST, "\n");
	}
}

/////////////////////////////////////
//...testing output of energy matrix;
template <typename T>
void CSolver<T>::test_energy_matrix(FILE * TST, int k, int variant, double eps)
{
	if (k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! this->JR[k][this->JR_DIAG] || ! dim || ! TST) return;
	int j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT, j;

	fprintf(TST, "\nEnergy matrix for %d block (loc_%d):", k, this->p[this->N+k]);
	if (this->TL[k][j_diag].GetMatrix()) for (int i = 0; i < dim[k]; i++) {
		fprintf(TST, "\n");
		if (mode(MASKS_MODE)) { //...печать маски матрицы энергии;
			for (j = 0; j < dim[k]; j++) 
				if (fabs(this->TL[k][j_diag][i][j]) < eps) fprintf(TST, " _"); 
				else fprintf(TST, " X");
		}
		else { //...печать значений элементов матрицы энергии;
			for (j = 0; j < dim[k]; j++) 
				if (typeid(T) == typeid(complex))
					fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->TL[k][j_diag][i][j], eps), filtr_Im(this->TL[k][j_diag][i][j], eps)); else
					fprintf(TST, "%0.15lg ", filtr_Re(this->TL[k][j_diag][i][j], eps));
		}
	}	fprintf(TST, "\n");  

	if (0 <= variant && variant < n && this->hh && this->hh[k]) { //...печать одного варианта правых частей;
		fprintf(TST, "\nRow vectors:");
		if (this->hh[k][2].GetMatrix()) for (int i = 0; i < dim[k]; i++) {
			if (typeid(T) == typeid(complex))
				fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->hh[k][2][variant][i], eps), filtr_Im(this->hh[k][2][variant][i], eps)); else
				fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][2][variant][i], eps));
		}  fprintf(TST, "\n");
	}
	else
	if (variant >= n && this->hh && this->hh[k] && this->hh[k][2].GetMatrix()) {
		fprintf(TST, "\nRow vectors:");
		for (int i = 0; i < dim[k]; i++) {
			fprintf(TST, "\n");
			if (mode(MASKS_MODE)) { //...печать маски правых частей;
				for (j = 0; j < id_norm; j++)
					if (fabs(this->hh[k][2][j][i]) < eps) fprintf(TST, " _"); 
					else fprintf(TST, " X");
			}
			else { //...печать значений правых частей;
				for (j = 0; j < id_norm; j++)
					if (typeid(T) == typeid(complex))
						fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(this->hh[k][2][j][i], eps), filtr_Im(this->hh[k][2][j][i], eps)); else
						fprintf(TST, "%0.15lg ", filtr_Re(this->hh[k][2][j][i], eps));
			}
		}  fprintf(TST, "\n");
	}
}

///////////////////////////////////////
//...testing output of transfer matrix;
template <typename T>
void CSolver<T>::test_transfer_matrix(FILE * TST, int k, double eps)
{
	if (k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! this->JR[k][this->JR_DIAG] || ! dim || ! TST) return;
	int i, j, j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT, m, l;
	for (j = 0; j < this->JR[k][0]; j++) if (j != j_diag) {
		T ** tt;
		if (this->TR && this->TR[k]) {
			fprintf(TST, "\nRight transfer matrix for (%d, %d) blocks:", k, i = this->JR[k][j+this->JR_SHIFT]);
			for (tt = this->TR[k][j].GetMatrix(), m = 0; tt && m < dim[k]; m++) {
				fprintf(TST, "\n");
				if (mode(MASKS_MODE)) { //...печать маски подматриц-блоков;
					for (l = 0; l < dim[i]; l++)
						if (fabs(tt[m][l]) < eps) fprintf(TST, " _"); 
						else fprintf(TST, " X");
				}
				else { //...печать значений подматриц-блоков;
					for (l = 0; l < dim[i]; l++)
						if (typeid(T) == typeid(complex))
							fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(tt[m][l], eps), filtr_Im(tt[m][l], eps)); else
							fprintf(TST, "%0.15lg ", filtr_Re(tt[m][l], eps));
				}
			}  fprintf(TST, "\n");
		}
		if (this->TL && this->TL[k]) {
			fprintf(TST, "\nLeft transfer matrix for (%d, %d) blocks:", k, i = this->JR[k][j+this->JR_SHIFT]);
			for (tt = this->TL[k][j].GetMatrix(), m = 0; tt && m < dim[k]; m++) {
				fprintf(TST, "\n");
				if (mode(MASKS_MODE)) { //...печать маски подматриц-блоков;
					for (l = 0; l < dim[i]; l++)
						if (fabs(tt[m][l]) < eps) fprintf(TST, " _"); 
						else fprintf(TST, " X");
				}
				else { //...печать значений подматриц-блоков;
					for (l = 0; l < dim[i]; l++)
						if (typeid(T) == typeid(complex))
							fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re(tt[m][l], eps), filtr_Im(tt[m][l], eps)); else
							fprintf(TST, "%0.15lg ", filtr_Re(tt[m][l], eps));
				}
			}  fprintf(TST, "\n");
		}
	}
}

///////////////////////////////////////////////////
//...testing output of symmetry of diagonal matrix;
template <typename T>
void CSolver<T>::test_gram_symmetry(FILE * TST, int k, double eps)
{
	if (k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! this->JR[k][this->JR_DIAG] || ! dim || ! TST) return;
	int j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT, i, j;

	fprintf(TST, "\nSymmetry of gram matrix for %d block (loc_%d):", k, this->p[this->N+k]);
	T ** tt;
	for (tt = this->TR[k][j_diag].GetMatrix(), i = 0; tt && i < dim[k]; i++) {
		fprintf(TST, "\n");
		for (j = 0; j < dim[k]; j++) 
			if (typeid(T) == typeid(complex))
				fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re((tt[i][j]-tt[j][i])*.5, eps), filtr_Im((tt[i][j]-tt[j][i])*.5, eps)); else
				fprintf(TST, "%0.15lg ", filtr_Re((tt[i][j]-tt[j][i])*.5, eps));
	}  fprintf(TST, "\n");
}

///////////////////////////////////////////////////
//...testing output of symmetry of energy matrix;
template <typename T>
void CSolver<T>::test_energy_symmetry(FILE * TST, int k, double eps)
{
	if (k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! this->JR[k][this->JR_DIAG] || ! dim || ! TST) return;
	int j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT, i, j;

	fprintf(TST, "\nSymmetry of energy matrix for %d block (loc_%d):", k, this->p[this->N+k]);
	T ** tt;
	for (tt = this->TL[k][j_diag].GetMatrix(), i = 0; tt && i < dim[k]; i++) {
		fprintf(TST, "\n");
		for (j = 0; j < dim[k]; j++) 
			if (typeid(T) == typeid(complex))
				fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re((tt[i][j]-tt[j][i])*.5, eps), filtr_Im((tt[i][j]-tt[j][i])*.5, eps)); else
				fprintf(TST, "%0.15lg ", filtr_Re((tt[i][j]-tt[j][i])*.5, eps));
	}  fprintf(TST, "\n");
}

//////////////////////////////////////////////////////
//...testing output of symmetry off-diagonal matrixes;
template <typename T>
void CSolver<T>::test_transfer_symmetry(FILE * TST, int k, double eps)
{
	if (k < 0 || k >= this->N || ! this->JR || ! this->JR[k] || ! dim || ! TST) return;
	int i, j, j_diag = this->JR[k][this->JR_DIAG]-this->JR_SHIFT, m, l;
	for (j = 0; j < this->JR[k][0]; j++) if (j != j_diag) {
		T ** tt, ** td;
		if (this->TR && this->TR[k] && this->TL && this->TL[k]) {
			fprintf(TST, "\nTransfer matrix for (%d, %d) blocks:", k, i = this->JR[k][j+this->JR_SHIFT]);
			for (tt = this->TR[k][j].GetMatrix(), td = this->TL[k][j].GetMatrix(), m = 0; tt && td && m < dim[k]; m++) {
				fprintf(TST, "\n");
				for (l = 0; l < dim[i]; l++)
					if (typeid(T) == typeid(complex))
						fprintf(TST, "(%0.15lg,%0.15lg)", filtr_Re((tt[m][l]-td[m][l])*.5, eps), filtr_Im((tt[m][l]-td[m][l])*.5, eps)); else
						fprintf(TST, "%0.15lg ", filtr_Re((tt[m][l]-td[m][l])*.5, eps));
			}  fprintf(TST, "\n");
		}
	}
}
#endif
