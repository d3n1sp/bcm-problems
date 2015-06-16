/*============================================*/
/*                  CCOMPUT2D                 */
/*============================================*/
#ifndef ___CCOMPUT2D___
#define ___CCOMPUT2D___

#include "cdraft.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

/////////////////////////////////////////////////////
//...,базовый класс для счетных схем блочного метода;
template <typename T>
class CComput2D : public CDraft<T> {
protected:
		void comput1(int opt); //...вычислительная схема формирования блочной матрицы для аналитической модели;
		void comput2(int opt); //...вычислительная схема формирования блочной матрицы для конечно-элементной модели;
		void comput3(int opt); //...вычислительная схема (одноблочная, конечно-элементная) для периодической задачи;
		void comput4(int opt); //...счетная схема для задачи Эшелби -- аппроксимация по границе включения;
		void comput5(int opt, T * K, Num_Value _FMF, int id_variant); //...вычислительная схема интегрирования эффективных характеристик по границе блока;
		void comput6(int opt, T * K, Num_Value _FMF, int id_variant); //...вычислительная схема интегрирования эффективных характеристик по объему блока;
		void comput7(int opt, T * K);  //...вычислительная схема дополнительных коррекций интегрируемой величины;
		Num_State comput_kernel1(Num_Comput Num); //...вычислительная схема для блочной структуры; 
		Num_State comput_kernel2(Num_Comput Num); //...вычислительная схема для блочной структуры с продолжением по времени;
//...testing output of block matrix;
		void test_block_matrix (int variant = ERR_STATE);
};

/////////////////////////////////////////////////////////////
//          TEMPLATE VARIANT OF COUNTING SCHEMES           //
/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...основная счетная схема на основе анализа размерных элементов с учетом включения и периодического условия;
template <typename T> // (еще не сделана)==-
void CComput2D<T>::comput1(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), j, k, l, m;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(1);

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(4);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
	double pp[5], Po[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int m_ilist = max(this->N/5, 1), loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++) {
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);
			//this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]].bar->cells_out("CCells");
		}

	for (k = 0; k < this->N; k++) SkeletonBounding(this->B[k], par, NULL_STATE);

////////////////////////////////////////////////
//...discrete norm on stitching sides of blocks;
   if (! this->solver.mode(NO_MESSAGE) && this->volm >= BOUND_VOLUME && this->volm != BLOCK_VOLUME && 
		 ! this->solver.mode(RAPID_MODE)) Message("Discrete stitching sides...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BOUND_VOLUME && ! this->solver.mode(RAPID_MODE)) {

		LOOP_MASK(opt, k);
/////////////////////////////////
//...накапливаем граничные точки;
		if (this->volm != BLOCK_VOLUME)
		for (j = min (this->B[k].bar->graph[0], this->B[k].link[0])-1; j >= 0; j--)
		if (this->B[k].bar->ce[j]->cells_dim() == 1) {

			bnd->zero_grid();
			this->B[k].bar->ce[j]->grid_cells1(bnd, 0., N_max);

			if ((m = this->B[k].link[j+1]) >= 0 && bnd->geom && (! this->solver.mode(NO_TR) || m == this->solver.p[opt]))
			for (l = 0; l < bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num;

				if (bnd->geom[num] == GL_LINE_STRIP) {
					for (; num < num_f-1; num++) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = 0.;
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = 0.;

						gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[0] = gauss_bnd->get_param(0, lp);
							block_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
															 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
			}
//////////////////////////////////
//...формирование блочной матрицы;
			if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
				sprintf(msg, "block %4i: block_bnd->N = %i", k, block_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING)
			this->TransferBB (block_bnd, k, m, 0, E_BASIC_COMPUT); else
			this->TransferBB (block_bnd, k, m, 0,   BASIC_COMPUT);

			if (! this->solver.mode(ACCUMULATION)) block_bnd->add_buffer(block_bnd->N);
		}
   }

///////////////////////////////
//...discrete norm on boundary;
   if (! this->solver.mode(NO_MESSAGE) && 
		this->volm >= BOUND_VOLUME && this->volm != BLOCK_VOLUME) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BOUND_VOLUME) {

		LOOP_MASK(opt, k);
/////////////////////////////////
//...накапливаем граничные точки;
		if (this->volm != BLOCK_VOLUME)
		for (j = min (this->B[k].bar->graph[0], this->B[k].link[0])-1; j >= 0; j--)
		if (this->B[k].bar->ce[j]->cells_dim() == 1) {

			bnd->zero_grid();
			this->B[k].bar->ce[j]->grid_cells1(bnd, 0., N_max);

			if ((m = this->B[k].link[j+1]) < 0 && bnd->geom	&&	! this->solver.mode(NO_TR)) 
			for (l = 0; l < bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num;
				if (bnd->geom[num] == GL_LINE_STRIP) {
					for (; num < num_f-1; num++) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = 0.;
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = 0.;

//////////////////////////////
//...задаем граничные условия;
                  pp[0] = pp[1] = pp[2] = 0.;
						switch (m) {
							case BOUND1_STATE: pp[0] = 1.; pp[1] = 0.; pp[2] = MAX_HIT/*1.*/; break;
							case BOUND2_STATE: pp[0] = 0.; pp[1] = 0.; pp[2] = MAX_HIT; break;
						}
						gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[3] = gauss_bnd->get_param(0, lp);
							bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
															 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
			}
//////////////////////////
//...коррекция квадратуры;
			if (m <= SRF_STATE && this->bar && this->bar->graph && -m+SRF_STATE < this->bar->graph[0])
				bound_bnd->QG_curve(this->bar->ce[-m+SRF_STATE]->mp);

//////////////////////////////////
//...формирование блочной матрицы;
			if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
				sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (this->solv >= ENERGY_SOLVING)
			this->GramAll(bound_bnd, k, 1, E_BASIC_COMPUT); else
			this->GramAll(bound_bnd, k, 1,   BASIC_COMPUT);

			if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
		}
	}

//////////////////////////////////////////////////////////
//...discrete norm on initial volume (еще не сделано !!!);
	if (! this->solver.mode(NO_MESSAGE) && this->volm >= BLOCK_VOLUME) Message("Discrete initial volume...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BLOCK_VOLUME) {
		
		LOOP_MASK(opt, k);
	}
   delete block_bnd;
   delete bound_bnd;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...счетная схема для конечно-элементной модели с учетом включения и периодического условия;
template <typename T>
void CComput2D<T>::comput2(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), 
		  id_dir, id_first, i, j, k, l, m, mm, elem, id_fast = ZERO_STATE;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * bnd = CreateNodes();

	CGrid * block_phs = CreateNodes(GRID_QG_NODES);
			  block_phs->add_params(2);

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(1);

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params((this->solv%ENERGY_SOLVING) ? 5 : 4);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
	int m_ilist = max(this->N/5, 1), loop, N_ini = 4;
	double pp[5], Po[6], pm[3*4], par[6]; this->SetGeomBounding(par);

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++) {
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);
			//this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]].bar->cells_out("CCells");
		}

////////////////////////////////
//...discrete norm on inclusion;
	if (! this->solver.mode(NO_MESSAGE) && this->volm >= BOUND_VOLUME && this->volm != BLOCK_VOLUME && 
		 ! this->solver.mode(RAPID_MODE)) Message("Discrete inclusion...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BOUND_VOLUME && ! this->solver.mode(RAPID_MODE)) {
		
		LOOP_MASK(opt, k);
		if (this->volm != BLOCK_VOLUME)
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			if (this->NUM_PHASE && this->B[k].link[0] > this->NUM_PHASE)
			for (j = 0;  j < this->B[k].bar->ce[i]->graph[1]; j++) 
			for (m = this->NUM_PHASE; m < this->B[k].link[0]; m++) 
			if (this->B[k].link[j+1] == -this->B[k].link[m+1]+SRF_STATE && (! this->solver.mode(NO_TR) || this->B[k].link[m+1] == this->solver.p[opt])) {

				bnd->zero_grid();

				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

				if  (bnd->geom) 
				for (l = 0; l < bnd->geom[0]; l++) {
					int num = bnd->geom_element(l), num_n = bnd->geom[num+1],
					  num_f = num_n+num;

					if (bnd->geom[num] == GL_LINE_STRIP) {
						for (; num < num_f-1; num++) {
							Po[0] = bnd->X[bnd->geom[num+2]];
							Po[1] = bnd->Y[bnd->geom[num+2]];
							Po[2] = 0.;
							Po[3] = bnd->X[bnd->geom[num+3]];
							Po[4] = bnd->Y[bnd->geom[num+3]];
							Po[5] = 0.;
							pp[1] = this->B[k].link[m+1];

							gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								block_phs->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
																gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
							}
							gauss_bnd->add_buffer(gauss_bnd->N);
						}
					}
				}
//////////////////////////
//...коррекция квадратуры;
				if (this->bar && this->bar->graph && this->bar->graph[0])
					block_phs->QG_curve(this->bar->ce[0]->mp);

//////////////////////////////////
//...формирование блочной матрицы;
				if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
					sprintf(msg, "block %4i: block_phs->N = %i", k, block_phs->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				if (this->solv >= ENERGY_SOLVING)
				this->TransferBB (block_phs, k, this->B[k].link[m+1], 0, E_PERIOD_COMPUT); else
				this->TransferBB (block_phs, k, this->B[k].link[m+1], 0,   PERIOD_COMPUT);

				if (! this->solver.mode(ACCUMULATION)) block_phs->add_buffer(block_phs->N);
			}
		}
	}

////////////////////////////////////////////////
//...discrete norm on stitching sides of blocks;
	if (! this->solver.mode(NO_MESSAGE) && this->volm >= BOUND_VOLUME && this->volm != BLOCK_VOLUME && 
		 ! this->solver.mode(RAPID_MODE)) Message("Discrete stitching sides...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BOUND_VOLUME && ! this->solver.mode(RAPID_MODE)) {

		LOOP_MASK(opt, k);
		if (this->volm != BLOCK_VOLUME)
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < this->B[k].bar->ce[i]->graph[1]; j++) 
			if ((m = this->B[k].link[j+1]) >= 0 && this->B[k].link[this->NUM_PHASE] == this->B[m].link[this->NUM_PHASE] &&	(! this->solver.mode(NO_TR) || m == this->solver.p[opt])) {
				bnd->zero_grid();

				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

				if (bnd->geom) {
					for (l = 0; l < bnd->geom[0]; l++) {
						int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
							num_f = num_n+num;

						if (bnd->geom[num] == GL_LINE_STRIP) {
							for (; num < num_f-1; num++) {
								Po[0] = bnd->X[bnd->geom[num+2]];
								Po[1] = bnd->Y[bnd->geom[num+2]];
								Po[2] = 0.;
								Po[3] = bnd->X[bnd->geom[num+3]];
								Po[4] = bnd->Y[bnd->geom[num+3]];
								Po[5] = 0.;

								gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[0] = gauss_bnd->get_param(0, lp);
									block_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
																	 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
								}
								gauss_bnd->add_buffer(gauss_bnd->N);
							}
						}
					}
				}
//////////////////////////////////
//...формирование блочной матрицы;
				if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
					sprintf(msg, "block %4i: block_bnd->N = %i", k, block_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				if (this->solv >= ENERGY_SOLVING)
				this->TransferBB (block_bnd, k, m, 0, E_BASIC_COMPUT); else
				this->TransferBB (block_bnd, k, m, 0,	 BASIC_COMPUT);

				if (! this->solver.mode(ACCUMULATION)) block_bnd->add_buffer(block_bnd->N);
			}
		}
	}

///////////////////////////////////////////////////////////
//...discrete norm on boundary and jump boundary condition;
	if (! this->solver.mode(NO_MESSAGE) && 
		this->volm >= BOUND_VOLUME && this->volm != BLOCK_VOLUME) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BOUND_VOLUME) {

		LOOP_MASK(opt, k);
		if (this->volm != BLOCK_VOLUME)
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < this->B[k].bar->ce[i]->graph[1]; j++) if (this->B[k].link[j+1] < 0) {

				for (mm = 1, m = this->NUM_PHASE; mm && this->NUM_PHASE && m < this->B[k].link[0]; m++) 
				if  (this->B[k].link[j+1] == -this->B[k].link[m+1]+SRF_STATE) mm = 0;
				if (mm) { 
					bnd->zero_grid();

					this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
					this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

					elem = block_plink_2D(this->B[k], l = i, m = j, id_dir, par);

/////////////////////////
//...интегрируем границу;
					if  (bnd->geom && (! this->solver.mode(NO_TR) || elem == this->solver.p[opt])) 
					for (l = 0; l < bnd->geom[0]; l++) {
						int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
							num_f = num_n+num;

						if (bnd->geom[num] == GL_LINE_STRIP) {
							for (; num < num_f-1; num++) {
								Po[0] = bnd->X[bnd->geom[num+2]];
								Po[1] = bnd->Y[bnd->geom[num+2]];
								Po[2] = 0.;
								Po[3] = bnd->X[bnd->geom[num+3]];
								Po[4] = bnd->Y[bnd->geom[num+3]];
								Po[5] = 0.;

////////////////////////////////////////////////////////////////////////////////////
//...задаем граничные условия растяжения ячейки или периодические граничные условия;
								pp[0] = pp[1] = pp[2] = pp[4] = 0.;
								if (! (this->solv%ENERGY_SOLVING)) { //...одноосное растяжение;
									if (fabs(Po[0]-par[0]) < EE_dop && fabs(Po[3]-par[0]) < EE_dop ||
										 fabs(Po[0]-par[1]) < EE_dop && fabs(Po[3]-par[1]) < EE_dop)
										pp[2] = NUMS_BND;
								}
								if ((this->solv%ENERGY_SOLVING) && elem >= 0) { //...периодические условия;
									pp[0] = par[1]-par[0];
									pp[1] = par[3]-par[2];
									pp[2] = id_dir; //...нумерация границ в стандарте (X:1-2, Y:3-4);
									pp[4] = elem;
								}
								gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
								for (int lp = 0; lp < gauss_bnd->N; lp++) {
									pp[3] = gauss_bnd->get_param(0, lp);
									bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
																	 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
								}
								gauss_bnd->add_buffer(gauss_bnd->N);
							}
						}
					}
//////////////////////////
//...коррекция квадратуры;
					if ((this->solv%ENERGY_SOLVING) && ! id_dir && this->bar && this->bar->graph && this->bar->graph[0]) //...коррекция квадратур;
					bound_bnd->QG_curve(this->bar->ce[0]->mp);
				}
//////////////////////////////////
//...формирование блочной матрицы;
				if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				if (this->solv >= ENERGY_SOLVING)
				this->GramAll(bound_bnd, k, 1, (this->solv%ENERGY_SOLVING) ? E_PERIOD_COMPUT : E_BASIC_COMPUT); else
				this->GramAll(bound_bnd, k, 1, (this->solv%ENERGY_SOLVING) ?   PERIOD_COMPUT :   BASIC_COMPUT);

				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
		}
	}

/////////////////////////////////////
//...discrete norm on initial volume;
	if (! this->solver.mode(NO_MESSAGE) && 
		this->volm >= BLOCK_VOLUME) Message("Discrete initial volume...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->volm >= BLOCK_VOLUME) {
		
		LOOP_MASK(opt, k);
		if (! this->solver.mode(RAPID_MODE)) {
			for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i)
			if (! this->solver.mode(NO_TR) || elem == this->solver.p[opt]) {

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
				gauss_bnd->zero_grid();

///////////////////////////////////////////////////////////////////////////////////
//...определяем положение криволинейной границы в блоке и строим квадратуру Гаусса;
				if (this->B[k].link[0] > this->NUM_PHASE) {
					for (j = id_first = 0; ! id_first && j < this->B[k].bar->ce[i]->graph[1]; j++) 
					for (m = this->NUM_PHASE;    ! id_first && m < this->B[k].link[0]; m++) 
					if  (this->B[k].link[j+1] == -this->B[k].link[m+1]+SRF_STATE) id_first = 1; 
					if (id_first && N_arc == 3) {
						if (j == 1) { //...переворачиваем описание контура;
							swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
							swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
						}
						if (j == 3) {
							swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
						}
						gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
						if (this->bar && this->bar->graph && this->bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_tria_curve(this->bar->ce[0]->mp, pm);
					}
					else
					if (id_first && N_arc == 4) {
						if (j == 1) { //...переворачиваем описание контура;
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
							swap(pm[0], pm[6]); swap(pm[1], pm[7]);  swap(pm[2], pm[8]);
							swap(pm[3], pm[9]); swap(pm[4], pm[10]); swap(pm[5], pm[11]);
						}
						if (j == 2) {
							swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
							swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
							swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						if (j == 3) {
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						if (j == 4) {
							swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
							swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
							swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
						if (this->bar && this->bar->graph && this->bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_quad_curve(NULL, this->bar->ce[0]->mp, pm);
					}
				}

////////////////////////////////////////////////////
//...строим квадратуру Гаусса в прямолинейном блоке;
				if (this->B[k].link[0] == this->NUM_PHASE ) {
					if (N_arc == 3)
						gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
					else
					if (N_arc == 4) {
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
					}
				}
			}
			if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
				sprintf(msg, "block %4i: gauss_bnd->N = %i", k, gauss_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
		}
		this->GramAll(gauss_bnd, k, 1, BASICtCOMPUT);
		if (! this->solver.mode(ACCUMULATION)) gauss_bnd->add_buffer(gauss_bnd->N);
	}
	delete block_phs;
	delete block_bnd;
	delete bound_bnd;
	delete gauss_bnd;
	delete bnd;
}

/////////////////////////////////////////////////////////////
//...счетная схема для одного блока и периодического условия;
template <typename T>
void CComput2D<T>::comput3(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), i, j, k, l;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(4);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[5], Po[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int  m_ilist = max(this->N/5, 1), loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);

	for (k = 0; k < this->N; k++) SkeletonBounding(this->B[k], par, NULL_STATE);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->B[k].link[this->NUM_PHASE] == -1) {

		LOOP_MASK(opt, k);
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < this->B[k].bar->ce[i]->graph[1]; j++) {//...прямоугольная ячейка;
				bnd->zero_grid();

				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

/////////////////////////
//...интегрируем границу;
				if  (bnd->geom) 
				for (l = 0; l < bnd->geom[0]; l++) {
					int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
						num_f = num_n+num;

					if (bnd->geom[num] == GL_LINE_STRIP)
					for (; num < num_f-1; num++) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = 0.;
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = 0.;

////////////////////////////////////////////
//...задаем периодические граничные условия;
						pp[0] = par[1]-par[0];
						pp[1] = par[3]-par[2];
						pp[2] = block_iddir_2D(this->B[k], i, j, par);

						gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[3] = gauss_bnd->get_param(0, lp);
							bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
															 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
				if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				this->gram2peri(bound_bnd, k, 0);

				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
		}
	}
   delete bound_bnd;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////
//...счетная схема для задачи Эшелби;
template <typename T>
void CComput2D<T>::comput4(int opt)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), i, j, k, l;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes(GRID_EL_NODES);

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(4);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);


//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[5], Po[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int  m_ilist = max(this->N/5, 1), loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);

	for (k = 0; k < this->N; k++) this->SkeletonBounding(this->B[k], par, NULL_STATE);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link && this->B[k].link[this->NUM_PHASE] == -1) {

		LOOP_MASK(opt, k);
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < this->B[k].bar->ce[i]->graph[1]; j++) {//...прямоугольная ячейка;
				bnd->zero_grid();

				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

/////////////////////////
//...интегрируем границу;
				if  (bnd->geom) 
				for (l = 0; l < bnd->geom[0]; l++) {
					int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
						num_f = num_n+num;

					if (bnd->geom[num] == GL_LINE_STRIP)
					for (; num < num_f-1; num++) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = 0.;
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = 0.;

////////////////////////////////////////////
//...задаем периодические граничные условия;
						pp[0] = par[1]-par[0];
						pp[1] = par[3]-par[2];
						pp[2] = block_iddir_2D(this->B[k], i, j, par);

						gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[3] = gauss_bnd->get_param(0, lp);
							bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
															 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
				if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
					sprintf(msg, "block %4i: bound_bnd->N = %i", k, bound_bnd->N);
					if (! this->solver.mode(NO_MESSAGE)) Message(msg);
				}
				this->gram2peri(bound_bnd, k, 0);

				if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);
			}
		}
	}
   delete bound_bnd;
	delete gauss_bnd;
   delete bnd;
}

//////////////////////////////////////////////////////////////////////////////////
//...счетная схема для интегрирования эффективных характеристик по границе блоков;
template <typename T>
void CComput2D<T>::comput5(int opt, T * K, Num_Value _FMF, int id_variant)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), id_dir, i, j, k, l, m, mm;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes(GRID_EL_NODES);

	CGrid * bound_bnd = CreateNodes();
			  bound_bnd->add_params(1);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
   double pp[1], Po[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	int  m_ilist = max(this->N/5, 1), loop;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);

	for (k = 0; k < this->N; k++) SkeletonBounding(this->B[k], par, NULL_STATE);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete boundary...");

	LOOP_OPT(loop, opt, k) if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		for (i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < this->B[k].bar->ce[i]->graph[1]; j++) {
				bnd->zero_grid();

				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->line_correct();
				this->B[k].bar->ce[this->B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);

//////////////////////////////////////////////////////////////////////////////////////////////////
//...определяем внутреннюю границу, внешнюю границу, прямоугольную ячейку (mm == 0 и id_dir == 0);
				for (mm = 1, m = this->NUM_PHASE; mm && this->NUM_PHASE && m < this->B[k].link[0]; m++) 
				if ( this->B[k].link[j+1] == -this->B[k].link[m+1]+SRF_STATE) mm = 0;
				if (! mm) m = SRF_STATE; else
				if ( this->B[k].link[j+1] < 0) block_plink_2D(this->B[k], l = i, m = j, id_dir, par);

/////////////////////////
//...интегрируем границу;
				if  (bnd->geom) 
				for (l = 0; l < bnd->geom[0]; l++) {
					int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
						num_f = num_n+num;

					if (bnd->geom[num] == GL_LINE_STRIP)
					for (; num < num_f-1; num++) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = 0.;
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = 0.;

						gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[0] = gauss_bnd->get_param(0, lp);
							bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
															 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
//////////////////////////
//...коррекция квадратуры;
				if ((! mm || mm && ! id_dir && (m = this->B[k].link[j+1]) <= SRF_STATE) && this->bar && this->bar->graph && -m+SRF_STATE < this->bar->graph[0])
					bound_bnd->QG_curve(this->bar->ce[-m+SRF_STATE]->mp);

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
	}
   delete bound_bnd;
	delete gauss_bnd;
   delete bnd;
}

/////////////////////////////////////////////////////////////////////////////////
//...счетная схема для интегрирования эффективных характеристик по объему блоков;
template <typename T>
void CComput2D<T>::comput6(int opt, T * K, Num_Value _FMF, int id_variant)
{
	int  N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), id_first, j, k, l, m, id_fast = ZERO_STATE;
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

//////////////////////////////////////////////////////////
//...initialization of axilliary arrays for discrete norm;
	double pm[3*4];
	int m_ilist = max(this->N/5, 1), loop, N_ini = 4;

	LOOP_OPT(loop, opt, k) 
		for (int j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (int j = 0; j < this->solver.JR[k][0]; j++)
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);

/////////////////////////////////////
//...discrete norm on initial volume;
	if (! this->solver.mode(NO_MESSAGE)) Message("Discrete initial volume...");

	LOOP_OPT(loop, opt, k) 
	if (this->B[k].bar && this->B[k].link) {

		LOOP_MASK(opt, k);
		for (int i = 0; i < this->B[k].bar->graph[0]; i++) IF_ANY_FACET(this->B[k].bar, i)
		if (! this->solver.mode(NO_TR)) {

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

///////////////////////////////////////////////////////////////////////////////////
//...определяем положение криволинейной границы в блоке и строим квадратуру Гаусса;
			if (this->B[k].link[0] > this->NUM_PHASE) {
				for (j = id_first = 0; ! id_first && j < this->B[k].bar->ce[i]->graph[1]; j++) 
				for (m = this->NUM_PHASE;    ! id_first && m < this->B[k].link[0]; m++) 
				if  (this->B[k].link[j+1] == -this->B[k].link[m+1]+SRF_STATE) id_first = 1; 
				if (id_first && N_arc == 3) {
					if (j == 1) { //...переворачиваем описание контура;
						swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
						swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
					}
					if (j == 3) {
						swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
						swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
					}
					gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
					if (this->bar && this->bar->graph && this->bar->graph[0]) //...коррекция квадратур;
					gauss_bnd->QG_tria_curve(this->bar->ce[0]->mp, pm);
				}
				else
				if (id_first && N_arc == 4) {
					if (j == 1) { //...переворачиваем описание контура;
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						swap(pm[0], pm[6]); swap(pm[1], pm[7]);  swap(pm[2], pm[8]);
						swap(pm[3], pm[9]); swap(pm[4], pm[10]); swap(pm[5], pm[11]);
					}
					if (j == 2) {
						swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
						swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
						swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
					}
					if (j == 3) {
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
					}
					if (j == 4) {
						swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
						swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
						swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
					}
					gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
					if (this->bar && this->bar->graph && this->bar->graph[0]) //...коррекция квадратур;
					gauss_bnd->QG_quad_curve(NULL, this->bar->ce[0]->mp, pm);
				}
			}

////////////////////////////////////////////////////
//...строим квадратуру Гаусса в прямолинейном блоке;
			if (this->B[k].link[0] == this->NUM_PHASE ) {
				if (N_arc == 3)
					gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
				else
				if (N_arc == 4) {
					swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
					gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
				}
			}
			if (! this->solver.mode(REDUCED_MESSAGE) && k%m_ilist == 0) {
				sprintf(msg, "block %4i: gauss_bnd->N = %i", k, gauss_bnd->N);
				if (! this->solver.mode(NO_MESSAGE)) Message(msg);
			}
			if (_FMF == ERR_VALUE) RigidyAll( gauss_bnd, k, K, VOLUME_COMPUT); else
			for (l = 0; l <  gauss_bnd->N; l++)
			GetFuncAllValues(gauss_bnd->X[l], gauss_bnd->Y[l], gauss_bnd->Z[l], K, k, _FMF, id_variant);

			if (! this->solver.mode(ACCUMULATION)) gauss_bnd->add_buffer(gauss_bnd->N);
		}
	}
	delete gauss_bnd;
}

///////////////////////////////////////////////////////
//...счетная схема для дополнительного перебора блоков;
template <typename T>
void CComput2D<T>::comput7(int opt, T * K)
{
	int loop, k, j;

	LOOP_OPT(loop, opt, k) 
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->solver.struct_init(this->solver.JR[k][j+this->solver.JR_SHIFT], NULL_STATE);

	LOOP_OPT(loop, opt, k) //...активизируем геометрию блочной строки;
		for (j = 0; j < this->solver.JR[k][0]; j++)
			this->bar_activate(this->B[this->solver.JR[k][j+this->solver.JR_SHIFT]], NULL_STATE, NULL_STATE);

///////////////////////////////
//...discrete norm on boundary;
	if (! this->solver.mode(NO_MESSAGE)) Message("Additional blocks search...");

	LOOP_OPT(loop, opt, k) if (this->B[k].mp) {

		LOOP_MASK(opt,  k);
		RigidyAll(NULL, k, K, COVERING_COMPUT);
	}
}

/////////////////////////////////////////
//...счетная схема для блочной структуры;
template <typename T>
Num_State CComput2D<T>::comput_kernel1(Num_Comput Num)
{
	if (this->computing_header(Num) != OK_STATE) 
		return ERR_STATE;

	if (! this->solver.changed(EXTERN_STATE)) {
		this->computing(-1, Num);

/////////////////////////////////////////////////////
//...формируем диагональ, блочную матрицу и печатаем;
		this->solver.diagonal();

		if (this->solver.mode(TESTI_GRAM))   this->solver.lagrangian(1., 0.); else
		if (this->solver.mode(TESTI_ENERGY)) this->solver.lagrangian(0., 1.); else
			this->solver.lagrangian(max(1., this->get_param(this->size_of_param()-1)), 1.-1./max(1., this->get_param(this->size_of_param()-1)));

		if (! this->solver.mode(REDUCED_PRINT)) 
			test_block_matrix(this->solver.n); else test_block_matrix(0);

///////////////////////////////////
//...запускаем внутренний решатель;
		if (! this->solver.mode(NO_MESSAGE)) Message("Blocked GaussJ...");

		double f_count = 0., g_count, dim_N = 2*this->solver.N;
		char buff[1000];
		int  k;

		for (k = 0; k < this->solver.N; k++) {
			if (! this->solver.solver(k)) {
				return ERR_STATE;
			}
			if ((g_count = (int)(100.*(k+1.)/(double)dim_N+.5)/5) > f_count) {
				f_count = g_count;
				sprintf(buff, "Solver : %d", (int)(f_count*5.));
				if (! this->solver.mode(NO_MESSAGE)) Message(buff);
			}
		}
		for (k = this->solver.N-1; k >= 0; k--) {
			if (! this->solver.solver(-k-1)) {
				return ERR_STATE;
			}
			if ((g_count = (int)(100.*(2*this->solver.N-k)/(double)dim_N+.5)/5) > f_count) {
				f_count = g_count;
				sprintf(buff, "Solver : %d", (int)(f_count*5.));
				if (! this->solver.mode(NO_MESSAGE)) Message(buff);
			}
		}
		this->solver.reset_struct();

//////////////////////////////////////////
//...переносим результаты в функции формы;
		this->shapes_init(OK_STATE);
		if (this->solver.mode(REDUCED_PRINT) || this->solver.mode(PRINT_MODE)) { //...тестовая печать;
		if (this->solver.mode(FULLY_MODE)  && ! this->solver.mode(REDUCED_PRINT)) 
			for (int i = 0; i < this->N; i++)
				this->B[this->solver.p[i]].shape->TestShape("res", this->solver.n+1, this->solver.p[i], EE*1000000.);
		else 
			for (int i = 0; i < min(this->N, 5); i++)
				this->B[this->solver.p[i]].shape->TestShape("res", this->solver.n+1, this->solver.p[i], EE*1000000.);
		}
	}
	return OK_STATE;
}

/////////////////////////////////////////////////////////////////////////////
//...счетная схема для блочной структуры с параметром продолжения по времени;
template <typename T>
Num_State CComput2D<T>::comput_kernel2(Num_Comput Num)
{
	if (this->computing_header(Num) != OK_STATE) 
		return ERR_STATE;

	if (! this->solver.changed(EXTERN_STATE)) {
//////////////////////////////////////////////////////
//...готовим тестовую печать помежуточных результатов;
		if (this->solver.mode(TESTI_MODE) && ! this->solver.mode(NO_MESSAGE)) Message("Preparing visualization...");
		
		CGrid * nd = CreateNodes(GRID_EL_NODES);
		if (this->solver.mode(TESTI_MODE)) {
			double par[6]; this->SetGeomBounding(par);

			this->BlockActivate(NULL_STATE, OK_STATE);

			int NX = 100, NY = 100, i, j;
			for (i = 0; i <= 2*NX; i++) nd->add_new_point_X(.5*i/NX*(par[1]-par[0]));
			for (j = 0; j <= 2*NY; j++) nd->add_new_point_Y(.5*j/NY*(par[3]-par[2]));

			nd->hit = (int *)new_struct(nd->N*nd->N1*sizeof(int));

			for (i = 0; i < nd->N;  i++)
			for (j = 0; j < nd->N1; j++) {
				int hit = -1;
				//Poly_struc_in2D (this, hit, nd->X[i], nd->Y[j]);
				nd->hit[i+j*nd->N] = hit;
			}
		}

///////////////////////////////////
//...запускаем внутренний решатель;
		if (! this->solver.mode(NO_MESSAGE)) Message("Blocked GaussJ...");
		this->solver.clean_mode(RAPID_MODE | CONTI_MODE);

/////////////////////////
//...итерации по времени;
		while (this->param[this->NUM_TIME+1]+EE < this->param[this->NUM_TIME+3]) {

/////////////////////////////////////////////////////
//...модификация интервала времени на последнем шаге;
			if ( this->param[this->NUM_TIME+1]+this->param[this->NUM_TIME+2] > this->param[this->NUM_TIME+3]) {
				this->param[this->NUM_TIME+2] = this->param[this->NUM_TIME+3]-this->param[this->NUM_TIME+1];
				this->solver.clean_mode(RAPID_MODE);
			}

////////////////////////////////////////////////////////////
//...формирование или пересчет правых частей и всей матрицы;
			this->solver.struct_clean(this->solver.mode(RAPID_MODE) ? NULL_STATE : 
			                    this->solver.mode(CONTI_MODE) ? ZERO_STATE : OK_STATE);
			this->computing(-1, Num);

			this->solver.lagrangian(max(1., this->get_param(this->size_of_param()-1)), 1.-1./max(1., this->get_param(this->size_of_param()-1)), 
									this->solver.mode(RAPID_MODE) ? NULL_STATE : OK_STATE);


			if (! this->solver.mode(REDUCED_PRINT)) 
				test_block_matrix(this->solver.n); else test_block_matrix(0);

////////////////////////////////////////////
//...полное или частичное обращение матрицы;
			double f_count, g_count, dim_N = 2*this->solver.N;
			char buff[1000];
			int  k;
			for (f_count = 0., k = 0; k < this->solver.N; k++) {
				if (! this->solver.solver(k)) {
					return ERR_STATE;
				}
				if ((g_count = (int)(100.*(k+1.)/(double)dim_N+.5)/5) > f_count) {
					f_count = g_count;
					sprintf(buff, "Solver : %d", (int)(f_count*5.));
					if (! this->solver.mode(NO_MESSAGE)) Message(buff);
				}
			}
			for (k = this->solver.N-1; k >= 0; k--) {
				if (! this->solver.solver(-k-1)) {
					return ERR_STATE;
				}
				if ((g_count = (int)(100.*(2*this->solver.N-k)/(double)dim_N+.5)/5) > f_count) {
					f_count = g_count;
					sprintf(buff, "Solver : %d", (int)(f_count*5.));
					if (! this->solver.mode(NO_MESSAGE)) Message(buff);
				}
			}
			this->param[this->NUM_TIME+1] += this->param[this->NUM_TIME+2]; //...интегральный момент времени, соответствующий концу интервала;

/////////////////////////////////////////
//...переносим результат в функции формы;
			this->shapes_init(OK_STATE);

///////////////////////////////////////////////
//...тестовая печать промежуточных результатов;
			if (this->solver.mode(TESTI_MODE)) {
				int kkk;
				if (! this->solver.mode(CONTI_MODE)) kkk = system("del *.grd");

				sprintf(buff, "rr_%g", f_count = ((int)(this->param[this->NUM_TIME+1]*10000.+.5))*.0001);
				this->GetSurferFormat(buff, nd, TESTI_VALUE);

				sprintf(buff, "test_%g", f_count);
				this->GetSurferFormat(buff, nd, TESTI_ANALYT_VALUE, 20);
			}
			this->solver.set_mode(RAPID_MODE | CONTI_MODE);
		}
		this->solver.reset_struct();
		this->solver.clean_mode (RAPID_MODE | CONTI_MODE);
 
		delete nd;

///////////////////////////////////
//...тестовая печать функций формы;
		if (this->solver.mode(REDUCED_PRINT) || this->solver.mode(PRINT_MODE)) for (int i = 0; i < min(this->N, 5); i++) //...тестовая печать;
			this->B[this->solver.p[i]].shape->TestShape("res", this->solver.n+1, this->solver.p[i], EE*1000000.);
	}
	return OK_STATE;
}

////////////////////////////////////
//...testing output of block matrix;
template <typename T>
void CComput2D<T>::test_block_matrix(int variant)
{
	if (this->solver.mode(PRINT_MODE)) {
		FILE * TST_gram	= fopen("GramMatrix.sta", "w"),	 //NULL,//
			  * TST_energy = fopen("EnergyMatrix.sta","w"),  //NULL,//
			  * TST_trans	= fopen("TransMatrix.sta", "w"),	 //NULL,//
			  * TST_gram_s	= NULL,//fopen("GramSymmetry.sta", "w"), //
			  * TST_ener_s = NULL,//fopen("EnergySymmetry.sta","w"),//
			  * TST_tran_s = NULL;//fopen("TransSymmetry.sta", "w");//
		if (TST_gram || TST_energy || TST_trans || TST_gram_s || TST_ener_s || TST_tran_s) {
			if (! this->solver.mode(NO_MESSAGE)) Message("Block matrix output...");
			for (int k = 0; k < this->N; k++) {
				this->solver.test_transfer_matrix(TST_trans, k);
				this->solver.test_transfer_symmetry(TST_tran_s, k);
				this->solver.test_energy_matrix(TST_energy, k, variant);
				this->solver.test_energy_symmetry(TST_ener_s, k);
				this->solver.test_gram_matrix(TST_gram, k, variant);
				this->solver.test_gram_symmetry(TST_gram_s, k);
			}
		}
		if (TST_gram)	 fclose(TST_gram);
		if (TST_energy) fclose(TST_energy);
		if (TST_trans)  fclose(TST_trans);
		if (TST_gram_s) fclose(TST_gram_s);
		if (TST_ener_s) fclose(TST_ener_s);
		if (TST_tran_s) fclose(TST_tran_s);
	}
}
#undef  Message
#endif
