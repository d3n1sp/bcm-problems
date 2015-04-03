#include "stdafx.h"
#include "cheat2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CHeat2D::NUM_GEOMT = 3;
int CHeat2D::NUM_BASIC = 5;
int CHeat2D::NUM_SHIFT = 1;
int CHeat2D::MAX_PHASE = 3;

//////////////////////////////////
//...initialization of the blocks;
int  CHeat2D::block_shape_init(Block<double> & B, int id_free)
{
	int m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == CLAYER_BLOCK) B.shape->add_shape(CreateShapeD(MP2D_CLAYER_SHAPE));
		else												  B.shape->add_shape(CreateShapeD(MP2D_POLY_SHAPE));

////////////////////////
//...setting parameters;
		if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
			B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_BASIC+NUM_SHIFT), get_param(NUM_BASIC+NUM_SHIFT*2), get_param(NUM_BASIC), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)); else 
			B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			B.shape->init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 1);

/////////////////////////////////////////////////////////////////
//...setting local system of coordinate and init parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////
//...setting parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
				B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_BASIC+NUM_SHIFT), get_param(NUM_BASIC+NUM_SHIFT*2), get_param(NUM_BASIC), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)); else 
				B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		}
		else
       if (id_free == OK_STATE) {
			for (m = 0; m < solver.id_norm; m++) {
				B.shape->set_potential(solver.hh[(int)(&B-B.B)][0][m], m);
			}
       }
		 else
		 if (id_free == NO_STATE) {
			for (m = 0; m < solver.id_norm; m++) {
				B.shape->get_potential(solver.hh[(int)(&B-B.B)][0][solver.id_norm+m], m);
			}
       }
		 else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

//////////////////////////////
//...realization of potential;
void CHeat2D::jump1(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->cpy(solver.hh[i][0][m]);
}

///////////////////////////////////
//...realization of heat intensity;
void CHeat2D::jump2(double * P, int i, int m)
{
	int    shift  = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = get_param(NUM_BASIC+shift); m += solver.id_norm;

	B[i].shape->admittance(solver.hh[i][0][m], NULL, 0., 0.);
	B[i].shape->adm_x(solver.hh[i][0][m], -P[3]*lambda);
	B[i].shape->adm_y(solver.hh[i][0][m], -P[4]*lambda);
}

///////////////////////////////////
//...realization of heat intensity;
void CHeat2D::jump2_compos(double * P, int i, int m)
{
	double lambda = 0.;
	if (B[i].shape->inverse() > 0) lambda = B[i].shape->get_param(2); else  
	if (! B[i].shape->inverse())   lambda = B[i].shape->get_param(0);	else lambda = B[i].shape->get_param(1);
	m += solver.id_norm;

	B[i].shape->admittance(solver.hh[i][0][m], NULL, 0., 0.);
	B[i].shape->adm_x(solver.hh[i][0][m], -P[3]*lambda);
	B[i].shape->adm_y(solver.hh[i][0][m], -P[4]*lambda);
}

////////////////////////////////////////////////////////////////////////////
//...realization of common jump boundary condition for matrix and inclusion;
int CHeat2D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double hh, p4, f, P[6];
		int m = solver.id_norm;

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
//		nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			hh = nd->get_param(0, l);
			p4 = nd->get_param(2, l);	if (p4 == 0.) hh = 0.;
			f  = nd->get_param(3, l);

//p4 = MIN_HIT;
//if (nd->X[l] == 0) {
//	hh = 1.-nd->Y[l];
//}
//else
//if (nd->Y[l] == 0) {
//	hh = 1.-nd->X[l];
//}

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////
//...jump of all moments;
			B[i].shape->parametrization_grad(P);

			if (p4 == MIN_HIT || p4 == MAX_HIT) {
				jump1(P, i, 0);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump2(P, i, 0);
			}

////////////////////////////
//...composition functional;
			solver.to_equationDD(i, solver.hh[i][0][m], solver.hh[i][0][m], f);
			if (fabs(hh) > EE)
			  solver.to_equationHH(i, 0, solver.hh[i][0][m], hh*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
int CHeat2D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, f, P[6], 
				 g1 = .5, f1 = 1., g2 = .5, g0 = -1., TX, TY, hh;
      int id_isolated = 0, m = solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
//		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			f  = nd->get_param(3, l);

			for (k  = (int)nd->get_param(4, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = 0.;       P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hh = -AX; break;
					case 2: TX = -AX; hh =  AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					 first = 0; k0 = k; j0 = j;
					 memset(solver.hh[i][0][m+2], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][m+2], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P); 
				jump1(P, i, 0); if (! first) solver.admittance (i, 2, 1., 0, f); else solver.admittance (i, 2, 0., 0, 1.);
				jump2(P, i, 1); 
				solver.admittance(i, 0, g1, 1, g2); 
				solver.admittance(i, 1, g0, 0, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P); 
				jump1(P, k, 0); if (! first) solver.admittance (k, 2, 1., 0, -f); else solver.admittance (k, 2, 0., 0, -1.);
				jump2(P, k, 1); 
				solver.admittance(k, 0, g1, 1, g2); 
				solver.admittance(k, 1, g0, 0, f1); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				if (fabs(hh) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m  ],  g1*hh*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m+1], -g2*hh*f);
				}
				if (first && solver.mode(REGULARIZATION2) && (id_dir == 1 || id_dir == 2)) {//...регуляризация матрицы через граничное условие;
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
 			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
			solver.to_transferTT(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для периодической задачи (один блок);
int CHeat2D::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, f, P[6], g1 = .5, f1 = 1., g2 = .5, g0 = -1.;
      int id_isolated = 0, m = solver.id_norm, id_dir;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

///////////////////////////////////
//...inclusion data in gram matrix;
		for ( int l = 0; l < nd->N; l++) if (nd->hit[l] && 
			((id_dir = (int)nd->get_param(2, l)) == 1 || id_dir == 3)) {
			P[0] = nd->X[l]; P[3] = nd->nX[l];
			P[1] = nd->Y[l]; P[4] = nd->nY[l];
			P[2] = 0.;       P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			f  = nd->get_param(3, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

////////////////////////////////////////////////////
//...вычисляем функции для формирования функционала;
			B[i].shape->parametrization_grad(P); 
			jump1(P, i, 0); 
			jump2(P, i, 1); 
			
			B[i].shape->make_common(P);
			B[i].shape->norm_common(P+3);

			if (id_dir == 1) B[i].mp[1] -= AX; else
			if (id_dir == 3) B[i].mp[2] -= AY; 
			B[i].shape->set_local_P0(B[i].mp+1);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			B[i].shape->parametrization_grad(P);
			jump1(P, i, 2); solver.admittance (i, 0, 1., 2, -1.); solver.admittance (i, 2, 2., 0, 1.);
			jump2(P, i, 3); solver.admittance (i, 1, 1., 3, -1.);
			solver.admittance(i, 0, g1, 1, g2); 
			solver.admittance(i, 1, g0, 0, f1); 

///////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов;
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m], f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			if (id_dir == 1) {//...регуляризация матрицы и правая часть;
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
				solver.to_equationHH(i, 0,	solver.hh[i][0][m], -g1*AX*f);
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], -g2*AX*f);
			}
			if (id_dir == 1) B[i].mp[1] += AX; else
			if (id_dir == 3) B[i].mp[2] += AY; 
			B[i].shape->set_local_P0(B[i].mp+1);
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////
//...inclusion of the stitching data to the solver for all blocks;
int CHeat2D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double f, P[6], f1 = 1., g1 = .5, g2 = .5, g0 = -1.;
      int id_isolated = 0, m = solver.id_norm, j, l;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
//		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump2(P, i, 1); 
				solver.admittance(i, 0, g1, 1, g2); 
				solver.admittance(i, 1, g0, 0, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump2(P, k, 1); 
				solver.admittance(k, 0, g1, 1, g2); 
				solver.admittance(k, 1, g0, 0, f1); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for all blocks;
int CHeat2D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double f, P[6], f0 = -1., f1 = 1., f2 = .5, g1 = .5;
      int id_isolated = 0, id_flag = 1, m = solver.id_norm, j, l;
		if (id_isolated) {
			f0 = g1 = 1.;
			f1 = f2 = 0.;
			id_flag = B[i].link[NUM_PHASE] == -1;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
//		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump2(P, i, 1); 
				solver.admittance(i, 0, g1, 1, f2); 
				solver.admittance(i, 1, f0, 0, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump2(P, k, 1); 
				solver.admittance(k, 0, g1, 1, f2); 
				solver.admittance(k, 1, f0, 0, f1); 

///////////////////////////
//...composition functional;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m], f);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CHeat2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, f, P[6], TX, TY, hh;
      int m = solver.id_norm, id_dir, k, j, first = 1, k0, j0;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			f  = nd->get_param(3, l);

			for (k  = (int)nd->get_param(4, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = 0.;       P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hh = -AX; break;
					case 2: TX = -AX; hh =  AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					 first = 0; k0 = k; j0 = j;
					 memset(solver.hh[i][0][m+2], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][m+2], 0, solver.dim[k]*sizeof(double));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); if (! first) solver.admittance (i, 2, 1., 0, f); else solver.admittance (i, 2, 0., 0, 1.); 
				jump2(P, i, 1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); if (! first) solver.admittance (k, 2, 1., 0, -f); else solver.admittance (k, 2, 0., 0, -1.);
				jump2(P, k, 1); 

//////////////////////////////////////////////////////////////////
//...сшивка функций и условие скачка методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				if (fabs(hh) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],  hh*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m], -hh*f);
				}
				if (first && solver.mode(REGULARIZATION2) && (id_dir == 1 || id_dir == 2)) {//...регуляризация матрицы через граничное условие;
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}

/////////////////////////////
//...энергетиеские слагаемые;
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+1],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+1], -f);
				
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
      }
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
			solver.to_transferTT(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+2], solver.hh[k0][0][m+2], 1.);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CHeat2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double hh, p4, f, P[6];
		int m = solver.id_norm;

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
//		nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hh = nd->get_param(0, l); 
			p4 = nd->get_param(2, l); if (p4 == 0.) hh = 0.;
			f  = nd->get_param(3, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
			B[i].shape->parametrization_grad(P);
			jump1(P, i, 0); 
			jump2(P, i, 1);

////////////////////////////////////////////////////
//...граничное условие методом наименьших квадратов;
			if (p4 == MIN_HIT || p4 >= NUMS_BND && p4 < NUMS_BND+20 || p4 == MAX_HIT) {
				solver.to_equationDD(i, solver.hh[i][0][m], solver.hh[i][0][m], f);
				if (fabs(hh) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m], hh*f);
			}

//////////////////////////////
//...энергетические слагаемые;
			solver.to_equationEE(i, solver.hh[i][0][m], solver.hh[i][0][m+1], f);
			if (1. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND) && fabs(hh) > EE)
				solver.to_equationEH(i, 0, solver.hh[i][0][m], -hh*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
int CHeat2D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double f, P[6];
      int m = solver.id_norm, j, l;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump2(P, i, 1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump2(P, k, 1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

/////////////////////////////
//...энергетиеские слагаемые;
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+1],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+1], -f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////
//...интегрирование параметров потока по границе блоков;
int CHeat2D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, N_elem = UnPackInts(get_param(NUM_QUAD)), N_max = UnPackInts(get_param(NUM_QUAD), 1), 
				  m = solver.id_norm, l;
		double kk = get_param(NUM_BASIC+shift), f, P[6], TH,
				 R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			memset(solver.hh[0][0][m], 0, solver.dim[0]*sizeof(double));

			B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][m], 0);
			if (B[i].shape->inverse() > 0) kk = B[i].shape->get_param(2); else  
			if (! B[i].shape->inverse())   kk = B[i].shape->get_param(0); else kk = B[i].shape->get_param(1);
			
			K[0] -= kk*TH*P[3]*(f = nd->get_param(0, l));
			K[1] -= kk*TH*P[4]* f;
			K[2+(-B[i].link[NUM_PHASE]-1)] -= P[0]*P[3]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////
//...интегрирование параметров потока по объему блоков;
int CHeat2D::rigidy2(CGrid * nd, int i, double * K)
{
	if (nd) {
      int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm, l;
      double f, P[6], PX, PY;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = -1;
			P[1] = nd->Y[l];  P[4] = 0;
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			memset(solver.hh[i][0][m], 0, solver.dim[i]*sizeof(double));
			B[i].shape->parametrization_grad(P);
			
			if (B[i].shape->type() == HT2D_CIRCLE_SHAPE)	
				jump2_compos(P, i, 0); else jump2(P, i, 0);
			PX = B[i].shape->potential(solver.hh[i][0][m], 0);

			P[4] = -1.; P[3] = P[5] = 0.;
			B[i].shape->norm_local(P+3);

			memset(solver.hh[i][0][m], 0, solver.dim[i]*sizeof(double));

			if (B[i].shape->type() == HT2D_CIRCLE_SHAPE)	
				jump2_compos(P, i, 0); else jump2(P, i, 0);
			PY = B[i].shape->potential(solver.hh[i][0][m], 0);
			
			K[0] += PX*(f = nd->get_param(0, l));
			K[1] += PY* f;
			K[2+(-B[i].link[NUM_PHASE]-1)] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...дополнительное интегрирование параметров потока для блоков с усложненными функциями;
int CHeat2D::rigidy5(CGrid * nd, int i, double * K)
{
   int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, N_elem = UnPackInts(get_param(NUM_QUAD)), N_max = UnPackInts(get_param(NUM_QUAD), 1), 
			  m = solver.id_norm, l, k;
	double kk = get_param(NUM_BASIC+shift), K1 = get_param(NUM_BASIC+NUM_SHIFT), K2 = get_param(NUM_BASIC+NUM_SHIFT*2), 
			 R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1), f, pp[1], P[6], TH, sum1, sum2, sum3, theta = M_PI*.25;
	if (K[5] < R2) theta -= acos(K[5]/R2); //...коррекция дуги окружности для взаимопроникающих включений;
	if (B[i].shape->inverse() > 0) kk = B[i].shape->get_param(2); else  
	if (! B[i].shape->inverse())   kk = B[i].shape->get_param(0); else kk = B[i].shape->get_param(1);

//////////////////////////
//...коррекция интегралов;
	if (B[i].shape->type() == MP2D_CLAYER_SHAPE && (R1 || R2))	{
		CGrid * bnd = CreateNodes();

		CGrid * bound_bnd = CreateNodes();
				  bound_bnd->add_params(1);

		CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
				  gauss_bnd->add_params(1);

////////////////////////////////////////
//...образуем временную дугу окружности;
		CCells * ce = new(CCells);
		((CCeBasic *)ce)->get_arc(R1, -M_PI*.25, M_PI*.25);
		ce->mp[1] = B[i].mp[1];
		ce->mp[2] = B[i].mp[2];

/////////////////////////////////
//...накапливаем граничные точки;
		for (ce->mp[6] = M_PI*.25, k = 0; k < 4; k++) {
			bnd->zero_grid();
			ce->grid_cells1(bnd, 0., N_max);
			if  (bnd->geom) 
			for (l = 0; l < bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num;

				if (bnd->geom[num] == GL_LINE_STRIP)
				for (; num < num_f-1; num++) {
					P[0] = bnd->X[bnd->geom[num+2]];
					P[1] = bnd->Y[bnd->geom[num+2]];
					P[2] = 0.;
					P[3] = bnd->X[bnd->geom[num+3]];
					P[4] = bnd->Y[bnd->geom[num+3]];
					P[5] = 0.;

					gauss_bnd->facet_QG(P, N_elem, SPECIAL_STATE);
					for (int lp = 0; lp < gauss_bnd->N; lp++) {
						pp[0] = gauss_bnd->get_param(0, lp);
						bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
														 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
					}
					gauss_bnd->add_buffer(gauss_bnd->N);
				}
			}
			ce->mp[6] += M_PI_2;
		}
		bound_bnd->QG_curve(ce->mp); //...коррекция квадратуры;

////////////////////////////////
//...интегрирование температуры;
//		bound_bnd->TestGrid("nodes.bln", 0.005, 0., 0., 0., AXIS_Z, 1);
		for (sum1 = sum2 = sum3 = 0., l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = 0.;					P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			memset(solver.hh[i][0][m], 0, solver.dim[i]*sizeof(double));
			B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][m], 0);

			sum1 += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum2 += TH*P[4]* f;
			sum3 += P[0]*P[3]*f;
		}
		if (sum3 < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum1 = -sum1;
			sum2 = -sum2;
			sum3 = -sum3;
		}
		K[0] += (K1-K2)*sum1;
		K[1] += (K1-K2)*sum2;
		K[3] += sum3;
		K[4] -= sum3;
		if (! solver.mode(ACCUMULATION))  bound_bnd->add_buffer(bound_bnd->N);

//////////////////////////////////////////////////////////////////////////////////////
//...повторяем граничные точки для второй окружности (с учетом проникающих включений);
		for (ce->mp[6] = M_PI*.25, ce->mp[7] = R2, ce->mp[8] = theta, k = 0; k < 4; k++) {
			bnd->zero_grid();
			ce->grid_cells1(bnd, 0., N_max);
			if  (bnd->geom) 
			for (l = 0; l < bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num;

				if (bnd->geom[num] == GL_LINE_STRIP)
				for (; num < num_f-1; num++) {
					P[0] = bnd->X[bnd->geom[num+2]];
					P[1] = bnd->Y[bnd->geom[num+2]];
					P[2] = 0.;
					P[3] = bnd->X[bnd->geom[num+3]];
					P[4] = bnd->Y[bnd->geom[num+3]];
					P[5] = 0.;

					gauss_bnd->facet_QG(P, N_elem, SPECIAL_STATE);
					for (int lp = 0; lp < gauss_bnd->N; lp++) {
						pp[0] = gauss_bnd->get_param(0, lp);
						bound_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
														 gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
					}
					gauss_bnd->add_buffer(gauss_bnd->N);
				}
			}
			ce->mp[6] += M_PI_2;
		}
		bound_bnd->QG_curve(ce->mp); //...коррекция квадратуры;

////////////////////////////////
//...интегрирование температуры;
//		bound_bnd->TestGrid("nodes.bln", 0.005, 0., 0., 0., AXIS_Z, 0);
		for (sum1 = sum2 = sum3 = 0., l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = 0.;       P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			memset(solver.hh[i][0][m], 0, solver.dim[i]*sizeof(double));
			B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][m], 0);

			sum1 += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum2 += TH*P[4]* f;
			sum3 += P[0]*P[3]*f;
		}
		if (sum3 < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum1 = -sum1;
			sum2 = -sum2;
			sum3 = -sum3;
		}
		if (K[5] < R2) //...правим площадь промежуточного слоя в случае взаимопроникающих включений;
			sum3 += sqrt(R2*R2-K[5]*K[5])*4.*K[5];

		K[0] += (K2-kk)*sum1;
		K[1] += (K2-kk)*sum2;
		K[4] += sum3;
		K[2+shift] -= sum3;
		if (! solver.mode(ACCUMULATION))  bound_bnd->add_buffer(bound_bnd->N);

		delete ce;
		delete bound_bnd;
		delete gauss_bnd;
		delete bnd;
	}
	return(OK_STATE);
}

///////////////////////////////////////
//...auxiliary function for parameters;
void CHeat2D::set_fasa_hmg(double R1, double R2, double K3, double K1, double K2)
{
	if (size_of_param() > NUM_GEOMT+2) {
		param[NUM_GEOMT] = R1;
		param[NUM_GEOMT+1] = R2;
	}
	if (size_of_param() > NUM_BASIC+NUM_SHIFT*2) {
		param[NUM_BASIC] = K3;
		param[NUM_BASIC+NUM_SHIFT] = K1;
		param[NUM_BASIC+NUM_SHIFT*2] = K2;
	}
}

///////////////////////////////////////////////////////
//...counting header for solving electrostatic problem;
int CHeat2D::computing_header(Num_Comput Num)
{
	int N_elem = UnPackInts(get_param(NUM_QUAD)), i, j, k, elem, id_dir, n_rhs = 1;
	char msg[201];
	
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CHeat2D sample: N_sm = %d, N_mpl = %d, N_elem = %d", N, UnPackInts(get_param(NUM_MPLS)), N_elem);
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
			case  PERIOD_COMPUT: Message("Periodic Blocks (one block)...");  break;
		}
		Message(" ");
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials!!!
	solver.n += 4;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par);
		for (k = 0; k < N; k++) SkeletonBounding(B[k], par);
		for (k = 0; k < N; k++) if (B[k].link && B[k].link[NUM_PHASE] == -1) {
			i = 0; while ((elem = geom_plink_2D(B[k], i, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
			i = j = 0; while ((elem = block_plink_2D(B[k], i, j, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
		}
	}
	LinkPhase2D(MAX_PHASE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : /*NULL*/OK_STATE) || 
		 ! solver.inverse_index()) {
		return ERR_STATE;
	}
	for (k = 0; k < solver.N; k++)
		solver.set_dimension(k, freedom_block(k));

	solver.struct_init();

	if (solver.mode(FULLY_MODE)) { 
		solver.test_struct("1", 0);
		solver.test_struct("2", 1);
	}
   return OK_STATE;
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CHeat2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp) {
		int m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case HEAT_VALUE: {
////////////////////////////////
//...calculation heat potential;
				B[i].shape->parametrization(P, 1);
				jump1(P, i, 0); 

				F[0] = F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant); 
				if (solv ==	SPECIAL_SOLVING) F[0] -= X-.5;
			}	break;
			case FLUX_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P);
				jump2(P, i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				jump2(P, i, 0);

				F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}  break;
			case FLUX_COMPOS_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P);

				if (B[i].shape->type() == MP2D_CLAYER_SHAPE)	
					jump2_compos(P, i, 0); else jump2(P, i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
			
				if (B[i].shape->type() == MP2D_CLAYER_SHAPE)	
					jump2_compos(P, i, 0); else jump2(P, i, 0);

				F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}  break;
			case HEAT_ANALYT_VALUE: {//...calculation heat function for the layer structure;
				double KK = get_param(NUM_BASIC)/get_param(NUM_BASIC+NUM_SHIFT), R2 = get_param(NUM_GEOMT+1);
				if (     X-.5  < -R2) F[0] = F[1] = (X-.5+R2*(1.-KK))/(1.-2.*R2*(1.-KK)); else 
				if (fabs(X-.5) <= R2) F[0] = F[1] = KK*(X-.5)/(1.-2.*R2*(1.-KK)); else 
				if (		X-.5  >  R2) F[0] = F[1] = (X-.5-R2*(1.-KK))/(1.-2.*R2*(1.-KK));
				if (solv ==	SPECIAL_SOLVING) F[0] -= X-.5;
			}	break;
			case FLUX_ANALYT_VALUE: {//...calculation heat flux for the layer structure;
				F[0] = -get_param(NUM_BASIC);
				F[1] = 0.;
			}	break;
			case LINK_VALUE: { 
				F[0] = B[i].link[1]; 
			}	break;
			default : F[0] = i; F[1] = 0;
     }
  }
}

////////////////////////////////////////////////////
//...трехфазная модель для цилиндрических включений;
double CHeat2D::TakeEshelby_two(double ff)
{
	double K1 = get_param(NUM_BASIC+NUM_SHIFT), 
			 K3 = get_param(NUM_BASIC),
			 KH = K3*(1.+ff/((1.-ff)/2.+K3/(K1-K3)));
	return(KH);
}

////////////////////////////////////////////////////////////////
//...четырехфазная модель для цилиндрических включений со слоем;
double CHeat2D::TakeEshelby(double ff, double ff_l)
{
	double K1 = get_param(NUM_BASIC+NUM_SHIFT), 
			 K2 = get_param(NUM_BASIC+NUM_SHIFT*2),
			 K3 = get_param(NUM_BASIC), c0 = ff+ff_l, c1 = ff/c0, KH;
	KH = K3*(((1.-c0)*(1.+c1)+K2/K3*(1.+c0)*(1.-c1))+K1/K2*((1.-c0)*(1.-c1)+K2/K3*(1.+c0)*(1.+c1)))/
			  (((1.+c0)*(1.+c1)+K2/K3*(1.-c0)*(1.-c1))+K1/K2*((1.+c0)*(1.-c1)+K2/K3*(1.-c0)*(1.+c1)));
	return(KH);
}

///////////////////////////////////////////////////////////
//...обобщение самосогласованной модели на множество слоев;
double CHeat2D::TakeEshelby(int nn, double * RR, double * KK)
{
	int k, dim_N = 4+2*nn;

	double **  matrix = NULL, ** hh = NULL, KH;
	set_matrix(matrix, dim_N, dim_N);
	set_matrix(hh,  2, dim_N);

//////////////////////////////////////////////////////////////////////////////////
//...заполняем и решаем систему линейных уравнений A0...Ak, Bk...Ann+1, Bnn+1, KH;
	matrix[0][0] =  1.;
	matrix[0][1] = -1.;
	matrix[0][2] = -1./sqr(RR[0]);
//
	matrix[1][0] =  KK[0];
	matrix[1][1] = -KK[1];
	matrix[1][2] =  KK[1]/sqr(RR[0]);
//
	for (k = 1; k <= nn; k++) { 
		matrix[2*k][2*k-1] =  1.;
		matrix[2*k][2*k]	 =  1./sqr(RR[k]);
		matrix[2*k][2*k+1] = -1.;
		matrix[2*k][2*k+2] = -1./sqr(RR[k]);
//
		matrix[2*k+1][2*k-1] =  KK[k];
		matrix[2*k+1][2*k]	= -KK[k]/sqr(RR[k]);
		matrix[2*k+1][2*k+1] = -KK[k+1];
		matrix[2*k+1][2*k+2] =  KK[k+1]/sqr(RR[k]);
	}
//
	matrix[2*k][2*k-1] = 1.;
	matrix[2*k][2*k]	 = 1./sqr(RR[k]);
	hh [1][2*k] = 1.;
//
	matrix[2*k+1][2*k-1] =  KK[k];
	matrix[2*k+1][2*k]	= -KK[k]/sqr(RR[k]);
	matrix[2*k+1][2*k+1] = -1.;
//
	solver.pivot_init(dim_N);
	solver.GaussJ(matrix, dim_N);
	int i, l;
	for (               i = 0; i < dim_N; i++)
	for (hh[0][i] = 0., l = 0; l < dim_N; l++) hh[0][i] += matrix[i][l]*hh[1][l];

	solver.pivot_init(dim_N);
	KH = hh[0][2*k+1];

	delete_struct(matrix);
	delete_struct(hh);

	return(KH);
}

///////////////////////////////////////////////////////////////////////////////////////
//...функция пересчета поля температур в однородной модели модели (схема установления);
void CHeat2D::TakeStabStep(double * T, int NN, double alpha)
{
////////////////////////////////
//...эффективные характеристики;
	int nn = 12, id_num;
	double KH[] = {425.1968504,3156.089194,5960.264901,2956.40327,945.0881612,415.5209072,
						292.5088532,184.6791526,127.6476378,82.76858049,57.1875,42.67500978},//...layer;
			 //KH[] = {543.9897698,3813.188479,6237.553343,3039.792388,1046.129335,415.6090244,
				//	   298.0609153,197.4879125,142.2481426,96.09193408,68.27027027,51.82525288},//...cylinder;
			 TH[] = {4.,10.,20.,30.,80.,150.,200.,300.,400.,600.,800.,1000.}, 
			 CP[] = {0.028413559, 0.04672395,0.363605601, 1.35566986,3.924635225,46.08837141,
					   190.0657332,411.3700811,609.1330877,737.8997789,806.7217391,875.9773029}, 
			 TC[] = {4.,5.,10.,15.,20.,40.,80.,150.,250.,400.,600.,1000.}, Ro = 18998, CC, KK;

///////////////////////////////////////////////////////////////////////////
//...пересчет температуры с учетом нелинейной зависимости теплопроводности;
	for (int k = 1; k < NN; k++) {
		id_num  = 0;
		while (id_num < nn && T[k] >= TH[id_num]) id_num++;
		if (id_num <= 0) KK = KH[id_num];
		if (id_num > 0 && id_num  < nn) KK = KH[id_num-1]+(KH[id_num]-KH[id_num-1])*(T[k]-TH[id_num-1])/(TH[id_num]-TH[id_num-1]);
		if (id_num == nn) KK = KH[id_num-1];

		id_num  = 0;
		while (id_num < nn && T[k] >= TC[id_num]) id_num++;
		if (id_num <= 0) CC = CP[id_num];
		if (id_num > 0 && id_num  < nn) CC = CP[id_num-1]+(CP[id_num]-CP[id_num-1])*(T[k]-TC[id_num-1])/(TC[id_num]-TC[id_num-1]);
		if (id_num == nn) CC = CP[id_num-1];
		
		T[k] += alpha*KK/(CC*Ro)*(T[k-1]+T[k+1]-2.*T[k]);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...функция пересчета поля температур в слоистой модели (схема установления);
void CHeat2D::TakeStabStep_layer(double * T, int N_SC, int N_CU, int N_cells, double alpha)
{
///////////////////////////////////
//...характеристики кремния и меди;
	int nn = 12, id_num, mm, N1, N2;
	double K_SC[] = {300.,2300.,5000.,3500.,1340.,410.,260.,150.,99.,62.,42.,31.,},
			 K_CU[] = {16200.,24000.,10800.,2170.,560.,429.,413.,401.,393.,379.,366.,352.,},
			 TH[] = {4.,10.,20.,30.,80.,150.,200.,300.,400.,600.,800.,1000.,}, 
			 C_SC[] = {0.018,0.03,0.28,1.1,3.37,44.,188.,426.,648.,794.,871.,946.,},
			 C_CU[] = {0.0916,0.1482,0.8709,2.907,7.29,58.76,202.6,322.6,373.3,397.5,416.7,451.1,},
			 TC[] = {4.,5.,10.,15.,20.,40.,80.,150.,250.,400.,600.,1000.,}, 
			 Ro_SC = 23300., Ro_CU = 8960., CC, KK, 
			 KK_SC_prev, KK_SC_next, KK_CU_prev, KK_CU_next;

///////////////////////////////////////////////////////////////////////////
//...пересчет температуры с учетом нелинейной зависимости теплопроводности;
	mm = 0;
	for (int j = 0; j < N_cells; j++) {
		if (j == 0) N1 = N_SC/2; else N1 = N_SC; N2 = N_CU;
		for (int k = 1; k < N1; k++) {
			id_num  = 0;
			while (id_num < nn && T[k+mm] >= TH[id_num]) id_num++;
			if (id_num <= 0) KK = K_SC[id_num];
			if (id_num > 0 && id_num  < nn) KK = K_SC[id_num-1]+(K_SC[id_num]-K_SC[id_num-1])*(T[k+mm]-TH[id_num-1])/(TH[id_num]-TH[id_num-1]);
			if (id_num == nn) KK = K_SC[id_num-1];

			id_num  = 0;
			while (id_num < nn && T[k+mm] >= TC[id_num]) id_num++;
			if (id_num <= 0) CC = C_SC[id_num];
			if (id_num > 0 && id_num  < nn) CC = C_SC[id_num-1]+(C_SC[id_num]-C_SC[id_num-1])*(T[k+mm]-TC[id_num-1])/(TC[id_num]-TC[id_num-1]);
			if (id_num == nn) CC = C_SC[id_num-1];
			
			T[k+mm] += alpha*KK/(CC*Ro_SC)*(T[k+mm-1]+T[k+mm+1]-2.*T[k+mm]);
			if (k == 1) KK_SC_prev = KK;
			if (k == N1-1) KK_SC_next = KK;
		}
		if (j) T[mm] = (KK_CU_next*T[mm-1]+KK_SC_prev*T[mm+1])/(KK_CU_next+KK_SC_prev);
		mm += N1; 
		for (int k = 1; k < N2; k++) {
			id_num  = 0;
			while (id_num < nn && T[k+mm] >= TH[id_num]) id_num++;
			if (id_num <= 0) KK = K_CU[id_num];
			if (id_num > 0 && id_num  < nn) KK = K_CU[id_num-1]+(K_CU[id_num]-K_CU[id_num-1])*(T[k+mm]-TH[id_num-1])/(TH[id_num]-TH[id_num-1]);
			if (id_num == nn) KK = K_CU[id_num-1];

			id_num  = 0;
			while (id_num < nn && T[k+mm] >= TC[id_num]) id_num++;
			if (id_num <= 0) CC = C_CU[id_num];
			if (id_num > 0 && id_num  < nn) CC = C_CU[id_num-1]+(C_CU[id_num]-C_CU[id_num-1])*(T[k+mm]-TC[id_num-1])/(TC[id_num]-TC[id_num-1]);
			if (id_num == nn) CC = C_CU[id_num-1];
			
			T[k+mm] += alpha*KK/(CC*Ro_CU)*(T[k+mm-1]+T[k+mm+1]-2.*T[k+mm]);
			if (k == 1) KK_CU_prev = KK;
			if (k == N2-1) KK_CU_next = KK;
		}
		T[mm] = (KK_SC_next*T[mm-1]+KK_CU_prev*T[mm+1])/(KK_SC_next+KK_CU_prev);
		mm += N2; 
	}
}
#undef  Message
