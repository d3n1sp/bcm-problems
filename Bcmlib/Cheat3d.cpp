#include "stdafx.h"
#include "cheat3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CHeat3D::NUM_GEOMT = 5;
int CHeat3D::NUM_BASIC = 7;
int CHeat3D::NUM_SHIFT = 1;
int CHeat3D::MAX_PHASE = 3;

//////////////////////////////////
//...initialization of the blocks;
int  CHeat3D::block_shape_init(Block<double> & B, int id_free)
{
	int m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE), 1);
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE), 1);

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, fabs(B.mp[8]));

			if (get_param(NUM_GEOMT-1)) { //...подключаем внешнюю систему функций;
				B.shape->add_shape(CreateShapeD(SK3D_ZOOM_SHAPE, 1));

				B.shape->init1(2, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
				B.shape->set_shape(2, fabs(B.mp[8]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
		}
		else
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE), 1);
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE), 1);

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, fabs(B.mp[8]));

			if (get_param(NUM_GEOMT-1)) {
				B.shape->add_shape(CreateShapeD(SK3D_ZOOM_SHAPE, 0)); //...внутренняя система функций;
				B.shape->add_shape(CreateShapeD(SK3D_ZOOM_SHAPE, 1)); //...подключаем внешнюю систему функций;

				B.shape->init1(2, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
				B.shape->set_shape(2, fabs(B.mp[7]), sqrt(fabs(get_param(NUM_GEOMT-1))));

				B.shape->init1(3, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
				B.shape->set_shape(3, fabs(B.mp[8]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
		}
		else {
			if ((B.type & ERR_CODE) == CLAYER_BLOCK) B.shape->add_shape(CreateShapeD(MP3D_CLAYER_SHAPE), 1); else
			if ((B.type & ERR_CODE) ==   ELLI_BLOCK) B.shape->add_shape(CreateShapeD(MP3D_ELLI_SHAPE),	1); else
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->add_shape(CreateShapeS(MP3D_SPHEROID_SHAPE), 1); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE), 1);
			else													  B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE), 1);
			
			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
				B.shape->set_shape(0, fabs(B.mp[7]), get_param(NUM_BASIC+NUM_SHIFT), get_param(NUM_BASIC+NUM_SHIFT*2), get_param(NUM_BASIC), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)); else 
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)/get_param(NUM_GEOMT)); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8])); else
				B.shape->set_shape(0, fabs(B.mp[7]));

			if (get_param(NUM_GEOMT-1)) { //...внутренняя система функций;
				B.shape->add_shape(CreateShapeS(SK3D_ZOOM_SHAPE, 0));

				B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
				B.shape->set_shape(1, fabs(B.mp[7]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
		}
 
////////////////////////////////////////////////////
//...local system of coordinate and parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////
//...setting parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
				B.shape->set_shape(2, fabs(B.mp[8]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
			else
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
				B.shape->set_shape(2, fabs(B.mp[7]), sqrt(fabs(get_param(NUM_GEOMT-1))));
				B.shape->set_shape(3, fabs(B.mp[8]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
			else {
				if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
					B.shape->set_shape(0, fabs(B.mp[7]), get_param(NUM_BASIC+NUM_SHIFT), get_param(NUM_BASIC+NUM_SHIFT*2), get_param(NUM_BASIC), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)); else 
				if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), get_param(NUM_GEOMT), get_param(NUM_GEOMT+1)/get_param(NUM_GEOMT)); else
				if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8])); else
					B.shape->set_shape(0, fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[7]), sqrt(fabs(get_param(NUM_GEOMT-1))));
			}
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
void CHeat3D::jump1(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->cpy(solver.hh[i][0][m]);
}

void CHeat3D::jump1_classic(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->cpy(solver.hh[i][0][m]);
	for (l = nm; l < B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		B[i].shape->admittance(l, solver.hh[i][0][m], NULL, 0., 0.);
}

//////////////////////////////////////
//...realization of normal derivative;
void CHeat3D::jump2(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->admittance(solver.hh[i][0][m], NULL, 0., 0.);
	B[i].shape->adm_x(solver.hh[i][0][m], P[3]);
	B[i].shape->adm_y(solver.hh[i][0][m], P[4]);
	B[i].shape->adm_z(solver.hh[i][0][m], P[5]);
}

//////////////////////////////////
//...realization of cohesion term;
void CHeat3D::jump3(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link [NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = get_param(NUM_BASIC+shift); m += solver.id_norm;

	B[i].shape->cpy		 (solver.hh[i][0][m]);
	B[i].shape->admittance(solver.hh[i][0][m], NULL, lambda, 0.);
	for (l = 0; l < nm; l++) //...обнуляем классическое поле;
		B[i].shape->admittance(l, solver.hh[i][0][m], NULL, 0., 0.);
}

///////////////////////////////////
//...realization of heat intensity;
void CHeat3D::jump4(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link [NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = get_param(NUM_BASIC+shift); m += solver.id_norm;

	B[i].shape->admittance(solver.hh[i][0][m], NULL, 0., 0.);
	B[i].shape->adm_x(solver.hh[i][0][m], -P[3]*lambda);
	B[i].shape->adm_y(solver.hh[i][0][m], -P[4]*lambda);
	B[i].shape->adm_z(solver.hh[i][0][m], -P[5]*lambda);
	for (l = nm; l < B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		B[i].shape->admittance(l, solver.hh[i][0][m], NULL, 0., 0.);
}

////////////////////////////////////////////////////////
//...realization of heat intensity for clayer functions;
void CHeat3D::jump4_compos(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link [NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = 0.;	m += solver.id_norm;

	if (B[i].shape->inverse() > 0) lambda = B[i].shape->get_param(2); else  
	if (! B[i].shape->inverse())   lambda = B[i].shape->get_param(0);	else lambda = B[i].shape->get_param(1);

	B[i].shape->admittance(solver.hh[i][0][m], NULL, 0., 0.);
	B[i].shape->adm_x(solver.hh[i][0][m], -P[3]*lambda);
	B[i].shape->adm_y(solver.hh[i][0][m], -P[4]*lambda);
	B[i].shape->adm_z(solver.hh[i][0][m], -P[5]*lambda);
	for (l = nm; l < B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		B[i].shape->admittance(l, solver.hh[i][0][m], NULL, 0., 0.);
}

////////////////////////////////////////////////////////////////////////////
//...realization of common jump boundary condition for matrix and inclusion;
int CHeat3D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double hh, p4, f, P[6];
		int m = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hh = nd->get_param(0, l); 
			p4 = nd->get_param(3, l); if (p4 == 0.) hh = 0.;
			f  = nd->get_param(4, l);

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
				jump4(P, i, 0);
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
int CHeat3D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, AZ, f, P[6], 
				 g1 = .5, f1 = 1., g2 = .5, g0 = -1., TX, TY, TZ, hh;
      int id_isolated = 0, m = solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY;
				B[k].mp[3] -= TZ; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					 first = 0; k0 = k; j0 = j;
					 memset(solver.hh[i][0][m+2], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][m+2], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); if (! first) solver.admittance (i, 2, 1., 0, f); else solver.admittance (i, 2, 0., 0, 1.);
				jump4(P, i, 1); 
				solver.admittance(i, 0, g1, 1, g2); 
				solver.admittance(i, 1, g0, 0, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); if (! first) solver.admittance (k, 2, 1., 0, -f); else solver.admittance (k, 2, 0., 0, -1.);
				jump4(P, k, 1); 
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
				if (first && solver.mode(REGULARIZATION2) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
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
int CHeat3D::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, AZ, f, P[6], g1 = .5, f1 = 1., g2 = .5, g0 = -1.;
      int id_isolated = 0, m = solver.id_norm, id_dir;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

///////////////////////////////////
//...inclusion data in gram matrix;
		for ( int l = 0; l < nd->N; l++) if (nd->hit[l] && 
			((id_dir = (int)nd->get_param(3, l)) == 1 || id_dir == 3 || id_dir == 5)) {
			P[0] = nd->X[l]; P[3] = nd->nX[l];
			P[1] = nd->Y[l]; P[4] = nd->nY[l];
			P[2] = nd->Z[l]; P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

////////////////////////////////////////////////////
//...вычисляем функции для формирования функционала;
			B[i].shape->parametrization_grad(P); 
			jump1(P, i, 0); 
			jump4(P, i, 1); 
			if (get_param(NUM_GEOMT-1)) {
				jump2(P, i, 4); 
				jump3(P, i, 5); 
			}
			B[i].shape->make_common(P);
			B[i].shape->norm_common(P+3);

			if (id_dir == 1) B[i].mp[1] -= AX; else
			if (id_dir == 3) B[i].mp[2] -= AY; else
			if (id_dir == 5) B[i].mp[3] -= AZ; 
			B[i].shape->set_local_P0(B[i].mp+1);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			B[i].shape->parametrization_grad(P);
			jump1(P, i, 2); solver.admittance (i, 0, 1., 2, -1.); solver.admittance (i, 2, 2., 0, 1.);
			jump4(P, i, 3); solver.admittance (i, 1, 1., 3, -1.);
			if (get_param(NUM_GEOMT-1)) {
				jump2(P, i, 6); solver.admittance (i, 4, 1., 6, -1.);
				jump3(P, i, 7); solver.admittance (i, 5, 1., 7, -1.);
			}
			solver.admittance(i, 0, g1, 1, g2); 
			solver.admittance(i, 1, g0, 0, f1); 
			solver.admittance(i, 4, g1, 5, g2); 
			solver.admittance(i, 5, g0, 4, f1); 

///////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов;
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m], f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			if (get_param(NUM_GEOMT-1)) {
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}
			if (id_dir == 5) {//...регуляризация матрицы и правая часть;
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
				solver.to_equationHH(i, 0,	solver.hh[i][0][m], -g1*AZ*f);
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], -g2*AZ*f);
			}
			if (id_dir == 1) B[i].mp[1] += AX; else
			if (id_dir == 3) B[i].mp[2] += AY; else
			if (id_dir == 5) B[i].mp[3] += AZ; 
			B[i].shape->set_local_P0(B[i].mp+1);
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////
//...inclusion of the stitching data to the solver for all blocks;
int CHeat3D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double f, P[6], f0 = 1., g1 = .5, g2 = .5, g0 = -1.;
      int id_isolated = 0, m = solver.id_norm, j, l;
		if (id_isolated) {
			g0 = g1 = 1.;
			f0 = g2 = 0.;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
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
				jump4(P, i, 1); 
				solver.admittance(i, 0, g1, 1, g2); 
				solver.admittance(i, 1, g0, 0, f0);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				solver.admittance(k, 0, g1, 1, g2); 
				solver.admittance(k, 1, g0, 0, f0); 

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
int CHeat3D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double f, P[6], f0 = -1., f1 = 1., f2 = .5, g1 = .5;
      int id_isolated = 1, id_flag = 1, m = solver.id_norm, j, l;
		if (id_isolated) {
			f0 = g1 = 1.;
			f1 = f2 = 0.;
			id_flag = B[i].link[NUM_PHASE] == -2;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 
				if (get_param(NUM_GEOMT-1)) {
					jump2(P, i, 2); 
					jump3(P, i, 3); 
				}
				solver.admittance(i, 0, g1, 1, f2); 
				solver.admittance(i, 1, f0, 0, f1); 
				solver.admittance(i, 2, g1, 3, f2); 
				solver.admittance(i, 3, f0, 2, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				if (get_param(NUM_GEOMT-1)) {
					jump2(P, k, 2); 
					jump3(P, k, 3); 
				}
				solver.admittance(k, 0, g1, 1, f2); 
				solver.admittance(k, 1, f0, 0, f1); 
				solver.admittance(k, 2, g1, 3, f2); 
				solver.admittance(k, 3, f0, 2, f1); 

///////////////////////////
//...composition functional;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					if (get_param(NUM_GEOMT-1)) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
						solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+3], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					}
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m], f);
					if (get_param(NUM_GEOMT-1)) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
						solver.to_transferTT(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+2], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					}
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for Eselby problem;
int CHeat3D::trans_esh(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double f, P[6], f0 = -1., f1 = 1., f2 = .5, g1 = .5, hh, pp;
      int id_isolated = 0, id_flag = 1, m = solver.id_norm, j, l;
		if (id_isolated) {
			f0 = g1 = 1.;
			f1 = f2 = 0.;
			id_flag = B[i].link[NUM_PHASE] == -2;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

/////////////////////////
//...jump of all moments;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 
				solver.admittance(i, 0, g1, 1, f2); 
				solver.admittance(i, 1, f0, 0, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				solver.admittance(k, 0, g1, 1, f2); 
				solver.admittance(k, 1, f0, 0, f1); 

///////////////////////////
//...composition functional;
				double lambda = get_param(NUM_BASIC);
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
/////////////////////////////////////
//...однородное поле -- правая часть;
				if (B[i].link[NUM_PHASE] == -2) f = -f;
				hh = -nd->Z[l]/lambda; pp = nd->nZ[l]; hh = hh*g1+pp*f2; pp = pp*f0+hh*f1;
				solver.to_equationHH(i, 0, solver.hh[i][0][m],  -hh*f);
				solver.to_equationHH(k, 0, solver.hh[k][0][m+1], pp*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CHeat3D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double AX, AY, AZ, f, P[6], TX, TY, TZ, hh;
      int m = solver.id_norm, id_dir, k, j, first = 1, k0, j0;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY;
				B[k].mp[3] -= TZ; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					 first = 0; k0 = k; j0 = j;
					 memset(solver.hh[i][0][m+2], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][m+2], 0, solver.dim[k]*sizeof(double));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); if (! first) solver.admittance (i, 2, 1., 0, f); else solver.admittance (i, 2, 0., 0, 1.);
				jump4(P, i, 1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); if (! first) solver.admittance (k, 2, 1., 0, -f); else solver.admittance (k, 2, 0., 0, -1.);
				jump4(P, k, 1); 

//////////////////////////////////////////////////////////////////
//...сшивка функций и условие скачка методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				if (fabs(hh) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],  hh*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m], -hh*f);
				}
				if (first && solver.mode(REGULARIZATION2) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}

/////////////////////////////
//...энергетиеские слагаемые;
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+1],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+1], -f);
				
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
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
int CHeat3D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double hh, p4, f, P[6];
		int m = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hh = nd->get_param(0, l); 
			p4 = nd->get_param(3, l); if (p4 == 0.) hh = 0.;
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
			B[i].shape->parametrization_grad(P);
			jump1(P, i, 0); 
			jump4(P, i, 1);

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
				solver.to_equationEH(i, 0, solver.hh[i][0][m], hh*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
int CHeat3D::transfer4(CGrid * nd, int i, int k, int id_local)
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
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
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
				jump4(P, i, 1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 

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
int CHeat3D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm, l;
      double kk = get_param(NUM_BASIC+shift), f, P[6], TH;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			B[i].shape->parametrization(P, 1);
			jump1_classic(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][m], 0);
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) { 
				if (B[i].shape->inverse() > 0) kk = B[i].shape->get_param(2); else  
				if (! B[i].shape->inverse())   kk = B[i].shape->get_param(0); else kk = B[i].shape->get_param(1);
			}		
			K[0] += kk*TH*P[3]*f;
			K[1] += kk*TH*P[4]*f;
			K[2] += kk*TH*P[5]*f;
			K[3+(-B[i].link[NUM_PHASE]-1)] += P[2]*P[5]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////
//...интегрирование параметров потока по объему блоков;
int CHeat3D::rigidy2(CGrid * nd, int i, double * K)
{
	if (nd) {
      int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm, l;
      double f, P[6], PX, PY, PZ;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = -1;
			P[1] = nd->Y[l];  P[4] = 0.;
			P[2] = nd->Z[l];  P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			B[i].shape->parametrization_grad(P);
			jump4(P, i, 0);

			PX = B[i].shape->potential(solver.hh[i][0][m], 0);

			P[4] = -1.; P[3] = P[5] = 0.;
			B[i].shape->norm_local(P+3);
			jump4(P, i, 0);

			PY = B[i].shape->potential(solver.hh[i][0][m], 0);

			P[5] = -1.; P[3] = P[4] = 0.;
			B[i].shape->norm_local(P+3);
			jump4(P, i, 0);

			PZ = B[i].shape->potential(solver.hh[i][0][m], 0);
			
			K[0] += PX*f;
			K[1] += PY*f;
			K[2] += PZ*f;
			K[3+(-B[i].link[NUM_PHASE]-1)] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...дополнительное интегрирование параметров потока для блоков с усложненными функциями;
int CHeat3D::rigidy5(CGrid * nd, int i, double * K)
{
   int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, 
		N_elem = UnPackInts(get_param(NUM_QUAD)), N_max = UnPackInts(get_param(NUM_QUAD), 1), m, k, l, l0;
	double kk = get_param(NUM_BASIC+shift), K1 = get_param(NUM_BASIC+NUM_SHIFT), K2 = get_param(NUM_BASIC+NUM_SHIFT*2), 
			 R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1), f, P[6], TH, sum[4];

	if (B[i].shape->type() == MP3D_CLAYER_SHAPE && (R1 || R2))	{
		CGrid * bound_bnd = CreateNodes();
				  bound_bnd->add_params(1);

//////////////////////////////////////////
//...образуем временную поверхность сферы;
		CCells * ce = new(CCells);
		ce->cells_new(1, 2, (l = size_of_map(2, SPHERE_GENUS))+1+size_of_dop(SPH_SEGMENT));
		ce->mp[0] = (CMap)ID_MAP(2, SPHERE_GENUS);
		ce->mp[1] = B[i].mp[1];
		ce->mp[2] = B[i].mp[2];
		ce->mp[3] = B[i].mp[3];
		ce->mp[7] = R1;
		ce->mp[l] = (CMap)SPH_SEGMENT;
      ce->mp[++l] = 2.*M_PI/N_max;
      ce->mp[++l] = 0.;
      ce->mp[++l] = M_PI/N_max; l0 = l;

/////////////////////////////////
//...накапливаем граничные точки;
		for (m = 0; m < N_max; m++, ce->mp[l0-1] = ce->mp[l0], ce->mp[l0] *= (m+1.)/m)
		for (k = 0; k < N_max; k++) {
			ce->mp[6] = ce->mp[l0-2]*k;
			ce->segms_QG(bound_bnd, N_elem);
		}

////////////////////////////////
//...интегрирование температуры;
		for (memset(sum, 0, 4*sizeof(double)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			memset(solver.hh[i][0][solver.id_norm], 0, solver.dim[i]*sizeof(double));
			B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][solver.id_norm], 0);

			sum[0] += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum[1] += TH*P[4]* f;
			sum[2] += TH*P[5]* f;
			sum[3] += P[2]*P[5]*f;
		}
		if (sum[3] < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum[0] = -sum[0];
			sum[1] = -sum[1];
			sum[2] = -sum[2];
			sum[3] = -sum[3];
		}
		K[0] += (K1-K2)*sum[0];
		K[1] += (K1-K2)*sum[1];
		K[2] += (K1-K2)*sum[2];
		K[4] += sum[3];
		K[5] -= sum[3];
		if (! solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

//////////////////////////////////////////////////////
//...повторяем граничные точки для второй поверхности;
		ce->mp[7] = R2;
      ce->mp[l0-1] = 0.;
      ce->mp[l0] = M_PI/N_max;
		
		for (m = 0; m < N_max; m++, ce->mp[l0-1] = ce->mp[l0], ce->mp[l0] *= (m+1.)/m)
		for (k = 0; k < N_max; k++) {
			ce->mp[6] = ce->mp[l0-2]*k;
			bound_bnd->sphere_intrusion_QG(ce->mp, N_elem, K[6]);
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			bound_bnd->TestGrid("nodes.bln", 0.005, 10., 20., 30., AXIS_X, 1);

////////////////////////////////
//...интегрирование температуры;
		for (memset(sum, 0, 4*sizeof(double)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			memset(solver.hh[i][0][solver.id_norm], 0, solver.dim[i]*sizeof(double));
			B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = B[i].shape->potential(solver.hh[i][0][solver.id_norm], 0);

			sum[0] += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum[1] += TH*P[4]* f;
			sum[2] += TH*P[5]* f;
			sum[3] += P[2]*P[5]*f;
		}
		if (sum[3] < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum[0] = -sum[0];
			sum[1] = -sum[1];
			sum[2] = -sum[2];
			sum[3] = -sum[3];
		}
		if (K[6] < R2) //...правим площадь промежуточного слоя в случае взаимопроникающих включений;
			sum[3] += (R2*R2-K[6]*K[6])*M_PI*2.*K[6];

		K[0] += (K2-kk)*sum[0];
		K[1] += (K2-kk)*sum[1];
		K[2] += (K2-kk)*sum[2];
		K[5] += sum[3];
		K[3+(-B[i].link[NUM_PHASE]-1)] -= sum[3];
		if (! solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

		delete ce;
		delete bound_bnd;
	}
	return(OK_STATE);
}

///////////////////////////////////////
//...auxiliary function for parameters;
void CHeat3D::set_fasa_hmg(double R1, double R2, double K3, double K1, double K2, double C)
{
	if (size_of_param() > NUM_GEOMT+2) {
		param[NUM_GEOMT-1] = C;
		param[NUM_GEOMT]   = R1;
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
int CHeat3D::computing_header(Num_Comput Num)
{
	int N_elem = UnPackInts(get_param(NUM_QUAD)), i, k, elem, id_dir, n_rhs = 2;
	char msg[201];
	
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CHeat3D sample: N_sm = %d, N_mpl = %d, N_elem = %d, Q_facet = %g", N, UnPackInts(get_param(NUM_MPLS)), N_elem, get_param(NUM_QUAD+1));
		Message(msg);

		Message(" ");
		switch (Num){
			case	 BASIC_COMPUT: Message("Junction Blocks..."); break;
			case MAPPING_COMPUT: Message("Mapping  Blocks..."); break;
			case  PERIOD_COMPUT: Message("Periodic Blocks..."); break;
		}
		Message(" ");
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += 8;//<==== number of additional auxilliary arrays!!!
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
			i = 0; while ((elem =  geom_plink_3D(B[k], i, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
			i = 0; while ((elem = block_plink_3D(B[k], i, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
		}
	}
	LinkPhase3D(MAX_PHASE);

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
void CHeat3D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
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
///////////////////////////
//...calculation potential;
				B[i].shape->parametrization(P, 1);
				jump1(P, i, 0); 

				F[0] = F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant); 

//if (sqr(P[0])+sqr(P[1])+sqr(P[2]) < 0.0025) F[0] = F[1] = 0.;

				if (solv ==	SPECIAL_SOLVING) F[0] -= Z;
			}	break;
			case HEAT_ESHE_VALUE: {
///////////////////////////
//...calculation potential;
				B[i].shape->parametrization(P, 1);
				jump1(P, i, 0); 

				F[0] = F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant); 
				if (B[i].link[NUM_PHASE] == -1) {
					double lambda = get_param(NUM_BASIC);
					F[0] = (F[1] -= Z/lambda);
				}
				if (solv ==	SPECIAL_SOLVING) F[0] -= Z;
			}	break;
			case FLUX_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P);
				jump4(P, i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[2] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}  break;
			case FLUX_ESHE_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P);
				jump4(P, i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[2] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
				if (B[i].link[NUM_PHASE] == -1)
					F[2] += P[5];
			}  break;
			case FLUX_COMPOS_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P);

				if (B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				if (B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[1] = B[i].shape->potential(solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);

				if (B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[2] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}  break;
        default : F[0] = i; F[1] = F[2] = 0;
     }
  }
}

/////////////////////////////////////////////////
//...трехфазная модель для сферических включений;
double CHeat3D::TakeEshelby_two(double ff)
{
	double K1 = get_param(NUM_BASIC+NUM_SHIFT), 
			 K3 = get_param(NUM_BASIC), 
			 KH = K3*(1.+ff/((1.-ff)/3.+K3/(K1-K3)));
	return(KH);
}

/////////////////////////////////////////////////////////////
//...четырехфазная модель для сферических включений со слоем;
double CHeat3D::TakeEshelby(double ff, double ff_l)
{
	double K1 = get_param(NUM_BASIC+NUM_SHIFT), 
			 K2 = get_param(NUM_BASIC+NUM_SHIFT*2),
			 K3 = get_param(NUM_BASIC), c0 = ff+ff_l, c1 = ff/c0, KH;
	KH = K3*(((1.-c0)*(2.+c1)+K2/K3*(1.+2.*c0)*(1.-c1))+K1/K2*((1.-c0)*(1.-c1)+K2/K3*(1.+2.*c0)*(.5+c1)))/
			  (((1.+c0*.5)*(2.+c1)+K2/K3*(1.-c0)*(1.-c1))+K1/K2*((1.+c0*.5)*(1.-c1)+K2/K3*(1.-c0)*(.5+c1)));
	return(KH);
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (прямой алгоритм);
double CHeat3D::TakeEshelby_grad(double ff)
{
	double K_I = get_param(NUM_BASIC+NUM_SHIFT), K_M = get_param(NUM_BASIC),
			 C_I = get_param(NUM_GEOMT-1), C_M = C_I;
	if (C_I == 0. || C_M == 0.) { //...классическая трехфазная модель;
		return(K_M*(1.+ff/((1.-ff)/3.+K_M/(K_I-K_M))));
	}
	double RR1 = pow(ff, 1./3.), fR3 = 1./ff, 
			 kk_I = sqrt(C_I), tt_I = exp(-2.*kk_I), HH1 = ((1.+tt_I)*kk_I-(1.-tt_I))*fR3*.5, HH2 = ((1.-tt_I)*(sqr(kk_I)+2.)-2.*(1.+tt_I)*kk_I)*fR3*.5,
			 kk_M = sqrt(C_M), tt_D = exp(kk_M*(1./RR1-1.)),
			 JJP1 = (kk_M-1.)*fR3, JJP2 = (sqr(kk_M)+2.*(1.-kk_M))*fR3, JJM1 = -(kk_M+1.)*fR3, JJM2 = (sqr(kk_M)+2.*(1.+kk_M))*fR3, 
			 matr[7][8] = {
					{ 1., HH1, -1.,   -fR3, -JJP1, -JJM1, 0., 0.},
					{ 1., HH2, -1., 2.*fR3, -JJP2, -JJM2, 0., 0.},
					{ 0., K_I*HH1, 0.,0., -K_M*JJP1, -K_M*JJM1, 0., 0.},
					{ K_I, 0., -K_M, 2.*K_M*fR3, 0., 0., 0., 0.},
					{ 0., 0., 0., 0., (kk_M/RR1-1.)*tt_D, -(kk_M/RR1+1.)/tt_D, 0., 0.},
					{ 0., 0., K_M, -2.*K_M, 0., 0., -1., 0.},
					{ 0., 0., 1., 1., 0., 0., 0., 1.}, //...равенство температур дает единицу в правую часть;
	};

//////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 7, ii[7] = {0, 0, 0, 0, 0, 0, 0}, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matr[k][l0]; matr[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
			}
	}

/////////////////////////////////////////////////////////
//...распределение температурного поля в среднем сечении;
	int id_visual = 0;
	if (id_visual) {
		FILE * TST = fopen("CCohes3D_sphere_eshelby.dat", "w");
		double A = .7, X, Y = 0., Z, RR, func; int M = 400;
		for (i = 0; i <= M; i++)
		for (l = 0; l <= M; l++) {
			X = (-1.+i*2./M)*A; Z = (-1.+l*2./M)*A; RR = sqrt(X*X+Y*Y+Z*Z); func = Z;
			if (RR < RR1) func *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR)); else
			if (RR < 1.0) func *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
			fprintf(TST, "%g   %g   %g\n", Z, X, func-Z);
		}
		fclose(TST);
		if (NULL_STATE) { //...проверка сшивки в одной точке;
			X = 0.; Z = RR1; RR = sqrt(X*X+Y*Y+Z*Z); double func1 = Z, func2 = Z;
			func1 *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR));
			func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 

			X = 0.; Z = 1.; RR = sqrt(X*X+Y*Y+Z*Z); func1 = Z; func2 = Z;
			func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
		}
	}
	return(matr[6][7]);
}

//////////////////////////////////////////////////////
//...двухфазная модель для эллипсоидального включения;
double CHeat3D::TakeEllipsoidEshelby(double ff, double eps, double phi_stream, int NX, int NY)
{
	double KI = get_param(NUM_BASIC+NUM_SHIFT), 
			 KM = get_param(NUM_BASIC), KK = KM/KI, KH = 0., z0, ch_alpha, sh_alpha, h1, h2, B1, B2, A1, A2;
	if (0. < eps && eps < 1.) {
		ch_alpha = 1./sqrt(1.-eps*eps); sh_alpha = eps*ch_alpha;
		h1 =  1./ch_alpha+(h2 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.))); h2 += 1./(eps*sh_alpha);
		B1 = (1.-KK)/(KK*ch_alpha*sqr(eps-1./eps)-(1.-KK)*h1); A1 = 1.+h1*B1;
		B2 = (KK-1.)/(KK*ch_alpha*sqr(eps-1./eps)*2.+(1.-KK)*h2); A2 = 1.+h2*B2;
		z0 = ch_alpha;

//////////////////////////////////////////
//...вычисление температуры и теплопотока;
 		int id_visual = 1;
		if (id_visual) {
			CGrid * nd = CreateNodes();
			double cos_phi = 1., sin_phi = 0., R0 = 1., f0 = R0/z0, 
				cos_stream = cos(phi_stream), sin_stream = sin(phi_stream), A = 2.5, B = 2.5;

///////////////////////////////////////////////////////////
//...разбиваем плоскость сечения в пространстве параметров;
			nd->grid_box2D(NX, NY); 
			for (int i = 0; i < nd->N; i++) {
				nd->Z[i] = (nd->X[i]-.5)*A;
				nd->X[i] = (nd->Y[i]-.5)*B;
				nd->Y[i] = 0.; 
			}

////////////////////////////////////////////////////////////////////
//...распределение поля температуры и теплопотока в среднем сечении;
			FILE * TST = fopen("CHeat3D_spheroid.dat", "w");
			double func, stream1, stream2, X, Y, Z, f, ff, f1, f2, FN2;
			for (int i = 0; i < nd->N; i++) {
				Z = nd->Z[i]*cos_stream+nd->X[i]*sin_stream;
				X = nd->X[i]*cos_stream-nd->Z[i]*sin_stream;
				Y = nd->Y[i]; 
				func  = Z*cos_stream-X*sin_stream;
				stream1 = KM*cos_stream;
				stream2 = -KM*sin_stream;
				f = sqr(X)+sqr(f0);	f1 = sqr(Z);
				if ((ff = sqr(Z)+sqr(X/eps)) < sqr(R0)) {
					func = A1*Z*cos_stream-A2*X*sin_stream;
					stream1 =  KI*A1*cos_stream;            
					stream2 = -KI*A2*sin_stream;
				}
				else {
					f2 = (f1+2.*f-4.*sqr(f0))/sqr(f);
					if (f1 < 1e-5) 
						ch_alpha = (1.-f2*f*.5*(1.-f1*f2*.25))*.5;
					else
						ch_alpha = (1.-f/f1*(sqrt(1.+f1*f2)-1.))*.5;
					ch_alpha = 1./sqrt(ch_alpha); 
					sh_alpha = sqrt(sqr(ch_alpha)-1.);
					FN2 = sqr(f0*ch_alpha)-sqr(Z/ch_alpha);
					h1 = 1./ch_alpha+(h2 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.))); h2 += ch_alpha/(sqr(ch_alpha)-1.);
					func += B1*h1*Z*cos_stream-B2*h2*X*sin_stream;
					stream1 = KM*( B1*h1*cos_stream+Z/(f0*FN2)*(B1*Z/(ch_alpha*sqr(ch_alpha))*cos_stream+2.*B2*X/(ch_alpha*sqr(sh_alpha))*sin_stream));
					stream2 = KM*(-B2*h2*sin_stream+X/(f0*FN2)*(B1*Z/(ch_alpha*sqr(sh_alpha))*cos_stream+2.*B2*X*ch_alpha/sqr(sqr(sh_alpha))*sin_stream));
				}
				fprintf(TST, "%g   %g   %g   %g   %g   %g   %g\n", nd->X[i], nd->Y[i], nd->Z[i],	func, func-nd->Z[i], 
					stream1*cos_stream-stream2*sin_stream, stream1*sin_stream+stream2*cos_stream);
			}
			fclose(TST);
		}
	}
	return(KH);
}

#undef  Message
