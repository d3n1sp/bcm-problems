#include "stdafx.h"

#include "shapes.h"
#include "clame2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CLame2D::NUM_SHEAR = 4;
int CLame2D::NUM_SHIFT = 2;

//////////////////////////////////
//...initialization of the blocks;
int CLame2D::block_shape_init(Block<double> & B, Num_State id_free)
{
	int k, m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));

////////////////////////
//...setting parameters;
      B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
      B.shape->init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));

		if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for inclusion!!!
		B.shape->init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));

/////////////////////////////////////////////////////////////////
//...setting local system of coordinate and init parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////////////
//...setting cohesion parameter and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
         B.shape->set_local(B.mp+1);
			B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		}
		else
		if (id_free == OK_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k = (int)(&B-B.B)][0][m], m);
		else
		if (id_free == NO_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k = (int)(&B-B.B)][0][solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

/////////////////////////////////////////////////
//...realization of classical displacements (Ux);
void CLame2D::jump1_classic_x(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_x     (B[i].shape->deriv);
  B[i].shape->cpy       (B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[0]*alpha*G1);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*alpha*G1, 0.);

  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////
//...realization of classical displacements (Uy);
void CLame2D::jump1_classic_y(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_y     (B[i].shape->deriv);
  B[i].shape->cpy       (B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[1]*alpha*G1);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Ux);
void CLame2D::jump2_classic_x(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_x		();
  B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
  B[i].shape->deriv_X	(B[i].shape->deriv, P[3]);
  B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[4]);
  B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
  B[i].shape->cpy_x		();

  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.-2.*alpha)*G1);
  B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.-alpha)*G1);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->p_cpy, -P[4]*alpha*G1);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uy);
void CLame2D::jump2_classic_y(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_y		();
  B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
  B[i].shape->deriv_X	(B[i].shape->deriv, P[3]);
  B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[4]);
  B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
  B[i].shape->cpy_y		();

  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
  B[i].shape->adm_y     (B[i].shape->deriv, -P[3]*alpha*G1);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->p_cpy, P[3]*(1.-alpha)*G1);
  B[i].shape->adm_y     (B[i].shape->p_cpy, P[4]*(1.-2.*alpha)*G1);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////////////////////////
//...realization shear derivative of classical displacements (Ux);
void CLame2D::jump3_classic_x(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_x		();
  B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
  B[i].shape->deriv_X	(B[i].shape->deriv, -P[4]);
  B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[3]);
  B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
  B[i].shape->cpy_x		();

  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->deriv,-P[4]*(1.-2.*alpha)*G1);
  B[i].shape->adm_y     (B[i].shape->deriv, P[3]*(1.-alpha)*G1);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->p_cpy, -P[3]*alpha*G1);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////////////////////////
//...realization shear derivative of classical displacements (Uy);
void CLame2D::jump3_classic_y(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double G1    = 1./get_param(NUM_SHEAR+shift), 
         alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_y		();
  B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
  B[i].shape->deriv_X	(B[i].shape->deriv, -P[4]);
  B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[3]);
  B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
  B[i].shape->cpy_y		();

  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
  B[i].shape->adm_y     (B[i].shape->deriv, P[4]*alpha*G1);
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
  B[i].shape->adm_x     (B[i].shape->p_cpy,-P[4]*(1.-alpha)*G1);
  B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*(1.-2.*alpha)*G1);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Px);
void CLame2D::jump4_classic_x(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double alpha = .5/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_x  ();
  B[i].shape->deriv_N();
  B[i].shape->cpy_x  ();
  B[i].shape->admittance(B[i].shape->deriv, NULL, -alpha, 0.);
  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, P[0], 0.);
  B[i].shape->adm_x     (B[i].shape->deriv, P[3]);
  B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.-alpha));
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1], 0.);
  B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*(2.*alpha-1.));
  B[i].shape->adm_x     (B[i].shape->p_cpy, P[4]*(1.-alpha));

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Py);
void CLame2D::jump4_classic_y(double * P, int i, int m)
{
  int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  double alpha = .5/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
  B[i].shape->cpy_y  ();
  B[i].shape->deriv_N();
  B[i].shape->cpy_y  ();
  B[i].shape->admittance(B[i].shape->deriv, NULL, -alpha, 0.);
  B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
  B[i].shape->admittance(B[i].shape->deriv, NULL, P[0], 0.);
  B[i].shape->adm_x     (B[i].shape->deriv, P[4]*(2.*alpha-1.));
  B[i].shape->adm_y     (B[i].shape->deriv, P[3]*(1.-alpha));
  B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1], 0.);
  B[i].shape->adm_x     (B[i].shape->p_cpy, P[3]*(1.-alpha));
  B[i].shape->adm_y     (B[i].shape->p_cpy, P[4]);

  B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
  B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////////////////////////////
//...realization double shear derivative of classical displacements (Ux);
void CLame2D::jump5_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double G1    = 1./get_param(NUM_SHEAR+shift), 
			 alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->cpy_xx	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, sqr(P[4]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[4]));
	B[i].shape->cpy_xx	 ();

	B[i].shape->cpy_yy	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, sqr(P[3]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[3]));
	B[i].shape->cpy_yy	 ();

	B[i].shape->cpy_xy	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, -P[4]*P[3]*2.);
	B[i].shape->adm		 (B[i].shape->p_cpy, -P[4]*P[3]*2.);
	B[i].shape->cpy_xy	 ();

	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[0]*alpha*G1);
	B[i].shape->adm_xx    (B[i].shape->p_cpy, -sqr(P[4])*2.*alpha*G1);
	B[i].shape->adm_xy    (B[i].shape->p_cpy,  P[4]*P[3]*2.*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, -sqr(P[3])*2.*alpha*G1);
	B[i].shape->adm_xx    (B[i].shape->deriv,  P[4]*P[3]*2.*alpha*G1);

	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////////////////////////////
//...realization double shear derivative of classical displacements (Uy);
void CLame2D::jump5_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double G1    = 1./get_param(NUM_SHEAR+shift), 
			 alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->cpy_xx	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, sqr(P[4]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[4]));
	B[i].shape->cpy_xx	 ();

	B[i].shape->cpy_yy	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, sqr(P[3]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[3]));
	B[i].shape->cpy_yy	 ();

	B[i].shape->cpy_xy	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, -P[4]*P[3]*2.);
	B[i].shape->adm		 (B[i].shape->p_cpy, -P[4]*P[3]*2.);
	B[i].shape->cpy_xy	 ();

	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[1]*alpha*G1);
	B[i].shape->adm_yy    (B[i].shape->p_cpy, -sqr(P[3])*2.*alpha*G1);
	B[i].shape->adm_xy    (B[i].shape->p_cpy,  P[4]*P[3]*2.*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, -sqr(P[4])*2.*alpha*G1);
	B[i].shape->adm_yy    (B[i].shape->deriv,  P[4]*P[3]*2.*alpha*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CLame2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CLame2D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_common(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CLame2D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			hx = nd->get_param(0, l);
			hy = nd->get_param(1, l);
			p3 = nd->get_param(2, l);
			f  = nd->get_param(3, l);
			 
			if (p3 == MIN_HIT || p3 == 2.) {
			  hy  = nd->nY[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p3 == NUMS_BND) { //...специальный случай -- одноосное растяжение;
				double nju = get_param(NUM_SHEAR+1), G1 = .5/get_param(NUM_SHEAR), 
						 AAA = -nju*G1, BBB = (1.-nju)*G1;
				hx = nd->X[l]*BBB;
				hy = nd->Y[l]*AAA; 
			}
			else 
			if (p3 == 0.) hx = hy = 0.;

/////////////////////////////
//...reset auxilliary arrays;
          for (int num = m; num < solver.n; num++)
               memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

//////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_grad(P, 2);

			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				jump1_classic_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_classic_y(P, i, 1); solver.admittance(i, 1, G1);
			}
			else
			if (0. <= p3 && p3 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
			}
			else
			if (p3 == (double)(NORMS_BND-SPECIAL_BND)) {
				jump2_classic_x(P, i, 0); solver.admittance(i, 0, G1);
				jump2_classic_y(P, i, 1); solver.admittance(i, 1, G1);
			}
			jump_make_common(i, 0);

////////////////////////////
//...composition functional;
			solver.to_equationDD(i, solver.hh[i][0][m],	 solver.hh[i][0][m], f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
		
			if (p3 == MIN_HIT || p3 == (double)(NORMS_BND-SPECIAL_BND) || p3 == NUMS_BND || p3 == MAX_HIT) f *= G1;
			if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m], hx*f);
			if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
Num_State CLame2D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx, 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1,  requl = get_param(NUM_GEOMT);
      int id_isolated = 0, m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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
				TX = TY = hx = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hx = -g1*AX; break;
					case 2: TX = -AX; hx =  g1*AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+3; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					for (int num = m+4; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2); 

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
				jump2_classic_x(P, i, 2); jump2_classic_y(P, i, 3);
				jump_make_common(i, 0);	  jump_make_common(i, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 4, 1., 0, f*G1); 
					solver.admittance (i, 5, 1., 1, f*G1);
				}
				else if (requl) {
					solver.admittance (i, 4, 0., 0, G1); 
					solver.admittance (i, 5, 0., 1, G1);
				}
				else {
					solver.admittance (i, 4, 0., 1, G1); 
					solver.admittance (i, 5, 0., 0, G1);
				}
				solver.admittance(i, 0, g1, 2, g2); solver.admittance(i, 2, g0, 0, f1); 
				solver.admittance(i, 1, g1, 3, g2); solver.admittance(i, 3, g0, 1, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump2_classic_x(P, k, 2); jump2_classic_y(P, k, 3);
				jump_make_common(k, 0);	  jump_make_common(k, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 4, 1., 0, -f*G1); 
					solver.admittance (k, 5, 1., 1, -f*G1);
				}
				else if (requl) {
					solver.admittance (k, 4, 0., 0, -G1); 
					solver.admittance (k, 5, 0., 1, -G1);
				}
				else {
					solver.admittance (k, 4, 0., 1, -G1); 
					solver.admittance (k, 5, 0., 0, -G1);
				}
				solver.admittance(k, 0, g1, 2, g2); solver.admittance(k, 2, g0, 0, f1); 
				solver.admittance(k, 1, g1, 3, g2); solver.admittance(k, 3, g0, 1, f1); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],   hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+1], hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m+2], -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+3], -hx*f);
				}
				if (first && solver.mode(REGUL_BOUNDARY)) {//...регуляризация матрицы через граничное условие;
					if (id_dir == 1 || id_dir == 2) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					}
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////////////
//...inclusion of the joining data to the solver for all blocks;
Num_State CLame2D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f1 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int id_isolated = 0, m = solver.id_norm;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump2_classic_x(P, i, 2); 
				solver.admittance (i, 0, g1, 2, g2); solver.admittance(i, 2, g0, 0, f1); 

				jump1_classic_y(P, i, 1); 
				jump2_classic_y(P, i, 3); 
				solver.admittance (i, 1, g1, 3, g2); solver.admittance(i, 3, g0, 1, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); 
				jump2_classic_x(P, k, 2); 
				solver.admittance (k, 0, g1, 2, g2); solver.admittance(k, 2, g0, 0, f1); 

				jump1_classic_y(P, k, 1); 
				jump2_classic_y(P, k, 3); 
				solver.admittance (k, 1, g1, 3, g2); solver.admittance(k, 3, g0, 1, f1); 
				
				jump_make_common(k, 0);
				jump_make_common(k, 2);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for all blocks;
Num_State CLame2D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f0 = -1., f1 = 1., f2 = .5, g1 = G1*.5;
      int id_isolated = 1, id_flag = 1, m = solver.id_norm;
		if (id_isolated) {
			f1 = f2 = 0.;
			g1 = G1; f0 = 1.;	
			id_flag = B[i].link[NUM_PHASE] == -1;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump4_classic_x(P, i, 2); 
				solver.admittance (i, 0, g1, 2, f2); solver.admittance(i, 2, f0, 0, f1); 

				jump1_classic_y(P, i, 1); 
				jump4_classic_y(P, i, 3); 
				solver.admittance (i, 1, g1, 3, f2); solver.admittance(i, 3, f0, 1, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); 
				jump4_classic_x(P, k, 2); 
				solver.admittance (k, 0, g1, 2, f2); solver.admittance(k, 2, f0, 0, f1); 

				jump1_classic_y(P, k, 1); 
				jump4_classic_y(P, k, 3); 
				solver.admittance (k, 1, g1, 3, f2); solver.admittance(k, 3, f0, 1, f1); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

///////////////////////////
//...composition functional;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
					solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m],   f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m],	  f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
Num_State CLame2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx,  requl = get_param(NUM_GEOMT);
      int	 m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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
				TX = TY = hx = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hx = -G1*AX; break;
					case 2: TX = -AX; hx =  G1*AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+3; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					for (int num = m+4; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1);
				jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3);
				jump_make_common(i, 0);	  jump_make_common(i, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 4, 1., 0, f*G1); 
					solver.admittance (i, 5, 1., 1, f*G1);
				}
				else if (requl) {
					solver.admittance (i, 4, 0., 0, G1); 
					solver.admittance (i, 5, 0., 1, G1);
				}
				else {
					solver.admittance (i, 4, 0., 1, G1); 
					solver.admittance (i, 5, 0., 0, G1);
				} 
				solver.admittance(i, 0, G1); solver.admittance(i, 1, G1);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump4_classic_x(P, k, 2); jump4_classic_y(P, k, 3);
				jump_make_common(k, 0);	  jump_make_common(k, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 4, 1., 0, -f*G1); 
					solver.admittance (k, 5, 1., 1, -f*G1);
				}
				else if (requl) {
					solver.admittance (k, 4, 0., 0, -G1); 
					solver.admittance (k, 5, 0., 1, -G1);
				}
				else {
					solver.admittance (k, 4, 0., 1, -G1); 
					solver.admittance (k, 5, 0., 0, -G1);
				}
				solver.admittance(k, 0, G1); solver.admittance(k, 1, G1);

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
			
				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],   hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+1], hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m],   -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+1], -hx*f);
				}
				if (first && solver.mode(REGUL_BOUNDARY)) {//...регуляризация матрицы через граничное условие;
					if (id_dir == 1 || id_dir == 2) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					}
				}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+2], -f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+2],  f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);

				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
Num_State CLame2D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm, shift;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
				jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3); 
				jump_make_common(i, 0);	  jump_make_common(i, 2);
				solver.admittance(i, 0, G1); 
				solver.admittance(i, 1, G1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump4_classic_x(P, k, 2); jump4_classic_y(P, k, 3); 
				jump_make_common(k, 0);   jump_make_common(k, 2);
				solver.admittance(k, 0, G1); 
				solver.admittance(k, 1, G1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && B[i].link[NUM_PHASE] < -1)  f = -f;
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+2], -f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+2],  f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);

////////////////////////////////////////////////////////////////
//...добавляем поверхностную энергию адгезии (нормировка ?????);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && 0) {
					double K1, K2, d1, d2;
					K1 = get_param(NUM_SHEAR-2+(shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT));
					d1 = get_param(NUM_SHEAR-1+ shift); 
					K2 = get_param(NUM_SHEAR-2+(shift = (-B[k].link[NUM_PHASE]-1)*NUM_SHIFT));
					d2 = get_param(NUM_SHEAR-1+ shift);

					jump3_classic_x(P, k, 4); 
					jump3_classic_y(P, k, 5); 
					solver.admittance(k, 6,  0.,   4, P[3]);
					solver.admittance(k, 4, -P[4], 5, P[3]);
					solver.admittance(k, 5,  P[4], 6, 1.);

					B[k].shape->make_local(P);
					B[k].shape->norm_local(P+3);

					B[i].shape->make_common(P);
					B[i].shape->norm_common(P+3);

					jump3_classic_x(P, i, 4); 
					jump3_classic_y(P, i, 5); 
					solver.admittance (i, 6,  0.,   4, P[3]);
					solver.admittance (i, 4, -P[4], 5, P[3]);
					solver.admittance (i, 5,  P[4], 6, 1.);
					
					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4],  f*K1);
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], -f*K2);

					solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5],  f*d1);
					solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], -f*d2);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CLame2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			hx = nd->get_param(0, l);
			hy = nd->get_param(1, l);
			p3 = nd->get_param(2, l);
			f  = nd->get_param(3, l);

			if (p3 == MIN_HIT || p3 == 2.) {
			  hy  = nd->nY[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p3 == NUMS_BND) { //...специальный случай -- одноосное растяжение;
				double nju = get_param(NUM_SHEAR+1), G1 = .5/get_param(NUM_SHEAR), 
						 AAA = -nju*G1, BBB = (1.-nju)*G1;
				hx = nd->X[l]*BBB;
				hy = nd->Y[l]*AAA; 
			}
			else 
			if (p3 == 0.) hx = hy = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

//////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_grad(P, 2);

			jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
			jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3); 	
			jump_make_common(i, 0);	  jump_make_common(i, 2);
			solver.admittance(i, 0, G1); 
			solver.admittance(i, 1, G1); 

////////////////////////////////////
//...composition collocation vector;
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 4, 0., 0, 1.);
				solver.admittance(i, 5, 0., 1, 1.);
			}
			else
			if (p3 == (double)(NORMS_BND-SPECIAL_BND)) {
				jump2_classic_x(P, i, 4); solver.admittance(i, 4, G1);
				jump2_classic_y(P, i, 5); solver.admittance(i, 5, G1);
				jump_make_common(i, 4);
			}

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p3 == MIN_HIT || p3 == (double)(NORMS_BND-SPECIAL_BND) || p3 == NUMS_BND || p3 == MAX_HIT) {
				if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
				
				if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], hy*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

			if (1. <= p3 && p3 <= (double)(NUMS_BND-SPECIAL_BND)) {
				if (fabs(hx) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m], -hx*f);
				if (fabs(hy) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+1], -hy*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////
//...интегрирование НДС на заданном наборе узлов для периодической задачи;
Num_State CLame2D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6], UX, UY;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[0][0][num], 0, solver.dim[0]*sizeof(double));

/////////////////////////////////////////////////////////////////////////////
//...вычисляем функции для формирования интеграла от деформаций и напряжений;
			B[i].shape->parametrization_grad(P, 2);
			jump1_classic_x(P, i, 0);
			jump1_classic_y(P, i, 1);

			UX = B[i].shape->potential(solver.hh[i][0][m],   0);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[0] -= G0*(UX*P[3]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			K[1] -= G0*(UX*P[4]+UY*P[3])*f;
			K[2] -= G0*(UY*P[4]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			K[3] -= UX*P[3]*f;
			K[4] -=(UX*P[4]+UY*P[3])*f*.5;
			K[5] -= UY*P[4]*f;

			UX = B[i].shape->potential(solver.hh[i][0][m],   1);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[6]  -= G0*(UX*P[3]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			K[7]  -= G0*(UX*P[4]+UY*P[3])*f;
			K[8]  -= G0*(UY*P[4]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			K[9]  -= UX*P[3]*f;
			K[10] -=(UX*P[4]+UY*P[3])*f*.5;
			K[11] -= UY*P[4]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CLame2D::set_fasa_hmg(double nju1, double nju2, double nju3, double G1, double G2, double G3)
{
	if (size_of_param() > NUM_SHEAR+NUM_SHIFT*2) {
		param[NUM_SHEAR] = G1;
		param[NUM_SHEAR+NUM_SHIFT] = G2;
		param[NUM_SHEAR+NUM_SHIFT*2] = G3;
		param[NUM_SHEAR+1] = nju1;
		param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
		param[NUM_SHEAR+1+NUM_SHIFT*2] = nju3;
	}
}

//////////////////////////////////////////////////////////
//...counting header for solving plane elasticity problem;
Num_State CLame2D::computing_header(Num_Comput Num)
{
	int i, j, k, elem, id_dir, n_rhs = 2;
	char msg[201];

	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CLame2D_draft: N_sm = %d, N_mpl = %d, G1 = %g, G2 = %g", N,
				  UnPackInts(get_param(NUM_MPLS)), get_param(NUM_SHEAR), get_param(NUM_SHEAR+NUM_SHIFT));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials!!!
	solver.n += 9;//<==== number of additional auxilliary arrays!!!
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
			i = 0; while ((elem =  geom_plink_2D(B[k], i, id_dir, par)) >= 0)			  
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
void CLame2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case DISPL_VALUE: {
///////////////////////////////
//...calculation displacements;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump1_classic_y(P, i, 1); 
				jump_make_common(i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);

				if ((solv%ENERGY_SOLVING) && 0) F[solv%ENERGY_SOLVING-1] -= X;
			}  break;
			case STRESS_X_VALUE: {
/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
				jump_make_common(i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case STRESS_Y_VALUE: {
/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
				jump_make_common(i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case DILAT_VALUE: {
/////////////////////////////////
//...dilatation of displacements;
				double alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)), G = 1./get_param(NUM_SHEAR+shift);
				B[i].shape->parametrization_hess(P, 1);

				B[i].shape->cpy_x(B[i].shape->FULL(solver.hh[i][0][m+1], 0));
				B[i].shape->cpy_y(B[i].shape->FULL(solver.hh[i][0][m+1], 0, 1));

				F[0] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant)*(1.-2*alpha)*G;
        }   break;
        case ROTATION_VALUE: {
///////////////////////////////
//...rotation of displacements;
				double G = .5/get_param(NUM_SHEAR+shift);
				B[i].shape->parametrization_hess(P, 1);

				B[i].shape->adm_y(B[i].shape->FULL(solver.hh[i][0][m+1], 0), -1.);
				B[i].shape->adm_x(B[i].shape->FULL(solver.hh[i][0][m+1], 0, 1), 1.);

				F[0] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant)*G;
        }   break;
        case ENERGY_VALUE: {
//////////////////////////
//...local energy in cell;				 
			  double E1, E2, e11, e12, e22, 
						  alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)), G = 1./get_param(NUM_SHEAR+shift);
             B[i].shape->parametrization_hess(P, 1);

				 B[i].shape->cpy_x(B[i].shape->FULL(solver.hh[i][0][m+1], 0));
				 B[i].shape->cpy_y(B[i].shape->FULL(solver.hh[i][0][m+1], 0, 1));

             E1 = B[i].shape->potential(solver.hh[i][0][m+1], 0);

             e11 = (1.-2.*alpha)*B[i].shape->potential(0, solver.hh[i][0][m+1], 0);
             e22 = (1.-2.*alpha)*B[i].shape->potential(1, solver.hh[i][0][m+1], 0);
             e12 = (1.-2.*alpha)*(B[i].shape->potential(1, solver.hh[i][0][m+1]-B[i].shape->get_NN(), 0)+
											 B[i].shape->potential(0, solver.hh[i][0][m+1]+B[i].shape->get_NN(), 0))*.5;

				 B[i].shape->cpy_xx(B[i].shape->FULL(solver.hh[i][0][m+1], 0));
				 B[i].shape->cpy_yy(B[i].shape->FULL(solver.hh[i][0][m+2], 0));
				 B[i].shape->cpy_xy(B[i].shape->FULL(solver.hh[i][0][m+3], 0));

             e11 -= (X-B[i].mp[1])*alpha*B[i].shape->potential(0, solver.hh[i][0][m+1], 0)+
						  (Y-B[i].mp[2])*alpha*B[i].shape->potential(1, solver.hh[i][0][m+1]-B[i].shape->get_NN(), 0);
             e22 -= (X-B[i].mp[1])*alpha*B[i].shape->potential(0, solver.hh[i][0][m+2], 0)+
						  (Y-B[i].mp[2])*alpha*B[i].shape->potential(1, solver.hh[i][0][m+2]-B[i].shape->get_NN(), 0);
             e12 -= (X-B[i].mp[1])*alpha*B[i].shape->potential(0, solver.hh[i][0][m+3], 0)+
						  (Y-B[i].mp[2])*alpha*B[i].shape->potential(1, solver.hh[i][0][m+3]-B[i].shape->get_NN(), 0);

//////////////////
//...local energy;
				 E1 = (4.*alpha-1.)*(1.-2.*alpha)*G*sqr(E1);

				 e11 = e11*G;
				 e22 = e22*G;
				 e12 = e12*G;

				 E2 = 2./G*(sqr(e11)+sqr(e22)+sqr(e12)*2.);

             F[0] = E1+E2;
        }     break;
        case NORMAL_X_VALUE: {
////////////////////////
//...normal derivatives;
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump2_classic_x(P, i, 0); 
				jump2_classic_y(P, i, 1); 
				jump_make_common(i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case NORMAL_Y_VALUE: {
////////////////////////
//...normal derivatives;
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump2_classic_x(P, i, 0); 
				jump2_classic_y(P, i, 1); 
				jump_make_common(i, 0);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case THIN_R_VALUE: {
/////////////////////////
//...thin surface forces;             
				P[3] = F[0]; P[4] = F[1]; P[5] = F[2];
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);

				jump5_classic_x(P, i, 0); 
				jump5_classic_y(P, i, 1); 
				solver.admittance(i, 2,  0.,   0, P[3]);
				solver.admittance(i, 0, -P[4], 1, P[3]);
				solver.admittance(i, 1,  P[4], 2, 1.);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+0], id_variant);
        }   break;
        case POTENTIAL_VALUE: {
/////////////////////////
//...potential fX and fY;
              B[i].shape->parametrization(P, 2);//...calculation with grad;
              F[0]  = B[i].shape->potential(0, B[i].shape->p_cpy, id_variant);
              F[1]  = B[i].shape->potential(1, B[i].shape->p_cpy-B[i].shape->get_NN(), id_variant);
        }     break;
        default : F[0] = i; F[1] = 0.;
     }
  }
}

//////////////////////////////////////////
//...calculaion energy by domain integral;
void CLame2D::GetEnergy(double * energy, Num_Value _FMF)
{
   if (! solver.mode(NO_MESSAGE)) Message("Energy calculation...");
	int N_elem = UnPackInts(get_param(3)), id_fast = ZERO_STATE;

	CGrid * gauss_bnd = CreateNodes(GRID_EL_NODES);
			  gauss_bnd->add_params(1);

	int i, j, k, m, l, N_ini = 4, id_first;
	double pm[3*4], f, F[3] = {0., 0., 0.};

   for (k = 0; k < N; k++) if (B[k].bar && B[k].link) {
       gauss_bnd->zero_grid();

       for (i = 0; i < B[k].bar->graph[0]; i++) 
       if  (B[k].bar->ce[i]->cells_dim() == 2 && 
            B[k].bar->ce[i]->mp[size_of_map (B[k].bar->ce[i]->mp)] == FACET_CELL) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
				int N_arc, arc, prev = B[k].bar->ce[i]->graph[(N_arc = B[k].bar->ce[i]->graph[1])+1], 
					 N_min = min(N_arc, N_ini), m1, m2;
				for (j = 1; j <= N_min; j++, prev = arc) {
					  arc = B[k].bar->ce[i]->graph[j+1];
					  if (id_fast == OK_STATE) m1 = arc; else {
						  m1  = get_num(B[k].bar->ce[i]->ce[arc]->graph, 0),
						  m2  = get_num(B[k].bar->ce[i]->ce[arc]->graph, 1);
						  if (! B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);
					  }
					  pm[3*(j-1)]   = B[k].bar->ce[i]->ce[m1]->mp[1];
					  pm[3*(j-1)+1] = B[k].bar->ce[i]->ce[m1]->mp[2];
					  pm[3*(j-1)+2] = B[k].bar->ce[i]->ce[m1]->mp[3];
				}

///////////////////////////////////////////////////////////////////////////
//...определяем положение криволинейной границы и строим квадратуру Гаусса;
            if (B[k].link[0] > NUM_PHASE) { //...криволинейный блок;
					for (j = id_first = 0; ! id_first && j < B[k].bar->ce[i]->graph[1]; j++) 
					for (m = NUM_PHASE;    ! id_first && m < B[k].link[0]; m++) 
					if  (B[k].link[j+1] == -B[k].link[m+1]+SRF_STATE) id_first = 1;
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
						if (bar && bar->graph && bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_tria_curve(bar->ce[0]->mp, pm);
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
						if (bar && bar->graph && bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_quad_curve(NULL, bar->ce[0]->mp, pm);
					}
				}
            if (B[k].link[0] == NUM_PHASE ) { //...прямолинейный блок;
					if (N_arc == 3)
						gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
					else
					if (N_arc == 4) {
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
					}
				}

/////////////////////////
//...интегрируем энергию;
				for (l = 0; l < gauss_bnd->N; l++) {
					f = gauss_bnd->get_param(0, l);
					GetFuncAllValues(gauss_bnd->X[l], gauss_bnd->Y[l], gauss_bnd->Z[l], F, k, _FMF);

					energy[-B[k].link[NUM_PHASE]*3-3] += F[0]*f;
					energy[-B[k].link[NUM_PHASE]*3-2] += F[1]*f;
					energy[-B[k].link[NUM_PHASE]*3-1] += F[2]*f;
				}
		 }
	}
	delete gauss_bnd;
   if (! solver.mode(NO_MESSAGE)) Message("O'K");
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CLame2D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k][0].GetMatrix() && 
							solver.TL && solver.TL[k]) {
		double ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]-solver.JR_SHIFT].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
		int i, j, NN = solver.dim[k];

///////////////////////
//...вычисляем энергию;
		block_shape_init(B[k], NO_STATE);
		if (h && TL) { 
			double measure = 1./get_param(NUM_SHEAR);
			for (i = 0; i < NN; i++)
			for (j = 0; j < NN; j++)
				energy[-B[k].link[NUM_PHASE]*3-3] += TL[i][j]*h[i]*h[j]*measure;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление продольного эффективного модуля растяжения/сжатия в одномерной слоистой модели цилиндрической симметрии;
double CLame2D::TakeLayer_kk(int N, double * ff, double * kk)
{
	double sum = 0.;
	int m;

/////////////////////////////////////////////////////
//...вычисляем эффективный модуль (по формуле смеси);
	for (m = 0; m < N; m++) sum += ff[m]*kk[m];
	return(sum);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного поперечного модуля Ляме в одномерной слоистой модели в цилиндрической симметрии;
double CLame2D::TakeLayer_kk(int N, double * ff, double * lm, double * mu, double * nj)
{
	double ** matr = NULL, s0, s1, sum = 0.; set_matrix(matr, 2*N, 2*N+1);
	int i1, m1, m = 2*N-1, mm;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] = nj[0]; matr[1][0] =  lm[0];
		matr[m][m] = -1.; matr[m][m+1] = -lm[0]; s0 = ff[0];
		if (N > 1) {
			matr[0][1] = -nj[1]; matr[0][2] = -1.;	matr[1][1] = -lm[1]; matr[1][2] = 2.*mu[1]; matr[1][m+1] = lm[1]-lm[0]; s1 = s0+ff[1];
			matr[m][m-2] = lm[N-1]; matr[m][m-1] = -2.*mu[N-1]*(1.-ff[N-1]); matr[m][m+1] = -lm[N-1]; matr[m-1][m-2] = nj[N-1]; matr[m-1][m-1] = (1.-ff[N-1]);
			if (N > 2) {
				for (mm = 2, m1 = 1, i1 = 1; i1 < N-1; i1++, m1 += 2) {
					matr[mm][m1] = nj[i1]; matr[mm][m1+1] =  s0/s1; matr[mm][m1+2] = -nj[i1+1]; matr[mm][m1+3] = -1.; mm++;
					matr[mm][m1] = lm[i1]; matr[mm][m1+1] = -s0/s1*mu[i1]*2.; matr[mm][m1+2] = -lm[i1+1]; matr[mm][m1+3] = 2.*mu[i1+1]; matr[mm][m+1] = lm[i1+1]-lm[i1]; mm++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}

///////////////////////////////////////
//...решаем систему линейных уравнений;
	int dim_N = 2*N, * ii, i, k, l, k0, l0; ii = (int *)new_struct(2*N*sizeof(int));
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
					else if (ii[l] > 1) goto err;
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) goto err;
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

//////////////////////////////////
//...вычисляем эффективный модуль;
	sum = matr[dim_N-1][dim_N];
err:
	delete_struct(matr); delete_struct(ii);
	return(sum);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного продольного модуля сдвига в одномерной слоистой модели в цилиндрической симметрии;
double CLame2D::TakeLayer_GG(int N, double * ff, double * mu, double * nj)
{
	double ** matr = NULL, s0, s1, sum = 0.; set_matrix(matr, 2*N, 2*N+1);
	int i1, m1, m = 2*N-1, mm;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] = 1.-2.*nj[0]; matr[1][0] = mu[0];
		matr[m][m] = 1.; matr[m][m+1] = mu[0]; s0 = ff[0];
		if (N > 1) {
			matr[0][1] = -1.+2.*nj[1]; matr[0][2] = -1.;	matr[1][1] = -mu[1]; matr[1][2] = mu[1]; matr[1][m+1] = mu[0]-mu[1]; s1 = s0+ff[1];
			matr[m][m-2] = mu[N-1]; matr[m][m-1] = -mu[N-1]*(1.-ff[N-1]); matr[m][m+1] = mu[N-1]; matr[m-1][m-2] = 1.-2.*nj[N-1]; matr[m-1][m-1] = (1.-ff[N-1]);
			if (N > 2) {
				for (mm = 2, m1 = 1, i1 = 1; i1 < N-1; i1++, m1 += 2) {
					matr[mm][m1] = 1.-2.*nj[i1]; matr[mm][m1+1] = s0/s1; matr[mm][m1+2] = -1.+2.*nj[i1+1]; matr[mm][m1+3] = -1.; mm++;
					matr[mm][m1] = mu[i1]; matr[mm][m1+1] = -s0/s1*mu[i1]; matr[mm][m1+2] = -mu[i1+1]; matr[mm][m1+3] = mu[i1+1]; matr[mm][m+1] = mu[i1]-mu[i1+1]; mm++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}

///////////////////////////////////////
//...решаем систему линейных уравнений;
	int dim_N = 2*N, * ii, i, k, l, k0, l0; ii = (int *)new_struct(2*N*sizeof(int));
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
					else if (ii[l] > 1) goto err;
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) goto err;
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

//////////////////////////////////
//...вычисляем эффективный модуль;
	sum = matr[dim_N-1][dim_N];
err:
	delete_struct(matr); delete_struct(ii);
	return(sum);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективных характеристик в одномерной слоистой модели в цилиндрической симметрии;
double CLame2D::TakeLayer_kk(int N, double * ff, double * kp, double * mu)
{
	double ** matr = NULL, s0, s1, sum = 0.; set_matrix(matr, 2*N, 2*N+1);
	int i1, m1, m = 2*N-1;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] =  1.; matr[1][0] = kp[0];
		matr[m][m] = -1.; matr[m-1][m+1] = 1.; s0 = ff[0];
		if (N > 1) {
			matr[0][1]   = -1.;     matr[0][2]   = -1.;						  matr[1][1] = -kp[1]; matr[1][2] = mu[1]; s1 = s0+ff[1];
			matr[m][m-2] = kp[N-1]; matr[m][m-1] = -mu[N-1]*(1.-ff[N-1]); matr[m-1][m-2] = 1.; matr[m-1][m-1] = (1.-ff[N-1]);
			if (N > 2) {
				for ( m = 2, m1 = 1, i1 = 1; i1 < N-1; i1++, m1 += 2) {
					matr[m][m1] = 1.;     matr[m][m1+1] =  s0/s1;        matr[m][m1+2] = -1.0;       matr[m][m1+3] = -1.0;      m++;
					matr[m][m1] = kp[i1]; matr[m][m1+1] = -s0/s1*mu[i1]; matr[m][m1+2] = -kp[i1+1];  matr[m][m1+3] =  mu[i1+1]; m++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}

///////////////////////////////////////
//...решаем систему линейных уравнений;
	int dim_N = 2*N, * ii, i, k, l, k0, l0; ii = (int *)new_struct(2*N*sizeof(int));
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
					else if (ii[l] > 1) goto err;
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) goto err;
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

//////////////////////////////////
//...вычисляем эффективный модуль;
	sum = matr[dim_N-1][dim_N];
err:
	delete_struct(matr); delete_struct(ii);
	return(sum);
}

//////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для многослойной модели;
double take_system_cyl(double ** matrix, int * ii, int dim_N, double mu_H, double nu_H)
{
	int i, k, l, k0, l0, m = dim_N-1; memset(ii, 0, dim_N*sizeof(int));
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[m-3][m] *= 1.-nu_H;
	matrix[m-2][m] *= .5-nu_H;

	matrix[m-1][m-1] *= mu_H;
	matrix[m-1][m]   *= mu_H;
	matrix[m-1][m+1] *= mu_H;

	matrix[m][m-1] *= mu_H;
	matrix[m][m]   *= mu_H;
	matrix[m][m+1] *= mu_H;

///////////////////////////////////////
//...решаем систему линейных уравнений;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matrix[k][l]) >= f) f = fabs(matrix[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matrix[k0][l]; matrix[k0][l] = matrix[l0][l]; matrix[l0][l] = f; 
			}
		if (matrix[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matrix[l0][l0]; matrix[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matrix[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matrix[k][l0]; matrix[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matrix[k][l] -= matrix[l0][l]*finv;
			}
	}
	return(matrix[m][m+1]); //...чистый сдвиг;
	//return((0.5-nu_H)*matrix[m][m+1]-4.*matrix[m-1][m+1]); //...простой сдвиг;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного поперечного модуля сдвига в слоистой модели цилиндрической симметрии;
double CLame2D::TakeLayer_GG(int N, double * ff, double * kp, double * mu, double * nj, double eps, int max_iter)
{
	double ** matr = NULL, mu_M = mu[N-1], s0 = ff[0], s1; set_matrix(matr, 4*N, 4*N+1);
	int i1, m1, m = 4*N-1;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] = 1.; matr[0][1] = -4.*nj[0]*s0; matr[m-3][m-1] = -2.; matr[m-3][m] = -1.; matr[m-3][m+1] = 1.;
		matr[1][0] = 1.; matr[1][1] = -2.*(3.-2.*nj[0])*s0; matr[m-2][m-1] =  2.; matr[m-2][m] = -1.; matr[m-2][m+1] = 1.;
		matr[2][0] = mu[0]; matr[m-1][m-1] = 6.; matr[m-1][m] = 1.; matr[m-1][m+1] = 1.;
		matr[3][0] = mu[0]; matr[3][1] = -6.*mu[0]*s0; matr[m][m-1] = -6.; matr[m][m] = -0.5; matr[m][m+1] = 1.;
		if (N > 1) {
			matr[0][2] = -1.; matr[0][3] = -(1.-nj[1]); matr[0][4] = -2.; matr[0][5] = 4.*nj[1]*s0; s1 = s0+ff[1];
			matr[m-3][m-5] = 1.; matr[m-3][m-4] = (1.-ff[N-1])*(1.-nj[N-1]); matr[m-3][m-3] = 2.*sqr(1.-ff[N-1]);  matr[m-3][m-2] = -4.*nj[N-1];

			matr[1][2] = -1.; matr[1][3] = -(0.5-nj[1]); matr[1][4] = 2.; matr[1][5] = 2.*(3.-2.*nj[1])*s0;            
			matr[m-2][m-5] = 1.; matr[m-2][m-4] = (1.-ff[N-1])*(0.5-nj[N-1]); matr[m-2][m-3] = -2.*sqr(1.-ff[N-1]);  matr[m-2][m-2] = -2.*(3.-2.*nj[N-1]);

			matr[2][2]  = -mu[1]; matr[2][3] = mu[1]; matr[2][4] = 6.*mu[1]; 
			matr[m-1][m-5] = mu_M; matr[m-1][m-4] = -mu_M*(1.-ff[N-1]); matr[m-1][m-3] = -6.*mu_M*sqr(1.-ff[N-1]);

			matr[3][2]  = -mu[1]; matr[3][3] = -0.5*mu[1]; matr[3][4] = -6.*mu[1]; matr[3][5] = 6.*mu[1]*s0; 
			matr[m][m-5] = mu_M; matr[m][m-4] = 0.5*mu_M*(1.-ff[N-1]); matr[m][m-3] = 6.*mu_M*sqr(1.-ff[N-1]); matr[m][m-2] = -6.*mu_M;
			
			if (N > 2) {
				for ( m = 4, m1 = 2, i1 = 1; i1 < N-1; i1++, m1 += 4) {
					matr[m][m1] = 1.; matr[m][m1+1] = s0/s1*(1.-nj[i1]); matr[m][m1+2] = 2.*sqr(s0/s1); matr[m][m1+3] = -4.*s1*nj[i1];      
					matr[m][m1+4] = -1.; matr[m][m1+5] = -(1.-nj[i1+1]); matr[m][m1+6] = -2.; matr[m][m1+7] = 4.*nj[i1+1]*s1; m++;

					matr[m][m1] = 1.; matr[m][m1+1] = s0/s1*(0.5-nj[i1]); matr[m][m1+2] = -2.*sqr(s0/s1); matr[m][m1+3] = -2.*s1*(3.-2.*nj[i1]);      
					matr[m][m1+4] = -1.; matr[m][m1+5] = -(0.5-nj[i1+1]); matr[m][m1+6] = 2.; matr[m][m1+7] = 2.*(3.-2.*nj[i1+1])*s1; m++;

					matr[m][m1] = mu[i1]; matr[m][m1+1] = -mu[i1]*s0/s1; matr[m][m1+2] = -6.*mu[i1]*sqr(s0/s1);
					matr[m][m1+4] = -mu[i1+1]; matr[m][m1+5] = mu[i1+1]; matr[m][m1+6] = 6.*mu[i1+1]; m++;

					matr[m][m1] = mu[i1]; matr[m][m1+1] = 0.5*mu[i1]*s0/s1; matr[m][m1+2] = 6.*mu[i1]*sqr(s0/s1); matr[m][m1+3] = -6.*mu[i1]*s1;
					matr[m][m1+4] = -mu[i1+1]; matr[m][m1+5] = -0.5*mu[i1+1]; matr[m][m1+6] = -6.*mu[i1+1]; matr[m][m1+7] = 6.*mu[i1+1]*s1; m++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1, KH = TakeLayer_kk(N, ff, kp, mu), 
		** matrix = NULL; set_matrix(matrix, m = 4*N, 4*N+1); int * ii = (int *)new_struct(m*sizeof(int)), k_iter = 0, k, l; 

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	for (mu_H1 = 0., l = 0; l < N; l++) mu_H1 += mu[l]*ff[l]; mu_H1 *= 10.;
	nu_H = 0.5*(1.-mu_H1/KH);
	for (k = 0; k < m; k++)
	for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do { 
		mu_H0 /= 2.; k_iter++;
		nu_H = 0.5*(1.-mu_H0/KH);
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	}
	while(optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = 0.5*(1.-mu_H/KH);
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
		if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	}
	while(fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5; 

//////////////////////////////////
//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6 ) mu_H = -mu_H;
	delete_struct(matr); delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

//////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для многослойной модели;
double take_system_cyl_sh(double ** matrix, int * ii, int dim_N, double mu_MH, double nu_H)
{
	int i, k, l, k0, l0, m = dim_N-1; memset(ii, 0, dim_N*sizeof(int));
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[m-3][m-1] *= mu_MH;
	matrix[m-3][m]   *= mu_MH*(1.-nu_H);
	matrix[m-3][m+1] *= mu_MH;
	matrix[m-2][m-1] *= mu_MH;
	matrix[m-2][m]   *= mu_MH*(.5-nu_H);
	matrix[m-2][m+1] *= mu_MH;
///////////////////////////////////////
//...решаем систему линейных уравнений;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matrix[k][l]) >= f) f = fabs(matrix[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matrix[k0][l]; matrix[k0][l] = matrix[l0][l]; matrix[l0][l] = f; 
			}
		if (matrix[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matrix[l0][l0]; matrix[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matrix[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matrix[k][l0]; matrix[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matrix[k][l] -= matrix[l0][l]*finv;
			}
	}
	return(matrix[m][m+1]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного поперечного модуля сдвига в слоистой модели цилиндрической симметрии;
double CLame2D::TakeLayer_sh(int N, double * ff, double * kp, double * mu, double * nj, double eps, int max_iter)
{
	double ** matr = NULL, mu_M = mu[N-1], s0 = ff[0], s1; set_matrix(matr, 4*N, 4*N+1);
	int i1, m1, m = 4*N-1;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] = 1.; matr[0][1] = -4.*nj[0]*s0; matr[m-3][m-1] = -2.; matr[m-3][m] = -1.; matr[m-3][m+1] = 1.;
		matr[1][0] = 1.; matr[1][1] = -2.*(3.-2.*nj[0])*s0; matr[m-2][m-1] =  2.; matr[m-2][m] = -1.; matr[m-2][m+1] = 1.;
		matr[2][0] = 1.; matr[m-1][m-1] = 6.; matr[m-1][m] = 1.; matr[m-1][m+1] = 1.;
		matr[3][0] = 1.; matr[3][1] = -6.*s0; matr[m][m-1] = -6.; matr[m][m] = -0.5; matr[m][m+1] = 1.;
		if (N > 1) {
			matr[0][0] *= mu[1]; matr[0][1] *= mu[1];	
			matr[1][0] *= mu[1]; matr[1][1] *= mu[1]; 

			matr[0][2] = -mu[0]; matr[0][3] = -mu[0]*(1.-nj[1]); matr[0][4] = -2.*mu[0]; matr[0][5] = 4.*mu[0]*nj[1]*s0; s1 = s0+ff[1];
			matr[m-3][m-5] = 1.; matr[m-3][m-4] = (1.-ff[N-1])*(1.-nj[N-1]); matr[m-3][m-3] = 2.*sqr(1.-ff[N-1]);  matr[m-3][m-2] = -4.*nj[N-1];

			matr[1][2] = -mu[0]; matr[1][3] = -mu[0]*(0.5-nj[1]); matr[1][4] = 2.*mu[0]; matr[1][5] = 2.*mu[0]*(3.-2.*nj[1])*s0;            
			matr[m-2][m-5] = 1.; matr[m-2][m-4] = (1.-ff[N-1])*(0.5-nj[N-1]); matr[m-2][m-3] = -2.*sqr(1.-ff[N-1]);  matr[m-2][m-2] = -2.*(3.-2.*nj[N-1]);

			matr[2][2]  = -1.; matr[2][3] = 1.; matr[2][4] = 6.; 
			matr[m-1][m-5] = 1.; matr[m-1][m-4] = -(1.-ff[N-1]); matr[m-1][m-3] = -6.*sqr(1.-ff[N-1]);

			matr[3][2]  = -1.; matr[3][3] = -0.5; matr[3][4] = -6.; matr[3][5] = 6.*s0; 
			matr[m][m-5] = 1.; matr[m][m-4] = 0.5*(1.-ff[N-1]); matr[m][m-3] = 6.*sqr(1.-ff[N-1]); matr[m][m-2] = -6.;
			
			if (N > 2) {
				for ( m = 4, m1 = 2, i1 = 1; i1 < N-1; i1++, m1 += 4) {
					matr[m][m1] = mu[i1+1]; matr[m][m1+1] = s0/s1*mu[i1+1]*(1.-nj[i1]); matr[m][m1+2] = 2.*sqr(s0/s1)*mu[i1+1]; matr[m][m1+3] = -4.*s1*mu[i1+1]*nj[i1];      
					matr[m][m1+4] = -mu[i1]; matr[m][m1+5] = -mu[i1]*(1.-nj[i1+1]); matr[m][m1+6] = -2.*mu[i1]; matr[m][m1+7] = 4.*mu[i1]*nj[i1+1]*s1; m++;

					matr[m][m1] = mu[i1+1]; matr[m][m1+1] = s0/s1*mu[i1+1]*(0.5-nj[i1]); matr[m][m1+2] = -2.*sqr(s0/s1)*mu[i1+1]; matr[m][m1+3] = -2.*s1*mu[i1+1]*(3.-2.*nj[i1]);      
					matr[m][m1+4] = -mu[i1]; matr[m][m1+5] = -mu[i1]*(0.5-nj[i1+1]); matr[m][m1+6] = 2.*mu[i1]; matr[m][m1+7] = 2.*mu[i1]*(3.-2.*nj[i1+1])*s1; m++;

					matr[m][m1] = 1.; matr[m][m1+1] = -s0/s1; matr[m][m1+2] = -6.*sqr(s0/s1);
					matr[m][m1+4] = -1.; matr[m][m1+5] = 1.; matr[m][m1+6] = 6.; m++;

					matr[m][m1] = 1.; matr[m][m1+1] = 0.5*s0/s1; matr[m][m1+2] = 6.*sqr(s0/s1); matr[m][m1+3] = -6.*s1;
					matr[m][m1+4] = -1.; matr[m][m1+5] = -0.5; matr[m][m1+6] = -6.; matr[m][m1+7] = 6.*s1; m++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}
	double optim, sgn0, sgn1, nu_H, mu_H, mu_MH, mu_MH0, mu_MH1, KH = TakeLayer_kk(N, ff, kp, mu), eps0 = 1e-6, 
		** matrix = NULL; set_matrix(matrix, m = 4*N, 4*N+1); int * ii = (int *)new_struct(m*sizeof(int)), k_iter = 0, k, l;
/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.; nu_H = -1000.;
	for (k = 0; k < m; k++)
	for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl_sh(matrix, ii, m, mu_MH0, nu_H);
	if (optim < 0.) sgn0 = -1.; else sgn0 = 1.;

	mu_MH1 = mu_M/KH;
	do { 
		mu_MH1 *= 10.; k_iter++; 
		nu_H = 0.5*(1.-mu_M/(mu_MH1*KH));
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl_sh(matrix, ii, m, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H = 0.5*(1.-mu_M/(mu_MH*KH));
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl_sh(matrix, ii, m, mu_MH, nu_H);
		if (optim*sgn0 > 0.) mu_MH0 = mu_MH; else
		if (optim*sgn0 < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);
	mu_H = 2.*mu_M/(mu_MH0+mu_MH1);

//////////////////////////////////
//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > eps0 ) mu_H = -mu_H;
	delete_struct(matr); delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

//////////////////////////////////////////////////////////////////////////////////////
//...эффективный модуль Юнга вдоль слоев на основе модели симметричной слоистой среды;
double CLame2D::TakeLayer_E1(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 ku1 = (1.-nj1)/(.5-nj1)*G1, lm1 = nj1/(1.-nj1)*ku1,
			 ku2 = (1.-nj2)/(.5-nj2)*G2, lm2 = nj2/(1.-nj2)*ku2,
			 alpha = ff*(1.-ff)/(ff*ku1+(1.-ff)*ku2);
	double kk_1 = ff*ku2+(1.-ff)*ku1-alpha*sqr(lm2-lm1), kk_2 = 1./(ff/ku2+(1.-ff)/ku1),
			 lm_2 = ff*lm2+(1.-ff)*lm1-alpha*(lm2-lm1)*(ku2-ku1);
	return(4.*sqr(ff*G2+(1.-ff)*G1)*kk_2/(kk_1*kk_2-sqr(lm_2)));
}

////////////////////////////////////////////////////////////////////////////////////////
//...эффективный модуль Юнга поперек слоев на основе модели симметричной слоистой среды;
double CLame2D::TakeLayer_E2(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 ku1 = (1.-nj1)/(.5-nj1)*G1, lm1 = nj1/(1.-nj1)*ku1,
			 ku2 = (1.-nj2)/(.5-nj2)*G2, lm2 = nj2/(1.-nj2)*ku2,
			 alpha = ff*(1.-ff)/(ff*ku1+(1.-ff)*ku2);
	double kk_1 = ff*ku2+(1.-ff)*ku1-alpha*sqr(lm2-lm1), kk_2 = 1./(ff/ku2+(1.-ff)/ku1),
			 lm_1 = ff*lm2+(1.-ff)*lm1-alpha*sqr(lm2-lm1),
			 lm_2 = lm_1-2.*alpha*(lm2-lm1)*(G2-G1);
	return(kk_2-2.*sqr(lm_2)/(kk_1+lm_1));
}

////////////////////////////////////////////////////////////////////////////////////////
//...эффективный модуль сдвига вдоль слоев на основе модели симметричной слоистой среды;
double CLame2D::TakeLayer_G1(double ff)
{
	double G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT);
	return(ff*G2+(1.-ff)*G1);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...эффективный модуль сдвига поперек слоев на основе модели симметричной слоистой среды;
double CLame2D::TakeLayer_G2(double ff)
{
	double G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT);
	return(1./(ff/G2+(1.-ff)/G1));
}
#undef  Message
