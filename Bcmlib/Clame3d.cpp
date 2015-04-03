#include "stdafx.h"

#include "cgrid_el.h"
#include "clame3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}


int CLame3D::NUM_GEOMT = 4;
int CLame3D::NUM_SHEAR = 8;
int CLame3D::NUM_SHIFT = 6;
int CLame3D::MAX_PHASE = 3;
int CLame3D::NUM_HESS = 14;
extern void (* FUNC)(double * pp);

/////////////////////////////////
int forced_phase = -777;       //...принудительное значение фазы;
/////////////////////////////////

//////////////////////////////////
//...destructor of CHydro3D class;
CLame3D::~CLame3D (void)
{
	for (int k = 0; TT && TT[k]; k++) {
		delete_struct(TT[k]);
		delete_struct(TH[k]);
	}
   delete_struct(TT);
   delete_struct(TH);
   delete_struct(C0);
   delete_struct(C1);
   delete_struct(C2);
   delete_struct(A1);
   delete_struct(B1);
   delete_struct(A2);
   delete_struct(B2);
}

//////////////////////////////////
//...initialization of the blocks;
int  CLame3D::block_shape_init(Block<double> & B, int id_free)
{
	int k, m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, fabs(B.mp[8]));
		}
		else
		if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE), NULL_STATE);

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 0);
			B.shape->set_shape(1, sqr(B.mp[8] = get_param(NUM_GEOMT+1))/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));
		}
		else {
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->add_shape(CreateShapeD(MP3D_ELLI_SHAPE)); else
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->add_shape(CreateShapeS(MP3D_SPHEROID_FULL_SHAPE)); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE));
			else													  B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			
			B.shape->init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), get_param(NUM_GEOMT+1)/get_param(NUM_GEOMT), get_param(NUM_GEOMT)); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8])); else
			B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for inclusion!!!
			B.shape->init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
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
			}
			else
			if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, sqr(B.mp[8])/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));
			}
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), get_param(NUM_GEOMT+1)/get_param(NUM_GEOMT), get_param(NUM_GEOMT)); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8]));
			else B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));
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
void CLame3D::jump1_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_x     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[0]*alpha);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[1]*alpha, 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Ux) for sphere inclusion (2body conditions);
void CLame3D::jump1_sphere_x_2body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_x     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[0]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xx(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_xy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_xz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Ux) for sphere inclusion (3body conditions);
void CLame3D::jump1_sphere_x_3body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_x     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[0]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xx(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_xy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_xz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		B[i].shape->cpy_xx(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
		B[i].shape->cpy_xy(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		B[i].shape->cpy_xz(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -sqr(RR)*alpha*G1*(6.+k-4.*nju);

		B[i].shape->cpy_x		 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->deriv, 0., P[0]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->deriv, 0., P[1]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->deriv, 0., P[2]);
		B[i].shape->adm		 (0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(6.+k-4.*nju)*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[0], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[0], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[0], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*RR*alpha*G1*(k+2.-(1.-nju)*(2.*k+1.));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[0]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[0]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[0]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(k+2.-2.*(1.-nju)*(2.*k+3.))*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

/////////////////////////////////////////////////
//...realization of classical displacements (Uy);
void CLame3D::jump1_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_y     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[1]*alpha);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0]*alpha, 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Uy) for sphere inclusion (2body);
void CLame3D::jump1_sphere_y_2body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_y     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[1]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0,  -l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_yy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_yz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Uy) for sphere inclusion (3body);
void CLame3D::jump1_sphere_y_3body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_y     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[1]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[2]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0,  -l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_yy(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_yz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		B[i].shape->cpy_xy(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
		B[i].shape->cpy_yy(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		B[i].shape->cpy_yz(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -sqr(RR)*alpha*G1*(6.+k-4.*nju);

		B[i].shape->cpy_y		 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->deriv, 0., P[0]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->deriv, 0., P[1]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->deriv, 0., P[2]);
		B[i].shape->adm		 (0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(6.+k-4.*nju)*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[1], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[1], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[1], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*RR*alpha*G1*(k+2.-(1.-nju)*(2.*k+1.));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[1]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[1]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[1]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(k+2.-2.*(1.-nju)*(2.*k+3.))*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

/////////////////////////////////////////////////
//...realization of classical displacements (Uz);
void CLame3D::jump1_classic_z(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_z     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[2]*alpha);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0]*alpha, 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1]*alpha, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Uz) for sphere inclusion (2body);
void CLame3D::jump1_sphere_z_2body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_z     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[2]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0,  -l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_yz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_zz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical displacements (Uz) for sphere inclusion (3body);
void CLame3D::jump1_sphere_z_3body(double * P, int i, int m) 
{ 
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, j = NUM_HESS+4+solver.id_norm, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_z     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, alpha+1., P[2]*alpha);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, G1, 0.);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, G1, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha, 0.);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0,  -l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		B[i].shape->cpy_xz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, -1));
		B[i].shape->cpy_yz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0));
		B[i].shape->cpy_zz(1, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		B[i].shape->cpy_xz(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
		B[i].shape->cpy_yz(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		B[i].shape->cpy_zz(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -sqr(RR)*alpha*G1*(6.+k-4.*nju);

		B[i].shape->cpy_z		 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->deriv, 0., P[0]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->deriv, 0., P[1]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->deriv, 0., P[2]);
		B[i].shape->adm		 (0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(6.+k-4.*nju)*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[2], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[2], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[2], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*RR*alpha*G1*(k+2.-(1.-nju)*(2.*k+1.));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[2]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[2]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[2]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(k+2.-2.*(1.-nju)*(2.*k+3.))*(2.*k+1.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Ux);
void CLame3D::jump2_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_xx	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[4]);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[5]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(2.*alpha+1.)*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.+alpha)*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*(1.+alpha)*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[4]*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[5]*alpha*G1);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Ux) for sphere inclusion (2body conditions);
void CLame3D::jump2_sphere_x_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[5]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Ux) for sphere inclusion (3body conditions);
void CLame3D::jump2_sphere_x_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[5]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
 		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[0]*P[3]+nr, P[0]*P[0]);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[0]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0]*P[1], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[0]*P[4]);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[1]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0]*P[2], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[0]*P[5]);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);
 
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[0]*P[3], P[0]*P[0]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[0]*P[4], P[0]*P[1]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[0]*P[5], P[0]*P[2]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+5.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*4.*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_x     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[0]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (0, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*RR*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[0]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[0]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[0]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -2.*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[3]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[3]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[3]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+3.-(k+2.)*.5/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[3], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[3], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[3], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uy);
void CLame3D::jump2_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_yy	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[5]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.+alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(2.*alpha+1.)*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*(1.+alpha)*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, P[5]*alpha*G1);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uy) for sphere inclusion (2body conditions);
void CLame3D::jump2_sphere_y_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uy) for sphere inclusion (3body conditions);
void CLame3D::jump2_sphere_y_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[1]*P[4]+nr, P[1]*P[1]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[1]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[0], 0.);
		B[i].shape->adm       (0, B[i].shape->deriv, P[1]*P[3]);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[0]*nr);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1]*P[2], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[1]*P[5]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[1]*P[3], P[1]*P[0]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[1]*P[4], P[1]*P[1]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[1]*P[5], P[1]*P[2]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+5.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[1]*P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[1]*P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[1]*P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_yy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_yz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*4.*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_y     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[4], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[3]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[1]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[3]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2], 0.);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*RR*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[1]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[1]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[1]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -2.*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[4]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[4]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[4]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+3.-(k+2.)*.5/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[4], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[4], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[4], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uz);
void CLame3D::jump2_classic_z(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift), 
			 G1 = 1./get_param(NUM_SHEAR+shift); m += solver.id_norm;
	B[i].shape->cpy_zz	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[5], 0.);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.+alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.+alpha)*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*(2.*alpha+1.)*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->p_cpy, P[3]*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, P[4]*alpha*G1);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uz) for sphere inclusion (2body conditions);
void CLame3D::jump2_sphere_z_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[5], 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(2.*alpha+1.)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_z     (l, B[i].shape->p_cpy, P[3]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uz) for sphere inclusion (3body conditions);
void CLame3D::jump2_sphere_z_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift), G1 = 1./get_param(NUM_SHEAR+shift), nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[5], 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2]*alpha*G1, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G1);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G1);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(2.*alpha+1.)*G1);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
		B[i].shape->adm_z     (l, B[i].shape->p_cpy, P[3]*alpha*G1);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*alpha*G1);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*G1*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[2]*P[5]+nr, P[2]*P[2]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[0], 0.);
		B[i].shape->adm       (0, B[i].shape->deriv, P[2]*P[3]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[0]*nr);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2]*P[1], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[2]*P[4]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[1]*nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= 2.*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[2]*P[3], P[2]*P[0]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[2]*P[4], P[2]*P[1]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[2]*P[5], P[2]*P[2]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+5.);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[2]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[2]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[2]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_yz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_zz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, nr, 0.);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*4.*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_z     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[5], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[3]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[2]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[3]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*alpha*G1*(2.*k+1.)*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*RR*alpha*G1*(6.+k-4.*nju);
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[2]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[2]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[2]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -2.*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[5]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[5]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[5]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= G1*(2.*k+1.)*(2.*k+3.-(k+2.)*.5/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[5], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[5], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[5], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*G1*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Px);
void CLame3D::jump4_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift)*2.; m += solver.id_norm;
	B[i].shape->cpy_xx	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[4]);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[5]);

	B[i].shape->admittance(B[i].shape->deriv, NULL, alpha, 0.);
	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.+alpha));
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*(1.+alpha));

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1], 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[4]*( 1.+alpha));
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2], 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, P[3]*(-2.*alpha-1.));
	B[i].shape->adm_x     (B[i].shape->deriv, P[5]*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Px) for sphere inclusion (2body);
void CLame3D::jump4_sphere_x_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[5]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Px) for sphere inclusion (3body);
void CLame3D::jump4_sphere_x_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[5]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
 		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[0]*P[3]+nr, P[0]*P[0]);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[0]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[1], 0.);
		B[i].shape->adm       (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_x     (0, B[i].shape->deriv, P[1]*nr);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0]*P[2], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[0]*P[5]);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(2.*k+1.-(k+2.)/(1.-nju));

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3]+nr, 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3]+nr, 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*P[3]+nr, 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[0]*P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[0]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k-3.-2.*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_x     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[0]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (0, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k+1.)*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0); //...deriv_N_xx;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*RR*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[0]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[0]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[0]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[3]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[3]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[3]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(3.*k+4.+2.*nju*sqr(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[3], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[3], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[3], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*(k-1.)+nju*(k+2.)*(2.*k+3.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Py);
void CLame3D::jump4_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift)*2.; m += solver.id_norm;
	B[i].shape->cpy_yy	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[5]);

	B[i].shape->admittance(B[i].shape->deriv, NULL, alpha, 0.);
	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[1], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.+alpha));
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*(1.+alpha));

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*( 1.+alpha));
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2], 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, P[4]*(-2.*alpha-1.));
	B[i].shape->adm_y     (B[i].shape->deriv, P[5]*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Py) for sphere inclusion (2body);
void CLame3D::jump4_sphere_y_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Py) for sphere inclusion (3body);
void CLame3D::jump4_sphere_y_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]*(1.+alpha));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[1]*P[4]+nr, P[1]*P[1]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[1]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[0], 0.);
		B[i].shape->adm       (0, B[i].shape->deriv, P[1]*P[3]);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[0]*nr);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1]*P[2], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[1]*P[5]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(2.*k+1.-(k+2.)/(1.-nju));

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[1]*P[4]+nr);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[1]*P[4]+nr);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[1]*P[4]+nr);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[1]*P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k-3.-2.*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_y     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[4], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[3]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[1]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[3]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2], 0.);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[5]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k+1.)*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2); //...deriv_N_yy;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*RR*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[1]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[1]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[1]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[4]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[4]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[4]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(3.*k+4.+2.*nju*sqr(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[4], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[4], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[4], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*(k-1.)+nju*(k+2.)*(2.*k+3.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Pz);
void CLame3D::jump4_classic_z(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = get_param(NUM_SHEAR+2+shift)*2.; m += solver.id_norm;
	B[i].shape->cpy_zz	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[5], 0.);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[4]);

	B[i].shape->admittance(B[i].shape->deriv, NULL, alpha, 0.);
	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[2], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.+alpha));
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.+alpha));
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 2));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[5]*(-2.*alpha-1.));
	B[i].shape->adm_z     (B[i].shape->p_cpy, P[3]*( 1.+alpha));
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[1], 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, P[5]*(-2.*alpha-1.));
	B[i].shape->adm_z     (B[i].shape->deriv, P[4]*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 1, 1));
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Pz) for sphere inclusion (2body);
void CLame3D::jump4_sphere_z_2body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[5], 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]);
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[5]*(-2.*alpha-1.));
		B[i].shape->adm_z     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*(-2.*alpha-1.));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 1, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		add_collocation(i, m, j+1, 1, B1);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...realization surface forces of classical displacements (Pz) for sphere inclusion (3body);
void CLame3D::jump4_sphere_z_3body(double * P, int i, int m)
{
	double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]); m += solver.id_norm;
	int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;

	if (forced_phase != -777) phase = forced_phase;

	int shift = (1-phase)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, j = num_hess+4, N_mpl = B[i].shape->get_N(0), k, l, num;
	double alpha = get_param(NUM_SHEAR+2+shift)*2., nju = get_param(NUM_SHEAR+1+shift), * p;

///////////////////////
//...простой потенциал;
	for (l = 0; l < 2; l++) if (! l || phase) {
		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[5], 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[2], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[5]);
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 2-l));

		B[i].shape->cpy       (l, B[i].shape->p_cpy, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[5]*(-2.*alpha-1.));
		B[i].shape->adm_z     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[5]*(-2.*alpha-1.));
		B[i].shape->adm_z     (l, B[i].shape->deriv, P[4]*( 1.+alpha));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+l], 0, -l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+l], 0, 1-l));
	}
	if (phase > 0) {
		solver.admittance(i, m-solver.id_norm, 0., j-solver.id_norm, 1.);
		add_collocation  (i, m, j+1, 2, C2);
		add_collocation  (i, m, j+1, 3, A2);
	}
	else 
	if (phase < 0) {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j,   1, NULL);
		add_collocation  (i, m, j+1, 2, A1);
	}
	else {
		solver.admittance(i, m-solver.id_norm, 0.);
		add_collocation  (i, m, j, 0, NULL);
	}
///////////////////////////////////////
//...сингулярный градиентный потенциал;
	if (phase) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);

		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -alpha*(6.+k-4.*nju);
		if (phase > 0)	add_collocation(i, m, j+1, 3, B2);
		else				add_collocation(i, m, j+1, 2, B1);
	}
//////////////////////////////////////
//...регулярный градиентный потенциал;
	if (phase <= 0) {
		double nr = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
		B[i].shape->cpy_x     (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, P[2]*P[5]+nr, P[2]*P[2]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[2]*nr);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 2));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[0], 0.);
		B[i].shape->adm       (0, B[i].shape->deriv, P[2]*P[3]);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[0]*nr);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[2]*P[1], 0.);
		B[i].shape->adm       (0, B[i].shape->p_cpy, P[2]*P[4]);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[1]*nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(2.*k+1.-(k+2.)/(1.-nju));

		B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_xy    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_xz    (0, B[i].shape->deriv, P[2]*P[5]+nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));

		B[i].shape->cpy_xy	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_yy    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[2]*P[5]+nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[2]*P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[2]*P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[2]*P[5]+nr);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k-3.-2.*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_xz	 (0, B[i].shape->deriv);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_yz    (0, B[i].shape->deriv, P[4]);
		B[i].shape->adm_zz    (0, B[i].shape->deriv, P[5]);
		B[i].shape->cpy_z     (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, 2.*P[5], 0.);
		B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[3]);
		B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, 1., P[2]);
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2));

		B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_z     (0, B[i].shape->deriv, P[3]);
		B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_z     (0, B[i].shape->p_cpy, P[4]);
		B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j+1], 0));
		B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1));
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*k+1.)*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0); //...deriv_N_xz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), p, 0., 1.);
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2); //...deriv_N_zz;
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), p, 0., 1.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= RR*RR*(2.+.5*(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[2]*nr, 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[2]*nr, 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[2]*nr, 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -4.*(2.*k+1.-(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy		 (0, B[i].shape->p_cpy);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), B[i].shape->p_cpy, 0., P[0]*P[5]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), B[i].shape->p_cpy, 0., P[1]*P[5]);
		B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), B[i].shape->p_cpy, 0., P[2]*P[5]);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= (2.*k+1.)*(3.*k+4.+2.*nju*sqr(k+2.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);

		B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 0), NULL, P[5], 0.);
		B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 1), NULL, P[5], 0.);
		B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2)); B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j+1], 0, 2), NULL, P[5], 0.);
		for (num = 0; num < 3; num++)
		for (p = B[i].shape->FULL(solver.hh[i][0][j+1], 0, num), k = 0; k <= N_mpl; k++)
		for (l = 0; l <= k*2; l++) p[k*k+l] *= -RR*(2.*(k-1.)+nju*(k+2.)*(2.*k+3.)/(1.-nju));
		solver.admittance(i, j-solver.id_norm, 1., j+1-solver.id_norm, 1.);
		if (! phase) add_collocation(i, m, j, 2, NULL, C0);
		else			 add_collocation(i, m, j, 2, NULL, C1);
	}
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
void CLame3D::hessian_deriv_N(double * P, int i)
{
	int l, num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->cpy_xx(); B[i].shape->deriv_N(); 
	B[i].shape->cpy_xx();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+2*l], 0, -l));
	B[i].shape->cpy_xy(); B[i].shape->deriv_N();
	B[i].shape->cpy_xy();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+2*l], 0, 1-l));
	B[i].shape->cpy_yy(); B[i].shape->deriv_N();
	B[i].shape->cpy_yy();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+2*l], 0, 2-l));
	B[i].shape->cpy_xz(); B[i].shape->deriv_N();
	B[i].shape->cpy_xz();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1+2*l], 0, -l));
	B[i].shape->cpy_yz(); B[i].shape->deriv_N();
	B[i].shape->cpy_yz();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1+2*l], 0, 1-l));
	B[i].shape->cpy_zz(); B[i].shape->deriv_N();
	B[i].shape->cpy_zz();
	for (l = 0; l < 2; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1+2*l], 0, 2-l));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CLame3D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 
							solver.hh[i][0][m+2][j] };
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
		solver.hh[i][0][m+2][j] = P[2];
	}
}

void CLame3D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 
							solver.hh[i][0][m+2][j] };
		B[i].shape->norm_common(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
		solver.hh[i][0][m+2][j] = P[2];
	}
}

////////////////////////////////////////////////////////////////////////////
//...realization of common jump boundary condition for matrix and inclusion;
int CLame3D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = solver.id_norm, n_rhs = 0;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);
		
////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == MIN_HIT || p4 == 2.) {
				hy = nd->nY[l]*hx;
				hz = nd->nZ[l]*hx;
				hx *= nd->nX[l];
			}
			else 
			if (p4 == NUMS_BND) { //...специальный случай -- одноосное растяжение;
				double nju = get_param(NUM_SHEAR+1), G0 = .5/get_param(NUM_SHEAR), 
						 BBB = G0/(1.+nju), AAA = -nju*BBB;
				hx = nd->X[l]*AAA;
				hy = nd->Y[l]*AAA; 
				hz = nd->Z[l]*BBB;
			}
			else 
			if (p4 == (double)(NORMS_BND-SPECIAL_BND) && hx != 0.) { //...антиплоский сдвиг в плоскости XZ;
				hx = nd->Z[l];
				hy = hz = 0.;
			}
			if (p4 == (double)(NORMS_BND-SPECIAL_BND) && hy != 0.) { //...антиплоский сдвиг в плоскости YZ;
				hy = nd->Z[l];
				hx = hz = 0.;
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_hess(P, 1);

			if (p4 == (double)(NORMS_BND-SPECIAL_BND) || p4 == (double)(SKEWS_BND-SPECIAL_BND)) { //...смешанные краевые условия; 
				jump1_classic_x (P, i, 0); solver.admittance(i, 0, G1); solver.admittance (i, 7, 0., 0, P[3]);
				jump1_classic_y (P, i, 1); solver.admittance(i, 1, G1); solver.admittance (i, 7, 1., 1, P[4]);
				jump1_classic_z (P, i, 2); solver.admittance(i, 2, G1); solver.admittance (i, 7, 1., 2, P[5]);
				jump_make_common(i, 0);
				solver.admittance (i, 0, 1., 7, -nd->nX[l]);
				solver.admittance (i, 1, 1., 7, -nd->nY[l]);
				solver.admittance (i, 2, 1., 7, -nd->nZ[l]);

				jump4_classic_x (P, i, 4); solver.admittance (i, 3, 0., 4, P[3]);
				jump4_classic_y (P, i, 5); solver.admittance (i, 3, 1., 5, P[4]);
				jump4_classic_z (P, i, 6); solver.admittance (i, 3, 1., 6, P[5]);
				jump_make_common(i, 4);
				solver.admittance (i, 4, 1., 3, -nd->nX[l]);
				solver.admittance (i, 5, 1., 3, -nd->nY[l]);
				solver.admittance (i, 6, 1., 3, -nd->nZ[l]); n_rhs = 1;

				if (p4 == (double)(SKEWS_BND-SPECIAL_BND)) { //...меняем компоненты;
					solver.admittance (i, 0, 0., 4, 1.);
					solver.admittance (i, 1, 0., 5, 1.);
					solver.admittance (i, 2, 0., 6, 1.);
					solver.admittance (i, 3, 0., 7, 1.); n_rhs = 0;
				}
			}
			else { //...традиционные краевые условия;
				if (p4 == MIN_HIT || p4 == NUMS_BND || p4 == MAX_HIT) {
					jump1_classic_x(P, i, 0); solver.admittance(i, 0, G1);
					jump1_classic_y(P, i, 1); solver.admittance(i, 1, G1);
					jump1_classic_z(P, i, 2); solver.admittance(i, 2, G1);
				}
				else
				if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
					jump4_classic_x(P, i, 0);
					jump4_classic_y(P, i, 1);
					jump4_classic_z(P, i, 2);
				}
				jump_make_common(i, 0); n_rhs = 0;
			}

////////////////////////////
//...composition functional;
			if (p4 == (double)(NORMS_BND-SPECIAL_BND) || p4 == (double)(SKEWS_BND-SPECIAL_BND))
			solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],	  f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
			if (p4 != (double)(SKEWS_BND-SPECIAL_BND)) {
				if (p4 == MIN_HIT || p4 == (double)(NORMS_BND-SPECIAL_BND) || p4 == NUMS_BND || p4 == MAX_HIT) f *= G1;
				if (fabs(hx) > EE)
					solver.to_equationHH(i, n_rhs, solver.hh[i][0][m],   hx*f);
				if (fabs(hy) > EE)
					solver.to_equationHH(i, n_rhs, solver.hh[i][0][m+1], hy*f);
				if (fabs(hz) > EE)
					solver.to_equationHH(i, n_rhs, solver.hh[i][0][m+2], hz*f);
			}
			else 
			if (fabs(hx) > EE)
				solver.to_equationHH(i, n_rhs, solver.hh[i][0][m+3], hx*f*G1);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
int CLame3D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz, 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1,  requl = 1.;
      int id_isolated = 0, m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);

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
				TX = TY= TZ = hz = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hz = -g1*AZ; break;
					case 6: TZ = -AZ; hz =  g1*AZ; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY;
				B[k].mp[3] -= TZ; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+5; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					for (int num = m+6; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);
				
				jump1_classic_x (P, i, 0); jump1_classic_y (P, i, 1); jump1_classic_z(P, i, 2); 
				jump2_classic_x (P, i, 3); jump2_classic_y (P, i, 4); jump2_classic_z(P, i, 5); 
				jump_make_common(i, 0);	   jump_make_common(i, 3);	

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 6, 1., 0, f*G1); 
					solver.admittance (i, 7, 1., 1, f*G1);
					solver.admittance (i, 8, 1., 2, f*G1);
				}
				solver.admittance(i, 0, g1, 3, g2); solver.admittance(i, 3, g0, 0, f1); 
				solver.admittance(i, 1, g1, 4, g2); solver.admittance(i, 4, g0, 1, f1); 
				solver.admittance(i, 2, g1, 5, g2); solver.admittance(i, 5, g0, 2, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				
				B[k].shape->parametrization_hess(P, 1);

				jump1_classic_x (P, k, 0); jump1_classic_y (P, k, 1); jump1_classic_z(P, k, 2); 
				jump2_classic_x (P, k, 3); jump2_classic_y (P, k, 4); jump2_classic_z(P, k, 5); 
				jump_make_common(k, 0);	   jump_make_common(k, 3);	

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 6, 1., 0, -f*G1); 
					solver.admittance (k, 7, 1., 1, -f*G1);
					solver.admittance (k, 8, 1., 2, -f*G1);
				}
				solver.admittance(k, 0, g1, 3, g2); solver.admittance(k, 3, g0, 0, f1); 
				solver.admittance(k, 1, g1, 4, g2); solver.admittance(k, 4, g0, 1, f1); 
				solver.admittance(k, 2, g1, 5, g2); solver.admittance(k, 5, g0, 2, f1); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				if (fabs(hz) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+2],  hz*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m],    hz*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m+5], -hz*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+3], -hz*f);
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
			}
      }
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
		}
/////////////////////////////////////////////////////////////////
//memset(solver.hh[i][0][0], 0, solver.dim[i]*sizeof(double));
//B[i].shape->FULL(solver.hh[i][0][0], 0, 2)[1] = 1.;	
/////////////////////////////////////////////////////////////////
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для периодической задачи (один блок);
int CLame3D::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, AZ, f, P[6], 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1;
      int id_isolated = 1, m  = solver.id_norm, id_dir;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l] && 
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
			B[i].shape->parametrization_hess(P, 1);
			
			jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); 
			jump2_x(P, i, 3); jump2_y(P, i, 4); jump2_z(P, i, 5); 
			jump_make_common(i, 0);	   jump_make_common(i, 3);	

			B[i].shape->make_common(P);
			B[i].shape->norm_common(P+3);

			if (id_dir == 1) B[i].mp[1] -= AX; else
			if (id_dir == 3) B[i].mp[2] -= AY; else
			if (id_dir == 5) B[i].mp[3] -= AZ; 
			B[i].shape->set_local_P0(B[i].mp+1);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			B[i].shape->parametrization_hess(P, 1);

			jump1_x(P, i, 6); jump1_y(P, i, 7);  jump1_z(P, i, 8); 
			jump2_x(P, i, 9); jump2_y(P, i, 10); jump2_z(P, i, 11); 
			jump_make_common(i, 6);	    jump_make_common(i, 9);	

			solver.admittance(i, 0, 1., 6, -1.); solver.admittance (i, 6, 2., 0, 1.); solver.admittance (i, 3, 1.,  9, -1.);
			solver.admittance(i, 1, 1., 7, -1.); solver.admittance (i, 7, 2., 1, 1.); solver.admittance (i, 4, 1., 10, -1.);
			solver.admittance(i, 2, 1., 8, -1.); solver.admittance (i, 8, 2., 2, 1.); solver.admittance (i, 5, 1., 11, -1.);

			solver.admittance(i, 0, g1, 3, g2); solver.admittance(i, 3, g0, 0, f1); 
			solver.admittance(i, 1, g1, 4, g2); solver.admittance(i, 4, g0, 1, f1); 
			solver.admittance(i, 2, g1, 5, g2); solver.admittance(i, 5, g0, 2, f1); 

///////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов;
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);

			solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);

			if (id_dir == 5) {//...регуляризация матрицы и правая часть;
				solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
				solver.to_equationDD(i, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
				solver.to_equationDD(i, solver.hh[i][0][m+8], solver.hh[i][0][m+8], f);

				solver.to_equationHH(i, 0,	solver.hh[i][0][m+2], -g1*AZ*f);
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], -g2*AZ*f);
				solver.to_equationHH(i, 1,	solver.hh[i][0][m],   -g1*AZ*f);
				solver.to_equationHH(i, 1, solver.hh[i][0][m+3], -g2*AZ*f);
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
int CLame3D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f0 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int id_isolated = 0, m = solver.id_norm;
		if (id_isolated) {
			g0 = g1 = G1;
			f0 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);
				
				jump1_classic_x(P, i, 0); 
				jump2_classic_x(P, i, 3); 
				solver.admittance (i, 0, g1, 3, g2); solver.admittance(i, 3, g0, 0, f0); 
				
				jump1_classic_y(P, i, 1); 
				jump2_classic_y(P, i, 4); 
				solver.admittance (i, 1, g1, 4, g2); solver.admittance(i, 4, g0, 1, f0); 

				jump1_classic_z(P, i, 2); 
				jump2_classic_z(P, i, 5); 
				solver.admittance (i, 2, g1, 5, g2); solver.admittance(i, 5, g0, 2, f0); 
				
				jump_make_common(i, 0);
				jump_make_common(i, 3);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				
				B[k].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, k, 0); 
				jump2_classic_x(P, k, 3); 
				solver.admittance (k, 0, g1, 3, g2); solver.admittance(k, 3, g0, 0, f0); 

				jump1_classic_y(P, k, 1); 
				jump2_classic_y(P, k, 4); 
				solver.admittance (k, 1, g1, 4, g2); solver.admittance(k, 4, g0, 1, f0); 

				jump1_classic_z(P, k, 2); 
				jump2_classic_z(P, k, 5); 
				solver.admittance (k, 2, g1, 5, g2); solver.admittance(k, 5, g0, 2, f0); 
				
				jump_make_common(k, 0);
				jump_make_common(k, 3);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for all blocks;
int CLame3D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f0 = -1., f1 = 1., f2 = .5, g1 = G1*.5;
      int id_isolated = 1, id_flag = 1, m = solver.id_norm;
		if (id_isolated) {
			f1 = f2 = 0.;
			g1 = G1; f0 = 1.;	
			id_flag = B[i].link[NUM_PHASE] == -2;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);

				jump1_x(P, i, 0); 
				jump4_x(P, i, 3); 
				solver.admittance (i, 0, g1, 3, f2); solver.admittance(i, 3, f0, 0, f1); 

				jump1_y(P, i, 1); 
				jump4_y(P, i, 4); 
				solver.admittance (i, 1, g1, 4, f2); solver.admittance(i, 4, f0, 1, f1); 

				jump1_z(P, i, 2); 
				jump4_z(P, i, 5); 
				solver.admittance (i, 2, g1, 5, f2); solver.admittance(i, 5, f0, 2, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 3);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);

				jump1_x(P, k, 0); 
				jump4_x(P, k, 3); 
				solver.admittance (k, 0, g1, 3, f2); solver.admittance(k, 3, f0, 0, f1); 

				jump1_y(P, k, 1); 
				jump4_y(P, k, 4); 
				solver.admittance (k, 1, g1, 4, f2); solver.admittance(k, 4, f0, 1, f1); 

				jump1_z(P, k, 2); 
				jump4_z(P, k, 5); 
				solver.admittance (k, 2, g1, 5, f2); solver.admittance(k, 5, f0, 2, f1); 

				jump_make_common(k, 0);
				jump_make_common(k, 3);

///////////////////////////
//...composition functional;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
					solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m],   f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m],   f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for Eshrlby problem;
int CLame3D::trans_esh(CGrid * nd, int i, int k, int id_local)
{
	if (1 && nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_SHEAR), f, P[6], hx, hy, hz, px, py, pz, 
				 f0 = -1., f1 = 1., f2 = .5, g1 = G1*.5;
      int id_isolated = 0, id_flag = 1, m = solver.id_norm;
		if (id_isolated) {
			f1 = f2 = 0.;
			g1 = G1; f0 = 1.;	
			id_flag = B[i].link[NUM_PHASE] == -2;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

//////////////////////////////////
//...тестовая проверка квадратуры;
#ifdef ___CONSTRUCTION___
		double sum = 0., sum_X = 0., sum_Y = 0., sum_Z = 0., surf, volm, 
			RR = fabs(B[i].bar->ce[0]->mp[7]), rr = fabs(B[i].bar->ce[0]->mp[8]), p = acos(RR < rr ? RR/rr : rr/RR);
		for (int k = 0; k < nd->N; k++) {
			sum_X += nd->X[k]*nd->nX[k]*nd->get_param(0, k);
			sum_Y += nd->Y[k]*nd->nY[k]*nd->get_param(0, k);
			sum_Z += nd->Z[k]*nd->nZ[k]*nd->get_param(0, k);
			sum += nd->get_param(0, k);
		}
		volm = 4./3.*M_PI*sqr(rr)*RR;
		surf = 2.*M_PI*(RR < rr ? sqr(rr)+sqr(RR)/sin(p)*log((1.+sin(p))/cos(p)) :
										  sqr(rr)+fabs(RR*rr)*p/sin(p));
#endif
////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);

				jump1_x(P, i, 0); 
				jump4_x(P, i, 3); 
				solver.admittance (i, 0, g1, 3, f2); solver.admittance(i, 3, f0, 0, f1); 

				jump1_y(P, i, 1); 
				jump4_y(P, i, 4); 
				solver.admittance (i, 1, g1, 4, f2); solver.admittance(i, 4, f0, 1, f1); 

				jump1_z(P, i, 2); 
				jump4_z(P, i, 5); 
				solver.admittance (i, 2, g1, 5, f2); solver.admittance(i, 5, f0, 2, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 3);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);

				jump1_x(P, k, 0); 
				jump4_x(P, k, 3); 
				solver.admittance (k, 0, g1, 3, f2); solver.admittance(k, 3, f0, 0, f1); 

				jump1_y(P, k, 1); 
				jump4_y(P, k, 4); 
				solver.admittance (k, 1, g1, 4, f2); solver.admittance(k, 4, f0, 1, f1); 

				jump1_z(P, k, 2); 
				jump4_z(P, k, 5); 
				solver.admittance (k, 2, g1, 5, f2); solver.admittance(k, 5, f0, 2, f1); 

				jump_make_common(k, 0);
				jump_make_common(k, 3);

///////////////////////////
//...composition functional;
				double nju = get_param(NUM_SHEAR+1), G0 = .5/get_param(NUM_SHEAR), BBB = G0/(1.+nju), AAA = -nju*BBB;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
					solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m],   f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m],   f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				}
/////////////////////////////////////
//...однородное поле -- правая часть;
				if (B[i].link[NUM_PHASE] == -2) f = -f;
				hx = nd->X[l]*AAA; px = 0.;			hx = hx*g1+px*f2; px = px*f0+hx*f1;
				hy = nd->Y[l]*AAA; py = 0.;			hy = hy*g1+py*f2; py = py*f0+hy*f1;
				hz = nd->Z[l]*BBB; pz = nd->nZ[l];	hz = hz*g1+pz*f2; pz = pz*f0+hz*f1;
				solver.to_equationHH(i, 0, solver.hh[i][0][m],   -hx*f);
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], -hy*f);
				solver.to_equationHH(i, 0, solver.hh[i][0][m+2], -hz*f);
				solver.to_equationHH(k, 0, solver.hh[k][0][m+3],  px*f);
				solver.to_equationHH(k, 0, solver.hh[k][0][m+4],  py*f);
				solver.to_equationHH(k, 0, solver.hh[k][0][m+5],  pz*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CLame3D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz, requl = 1.;
      int	 m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hz = -G1*AZ; break;
					case 6: TZ = -AZ; hz =  G1*AZ; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY;
				B[k].mp[2] -= TZ; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+5; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					for (int num = m+6; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); jump1_classic_z(P, i, 2); 
				jump4_classic_x(P, i, 3); jump4_classic_y(P, i, 4); jump4_classic_z(P, i, 5); 
				jump_make_common(i, 0);	  jump_make_common(i, 3);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 6, 1., 0, f*G1); 
					solver.admittance (i, 7, 1., 1, f*G1);
					solver.admittance (i, 8, 1., 2, f*G1);
				}
				solver.admittance(i, 0, G1); solver.admittance(i, 1, G1); solver.admittance(i, 2, G1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); jump1_classic_z(P, k, 2); 
				jump4_classic_x(P, k, 3); jump4_classic_y(P, k, 4); jump4_classic_z(P, k, 5); 
				jump_make_common(k, 0);	  jump_make_common(k, 3);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 6, 1., 0, f*G1); 
					solver.admittance (k, 7, 1., 1, f*G1);
					solver.admittance (k, 8, 1., 2, f*G1);
				}
				solver.admittance(k, 0, G1); solver.admittance(k, 1, G1); solver.admittance(k, 2, G1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

				if (fabs(hz) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+2],  hz*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m],    hz*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m+2], -hz*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m],   -hz*f);
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+3],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+3], -f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+4], -f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+5], -f);
				
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
			}
      }
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+6], solver.hh[k0][0][m+6], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+7], solver.hh[k0][0][m+7], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+8], solver.hh[k0][0][m+8], requl);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CLame3D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == MIN_HIT || p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p4 == NUMS_BND) { //...специальный случай -- одноосное растяжение;
				double nju = get_param(NUM_SHEAR+1), G1 = .5/get_param(NUM_SHEAR), 
						 BBB = G1/(1.+nju), AAA = -nju*BBB;
				hx = nd->X[l]*AAA;
				hy = nd->Y[l]*AAA; 
				hz = nd->Z[l]*BBB;
			}
			else
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

//////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_hess(P, 1);

			jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); jump1_classic_z(P, i, 2); 
			jump4_classic_x(P, i, 3); jump4_classic_y(P, i, 4); jump4_classic_z(P, i, 5);
			jump_make_common(i, 0);	  jump_make_common(i, 3);
			solver.admittance(i, 0, G1); 
			solver.admittance(i, 1, G1); 
			solver.admittance(i, 2, G1); 

////////////////////////////////////
//...composition collocation vector;
			if (p4 == (double)(SKEWS_BND-SPECIAL_BND)) { //...специальные краевые условия на ячейке периодичности; 
				jump2_classic_x (P, i, 6); solver.admittance(i, 6, G1);
				jump2_classic_y (P, i, 7); solver.admittance(i, 7, G1);
				jump2_classic_z (P, i, 8); solver.admittance(i, 8, G1);
				jump_make_common(i, 6);

				if (fabs(fabs(P[3])-1.) < EE_ker) solver.admittance (i, 6, 0., 0, 1.); else
				if (fabs(fabs(P[4])-1.) < EE_ker) solver.admittance (i, 7, 0., 1, 1.); else
				if (fabs(fabs(P[5])-1.) < EE_ker) solver.admittance (i, 8, 0., 2, 1.);
			}
			else
			if (p4 == MIN_HIT || p4 == NUMS_BND || p4 == MAX_HIT) {
				solver.admittance(i, 6, 0., 0, 1.);
				solver.admittance(i, 7, 0., 1, 1.);
				solver.admittance(i, 8, 0., 2, 1.);
			}
			else
			if (p4 == (double)(NORMS_BND-SPECIAL_BND)) {
				jump2_classic_x(P, i, 6); solver.admittance(i, 6, G1);
				jump2_classic_y(P, i, 7); solver.admittance(i, 7, G1);
				jump2_classic_z(P, i, 8); solver.admittance(i, 8, G1);
				jump_make_common(i, 6);
			}

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p4 == MIN_HIT || p4 == (double)(SKEWS_BND-SPECIAL_BND) || p4 == (double)(NORMS_BND-SPECIAL_BND) || p4 == NUMS_BND || p4 == MAX_HIT) {
				solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
				solver.to_equationDD(i, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
				solver.to_equationDD(i, solver.hh[i][0][m+8], solver.hh[i][0][m+8], f);

				if (fabs(hx) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+6], hx*G1*f);
				if (fabs(hy) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+7], hy*G1*f);
				if (fabs(hz) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+8], hz*G1*f);
			}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+3], f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4], f);
			solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5], f);

			if (1. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				if (fabs(hx) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m],   -hx*f);
				if (fabs(hy) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+1], -hy*f);
				if (fabs(hz) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+2], -hz*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
int CLame3D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
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

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); jump1_classic_z(P, i, 2); 
				jump4_classic_x(P, i, 3); jump4_classic_y(P, i, 4); jump4_classic_z(P, i, 5); 
				jump_make_common(i, 0);	  jump_make_common(i, 3);
				solver.admittance(i, 0, G1); 
				solver.admittance(i, 1, G1); 
				solver.admittance(i, 2, G1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); jump1_classic_z(P, k, 2); 
				jump4_classic_x(P, k, 3); jump4_classic_y(P, k, 4); jump4_classic_z(P, k, 5); 
				jump_make_common(k, 0);	  jump_make_common(k, 3);
				solver.admittance(k, 0, G1); 
				solver.admittance(k, 1, G1); 
				solver.admittance(k, 2, G1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+3],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+3], -f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+4], -f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+5], -f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////
//...интегрирование НДС по границе блоков;
int CLame3D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6], U[3];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];	f = nd->get_param(0, l);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = solver.id_norm; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

////////////////////////////////////////////////////////////////////////////
//...вычисляем перемещения для формирования интеграла от тензора напряжений;
			B[i].shape->parametrization_hess(P, 1);

			jump1_x(P, i, 0); 
			jump1_y(P, i, 1); 
			jump1_z(P, i, 2); 

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   0);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 0);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 0);
			B[i].shape->norm_common(U);
			if ((B[i].type & ERR_CODE)	== CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS)) { 
				double RR = sqr(P[0])+sqr(P[1])+sqr(P[2]);
				int phase = RR > sqr(get_param(NUM_GEOMT+1)); if (! phase && RR > sqr(get_param(NUM_GEOMT))) phase = -1;
				shift = (1-phase)*NUM_SHIFT; G0 = get_param(NUM_SHEAR+shift);
				alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)); 
			}		
			K[0] += G0*(U[0]*nd->nX[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;
			K[1] += G0*(U[0]*nd->nY[l]+U[1]*nd->nX[l])*f;
			K[2] += G0*(U[1]*nd->nY[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;
			K[3] += G0*(U[0]*nd->nZ[l]+U[2]*nd->nX[l])*f;
			K[4] += G0*(U[1]*nd->nZ[l]+U[2]*nd->nY[l])*f;
			K[5] += G0*(U[2]*nd->nZ[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;

			K[6] +=  U[0]*nd->nX[l]*f;
			K[7] += (U[0]*nd->nY[l]+U[1]*nd->nX[l])*.5*f;
			K[8] +=  U[1]*nd->nY[l]*f;
			K[9] += (U[0]*nd->nZ[l]+U[2]*nd->nX[l])*.5*f;
			K[10]+= (U[1]*nd->nZ[l]+U[2]*nd->nY[l])*.5*f;
			K[11]+=  U[2]*nd->nZ[l]*f;

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   1);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 1);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 1);
			B[i].shape->norm_common(U);
			
			K[12] += G0*(U[0]*nd->nX[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;
			K[13] += G0*(U[0]*nd->nY[l]+U[1]*nd->nX[l])*f;
			K[14] += G0*(U[1]*nd->nY[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;
			K[15] += G0*(U[0]*nd->nZ[l]+U[2]*nd->nX[l])*f;
			K[16] += G0*(U[1]*nd->nZ[l]+U[2]*nd->nY[l])*f;
			K[17] += G0*(U[2]*nd->nZ[l]*2.+alpha*(U[0]*nd->nX[l]+U[1]*nd->nY[l]+U[2]*nd->nZ[l]))*f;
			K[18] +=  U[0]*nd->nX[l]*f;
			K[19] += (U[0]*nd->nY[l]+U[1]*nd->nX[l])*.5*f;
			K[20] +=  U[1]*nd->nY[l]*f;
			K[21] += (U[0]*nd->nZ[l]+U[2]*nd->nX[l])*.5*f;
			K[22] += (U[1]*nd->nZ[l]+U[2]*nd->nY[l])*.5*f;
			K[23] +=  U[2]*nd->nZ[l]*f;
			K[24+(-B[i].link[NUM_PHASE]-1)] += nd->Z[l]*nd->nZ[l]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////
//...интегрирование напряжений по объему блоков;
int CLame3D::rigidy2(CGrid * nd, int i, double * K)
{
	if (nd) {
      int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm, l;
      double f, P[6], F[3];

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = 1;
			P[1] = nd->Y[l];  P[4] = 0.;
			P[2] = nd->Z[l];  P[5] = 0.; f = nd->get_param(0, l);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

//////////////////////////
//...вычисляем напряжения;
			B[i].shape->parametrization_hess(P, 1);
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

			jump4_x(P, i, 0); 
			jump4_y(P, i, 1); 
			jump4_z(P, i, 2); 

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   0);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 0);
			B[i].shape->norm_common(F);
			K[0] += F[0]*f;
			K[1] += F[1]*f;
			K[3] += F[2]*f;

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   1);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 1);
			B[i].shape->norm_common(F);
			K[12] += F[0]*f;
			K[13] += F[1]*f;
			K[15] += F[2]*f;

			P[4] = 1.; P[3] = P[5] = 0.;
			B[i].shape->norm_local(P+3);
			B[i].shape->parametrization_hess(P, 1);
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

			jump4_x(P, i, 0); 
			jump4_y(P, i, 1); 
			jump4_z(P, i, 2); 

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   0);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 0);
			B[i].shape->norm_common(F);
			K[2] += F[1]*f;
			K[4] += F[2]*f;

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   1);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 1);
			B[i].shape->norm_common(F);
			K[14] += F[1]*f;
			K[16] += F[2]*f;

			P[5] = 1.; P[3] = P[4] = 0.;
			B[i].shape->norm_local(P+3);
			B[i].shape->parametrization_hess(P, 1);
			if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

			jump4_x(P, i, 0); 
			jump4_y(P, i, 1); 
			jump4_z(P, i, 2); 

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   0);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 0);
			B[i].shape->norm_common(F);
			K[5] += F[2]*f;

			F[0] = B[i].shape->potential(solver.hh[i][0][m],   1);
			F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 1);
			B[i].shape->norm_common(F);
			K[17] += F[2]*f;
			K[24+(-B[i].link[NUM_PHASE]-1)] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...дополнительное интегрирование параметров потока для блоков с усложненными функциями;
int CLame3D::rigidy5(CGrid * nd, int i, double * K)
{
   int	N_elem = UnPackInts(get_param(NUM_QUAD)), N_max = UnPackInts(get_param(NUM_QUAD), 1), m, k, l, l0;
   double 	 R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1), f, P[6], U[3], sum[25],
			alphaM = get_param(NUM_SHEAR+1)/(.5-get_param(NUM_SHEAR+1)), GM = get_param(NUM_SHEAR), 
			alphaI = get_param(NUM_SHEAR+NUM_SHIFT+1)/(.5-get_param(NUM_SHEAR+NUM_SHIFT+1)), GI = get_param(NUM_SHEAR+NUM_SHIFT), 
			alphaL = get_param(NUM_SHEAR+NUM_SHIFT*2+1)/(.5-get_param(NUM_SHEAR+NUM_SHIFT*2+1)), GL = get_param(NUM_SHEAR+NUM_SHIFT*2);
	
///////////////////////////////////////////////////////////////////////
//...коррекция интегральных характеристик в случае композитной функции;
	if ((B[i].type & ERR_CODE) == CLAYER_BLOCK && B[i].mp[0] == ID_MAP(2, SPHEROID_GENUS) && (R1 || R2))	{
		CGrid * bound_bnd = CreateNodes(GRID_QG_NODES);
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

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			bound_bnd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_X, 1);

/////////////////////////////
//...интегрируем перемещения;
		for (memset(sum, 0, 25*sizeof(double)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l]; f = bound_bnd->get_param(0, l);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			for (int num = solver.id_norm; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

			B[i].shape->parametrization_hess(P, 1);
			jump1_x(P, i, 0); 
			jump1_y(P, i, 1); 
			jump1_z(P, i, 2); 

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   0);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 0);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 0);
			B[i].shape->norm_common(U);

			sum[0]  += ((GI-GL)* U[0]*bound_bnd->nX[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[1]  += ((GI-GL)*(U[0]*bound_bnd->nY[l]+U[1]*bound_bnd->nX[l]))*f;
			sum[2]  += ((GI-GL)* U[1]*bound_bnd->nY[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[3]  += ((GI-GL)*(U[0]*bound_bnd->nZ[l]+U[2]*bound_bnd->nX[l]))*f;
			sum[4]  += ((GI-GL)*(U[1]*bound_bnd->nZ[l]+U[2]*bound_bnd->nY[l]))*f;
			sum[5]  += ((GI-GL)* U[2]*bound_bnd->nZ[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   1);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 1);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 1);
			B[i].shape->norm_common(U);
			
			sum[12] += ((GI-GL)* U[0]*bound_bnd->nX[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[13] += ((GI-GL)*(U[0]*bound_bnd->nY[l]+U[1]*bound_bnd->nX[l]))*f;
			sum[14] += ((GI-GL)* U[1]*bound_bnd->nY[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[15] += ((GI-GL)*(U[0]*bound_bnd->nZ[l]+U[2]*bound_bnd->nX[l]))*f;
			sum[16] += ((GI-GL)*(U[1]*bound_bnd->nZ[l]+U[2]*bound_bnd->nY[l]))*f;
			sum[17] += ((GI-GL)* U[2]*bound_bnd->nZ[l]*2.+(GI*alphaI-GL*alphaL)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[24] += P[2]*P[5]*f;
		}
		if (sum[24] < 0.)  //...правим знак в случае ошибки в направлении нормали;
			for (m = 0; m < 25; m++) sum[m] -= sum[m];
		K[25] += sum[24];
		K[26] -= sum[24];
		for (m = 0; m < 24; m++) K[m] += sum[m];
		if (! solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

//////////////////////////////////////////////////////
//...повторяем граничные точки для второй поверхности;
		ce->mp[7] = R2;
      ce->mp[l0-1] = 0.;
      ce->mp[l0] = M_PI/N_max;
		
		for (m = 0; m < N_max; m++, ce->mp[l0-1] = ce->mp[l0], ce->mp[l0] *= (m+1.)/m)
		for (k = 0; k < N_max; k++) {
			ce->mp[6] = ce->mp[l0-2]*k;
			bound_bnd->sphere_intrusion_QG(ce->mp, N_elem, K[27]);
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			bound_bnd->TestGrid("nodes.bln", 0.005, 10., 20., 30., AXIS_X, 1);

////////////////////////////////
//...интегрирование перемещений;
		for (memset(sum, 0, 25*sizeof(double)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l]; f = bound_bnd->get_param(0, l);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			for (int num = solver.id_norm; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

			B[i].shape->parametrization_hess(P, 1);
			jump1_x(P, i, 0); 
			jump1_y(P, i, 1); 
			jump1_z(P, i, 2); 

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   0);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 0);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 0);
			B[i].shape->norm_common(U);

			sum[0]  += ((GL-GM)* U[0]*bound_bnd->nX[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[1]  += ((GL-GM)*(U[0]*bound_bnd->nY[l]+U[1]*bound_bnd->nX[l]))*f;
			sum[2]  += ((GL-GM)* U[1]*bound_bnd->nY[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[3]  += ((GL-GM)*(U[0]*bound_bnd->nZ[l]+U[2]*bound_bnd->nX[l]))*f;
			sum[4]  += ((GL-GM)*(U[1]*bound_bnd->nZ[l]+U[2]*bound_bnd->nY[l]))*f;
			sum[5]  += ((GL-GM)* U[2]*bound_bnd->nZ[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;

			U[0] = B[i].shape->potential(solver.hh[i][0][solver.id_norm],   1);
			U[1] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+1], 1);
			U[2] = B[i].shape->potential(solver.hh[i][0][solver.id_norm+2], 1);
			B[i].shape->norm_common(U);
			
			sum[12] += ((GL-GM)* U[0]*bound_bnd->nX[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[13] += ((GL-GM)*(U[0]*bound_bnd->nY[l]+U[1]*bound_bnd->nX[l]))*f;
			sum[14] += ((GL-GM)* U[1]*bound_bnd->nY[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[15] += ((GL-GM)*(U[0]*bound_bnd->nZ[l]+U[2]*bound_bnd->nX[l]))*f;
			sum[16] += ((GL-GM)*(U[1]*bound_bnd->nZ[l]+U[2]*bound_bnd->nY[l]))*f;
			sum[17] += ((GL-GM)* U[2]*bound_bnd->nZ[l]*2.+(GL*alphaL-GM*alphaM)*(U[0]*bound_bnd->nX[l]+U[1]*bound_bnd->nY[l]+U[2]*bound_bnd->nZ[l]))*f;
			sum[24] += P[2]*P[5]*f;
		}
		if (sum[24] < 0.)  //...правим знак в случае ошибки в направлении нормали;
			for (m = 0; m < 25; m++) sum[m] -= sum[m];
		if ( K[24] < R2) //...правим площадь промежуточного слоя в случае взаимопроникающих включений;
			sum[24] += (R2*R2-K[24]*K[24])*M_PI*2.*K[24];

		K[26] += sum[24];
		K[24] -= sum[24];
		for (m = 0; m < 24; m++)  K[m] += sum[m];
		if (! solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

		delete ce;
		delete bound_bnd;
	}
	return(OK_STATE);
}

/////////////////////////////////////////////////////
//...installing parameters of doubly connected media;
void CLame3D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double nju3, double G3)
{
  if (size_of_param() > NUM_SHEAR+2+NUM_SHIFT*2) {
      param[NUM_SHEAR] = G1;
      param[NUM_SHEAR+NUM_SHIFT] = G2;
      param[NUM_SHEAR+NUM_SHIFT*2] = G3;
      param[NUM_SHEAR+1] = nju1;
      param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
      param[NUM_SHEAR+1+NUM_SHIFT*2] = nju3;
      param[NUM_SHEAR+2] = .25/(nju1-1.);
      param[NUM_SHEAR+2+NUM_SHIFT] = .25/(nju2-1.);
      param[NUM_SHEAR+2+NUM_SHIFT*2] = .25/(nju3-1.);
  }
  return;
}

void CLame3D::set_fasa_hmg(double R1, double R2, double nju1, double nju2, double nju3, double G1, double G2, double G3, double alpha)
{
	if (size_of_param() > NUM_GEOMT+2) {
		param[NUM_GEOMT] = R1;
		param[NUM_GEOMT+1] = R2;
		param[NUM_GEOMT+2] = alpha;
	}
	set_fasa_hmg(nju3, nju1, G3, G1, nju2, G2);
}

///////////////////////////////////////////////////////////
//...counting header for solving linear elasticity problem;
int CLame3D::computing_header(Num_Comput Num)
{
	int  N_elem = UnPackInts(get_param(NUM_QUAD)), id_dir, i, k, elem, n_rhs = 6;
	char msg[201];
	
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CLame3D sample: N_sm = %d, N_mpl = %d, N_elem = %d, Q_facet = %g", N, UnPackInts(get_param(NUM_MPLS)), N_elem, get_param(NUM_QUAD+1));
		Message(msg);

		Message(" ");
		switch (Num){
			case	 BASIC_COMPUT: Message("Junction Blocks..."); break;
			case MAPPING_COMPUT: Message("Mapping  Blocks..."); break;
			case  PERIOD_COMPUT: Message("Periodic Blocks..."); break;
			case ESHELBY_COMPUT: Message("Eshelby Problem..."); break;
		}
		Message(" ");
	}
	if (N == 1 && (B[0].type & ERR_CODE) == CLAYER_BLOCK)
		set_eshelby_matrix (UnPackInts(get_param(NUM_MPLS)));

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += NUM_HESS+6;//<==== number of additional auxilliary arrays!!!
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

///////////////////////////////////////////////////////////////////////////
//...расчет функций формы и матрицы жесткости с помощью нового функционала;
extern "C" {//...lapack procedure;
int dgesvd_ (char * jobu, char * jobvt, int * m, int * n, 
				double * a,  int * lda,  double * s, double * u, int * ldu, 
				double * vt, int * ldvt, double * work, int * lwork, int * info);
};
int CLame3D::rigidy1_collocat(CGrid * nd, int i, double * K)
{
	if (nd && 0 <= i && i < solver.N) {
		double G1 = get_param(NUM_SHEAR), P[6] = {0., 0., 0., 0., 0., 1.}, Po[12], F[3], h1, h2, h3, h4, f;
		int 	 m  = solver.id_norm, m1, m2, m3, m4, j, k, l, shift, * mask = NULL, N_elem = 10;

//////////////////////////////////////////
//...формируем матрицу узловых коллокаций;
		int id_collocat = 0;
		if (id_collocat) {
			m1 = m/3; k = shift = 0;
			while (shift < m1-nd->N) {
				for ( ; k < nd->N+shift; k++) {
					P[0] = nd->X[k-shift];				
					P[1] = nd->Y[k-shift];				
					P[2] = nd->Z[k-shift];		
					B[i].shape->make_local(P);
					B[i].shape->parametrization_hess(P, 1);

					jump1_classic_x(P, i, 0);
					jump1_classic_y(P, i, 1);
					jump1_classic_z(P, i, 2);
					jump_make_common(i, 0); //...make_common displacement {Ux, Uy, Uz};

					solver.to_equationEH(i, k*3,	 solver.hh[i][0][m],   1.);
					solver.to_equationEH(i, k*3+1, solver.hh[i][0][m+1], 1.);
					solver.to_equationEH(i, k*3+2, solver.hh[i][0][m+2], 1.);
				}
				shift += nd->N;
			}
			for ( ; k < m1; k++) {
				P[0] = nd->X[k-shift]; P[3] = nd->nX[k-shift];				
				P[1] = nd->Y[k-shift]; P[4] = nd->nY[k-shift];				
				P[2] = nd->Z[k-shift]; P[5] = nd->nZ[k-shift];		
				B[i].shape->make_local(P);
				B[i].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, i, 0);
				jump1_classic_y(P, i, 1);
				jump1_classic_z(P, i, 2);
				jump_make_common(i, 0); //...make_common displacement {Ux, Uy, Uz};

				solver.to_equationEH(i, k*3,	 solver.hh[i][0][m],	  1.);
				solver.to_equationEH(i, k*3+1, solver.hh[i][0][m+1], 1.);
				solver.to_equationEH(i, k*3+2, solver.hh[i][0][m+2], 1.);
			}

/////////////////////////
//...SVD-decomposition!!!
			int id_action = 0;
			if (id_action)
			memcpy(K, solver.hh[i][2][0], solver.dim[i]*solver.dim[i]*sizeof(double));
#ifdef ___MPI_INIT___
			dgesvd_("A", "O", & solver.dim[i], & solver.dim[i], solver.hh[i][2][0],
									& solver.dim[i], solver.hh[i][0][m], solver.hh[i][0][0],
									& solver.dim[i], NULL, & solver.dim[i], solver.hh[i][0][m+1],
									& (m2 = (solver.n-solver.id_norm)*solver.dim[i]), & (k = 0));
			for (j = nd->N*3; j < m; j++) solver.hh[i][0][m][j] = 0.;
			for (j = min(m, nd->N*3)-1; j >= 0; j--) solver.hh[i][0][m][j] = fabs(solver.hh[i][0][m][j]) > EE ? 1./solver.hh[i][0][m][j] : 0.;
#endif
////////////////////////////////////////////
//...формируем коллокационные функции формы;
			for (k = 0; k < solver.dim[i]; k++) {
				for (l = 0; l < solver.dim[i]; l++)
				for (	 solver.hh[i][0][m+1][l] = 0., j = 0; j < solver.dim[i]; j++) //...VT*L^(-1)*U
						 solver.hh[i][0][m+1][l] += solver.hh[i][2][k][j]*solver.hh[i][0][m][j]*solver.hh[i][0][j][l];
				memcpy(solver.hh[i][2][k], solver.hh[i][0][m+1], solver.dim[i]*sizeof(double)); 
			}

//////////////////////////////////////////////////////////////
//...склеиваем коллокационные функции формы для кратных узлов;
			k = shift = nd->N;
			while (shift < m1-nd->N) {
				for ( ; k < nd->N+shift; k++) {
					solver.admittance(solver.hh[i][2][(k-shift)*3],	solver.hh[i][2][k*3],	 solver.dim[i], 1., 1.);
					solver.admittance(solver.hh[i][2][(k-shift)*3+1], solver.hh[i][2][k*3+1], solver.dim[i], 1., 1.);
					solver.admittance(solver.hh[i][2][(k-shift)*3+2], solver.hh[i][2][k*3+2], solver.dim[i], 1., 1.);
				}
				shift += nd->N;
			}
			for ( ; k < m1; k++) {
				solver.admittance(solver.hh[i][2][(k-shift)*3],	solver.hh[i][2][k*3],	 solver.dim[i], 1., 1.);
				solver.admittance(solver.hh[i][2][(k-shift)*3+1], solver.hh[i][2][k*3+1], solver.dim[i], 1., 1.);
				solver.admittance(solver.hh[i][2][(k-shift)*3+2], solver.hh[i][2][k*3+2], solver.dim[i], 1., 1.);
			}

////////////////////////////////////////////////////////////
//...тестовое формирование и печать коллокационных значений;
			if (id_action) {
				for (k = 0; k < solver.dim[i]; k++) {
					for (l = 0; l < solver.dim[i]; l++) 
					for (	 solver.hh[i][0][m+1][l] = 0., j = 0; j < solver.dim[i]; j++)
							 solver.hh[i][0][m+1][l] += solver.hh[i][2][k][j]*K[l*solver.dim[i]+j];
					memcpy(solver.hh[i][2][k], solver.hh[i][0][m+1], solver.dim[i]*sizeof(double)); 
				}
				
				FILE  * TST = fopen("StiffnessMatrix.sta", "w");
				fprintf(TST, "\nCollocation values (NN = %i, NDP = %i):\n", m3 = min(m, nd->N*3), nd->N);
				for (l = 0; l < m3; l++) {
					fprintf(TST, "\n");
					for (k = 0; k < m3; k++) 
					fprintf(TST, "%g ", filtr_(solver.hh[i][2][k][l], EE_dop));
				}  fprintf(TST, "\n");
				fclose(TST);

				return(1);
			}
		}

///////////////////////////////////////////////////////
//...строим функции формы методом наименьших квадратов;
      if (1 && nd->geom && solver.id_norm >= nd->N*3) {

//////////////////////////////////
//...строим вспомогательную маску;
			for (j = k = 0; k < nd->N; k++) j = max(j, nd->hit[k]);
			if ((mask = (int *)new_struct((j+1)*sizeof(int))) != NULL)
			for (k = 0; k < nd->N; k++) mask[nd->hit[k]] = k;

			CGrid_el * bnd = (CGrid_el * )CreateNodes(GRID_EL_NODES);
						  bnd->add_params(1);

///////////////////////////////////////////////////////////////
//...формируем систему уравнений для определения функций формы;
			for (l = 0, k = 0; k < nd->geom[0]; k++) {
				switch (nd->geom[(++l)++]) {
					case GL_TRIANGLES: if (nd->geom[l] == 5) {
						j   = -nd->geom[l+2];
						m1  = mask[nd->geom[l+3]];
						m2  = mask[nd->geom[l+4]];
						m3  = mask[nd->geom[l+5]];

////////////////
//...tria facet;
						Po[0] = nd->X[m1]; Po[1]  = nd->Y[m1]; Po[2]  = nd->Z[m1];
						Po[3] = nd->X[m2]; Po[4]  = nd->Y[m2]; Po[5]  = nd->Z[m2];
						Po[6] = nd->X[m3]; Po[7]  = nd->Y[m3]; Po[8]  = nd->Z[m3];
					}  break;
					case GL_QUADS: if (nd->geom[l] == 6) {
						j   = -nd->geom[l+2]; //...свойство указывает на криволинейную поверхность!!!
						m1  = mask[nd->geom[l+3]];
						m2  = mask[nd->geom[l+4]];
						m3  = mask[nd->geom[l+5]];
						m4  = mask[nd->geom[l+6]];

////////////////
//...quad facet;
						Po[0] = nd->X[m1]; Po[1]  = nd->Y[m1]; Po[2]  = nd->Z[m1];
						Po[3] = nd->X[m4]; Po[4]  = nd->Y[m4]; Po[5]  = nd->Z[m4];
						Po[6] = nd->X[m2]; Po[7]  = nd->Y[m2]; Po[8]  = nd->Z[m2];
						Po[9] = nd->X[m3]; Po[10] = nd->Y[m3]; Po[11] = nd->Z[m3];

						bnd->facet_QG(Po, N_elem, OK_STATE, NULL_STATE);
						bnd->QG_quad_bi(Po);
						if (bar && bar->graph && 0 < j && j <= bar->graph[0])
						bnd->QG_surface(bar->ce[j-1]->mp);

						for (int lp = 0; lp < bnd->N; lp++) {
							P[0] = bnd->X[lp]; P[3] = bnd->nX[lp];
							P[1] = bnd->Y[lp]; P[4] = bnd->nY[lp];
							P[2] = bnd->Z[lp]; P[5] = bnd->nZ[lp];
							h1 = (1.-bnd->get_param(1, lp))*(1.-bnd->get_param(2, lp));
							h2 = bnd->get_param(1, lp)*(1.-bnd->get_param(2, lp));
							h3 = bnd->get_param(1, lp)*bnd->get_param(2, lp);
							h4 = (1.-bnd->get_param(1, lp))*bnd->get_param(2, lp);
							f  = bnd->get_param(0, lp);

							for (int num = m; num < solver.n; num++)
								memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

							B[i].shape->make_local(P);
							B[i].shape->norm_local(P+3);
							B[i].shape->parametrization_hess(P, 1);

							jump1_classic_x(P, i, 0);
							jump1_classic_y(P, i, 1);
							jump1_classic_z(P, i, 2);

							jump4_classic_x(P, i, 3);
							jump4_classic_y(P, i, 4);
							jump4_classic_z(P, i, 5);

							solver.admittance(i, 0, G1);
							solver.admittance(i, 1, G1);
							solver.admittance(i, 2, G1);

							solver.to_equationDD(i, solver.hh[i][0][m], solver.hh[i][0][m], f);
							solver.to_equationHH(i, m1*3, solver.hh[i][0][m], h1*f);
							solver.to_equationHH(i, m2*3, solver.hh[i][0][m], h2*f);
							solver.to_equationHH(i, m3*3, solver.hh[i][0][m], h3*f);
							solver.to_equationHH(i, m4*3, solver.hh[i][0][m], h4*f);

							solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
							solver.to_equationHH(i, m1*3+1, solver.hh[i][0][m+1], h1*f);
							solver.to_equationHH(i, m2*3+1, solver.hh[i][0][m+1], h2*f);
							solver.to_equationHH(i, m3*3+1, solver.hh[i][0][m+1], h3*f);
							solver.to_equationHH(i, m4*3+1, solver.hh[i][0][m+1], h4*f);

							solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
							solver.to_equationHH(i, m1*3+2, solver.hh[i][0][m+2], h1*f);
							solver.to_equationHH(i, m2*3+2, solver.hh[i][0][m+2], h2*f);
							solver.to_equationHH(i, m3*3+2, solver.hh[i][0][m+2], h3*f);
							solver.to_equationHH(i, m4*3+2, solver.hh[i][0][m+2], h4*f);

//////////////////////////////
//...энергетические слагаемые;
							solver.to_equationEE(i, solver.hh[i][0][m],   solver.hh[i][0][m+3], f);
							solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4], f);
							solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5], f);
						}
						bnd->add_buffer(bnd->N);
					}  break;
				}
				l += nd->geom[l];
			}
			delete_struct(mask);
			delete bnd;

//////////////////////////////////////////
//...test output transfer and gram matrix;
			FILE * TST_gram = NULL,//fopen("GramMatrix.sta", "w"),//
				  * TST_trans = NULL;//fopen("TransMatrix.sta", "w");//
			if (TST_gram || TST_trans) {

				Message("Block matrix output...");

				for (k = 0; k < N; k++) {
					solver.test_transfer_matrix(TST_trans, k);
					solver.test_gram_matrix(TST_gram, k, 6);
				}
			}
			if (TST_gram) fclose(TST_gram);
			if (TST_trans) fclose(TST_trans);

///////////////////////
//...Blocked GaussJ !!!
			Message("Blocked GaussJ...");

			double f_count = 0., g_count, dim_N = 2*solver.N;
			char buff[1000];

			//solver.lagrangian(get_param(size_of_param()-1));
			for (k = 0; k < solver.N; k++) {
				if (! solver.solver(k)) {
					return ERR_STATE;
				}
				if ((g_count = (int)(100.*(k+1.)/(double)dim_N+.5)/5) > f_count) {
					f_count = g_count;
					sprintf(buff, "Solver : %d", (int)(f_count*5.));
					Message(buff);
				}
			}
			for (k = solver.N-1; k >= 0; k--) {
				if (! solver.solver(-k-1)) {
					return ERR_STATE;
				}
				if ((g_count = (int)(100.*(2*solver.N-k)/(double)dim_N+.5)/5) > f_count) {
					f_count = g_count;
					sprintf(buff, "Solver : %d", (int)(f_count*5.));
					Message(buff);
				}
			}
			solver.reset_struct();
			shapes_init(OK_STATE);
		}

//////////////////////////////////////////////////////////////////////////////
//...правим узловые значения функций формы при помощи коллокационных значений;
		for (m3 = min(m/3, nd->N), k = 0; k < m3; k++) {
			P[0] = nd->X[k];
			P[1] = nd->Y[k];
			P[2] = nd->Z[k];

			solver.admittance(solver.hh[i][0][k*3],	  solver.hh[i][2][k*3], solver.dim[i], 0., 1./*, 0.*/); //...??? -- было 0. на hEE!!!
			solver.admittance(solver.hh[i][0][k*3+1], solver.hh[i][2][k*3+1], solver.dim[i], 0., 1./*, 0.*/); //...???
			solver.admittance(solver.hh[i][0][k*3+2], solver.hh[i][2][k*3+2], solver.dim[i], 0., 1./*, 0.*/); //...???

			B[i].shape->make_local(P);
			B[i].shape->parametrization_hess(P, 1);

			jump1_classic_x(P, i, 0); 
			jump1_classic_y(P, i, 1); 
			jump1_classic_z(P, i, 2); 
//			jump_make_common(i, 0); //...make_common displacement {Ux, Uy, Uz};

/////////////////////////////////////////////////
//...тестовое вычисление коллокационных значений;
			int id_test_collocat = 0;
			if (id_test_collocat) {
				shapes_init(OK_STATE);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],                    k*9+0)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN,     k*9+1)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN*2,   k*9+2);

				F[1] = B[i].shape->potential(solver.hh[i][0][m+1],                  k*9+0)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN,   k*9+1)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN*2, k*9+2);

				F[2] = B[i].shape->potential(solver.hh[i][0][m+2],                  k*9+0)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN,   k*9+1)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN*2, k*9+2);
				B[i].shape->norm_common(F);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],                    k*9+3)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN,     k*9+4)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN*2,   k*9+5);

				F[1] = B[i].shape->potential(solver.hh[i][0][m+1],                  k*9+3)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN,   k*9+4)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN*2, k*9+5);

				F[2] = B[i].shape->potential(solver.hh[i][0][m+2],                  k*9+3)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN,   k*9+4)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN*2, k*9+5);
				B[i].shape->norm_common(F);

				F[0] = B[i].shape->potential(solver.hh[i][0][m],                    k*9+6)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN,     k*9+7)+
						 B[i].shape->potential(solver.hh[i][0][m]+B[i].shape->NN*2,   k*9+8);

				F[1] = B[i].shape->potential(solver.hh[i][0][m+1],                  k*9+6)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN,   k*9+7)+
						 B[i].shape->potential(solver.hh[i][0][m+1]+B[i].shape->NN*2, k*9+8);

				F[2] = B[i].shape->potential(solver.hh[i][0][m+2],                  k*9+6)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN,   k*9+7)+
						 B[i].shape->potential(solver.hh[i][0][m+2]+B[i].shape->NN*2, k*9+8);
				B[i].shape->norm_common(F);
			}
		}
		shapes_init(OK_STATE);

		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////////////////////////
//...формирование матрицы переразложения для оператора vec{r}grad{vec{f}_0};
void CLame3D::rgradf_matrix(double ** T, int n, int shift_m, int shift_n, double f)
{
	if (T && n > 0) {
		int nn = 2*n+1, m, shifn, shifm;
		double dd; f /= (nn-2.);
//...Re fw;
		T[shift_m][shift_n] += (dd = (n-1.)*n*.25*f); 
		T[shift_m+nn][shift_n+nn] -= dd;
		if (n >= 2) {
			T[shift_m+nn][shift_n+3] += f; 
			T[shift_m][shift_n+4] -= f; T[shift_m][3+nn+shift_n] -= f; T[shift_m+nn][shift_n+4+nn] -= f;
		}
		if (n >= 1) {
			T[shift_m+nn][shift_n+1+nn*2] += (dd = -n*.5*f); 
			T[shift_m][shift_n+2+nn*2] -= dd;
		}
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2;
			T[shifm][shifn] += (dd = (n+m)*(n+m-1.)*.125*f); 
			T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; T[shifm-1][shifn+nn] -= dd;
			if (m+2 <= n) {
				T[shifm-1][shifn+3] += (dd = (m+1.)*(m+2.)*.5*f); 
				T[shifm][shifn+4] -= dd; T[shifm][shifn+3+nn] -= dd; T[shifm-1][shifn+4+nn] -= dd;
			}
			if (m+1 <= n) {
				T[shifm-1][shifn+1+nn*2] += (dd = -(m+1.)*(n+m)*.5*f); 
				T[shifm][shifn+2+nn*2] -= dd;
			}
		}
//...Im fw;
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn;
			T[shifm-1][shifn-1] += (dd = -(n-m-1.)*(n-m)*.125*f); 
			T[shifm][shifn] -= dd; T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn+nn] -= dd;
			if (m == 2) {
				T[shifm-1][shift_n+nn] += (dd = (n-1.)*n*(n+1.)*(n+2.)*.03125*f); 
				T[shifm][shift_n] -= dd;
			}
			if (m == 1)	{	
				T[shifm-1][shifn+nn] += (dd = -((n-2.)*(n-1.)*.5+(nn-2.))*.25*f);
				T[shifm][shifn] -= dd; T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; 
				T[shifm][shift_n+nn*2] += (1.-n*n)*n*.25*f;
			}
			if (m > 2) {
				T[shifm][shifn-4] += (dd = -(n-m+1.)*(n-m+2.)*(n+m-1.)*(n+m)/((m-1.)*m)*.03125*f);
				T[shifm][shifn-5+nn] -= dd; T[shifm-1][shifn-5] -= dd; T[shifm-1][shifn-4+nn] -= dd;
			}
			if (m > 1) {
				T[shifm-1][shifn-3+nn*2] += (dd = (n-m)*(n-m+1.)*(n+m)/m*.125*f); 
				T[shifm][shifn-2+nn*2] -= dd;
			}
		}
//...fz;
		T[shift_m+nn*2][shift_n+nn*2] += n*n*f;
		if (n >= 1) {
			T[shift_m+nn*2][shift_n+1+nn] += (dd = -(n-1.)*f);
			T[shift_m+nn*2][shift_n+2] += dd;
		}
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn*2;
			T[shifm-1][shifn-1+nn*2] += (dd = -(n-m)*(n+m)*.5*f);
			T[shifm][shifn+nn*2] -= dd;
			if (m == 1)	{
				T[shifm-1][shift_n+nn] += (dd = -n*n*(n+1.)*.25*f);
				T[shifm][shift_n] -= dd;
			}
			if (n >= 1+m) {
				T[shifm-1][shifn+1] += (dd = (m+1.)*(n-m-1.)*.5*f);
				T[shifm][shifn+2] -= dd; T[shifm][shifn+1+nn] -= dd; T[shifm-1][shifn+2+nn] -= dd; 
			}
			if (m > 1) {
				T[shifm][shifn-2] += (dd = (n-m+1.)*(n+m-1.)*(n+m)/m*.125*f);
				T[shifm][shifn-3+nn] -= dd; T[shifm-1][shifn-3] -= dd; T[shifm-1][shifn-2+nn] -= dd;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////
//...формирование матрицы переразложения для оператора vec{r}div{vec{f}_0};
void CLame3D::rdivrf_matrix(double ** T, int n, int shift_m, int shift_n, double f)
{
	if (T && n > 0) {
		int nn = 2*n+1, m, shifn, shifm;
		double dd; f /= (nn-2.);
//...Re fw;
		T[shift_m][shift_n] += (dd = (n-1.)*n*.25*f);
		T[shift_m+nn][shift_n+nn] -= dd;
		if (n >= 2) {
			T[shift_m+nn][shift_n+3] += f;
			T[shift_m][shift_n+4] -= f; T[shift_m][3+nn+shift_n] -= f; T[shift_m+nn][shift_n+4+nn] -= f;
		}
		if (n >= 1) {
			T[shift_m+nn][shift_n+1+nn*2] += (dd = (n-1.)*.5*f);
			T[shift_m][shift_n+2+nn*2] -= dd;
		}
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2;
			T[shifm][shifn] += (dd = (n-m-1.)*(n-m)*.125*f);
			T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; T[shifm-1][shifn+nn] -= dd;
			if (m+2 <= n) {
				T[shifm-1][shifn+3] += (dd = (m+1.)*(m+2.)*.5*f);
				T[shifm][shifn+4] -= dd; T[shifm][shifn+3+nn] -= dd; T[shifm-1][shifn+4+nn] -= dd;
			}
			if (m+1 <= n) {
				T[shifm-1][shifn+1+nn*2] += (dd = (m+1.)*(n-m-1.)*.5*f);
				T[shifm][shifn+2+nn*2] -= dd;
			}
		}
//...Im fw;
		for (m = 1; m <= n; m++) { 
			double d = (n+m-1.)*(n+m)*.125*f;
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn;
			T[shifm-1][shifn-1] -= d;
			T[shifm][shifn] += d; T[shifm][shifn-1+nn] += d; T[shifm-1][shifn+nn] += d;
			if (m == 2) {
				T[shifm-1][shift_n+nn] += (dd = (n-1.)*n*.25*d);
				T[shifm][shift_n] -= dd;
			}
			if (m == 1)	{	
				T[shifm-1][shifn+nn] -= d;
				T[shifm][shifn] += d; T[shifm][shifn-1+nn] += d; T[shifm-1][shifn-1] += d; 
				T[shifm][shift_n+nn*2] += 2.*n*d;
			}
			if (m > 2) {
				T[shifm][shifn-4] += (dd = -(n-m+1.)*(n-m+2.)/((m-1.)*m)*.25*d);
				T[shifm][shifn-5+nn] -= dd; T[shifm-1][shifn-5] -= dd; T[shifm-1][shifn-4+nn] -= dd;
			}
			if (m > 1) {
				T[shifm-1][shifn-3+nn*2] += (dd = -(n-m+1.)/m*d);
				T[shifm][shifn-2+nn*2] -= dd;
			}
		}
//...fz;
		T[shift_m+nn*2][shift_n+nn*2] += n*n*f;
		if (n >= 1) {
			T[shift_m+nn*2][shift_n+1+nn] += (dd = n*f);
			T[shift_m+nn*2][shift_n+2] += dd;
		}
		for (m = 1; m <= n; m++) { 
			double d = (n+m)*.5*f;
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn*2;
			T[shifm-1][shifn-1+nn*2] += (dd = -(n-m)*d);
			T[shifm][shifn+nn*2] -= dd;
			if (m == 1)	{
				T[shifm-1][shift_n+nn] += (dd = (n-1.)*n*.5*d);
				T[shifm][shift_n] -= dd;
			}
			if (n >= 1+m) {
				T[shifm-1][shifn+1] += (dd = -(m+1.)*d);
				T[shifm][shifn+2] -= dd; T[shifm][shifn+1+nn] -= dd; T[shifm-1][shifn+2+nn] -= dd;
			}
			if (m > 1) {
				T[shifm][shifn-2] += (dd = -(n-m)*(n-m+1.)/m*.25*d);
				T[shifm][shifn-3+nn] -= dd; T[shifm-1][shifn-3] -= dd; T[shifm-1][shifn-2+nn] -= dd;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////
//...формирование матрицы переразложения для тождественного оператора;
void CLame3D::unit_f_matrix(double ** T, int n, int shift_m, int shift_n, double f)
{
	if (T) {
		int nn = 2*n+1, m, shifn, shifm;
		double dd;
//...Re fw;
		T[shift_m][shift_n] += (dd = .5*f);
		T[shift_m+nn][shift_n+nn] -= dd;
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2;
			T[shifm][shifn] += (dd = .25*f);
			T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; T[shifm-1][shifn+nn] -= dd;
		}
//...Im fw;
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn;
			T[shifm-1][shifn-1] -= (dd = .25*f);
			T[shifm][shifn] += dd; T[shifm][shifn-1+nn] += dd; T[shifm-1][shifn+nn] += dd;		
		}
//...fz;
		T[shift_m+nn*2][shift_n+nn*2] += f;
		for (m = 1; m <= n; m++) { 
			shifn = shift_n+m*2; shifm = shift_m+m*2+nn*2;
			T[shifm-1][shifn-1+nn*2] -= (dd = .5*f);
			T[shifm][shifn+nn*2] += dd;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
//...раскомплексификация строчек матрицы переразложения граничного оператора;
void CLame3D::toreal_matrix(double ** T, int n, int shift_m, int shift_n)
{
	if (T) {
		int l, m, nn = 2*n+1, shifm;
		for (l = 3*nn-1+shift_n; l >= shift_n; l--) {				
			T[shift_m][l] *= 2.; 
			T[shift_m+nn][l] *= -2.;
			for (m = 1; m <= n; m++) {
				shifm = shift_m+m*2;
				T[shifm][l] = 2.*(T[shifm+nn][l]+T[shifm][l]); 
				T[shifm-1][l] = -2.*(T[shifm-1+nn][l]+T[shifm-1][l]);
				T[shifm+nn][l] = 4.*T[shifm+nn][l]-T[shifm][l];
				T[shifm-1+nn][l] = 4.*T[shifm-1+nn][l]+T[shifm-1][l];
				T[shifm+nn*2][l] *=  2.; 
				T[shifm-1+nn*2][l] *= -2.; 
				swap(T[shifm+nn][l], T[shifm-1+nn][l]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы перехода для оператора grad{div{vec{f}_0}};
void CLame3D::grdivf_transf(double ** T, int n, int shift_m, int shift_n, double f)
{ 
	if (T && n > 1) {
		int nn = 2*n+1, m, mm = nn-4, shifn, shifm; 
		double dd;
//...Re fw;
		T[shift_m+mm][shift_n+nn] += (dd = (n-1.)*n*.25*f);
		T[shift_m][shift_n] -= dd;
		T[shift_m+mm][shift_n+3] -= f;
		T[shift_m][shift_n+4] += f; T[shift_m][shift_n+3+nn] += f; T[shift_m+mm][shift_n+4+nn] += f;
		T[shift_m+mm][shift_n+1+nn*2] += (dd = -(n-1.)*.5*f);
		T[shift_m][shift_n+2+nn*2] -= dd;
		for (m = 1; m <= n-2; m++) { 
			shifn = shift_n+m*2;	shifm = shift_m+m*2;
			T[shifm][shifn] += (dd = -(n-m-1.)*(n-m)*.125*f);
			T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; T[shifm-1][shifn+nn] -= dd;
			T[shifm-1][shifn+3] += (dd = -(m+1.)*(m+2.)*.5*f);
			T[shifm][shifn+4] -= dd; T[shifm][shifn+3+nn] -= dd; T[shifm-1][shifn+4+nn] -= dd;
			T[shifm-1][shifn+1+nn*2] += (dd = -(m+1.)*(n-m-1.)*.5*f);
			T[shifm][shifn+2+nn*2] -= dd;
		}
//...Im fw;
		for (m = 1; m <= n-2; m++) { 
			shifn = shift_n+m*2;	shifm = shift_m+m*2+mm;
			T[shifm-1][shifn-1] += (dd = (n-m-1.)*(n-m)*.125*f);
			T[shifm][shifn] -= dd; T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn+nn] -= dd;
			if (m == 2) {
				T[shifm-1][shift_n+nn] += (dd = -(n-3.)*(n-2.)*(n-1.)*n*.03125*f);
				T[shifm][shift_n] -= dd;
			}
			if (m == 1)	{	
				T[shifm-1][shifn+nn] += (dd = (n-2.)*(n-1.)*.125*f);
				T[shifm][shifn] -= dd; T[shifm][shifn-1+nn] -= dd; T[shifm-1][shifn-1] -= dd; 
				T[shifm][shift_n+nn*2] -= 2.*n*dd;
			}
			if (m > 2) {
				T[shifm][shifn-4] += (dd = (n-m-1.)*(n-m)*(n-m+1.)*(n-m+2.)/((m-1.)*m)*.03125*f);
				T[shifm][shifn-5+nn] -= dd; T[shifm-1][shifn-4+nn] -= dd; T[shifm-1][shifn-5] -= dd;
			}
			if (m > 1) {
				T[shifm-1][shifn-3+nn*2] += (dd = (n-m-1.)*(n-m)*(n-m+1.)/m*.125*f);
				T[shifm][shifn-2+nn*2] -= dd;
			}
		}
//...fz;
		T[shift_m+mm*2][shift_n+nn*2] += (n-1.)*n*f;
		T[shift_m+mm*2][shift_n+1+nn] += (dd = (n-1.)*f);
		T[shift_m+mm*2][shift_n+2] += dd;
		for (m = 1; m <= n-2; m++) { 
			shifn = shift_n+m*2;	shifm = shift_m+m*2+mm*2;
			T[shifm-1][shifn-1+nn*2] += (dd = -(n-m-1.)*(n-m)*.5*f);
			T[shifm][shifn+nn*2] -= dd;
			if (m == 1)	{
				T[shifm-1][shift_n+nn] += (dd = (n-2.)*(n-1.)*n*.25*f);
				T[shifm][shift_n] -= dd;
			}
			T[shifm-1][shifn+1] += (dd = -(m+1.)*(n-m-1.)*.5*f);
			T[shifm][shifn+2] -= dd; T[shifm][shifn+1+nn] -= dd; T[shifm-1][shifn+2+nn] -= dd;
			if (m > 1) {
				T[shifm][shifn-2] += (dd = -(n-m-1.)*(n-m)*(n-m+1.)/m*.125*f);
				T[shifm][shifn-3+nn] -= dd; T[shifm-1][shifn-3] -= dd; T[shifm-1][shifn-2+nn] -= dd;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////
//...раскомплексификация строчек матрицы перехода граничного оператора;
void CLame3D::toreal_transf(double ** T, int n, int shift_m, int shift_n)
{
	if (T && n > 1) {
		int l, m, nn = 2*n+1, mm = nn-4, shifm; 
		for (l = 3*nn-1+shift_n; l >= shift_n; l--) {				
			T[shift_m][l] *= 2.; 
			T[shift_m+mm][l] *= -2.;
			for (m = 1; m <= n-2; m++) {
				shifm = shift_m+m*2;
				T[shifm][l] = 2.*(T[shifm+mm][l]+T[shifm][l]); 
				T[shifm-1][l] = -2.*(T[shifm-1+mm][l]+T[shifm-1][l]);
				T[shifm+mm][l] = 4.*T[shifm+mm][l]-T[shifm][l];
				T[shifm-1+mm][l] = 4.*T[shifm-1+mm][l]+T[shifm-1][l];
				T[shifm+mm*2][l] *=  2.; 
				T[shifm-1+mm*2][l] *= -2.; 
				swap(T[shifm+mm][l], T[shifm-1+mm][l]);
			}
		}
	}
}

/////////////////////////////////////////////////////
//...формироване матриц Эшелби для граничных функций;
void CLame3D::set_eshelby_matrix(int N_mpl)
{ 
	int num_bound = 4;
	double ** TX = NULL, ** TD = NULL, ** TD_dop = NULL, * H = (double *)new_struct((6*N_mpl+3)*num_bound*sizeof(double));
	TT = (double ***)new_struct((N_mpl+2)*sizeof(double **));
	TH = (double ***)new_struct((N_mpl+1)*sizeof(double **));
	C0 = (double *)new_struct((N_mpl+1)*sizeof(double));
	C1 = (double *)new_struct((N_mpl+1)*sizeof(double));
	C2 = (double *)new_struct((N_mpl+1)*sizeof(double));
	A1 = (double *)new_struct((N_mpl+1)*sizeof(double));
	B1 = (double *)new_struct((N_mpl+1)*sizeof(double));
	A2 = (double *)new_struct((N_mpl+1)*sizeof(double));
	B2 = (double *)new_struct((N_mpl+1)*sizeof(double));

/////////////////////////////////////////////// 
//...формирование решения граничного уравнения;
	if (TT && TH && H) {
		solver.pivot_init((6*N_mpl+3)*num_bound);
		int  k, j, l, num, nn, mm, mn, nm;
		for (k = N_mpl; k >= 0; k--) { 
			for (mm = (nn = 6*k+3)*num_bound, mn = (nm = nn+12)*(num_bound/2), num = 0; num < nm && TH[k]; num++) {
				for (l = 0; l < mm; l++) 
				for (j = 0; j < mn; j++) TH[k][l][num] -= TD[l][j]*TT[k+2][j][num]; 
			}
			eshelby_matrix_3body(TX, TT[k], k, TD, TH[k-2], num_bound); swap(TD, TD_dop);
			if (solver.GaussJ(TX, mm)) {
				for (num = 0; num < nn && TT[k]; num++) {
					for (           l = 0; l < mm; l++) 
					for (H[l] = 0., j = 0; j < mm; j++) H[l] += TX[l][j]*TT[k][j][num]; 
					for (           l = 0; l < mm; l++) TT[k][l][num] = H[l];
				}
				for (num = 0; num < nm && TH[k]; num++) {
					for (           l = 0; l < mm; l++) 
					for (H[l] = 0., j = 0; j < mm; j++) H[l] += TX[l][j]*TH[k][j][num]; 
					for (           l = 0; l < mm; l++) TH[k][l][num] = H[l]; 
				}
			}
		}
	}
	delete_struct(TX);
	delete_struct(TD);
	delete_struct(TD_dop);
	delete_struct(H);
}


/////////////////////////////////////////////////////////////////////////////////
//...формирование внутренних матриц в аналитическом методе для граничных функций;
void CLame3D::eshelby_matrix_3body(double **& TT, double **& TX, int n, double **& TD, double **& TH, int num_bound)
{
	int nn = 6*n+3, mm = nn-12;

///////////////////////////////////////////////////////////
//...распределяем матрицы для решения граничного уравнения;
	if (! TT) set_matrix(TT, num_bound*nn, num_bound*nn);
	else		clean_matrix(TT, num_bound*nn, num_bound*nn);

	if (! TX) set_matrix(TX, num_bound*nn, nn);
	else		clean_matrix(TX, num_bound*nn, nn);

	if (n > 1) {
	  	if (! TD) set_matrix(TD, num_bound*mm, nn*(num_bound/2));	
		else		clean_matrix(TD, num_bound*mm, nn*(num_bound/2));

		if (! TH) set_matrix(TH, num_bound*mm, nn);	
		else		clean_matrix(TH, num_bound*mm, nn);
	}

/////////////////////////////////////////////////////////////
//...находим сшивочные коэффициенты для граничного уравнения;
	C2[n] =  1.;
	A2[n] =  1.;
	double R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1), RR1 = R1*R1, RR2 = R2*R2, Rn1 = RR1/(2.*n-1.), Rn2 = RR2/(2.*n-1.),
			 GM = get_param(NUM_SHEAR), GI = get_param(NUM_SHEAR+NUM_SHIFT), GL = get_param(NUM_SHEAR+NUM_SHIFT*2),
			 nM = get_param(NUM_SHEAR+1), nI = get_param(NUM_SHEAR+1+NUM_SHIFT), nL = get_param(NUM_SHEAR+1+NUM_SHIFT*2);
	double Y2 = C2[n]*.25/(1.-nM), d = (n+2.)*(n+3.)/((1.-nI)*(2.*n+3.)),
			 Z0 = (GL/GM-1.)*RR1*Y2/RR2*((1.-d)/((2.*(n+3.)*GL/GI-1.)-(GL/GI-1.)*d)-(1.-(n+2.)*(n+3.)/((1.-nL)*(2.*n+3.)))/(2.*n+5.)),
			 Z1 = (1.-(GL/GM-1.)/(2.*n+5.)*(1.-(n+2.)*(n+3.)/((1.-nL)*(2.*n+3.))))*Y2,
			 dd = powI(R2/R1, 2*n+1), regul = (n ? 0. : 1.);
	C0[n] = (GL/GM-1.)*Y2/(RR2*((GL/GI-.5/(n+3.))*(2.*n+3.)-(GL/GI-1.)*.5*(n+2.)/(1.-nI)));
	C1[n] = (GL/GM-1.)*2.*(n+3.)*Y2/(RR2*(2.*n+3.)*(2.*n+5.));
	A1[n] = (Z1*RR2-Z0/dd*RR1)/(RR2-RR1)*4.*(1.-nL);
	B1[n] = (Z0/dd-Z1)*RR1*RR2/((RR2-RR1)*(2.*n+3.)*(6.+n-4.*nL))*4.*(1.-nL);
	B2[n] = -RR2/((2.*n+3.)*(6.+n-4.*nM))*A2[n];
 
////////////////////////////////////////////////////////
//...заполняем матрицы для решения граничного уравнения;
 	double f = .25/(GI*(1.-nI)), d0, d1, d2, d3, d4; //...перемещения на границе включений - со стороны включения;
	unit_f_matrix(TT, n, 0, 0, (3.-4.*nI)*f);
	rgradf_matrix(TT, n, 0, 0, -f);
	grdivf_transf(TD, n, 0, 0, -Rn1*f);
	unit_f_matrix(TT, n, 0, 2*nn, d1 = -(d0 = C0[n]*RR1*f)*(2.*n+1.)*(6.+n-4.*nI));
	rgradf_matrix(TT, n, 0, 2*nn, d1);
	rdivrf_matrix(TT, n, 0, 2*nn, d2 = d1+d0*(n+2.)*(2.*n+5.));

///////////////////////////////////////////////////////////////////
//...перемещения на границе включений - со стороны межфазного слоя;
	f = -.25/(GL*(1.-nL));
	unit_f_matrix(TT, n, 0, nn, d3 = (3.-4.*nL)*f);
	rgradf_matrix(TT, n, 0, nn, -f);
	grdivf_transf(TD, n, 0, nn, -Rn1*f);
	unit_f_matrix(TT, n, 0, 2*nn, d3*A1[n]*dd+(d1 = -(B1[n]*dd*f/RR1+(d0 = C1[n]*RR1*f))*(2.*n+1.)*(6.+n-4.*nL)));
	rgradf_matrix(TT, n, 0, 2*nn, -f*A1[n]*dd+ d1);
	rdivrf_matrix(TT, n, 0, 2*nn, d2 = d1+d0*(n+2.)*(2.*n+5.));

/////////////////////////////////////////////////////////////////////////
//...перемещения на границе межфазного слоя - со стороны межфазного слоя;
	f = -f; d3 = -d3;
	unit_f_matrix(TT, n, nn, nn, d3);
	rgradf_matrix(TT, n, nn, nn, -f);
	grdivf_transf(TD, n, mm, nn, -Rn2*f);
	unit_f_matrix(TT, n, nn, 2*nn, d3*A1[n]+(d1 = -(B1[n]*f/RR2+(d0 = C1[n]*RR2*f))*(2.*n+1.)*(6.+n-4.*nL)));
	rgradf_matrix(TT, n, nn, 2*nn, -f*A1[n]+ d1);
	rdivrf_matrix(TT, n, nn, 2*nn, d2 = d1+d0*(n+2.)*(2.*n+5.));

/////////////////////////////////////////////////////////////////
//...перемещения на границе межфазного слоя - со стороны матрицы;
	f = -.25/(GM*(1.-nM));
	unit_f_matrix(TT, n, nn, 2*nn, (d3 = (3.-4.*nM)*f)*C2[n]);
	rgradf_matrix(TT, n, nn, 2*nn, -f*C2[n]);
	unit_f_matrix(TT, n, nn, 3*nn, d3*A2[n]);
	rgradf_matrix(TT, n, nn, 3*nn, -f*A2[n]);
	unit_f_matrix(TT, n, nn, 3*nn, d1 = -B2[n]*f/RR2*(6.+n-4.*nM)*(2.*n+1.));
	rgradf_matrix(TT, n, nn, 3*nn, d1);
	rdivrf_matrix(TT, n, nn, 3*nn, d1);

//////////////////
//...правая часть;
	f = -f; d3 = -d3;
	unit_f_matrix(TX, n, nn, 0, d3);
	rgradf_matrix(TX, n, nn, 0, -f);
	grdivf_transf(TH, n, mm, 0, -Rn2*f);

////////////////////////////////////////////////////////////
//...напряжения на границе включений - со стороны включения;
	f = .5/(1.-nI);
	unit_f_matrix(TT, n, 2*nn, 0, n*(1.-2.*nI)*f);
	rgradf_matrix(TT, n, 2*nn, 0, (2.-n-2.*nI)*f);
	rdivrf_matrix(TT, n, 2*nn, 0,  2.*nI*f);
	grdivf_transf(TD, n, 2*mm, 0, Rn1*(2.-n)*f);
	unit_f_matrix(TT, n, 2*nn, 2*nn, d1 = -(sqr(n+2.)*f-1.)*(d0 = C0[n]*RR1)*(2.*n+1.));
	rgradf_matrix(TT, n, 2*nn, 2*nn, d1);
	rdivrf_matrix(TT, n, 2*nn, 2*nn, d2 = d0*(2.*n+1.-2.*(n+2.)*f));

//////////////////////////////////////////////////////////////////
//...напряжения на границе включений - со стороны межфазного слоя;
	f = -.5/(1.-nL);
	unit_f_matrix(TT, n, 2*nn, nn, n*(1.-2.*nL)*f);
	rgradf_matrix(TT, n, 2*nn, nn, (2.-n-2.*nL)*f);
	rdivrf_matrix(TT, n, 2*nn, nn,  2.*nL*f);
	grdivf_transf(TD, n, 2*mm, nn, Rn1*(2.-n)*f);
	unit_f_matrix(TT, n, 2*nn, 2*nn, A1[n]*dd*(1.+n)*(2.*nL-1.)*f+(d1 = -((d4 = -B1[n]*dd*f/RR1*(n+3.)*(6.+n-4.*nL))+(d0 = C1[n]*RR1)*(sqr(n+2.)*f+1.))*(2.*n+1.))+regul);
	rgradf_matrix(TT, n, 2*nn, 2*nn, A1[n]*dd*(3.+n - 2.*nL)*f+d1);
	rdivrf_matrix(TT, n, 2*nn, 2*nn, A1[n]*dd*2.*nL*f+(d2 = -d4*(2.*n+1.)-d0*(2.*n+1.+2.*(n+2.)*f)));

////////////////////////////////////////////////////////////////////////
//...напряжения на границе межфазного слоя - со стороны межфазного слоя;
	f = -f;
	unit_f_matrix(TT, n, 3*nn, nn, n*(1.-2.*nL)*f);
	rgradf_matrix(TT, n, 3*nn, nn, (2.-n-2.*nL)*f);
	rdivrf_matrix(TT, n, 3*nn, nn,  2.*nL*f);
	grdivf_transf(TD, n, 3*mm, nn, Rn2*(2.-n)*f);
	unit_f_matrix(TT, n, 3*nn, 2*nn, A1[n]*(1.+n)*(2.*nL-1.)*f+(d1 = -((d4 = -B1[n]*f/RR2*(n+3.)*(6.+n-4.*nL))+(d0 = C1[n]*RR2)*(sqr(n+2.)*f-1.))*(2.*n+1.)));
	rgradf_matrix(TT, n, 3*nn, 2*nn, A1[n]*(3.+n - 2.*nL)*f+d1);
	rdivrf_matrix(TT, n, 3*nn, 2*nn, A1[n]*2.*nL*f+(d2 = -d4*(2.*n+1.)+d0*(2.*n+1.-2.*(n+2.)*f)));

////////////////////////////////////////////////////////////////
//...напряжения на границе межфазного слоя - со стороны матрицы;
	f = -.5/(1.-nM);
	unit_f_matrix(TT, n, 3*nn, 2*nn, C2[n]*(d3 = (1.+n)*(2.*nM-1.)*f));
	rgradf_matrix(TT, n, 3*nn, 2*nn, C2[n]*(3.+n-2.*nM)*f);
	rdivrf_matrix(TT, n, 3*nn, 2*nn, C2[n]* 2.*nM*f);
	unit_f_matrix(TT, n, 3*nn, 3*nn, A2[n]*d3);
	rgradf_matrix(TT, n, 3*nn, 3*nn, A2[n]*(3.+n-2.*nM)*f);
	rdivrf_matrix(TT, n, 3*nn, 3*nn, A2[n]* 2.*nM*f);
	unit_f_matrix(TT, n, 3*nn, 3*nn, d1 = B2[n]*f/RR2*(n+3.)*(6.+n-4.*nM)*(2.*n+1.)+regul);
	rgradf_matrix(TT, n, 3*nn, 3*nn, d1);
	rdivrf_matrix(TT, n, 3*nn, 3*nn, d1);

//////////////////
//...правая часть;
	f = -f;
	unit_f_matrix(TX, n, 3*nn, 0, n*(1.-2.*nM)*f);
	rgradf_matrix(TX, n, 3*nn, 0, (2.-n-2.*nM)*f);
	rdivrf_matrix(TX, n, 3*nn, 0,  2.*nM*f);
	grdivf_transf(TH, n, 3*mm, 0, Rn2*(2.-n)*f);
}

/////////////////////////////////////////////////////////////////////////////////
//...формирование внутренних матриц в аналитическом методе для граничных функций;
void CLame3D::eshelby_matrix_2body(double **& TT, double **& TX, int n, double **& TD, double **& TH, int num_bound)
{
	int nn = 6*n+3, mm = nn-12;

///////////////////////////////////////////////////////////
//...распределяем матрицы для решения граничного уравнения;
	if (! TT) set_matrix(TT, num_bound*nn, num_bound*nn);
	else		clean_matrix(TT, num_bound*nn, num_bound*nn);

	if (! TX) set_matrix(TX, num_bound*nn, nn);
	else		clean_matrix(TX, num_bound*nn, nn);

	if (n > 1) {
	  	if (! TD) set_matrix(TD, num_bound*mm, nn*(num_bound/2));	
		else		clean_matrix(TD, num_bound*mm, nn*(num_bound/2));

		if (! TH) set_matrix(TH, num_bound*mm, nn);	
		else		clean_matrix(TH, num_bound*mm, nn);
	}

/////////////////////////////////////////////////////////////
//...находим сшивочные коэффициенты для граничного уравнения;
	double R1 = get_param(NUM_GEOMT), R2 = get_param(NUM_GEOMT+1), RR1 = R1*R1, Rn1 = RR1/(2.*n-1.),
			 GI = get_param(NUM_SHEAR+NUM_SHIFT), GL = get_param(NUM_SHEAR+NUM_SHIFT*2),
			 nI = get_param(NUM_SHEAR+1+NUM_SHIFT), nL = get_param(NUM_SHEAR+1+NUM_SHIFT*2);
	double dd = powI(R2/R1, 2*n+1), regul = (n ? 0. : 1.);
	A1[n] = 1.;
	B1[n] = -RR1/((2.*n+3.)*(6.+n-4.*nL))*A1[n];

////////////////////////////////////////////////////////
//...заполняем матрицы для решения граничного уравнения;
 	double f = .25/(GI*(1.-nI)), d1, d2, d3, d4; //...перемещения на границе включений - со стороны включения;
	unit_f_matrix(TT, n, 0, 0, (3.-4.*nI)*f);
	rgradf_matrix(TT, n, 0, 0, -f);
	grdivf_transf(TD, n, 0, 0, -Rn1*f);

///////////////////////////////////////////////////////////////////
//...перемещения на границе включений - со стороны межфазного слоя;
	f = -.25/(GL*(1.-nL)); d3 = (3.-4.*nL)*f;
	unit_f_matrix(TT, n, 0, nn, d3*A1[n]*dd+(d1 = -(d4 = B1[n]*dd*f/RR1)*(2.*n+1.)*(6.+n-4.*nL)));
	rgradf_matrix(TT, n, 0, nn, -f*A1[n]*dd+ d1);
	rdivrf_matrix(TT, n, 0, nn, d2 = d1);

//////////////////
//...правая часть;
	f = -f; d3 = -d3;
	unit_f_matrix(TX, n, 0, 0, d3);
	rgradf_matrix(TX, n, 0, 0, -f);
	grdivf_transf(TH, n, 0, 0, -Rn1*f);

////////////////////////////////////////////////////////////
//...напряжения на границе включений - со стороны включения;
	f = .5/(1.-nI);
	unit_f_matrix(TT, n, nn, 0, n*(1.-2.*nI)*f);
	rgradf_matrix(TT, n, nn, 0, (2.-n-2.*nI)*f);
	rdivrf_matrix(TT, n, nn, 0,  2.*nI*f);
	grdivf_transf(TD, n, mm, 0, Rn1*(2.-n)*f);

//////////////////////////////////////////////////////////////////
//...напряжения на границе включений - со стороны межфазного слоя;
	f = -.5/(1.-nL);
	unit_f_matrix(TT, n, nn, nn, A1[n]*dd*(1.+n)*(2.*nL-1.)*f+(d1 = -(d4 = -B1[n]*dd*f/RR1*(n+3.)*(6.+n-4.*nL))*(2.*n+1.))+regul);
	rgradf_matrix(TT, n, nn, nn, A1[n]*dd*(3.+n - 2.*nL)*f+d1);
	rdivrf_matrix(TT, n, nn, nn, A1[n]*dd*2.*nL*f+(d2 = -d4*(2.*n+1.)));

//////////////////
//...правая часть;
	f = -f;
	unit_f_matrix(TX, n, nn, 0, n*(1.-2.*nL)*f);
	rgradf_matrix(TX, n, nn, 0, (2.-n-2.*nL)*f);
	rdivrf_matrix(TX, n, nn, 0,  2.*nL*f);
	grdivf_transf(TH, n, mm, 0, Rn1*(2.-n)*f);
}

//////////////////////////////////////////////////////////////
//...сложение коллокационных векторов с учетом матрицы Эшелби;
void CLame3D::add_collocation(int i, int m, int j, int n, double * ff, double * ff_reg)
{
	if (TT && TH) {
		double * px = B[i].shape->FULL(solver.hh[i][0][j], 0, 0), 
				 * py = B[i].shape->FULL(solver.hh[i][0][j], 0, 1), 
				 * pz = B[i].shape->FULL(solver.hh[i][0][j], 0, 2), * p, 
					fd = B[i].shape->get_R(0)/fabs(B[i].mp[8]), RR0_inv = sqr(B[i].shape->get_R_inv(0)), fd_inv = fd*RR0_inv, dd;
		int   N_mpl = B[i].shape->get_N(0), nn, mm, kk, k, l, ll, shift, shifn, num;
		for (k = N_mpl; k >= 0; k--) {
			nn = 2*k+1; kk = k*k; mm = nn; ll = kk; shift = n*nn*3;
			if (TT[k]) { //...вычисляем текущую степень;
				if (ff) dd = ff[k]*fd; else dd = 1.;
				if (ff_reg) dd *= ff_reg[k];
				for (shifn = num = 0; num < 3; num++, shifn += nn)
				for (p = B[i].shape->FULL(solver.hh[i][0][m], 0, num), l = 0; l < nn; l++)
				for (j = 0; j < nn; j++) p[kk+l] += dd*(px[kk+j]*TT[k][shift+j][shifn+l]+py[kk+j]*TT[k][shift+j+nn][shifn+l]+pz[kk+j]*TT[k][shift+j+nn*2][shifn+l]);
			}
			if (TH[k] && k+2 <= N_mpl) { //...правим старшую степень;
				mm += 4; ll += 4*(k+1);
				if (ff) dd = ff[k]*fd_inv; else dd = RR0_inv;
				if (ff_reg) dd *= ff_reg[k];
				for (shifn = num = 0; num < 3; num++, shifn += mm)
				for (p = B[i].shape->FULL(solver.hh[i][0][m], 0, num), l = 0; l < mm; l++)
				for (j = 0; j < nn; j++) p[ll+l] += dd*(px[kk+j]*TH[k][shift+j][shifn+l]+py[kk+j]*TH[k][shift+j+nn][shifn+l]+pz[kk+j]*TH[k][shift+j+nn*2][shifn+l]);
			}
		}
	}
}

//////////////////////////////////////
//...вычисление невязки между блоками;
void CLame3D::block_descrap(char * OUT_FILE)
{
   int N_elem = UnPackInts(get_param(NUM_QUAD));
   char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid_el * bnd = (CGrid_el *)CreateNodes(GRID_EL_NODES);

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(3);

	CGrid_el * gauss_bnd = (CGrid_el * )CreateNodes(GRID_EL_NODES);
				gauss_bnd->add_params(1);

///////////////////////////////////////////
//...вычисление среднеквадратичной невязки;
	FILE * OUT = OUT_FILE ? fopen(OUT_FILE, "w") : NULL;
	sprintf(msg, "Block descrapency...");
	Message(msg); fprintf(OUT, "%s\n", msg);

   double pp[6], P0[6], Po[12];
	int  k, l, i, j, j_surf, m;

	for (k = 0; k < N; k++) if (B[k].bar && B[k].link) {
		sprintf(msg, "block %4i: ", k);
		Message(msg); fprintf(OUT, "%s\n", msg);

		for (i = 0; i < B[k].link[0]; i++) if ((j = B[k].link[i+1]) >= 0) {
			bnd->zero_grid(); 
			m = block_comput(bnd, k, j, sqr(get_param(4)), j_surf, 1);
       
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
					gauss_bnd->facet_QG(Po, N_elem, NULL_STATE, NULL_STATE);
					if (NUM_PHASE >= i+1) { //...коррекция квадратур;
						if (m && bar && bar->graph && -j_surf+SRF_STATE < bar->graph[0]) gauss_bnd->QG_tria_surface(bar->ce[-j_surf+SRF_STATE]->mp, Po);
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

						gauss_bnd->facet_QG(Po, N_elem, OK_STATE, NULL_STATE);
						if (NUM_PHASE >= i+1) { //...коррекция квадратур;
							if (m && bar && bar->graph && -j_surf+SRF_STATE < bar->graph[0]) gauss_bnd->QG_quad_surface(bar->ce[-j_surf+SRF_STATE]->mp, NULL, Po);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								block_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
							}
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
			}

/////////////////////////////////////
//...выисление невязки между блоками;
			if (NUM_PHASE >= i) {
				double F[6], sum = 0., norm = 0.;
				for (int lp = 0; lp < block_bnd->N; lp++) if (block_bnd->hit[lp]) {

					memset(F, 0, 6*sizeof(double));
					GetFuncAllValues(block_bnd->X[lp], block_bnd->Y[lp], block_bnd->Z[lp], F,   k, DISPL_VALUE);
					GetFuncAllValues(block_bnd->X[lp], block_bnd->Y[lp], block_bnd->Z[lp], F+3, j, DISPL_VALUE);

					sum  += (sqr(F[0]-F[3])+sqr(F[1]-F[4])+sqr(F[2]-F[5]))*block_bnd->get_param(0, lp);
					norm += (sqr(F[0])+sqr(F[1])+sqr(F[2]))*block_bnd->get_param(0, lp);
				}
				sum  = sqrt(sum);
				norm = sqrt(norm);

				sprintf(msg, "            block_bnd->N = %i  sum = %g", block_bnd->N, sum/(1.+norm));
				Message(msg); fprintf(OUT, "%s\n", msg);
			}
			block_bnd->add_buffer(block_bnd->N);
		}
	}
	sprintf(msg, "");
	Message(msg); fprintf(OUT, "%s\n", msg);
	if (OUT) fclose(OUT);

   delete block_bnd;
	delete gauss_bnd;
   delete bnd;
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CLame3D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 0., 0., 1.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp) {
		int m = solver.id_norm;
//////////////////////////////////////////////////////////////////
//memset(solver.hh[i][0][0], 0, solver.dim[i]*sizeof(double));
//for (int kk = 13; kk < 14; kk++) {
//	B[i].shape->FULL(solver.hh[i][0][0], 0, 0)[kk] = 1.;	
//	B[i].shape->FULL(solver.hh[i][0][0], 0, 1)[kk] = 1.;	
//	B[i].shape->FULL(solver.hh[i][0][0], 0, 2)[kk] = 1.;	
//}
//B[i].shape->set_potential(solver.hh[i][0][0], 0);
//////////////////////////////////////////////////////////////////
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
				B[i].shape->parametrization_hess(P, 1);

				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   0);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 0);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 0);
				B[i].shape->norm_common(F);

				if (solv == SPECIAL_SOLVING && id_variant == 0) F[2] -= Z;
				if (solv == SPECIAL_SOLVING && id_variant == 1) F[0] -= Z;
		}		break;
		case DISPL_ESHE_VALUE: {
///////////////////////////////
//...calculation displacements;
				B[i].shape->parametrization_hess(P, 1);

				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   0);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], 0);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], 0);
				B[i].shape->norm_common(F);
				if (B[i].link[NUM_PHASE] == -1) {
					double nju = get_param(NUM_SHEAR+1), G0 = .5/get_param(NUM_SHEAR), BBB = G0/(1.+nju), AAA = -nju*BBB;
					F[0] += X*AAA;
					F[1] += Y*AAA;
					F[2] += Z*BBB;
				}
				if (/*solv == SPECIAL_SOLVING && */id_variant == 0) F[2] -= Z;
				if (solv == SPECIAL_SOLVING && id_variant == 1) F[0] -= Z;
		}		break;
		case STRESS_X_VALUE: { 
//////////////////////////////////////////////////
//...calculation stress tensor (txx, txy and txz);
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

				jump4_x(P, i, 0); 
				jump4_y(P, i, 1); 
				jump4_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				B[i].shape->norm_common(F);
		}    break;
		case STRESS_Y_VALUE: { 
//////////////////////////////////////////////////
//...calculation stress tensor (tyx, tyy and tyz);
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

				jump4_x(P, i, 0); 
				jump4_y(P, i, 1); 
				jump4_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				B[i].shape->norm_common(F);
		}    break;
		case STRESS_Z_VALUE: { 
//////////////////////////////////////////////////
//...calculation stress tensor (tzx, tzy and tzz);
				P[5] = 1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

				jump4_x(P, i, 0); 
				jump4_y(P, i, 1); 
				jump4_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				B[i].shape->norm_common(F);
		}    break;
		case STRESS_Z_ESHE_VALUE: { 
//////////////////////////////////////////////////
//...calculation stress tensor (tzx, tzy and tzz);
				P[5] = 1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);

				jump4_x(P, i, 0); 
				jump4_y(P, i, 1); 
				jump4_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				B[i].shape->norm_common(F);
				if (B[i].link[NUM_PHASE] == -1) {
					F[2] += P[5];
				}
		}    break;
		case STRESS_R_VALUE: { 
////////////////////////////////////////////////
//...calculation surface forces (px, py and pz);
				double RR = sqrt(sqr(X)+sqr(Y)+sqr(Z)); 
				if (RR < 1e-16) { P[3] = P[4] = 0.; P[5] = 1.; }
				else {				P[3] = X/RR; P[4] = Y/RR; P[5] = Z/RR;	}
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if ((B[i].type & ERR_CODE) == CLAYER_BLOCK) hessian_deriv_N(P, i);

				jump4_x(P, i, 0); 
				jump4_y(P, i, 1); 
				jump4_z(P, i, 2); 

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				F[2] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				B[i].shape->norm_common(F);
		}    break;
		case POTENTIAL_VALUE: {
//////////////////////////////////////////
//...calculation potentials fX, fY and fZ;
				B[i].shape->parametrization_hess(P, 1);

				F[0]  = B[i].shape->potential(0, B[i].shape->p_cpy, id_variant);
				F[1]  = B[i].shape->potential(1, B[i].shape->p_cpy-B[i].shape->get_NN(), id_variant);
				F[2]  = B[i].shape->potential(2, B[i].shape->p_cpy-B[i].shape->get_NN()*2, id_variant);
				B[i].shape->norm_common(F);
		}		break;
		default: F[0] = i; F[1] = F[2] = 0.;
		}
	}
}

////////////////////////////////////////////////
//...calculaion additional member of asymptotic;
void CLame3D::GetEnergy(double * energy, int id_variant, Num_Comput Num)
{
	int N_elem = UnPackInts(get_param(NUM_QUAD));

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

////////////////////////////////////////////////////////////
//...discrete norm on inclusion and jump boundary condition;
	for (int k = 0; k < N; k++) 
	if (B[k].bar && B[k].link && B[k].link[NUM_PHASE] == -1) {
		solver.struct_init(k, NULL_STATE);

		int nums_facet = B[k].bar->arcs_number();
      for (int i = 1; i < nums_facet; i++) if (B[k].bar->ce[i]->segms_id()) {

			B[k].bar->ce[i]->segms_QG(gauss_bnd, N_elem);

			EnergyAll(gauss_bnd, k, energy, id_variant, Num);
			gauss_bnd->add_buffer(gauss_bnd->N);
		}
	}
	delete gauss_bnd;
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CLame3D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k] && 
							solver.TL && solver.TL[k] && solver.hh[k][0].GetMatrix()) {
		double ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
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

/////////////////////////////////////////////////////////////////////////
//...вычисление эффективных характеристик в классической слоистой модели;
double CLame3D::TakeLayer_kk(int N, double * ff, double * kk)
{
	double ** matr = NULL; set_matrix(matr, 2*N, 2*N+1);
	int i1, i2, m1, m2, m;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	for ( m2 = 2, i2 = 1, m = m1 = i1 = 0; i1 < N; i1++, i2 = (i1+1)%N, m1 = 2*i1, m2 = 2*i2) {
		matr[m][m1] = ff[i1]; matr[m][m1+1] =  1.0; matr[m][m2+1] = -1.0;  matr[m][2*N] = i2 ? 0. : 1.; m++;
		if (i2) {
			matr[m][m1] = kk[i1]; matr[m][m2] = -kk[i2]; m++;
		}
	}
	for ( m1 = i1 = 0; i1 < N; i1++, m1 = 2*i1) { //...добавляем условие среднего по периоду;
		matr[m][m1] = 0.5*ff[i1]*ff[i1];	matr[m][m1+1] = ff[i1];
	}

///////////////////////////////////////
//...решаем систему линейных уравнений;
	double sum = 0.;
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
//...вычисляем эффективные модули;
	for (i = 0; i < N; i++) sum += ff[i]*kk[i]*matr[2*i][dim_N];
err:
	delete_struct(matr); delete_struct(ii);
	return(sum);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективных характеристик в одномерной слоистой модели в сферической симметрии;
double CLame3D::TakeLayer_kk(int N, double * ff, double * kv, double * mu)
{
	double ** matr = NULL, s0, s1, sum = 0.; set_matrix(matr, 2*N, 2*N+1);
	int i1, m1, m = 2*N-1;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] =  1.; matr[1][0] = kv[0];
		matr[m][m] = -1.; matr[m-1][m+1] = 1.; s0 = ff[0];
		if (N > 1) {
			matr[0][1]   = -1.;     matr[0][2]   = -1.;								  matr[1][1] = -kv[1]; matr[1][2] = 4./3.*mu[1]; s1 = s0+ff[1];
			matr[m][m-2] = kv[N-1]; matr[m][m-1] = -4./3.*mu[N-1]*(1.-ff[N-1]); matr[m-1][m-2] = 1.; matr[m-1][m-1] = (1.-ff[N-1]);
			if (N > 2) {
				for ( m = 2, m1 = 1, i1 = 1; i1 < N-1; i1++, m1 += 2) {
					matr[m][m1] = 1.;     matr[m][m1+1] =  s0/s1;				  matr[m][m1+2] = -1.0;       matr[m][m1+3] = -1.0;			  m++;
					matr[m][m1] = kv[i1]; matr[m][m1+1] = -s0/s1*4./3.*mu[i1]; matr[m][m1+2] = -kv[i1+1];  matr[m][m1+3] = 4./3.*mu[i1+1]; m++; s0 = s1; s1 = s0+ff[i1+1];
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
//
////////////////////////////////////////////////////////////////////////////////////////
////...вычисление эффективного КТР в одномерной слоистой модели в сферической симметрии;
//double CLame3D::TakeLayer_kk(int N, double * ff, double * kv, double * mu, double * alpha)
//{
//	double ** matr = NULL, s0 = ff[0], s1, f1, sum = 0., KH = TakeLayer_kk(N, ff, kv, mu); set_matrix(matr, 2*N, 2*N+1);
//	int i1, m1, m = 2*N-1, mm;
//
////////////////////////////////////////////
////...заполняем систему линейных уравнений;
//	if (N > 0) {
//		matr[0][0] = 1.; matr[0][m+1] = alpha[0]*(f1 = pow(s0, 1./3.)); 
//		matr[1][0] = kv[0]; matr[1][m+1] = kv[0]*alpha[0]*f1;
//		
//		matr[m-1][m] = 1.; 
//		matr[m][m]   = KH; 
//		if (N > 1) {
//			matr[0][1] = -1.; matr[0][2] = -1.; matr[0][m+1] = -alpha[1]*f1; 
//			matr[1][1] = -kv[1]; matr[1][2] = 4./3.*mu[1]; matr[1][m+1] = -kv[1]*alpha[1]*f1; s1 = s0+ff[1];
//			
//			matr[m-1][m-2] = 1.; matr[m-1][m-1] = (1.-ff[N-1]); matr[m-1][m+1] = alpha[N-1]; 
//			matr[m][m-2]   = kv[N-1]; matr[m][m-1] = -4./3.*mu[N-1]*(1.-ff[N-1]); matr[m][m+1] = kv[N-1]*alpha[N-1]; 
//			if (N > 2) {
//				for (mm = 2, m1 = 1, i1 = 1; i1 < N-1; i1++, m1 += 2) {
//					matr[mm][m1] = 1.; matr[mm][m1+1] = s0/s1; matr[mm][m1+2] = -1.0; matr[mm][m1+3] = -1.0; matr[mm][m+1] = (alpha[i1]-alpha[i1+1])*(f1 = pow(s1, 1./3.)); mm++;
//					matr[mm][m1] = kv[i1]; matr[mm][m1+1] = -s0/s1*4./3.*mu[i1]; matr[mm][m1+2] = -kv[i1+1]; matr[mm][m1+3] = 4./3.*mu[i1+1]; matr[mm][m+1] = (kv[i1]*alpha[i1]-kv[i1+1]*alpha[i1+1])*f1; mm++; s0 = s1; s1 = s0+ff[i1+1];
//				}
//			}
//		}
//	}
//
/////////////////////////////////////////
////...решаем систему линейных уравнений;
//	int dim_N = 2*N, * ii, i, k, l, k0, l0; ii = (int *)new_struct(2*N*sizeof(int));
//	for (i = 0; i < dim_N; i++) {
//		double f = 0.;
/////////////////////////////////////////
////...look for position maximal element;
//		for (k = 0; k < dim_N; k++)
//			if (ii[k] != 1) 
//				for (l = 0; l < dim_N; l++) 
//					if (! ii[l]) {
//						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]); 
//					}
//					else if (ii[l] > 1) goto err;
//		++(ii[l0]);
/////////////////////////////////////////////////////////////
////...swapping row for diagonal position of maximal element;
//		if (k0 != l0) 
//			for (l = 0; l <= dim_N; l++) {
//				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
//			}
//		if (matr[l0][l0] == 0.) goto err;
//////////////////////////////////
////...diagonal row normalization;
//		double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
//		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
///////////////////////////////////
////...elimination all outher rows;
//		for (k = 0; k < dim_N; k++)
//			if ( k != l0) {
//				finv = matr[k][l0]; matr[k][l0] = 0.;
//				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
//			}
//	}
//
////////////////////////////////////
////...вычисляем эффективный модуль;
//	sum = matr[dim_N-1][dim_N];
//err:
//	delete_struct(matr); delete_struct(ii);
//	return(sum);
//}

//////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для многослойной модели;
double take_system_sph(double ** matrix, int * ii, int dim_N, double mu_H, double nu_H)
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
}

///////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного модуля сдвига в слоистой модели сферической симметрии;
double CLame3D::TakeLayer_GG(int N, double * ff, double * kv, double * mu, double * nj, double eps, int max_iter)
{
	double ** matr = NULL, mu_M = mu[N-1], s0 = ff[0], s1, f0, f1; set_matrix(matr, 4*N, 4*N+1);
	int i1, m1, m = 4*N-1;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (N > 0) {
		matr[0][0] = 1.; matr[0][1] = -3.*nj[0]*(f0 = pow(s0, 2./3.)); matr[m-3][m-1] = -.75; matr[m-3][m] = -1.; matr[m-3][m+1] = 1.;
		matr[1][0] = 1.; matr[1][1] = -(3.5-2.*nj[0])*f0; matr[m-2][m-1] =  .5; matr[m-2][m] = -1.; matr[m-2][m+1] = 1.;
		matr[2][0] = mu[0]; matr[2][1] = 1.5*mu[0]*nj[0]*f0; matr[m-1][m-1] = 3.; matr[m-1][m] = .5; matr[m-1][m+1] = 1.;
		matr[3][0] = mu[0]; matr[3][1] = -mu[0]*(3.5+nj[0])*f0; matr[m][m-1] = -2.; matr[m][m] = -0.5; matr[m][m+1] = 1.;
		if (N > 1) {
			matr[0][2] = -1.; matr[0][3] = -(1.25-nj[1]); matr[0][4] = -.75; matr[0][5] = 3.*nj[1]*f0; s1 = s0+ff[1];
			matr[m-3][m-5] = 1.; matr[m-3][m-4] = (1.-ff[N-1])*(1.25-nj[N-1]); matr[m-3][m-3] = .75*(f1 = pow(1.-ff[N-1], 5./3.));  matr[m-3][m-2] = -3.*nj[N-1];

			matr[1][2] = -1.; matr[1][3] = -(0.5-nj[1]); matr[1][4] = .5; matr[1][5] = (3.5-2.*nj[1])*f0;            
			matr[m-2][m-5] = 1.; matr[m-2][m-4] = (1.-ff[N-1])*(0.5-nj[N-1]); matr[m-2][m-3] = -.5*f1;  matr[m-2][m-2] = -(3.5-2.*nj[N-1]);

			matr[2][2]  = -mu[1]; matr[2][3] = mu[1]*(5.-nj[1])*.5; matr[2][4] = 3.*mu[1]; matr[2][5] = -1.5*mu[1]*nj[1]*f0; 
			matr[m-1][m-5] = mu_M; matr[m-1][m-4] = -mu_M*(5.-nj[N-1])*.5*(1.-ff[N-1]); matr[m-1][m-3] = -3.*mu_M*f1; matr[m-1][m-2] = 1.5*mu_M*nj[N-1];

			matr[3][2]  = -mu[1]; matr[3][3] = -0.5*mu[1]*(1.+nj[1]); matr[3][4] = -2.*mu[1]; matr[3][5] = mu[1]*(3.5+nj[1])*f0; 
			matr[m][m-5] = mu_M; matr[m][m-4] = 0.5*mu_M*(1.+nj[N-1])*(1.-ff[N-1]); matr[m][m-3] = 2.*mu_M*f1; matr[m][m-2] = -mu_M*(3.5+nj[N-1]);
			if (N > 2) {
				for ( m = 4, m1 = 2, i1 = 1; i1 < N-1; i1++, m1 += 4) {
					matr[m][m1] = 1.; matr[m][m1+1] = s0/s1*(1.25-nj[i1]); matr[m][m1+2] = .75*(f0 = pow(s0/s1, 5./3.)); matr[m][m1+3] = -3.*nj[i1]*(f1 = pow(s1, 2./3.));      
					matr[m][m1+4] = -1.; matr[m][m1+5] = -(1.25-nj[i1+1]); matr[m][m1+6] = -0.75; matr[m][m1+7] = 3.*nj[i1+1]*f1; m++;

					matr[m][m1] = 1.; matr[m][m1+1] = s0/s1*(0.5-nj[i1]); matr[m][m1+2] = -.5*f0; matr[m][m1+3] = -(3.5-2.*nj[i1])*f1;      
					matr[m][m1+4] = -1.; matr[m][m1+5] = -(0.5-nj[i1+1]); matr[m][m1+6] = .5; matr[m][m1+7] = (3.5-2.*nj[i1+1])*f1; m++;

					matr[m][m1] = mu[i1]; matr[m][m1+1] = -mu[i1]*(5.-nj[i1])*.5*s0/s1; matr[m][m1+2] = -3.*mu[i1]*f0; matr[m][m1+3] = 1.5*mu[i1]*nj[i1]*f1;
					matr[m][m1+4] = -mu[i1+1]; matr[m][m1+5] = mu[i1+1]*(5.-nj[i1+1])*.5; matr[m][m1+6] = 3.*mu[i1+1]; matr[m][m1+7] = -1.5*mu[i1+1]*nj[i1+1]*f1; m++;

					matr[m][m1] = mu[i1]; matr[m][m1+1] = 0.5*mu[i1]*(1.+nj[i1])*s0/s1; matr[m][m1+2] = 2.*mu[i1]*f0; matr[m][m1+3] = -mu[i1]*(3.5+nj[i1])*s1;
					matr[m][m1+4] = -mu[i1+1]; matr[m][m1+5] = -0.5*mu[i1+1]*(1.+nj[i1+1]); matr[m][m1+6] = -2.*mu[i1+1]; matr[m][m1+7] = mu[i1+1]*(3.5+nj[i1+1])*s1; m++; s0 = s1; s1 = s0+ff[i1+1];
				}
			}
		}
	}
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1, KH = TakeLayer_kk(N, ff, kv, mu), 
		** matrix = NULL; set_matrix(matrix, m = 4*N, 4*N+1); int * ii = (int *)new_struct(m*sizeof(int)), k_iter = 0, k, l; 

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	for (mu_H1 = 0., l = 0; l < N; l++) mu_H1 += mu[l]*ff[l]; mu_H1 *= 1000.;
	nu_H = 0.5*(1.-mu_H1/KH);
	for (k = 0; k < m; k++)
	for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_sph(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do { 
		mu_H0 /= 2.; k_iter++;
		nu_H = 0.5*(1.-mu_H0/KH);
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_sph(matrix, ii, m, mu_H0, nu_H);
	}
	while(optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = 0.5*(1.-mu_H/KH);
		for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_sph(matrix, ii, m, mu_H, nu_H);
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

/////////////////////////////////////////////////
//...трехфазная модель для сферических включений;
double CLame3D::TakeEshelby_volm_two(double ff)
{
	double K1 = 2.*get_param(NUM_SHEAR+NUM_SHIFT)*(get_param(NUM_SHEAR+1+NUM_SHIFT)+1.)/(3.-6.*get_param(NUM_SHEAR+1+NUM_SHIFT)), 
			 K3 = 2.*get_param(NUM_SHEAR)*(get_param(NUM_SHEAR+1)+1.)/(3.-6.*get_param(NUM_SHEAR+1)),
 			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*get_param(NUM_SHEAR)));
	return(KH);
}

////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для трехфазной модели;
double take_system(double matrix[8][9], double mu_MH, double nu_H)
{
	int dim_N = 8, ii[8] = {0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[4][6] *= mu_MH;
	matrix[4][7] *= mu_MH;
	matrix[4][8] *= mu_MH;
	matrix[5][6] *= mu_MH;
	matrix[5][8] *= mu_MH;
	matrix[5][7] *= mu_MH*2.*(1.-2.*nu_H)/(5.-4.*nu_H);
	matrix[6][7] *= 2.*(5.-nu_H)/(5.-4.*nu_H);
	matrix[7][7] *= 2.*(1.+nu_H)/(5.-4.*nu_H);
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
	return(matrix[7][8]);
}

//////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (итерационный алгоритм);
double CLame3D::TakeEshelby_shear_two(double ff, double eps, int max_iter)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR1 = pow(ff, 1./3.), fR2 = RR1*RR1, fR3 = 1./(RR1*fR2), fR5 = fR3/fR2, 
			 mu_MI = mu_M/mu_I, 
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, nu_H, mu_MH, mu_MH0, mu_MH1, matrix[8][9],
			 matr[8][9] = {
					{ mu_MI, -2.*mu_MI*QI1*fR2, -1., -fR3, -3.*fR5, 2.*QM1*fR2, 0., 0., 0.},
					{ mu_MI, -mu_MI*fR2, -1.,  -QM3*fR3, 2.*fR5, fR2, 0., 0., 0.},
					{ 1.,  QI1*fR2, -1.,  QM4*fR3, 12.*fR5, -QM1*fR2, 0., 0., 0.},
					{ 1., -QI2*fR2, -1., -QM5*fR3, -8.*fR5,  QM2*fR2, 0., 0., 0.},
					{ 0., 0., 1., 1., 3., -2.*QM1, -3., -1., 1.},
					{ 0., 0., 1.,  QM3, -2., -1., 2., -1., 1.},
					{ 0., 0., 1., -QM4, -12., QM1, 12., 1., 1.},
					{ 0., 0., 1.,  QM5, 8.,  -QM2, -8., -1., 1.},
	};
	double alpha = -1.;
	if (alpha >= 0.) {
		double alpha_I = alpha/mu_I*.5, alpha_M = alpha/mu_M*.5;
		matr[1][0] = 1.-alpha_I; matr[1][1] = -(QI2-alpha_I)*fR2; matr[1][2] = alpha_M;
		matr[1][3] = alpha_M*QM3*fR3;	matr[1][4] = -2.*alpha_M*fR5;	matr[1][5] = -alpha_M*fR2;
		matr[3][0] = alpha_I; matr[3][1] = -alpha_I*fR2; matr[3][2] = -(1.+alpha_M);
		matr[3][3] = -(QM5+QM3*alpha_M)*fR3; matr[3][4] = -8.*(1.-alpha_M*.25)*fR5; matr[3][5] = (QM2+alpha_M)*fR2;
	}
	int k_iter = 0, k, l;
/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.;
	mu_MH1 = 2.*mu_M/(3.*KH)/10.;
	do { 
		mu_MH1 *= 10.; k_iter++;
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 8; k++)
		for (l = 0; l < 9; l++)	matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH1, nu_H);
	}
	while(optim > 0. && k_iter < max_iter); k_iter = 0;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H  = (1.5*KH*mu_MH/get_param(NUM_SHEAR)-1.)/(3.*KH*mu_MH/get_param(NUM_SHEAR)+1.);
		for (k = 0; k < 8; k++)
		for (l = 0; l < 9; l++)	matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH, nu_H);
		if (optim > 0.) mu_MH0 = mu_MH; else
		if (optim < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);

	return(2.*mu_M/(mu_MH0+mu_MH1));
}

/////////////////////////////////////////////////////////////
//...четырехфазная модель для сферических включений со слоем;
double CLame3D::TakeEshelby_volm(double ff, double ff_l)
{
	double K1 = 2.*get_param(NUM_SHEAR+NUM_SHIFT)*(get_param(NUM_SHEAR+1+NUM_SHIFT)+1.)/(3.-6.*get_param(NUM_SHEAR+1+NUM_SHIFT)), 
			 K2 = 2.*get_param(NUM_SHEAR+NUM_SHIFT*2)*(get_param(NUM_SHEAR+1+NUM_SHIFT*2)+1.)/(3.-6.*get_param(NUM_SHEAR+1+NUM_SHIFT*2)),
			 K3 = 2.*get_param(NUM_SHEAR)*(get_param(NUM_SHEAR+1)+1.)/(3.-6.*get_param(NUM_SHEAR+1)),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT*2), c0 = ff+ff_l, c1 = ff/c0, DD, KH;
	DD = ((K1-K3)-(1.-c1)*(K1-K2)*(1.+(K3-K2)/(K2+4./3.*G2)))/(1.+(1.-c1)*(K1-K2)/(K2+4./3.*G2));
	KH = K3+c0*DD/(1.+(1.-c0)*DD/(K3+4./3.*get_param(NUM_SHEAR)));
	return(KH);
}

///////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для четырехфазной модели;
double take_system(double matrix[12][13], double mu_MH, double nu_H)
{
	int dim_N = 12, ii[12] = {0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[8][10] *= mu_MH;
	matrix[8][11] *= mu_MH;
	matrix[8][12] *= mu_MH;
	matrix[9][10] *= mu_MH;
	matrix[9][12] *= mu_MH;
	matrix[9][11] *= mu_MH*2.*(1.-2.*nu_H)/(5.-4.*nu_H);
	matrix[10][11] *= 2.*(5.-nu_H)/(5.-4.*nu_H);
	matrix[11][11] *= 2.*(1.+nu_H)/(5.-4.*nu_H);
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
	return(matrix[11][12]);
}

/////////////////////////////////////////////////////////////////////////////
//...четырехфазная модель для сферического включения (итерационный алгоритм);
double CLame3D::TakeEshelby_shear(double ff, double ff_l, double eps, int max_iter)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K2 = 2.*mu_L*(1.+nu_L)/(3.-6.*nu_L),	
			 K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M), c0 = ff+ff_l, c1 = ff/c0,
			 DD = ((K1-K3)-(1.-c1)*(K1-K2)*(1.+(K3-K2)/(K2+4./3.*mu_L)))/(1.+(1.-c1)*(K1-K2)/(K2+4./3.*mu_L)),
			 KH = K3+c0*DD/(1.+(1.-c0)*DD/(K3+4./3.*mu_M)),
			 RR1 = pow(ff, 1./3.), f1R2 = RR1*RR1, f1R3 = 1./(RR1*f1R2), f1R5 = f1R3/f1R2, 
			 RR2 = pow(ff+ff_l, 1./3.), f2R2 = RR2*RR2, f2R3 = 1./(RR2*f2R2), f2R5 = f2R3/f2R2,
			 mu_LI = mu_L/mu_I, mu_LM = mu_L/mu_M, 
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QL1 = 3.*nu_L/(7.-4.*nu_L), 
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QL2 = (7.+2.*nu_L)/(7.-4.*nu_L), 
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QL3 = 2.*(1.-2.*nu_L)/(5.-4.*nu_L), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QL4 = 2.*(5.-nu_L)/(5.-4.*nu_L),
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QL5 = 2.*(1.+nu_L)/(5.-4.*nu_L),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, nu_H, mu_MH, mu_MH0, mu_MH1, matrix[12][13],
			 matr[12][13] = {
					{ mu_LI, -2.*mu_LI*QI1*f1R2, -1., -f1R3, -3.*f1R5, 2.*QL1*f1R2, 0., 0., 0., 0., 0., 0., 0.},
					{ mu_LI, -mu_LI*f1R2, -1.,  -QL3*f1R3, 2.*f1R5, f1R2, 0., 0., 0., 0., 0., 0., 0.},
					{ 1.,  QI1*f1R2, -1.,  QL4*f1R3, 12.*f1R5, -QL1*f1R2, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -QI2*f1R2, -1., -QL5*f1R3, -8.*f1R5,  QL2*f1R2, 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 1., f2R3, 3.*f2R5, -2.*QL1*f2R2, -mu_LM, -mu_LM*f2R3, -3.*mu_LM*f2R5, 2.*mu_LM*QM1*f2R2, 0., 0., 0.},
					{ 0., 0., 1., QL3*f2R3, -2.*f2R5, -f2R2, -mu_LM, -mu_LM*QM3*f2R3, 2.*mu_LM*f2R5, mu_LM*f2R2, 0., 0., 0.},
					{ 0., 0., 1., -QL4*f2R3, -12.*f2R5, QL1*f2R2, -1.,  QM4*f2R3, 12.*f2R5, -QM1*f2R2, 0., 0., 0.},
					{ 0., 0., 1., QL5*f2R3, 8.*f2R5, -QL2*f2R2,   -1., -QM5*f2R3, -8.*f2R5,  QM2*f2R2, 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 1., 1., 3., -2.*QM1, -3., -1., 1.},
					{ 0., 0., 0., 0., 0., 0., 1.,  QM3,  -2., -1.,  2., -1., 1.},
					{ 0., 0., 0., 0., 0., 0., 1., -QM4, -12., QM1, 12.,  1., 1.},
					{ 0., 0., 0., 0., 0., 0., 1.,  QM5,  8., -QM2, -8., -1., 1.},
	};
	int k_iter = 0, k, l;
/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.;
	mu_MH1 = 2.*mu_M/(3.*KH);
	do { 
		mu_MH1 *= 10.; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 12; k++)
		for (l = 0; l < 13; l++) matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH1, nu_H);
	}
	while(optim > 0.);
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H  = (1.5*KH*mu_MH/mu_M-1.)/(3.*KH*mu_MH/mu_M+1.);
		for (k = 0; k < 12; k++)
		for (l = 0; l < 13; l++) matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH, nu_H);
		if (optim > 0.) mu_MH0 = mu_MH; else
		if (optim < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);



///////////////////////////////////////////////////////////////////////////////////////////////////
//...распределение инвариантов тензора напряжений поля в среднем сечении (и на межфазных границах);
	//int id_visual = 0;
	//if (id_visual) {
	//	FILE * TST = fopen("CCohes3D_sphere_eshelby.dat", "w");
	//	double A = .7, X, Y = 0., Z, RR, func; int M = 400;
	//	for (i = 0; i <= M; i++)
	//	for (l = 0; l <= M; l++) {
	//		X = (-1.+i*2./M)*A; Z = (-1.+l*2./M)*A; RR = sqrt(X*X+Y*Y+Z*Z); func = Z;
	//		if (RR < RR1) func *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR)); else
	//		if (RR < 1.0) func *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
	//		fprintf(TST, "%g   %g   %g\n", Z, X, func-Z);
	//	}
	//	fclose(TST);
	//	if (NULL_STATE) { //...проверка сшивки в одной точке;
	//		X = 0.; Z = RR1; RR = sqrt(X*X+Y*Y+Z*Z); double func1 = Z, func2 = Z;
	//		func1 *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR));
	//		func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 

	//		X = 0.; Z = 1.; RR = sqrt(X*X+Y*Y+Z*Z); func1 = Z; func2 = Z;
	//		func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
	//	}



	return(2.*mu_M/(mu_MH0+mu_MH1));
}

//////////////////////////////////////////////////////////////////
//...явное решение системы уравнений Эшелби для трехфазной модели;
inline double take_system(double RR2, double mu_MI, double mu_MH, double nu_I, double nu_M, double nu_H, 
								  double QQM, double QRM, double QRI, double QM, double QI)
{
	double fR = 1.-RR2*sqr(sqr(RR2)), ffR = 1.-RR2*sqr(RR2),
	C1 = (mu_MH-1.)/(1.+QM),
	C0 = (mu_MH-1.)/(mu_MI+QI),
	B1 = (mu_MH-C1-1./(1.-fR)*(mu_MI*C0-C1))/(1.-sqr(RR2)),
	A1 =  mu_MH-C1-B1,
	UI = (mu_MI*(1.+nu_M)-(1.-2.*nu_M))*A1*ffR+(4.*mu_MI+1.)*.6*B1*fR-QRI*C0-QRM*C1-mu_MI*(1.+nu_H)+mu_MH*(1.-2.*nu_H),
	FI =  mu_MI*(1.4-nu_H)+mu_MH*(1.6-2.*nu_H),
	UM =  3.*(nu_M*A1+B1)-QQM*C1-((1.+nu_H)-mu_MH*(1.-2.*nu_H)),
	FM = (1.4-nu_H)+mu_MH*(1.6-2.*nu_H),
	UH = mu_MI-mu_MH,
	FH = 1.-mu_MH,
	QQ = UM*FI-UI*FM,
	hH = (UH*UM-FH*UI)/QQ,
	hM = (FH*FI-UH*FM)/QQ,
	optim = hM+hH;
//////////////
//...проверки:
	//double f1, f2, f3, f4, f5, f6, fM, fI;
	//f1 = A1+B1+C1-mu_MH;
	//f2 = -4.*(A1+B1)+C1*(7.+5.*nu_M)/(7.-10.*nu_M)+4.;
	//f3 = (1.-fR)*(A1+B1*sqr(RR2))+C1-mu_MI*C0;
	//f4 = 4.*(1.-fR)*(A1+B1*sqr(RR2))+C0*(7.+5.*nu_I)/(7.-10.*nu_I)-C1*(7.+5.*nu_M)/(7.-10.*nu_M);
	//fM = 1.-(1.4-nu_H)*hH-((1.+nu_M)*A1+2.4*B1-1.5*(7.+2.*nu_M)/(7.-10.*nu_M)*C1-(1.+nu_H))*hM;
	//f5 = fM-mu_MH-mu_MH*(1.6-2.*nu_H)*hH+((1.-2.*nu_M)*A1-0.6*B1-1.5*(7.-4.*nu_M)/(7.-10.*nu_M)*C1-mu_MH*(1.-2.*nu_H))*hM;
	//fI = fM+((1.+nu_M)*A1*(1.-ffR)+2.4*B1*(1.-fR)+1.5/sqr(RR2)*((7.+2.*nu_I)/(7.-10.*nu_I)*C0-(7.+2.*nu_M)/(7.-10.*nu_M)*C1))*hM;
	//f6 = fM-mu_MI*fI+((1.-2.*nu_M)*A1*(1.-ffR)-0.6*B1*(1.-fR)+1.5/sqr(RR2)*(mu_MI*(7.-4.*nu_I)/(7.-10.*nu_I)*C0-(7.-4.*nu_M)/(7.-10.*nu_M)*C1))*hM;
//////////////
	return(optim);
}

///////////////////////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (явное решение системы уравнений Эшелби);
double CLame3D::TakeEshelby_shear_sys(double ff, double eps, int max_iter)
{
	double K1 = 2.*get_param(NUM_SHEAR+NUM_SHIFT)*(get_param(NUM_SHEAR+1+NUM_SHIFT)+1.)/(3.-6.*get_param(NUM_SHEAR+1+NUM_SHIFT)), 
			 K3 = 2.*get_param(NUM_SHEAR)*(get_param(NUM_SHEAR+1)+1.)/(3.-6.*get_param(NUM_SHEAR+1)),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*get_param(NUM_SHEAR))),
			 QM = (1.75+1.25*get_param(NUM_SHEAR+1))/(7.-10.*get_param(NUM_SHEAR+1)),
			 QI = (1.75+1.25*get_param(NUM_SHEAR+NUM_SHIFT+1))/(7.-10.*get_param(NUM_SHEAR+NUM_SHIFT+1)),
			 nu_M = get_param(NUM_SHEAR+1), nu_I = get_param(NUM_SHEAR+NUM_SHIFT+1), 
			 mu_MI = get_param(NUM_SHEAR)/get_param(NUM_SHEAR+NUM_SHIFT), 
			 RR1 = pow(ff, 1./3.),
			 QQM = 9.*nu_M/(7.-10.*nu_M), 
			 QRI = 9.*nu_I/(7.-10.*nu_I)*sqr(RR1)*mu_MI,
			 QRM = 1.5*(mu_MI*(7.+2.*nu_M)-(7.-4.*nu_M))/(7.-10.*nu_M)*(1.-sqr(RR1)), 
			 RR2 = 1./RR1, nu_H, optim, mu_MH, mu_MH0, mu_MH1;
	int k_iter = 0;
/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.;
	mu_MH1 = 2.*get_param(NUM_SHEAR)/(3.*KH);
	do { 
		mu_MH1 *= 10.; 
		nu_H = (1.5*KH*mu_MH1/get_param(NUM_SHEAR)-1.)/(3.*KH*mu_MH1/get_param(NUM_SHEAR)+1.);
		optim = take_system (RR2, mu_MI, mu_MH1, nu_I, nu_M, nu_H, QQM, QRM, QRI, QM, QI);
	}
	while(optim > 0.);
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H  = (1.5*KH*mu_MH/get_param(NUM_SHEAR)-1.)/(3.*KH*mu_MH/get_param(NUM_SHEAR)+1.);
		optim = take_system (RR2, mu_MI, mu_MH, nu_I, nu_M, nu_H, QQM, QRM, QRI, QM, QI);
		if (optim > 0.) mu_MH0 = mu_MH; else
		if (optim < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);
//////////////////////////////////////////////////////////
//...формула из Кристенсена;
	//double mu_HM = 1.-15.*ff*(1.-nu_M)*(1.-1./mu_MI)/(7.-5.*nu_M+2.*(4.-5.*nu_M)/mu_MI);
	//mu_MH0 = mu_MH1 = 1./mu_HM;
//////////////////////////////////////////////////////////
	return(2.*get_param(NUM_SHEAR)/(mu_MH0+mu_MH1));
}

/////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного модуля сдвига через определитель (алгоритм Кристенсена);
inline double DETER3(double matr[8][8], int i1, int i2, int i3, int j1, int j2, int j3)
{	return matr[i1-1][j1-1]*matr[i2-1][j2-1]*matr[i3-1][j3-1]+
			 matr[i1-1][j2-1]*matr[i2-1][j3-1]*matr[i3-1][j1-1]+ 
			 matr[i2-1][j1-1]*matr[i3-1][j2-1]*matr[i1-1][j3-1]- 
			 matr[i3-1][j1-1]*matr[i2-1][j2-1]*matr[i1-1][j3-1]- 
			 matr[i2-1][j1-1]*matr[i1-1][j2-1]*matr[i3-1][j3-1]- 
			 matr[i3-1][j2-1]*matr[i2-1][j3-1]*matr[i1-1][j1-1];
}
double CLame3D::TakeEshelby_shear_det(double ff)
{
	double RR1 = pow(ff, 1./3.), fR2 = RR1*RR1, fR3 = 1./(RR1*fR2), fR5 = fR3/fR2, alpha = -1., 
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), 
			 mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+NUM_SHIFT+1), 
			 mu_MI = mu_M/mu_I, 
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), mu_MH, A, B, C, D,
			 matr[8][8] = {
					{ mu_MI, -2.*mu_MI*QI1*fR2, -1., -fR3, -3.*fR5, 2.*QM1*fR2, 0., 0.},
					{ mu_MI, -mu_MI*fR2, -1.,  -QM3*fR3, 2.*fR5, fR2, 0., 0.},
					{ 1.,  QI1*fR2, -1.,  QM4*fR3, 12.*fR5, -QM1*fR2, 0., 0.},
					{ 1., -QI2*fR2, -1., -QM5*fR3, -8.*fR5,  QM2*fR2, 0., 0.},
					{ 0., 0., 1., 1., 3., -2.*QM1, -3., -1.},
					{ 0., 0., 1.,  QM3, -2.,   -1., 2., -1.},
					{ 0., 0., 1., -QM4, -12., QM1, 12., -1.},
					{ 0., 0., 1.,  QM5, 8.,  -QM2, -8., -1.},
	};
	if (alpha >= 0.) {
		double alpha_I = alpha/mu_I*.5, alpha_M = alpha/mu_M*.5;
		matr[1][0] = 1.-alpha_I; matr[1][1] = -(QI2-alpha_I)*fR2; matr[1][2] = alpha_M;
		matr[1][3] = alpha_M*QM3*fR3;	matr[1][4] = -2.*alpha_M*fR5;	matr[1][5] = -alpha_M*fR2;
		matr[3][0] = alpha_I; matr[3][1] = -alpha_I*fR2; matr[3][2] = -(1.+alpha_M);
		matr[3][3] = -(QM5+QM3*alpha_M)*fR3; matr[3][4] = -8.*(1.-alpha_M*.25)*fR5; matr[3][5] = (QM2+alpha_M)*fR2;
	}

/////////////////////////////
//...вычисляем олпределители;
	A = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 7, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 7, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 7, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 7, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6))*5.;
	B = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 6, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 6, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 6, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 6, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6))*(-15.)+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 6, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 6, 7, 4, 5, 6)+
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 6, 7, 4, 5, 6)-
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 6, 7, 4, 5, 6)+
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6)-
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6))*(-5.)+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*10.+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 7, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 7, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 7, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6)-
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*10.;
	C = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 6, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 6, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)-
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 6, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 6, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6)-
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*(-20.);

/////////////////////////////////
//...решаем квадратное уравнение;
	D = B*B-4.*A*C;
	if (A < 0.) mu_MH = (-B-sqrt(D))/(2.*A); else
	if (A > 0.) mu_MH = (-B+sqrt(D))/(2.*A); else mu_MH = 1.;

	return(mu_M/mu_MH);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для трехфазной модели (модуль сдвига);
double take_system_shear(double matrix[8][9], double mu_MH, double nu_H)
{
	int dim_N = 8, ii[8] = {0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[4][6] *= mu_MH;
	matrix[4][7] *= mu_MH*(1.25-nu_H);
	matrix[4][8] *= mu_MH;
	matrix[5][6] *= mu_MH;
	matrix[5][7] *= mu_MH*(0.5-nu_H);
	matrix[5][8] *= mu_MH;
	matrix[6][7] *= (5.-nu_H)*.5;
	matrix[7][7] *= (1.+nu_H)*.5;
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
	return(matrix[7][8]);
}

//////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (итерационный алгоритм);
double CLame3D::TakeEshelby_shear(double ff, double nju1, double nju2, double E1, double E2)
{
	double mu_I = E1/(2.*(1.+nju1)), nu_I = nju1,
			 mu_M = E2/(2.*(1.+nju2)), nu_M = nju2,
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR1 = pow(ff, 1./3.), fR2 = RR1*RR1, c0 = ff, s0 = ff*fR2, 
			 QI1 = 1.5*(7.-4.*nu_I),
			 QM1 = 1.5*(7.-4.*nu_M), 
			 QI2 = 1.5*(7.+2.*nu_I),
			 QM2 = 1.5*(7.+2.*nu_M),
			 QI3 = 9.*nu_I, 
			 QM3 = 9.*nu_M, 
			 QM4 = (5.-nu_M)*.5,
			 QM5 = (1.+nu_M)*.5,
			 QM6 = 2.5-2.*nu_M,
			 QM7 = 0.5-nu_M,
			 QM8 = 1.25-nu_M,
			 QM9 = 1.-2.*nu_M, optim, sgn0, sgn1, nu_H, mu_H, mu_MH, mu_MH0, mu_MH1, matrix[8][9], eps = EE, max_iter = 100, eps0 = 1e-6,
			 matr[8][9] = {
					{ mu_M, -mu_M*QI3*fR2,    -mu_I, -mu_I*QM8, -mu_I*2.25, mu_I*QM3*fR2,       0., 0., 0.},
					{ mu_M, -mu_M*QI1*fR2,    -mu_I, -mu_I*QM7,  mu_I*1.5,  mu_I*QM1*fR2,       0., 0., 0.},
					{   1.,       QI3*fR2*.5,   -1.,       QM4,       9.,       -QM3*fR2*.5,    0., 0., 0.},
					{   1.,      -QI2*fR2,      -1.,      -QM5,      -6.,        QM2*fR2,       0., 0., 0.},
					{   0.,            0.,       1.,    QM8*c0,  2.25*s0,           -QM3,    -2.25, -1., 1.},
					{   0.,            0.,       1.,    QM7*c0,  -1.5*s0,           -QM1,      1.5, -1., 1.},
					{   0.,            0.,       1.,   -QM4*c0,   -9.*s0,            QM3*.5,    9.,  1., 1.},
					{   0.,            0.,       1.,    QM5*c0,    6.*s0,           -QM2,      -6., -1., 1.},
	};
	int k_iter = 0, k, l;
/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.; nu_H = -1.;
	for (k = 0; k < 8; k++)
	for (l = 0; l < 9; l++) matrix[k][l] = matr[k][l];
	optim = take_system_shear(matrix, mu_MH0, nu_H);
	if (optim < 0.) sgn0 = -1.; else sgn0 = 1.;

	mu_MH1 = 2.*mu_M/(3.*KH);
	do { 
		mu_MH1 *= 10.; k_iter++; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 8; k++)
		for (l = 0; l < 9; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_shear(matrix, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H = (1.5*KH*mu_MH/mu_M-1.)/(3.*KH*mu_MH/mu_M+1.);
		for (k = 0; k < 8; k++)
		for (l = 0; l < 9; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_shear(matrix, mu_MH, nu_H);
		if (optim*sgn0 > 0.) mu_MH0 = mu_MH; else
		if (optim*sgn0 < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);

	mu_H = 2.*mu_M/(mu_MH0+mu_MH1);
	if (sgn0*sgn1 > 0. || fabs(optim) > eps0 ) mu_H = -mu_H;
	return(mu_H);
}
#undef  Message
