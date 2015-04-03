#include "stdafx.h"

#include "cgrid_el.h"
#include "chydro3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CHydro3D::NUM_HESS = 14;
int CHydro3D::NUM_SHEAR = 6;
int CHydro3D::NUM_GEOMT = 4;
int CHydro3D::MAX_PHASE = 2;
int regul = 1;

//////////////////////////////////
//...destructor of CHydro3D class;
CHydro3D::~CHydro3D (void)
{
	for (int k = 0; TT && TT[k]; k++) {
		delete_struct(TT[k]);
		delete_struct(TH[k]);
	}
   delete_struct(TT);
   delete_struct(TH);
   delete_struct(TN);
}

//////////////////////////////////
//...initialization of the blocks;
int  CHydro3D::block_shape_init(Block<real_T> & B, int id_free)
{
	int k;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<real_T>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, sqr(B.mp[8])/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));

			if (get_param(NUM_SHEAR+1)) { //...подключаем две дополнительные экспоненциальные системы функций;
				B.shape->add_shape(CreateShapeD(SK3D_EXPP_SHAPE),	 NULL_STATE);
				B.shape->add_shape(CreateShapeD(SK3D_EXPP_SHAPE, 1), NULL_STATE);

				B.shape->init1(2, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 0);
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)));

				B.shape->init1(3, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, 0);
				B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)));
			}
		}
		else
		if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE), NULL_STATE);

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, 0);
			B.shape->set_shape(1, sqr(B.mp[8] = get_param(NUM_GEOMT))/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));

			if (get_param(NUM_SHEAR+1)) { //...подключаем две дополнительные экспоненциальные системы функций;
				B.shape->add_shape(CreateShapeD(SK3D_EXPP_SHAPE),	  NULL_STATE);
				B.shape->add_shape(CreateShapeD(SK3D_EXPP_SHAPE, 1), NULL_STATE);

				B.shape->init1(2, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 0);
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), fabs(B.mp[8])); 

				B.shape->init1(3, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, 0);
				B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), fabs(B.mp[8])); 
			}
		}
		else
		if ((B.type & ERR_CODE) == POLY_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_SPHERE_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 1/*draft_dim(type())*/);
 			B.shape->set_shape(0, fabs(B.mp[7]), fabs(B.mp[8]/B.mp[7]));
		}
		else {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			if (get_param(NUM_SHEAR+1)) { //...подключаем дополнительную экспоненциальную систему функций;
				B.shape->add_shape(CreateShapeD(SK3D_ZOOM_SHAPE), NULL_STATE);

				B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 0);
				B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			}
		}

/////////////////////////////////////////////////////
//...local system of coordinates and parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release();
   }

///////////////////////////////////////
//...setting parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, sqr(B.mp[8])/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), fabs(B.mp[8])); 
				B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), fabs(B.mp[8])); 
			}
			else 
			if ((B.type & ERR_CODE) == CLAYER_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, sqr(B.mp[8])/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));
			}
 			else {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), sqrt(get_param(NUM_SHEAR+1)), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			}
		}
		else
			if (id_free == OK_STATE && solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k][0][m], m);
		else
		if (id_free == NO_STATE && solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k][0][solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

/////////////////////////////////////////////////////////////
//...realization pressure (p) for fluid and brinkman current;
void CHydro3D::jump0_pressure(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, -1.);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));

	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, -1.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, -1.);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
}

////////////////////////////////////////////////////////////
//...realization of fluid pressure (p) for sphere inclusion;
void CHydro3D::jump0_sphere(double * P, int i, int m)
{
	int j = NUM_HESS+solver.id_norm; m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, -1.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));

	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, -1.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, -1.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization pressure for brinkman current arround sphere (p);
void CHydro3D::jump0_current(double * P, int i, int m)
{
	int nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(0); m += solver.id_norm;
	double RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
/////////////////////////////////////////
//...главная лапласовская часть давления;
	B[i].shape->cpy_x(0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), B[i].shape->deriv, 0., -1.);

	B[i].shape->cpy_y(0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->deriv, 0., -1.);

	B[i].shape->cpy_z(0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->deriv, 0., -1.);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), B[i].shape->FULL(B[i].shape->deriv, 0, nm), 1., RR);

	B[i].shape->cpy_y(nm, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(B[i].shape->deriv, 0, nm), 1., RR);

	B[i].shape->cpy_z(nm, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(B[i].shape->deriv, 0, nm), 1., RR);
}

////////////////////////////////////////////////
//...realization pressure (p) for darcy current;
void CHydro3D::jump0_darcy(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->adm		 (B[i].shape->p_cpy, -.5/get_param(NUM_SHEAR+2));
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vx);
void CHydro3D::jump1_fluid_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_x     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[0]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2], 0.);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[0]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[0]*G1*B[i].shape->get_R_inv();
	}
}

///////////////////////////////////////////////////////////// 
//...realization of fluid velocity (Vx) for sphere inclusion;
void CHydro3D::jump1_sphere_x(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_x     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[0]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2], 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[0]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[0]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xx(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	B[i].shape->cpy_xy(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	B[i].shape->cpy_xz(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vx);
void CHydro3D::jump1_brinkmn_x(double * P, int i, int m)
{                                                                                         
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_xx(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_xz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm	(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), G1);
		B[i].shape->adm_xx(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_xy(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_xz(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vx);
void CHydro3D::jump1_current_x(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
///////////////////////////////// 
//...лапласовская часть скорости;
	B[i].shape->cpy_xx(0, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy_xy(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy_xz(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xx(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	B[i].shape->cpy_xy(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	B[i].shape->cpy_xz(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., G1);

	B[i].shape->cpy_xx(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_xy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_xz(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., G1);

	B[i].shape->cpy_xx(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_xy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_xz(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vx);
void CHydro3D::jump1_darcy_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_x		 (B[i].shape->deriv, .5/get_param(NUM_SHEAR));
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vy);
void CHydro3D::jump1_fluid_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_y     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[1]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2], 0.);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[1]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[1]*G1*B[i].shape->get_R_inv();
	}
}

/////////////////////////////////////////////////////////////
//...realization of fluid velocity (Vy) for sphere inclusion;
void CHydro3D::jump1_sphere_y(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_y     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[1]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2], 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[1]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[1]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xy(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	B[i].shape->cpy_yy(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	B[i].shape->cpy_yz(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vy);
void CHydro3D::jump1_brinkmn_y (double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_yy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm	(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), G1);
		B[i].shape->adm_xy(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_yy(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_yz(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vy);
void CHydro3D::jump1_current_y(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	B[i].shape->cpy_xy(0, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy_yy(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy_yz(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xy(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	B[i].shape->cpy_yy(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	B[i].shape->cpy_yz(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., G1);

	B[i].shape->cpy_xy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_yy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_yz(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., G1);

	B[i].shape->cpy_xy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_yy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_yz(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vy);
void CHydro3D::jump1_darcy_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_y		 (B[i].shape->deriv, .5/get_param(NUM_SHEAR));
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
}

////////////////////////////////////////
//...realization of fluid velocity (Vz);
void CHydro3D::jump1_fluid_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_z     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[2]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 2));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1], 0.);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[2]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[2]*G1*B[i].shape->get_R_inv();
	}
}

/////////////////////////////////////////////////////////////
//...realization of fluid velocity (Vz) for sphere inclusion;
void CHydro3D::jump1_sphere_z(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_z     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, 1., -P[2]);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, G1, 0.);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0], 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1], 0.);
	
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[2]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[2]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xz(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	B[i].shape->cpy_yz(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	B[i].shape->cpy_zz(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

///////////////////////////////////////////
//...realization of brinkman velocity (Vz);
void CHydro3D::jump1_brinkmn_z (double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_xz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_zz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., alpha);
	}

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm	(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), G1);
		B[i].shape->adm_xz(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_yz(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_zz(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

///////////////////////////////////////////////////////////////////
//...realization velocity for brinkman current arround sphere (Vz);
void CHydro3D::jump1_current_z(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	B[i].shape->cpy_xz(0, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy_yz(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy_zz(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_xz(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	B[i].shape->cpy_yz(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	B[i].shape->cpy_zz(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., G1);

	B[i].shape->cpy_xz(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_yz(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_zz(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., G1);

	B[i].shape->cpy_xz(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_yz(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_zz(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////
//...realization of darcy velocity (Vz);
void CHydro3D::jump1_darcy_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->adm_z		 (B[i].shape->deriv, .5/get_param(NUM_SHEAR));
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVx);
void CHydro3D::jump2_fluid_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_xx	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[4]);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[5]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*G1);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, -P[4]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, -P[5]*G1);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[3]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[3]*G1*B[i].shape->get_R_inv();
	}
}

//////////////////////////////////////////////////////////////////////////////////// 
//...realization of normal derivative of fluid velocity (dnVx) for sphere inclusion;
void CHydro3D::jump2_sphere_x(double * P, int i, int m)
{
	int num_hess = NUM_HESS+solver.id_norm, j = num_hess+2, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_xx	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[4]);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[5]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, -P[4]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, -P[5]*G1);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[3]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[3]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0); //...deriv_N_xx;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0); //...deriv_N_xz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVx);
void CHydro3D::jump2_brinkmn_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1); real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[3]*G1);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[4]*G1);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[5]*G1);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),		B[i].shape->FULL(solver.hh[i][0][num_hess+2], l),	 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l),	 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVx);
void CHydro3D::jump2_current_x(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), num_hess = NUM_HESS+solver.id_norm, k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
///////////////////////////////// 
//...лапласовская часть скорости;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0),    B[i].shape->FULL(solver.hh[i][0][num_hess], 0), 0., 1.); //...deriv_N_xx;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., 1.); //...deriv_N_xy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0), 0., 1.); //...deriv_N_xz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0),    B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0), alpha, -alpha*RR); //...deriv_N_xx;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1), alpha, -alpha*RR); //...deriv_N_xy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0), alpha, -alpha*RR); //...deriv_N_xz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy_x(nm+1, B[i].shape->deriv); B[i].shape->admittance(nm+1, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+1, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+1, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+4], 0); //...deriv_N_xx;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+4], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+5], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm+2, B[i].shape->deriv); B[i].shape->admittance(nm+2, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+2, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+2, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+6], 0); //...deriv_N_xx;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+6], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+7], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVy);
void CHydro3D::jump2_fluid_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_yy	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[5]); 

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.); 
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*G1);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, -P[3]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, -P[5]*G1);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[4]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[4]*G1*B[i].shape->get_R_inv();
	}
}

////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVy) for sphere inclusion;
void CHydro3D::jump2_sphere_y(double * P, int i, int m)
{
	int num_hess = NUM_HESS+solver.id_norm, j = num_hess+2, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_yy	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[5]); 

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.); 
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*G1);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, -P[3]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, -P[5]*G1);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[4]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[4]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1); //...deriv_N_xy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2); //...deriv_N_yy;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVy);
void CHydro3D::jump2_brinkmn_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1); real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[3]*G1);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[4]*G1);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[5]*G1);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),		B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 2), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVy);
void CHydro3D::jump2_current_y(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), num_hess = NUM_HESS+solver.id_norm, k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0),    B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., 1.); //...deriv_N_xy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2), 0., 1.); //...deriv_N_yy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., 1.); //...deriv_N_yz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 1), alpha, -alpha*RR); //...deriv_N_xy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], 0, 2), alpha, -alpha*RR); //...deriv_N_yy;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1), alpha, -alpha*RR); //...deriv_N_yz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy_x(nm+1, B[i].shape->deriv); B[i].shape->admittance(nm+1, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+1, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+1, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+4], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+4], 0, 2); //...deriv_N_yy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+5], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm+2, B[i].shape->deriv); B[i].shape->admittance(nm+2, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+2, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+2, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+6], 0, 1); //...deriv_N_xy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+6], 0, 2); //...deriv_N_yy;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+7], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

///////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVz);
void CHydro3D::jump2_fluid_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy_zz	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[5], 0.);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*G1);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->p_cpy, -P[3]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, -P[4]*G1);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
	}
	if (B[i].shape->get_N() && regul) {
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[5]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[5]*G1*B[i].shape->get_R_inv();
	}
}

////////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of fluid velocity (dnVz) for sphere inclusion;
void CHydro3D::jump2_sphere_z(double * P, int i, int m)
{
	int num_hess = NUM_HESS+solver.id_norm, j = num_hess+2, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy_zz	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[5], 0.);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2]*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->p_cpy, -P[3]*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, -P[4]*G1);
	
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[5]*G1*B[i].shape->get_R_inv();
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[5]*G1*B[i].shape->get_R_inv();
	}

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0); //...deriv_N_xz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1); //...deriv_N_yz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2); //...deriv_N_zz;
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/(2.*k+3.);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////
//...realization of normal derivative of brinkman velocity (dnVz);
void CHydro3D::jump2_brinkmn_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1); real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[3]*G1);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[4]*G1);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[5]*G1);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),		B[i].shape->FULL(solver.hh[i][0][num_hess+3], l),	 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 2), 1., -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization normal derivative for brinkman current arround sphere (dnVz);
void CHydro3D::jump2_current_z(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), num_hess = NUM_HESS+solver.id_norm, k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0),    B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0), 0., 1.); //...deriv_N_xz;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., 1.); //...deriv_N_yz;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2), 0., 1.); //...deriv_N_zz;

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0), alpha, -alpha*RR); //...deriv_N_xz;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 1), alpha, -alpha*RR); //...deriv_N_yz;
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], 0, 2), alpha, -alpha*RR); //...deriv_N_zz;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy_x(nm+1, B[i].shape->deriv); B[i].shape->admittance(nm+1, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+1, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+1, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+5], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+5], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+5], 0, 2); //...deriv_N_zz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm+2, B[i].shape->deriv); B[i].shape->admittance(nm+2, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_y(nm+2, B[i].shape->deriv, P[4]);
	B[i].shape->adm_z(nm+2, B[i].shape->deriv, P[5]); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., G1);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+7], 0); //...deriv_N_xz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+7], 0, 1); //...deriv_N_yz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	p = B[i].shape->FULL(solver.hh[i][0][num_hess+7], 0, 2); //...deriv_N_zz;
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////////////////////
//...realization biharmonic potential of fluid velocity;
void CHydro3D::jump3_fluid(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));

	B[i].shape->cpy       (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2]*G1, 0.);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
}

/////////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential of fluid velocity for sphere inclusion;
void CHydro3D::jump3_sphere(double * P, int i, int m)
{
	int num_hess = NUM_HESS+solver.id_norm, j = num_hess+2, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = .5/get_param(NUM_SHEAR), RR = sqr(B[i].mp[8]); real_T * p;
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0]*G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));

	B[i].shape->cpy       (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*G1, 0.);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[2]*G1, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));
	B[i].shape->cpy(1, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0), p, 1., -G1);

	B[i].shape->cpy_y(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 1), p, 1., -G1);

	B[i].shape->cpy_z(1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, 1);
	for (int k = 0; k <= N_mpl; k++)
	for (int l = 0; l <= k*2; l++) p[k*k+l] *= RR/((k+4.)*(2.*k+3.));
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][j], 0, 2), p, 1., -G1);

	add_collocation(i, m, j);
}

////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman velocity;
void CHydro3D::jump3_brinkmn(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1);
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_x(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_y(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., alpha);

		B[i].shape->cpy_z(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., alpha);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), -alpha);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), -alpha);
	}
}

//////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman current arround sphere;
void CHydro3D::jump3_current(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1), 
			 RR = get_param(NUM_MPLS+1)*fabs(B[i].mp[7])/fabs(B[i].mp[8]);
	real_T * p;

///////////////////////////////// 
//...лапласовская часть скорости;
	B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy_z(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, alpha, -alpha*RR);

	B[i].shape->cpy_y(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, alpha, -alpha*RR);

	B[i].shape->cpy_z(nm, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm);
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, alpha, -alpha*RR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->cpy_x(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_y(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_z(nm+1, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->cpy_x(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., -alpha);

	B[i].shape->cpy_y(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., -alpha);

	B[i].shape->cpy_z(nm+2, B[i].shape->deriv); p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., -alpha);
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (px);
void CHydro3D::jump4_fluid_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->cpy_xx	 (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[4]);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[5]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv,	P[3]);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1], 0.);
	//B[i].shape->adm_y(B[i].shape->p_cpy, P[3]);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2], 0.);
	//B[i].shape->adm_z(B[i].shape->deriv, P[3]);
	if (B[i].shape->get_N() && regul) {
		B[i].shape->p_cpy[3] = 0.;
		B[i].shape->deriv[3] = 0.;

		B[i].shape->p_cpy[2] = -P[3]*.5*B[i].shape->get_R_inv();
		B[i].shape->deriv[1] =  P[3]*.5*B[i].shape->get_R_inv();
	}
	B[i].shape->adm_y(B[i].shape->p_cpy, P[3]);
	B[i].shape->adm_z(B[i].shape->deriv, P[3]);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (px);
void CHydro3D::jump4_brinkmn_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double alpha = 1./get_param(NUM_SHEAR+1);	real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);
		B[i].shape->adm_x		 (l, B[i].shape->FULL(solver.hh[i][0][m], l), P[3]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);
		B[i].shape->adm_y		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), P[3]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
		B[i].shape->adm_z		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), P[3]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[3]);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[4]);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l), nm, 0, -nm+l), P[5]);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),	 B[i].shape->FULL(solver.hh[i][0][num_hess+2], l),	 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l),	 1., -alpha);
	}
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (py);
void CHydro3D::jump4_fluid_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->cpy_yy	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[5]); 

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1], 0.); 
	B[i].shape->adm_y     (B[i].shape->deriv,	P[4]);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0], 0.);
	//B[i].shape->adm_x(B[i].shape->p_cpy, P[4]);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2], 0.);
	//B[i].shape->adm_z(B[i].shape->deriv, P[4]);
	if (B[i].shape->get_N() && regul) {
		B[i].shape->p_cpy[2] = 0.;
		B[i].shape->deriv[2] = 0.;

		B[i].shape->p_cpy[3] =  P[4]*.5*B[i].shape->get_R_inv();
		B[i].shape->deriv[1] = -P[4]*.5*B[i].shape->get_R_inv();
	}
	B[i].shape->adm_x(B[i].shape->p_cpy, P[4]);
	B[i].shape->adm_z(B[i].shape->deriv, P[4]);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (py);
void CHydro3D::jump4_brinkmn_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1); real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);
		B[i].shape->adm_x		 (l, B[i].shape->FULL(solver.hh[i][0][m], l), P[4]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);
		B[i].shape->adm_y		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), P[4]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
		B[i].shape->adm_z		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), P[4]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[3]*G1);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[4]*G1);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), P[5]*G1);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),	 B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 2), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
	}
}

////////////////////////////////////////////////////////////
//...realization of surface tension for fluid velocity (pz);
void CHydro3D::jump4_fluid_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	B[i].shape->cpy_zz	(B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[5], 0.);
	B[i].shape->adm_xz    (B[i].shape->deriv, P[3]);
	B[i].shape->adm_yz    (B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[2], 0.);
	B[i].shape->adm_z     (B[i].shape->deriv, P[5]);
	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));

	B[i].shape->cpy       (B[i].shape->p_cpy, B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[0], 0.);
	//B[i].shape->adm_x(B[i].shape->p_cpy, P[5]);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1], 0.);
	//B[i].shape->adm_y(B[i].shape->deriv, P[5]);
	if (B[i].shape->get_N() && regul) {
		B[i].shape->p_cpy[1] = 0.;
		B[i].shape->deriv[1] = 0.;

		B[i].shape->p_cpy[3] = -P[5]*.5*B[i].shape->get_R_inv();
		B[i].shape->deriv[2] =  P[5]*.5*B[i].shape->get_R_inv();
	}
	B[i].shape->adm_x(B[i].shape->p_cpy, P[5]);
	B[i].shape->adm_y(B[i].shape->deriv, P[5]);
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));
	}
}

///////////////////////////////////////////////////////////////
//...realization of surface tension for brinkman velocity (pz);
void CHydro3D::jump4_brinkmn_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR), alpha = G1/get_param(NUM_SHEAR+1); real_T * p;
/////////////////////////////////
//...лапласовская часть скорости;
	for (l = 0; l < nm; l++) {
		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), p, 0., alpha);
		B[i].shape->adm_x		 (l, B[i].shape->FULL(solver.hh[i][0][m], l), P[5]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), p, 0., alpha);
		B[i].shape->adm_y		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), P[5]);

		p = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), p, 0., alpha);
		B[i].shape->adm_z		 (l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), P[5]);
	}
/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->adm_x(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[3]*G1);
		B[i].shape->adm_y(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[4]*G1);
		B[i].shape->adm_z(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), P[5]*G1);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),		B[i].shape->FULL(solver.hh[i][0][num_hess+3], l),	 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 1), 1., -alpha);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l, 2), 1., -alpha);
	}
}

////////////////////////////////////////////////////
//...realization primitive of vector potential (fx);
void CHydro3D::jump5_fluid_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));

	if (B[i].shape->get_N() && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[0]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[0]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for sphere inclusion;
void CHydro3D::jump5_sphere_x(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);
	
	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, -1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[3] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[3] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] = -P[0]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] =  P[0]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for brinkman velocity;
void CHydro3D::jump5_brinkmn_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->primitive (nm+l, B[i].shape->deriv);
		B[i].shape->admittance(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 0), nm, 0, -nm+l), B[i].shape->deriv, 0., G1);
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fx) for brinkman current arround sphere;
void CHydro3D::jump5_current_x(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR); real_T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->primitive(nm+1, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->primitive(nm+2, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0), p, 1., G1);
}

////////////////////////////////////////////////////
//...realization primitive of vector potential (fy);
void CHydro3D::jump5_fluid_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

	if (B[i].shape->get_N(0) && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[1]*.5*B[i].shape->get_R_inv(0)*P[2]*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[1]*.5*B[i].shape->get_R_inv(0)*P[2]*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for sphere inclusion;
void CHydro3D::jump5_sphere_y(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 0));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[2] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[2] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] =  P[1]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 2)[1] = -P[1]*.5*B[i].shape->get_R_inv()*P[2]*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for brinkman velocity;
void CHydro3D::jump5_brinkmn_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->primitive (nm+l, B[i].shape->deriv);
		B[i].shape->admittance(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 1), nm, 0, -nm+l), B[i].shape->deriv, 0., G1);
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fy) for brinkman current arround sphere;
void CHydro3D::jump5_current_y(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR); real_T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->primitive(nm+1, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->primitive(nm+2, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), p, 1., G1);
}

///////////////////////////////////////////////////
//...realization primitiveof vector potential (fz);
void CHydro3D::jump5_fluid_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	for (l = 0; l < nm; l++)
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 2));

	if (B[i].shape->get_N() && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[2]*.5*B[i].shape->get_R_inv()*P[2]*.5*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[2]*.5*B[i].shape->get_R_inv()*P[2]*.5*P[5];
	}
}

/////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fz) for sphere inclusion;
void CHydro3D::jump5_sphere_z(double * P, int i, int m)
{
	int j = NUM_HESS+2+solver.id_norm, N_mpl = B[i].shape->get_N(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

	B[i].shape->primitive (B[i].shape->deriv);
	B[i].shape->admittance(B[i].shape->deriv, NULL, G1, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 2));
	B[i].shape->cpy(1, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][j], 0, 1));

	if (N_mpl && regul) { //...регуляризация скорости;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[1] = 0.;
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[1] = 0.;

		B[i].shape->FULL(solver.hh[i][0][m], 0, 0)[3] = -P[2]*.5*B[i].shape->get_R_inv()*P[2]*.5*P[5];
		B[i].shape->FULL(solver.hh[i][0][m], 0, 1)[2] =  P[2]*.5*B[i].shape->get_R_inv()*P[2]*.5*P[5];
	}
	add_collocation(i, m, j);
}

//////////////////////////////////////////////////////////////////////////
//...realization primitive of vector potential (fz) for brinkman velocity;
void CHydro3D::jump5_brinkmn_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR);

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	for (l = 0; l < nm; l++) {
		B[i].shape->primitive (nm+l, B[i].shape->deriv);
		B[i].shape->admittance(nm+l, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][m], nm+l, 2), nm, 0, -nm+l), B[i].shape->deriv, 0., G1);
	}
}

//////////////////////////////////////////////////////////////////////////
//...realization biharmonic potential for brinkman current arround sphere;
void CHydro3D::jump5_current_z(double * P, int i, int m)
{
	int    nm = B[i].shape->num_usual(), N_mpl = B[i].shape->get_N(), k, j; m += solver.id_norm;
	double G1 = 1./get_param(NUM_SHEAR); real_T * p;

/////////////////////////////////////
//...гельмгольцевская часть скорости;
	B[i].shape->primitive(nm+1, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+1);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[0][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 0., G1);

/////////////////////////////////////////// 
//...дополнительное слагаемое к потенциалу;
	B[i].shape->primitive(nm+2, B[i].shape->deriv);	p = B[i].shape->FULL(B[i].shape->deriv, 0, nm+2);
	for (k = 0; k <= N_mpl; k++)
	for (j = 0; j <= k*2;   j++) p[k*k+j] *= TN[1][k];
	B[i].shape->admittance(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 2), p, 1., G1);
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
void CHydro3D::hessian_deriv_N(double * P, int i)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->cpy_xx(); B[i].shape->deriv_N(); 
	B[i].shape->cpy_xx();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], l));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+2], nm+l), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+4], nm*2+l), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+6], nm*3+l), nm*3, 0, -nm*3+l));
	}
	B[i].shape->cpy_xy(); B[i].shape->deriv_N();
	B[i].shape->cpy_xy();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+2], nm+l, 1), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+4], nm*2+l, 1), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+6], nm*3+l, 1), nm*3, 0, -nm*3+l));
	}
	B[i].shape->cpy_yy(); B[i].shape->deriv_N();
	B[i].shape->cpy_yy();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+2], nm+l, 2), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+4], nm*2+l, 2), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+6], nm*3+l, 2), nm*3, 0, -nm*3+l));
	}
	B[i].shape->cpy_xz(); B[i].shape->deriv_N();
	B[i].shape->cpy_xz();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], l));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+3], nm+l), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+5], nm*2+l), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+7], nm*3+l), nm*3, 0, -nm*3+l));
	}
	B[i].shape->cpy_yz(); B[i].shape->deriv_N();
	B[i].shape->cpy_yz();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+3], nm+l, 1), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+5], nm*2+l, 1), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+7], nm*3+l, 1), nm*3, 0, -nm*3+l));
	}
	B[i].shape->cpy_zz(); B[i].shape->deriv_N();
	B[i].shape->cpy_zz();
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2));
		B[i].shape->cpy(nm+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+3], nm+l, 2), nm, 0, -nm+l));
		B[i].shape->cpy(nm*2+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+5], nm*2+l, 2), nm*2, 0, -nm*2+l));
		B[i].shape->cpy(nm*3+l, B[i].shape->deriv, B[i].shape->FULL(B[i].shape->FULL(solver.hh[i][0][num_hess+7], nm*3+l, 2), nm*3, 0, -nm*3+l));
	}
}

///////////////////////////////////////////////////////////
//...realization normal derivatives of hessian for dim = 0;
void CHydro3D::hessian_deriv_N(int k, double * P, int i)
{
	int num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->cpy_xx(); B[i].shape->deriv_N(k); 
	B[i].shape->cpy_xx();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, -k));

	B[i].shape->cpy_xy(); B[i].shape->deriv_N(k);
	B[i].shape->cpy_xy();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1-k));

	B[i].shape->cpy_yy(); B[i].shape->deriv_N(k);
	B[i].shape->cpy_yy();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2-k));

	B[i].shape->cpy_xz(); B[i].shape->deriv_N(k);
	B[i].shape->cpy_xz();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, -k));

	B[i].shape->cpy_yz(); B[i].shape->deriv_N(k);
	B[i].shape->cpy_yz();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1-k));

	B[i].shape->cpy_zz(); B[i].shape->deriv_N(k);
	B[i].shape->cpy_zz();
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2-k));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CHydro3D::jump_make_local(int i, int m)
{
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		real_T P[] = { solver.hh[i][0][m  ][j], 
							 solver.hh[i][0][m+1][j], 
							 solver.hh[i][0][m+2][j] };
		B[i].shape->norm_local_T(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
		solver.hh[i][0][m+2][j] = P[2];
	}
}

void CHydro3D::jump_make_common(int i, int m)
{
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		real_T P[] = { solver.hh[i][0][m  ][j], 
							 solver.hh[i][0][m+1][j], 
							 solver.hh[i][0][m+2][j] };
		B[i].shape->norm_common_T(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
		solver.hh[i][0][m+2][j] = P[2];
	}
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для фазы свободного течения жидкости;
int CHydro3D::gram1_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, hz, p4, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);
		
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

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
				hy = nd->nY[l]*hx;
				hz = nd->nZ[l]*hx;
				hx *= nd->nX[l];
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////////////////////////
//...jump of all displacement moments and composition functional;
			B[i].shape->parametrization_hess(P, 1);
			if (p4 == SKEWS_BND-SPECIAL_BND) { //...краевые условия для скорости поперек потока на ячейке;
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				jump2_x(P, i, 0); solver.admittance(i, 0, G1); solver.admittance (i, 3, 0., 0, P[3]);
				jump2_y(P, i, 1); solver.admittance(i, 1, G1); solver.admittance (i, 3, 1., 1, P[4]);
				jump2_z(P, i, 2); solver.admittance(i, 2, G1); solver.admittance (i, 3, 1., 2, P[5]);
				jump_make_common(i, 0);
				solver.admittance (i, 0, 1., 3, -nd->nX[l]);
				solver.admittance (i, 1, 1., 3, -nd->nY[l]);
				solver.admittance (i, 2, 1., 3, -nd->nZ[l]);
				jump1_x(P, i, 4); solver.admittance (i, 3, 0., 4, P[3]);
				jump1_y(P, i, 5); solver.admittance (i, 3, 1., 5, P[4]);
				jump1_z(P, i, 6); solver.admittance (i, 3, 1., 6, P[5]);

				solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
				solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
				solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
			}
			else
			if (p4 == MAX_HIT || p4 == NUMS_BND) { //...краевые условия для давления вдоль потока на ячейке;
				jump0(P, i, 3);
 				solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
				
				if (fabs(hx) > EE) 
					solver.to_equationHH(i, 0, solver.hh[i][0][m+3], hx*f);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump1_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_y(P, i, 1); solver.admittance(i, 1, G1);
				jump1_z(P, i, 2); solver.admittance(i, 2, G1); jump_make_common(i, 0);

				solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
				solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
				
				if (fabs(hx) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m], hx*(f *= G1));
				if (fabs(hy) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*f);
				if (fabs(hz) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+2], hz*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


//////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для фазы фильтрационного течения жидкости;
int CHydro3D::gram1_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
		double hx, p4, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);
		
////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == 0.) hx = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////////////////////////
//...jump of all displacement moments and composition functional;
			B[i].shape->parametrization_hess(P, 1);
			jump1_x(P, i, 0); solver.admittance (i, 3, 0., 0, P[3]);
			jump1_y(P, i, 1); solver.admittance (i, 3, 1., 1, P[4]);
			jump1_z(P, i, 2); solver.admittance (i, 3, 1., 2, P[5]);
			solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
			if (fabs(hx) > EE) 
				solver.to_equationHH(i, 0, solver.hh[i][0][m+3], hx*f);

		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
int CHydro3D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz, 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1;
      int id_isolated = 1, m  = solver.id_norm, id_dir, k; int j;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 0);

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
				TX = TY = TZ = hz = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hz = -AZ; break;
					case 6: TZ = -AZ; hz =  AZ; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY;
				B[k].mp[3] -= TZ; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				
				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);

				solver.admittance(i, 0, g1, 4, g2); solver.admittance(i, 4, g0, 0, f1); 
				solver.admittance(i, 1, g1, 5, g2); solver.admittance(i, 5, g0, 1, f1); 
				solver.admittance(i, 2, g1, 6, g2); solver.admittance(i, 6, g0, 2, f1);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				
				B[k].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6); jump_make_common(k, 4);

				solver.admittance(k, 0, g1, 4, g2); solver.admittance(k, 4, g0, 0, f1); 
				solver.admittance(k, 1, g1, 5, g2); solver.admittance(k, 5, g0, 1, f1); 
				solver.admittance(k, 2, g1, 6, g2); solver.admittance(k, 6, g0, 2, f1);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+6], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+6], solver.hh[k][0][m+6], f);

				if (abs(id_dir) < 5) { //...давление с учетом регуляризации;
					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				}
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				if (fabs(hz) > EE) {
  				  solver.to_equationHH(i, 0, solver.hh[i][0][m+3], -hz*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m+3],  hz*f);
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
			}
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для периодической задачи (один блок);
int CHydro3D::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
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
			if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
			
			jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
			jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);

			B[i].shape->make_common(P);
			B[i].shape->norm_common(P+3);

			if (id_dir == 1) B[i].mp[1] -= AX; else
			if (id_dir == 3) B[i].mp[2] -= AY; else
			if (id_dir == 5) B[i].mp[3] -= AZ; 
			B[i].shape->set_local_P0(B[i].mp+1);
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			B[i].shape->parametrization_hess(P, 1);
			if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			jump1_x(P, i,  7); jump1_y(P, i,  8); jump1_z(P, i,  9); jump_make_common(i,  7); jump0(P, i, 10);
			jump2_x(P, i, 11); jump2_y(P, i, 12); jump2_z(P, i, 13); jump_make_common(i, 11);

			solver.admittance(i, 0, 1.,  7, -1.); solver.admittance(i,  4, 1., 11, -1.);
			solver.admittance(i, 1, 1.,  8, -1.); solver.admittance(i,  5, 1., 12, -1.);
			solver.admittance(i, 2, 1.,  9, -1.); solver.admittance(i,  6, 1., 13, -1.); 
			solver.admittance(i, 3, 1., 10, -1.); solver.admittance(i, 10, 2.,  3,  1.);

			solver.admittance(i, 0, g1,  4, g2);  solver.admittance(i,  4, g0, 0, f1); 
			solver.admittance(i, 1, g1,  5, g2);  solver.admittance(i,  5, g0, 1, f1); 
			solver.admittance(i, 2, g1,  6, g2);  solver.admittance(i,  6, g0, 2, f1); 

/////////////////////////////////////////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов (скорость, давление, производная);
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
			solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);

			if (id_dir == 5) {//...правая часть -- скачок давления и регуляризация матрицы;
				solver.to_equationHH(i, 0,	solver.hh[i][0][m+3], AZ*f);
				solver.to_equationDD(i, solver.hh[i][0][m+10], solver.hh[i][0][m+10], f);
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

////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе свободное - свободное течение жидкости;
int CHydro3D::transfer1_phase1(CGrid * nd, int i, int k, int id_local)
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
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 1);

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
					//solver.hh[i][0].clear_row_matrix(num, solver.dim[i]);
					//solver.hh[k][0].clear_row_matrix(num, solver.dim[k]);
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);
				
				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6); jump_make_common(i, 4);	

				solver.admittance (i, 0, g1, 4, g2); solver.admittance(i, 4, g0, 0, f0); 
				solver.admittance (i, 1, g1, 5, g2); solver.admittance(i, 5, g0, 1, f0); 
				solver.admittance (i, 2, g1, 6, g2); solver.admittance(i, 6, g0, 2, f0);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				
				B[k].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0);	jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6); jump_make_common(k, 4);

				solver.admittance (k, 0, g1, 4, g2); solver.admittance(k, 4, g0, 0, f0); 
				solver.admittance (k, 1, g1, 5, g2); solver.admittance(k, 5, g0, 1, f0); 
				solver.admittance (k, 2, g1, 6, g2); solver.admittance(k, 6, g0, 2, f0); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+6], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+6], solver.hh[k][0][m+6], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


//////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе фильтрационное - фильтрационное течение жидкости;
int CHydro3D::transfer1_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      int m = solver.id_norm;
      double f, P[6];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.001, 10., 20., 30., AXIS_Z, 1);

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
				jump1_x(P, i, 0); solver.admittance (i, 3, 0., 0, P[3]);
				jump1_y(P, i, 1); solver.admittance (i, 3, 1., 1, P[4]);
				jump1_z(P, i, 2); solver.admittance (i, 3, 1., 2, P[5]);
				
				B[k].shape->parametrization_hess(P, 1);
				jump0(P, k, 3);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе свободное - фильтрационное течение жидкости;
int CHydro3D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double f, P[6], beta = 0.1;
      int m = solver.id_norm;

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

///////////////////////////////////////////////
//...jump of all neaded moments for this block;
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				if (B[i].link[NUM_PHASE] == -1) {
					jump1_x(P, i, 0); solver.admittance (i, 3, 0., 0, P[3]);
					jump1_y(P, i, 1); solver.admittance (i, 3, 1., 1, P[4]);
					jump1_z(P, i, 2); solver.admittance (i, 3, 1., 2, P[5]);
					solver.admittance (i, 0, 1., 3, -P[3]);
					solver.admittance (i, 1, 1., 3, -P[4]);
					solver.admittance (i, 2, 1., 3, -P[5]);

					jump4_x(P, i, 4); solver.admittance (i, 7, 0., 4, P[3]);
					jump4_y(P, i, 5); solver.admittance (i, 7, 1., 5, P[4]);
					jump4_z(P, i, 6); solver.admittance (i, 7, 1., 6, P[5]);
					solver.admittance (i, 4, 1., 7, -P[3]); solver.admittance (i, 4, 1., 0, beta);
					solver.admittance (i, 5, 1., 7, -P[4]); solver.admittance (i, 5, 1., 1, beta);
					solver.admittance (i, 6, 1., 7, -P[5]); solver.admittance (i, 6, 1., 2, beta);
					jump0(P, i, 7);
				}
				if (B[i].link[NUM_PHASE] == -2) {
					jump1_x(P, i, 0); solver.admittance (i, 3, 0., 0, P[3]);
					jump1_y(P, i, 1); solver.admittance (i, 3, 1., 1, P[4]);
					jump1_z(P, i, 2); solver.admittance (i, 3, 1., 2, P[5]);
					jump0  (P, i, 7); 
					solver.admittance (i, 0, beta, 3, -P[3]*beta);
					solver.admittance (i, 1, beta, 3, -P[4]*beta);
					solver.admittance (i, 2, beta, 3, -P[5]*beta);
				}

///////////////////////////////////////////////////////
//...jump of all neaded moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
	
				B[k].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				if (B[k].link[NUM_PHASE] == -1) {
					jump1_x(P, k, 0); solver.admittance (k, 3, 0., 0, P[3]);
					jump1_y(P, k, 1); solver.admittance (k, 3, 1., 1, P[4]);
					jump1_z(P, k, 2); solver.admittance (k, 3, 1., 2, P[5]);
					solver.admittance (k, 0, 1., 3, -P[3]);
					solver.admittance (k, 1, 1., 3, -P[4]);
					solver.admittance (k, 2, 1., 3, -P[5]);

					jump4_x(P, k, 4); solver.admittance (k, 7, 0., 4, P[3]);
					jump4_y(P, k, 5); solver.admittance (k, 7, 1., 5, P[4]);
					jump4_z(P, k, 6); solver.admittance (k, 7, 1., 6, P[5]);
					solver.admittance (k, 4, 1., 7, -P[3]); solver.admittance (k, 4, 1., 0, beta);
					solver.admittance (k, 5, 1., 7, -P[4]); solver.admittance (k, 5, 1., 1, beta);
					solver.admittance (k, 6, 1., 7, -P[5]); solver.admittance (k, 6, 1., 2, beta);
					jump0(P, k, 7);
				}
				if (B[k].link[NUM_PHASE] == -2) {
					jump1_x(P, k, 0); solver.admittance (k, 3, 0., 0, P[3]);
					jump1_y(P, k, 1); solver.admittance (k, 3, 1., 1, P[4]);
					jump1_z(P, k, 2); solver.admittance (k, 3, 1., 2, P[5]);
					jump0	 (P, k, 7); 
					solver.admittance (k, 0, beta, 3, -P[3]*beta);
					solver.admittance (k, 1, beta, 3, -P[4]*beta);
					solver.admittance (k, 2, beta, 3, -P[5]*beta);
				}

//////////////////////////////////////////////////////////////
//...сшивка скоростей и давлений методом наименьших квадратов;
				if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -2) {
					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m],   f);
					solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+1], f);
					solver.to_transferTR(i, j, solver.hh[i][0][m+6], solver.hh[k][0][m+2], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+7], solver.hh[k][0][m+7], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+7], solver.hh[k][0][m+7], f);
				}
				if (B[k].link[NUM_PHASE] == -1 && B[i].link[NUM_PHASE] == -2) {
					solver.to_transferTR(i, j, solver.hh[i][0][m+7], solver.hh[k][0][m+7], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+5], solver.hh[k][0][m+5], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+6], solver.hh[k][0][m+6], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+4], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+5], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+6], f);
				}
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
int CHydro3D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, AZ, f, P[6], TX, TY, TZ, hz;
      int	 m  = solver.id_norm, id_dir, k; int j;

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
				TX = TY = TZ = hz = 0.;
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
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
				solver.admittance(i, 0, G1); solver.admittance(i, 1, G1); solver.admittance(i, 2, G1); 
				solver.admittance(i, 7, 0., 0, nd->nX[l]); 
				solver.admittance(i, 7, 1., 1, nd->nY[l]); 
				solver.admittance(i, 7, 1., 2, nd->nZ[l]);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6);	jump_make_common(k, 4);
				solver.admittance(k, 0, G1); solver.admittance(k, 1, G1); solver.admittance(k, 2, G1); 
				solver.admittance(k, 7, 0., 0, nd->nX[l]); 
				solver.admittance(k, 7, 1., 1, nd->nY[l]); 
				solver.admittance(k, 7, 1., 2, nd->nZ[l]);

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

				solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				if (fabs(hz) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+3], -hz*f);
				  solver.to_equationHL(k, 0, solver.hh[k][0][m+3],  hz*f);
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+4],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+4], -f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+5],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+5], -f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+6],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+6], -f);

				solver.to_equationER(i, solver.hh[i][0][m+7], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+7], solver.hh[k][0][m+3],  f);

////////////////////////////////////
//...естественные граничные условия;
				if (fabs(hz) > EE) {
					solver.to_equationEH(i, 0, solver.hh[i][0][m+7],  hz*f);
					solver.to_equationEH(k, 0, solver.hh[k][0][m+7], -hz*f);
				}				
				B[k].mp[1] += TX;
				B[k].mp[2] += TY;
				B[k].mp[3] += TZ; B[k].shape->set_local_P0(B[k].mp+1);
			}
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для фазы свободного течения жидкости;
int CHydro3D::gram4_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
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

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
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
			if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
			jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
			solver.admittance(i, 0, G1); solver.admittance(i, 1, G1); solver.admittance(i, 2, G1); 
			solver.admittance(i, 7, 0., 0, nd->nX[l]); 
			solver.admittance(i, 7, 1., 1, nd->nY[l]); 
			solver.admittance(i, 7, 1., 2, nd->nZ[l]);

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p4 == MAX_HIT || p4 == NUMS_BND) {
				solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+3], f);
				
				if (fabs(hx) > EE) 
					solver.to_equationHH(i, 0, solver.hh[i][0][m+3], hx*f);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
				solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);

				if (fabs(hx) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m],   hx*G1*f);
				if (fabs(hy) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*G1*f);
				if (fabs(hz) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+2], hz*G1*f);
			}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+4],  f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+5],  f);
			solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+6],  f);
			solver.to_equationEE(i, solver.hh[i][0][m+7], solver.hh[i][0][m+3], -f);

////////////////////////////////////
//...естественные граничные условия;
			if ((p4 == MAX_HIT || p4 == NUMS_BND) && fabs(hx) > EE)
				solver.to_equationEH(i, 0, solver.hh[i][0][m+7], -hx*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для фазы фильтрационного течения жидкости;
int CHydro3D::gram4_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < (int)solver.N && B[i].shape && B[i].mp) {
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

			if (p4 == NUMS_BND) hx = -nd->Z[l];
			else 
			if (p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
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
			if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе свободное - свободное течение жидкости;
int CHydro3D::transfer4_phase1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < (int)solver.N) {
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
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump1_x(P, i, 0); jump1_y(P, i, 1); jump1_z(P, i, 2); jump_make_common(i, 0); jump0(P, i, 3);
				jump2_x(P, i, 4); jump2_y(P, i, 5); jump2_z(P, i, 6);	jump_make_common(i, 4);
				solver.admittance(i, 0, G1); solver.admittance(i, 1, G1); solver.admittance(i, 2, G1); 
				solver.admittance(i, 7, 0., 0, nd->nX[l]); 
				solver.admittance(i, 7, 1., 1, nd->nY[l]); 
				solver.admittance(i, 7, 1., 2, nd->nZ[l]);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, k); else hessian_deriv_N(1, P, k);

				jump1_x(P, k, 0); jump1_y(P, k, 1); jump1_z(P, k, 2); jump_make_common(k, 0); jump0(P, k, 3);
				jump2_x(P, k, 4); jump2_y(P, k, 5); jump2_z(P, k, 6);	jump_make_common(k, 4);
				solver.admittance(k, 0, G1); solver.admittance(k, 1, G1); solver.admittance(k, 2, G1); 
				solver.admittance(k, 7, 0., 0, nd->nX[l]); 
				solver.admittance(k, 7, 1., 1, nd->nY[l]); 
				solver.admittance(k, 7, 1., 2, nd->nZ[l]);

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

				solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m],   solver.hh[i][0][m+4],  f);
				solver.to_equationEL(k, solver.hh[k][0][m],   solver.hh[k][0][m+4], -f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+5],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+5], -f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+6],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+6], -f);

				solver.to_equationER(i, solver.hh[i][0][m+7], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+7], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе фильтрационное - фильтрационное течение жидкости;
int CHydro3D::transfer4_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < (int)solver.N) {
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
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////
//...интегрирование скоростей фильтрации по границе блоков;
int CHydro3D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, m = solver.id_norm;
      double f, P[6], V[3];

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

///////////////////////////////////////////////////////////////////////////////////
//...вычисляем скорости в сепарированном виде для формирования интеграла по объему;
			B[i].shape->parametrization_hess(P, 1);
			jump3  (P, i, 0); 
			jump5_x(P, i, 1); //...нужно ли приведение нормали к общей системе координат???
			jump5_y(P, i, 2); 
			jump5_z(P, i, 3); 

			V[0] = to_double(B[i].shape->potential(solver.hh[i][0][m], 0));
			K[0] += V[0]*nd->nX[l]*f;
			K[1] += V[0]*nd->nY[l]*f;
			K[2] += V[0]*nd->nZ[l]*f;

			V[0] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], 0));
			V[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], 0));
			V[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+3], 0));
			B[i].shape->norm_common(V);

			K[0] += V[0]*f;//...нормаль убрали внутрь функций;
			K[1] += V[1]*f;
			K[2] += V[2]*f;
			K[3+(-B[i].link[NUM_PHASE]-1)] += nd->Z[l]*nd->nZ[l]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);}

//////////////////////////////////////////////////////////
//...интегрирование скоростей фильтрации по объему блоков;
int CHydro3D::rigidy2(CGrid * nd, int i, double *  K)
{
	if (nd) {
      int m = solver.id_norm, l;
      double f, P[6], V[3];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 0., 0., 0., AXIS_Z, 1);
/////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = -1.;
			P[1] = nd->Y[l];  P[4] = 0.;
			P[2] = nd->Z[l];  P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

///////////////////////////////////
//...вычисляем скорости фильтрации;
			B[i].shape->parametrization_hess(P, 1);

			jump1_x(P, i, 0); 
			jump1_y(P, i, 1); 
			jump1_z(P, i, 2); 

			V[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   0));
			V[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], 0));
			V[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], 0));
			B[i].shape->norm_common(V);
			
			K[0] += V[0]*f;
			K[1] += V[1]*f;
			K[2] += V[2]*f;
			K[3] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////
//...istalling parameters of doubly connected media;
void CHydro3D::set_fasa_hmg(double R0, double ll, double G_dynamic, double alpha, double kk)
{
	if (size_of_param() > NUM_GEOMT+1) {
		param[NUM_GEOMT] = R0;
		param[NUM_GEOMT+1] = R0+ll;
	}
	if (size_of_param() > NUM_SHEAR+1) {
		param[NUM_SHEAR]	 = G_dynamic;
		param[NUM_SHEAR+1] = alpha;
		param[NUM_SHEAR+2] = kk;
	}
	return;
}

///////////////////////////////////////////////////////////
//...counting header for solving linear elasticity problem;
int CHydro3D::computing_header(Num_Comput Num)
{
	int  N_elem = UnPackInts(get_param(NUM_QUAD)), id_dir, i, k, l, elem, n_rhs = 2;
	char msg[201];
	
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CHydro3D_draft: N_sm = %d, N_mpl = %d, N_elem = %d, alpha = %g", N, UnPackInts(get_param(NUM_MPLS)), N_elem, get_param(NUM_SHEAR+1));
		Message(msg);

		Message(" ");
		switch (Num){
			case	 BASIC_COMPUT: Message("Junction Blocks...");	break;
			case MAPPING_COMPUT: Message("Mapping  Blocks..."); break;
			case  PERIOD_COMPUT: Message("Periodic Blocks..."); break;
		}
		Message(" ");
	}
	if (N == 1 && (B[0].type & ERR_CODE) == CLAYER_BLOCK) {
		set_eshelby_matrix (UnPackInts(get_param(NUM_MPLS)));
		set_normaliz_vector_old(UnPackInts(get_param(NUM_MPLS)), get_param(NUM_GEOMT), sqrt(get_param(NUM_SHEAR+1)));
	}

///////////////////////
//...блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += NUM_HESS+9;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par);
		for (k = 0; k < N; k++) SkeletonBounding(B[k], par);
		for (k = 0; k < N; k++) if (B[k].link) {
			for (i = 0; i < B[k].link[0]; i++) if ((elem = geom_plink_3D(B[k], l = i, id_dir, par)) >= 0) 
			solver.add_link(k, elem);
			for (i = 0; i < B[k].link[0]; i++) if ((elem = block_plink_3D(B[k], l = i, id_dir, par)) >= 0) 
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
	for (k = 0; k < (int)solver.N; k++)
		solver.set_dimension(k, freedom_block(k));

	solver.struct_init();

	if (solver.mode(FULLY_MODE)) { 
		solver.test_struct("1", 0);
		solver.test_struct("2", 1);
	}
   return OK_STATE;
}

////////////////////////////////////////////////////////////
//...счетная схема для явного аналитического решения задачи;
int CHydro3D::computing_kernel5()
{
	set_param(NUM_MPLS, PackInts(2, 2)); //...multipoles degree;
	if (computing_header(SPECIAL3COMPUT) != OK_STATE) 
		return ERR_STATE;

///////////////////////////////////////////////////////////////
//...переносим результат в функции формы (задаем решение явно);
	shapes_init(OK_STATE);
	for (int i = 0; i < 1; i++) {
		memset(solver.hh[i][0][0], 0, solver.dim[i]*sizeof(double));
		B[i].shape->FULL(solver.hh[i][0][0], 0, 0)[6] = 0.5*sqr(B[i].shape->get_R());	
		B[i].shape->FULL(solver.hh[i][0][0], 0, 1)[5] = 0.5*sqr(B[i].shape->get_R());	
		B[i].shape->FULL(solver.hh[i][0][0], 0, 2)[0] = 0.5;	
		B[i].shape->set_potential(solver.hh[i][0][0], 0);
	}
   return OK_STATE;
}

//////////////////////////////////////////////////////////////
//...сложение коллокационных векторов с учетом матрицы Эшелби;
void CHydro3D::add_collocation(int i, int m, int j)
{
	if (TT && TH) {
		real_T * px = B[i].shape->FULL(solver.hh[i][0][j], 0, 0), 
				  * py = B[i].shape->FULL(solver.hh[i][0][j], 0, 1), 
				  * pz = B[i].shape->FULL(solver.hh[i][0][j], 0, 2), * p, 
					ff = fabs(B[i].mp[8])/B[i].shape->get_R(0), ff_inv = 1./ff;
		int   N_mpl = B[i].shape->get_N(0), nn, mm, kk, k, l, ll;
		for (k = N_mpl; k >= 0; k--) {
			nn = 2*k+1; kk = k*k; mm = nn+4; ll = (k+2)*(k+2);
			if (TT[k]) //...вычисляем текущую степень;
			for (l = 0; l < nn; l++) {				
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 0), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l]+py[kk+j]*TT[k][j+nn][l]+pz[kk+j]*TT[k][j+nn*2][l]);
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 1), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l+nn]+py[kk+j]*TT[k][j+nn][l+nn]+pz[kk+j]*TT[k][j+nn*2][l+nn]);
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 2), j = 0; j < nn; j++) p[kk+l] += ff_inv*(px[kk+j]*TT[k][j][l+nn*2]+py[kk+j]*TT[k][j+nn][l+nn*2]+pz[kk+j]*TT[k][j+nn*2][l+nn*2]);
			}
			if (TH[k] && k+2 <= N_mpl) //...правим более высокую степень;
			for (l = 0; l < mm; l++) {
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 0), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l]+py[kk+j]*TH[k][j+nn][l]+pz[kk+j]*TH[k][j+nn*2][l]);
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 1), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l+mm]+py[kk+j]*TH[k][j+nn][l+mm]+pz[kk+j]*TH[k][j+nn*2][l+mm]);
				for ( p = B[i].shape->FULL(solver.hh[i][0][m], 0, 2), j = 0; j < nn; j++) p[ll+l] += ff*(px[kk+j]*TH[k][j][l+mm*2]+py[kk+j]*TH[k][j+nn][l+mm*2]+pz[kk+j]*TH[k][j+nn*2][l+mm*2]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
//...формирование внутренних матриц в аналитическом методе для граничных функций;
void CHydro3D::intern_matrix(real_T **& T, int N, real_T **& TH)
{
	int l, m, n = N, nn = 2*n+1, mm = nn-4; 
	if (! T)   set_matrix(T, 3*nn, 3*nn);
	else		clean_matrix(T, 3*nn, 3*nn);
//...Re fw;
	double ff = (2.*n+3.)/(2.*n-1.)*.5;
	T[nn][nn] = -(T[0][0] = ff*((n-1.)*n*.5-(2.*n-1.))*.5);
	if (n >= 2) T[0][4] = T[0][3+nn] = T[nn][4+nn] = -(T[nn][3] = ff);
	if (n >= 1) T[0][2+nn*2] = -(T[nn][1+nn*2] = -ff*n*.5);
	for (m = 1; m <= N; m++) { 
		T[m*2][m*2-1+nn] = T[m*2-1][m*2-1] = T[m*2-1][m*2+nn] = -(T[m*2][m*2] = ff*((m-1.)*(2.*n-1.)+(n-m-1.)*(n-m)*.5)*.25);
		if (m+2 <= n) T[m*2][m*2+4] = T[m*2][m*2+3+nn] = T[m*2-1][m*2+4+nn] = -(T[m*2-1][m*2+3] = ff*(m+1.)*(m+2.)*.5);
		if (m+1 <= n) T[m*2][m*2+2+nn*2] = -(T[m*2-1][m*2+1+nn*2] = -ff*(m+1.)*(n+m)*.5);
	}
//...Im fw;
	for (m = 1; m <= N; m++) { 
		T[m*2+nn][m*2] = T[m*2+nn][m*2-1+nn] = T[m*2-1+nn][m*2+nn] = -(T[m*2-1+nn][m*2-1] = -ff*((n-m-1.)*(n-m)*.5-(2.*n-1.))*.25);
		if (m == 2) T[m*2+nn][0] = -(T[m*2-1+nn][nn] = ff*(2.*n-1.+(n-3.)*(n-2.)*.25)*(n-1.)*n*.125);
		if (m == 1)	{	
			double dd = ff*((n-2.)*(n-1.)*.5+(2.*n-1.))*.25;
			T[m*2+nn][m*2] += dd; T[m*2+nn][m*2-1+nn] += dd; T[m*2-1+nn][m*2-1] += dd; T[m*2-1+nn][m*2+nn] -= dd;
			T[m*2+nn][nn*2] = ff*(1.-n*n)*n*.25;
		}
		if (m > 2) T[m*2+nn][m*2-5+nn] = T[m*2-1+nn][m*2-5] = T[m*2-1+nn][m*2-4+nn] = 
					-(T[m*2+nn][m*2-4] = -ff*(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*(n-m+1.)*(n-m+2.)*.0625/(m-1.));
		if (m > 1) T[m*2+nn][m*2-2+nn*2] = -(T[m*2-1+nn][m*2-3+nn*2] = ff*(n+m)*(n-m)*(n-m+1.)*.125/m);
	}
//...fz;
	T[nn*2][nn*2] = ff*(n-1.)*(n-1.);
	if (n >= 1) T[nn*2][2] = T[nn*2][1+nn] = -ff*(n-1.);
	for (m = 1; m <= N; m++) { 
		T[m*2+nn*2][m*2+nn*2] = -(T[m*2-1+nn*2][m*2-1+nn*2] = -ff*(n+m-1.)*(n-m-1.)*.5);
		if (m == 1)	T[m*2+nn*2][0] = -(T[m*2-1+nn*2][nn] = -ff*n*(2.*n-1.+(n-2.)*(n-1.)*.5)*.5);
		if (n >= 1+m) T[m*2+nn*2][m*2+2] = T[m*2+nn*2][m*2+1+nn] = T[m*2-1+nn*2][m*2+2+nn] = 
						-(T[m*2-1+nn*2][m*2+1] = ff*(m+1.)*(n-m-1.)*.5);
		if (m > 1) T[m*2+nn*2][m*2-3+nn] = T[m*2-1+nn*2][m*2-3] = T[m*2-1+nn*2][m*2-2+nn] = 
					-(T[m*2+nn*2][m*2-2] = ff*(n-m+1.)*(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*.25);
	}
	if (N == 1 && regul) {
		double G1 = 1./get_param(NUM_SHEAR);
		T[8][0] -= 1.25*G1;
		T[1][1] += .625*G1; T[4][1] -= .625*G1;
		T[2][2] += .625*G1; T[5][2] -= .625*G1; T[6][2] += 2.5*G1;
		T[7][3] += 1.25*G1;
		T[2][4] += .625*G1; T[5][4] += .625*G1; T[6][4] -= 2.5*G1;
		T[1][5] += .625*G1; T[4][5] += .625*G1;
		T[2][6] -= 1.25*G1;
		T[3][7] += 1.25*G1;
		T[0][8] -= 1.25*G1;
	}

/////////////////////////////////
//...раскомплексификация строчек;
	for (l = 3*nn-1; l >= 0; l--) {				
		T[0][l] *= 2.; T[nn][l] *= -2.;
		for (m = 1; m <= N; m++) {
			T[m*2][l] = 2.*(T[m*2+nn][l]+T[m*2][l]); 
			T[m*2-1][l] = -2.*(T[m*2-1+nn][l]+T[m*2-1][l]);
			T[m*2+nn][l] = 4.*T[m*2+nn][l]-T[m*2][l];
			T[m*2-1+nn][l] = 4.*T[m*2-1+nn][l]+T[m*2-1][l];
			T[m*2+nn*2][l] *=  2.; T[m*2-1+nn*2][l] *= -2.; 
			swap(T[m*2+nn][l], T[m*2-1+nn][l]);
		}
	}

//////////////////////////////////////////////
//...повторяем процедуру для матрицы перехода;
	if (n > 1) {
	  	if (! TH) set_matrix(TH, 3*mm, 3*nn);	
		else		clean_matrix(TH, 3*mm, 3*nn); N -= 2;
//...Re fw;
		TH[0][0] = -(TH[mm][nn] = n*(n-1.)*.125);
		TH[0][4] = TH[0][3+nn] = TH[mm][4+nn] = -(TH[mm][3] = -.5);
		TH[0][2+nn*2] = -(TH[mm][1+nn*2] = -(n-1.)*.25);
		for (m = 1; m <= N; m++) { 
			TH[m*2][m*2-1+nn] = TH[m*2-1][m*2-1] = TH[m*2-1][m*2+nn] = -(TH[m*2][m*2] = -(n-m-1.)*(n-m)*.0625);
			TH[m*2][m*2+4] = TH[m*2][m*2+3+nn] = TH[m*2-1][m*2+4+nn] = -(TH[m*2-1][m*2+3] = -(m+1.)*(m+2.)*.25);
			TH[m*2][m*2+2+nn*2] = -(TH[m*2-1][m*2+1+nn*2] = -(m+1.)*(n-m-1.)*.25);
		}
//...Im fw;
		for (m = 1; m <= N; m++) { 
			TH[m*2+mm][m*2] = TH[m*2+mm][m*2-1+nn] = TH[m*2-1+mm][m*2+nn] = -(TH[m*2-1+mm][m*2-1] = (n-m-1.)*(n-m)*.0625);
			if (m == 2) TH[m*2+mm][0] = -(TH[m*2-1+mm][nn] = -(n-3.)*(n-2.)*(n-1.)*n*.015625);
			if (m == 1)	{	
				double dd = -(n-2.)*(n-1.)*.0625;
				TH[m*2+mm][m*2] += dd; TH[m*2+mm][m*2-1+nn] += dd; TH[m*2-1+mm][m*2-1] += dd; TH[m*2-1+mm][m*2+nn] -= dd;
				TH[m*2+mm][nn*2] = -(n-2.)*(n-1.)*n*.125;
			}
			if (m > 2) TH[m*2+mm][m*2-5+nn] = TH[m*2-1+mm][m*2-4+nn] = TH[m*2-1+mm][m*2-5] = 
						-(TH[m*2+mm][m*2-4] = (n-m-1.)*(n-m)*(n-m+1.)*(n-m+2.)*.015625/(m*(m-1.)));
			if (m > 1) TH[m*2+mm][m*2-2+nn*2] = -(TH[m*2-1+mm][m*2-3+nn*2] = (n-m-1.)*(n-m)*(n-m+1.)*.0625/m);
		}
//...fz;
		TH[mm*2][nn*2] = n*(n-1.)*.5;
		TH[mm*2][2] = TH[mm*2][1+nn] = (n-1.)*.5;
		for (m = 1; m <= N; m++) { 
			TH[m*2+mm*2][m*2+nn*2] = -(TH[m*2-1+mm*2][m*2-1+nn*2] = -(n-m-1.)*(n-m)*.25);
			if (m == 1)	TH[m*2+mm*2][0] = -(TH[m*2-1+mm*2][nn] = (n-2.)*(n-1.)*n*.125 );
			TH[m*2+mm*2][m*2+2] = TH[m*2+mm*2][m*2+1+nn] = TH[m*2-1+mm*2][m*2+2+nn] = 
									  -(TH[m*2-1+mm*2][m*2+1] = -(m+1.)*(n-m-1.)*.25);
			if (m > 1) TH[m*2+mm*2][m*2-3+nn] = TH[m*2-1+mm*2][m*2-3] = TH[m*2-1+mm*2][m*2-2+nn] = 
						-(TH[m*2+mm*2][m*2-2] = -(n-m-1.)*(n-m)*(n-m+1.)*.0625/m);
		}

/////////////////////////////////
//...раскомплексификация строчек;
		for (l = 3*nn-1; l >= 0; l--) {				
			TH[0][l] *= 2.; TH[mm][l] *= -2.;
			for (m = 1; m <= N; m++) {
				TH[m*2][l] = 2.*(TH[m*2+mm][l]+TH[m*2][l]); 
				TH[m*2-1][l] = -2.*(TH[m*2-1+mm][l]+TH[m*2-1][l]);
				TH[m*2+mm][l] = 4.*TH[m*2+mm][l]-TH[m*2][l];
				TH[m*2-1+mm][l] = 4.*TH[m*2-1+mm][l]+TH[m*2-1][l];
				TH[m*2+mm*2][l] *=  2.; TH[m*2-1+mm*2][l] *= -2.; 
				swap(TH[m*2+mm][l], TH[m*2-1+mm][l]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//...формирование внешней матрицы в аналитическом методе для граничных функций;
void CHydro3D::extern_matrix(real_T **& T, int N)
{
	int l, m, n = N, nn = 2*n+1; 
	if (! T)  set_matrix(T, 3*nn, 3*nn);	
	else		clean_matrix(T, 3*nn, 3*nn);
//...Re fw;
	T[nn][nn] = -(T[0][0] = (n+1.)+(n-1.)*n*.125);
	if (n >= 2) T[0][4] = T[0][3+nn] = T[nn][4+nn] = -(T[nn][3] = .5);
	if (n >= 1) T[0][2+nn*2] = -(T[nn][1+nn*2] = (n+1.)*.25);
	for (m = 1; m <= N; m++) { 
		T[m*2][m*2-1+nn] = T[m*2-1][m*2-1] = T[m*2-1][m*2+nn] = -(T[m*2][m*2] = (n+1.)*.5-m*.25+(n-m-1.)*(n-m)*.0625);
		if (m+2 <= n) T[m*2][m*2+4] = T[m*2][m*2+3+nn] = T[m*2-1][m*2+4+nn] = -(T[m*2-1][m*2+3] = (m+1.)*(m+2.)*.25);
		if (m+1 <= n) T[m*2][m*2+2+nn*2] = -(T[m*2-1][m*2+1+nn*2] = (m+1.)*(n-m+1.)*.25);
	}
//...Im fw;
	for (m = 1; m <= N; m++) { 
		T[m*2+nn][m*2] = T[m*2+nn][m*2-1+nn] = T[m*2-1+nn][m*2+nn] = -(T[m*2-1+nn][m*2-1] = -((n+1.)*.5+m*(n+.5)*.25+(n-m-1.)*(n-m)*.0625));
		if (m == 2) T[m*2+nn][0] = -(T[m*2-1+nn][nn] = (2.*n-1.+(n-3.)*(n-2.)*.25)*(n-1.)*n*.0625);
		if (m == 1)	{	
			double dd = (n-.5)*.25+(n-2.)*(n-1.)*.0625;
			T[m*2+nn][m*2] += dd; T[m*2+nn][m*2-1+nn] += dd; T[m*2-1+nn][m*2-1] += dd; T[m*2-1+nn][m*2+nn] -= dd;
			T[m*2+nn][nn*2] = (2.*n+1.+(n-1.)*n*.5)*n*.25;
		}
		if (m > 2) T[m*2+nn][m*2-5+nn] = T[m*2-1+nn][m*2-5] = T[m*2-1+nn][m*2-4+nn] = 
					-(T[m*2+nn][m*2-4] = -(2.*n-1.+(n-m-1.)*(n-m)*.5/m)*(n-m+1.)*(n-m+2.)*.03125/(m-1.));
		if (m > 1) T[m*2+nn][m*2-2+nn*2] = -(T[m*2-1+nn][m*2-3+nn*2] = -(2.*n+1.+(n-m)*(n-m+1.)*.5/m)*(n-m+1.)*.125);
	}
//...fz;
	T[nn*2][nn*2] = (n+2.)*(n+2.)*.5;
	if (n >= 1) T[nn*2][2] = T[nn*2][1+nn] = (n+2.)*.5;
	for (m = 1; m <= N; m++) { 
		T[m*2+nn*2][m*2+nn*2] = -(T[m*2-1+nn*2][m*2-1+nn*2] = -((n+1.)+(n-m)*(n+m)*.25));
		if (m == 1)	T[m*2+nn*2][0] = -(T[m*2-1+nn*2][nn] = n*((n-1.)*(n+3.)+4.)*.125);
		if (n >= 1+m) T[m*2+nn*2][m*2+2] = T[m*2+nn*2][m*2+1+nn] = T[m*2-1+nn*2][m*2+2+nn] = 
						-(T[m*2-1+nn*2][m*2+1] = -(m+1.)*(n+m+2.)*.25);
		if (m > 1) T[m*2+nn*2][m*2-3+nn] = T[m*2-1+nn*2][m*2-3] = T[m*2-1+nn*2][m*2-2+nn] = 
					-(T[m*2+nn*2][m*2-2] = -((n-m)*(n+m+2.)+4.*m)*(n-m+1.)*.0625/m);
	}

/////////////////////////////////
//...раскомплексификация строчек;
	for (l = 3*nn-1; l >= 0; l--) {				
		T[0][l] *= 2.; T[nn][l] *= -2.;
		for (m = 1; m <= N; m++) {
			T[m*2][l] = 2.*(T[m*2+nn][l]+T[m*2][l]); 
			T[m*2-1][l] = -2.*(T[m*2-1+nn][l]+T[m*2-1][l]);
			T[m*2+nn][l] = 4.*T[m*2+nn][l]-T[m*2][l];
			T[m*2-1+nn][l] = 4.*T[m*2-1+nn][l]+T[m*2-1][l];
			T[m*2+nn*2][l] *=  2.; T[m*2-1+nn*2][l] *= -2.; 
			swap(T[m*2+nn][l], T[m*2-1+nn][l]);
		}
	}
}

/////////////////////////////////////////////////////
//...формироване матриц Эшелби для граничных функций;
void CHydro3D::set_eshelby_matrix(int N_mpl)
{
	real_T ** TX = NULL, * H = (real_T *)new_struct((2*N_mpl+1)*3*sizeof(real_T));
	TT = (real_T ***)new_struct((N_mpl+2)*sizeof(real_T **));
	TH = (real_T ***)new_struct((N_mpl+1)*sizeof(real_T **));

/////////////////////////////////////////////// 
//...формирование решения граничного уравнения;
	if (TT && TH && H) {
		solver.pivot_init((2*N_mpl+1)*3);

		int  k, j, l, num, nn;
		for (k = N_mpl; k >= 0; k--) { 
			intern_matrix(TT[k], k, TH[k-2]);
			extern_matrix(TX, k);
			if (solver.GaussJ(TX,  nn = (2*k+1)*3)) {
				for (num = 0; num < nn && TT[k]; num++) {
					for (           l = 0; l < nn; l++) 
					for (H[l] = 0., j = 0; j < nn; j++) H[l] += TX[l][j]*TT[k][j][num]; 
					for (           l = 0; l < nn; l++) TT[k][l][num] = H[l]/*TT[k][l][num]/(2.*k+3.)*/; 
				}
				for (num = 0; num < nn+12 && TH[k]; num++) {
					for (           l = 0; l < nn; l++) 
					for (H[l] = 0., j = 0; j < nn; j++) H[l] += TX[l][j]*TH[k][j][num]; 
					for (           l = 0; l < nn; l++) TH[k][l][num] = H[l]/*TH[k][l][num]/(2.*k+3.)*/; 
				}
			}
		}
	}
	delete_struct(H);
	delete_struct(TX);
}

///////////////////////////////////////////////////////////////////////////////
//...вычисление нормирующих множителей в представлении для уравнения Бринкмана;
void CHydro3D::set_normaliz_vector(int N_mpl, double rad, double kk)
{
	if (rad > 0. && kk > 0.) {
		real_T h1 = 1., h2 = -1./(kk*rad), h0 = -h2, hh = sqr(kk*rad);
		//real_T he = exp(kk*rad), hi = 1./he, h0 = 1./(kk*rad), h1 = (he+hi)*.5, h2, hh = sqr(kk*rad);
		int m;
		set_matrix(TN, 2, N_mpl+2);

///////////////////////////////////
//...prepare of the skin functions;
		for ( TN[0][m = 0] = h1/*, h2 = -(he-hi)*.5*h0*/; m <= N_mpl; TN[0][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}
		for ( TN[1][m = 0] = (h1 = 1.), h2 = h0; m <= N_mpl; TN[1][++m] = h1) {
		//for ( TN[1][m = 0] = (h1 = (he-hi)*.5), h2 = -(he+hi)*.5*h0; m <= N_mpl; TN[1][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}

////////////////////////////////                                                                                                                                              
//...пересчитываем коэффициенты;
		for ( m = 0, h1 = sqr(hh = -kk); m <= N_mpl; m++, hh *= -h1/(2.*m+1.)) { 
			h2 = TN[1][m];
			TN[1][m] =  0.5*TN[0][m]/hh;
			TN[0][m] = -0.5*h2/hh;
			//TN[1][m] =  TN[0][m]/hh;
			//TN[0][m] = -h2/hh;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//...прежний способ вычисления нормирующих множителей в представлении для уравнения Бринкмана;
void CHydro3D::set_normaliz_vector_old(int N_mpl, double rad, double kk)
{
	if (rad > 0. && kk > 0.) {
		real_T h1 = 1., h2 = -1./(kk*rad), h0 = -h2, hh = sqr(kk*rad);
		int m;
		set_matrix(TN, 2, N_mpl+2);

///////////////////////////////////
//...prepare of the skin functions;
		for ( TN[0][m = 0] = h1; m <= N_mpl; TN[0][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1; 
			  swap(h1, h2);
		}
		for ( TN[1][m = 0] = (h1 = 1.), h2 = h0; m <= N_mpl; TN[1][++m] = h1) {
			  h2 = hh*h2/(4.*m*m-1.)-h1 ; 
			  swap(h1, h2);
		}

////////////////////////////////
//...пересчитываем коэффициенты;
		for ( m = 0, h1 = sqr(hh = rad); m <= N_mpl; m++, hh *= h1/(2.*m-1.)) { 
			TN[1][m] = -TN[0][m]/TN[1][m]; h2 = (TN[0][m+1]+TN[1][m]*TN[1][m+1]);
			if (to_double(h2)) TN[0][m] = hh/h2; else TN[0][m] = 0.; 
			TN[1][m] *= TN[0][m]; 
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CHydro3D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6] = { X, Y, Z, 0., 0., 1.};

/////////////////////////////////////
//...operation with all input points;
	if (0 <= i && i < N && B[i].shape && B[i].mp) {
		int m = solver.id_norm;
//////////////////////////////////////////////////////////////////
//memset(solver.hh[i][0][0], 0, solver.dim[i]*sizeof(double));
//B[i].shape->FULL(solver.hh[i][0][0], 0, 2)[4] = B[i].shape->get_R();	
//B[i].shape->set_potential(solver.hh[i][0][0], 0);
////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
		case PRESSURE_VALUE: {
/////////////////////////////////////
//...calculation filtration pressure;
				B[i].shape->parametrization_hess(P, 1);
				jump0(P, i, 0); 

				F[0] = F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m], id_variant));
				F[0] += Z;
		}		break;
		case VELOCITY_VALUE: {
/////////////////////////////////////
//...calculation filtration velocity;
				B[i].shape->parametrization_hess(P, 1);
				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   id_variant));
				F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], id_variant));
////////////////////////////
//...вычисляем сглаживатель;
//double rr = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])), r0 = get_param(NUM_GEOMT), r1 = r0 + 0.03, 
//		smooth = 1.;
//if (rr < r1) smooth = 1.-exp((rr-r0)/(rr-r1));
////////////////////////////
				F[2] = /*smooth**/to_double(B[i].shape->potential(solver.hh[i][0][m+2], id_variant));
				B[i].shape->norm_common(F);
		}		break;
		case NORMAL_R_VALUE: {
/////////////////////////////////////////////
//...normal component of filtration velocity;
			   P[5] = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
				P[3] = P[0]/P[5]; P[4] = P[1]/P[5]; P[5] = P[2]/P[5];

				B[i].shape->parametrization_hess(P, 1);
				jump1_x(P, i, 0); 
				jump1_y(P, i, 1); 
				jump1_z(P, i, 2); 

				F[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   id_variant));
				F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], id_variant));
				F[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], id_variant));
				B[i].shape->norm_common(F); 
				F[0] = F[0]*P[3]+F[1]*P[4]+F[2]*P[5];
		}		break;
		case NORMAL_X_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   id_variant));
				F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], id_variant));
				F[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], id_variant));
				B[i].shape->norm_common(F);
		}   break;
		case NORMAL_Y_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   id_variant));
				F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], id_variant));
				F[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], id_variant));
				B[i].shape->norm_common(F);
		}   break;
		case NORMAL_Z_VALUE: {
///////////////////////////////////////////////
//...normal derivatives of filtration velocity;
				P[5] = 1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_hess(P, 1);
				if (get_param(NUM_SHEAR+1)) hessian_deriv_N(P, i); else hessian_deriv_N(1, P, i);

				jump2_x(P, i, 0); 
				jump2_y(P, i, 1); 
				jump2_z(P, i, 2); 

				F[0] = to_double(B[i].shape->potential(solver.hh[i][0][m],   id_variant));
				F[1] = to_double(B[i].shape->potential(solver.hh[i][0][m+1], id_variant));
				F[2] = to_double(B[i].shape->potential(solver.hh[i][0][m+2], id_variant));
				B[i].shape->norm_common(F);
		}   break;
		case POTENTIAL_VALUE: {
//////////////////////////////////////////
//...calculation potentials fX, fY and fZ;
				B[i].shape->parametrization_hess(P, 1);

				F[0]  = to_double(B[i].shape->potential(0, B[i].shape->p_cpy, id_variant));
				F[1]  = to_double(B[i].shape->potential(1, B[i].shape->p_cpy-B[i].shape->get_NN(), id_variant));
				F[2]  = to_double(B[i].shape->potential(2, B[i].shape->p_cpy-B[i].shape->get_NN()*2, id_variant));
				B[i].shape->norm_common(F);
		}		break;
		case ANALYT_VALUE: {
/////////////////////////////////////
//...flow in the cylindrical channel;
				F[0] = F[1] = 0.;
				F[2] = CylinderVelocity(sqr(P[0])+sqr(P[1]), sqr(get_param(NUM_GEOMT)), get_param(NUM_SHEAR+1));
				B[i].shape->norm_common(F);
		}		break;
		case TESTI_VALUE: {
/////////////////////////////////
//...test calculation potentials;
				B[i].shape->parametrization_hess(P, 1);

				F[0]  = to_double(B[i].shape->potential(0, B[i].shape->FULL(B[i].shape->p_cpy, 0, iparam),	id_variant));
				F[1]  = to_double(B[i].shape->potential(0, B[i].shape->FULL(B[i].shape->p_cpy, 0, iparam+1), id_variant));
				F[2]  = to_double(B[i].shape->potential(0, B[i].shape->FULL(B[i].shape->p_cpy, 0, iparam+2), id_variant));
		}		break;
		default: F[0] = i; F[1] = F[2] = 0.;
		}
	}
}

//////////////////////////////////////////////////
//...проницаемость Бринкмана в трещиноватой среде;
double CHydro3D::TakeLayer(double ff)
{
	double ll = ff*.5, alpha = get_param(NUM_SHEAR+1);
	if (alpha == 0.) return(2.*sqr(ll)*ll/3.);
	else {
		double kappa = sqrt(fabs(alpha))*ll, EE = exp(-2.*kappa);
		return 2.*ll*(1.-(1.-EE)/((1.+EE)*kappa))/alpha;
	}
}

////////////////////////////////////////////////////
//...проницаемость Бринкмана в цилиндрических порах;
double CHydro3D::TakeCylinder(double ff, double eps)
{ //...rr - квадрат радиуса;
	double rr = fabs(ff/M_PI), alpha = get_param(NUM_SHEAR+1);
	if (alpha == 0.) return(M_PI*sqr(rr)*.125);
	else {
		double kappa = fabs(alpha)*rr*.25, f1 = .5, f2 = 1., sum1 = f1, sum2 = f2;
		int k = 0, k_max = 100;
		while ((f1 > eps || f2 > eps) && k++ < k_max ) {
			sum1 += (f1 *= kappa/(k*(k+2.)));
			sum2 += (f2 *= kappa/(k*(k+0.)));
		}
		return M_PI*sqr(rr)*sum1/(4.*sum2);
	}
}

//////////////////////////////////////////////
//...течение Бринкмана в цилиндрических порах;
double CHydro3D::CylinderVelocity(double rr, double RR, double eps)
{ //...rr, RR - квадраты радиусов;
	double alpha = get_param(NUM_SHEAR+1);
	if (alpha == 0.) return((RR-rr)*.25);
	else {
		double kappa1 = fabs(alpha)*rr*.25, kappa2 = fabs(alpha)*RR*.25, f1 = 1., f2 = 1., sum1 = f1, sum2 = f2;
		int k = 0, k_max = 100;
		while ((f1 > eps || f2 > eps) && k++ < k_max ) {
			sum1 += (f1 *= kappa1/(k*k));
			sum2 += (f2 *= kappa2/(k*k));
		}
		return((1.-sum1/sum2)/fabs(alpha));
	}
}
#undef  Message
