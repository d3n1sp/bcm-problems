#include "stdafx.h"

#include "shapes.h"
#include "cmindl3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CMindl3D::NUM_ADHES = 4;
int CMindl3D::NUM_SHEAR = 6;
int CMindl3D::NUM_SHIFT = 4;
int CMindl3D::MAX_PHASE = 2;
int CMindl3D::NUM_HESS  = 8;

//////////////////////////////////
//...initialization of the blocks;
int  CMindl3D::block_shape_init(Block<double> & B, Num_State id_free)
{
	int k;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<double>(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(MP3D_ZOOM_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
 			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, sqr(B.mp[8])/(get_param(NUM_MPLS+1)*fabs(B.mp[7])));

////////////////////////////////////////////////////////////////
//...подключаем дополнительные экспоненциальные системы функций;
			B.shape->add_shape(CreateShape<double>(SK3D_EXPP_SHAPE),	  NULL_STATE);
			B.shape->add_shape(CreateShape<double>(SK3D_EXPP_SHAPE, 1), NULL_STATE);

			B.shape->add_shape(CreateShape<double>(SK3D_EXPP_SHAPE),	  NULL_STATE);
			B.shape->add_shape(CreateShape<double>(SK3D_EXPP_SHAPE, 1), NULL_STATE);

			B.shape->init1(2, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3));

			B.shape->init1(3, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3));

			B.shape->init1(4, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2));

			B.shape->init1(5, UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2));
		}
		else {
			B.shape->add_shape(CreateShape<double>(MP3D_POLY_SHAPE));

			B.shape->init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

////////////////////////////////////////////////////////////
//...подключаем регулярные экспоненциальные системы функций;
			B.shape->add_shape(CreateShape<double>(SK3D_ZOOM_SHAPE), NULL_STATE);
			B.shape->add_shape(CreateShape<double>(SK3D_ZOOM_SHAPE), NULL_STATE);

			B.shape->init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3));

			B.shape->init1(2, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2));
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
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3), fabs(B.mp[8])); 
				B.shape->set_shape(3, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3), fabs(B.mp[8])); 
				B.shape->set_shape(4, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2), fabs(B.mp[8])); 
				B.shape->set_shape(5, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2), fabs(B.mp[8])); 
			}
 			else {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3), fabs(B.mp[7]));
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2), fabs(B.mp[7]));
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

//////////////////////////////////////////////
//...realization of common displacements (Rx);
void CMindl3D::jump1_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr;
////////////////////////////////////
//...лапласовская часть перемещений;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_x     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, (alpha+1.)*G0, P[0]*alpha*G0);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0.,	P[1]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0.,	P[2]*alpha*G0);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (; l < nm*2; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0, 1.);

		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., 1.);

		B[i].shape->cpy_xz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., 1.);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_xz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -1.);
	}
}

//////////////////////////////////////////////
//...realization of common displacements (Ry);
void CMindl3D::jump1_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr;
////////////////////////////////////
//...лапласовская часть перемещений;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_y     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, (alpha+1.)*G0, P[1]*alpha*G0);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0.,	P[0]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0.,	P[2]*alpha*G0);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (; l < nm*2; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., 1.);

		B[i].shape->cpy_yy(l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0, 1.);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., 1.);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_yy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -1.);
	}
}

//////////////////////////////////////////////
//...realization of common displacements (Rz);
void CMindl3D::jump1_common_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr;
////////////////////////////////////
//...лапласовская часть перемещений;
	for (l = 0; l < nm; l++) {
		B[i].shape->cpy_z     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, (alpha+1.)*G0, P[2]*alpha*G0);

		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0.,	P[0]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0.,	P[2]*alpha*G0);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (; l < nm*2; l++) {
		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., 1.);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., 1.);

		B[i].shape->cpy_zz(l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0, 1.);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_yz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -1.);

		B[i].shape->cpy_zz(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -1.);
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (DnRx);
void CMindl3D::jump2_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr; m += solver.id_norm;
///////////////////////////////////////////////
//...лапласовская часть нормальных производных;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xx    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*(2.*alpha+1.)*G0, P[0]*alpha*G0);
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.)*G0);
		B[i].shape->adm_z     (l, ptr,  P[5]*(alpha+1.)*G0);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*alpha*G0, P[1]*alpha*G0);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5]*alpha*G0, P[2]*alpha*G0);
	}
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l), ptr, 0., 1.);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
		B[i].shape->adm_z(l, ptr, -P[5]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., 1.);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -1.);
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (DnRy);
void CMindl3D::jump2_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr; m += solver.id_norm;
///////////////////////////////////////////////
//...лапласовская часть нормальных производных;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*(2.*alpha+1.)*G0, P[1]*alpha*G0);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.)*G0);
		B[i].shape->adm_z     (l, ptr,  P[5]*(alpha+1.)*G0);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*alpha*G0, P[0]*alpha*G0);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5]*alpha*G0, P[2]*alpha*G0);
	}
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., 1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
		B[i].shape->adm_z(l, ptr, -P[5]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., 1.);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -1.);
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (DnRz);
void CMindl3D::jump2_common_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = 1./get_param(NUM_SHEAR+shift), alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 =  get_param(NUM_SHEAR+3+shift), * ptr; m += solver.id_norm;
///////////////////////////////////////////////
//...лапласовская часть нормальных производных;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_zz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_z     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5]*(2.*alpha+1.)*G0, P[2]*alpha*G0);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.)*G0);
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.)*G0);

		B[i].shape->cpy_z     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*alpha*G0, P[0]*alpha*G0);

		B[i].shape->cpy_z     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*alpha*G0, P[1]*alpha*G0);
	}
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 0); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., 1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., 1.);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
		B[i].shape->adm_z(l, ptr, -P[5]*C0);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 0); //...deriv_N_xz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -1.);
	}
}

///////////////////////////////////////////
//...realization of cohesion moment M(n)_x;
void CMindl3D::jump3_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), 
			 C0 = get_param(NUM_SHEAR+3+shift), fnx = (G0+alpha*sqr(P[3])), * ptr;
/////////////////////////////////////////////////////////
//...лапласовская часть когезионного момента (обнуление);
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), NULL, 0., 0.);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*fnx, fnx);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_xz	 (l, ptr, alpha*P[3]*P[5]);

		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[4], fnx);
		B[i].shape->adm_yy	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_yz	 (l, ptr, alpha*P[3]*P[5]);

		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[5], fnx);
		B[i].shape->adm_yz	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_zz	 (l, ptr, alpha*P[3]*P[5]);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -fnx);
		B[i].shape->adm_xy	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_xz	 (l, ptr, -alpha*P[3]*P[5]);

		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -fnx);
		B[i].shape->adm_yy	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_yz	 (l, ptr, -alpha*P[3]*P[5]);

		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -fnx);
		B[i].shape->adm_yz	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_zz	 (l, ptr, -alpha*P[3]*P[5]);
	}
}

///////////////////////////////////////////
//...realization of cohesion moment M(n)_y;
void CMindl3D::jump3_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), 
			 C0 = get_param(NUM_SHEAR+3+shift), fny = (G0+alpha*sqr(P[4])), * ptr;
/////////////////////////////////////////////////////////
//...лапласовская часть когезионного момента (обнуление);
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), NULL, 0., 0.);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[4], fny);
		B[i].shape->adm_xx	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_xz	 (l, ptr, alpha*P[4]*P[5]);

		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*fny, fny);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_yz	 (l, ptr, alpha*P[4]*P[5]);

		B[i].shape->cpy_yz	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[4]*P[5], fny);
		B[i].shape->adm_xz	 (l, ptr, alpha*P[3]*P[4]);
		B[i].shape->adm_zz	 (l, ptr, alpha*P[4]*P[5]);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -fny);
		B[i].shape->adm_xx	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_xz	 (l, ptr, -alpha*P[4]*P[5]);

		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -fny);
		B[i].shape->adm_xy	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_yz	 (l, ptr, -alpha*P[4]*P[5]);

		B[i].shape->cpy_yz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -fny);
		B[i].shape->adm_xz	 (l, ptr, -alpha*P[3]*P[4]);
		B[i].shape->adm_zz	 (l, ptr, -alpha*P[4]*P[5]);
	}
}

///////////////////////////////////////////
//...realization of cohesion moment M(n)_z;
void CMindl3D::jump3_common_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), 
			 C0 = get_param(NUM_SHEAR+3+shift), fnz = (G0+alpha*sqr(P[5])), * ptr;
/////////////////////////////////////////////////////////
//...лапласовская часть когезионного момента (обнуление);
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 2), NULL, 0., 0.);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[5], fnz);
		B[i].shape->adm_xx	 (l, ptr, alpha*P[3]*P[5]);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[4]*P[5]);

		B[i].shape->cpy_yz	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[4]*P[5], fnz);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[3]*P[5]);
		B[i].shape->adm_yy	 (l, ptr, alpha*P[4]*P[5]);

		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*fnz, fnz);
		B[i].shape->adm_xz	 (l, ptr, alpha*P[3]*P[5]);
		B[i].shape->adm_yz	 (l, ptr, alpha*P[4]*P[5]);
	}
	for (; l < nm*3; l++) {
		B[i].shape->cpy_xz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., -fnz);
		B[i].shape->adm_xx	 (l, ptr, -alpha*P[3]*P[5]);
		B[i].shape->adm_xy	 (l, ptr, -alpha*P[4]*P[5]);

		B[i].shape->cpy_yz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., -fnz);
		B[i].shape->adm_xy	 (l, ptr, -alpha*P[3]*P[5]);
		B[i].shape->adm_yy	 (l, ptr, -alpha*P[4]*P[5]);

		B[i].shape->cpy_zz	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), B[i].shape->deriv, 0., -fnz);
		B[i].shape->adm_xz	 (l, ptr, -alpha*P[3]*P[5]);
		B[i].shape->adm_yz	 (l, ptr, -alpha*P[4]*P[5]);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical surface tensions for common displacements Px(R);
void CMindl3D::jump4_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), C0 = get_param(NUM_SHEAR+3+shift)*G0*.5,
			 C1 = get_param(NUM_SHEAR+2+shift)*G0*get_param(NUM_SHEAR+1+shift)/(1.-get_param(NUM_SHEAR+1+shift)*2.), * ptr; m += solver.id_norm;
/////////////////////////////////////////////////
//...лапласовская часть поверхностных напряжений;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xx    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3], P[0]*alpha);
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.));
		B[i].shape->adm_z     (l, ptr,  P[5]*(alpha+1.));

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*(alpha+1.), P[1]*alpha);
		B[i].shape->adm_y     (l, ptr,  P[3]*(-2.*alpha-1.));

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5]*(alpha+1.), P[2]*alpha);
		B[i].shape->adm_z     (l, ptr,  P[3]*(-2.*alpha-1.));
	}
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C0*2.);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
		B[i].shape->adm_z(l, ptr, -P[5]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[5]*C0);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[3]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_xz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[3]*C1);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical surface tensions for common displacements Py(R);
void CMindl3D::jump4_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), C0 = get_param(NUM_SHEAR+3+shift)*G0*.5,
			 C1 = get_param(NUM_SHEAR+2+shift)*G0*get_param(NUM_SHEAR+1+shift)/(1.-get_param(NUM_SHEAR+1+shift)*2.), * ptr; m += solver.id_norm;
/////////////////////////////////////////////////
//...лапласовская часть поверхностных напряжений;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yy    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*(alpha+1.), P[0]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[4]*(-2.*alpha-1.));

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4], P[1]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.));
		B[i].shape->adm_z     (l, ptr,  P[5]*(alpha+1.));

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5]*(alpha+1.), P[2]*alpha);
		B[i].shape->adm_z     (l, ptr,  P[4]*(-2.*alpha-1.));
	}
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C0*2.);
		B[i].shape->adm_x(l, ptr, -P[3]*C0);
		B[i].shape->adm_z(l, ptr, -P[5]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[5]*C0);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[4]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 2); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[4]*C1);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical surface tensions for common displacements Pz(R);
void CMindl3D::jump4_common_z(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), C0 = get_param(NUM_SHEAR+3+shift)*G0*.5,
			 C1 = get_param(NUM_SHEAR+2+shift)*G0*get_param(NUM_SHEAR+1+shift)/(1.-get_param(NUM_SHEAR+1+shift)*2.), * ptr; m += solver.id_norm;
/////////////////////////////////////////////////
//...лапласовская часть поверхностных напряжений;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xz    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yz    (l, B[i].shape->deriv, P[4]);
		B[i].shape->adm_zz    (l, B[i].shape->deriv, P[5]);

		B[i].shape->cpy_z     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*(alpha+1.), P[0]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[5]*(-2.*alpha-1.));

		B[i].shape->cpy_z     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*(alpha+1.), P[1]*alpha);
		B[i].shape->adm_y     (l, ptr,  P[5]*(-2.*alpha-1.));

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[5], P[2]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.));
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.));
	}
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 0); //...deriv_N_xz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[4]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[5]*C0*2.);
		B[i].shape->adm_x(l, ptr, -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
	}
	for (; l < nm*3; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 0); //...deriv_N_xz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[5]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 1); //...deriv_N_yz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[5]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l, 2); //...deriv_N_zz;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 2), ptr, 0., -G0);
		B[i].shape->adm_z(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 2), -P[5]*C1);
	}
}

////////////////////////////////////////////////////////
//...realization partial normal derivatives of hessian;
void CMindl3D::hessian_deriv_N(int k, double * P, int i)
{
	int num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->set_norm_cs(P);

	B[i].shape->cpy_xx(k);	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xx(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k));

	B[i].shape->cpy_xy(k); B[i].shape->deriv_N(k);
	B[i].shape->cpy_xy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k, 1));

	B[i].shape->cpy_yy(k);	B[i].shape->deriv_N(k);
	B[i].shape->cpy_yy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k, 2));

	B[i].shape->cpy_xz(k); B[i].shape->deriv_N(k);
	B[i].shape->cpy_xz(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k));

	B[i].shape->cpy_yz(k);	B[i].shape->deriv_N(k);
	B[i].shape->cpy_yz(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k, 1));

	B[i].shape->cpy_zz(k);	B[i].shape->deriv_N(k);
	B[i].shape->cpy_zz(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k, 2));
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CMindl3D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B)
{
	if (size_of_param() > NUM_ADHES+2+NUM_SHIFT) {
		param[NUM_ADHES] = A;
		param[NUM_ADHES+1] = B;
	}
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT) {
      param[NUM_SHEAR] = G1;
      param[NUM_SHEAR+NUM_SHIFT] = G2;
      param[NUM_SHEAR+1] = nju1;
      param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
		param[NUM_SHEAR+2] = 1./sqrt(l11*l11+l12*l12);
		param[NUM_SHEAR+2+NUM_SHIFT] = 1./sqrt(l21*l21+l22*l22);
		param[NUM_SHEAR+3] = 1./fabs(l11);
		param[NUM_SHEAR+3+NUM_SHIFT] = 1./fabs(l21);
	}
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (прямой алгоритм);
double CMindl3D::TakeEshelby_volm_two(double ff)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), 
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, 
			 K_I = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), 
			 K_M = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M), AA = get_param(NUM_ADHES);
	if (get_param(NUM_SHEAR+3) == 0. || 
		 get_param(NUM_SHEAR+3+NUM_SHIFT) == 0.) { //...классическая трехфазная модель;
		return(K_M+ff*(K_I-K_M)/(1.+(1.-ff)*(K_I-K_M)/(K_M+4./3.*mu_M)));
	}
	double kk_I = get_param(NUM_SHEAR+2+NUM_SHIFT), tt_I = exp(-2.*kk_I), RR2 = 1./pow(ff, 1./3.), 
			 kk_M = get_param(NUM_SHEAR+2), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), JHP1 =  (kk_M*RR2-1.)*ff*tt_D, 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), JHM1 = -(kk_M*RR2+1.)*ff/tt_D,
			 matr[7][8] = {
					{ 1., -HH1, -1., -1., JJP1, JJM1, 0., 0.},
					{ 0., -HH2,  0.,  3., JJP2, JJM2, 0., 0.},
					{ AA, -ku_I*HH1-AA*(HH2+HH1), 0.,0., ku_M*JJP1, ku_M*JJM1, 0., 0.},
//...почему индекс 2, а не 4 при mu_I и mu_M ???
					//{ 3.*K_I, (ku_I+2.*mu_I)*HH1, -3.*K_M, 4.*mu_M, -(ku_M+2.*mu_M)*JJP1, -(ku_M+2.*mu_M)*JJM1, 0., 0.},
					{ 3.*K_I, (ku_I+4.*mu_I)*HH1, -3.*K_M, 4.*mu_M, -(ku_M+4.*mu_M)*JJP1, -(ku_M+4.*mu_M)*JJM1, 0., 0.},
					{ 0., 0., 0., 0., JHP1, JHM1,  0., 0.},
					{ 0., 0., 3.*K_M, -4.*mu_M*ff, 0., 0., -3., 0.}, //...эффективный модуль с коэффициентом три;
					{ 0., 0., 1., ff, 0., 0., 0., 1.},					 //...перемещения дают единицу в правую часть;
	};

//////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, D0, A1, B1, C1, D1, KH;
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

/////////////////////////////////////////////////////////////////////
//...распределение поля перемещений в среднем сечении (проверить!!!);
	int id_visual = 0;
	if (id_visual) {
		FILE * TST = fopen("CCohes3D_sphere_eshelby.dat", "w");
		double A = .7, X, Y = 0., Z, RR, func; int M = 400;
		for (i = 0; i <= M; i++)
		for (l = 0; l <= M; l++) {
			X = (-1.+i*2./M)*A; Z = (-1.+l*2./M)*A; RR = sqrt(X*X+Y*Y+Z*Z); func = RR;
			//if (RR < 1.0) func *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR)); else
			//if (RR < RR2) func *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
			fprintf(TST, "%g   %g   %g\n", Z, X, func-Z);
		}
		fclose(TST);
		if (NULL_STATE) { //...проверка сшивки в одной точке;
			//X = 0.; Z = RR1; RR = sqrt(X*X+Y*Y+Z*Z); double func1 = RR, func2 = RR;
			//func1 *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR));
			//func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 

			//X = 0.; Z = 1.; RR = sqrt(X*X+Y*Y+Z*Z); func1 = Z; func2 = Z;
			//func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
		}
	}
	return(matr[6][7]);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для трехфазной модели (модуль сдвига);
double take_system2(double matrix[14][15], double mu_MH, double nu_H)
{
	int dim_N = 14, ii[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[10][12] *= mu_MH;
	matrix[10][13] *= mu_MH;
	matrix[10][14] *= mu_MH;
	matrix[11][12] *= mu_MH;
	matrix[11][14] *= mu_MH;
	matrix[11][13] *= mu_MH*2.*(1.-2.*nu_H)/(5.-4.*nu_H);
	matrix[12][13] *= 2.*(5.-nu_H)/(5.-4.*nu_H);
	matrix[13][13] *= 2.*(1.+nu_H)/(5.-4.*nu_H);
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
	return(matrix[13][14]);
}
//////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (итерационный алгоритм);
double CMindl3D::TakeEshelby_shear_two(double ff, double eps, int max_iter)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR2 = 1./pow(ff, 1./3.), fR2 = RR2*RR2, fR5 = ff/fR2, 
			 mu_MI = mu_M/mu_I, AA = get_param(NUM_ADHES)/mu_I, BB = get_param(NUM_ADHES+1)/mu_I,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M), 
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, sgn0, sgn1, nu_H, mu_MH, mu_MH0, mu_MH1, matrix[14][15], f1 = 0.5/*1.0*/, f2 = 2.0/*1.0*/,
			 matr[14][15] = {
					{ mu_MI, -2.*mu_MI*QI1, 0., 0., -1., -1., -3., 2.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -mu_MI,			0., 0., -1., -QM3, 2., 1., 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -6.*mu_MI*QI1, 0., 0., -1., 2., 12., 6.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -3.*mu_MI,		0., 0., -1., 2.*QM3, -8., 3., 0., 0., 0., 0., 0., 0., 0.},
					{ AA*f1, -6.*QI1*AA*f1,	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, //...надо уполовинить два уравнения (22.02.14)!
					{ BB*f1, -3.*BB*f1,		0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 1.,  QI1,					0., 0., -1.,  QM4, 12., -QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -QI2,					0., 0., -1., -QM5, -8.,  QM2, 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 1., ff, 3.*fR5, -2.*QM1*fR2, 0., 0., 0., 0., -3., -1., 1.}, //...нормировку по радиусу в последних коэффициентах можно убрать!!!
					{ 0., 0., 0., 0., 1.,  QM3*ff, -2.*fR5, -fR2,  0., 0., 0., 0.,  2., -1., 1.},
					{ 0., 0., 0., 0., 1., -QM4*ff, -12.*fR5, QM1*fR2, 0., 0., 0., 0., 12.,  1., 1.},
					{ 0., 0., 0., 0., 1.,  QM5*ff, 8.*fR5,  -QM2*fR2, 0., 0., 0., 0., -8., -1. , 1.},
	};
	int k_iter = 0, k = 0, l;
////////////////////////////////////////////
//...дописываем когезионную часть в матрице;
	double kk_I = get_param(NUM_SHEAR+2+NUM_SHIFT), tt_I = exp(-2.*kk_I),
			 kk_M = get_param(NUM_SHEAR+2), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), 
			 JHP1 =  (kk_M*RR2-1.)*ff*tt_D, JHP2 = ((sqr(kk_M*RR2)+3.)-3.*kk_M*RR2)*fR5*tt_D, 
			 JHM1 = -(kk_M*RR2+1.)*ff/tt_D, JHM2 = ((sqr(kk_M*RR2)+3.)+3.*kk_M*RR2)*fR5/tt_D,
			 kp_I = get_param(NUM_SHEAR+3+NUM_SHIFT), tp_I = exp(-2.*kp_I), 
			 kp_M = get_param(NUM_SHEAR+3), tp_D = exp(kp_M*(RR2-1.)),
			 HP1  = ((1.+tp_I)*kp_I-(1.-tp_I))*.5, 
			 HP2  = ((1.-tp_I)*(sqr(kp_I)+3.)-3.*(1.+tp_I)*kp_I)*.5,
			 JPP1 =  (kp_M-1.), JPP2 = ((sqr(kp_M)+3.)-3.*kp_M), 
			 JPM1 = -(kp_M+1.), JPM2 = ((sqr(kp_M)+3.)+3.*kp_M), 
			 JDP1 =  (kp_M*RR2-1.)*ff*tp_D, JDP2 = ((sqr(kp_M*RR2)+3.)-3.*kp_M*RR2)*fR5*tt_D, 
			 JDM1 = -(kp_M*RR2+1.)*ff/tp_D, JDM2 = ((sqr(kp_M*RR2)+3.)+3.*kp_M*RR2)*fR5/tt_D;
////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = -f2*mu_M*HP2*3.; //...надо удвоить четыре уравнения (22.02.14)!
	matr[0][3] = -f2*mu_M*(HH1*sqr(kk_I)-3.*HH2);
	matr[0][8] = f2*3.*mu_M*JPP2; 
	matr[0][9] = f2*3.*mu_M*JPM2; 
	matr[0][10] = f2*mu_M*(JJP1*sqr(kk_M)-3.*JJP2); 
	matr[0][11] = f2*mu_M*(JJM1*sqr(kk_M)-3.*JJM2);

	matr[1][2] = -f2*mu_M*(HP1*sqr(kp_I)-2.*HP2); 
	matr[1][3] = -f2*mu_M*HH2*2.;
	matr[1][8] = f2*mu_M*(JPP1*sqr(kp_M)-2.*JPP2); 
	matr[1][9] = f2*mu_M*(JPM1*sqr(kp_M)-2.*JPM2); 
	matr[1][10] = f2*mu_M*JJP2*2.; 
	matr[1][11] = f2*mu_M*JJM2*2.;

	matr[2][2] = -f2*mu_M*(HP1*sqr(kp_I)-4.*HP2)*3.; 
	matr[2][3] = -f2*mu_M*((HH2-2.*HH1)*sqr(kk_I)+HH2*12.);
	matr[2][8] = f2*mu_M*(JPP1*sqr(kp_M)-4.*JPP2)*3.; 
	matr[2][9] = f2*mu_M*(JPM1*sqr(kp_M)-4.*JPM2)*3.; 
	matr[2][10] = f2*mu_M*((JJP2-2.*JJP1)*sqr(kk_M)+JJP2*12.); 
	matr[2][11] = f2*mu_M*((JJM2-2.*JJM1)*sqr(kk_M)+JJM2*12.);

	matr[3][2] = -f2*mu_M*((HP2-HP1)*sqr(kp_I)+HP2*8.); 
	matr[3][3] = -f2*mu_M*(HH1*sqr(kk_I)-4.*HH2)*2.;
	matr[3][8] = f2*mu_M*((JPP2-JPP1)*sqr(kp_M)+JPP2*8.); 
	matr[3][9] = f2*mu_M*((JPM2-JPM1)*sqr(kp_M)+JPM2*8.); 
	matr[3][10] = f2*mu_M*(JJP1*sqr(kk_M)-4.*JJP2)*2.; 
	matr[3][11] = f2*mu_M*(JJM1*sqr(kk_M)-4.*JJM2)*2.;

	matr[4][2] = -ku_I*HP2*3.-AA*mu_I*(HP1*sqr(kp_I)-4.*HP2)*3.; 
	matr[4][3] = -ku_I*(HH1*sqr(kk_I)-3.*HH2)-AA*mu_I*((HH2-2.*HH1)*sqr(kk_I)+12.*HH2);
	matr[4][8] = ku_M*JPP2*3.; 
	matr[4][9] = ku_M*JPM2*3.;
	matr[4][10] = ku_M*(JJP1*sqr(kk_M)-3.*JJP2); 
	matr[4][11] = ku_M*(JJM1*sqr(kk_M)-3.*JJM2);

	matr[5][2] = -mu_I*(HP1*sqr(kp_I)-2.*HP2)-BB*mu_I*((HP2-HP1)*sqr(kp_I)+8.*HP2); 
	matr[5][3] = -mu_I*HH2*2.-BB*mu_I*(HH1*sqr(kk_I)-4.*HH2)*2.;
	matr[5][8] = mu_M*(JPP1*sqr(kp_M)-2.*JPP2); 
	matr[5][9] = mu_M*(JPM1*sqr(kp_M)-2.*JPM2);
	matr[5][10] = mu_M*JJP2*2.; 
	matr[5][11] = mu_M*JJM2*2.;

	matr[6][2] = -(mu_I*HP1*sqr(kp_I)-(6.*mu_I+ku_I)*HP2)*3.; 
	matr[6][3] = ((4.*mu_I+ku_I)*HH1*sqr(kk_I)-3.*(6.*mu_I+ku_I)*HH2);
	matr[6][8] = (mu_M*JPP1*sqr(kp_M)-(6.*mu_M+ku_M)*JPP2)*3.; 
	matr[6][9] = (mu_M*JPM1*sqr(kp_M)-(6.*mu_M+ku_M)*JPM2)*3.; 
	matr[6][10] = -((4.*mu_M+ku_M)*JJP1*sqr(kk_M)-3.*(6.*mu_M+ku_M)*JJP2); 
	matr[6][11] = -((4.*mu_M+ku_M)*JJM1*sqr(kk_M)-3.*(6.*mu_M+ku_M)*JJM2);

	matr[7][2] = (2.*mu_I*HP1*sqr(kp_I)-(10.*mu_I-3.*ku_I)*HP2)*2.; 
	matr[7][3] = ((ku_I-2.*mu_I)*HH1*sqr(kk_I)+(10.*mu_I-3.*ku_I)*HH2)*2.;
	matr[7][8] = -(2.*mu_M*JPP1*sqr(kp_M)-(10.*mu_M-3.*ku_M)*JPP2)*2.; 
	matr[7][9] = -(2.*mu_M*JPM1*sqr(kp_M)-(10.*mu_M-3.*ku_M)*JPM2)*2.; 
	matr[7][10] = -((ku_M-2.*mu_M)*JJP1*sqr(kk_M)+(10.*mu_M-3.*ku_M)*JJP2)*2.; 
	matr[7][11] = -((ku_M-2.*mu_M)*JJM1*sqr(kk_M)+(10.*mu_M-3.*ku_M)*JJM2)*2.;

	matr[8][8] = 3.*JDP2; 
	matr[8][9] = 3.*JDM2; 
	matr[8][10] = (JHP1*sqr(kk_M)-3.*JHP2); 
	matr[8][11] = (JHM1*sqr(kk_M)-3.*JHM2);

	matr[9][8] = (JDP1*sqr(kp_M)-2.*JDP2); 
	matr[9][9] = (JDM1*sqr(kp_M)-2.*JDM2); 
	matr[9][10] = 2.*JHP2; 
	matr[9][11] = 2.*JHM2;

//////////////////////////////////////////////////////////////////////////////////////////
//...упрощенная формулировка когезионных сил в граничном условии на поверхности включений;	
	//matr[6][2] = -(HP1*sqr(kp_I)-4.*HP2)*3.; 
	//matr[6][3] =  (HH1*sqr(kk_I)-6.*HH2)*2.;
	//matr[6][8] = (JPP1*sqr(kp_M)-4.*JPP2)*3.; 
	//matr[6][9] = (JPM1*sqr(kp_M)-4.*JPM2)*3.; 
	//matr[6][10] = -(JJP1*sqr(kk_M)-6.*JJP2)*2.; 
	//matr[6][11] = -(JJM1*sqr(kk_M)-6.*JJM2)*2.;

	//matr[7][2] = -(HP1*sqr(kp_I)-8.*HP2); 
	//matr[7][3] =  (HH1*sqr(kk_I)-4.*HH2)*2.;
	//matr[7][8] = (JPP1*sqr(kp_M)-8.*JPP2); 
	//matr[7][9] = (JPM1*sqr(kp_M)-8.*JPM2); 
	//matr[7][10] = -(JJP1*sqr(kk_M)-4.*JJP2)*2.; 
	//matr[7][11] = -(JJM1*sqr(kk_M)-4.*JJM2)*2.;

	//matr[8][8] = 3.*JPP2; 
	//matr[8][9] = 3.*JPM2; 
	//matr[8][10] = (JJP1*sqr(kk_M)-3.*JJP2); 
	//matr[8][11] = (JJM1*sqr(kk_M)-3.*JJM2);

	//matr[9][8] = (JPP1*sqr(kp_M)-2.*JPP2); 
	//matr[9][9] = (JPM1*sqr(kp_M)-2.*JPM2); 
	//matr[9][10] = 2.*JJP2; 
	//matr[9][11] = 2.*JJM2;
///////////////////////////////////////////////////////////

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.; 
	nu_H = (1.5*KH*mu_MH0/mu_M-1.)/(3.*KH*mu_MH0/mu_M+1.);
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system2(matrix, mu_MH0, nu_H);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.;
	mu_MH1 = 2.*mu_M/(3.*KH);

	do { 
		mu_MH1 *= 10.; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
		optim = take_system2(matrix, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0.);

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H  = (1.5*KH*mu_MH/get_param(NUM_SHEAR)-1.)/(3.*KH*mu_MH/get_param(NUM_SHEAR)+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system2(matrix, mu_MH, nu_H);
		if (optim*sgn0 > 0.) mu_MH0 = mu_MH; else
		if (optim*sgn0 < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);
	mu_MH = (mu_MH0+mu_MH1)*.5; 
	
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6 ) mu_MH = -mu_MH;
	return(mu_M/mu_MH);
}

/////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного модуля сдвига через определитель (алгоритм Кристенсена);
double DETER12(double matr[14][14], int i1, int i2)
{	
	int dim_N = 12, ii[12] = {0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
	double  matrix[12][12], deter = 1.;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	if (i1 == 10 && i2 == 11) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[12][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[13][l];
	} else
	if (i1 == 10 && i2 == 12) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[11][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[13][l];
	} else
	if (i1 == 10 && i2 == 13) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[11][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[12][l];
	} else
	if (i1 == 11 && i2 == 12) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[10][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[13][l];
	} else
	if (i1 == 11 && i2 == 13) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[10][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[12][l];
	} else
	if (i1 == 12 && i2 == 13) {
		for (l = 0; l < 12; l++) matrix[10][l] = matr[10][l];
		for (l = 0; l < 12; l++) matrix[11][l] = matr[11][l];
	}
	for (k = 0; k < 10; k++)
	for (l = 0; l < 12; l++) matrix[k][l] = matr[k][l];

////////////////////////////
//...вычисляем определитель;
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
		if (k0 != l0) {
			for (l = 0; l < dim_N; l++) {
				f = matrix[k0][l]; matrix[k0][l] = matrix[l0][l]; matrix[l0][l] = f; 
			}
			deter = -deter;
		}
		if (matrix[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matrix[l0][l0]; deter *= matrix[l0][l0]; matrix[l0][l0] = 1.;
		for (l = 0; l < dim_N; l++) matrix[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matrix[k][l0]; matrix[k][l0] = 0.;
				for (l = 0; l < dim_N; l++) matrix[k][l] -= matrix[l0][l]*finv;
			}
	}
	return deter;
}
double CMindl3D::TakeEshelby_shear_det(double ff)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR2 = 1./pow(ff, 1./3.), fR2 = RR2*RR2, fR5 = ff/fR2, 
			 mu_MI = mu_M/mu_I, AA = get_param(NUM_ADHES)/mu_I, BB = get_param(NUM_ADHES+1)/mu_I,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), A, B, C, D, mu_MH,
			 matr[14][14] = {
					{ mu_MI, -2.*mu_MI*QI1, 0., 0., -1., -1., -3., 2.*QM1, 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -mu_MI,			0., 0., -1., -QM3, 2., 1., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -6.*mu_MI*QI1, 0., 0., -1., 2., 12., 6.*QM1, 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -3.*mu_MI,		0., 0., -1., 2.*QM3, -8., 3., 0., 0., 0., 0., 0., 0.},
					{ AA, -6.*QI1*AA,			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ BB, -3.*BB,				0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 1.,  QI1,					0., 0., -1.,  QM4, 12., -QM1, 0., 0., 0., 0., 0., 0.},
					{ 1., -QI2,					0., 0., -1., -QM5, -8.,  QM2, 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 1.,    ff, 3.*fR5, -2.*QM1*fR2, 0., 0., 0., 0., -3., 1.}, //...нормировку по радиусу в последних коэффициентах можно убрать!!!
					{ 0., 0., 0., 0., 1.,  QM3*ff,  -2.*fR5,    -fR2, 0., 0., 0., 0.,  2., 1.},
					{ 0., 0., 0., 0., 1., -QM4*ff, -12.*fR5, QM1*fR2, 0., 0., 0., 0., 12., 1.},
					{ 0., 0., 0., 0., 1.,  QM5*ff, 8.*fR5,  -QM2*fR2, 0., 0., 0., 0., -8., 1.},
	};
////////////////////////////////////////////
//...дописываем когезионную часть в матрице;
	double kk_I = get_param(NUM_SHEAR+2+NUM_SHIFT), tt_I = exp(-2.*kk_I),
			 kk_M = get_param(NUM_SHEAR+2), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), 
			 JHP1 =  (kk_M*RR2-1.)*ff*tt_D, JHP2 = ((sqr(kk_M*RR2)+3.)-3.*kk_M*RR2)*fR5*tt_D, 
			 JHM1 = -(kk_M*RR2+1.)*ff/tt_D, JHM2 = ((sqr(kk_M*RR2)+3.)+3.*kk_M*RR2)*fR5/tt_D,
			 kp_I = get_param(NUM_SHEAR+3+NUM_SHIFT), tp_I = exp(-2.*kp_I), 
			 kp_M = get_param(NUM_SHEAR+3), tp_D = exp(kp_M*(RR2-1.)),
			 HP1  = ((1.+tp_I)*kp_I-(1.-tp_I))*.5, 
			 HP2  = ((1.-tp_I)*(sqr(kp_I)+3.)-3.*(1.+tp_I)*kp_I)*.5,
			 JPP1 =  (kp_M-1.), JPP2 = ((sqr(kp_M)+3.)-3.*kp_M), 
			 JPM1 = -(kp_M+1.), JPM2 = ((sqr(kp_M)+3.)+3.*kp_M), 
			 JDP1 =  (kp_M*RR2-1.)*ff*tp_D, JDP2 = ((sqr(kp_M*RR2)+3.)-3.*kp_M*RR2)*fR5*tt_D, 
			 JDM1 = -(kp_M*RR2+1.)*ff/tp_D, JDM2 = ((sqr(kp_M*RR2)+3.)+3.*kp_M*RR2)*fR5/tt_D;
////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = -mu_M*HP2*3.; 
	matr[0][3] = -mu_M*(HH1*sqr(kk_I)-3.*HH2);
	matr[0][8] = 3.*mu_M*JPP2; 
	matr[0][9] = 3.*mu_M*JPM2; 
	matr[0][10] = mu_M*(JJP1*sqr(kk_M)-3.*JJP2); 
	matr[0][11] = mu_M*(JJM1*sqr(kk_M)-3.*JJM2);

	matr[1][2] = -mu_M*(HP1*sqr(kp_I)-2.*HP2); 
	matr[1][3] = -mu_M*HH2*2.;
	matr[1][8] = mu_M*(JPP1*sqr(kp_M)-2.*JPP2); 
	matr[1][9] = mu_M*(JPM1*sqr(kp_M)-2.*JPM2); 
	matr[1][10] = mu_M*JJP2*2.; 
	matr[1][11] = mu_M*JJM2*2.;

	matr[2][2] = -mu_M*(HP1*sqr(kp_I)-4.*HP2)*3.; 
	matr[2][3] = -mu_M*((HH2-2.*HH1)*sqr(kk_I)+HH2*12.);
	matr[2][8] = mu_M*(JPP1*sqr(kp_M)-4.*JPP2)*3.; 
	matr[2][9] = mu_M*(JPM1*sqr(kp_M)-4.*JPM2)*3.; 
	matr[2][10] = mu_M*((JJP2-2.*JJP1)*sqr(kk_M)+JJP2*12.); 
	matr[2][11] = mu_M*((JJM2-2.*JJM1)*sqr(kk_M)+JJM2*12.);

	matr[3][2] = -mu_M*((HP2-HP1)*sqr(kp_I)+HP2*8.); 
	matr[3][3] = -mu_M*(HH1*sqr(kk_I)-4.*HH2)*2.;
	matr[3][8] = mu_M*((JPP2-JPP1)*sqr(kp_M)+JPP2*8.); 
	matr[3][9] = mu_M*((JPM2-JPM1)*sqr(kp_M)+JPM2*8.); 
	matr[3][10] = mu_M*(JJP1*sqr(kk_M)-4.*JJP2)*2.; 
	matr[3][11] = mu_M*(JJM1*sqr(kk_M)-4.*JJM2)*2.;

	matr[4][2] = -ku_I*HP2*3.-AA*mu_I*(HP1*sqr(kp_I)-4.*HP2)*3.; 
	matr[4][3] = -ku_I*(HH1*sqr(kk_I)-3.*HH2)-AA*mu_I*((HH2-2.*HH1)*sqr(kk_I)+12.*HH2);
	matr[4][8] = ku_M*JPP2*3.; 
	matr[4][9] = ku_M*JPM2*3.;
	matr[4][10] = ku_M*(JJP1*sqr(kk_M)-3.*JJP2); 
	matr[4][11] = ku_M*(JJM1*sqr(kk_M)-3.*JJM2);

	matr[5][2] = -mu_I*(HP1*sqr(kp_I)-2.*HP2)-BB*mu_I*((HP2-HP1)*sqr(kp_I)+8.*HP2); 
	matr[5][3] = -mu_I*HH2*2.-BB*mu_I*(HH1*sqr(kk_I)-4.*HH2)*2.;
	matr[5][8] = mu_M*(JPP1*sqr(kp_M)-2.*JPP2); 
	matr[5][9] = mu_M*(JPM1*sqr(kp_M)-2.*JPM2);
	matr[5][10] = mu_M*JJP2*2.; 
	matr[5][11] = mu_M*JJM2*2.;

	matr[6][2] = -(mu_I*HP1*sqr(kp_I)-(6.*mu_I+ku_I)*HP2)*3.; 
	matr[6][3] = ((4.*mu_I+ku_I)*HH1*sqr(kk_I)-3.*(6.*mu_I+ku_I)*HH2);
	matr[6][8] = (mu_M*JPP1*sqr(kp_M)-(6.*mu_M+ku_M)*JPP2)*3.; 
	matr[6][9] = (mu_M*JPM1*sqr(kp_M)-(6.*mu_M+ku_M)*JPM2)*3.; 
	matr[6][10] = -((4.*mu_M+ku_M)*JJP1*sqr(kk_M)-3.*(6.*mu_M+ku_M)*JJP2); 
	matr[6][11] = -((4.*mu_M+ku_M)*JJM1*sqr(kk_M)-3.*(6.*mu_M+ku_M)*JJM2);

	matr[7][2] = (2.*mu_I*HP1*sqr(kp_I)-(10.*mu_I-3.*ku_I)*HP2)*2.; 
	matr[7][3] = ((ku_I-2.*mu_I)*HH1*sqr(kk_I)+(10.*mu_I-3.*ku_I)*HH2)*2.;
	matr[7][8] = -(2.*mu_M*JPP1*sqr(kp_M)-(10.*mu_M-3.*ku_M)*JPP2)*2.; 
	matr[7][9] = -(2.*mu_M*JPM1*sqr(kp_M)-(10.*mu_M-3.*ku_M)*JPM2)*2.; 
	matr[7][10] = -((ku_M-2.*mu_M)*JJP1*sqr(kk_M)+(10.*mu_M-3.*ku_M)*JJP2)*2.; 
	matr[7][11] = -((ku_M-2.*mu_M)*JJM1*sqr(kk_M)+(10.*mu_M-3.*ku_M)*JJM2)*2.;

	matr[8][8] = 3.*JDP2; 
	matr[8][9] = 3.*JDM2; 
	matr[8][10] = (JHP1*sqr(kk_M)-3.*JHP2); 
	matr[8][11] = (JHM1*sqr(kk_M)-3.*JHM2);

	matr[9][8] = (JDP1*sqr(kp_M)-2.*JDP2); 
	matr[9][9] = (JDM1*sqr(kp_M)-2.*JDM2); 
	matr[9][10] = 2.*JHP2; 
	matr[9][11] = 2.*JHM2;

/////////////////////////////
//...вычисляем олпределители;
	A = DETER12(matr, 10, 11)*5.;
	B = DETER12(matr, 10, 12)*(-15.)+DETER12(matr, 10, 13)*(-5.)+DETER12(matr, 11, 12)*10.+DETER12(matr, 11, 13)*10.;
	C = DETER12(matr, 12, 13)*(-20.); 
	
/////////////////////////////////
//...решаем квадратное уравнение;
	if (A == 0.) mu_MH = B ? -C/B : -1.; else {	
		B /= A; C /= A;
		D = B*B*.25-C;
		if (D < 0.) return(0.); else	
		//if (B < 0. && C > 0.) mu_MH = -B*.5-sqrt(D);
		//else	mu_MH = -B*.5+sqrt(D);
		//mu_MH = -B*.5-sqrt(D);
		mu_MH = -B*.5+sqrt(D);
	}
	return(mu_M/mu_MH);
}
#undef  Message
