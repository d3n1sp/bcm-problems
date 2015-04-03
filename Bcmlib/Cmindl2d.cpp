#include "stdafx.h"

#include "cgrid_el.h"
#include "cmindl2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CMindl2D::NUM_ADHES = 3;
int CMindl2D::NUM_SHEAR = 5;
int CMindl2D::NUM_SHIFT = 4;
int CMindl2D::MAX_PHASE = 2;
int CMindl2D::NUM_HESS  = 8;
int CMindl2D::regul = 1;
int  gradient_model = 1;

//////////////////////////////////
//...initialization of the blocks;
int  CMindl2D::block_shape_init(Block<double> & B, int id_free)
{
	int k;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		B.shape->add_shape(CreateShapeD(MP2D_POLY_SHAPE));

		B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

////////////////////////////////////////////////////////////
//...подключаем регулярную экспоненциальную систему функций;
		if (gradient_model)	{//...using gradient displacements;
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->add_shape(CreateShapeD(SK2D_ELLI_SHAPE), NULL_STATE);
			else											   B.shape->add_shape(CreateShapeD(SK2D_POLY_SHAPE), NULL_STATE);
		}
		if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 

/////////////////////////////////////////////
//...задаем степень аппроксимирующих функций;
		B.shape->init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
		if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for inclusion!!!
		B.shape->init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, draft_dim(type()));

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
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
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
void CMindl2D::jump1_common_x(double * P, int i, int m)
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
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (; l < nm*2; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0, 1.);

		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 0., 1.);
	}
}

////////////////////////////////////////////////////////////////////
//...additional realization of common compression displacement (ux);
void CMindl2D::jump1_compress_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 1., -1.);

		B[i].shape->cpy_xy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 1., -1.);
	}
}

//////////////////////////////////////////////
//...realization of common displacements (Ry);
void CMindl2D::jump1_common_y(double * P, int i, int m)
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
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (; l < nm*2; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 0., 1.);

		B[i].shape->cpy_yy(l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0, 1.);
	}
}

////////////////////////////////////////////////////////////////////
//...additional realization of common compression displacement (uy);
void CMindl2D::jump1_compress_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(); m += solver.id_norm;
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), B[i].shape->deriv, 1., -1.);

		B[i].shape->cpy_yy(l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->deriv, 1., -1.);
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (DnRx);
void CMindl2D::jump2_common_x(double * P, int i, int m)
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

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*(2.*alpha+1.)*G0, P[0]*alpha*G0);
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.)*G0);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*alpha*G0, P[1]*alpha*G0);
	}
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l), ptr, 0., 1.);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);
	}
}

////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of common compression displacement (Dnux);
void CMindl2D::jump2_compress_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm;
	double * ptr; m += solver.id_norm;
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (l = nm; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 1., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 1., -1.);
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (DnRy);
void CMindl2D::jump2_common_y(double * P, int i, int m)
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

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*(2.*alpha+1.)*G0, P[1]*alpha*G0);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.)*G0);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*alpha*G0, P[0]*alpha*G0);
	}
///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., 1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[3]*C0);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);
	}
}

////////////////////////////////////////////////////////////////////////////////
//...realization of normal derivative of common compression displacement (Dnuy);
void CMindl2D::jump2_compress_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), num_hess = NUM_HESS+solver.id_norm;
	double * ptr; m += solver.id_norm;

///////////////////////////////////////////////////
//...гельмгольцевская часть нормальных производных;
	for (l = nm; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 1., -1.);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_yy;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 1., -1.);
	}
}

///////////////////////////////////////////
//...realization of cohesion moment M(n)_x;
void CMindl2D::jump3_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), 
			 C0 = get_param(NUM_SHEAR+3+shift), fnx = (G0+alpha*sqr(P[3])), * ptr;
/////////////////////////////////////////////////////////
//...лапласовская часть когезионного момента (обнуление);
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, 0., 0.);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*fnx, fnx);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[3]*P[4]);

		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[4], fnx);
		B[i].shape->adm_yy	 (l, ptr, alpha*P[3]*P[4]);
	}
}

///////////////////////////////////////////////////////////////
//...realization of compression part of cohesion moment M(n)_x;
void CMindl2D::jump3_compress_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), fnx = (G0+alpha*sqr(P[3])), * ptr;
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->adm_xx(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -fnx);
		B[i].shape->adm_xy(l, ptr, -alpha*P[3]*P[4]);

		B[i].shape->adm_xy(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -fnx);
		B[i].shape->adm_yy(l, ptr, -alpha*P[3]*P[4]);
	}
}

///////////////////////////////////////////
//...realization of cohesion moment M(n)_y;
void CMindl2D::jump3_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), 
			 C0 = get_param(NUM_SHEAR+3+shift), fny = (G0+alpha*sqr(P[4])), * ptr;
/////////////////////////////////////////////////////////
//...лапласовская часть когезионного момента (обнуление);
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 0), NULL, 0., 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, 0., 0.);
	}
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->cpy_xy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*alpha*P[3]*P[4], fny);
		B[i].shape->adm_xx	 (l, ptr, alpha*P[3]*P[4]);

		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, -C0*fny, fny);
		B[i].shape->adm_xy	 (l, ptr, alpha*P[3]*P[4]);
	}
}

///////////////////////////////////////////////////////////////
//...realization of compression part of cohesion moment M(n)_y;
void CMindl2D::jump3_compress_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift), alpha = G0/(1.-get_param(NUM_SHEAR+1+shift)*2.), fny = (G0+alpha*sqr(P[4])), * ptr;
////////////////////////////////////////
//...гельмгольцевская часть перемещений;
	for (l = nm; l < nm*2; l++) {
		B[i].shape->adm_xy(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -fny);
		B[i].shape->adm_xx(l, ptr, -alpha*P[3]*P[4]);

		B[i].shape->adm_yy(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -fny);
		B[i].shape->adm_xy(l, ptr, -alpha*P[3]*P[4]);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical surface tensions for common displacements Px(R);
void CMindl2D::jump4_common_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 = get_param(NUM_SHEAR+3+shift)*G0*.5, * ptr; m += solver.id_norm;
/////////////////////////////////////////////////
//...лапласовская часть поверхностных напряжений;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xx    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3], P[0]*alpha);
		B[i].shape->adm_y     (l, ptr,  P[4]*(alpha+1.));

		B[i].shape->cpy_x     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4]*(alpha+1.), P[1]*alpha);
		B[i].shape->adm_y     (l, ptr,  P[3]*(-2.*alpha-1.));
	}
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C0*2.);
		B[i].shape->adm_y(l, ptr, -P[4]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C0);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of compression part of classical surface tensions for common displacements Px(R);
void CMindl2D::jump4_compress_x(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., C1 = get_param(NUM_SHEAR+2+shift)*G0*get_param(NUM_SHEAR+1+shift)/(1.-get_param(NUM_SHEAR+1+shift)*2.), * ptr;
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (l = nm; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l); //...deriv_N_xx;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 1., -G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 1., -G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[3]*C1);
	}
}

//////////////////////////////////////////////////////////////////////////////
//...realization of classical surface tensions for common displacements Py(R);
void CMindl2D::jump4_common_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), 
			 C0 = get_param(NUM_SHEAR+3+shift)*G0*.5, * ptr; m += solver.id_norm;
/////////////////////////////////////////////////
//...лапласовская часть поверхностных напряжений;
	for (l = 0; l < nm; l++) {
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, 0., 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_yy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[3]*(alpha+1.), P[0]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[4]*(-2.*alpha-1.));

		B[i].shape->cpy_y     (l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1));
		B[i].shape->admittance(l, ptr,  B[i].shape->deriv, P[4], P[1]*alpha);
		B[i].shape->adm_x     (l, ptr,  P[3]*(alpha+1.));
	}
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[3]*C0);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., 1.);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C0*2.);
		B[i].shape->adm_x(l, ptr, -P[3]*C0);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...realization of compression part of classical surface tensions for common displacements Py(R);
void CMindl2D::jump4_compress_y(double * P, int i, int m)
{
	int l, nm = B[i].shape->num_usual(), shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm; m += solver.id_norm;
	double G0 = get_param(NUM_SHEAR+shift)*2., C1 = get_param(NUM_SHEAR+2+shift)*G0*get_param(NUM_SHEAR+1+shift)/(1.-get_param(NUM_SHEAR+1+shift)*2.), * ptr;
/////////////////////////////////////////////////////
//...гельмгольцевская часть поверхностных напряжений;
	for (l = nm; l < nm*2; l++) {
		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1); //...deriv_N_xy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 0), ptr, 0., -G0);
		B[i].shape->adm_x(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 0), -P[4]*C1);

		ptr = B[i].shape->FULL(solver.hh[i][0][num_hess+1], l); //...deriv_N_yy;
		B[i].shape->admittance(l,  B[i].shape->FULL(solver.hh[i][0][m], l, 1), ptr, 0., -G0);
		B[i].shape->adm_y(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -P[4]*C1);
	}
}

////////////////////////////////////////////////////////
//...realization partial normal derivatives of hessian;
void CMindl2D::hessian_deriv_N(int k, double * P, int i)
{
	int num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->set_norm_cs(P);

	B[i].shape->cpy_xx(k);	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xx(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k));

	B[i].shape->cpy_xy(k); B[i].shape->deriv_N(k);
	B[i].shape->cpy_xy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k, 1));

	B[i].shape->cpy_yy(k); B[i].shape->deriv_N(k);
	B[i].shape->cpy_yy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k));
}

///////////////////////////////////////
//...composition of collocation vector;
void CMindl2D::jump_admittance(int l, int i, int m, double adm_re, int k, double adm_im)
{
	m += solver.id_norm;
	if (m >= 0 && adm_re != 1.) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),    NULL, adm_re, 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, adm_re, 0.);
	}
	if (k >= 0) {
		k += solver.id_norm;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),	 B[i].shape->FULL(solver.hh[i][0][k], l),	   1., adm_im);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][k], l, 1), 1., adm_im);
	}
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CMindl2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0. };
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CMindl2D::jump_make_common(int i, int m)
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

/////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии (условие периодического скачка);
int CMindl2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx;
      int	 m  = solver.id_norm, id_dir, k, j, shift;

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
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

////////////////////////////////////////////////////////////////////////////////////
//...условие скачка для классической составляющей поля методом наименьших квадратов;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+4],  hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+5],  hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m+4], -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+5], -hx*f);
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...граничное условие для классической составляющей поля перемещений (регуляризация в случае симметричного включения);
				if (regul) {
					if (id_dir == 1 || id_dir == 2) {
						solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
						solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*f*.5);

						solver.to_equationDD(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_equationHH(k, 0, solver.hh[k][0][m+4], -hx*f*.5);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
						solver.to_equationDD(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], f);
					}
					//аналогичная регуляризация второй задачи???
				}

//////////////////////////////////////////////////////////
//...сшивка когезионого поля методом наименьших квадратов;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
	
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе фаз;
int CMindl2D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], Ai, Ak, Bi, Bk;
      int m = solver.id_norm, shift;

////////////////////////////////////////
//...коэффициенты поверхностной адгезии;
		Ai = Ak = get_param(NUM_ADHES)*.5-(Bi = Bk = get_param(NUM_ADHES+1)*.5);

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
				f = nd->get_param(0, l);
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, i, 4);
					jump2_common_y (P, i, 5);
				}
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, i, 4); solver.admittance(i, 4, G1);
					jump2_compress_y(P, i, 5); solver.admittance(i, 5, G1); 
				}
				jump_make_common(i, 0);
				jump_make_common(i, 2);
				jump_make_common(i, 4);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, k, 4);
					jump2_common_y (P, k, 5);
				}
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, k, 4); solver.admittance(k, 4, G1); 
					jump2_compress_y(P, k, 5); solver.admittance(k, 5, G1); 
				}
				jump_make_common(k, 0);
				jump_make_common(k, 2);
				jump_make_common(k, 4);

//////////////////////////////////////////////////////////////////////////
//...сшивка функций и нормальных производных методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////
//...добавляем поверхностную энергию адгезии (нормированную G1);
				if (gradient_model)	{//...using gradient displacements;
					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], -f*(Ai*sqr(P[3])+Bi*sqr(P[4])));
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], -f*(Ak*sqr(P[3])+Bk*sqr(P[4])));

					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+5], -f*Ai*Bi*P[3]*P[4]*2.);
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+5], -f*Ak*Bk*P[3]*P[4]*2.);

					solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], -f*(Ai*sqr(P[4])+Bi*sqr(P[3])));
					solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], -f*(Ak*sqr(P[4])+Bk*sqr(P[3])));
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && B[i].link[NUM_PHASE] < -1) f = -f;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
int CMindl2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;

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

/////////////////////////////////////////////
//...jump of all neaded displacement moments;
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
			B[i].shape->parametrization_hess(P, 1);

			hessian_deriv_N(1, P, i);
			jump1_common_x (P, i, 0);
			jump1_common_y (P, i, 1);

			jump4_common_x (P, i, 2);
			jump4_common_y (P, i, 3);

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_hess(1, P, 1);

			hessian_deriv_N (1, P, i);
			jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
			jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

			jump4_compress_x(P, i, 2); 
			jump4_compress_y(P, i, 3); 

			jump_make_common(i, 0);
			jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов (классическая составляющая);
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.); 

				if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);

				if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], hy*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}

///////////////////////////////////////////////////
//...обнуление когезионного поля (на всей границе);
			solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
			solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.); 
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
			solver.to_equationEE(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
int CMindl2D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
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

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

/////////////////////////////////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов (классическая составляющая);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов (когезионное поле);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////
//...интегрирование НДС на заданном наборе узлов для периодической задачи;
int CMindl2D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6], UX, UY, RX, RY;

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

///////////////////////////////////////////////////////////////////////////////////////////////
//...формирование интеграла от полных деформаций и классических напряжений (два состояния НДС);
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 			B[i].shape->parametrization_grad(0, P, 1);
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_common_x(P, i, 0); 
			jump1_common_y(P, i, 1); 

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_compress_x(P, i, 0);
			jump1_compress_y(P, i, 1);
			jump_make_common(i, 0); //...переходим к общей системе координат;

//////////////////////////////////
//...вычисляем полные перемещения;
			RX = B[i].shape->potential(solver.hh[i][0][m],   0);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[3] -= RX*nd->nX[l]*f;
			K[4] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[5] -= RY*nd->nY[l]*f;

			RX = B[i].shape->potential(solver.hh[i][0][m],   1);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[9]  -= RX*nd->nX[l]*f;
			K[10] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[11] -= RY*nd->nY[l]*f;

/////////////////////////////////////////////////////////////////
//...вычисляем классические перемещения и интегрируем напряжения;
			jump_admittance(1, i, 0, 0.);
			jump_admittance(1, i, 1, 0.);

			UX = B[i].shape->potential(solver.hh[i][0][m],   0);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[0] -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[1] -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[2] -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;

			K[12] -= G0* UX*nd->nX[l]*f*alpha;
			K[13] -= G0*(UX*nd->nX[l]+UY*nd->nY[l])*f*alpha;

			UX = B[i].shape->potential(solver.hh[i][0][m],   1);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[6]  -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[7]  -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[8]  -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;

/////////////////////////////////
//...вычисляем объем каждой фазы;
			K[14+(-B[i].link[NUM_PHASE]-1)] -= (P[0]*P[3]+P[1]*P[4])*f*.5;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CMindl2D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B)
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

/////////////////////////////////////////////////////////////
//...counting kernel for solving gradient elasticity problem;
int CMindl2D::computing_header(Num_Comput Num)
{
	int i, j, k, elem, id_dir, n_rhs = 2;
	char msg[201];

	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CMindl2D sample: N_sm = %d, N_mpl = %d, l1 = %g, l2 = %g", N,
				  UnPackInts(get_param(NUM_MPLS)), 1./get_param(NUM_SHEAR+3), 1./get_param(NUM_SHEAR+NUM_SHIFT+3));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

//////////////////////////////////////////////
//...устанавливаем мультиполи межфазного слоя;
	if (NUM_PHASE > 0)
	for (k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE)
	for (j = B[k].link[0]; j > 0; j--) 
	if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE])
	if (/*B[k].link[NUM_PHASE] == -1 &&*/ B[k].mp) { 
		int m, m1, m2;
		m  = stru.geom[(i = stru.geom_ptr[k])+1]-2;
		m1 = stru.geom[i+4+(j-1)%m];
		m2 = stru.geom[i+4+j%m];
		double X0 = stru.X[m1], Y0 = stru.Y[m1], 
				 X1 = stru.X[m2], Y1 = stru.Y[m2];
		B[k].type = ELLI_BLOCK;
		B[k].mp[4] = arg0(comp(X1-X0, Y1-Y0))-M_PI_2;
		B[k].mp[5] = 0.;
		B[k].mp[6] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials!!!
	solver.n += 10;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par); d_parm[0] = par[1]-par[0]; d_parm[1] = par[3]-par[2];
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
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : OK_STATE) || 
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
void CMindl2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp && B[i].link[0] >= NUM_PHASE) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case DISPL_VALUE: {
//////////////////////////
//...common displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);
//				jump_make_common(i, 0); //...вместо B[i].shape->norm_common(F);

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				if ((solv%ENERGY_SOLVING) && 1) F[(solv%ENERGY_SOLVING)-1] -=X;
			}  break;
			case DISPL_CLASSIC_VALUE: {
//////////////////////////////////
//...displacement (classic field);
				B[i].shape->parametrization_grad(0, P, 1);

				jump1_common_x(P, i, 0); jump_admittance(1, i, 0, 0.);
				jump1_common_y(P, i, 1); jump_admittance(1, i, 1, 0.);

////////////////////////////////////////
//...calculation classical displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				if ((solv%ENERGY_SOLVING) && 0) F[(solv%ENERGY_SOLVING)-1] -=X;
			}  break;
			case DISPL_COHESION_VALUE: {
///////////////////////////////////
//...displacement (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_common_x (P, i, 0); jump_admittance(0, i, 0, 0.);
				jump1_common_y (P, i, 1); jump_admittance(0, i, 1, 0.);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);

///////////////////////////////////////
//...calculation cohesion displacement;
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_CLASSIC_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////
//...surface stresses (classic field);
				B[i].shape->parametrization_hess(0, P, 1);

				jump4_common_x(P, i, 0); jump_admittance(1, i, 0, 0.);
				jump4_common_y(P, i, 1); jump_admittance(1, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_CLASSIC_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////
//...surface stresses (classic field);
				B[i].shape->parametrization_hess(0, P, 1);

				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_COHESION_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_COHESION_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case NORMAL_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case NORMAL_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case DILAT_VALUE: {//...divergency of displacements;
				double G = (.5-get_param(NUM_SHEAR+1+shift))/((1.-get_param(NUM_SHEAR+1+shift))*get_param(NUM_SHEAR+shift));
///////////////////////////////////////
//...classical and pressure potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), 1.); 
				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), -1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

/////////////////////////////////
//...divergency of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ROTATION_VALUE: {//...rotation of displacements;
				double G = .5/get_param(NUM_SHEAR+shift);
////////////////////////////////////
//...classical and shear potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), -1.); 
				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

///////////////////////////////
//...rotation of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ENERGY_VALUE: {//...local energy in cell;
				double e11, e12, e22, theta, u1, u2, G1 = get_param(NUM_SHEAR+shift), 
						 G2 = get_param(NUM_SHEAR+1+shift)*G1/(.5-get_param(NUM_SHEAR+1+shift)), C1 = get_param(NUM_SHEAR-1+shift);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				jump1_common_x (P, i, 4); jump_admittance(0, i, 4, 0.);
				jump1_common_y (P, i, 5); jump_admittance(0, i, 5, 0.);

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем дилатацию, компоненты тензора деформации и когезионное поле (с учетом поворота компонент);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				e11 =	F[0];
				e12 =	F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F);

				e22 = F[1];
				e12 =(F[0]+e12)*.5;
				theta = e11+e22;

				u1 = B[i].shape->potential(solver.hh[i][0][m+4], id_variant);
				u2 = B[i].shape->potential(solver.hh[i][0][m+5], id_variant);

//////////////////
//...local energy;
				F[0] = G1*(e11*e11+2.*e12*e12+e22*e22)+(G2*theta*theta+C1*(u1*u1+u2*u2))*.5;
			}	break;
			case STRAIN_X_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F); double exy = F[0];

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				F[1] = (F[1]+exy)*.5;
			}  break;
			case STRAIN_Y_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F); double exy = F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F);

				F[0] = (F[0]+exy)*.5;
			}  break;
			default : F[0] = i; F[1] = 0.;
     }
  }
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CMindl2D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k][0].GetMatrix() && 
							solver.TL && solver.TL[k]) {
		double ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]-solver.JR_SHIFT].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
		int i, j, NN = solver.dim[k];

///////////////////////
//...вычисляем энергию;
		block_shape_init(B[k], NO_STATE);
		if (h && TL) { 
			double measure = .5/get_param(NUM_SHEAR);
			for (i = 0; i < NN; i++)
			for (j = 0; j < NN; j++)
				energy[-B[k].link[NUM_PHASE]-1] += TL[i][j]*h[i]*h[j]*measure;
		}
	}
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (прямой алгоритм);
double CMindl2D::TakeEshelby_volm_two(double ff)
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
					{ 3.*K_I, (ku_I+2.*mu_I)*HH1, -3.*K_M, 4.*mu_M, -(ku_M+2.*mu_M)*JJP1, -(ku_M+2.*mu_M)*JJM1, 0., 0.},
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
extern double take_system(double matrix[14][15], double mu_MH, double nu_H);

//////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (итерационный алгоритм);
double CMindl2D::TakeEshelby_shear_two(double ff, double eps, int max_iter)
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
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, sgn0, sgn1, nu_H, mu_MH, mu_MH0, mu_MH1, matrix[14][15],
			 matr[14][15] = {
					{ mu_MI, -2.*mu_MI*QI1, 0., 0., -1., -1., -3., 2.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -mu_MI,			0., 0., -1., -QM3, 2., 1., 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -6.*mu_MI*QI1, 0., 0., -1., 2., 12., 6.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ mu_MI, -3.*mu_MI,		0., 0., -1., 2.*QM3, -8., 3., 0., 0., 0., 0., 0., 0., 0.},
					{ AA, -6.*QI1*AA,			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ BB, -3.*BB,				0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 1.,  QI1,					0., 0., -1.,  QM4, 12., -QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -QI2,					0., 0., -1., -QM5, -8.,  QM2, 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 1., ff, 3.*fR5, -2.*QM1*fR2, 0., 0., 0., 0., -3., -1., 1.}, //...нормировку по радиусу в последних коэффициентах можно убрать!!!
					{ 0., 0., 0., 0., 1.,  QM3*ff, -2.*fR5, -fR2,  0., 0., 0., 0.,  2., -1., 1.},
					{ 0., 0., 0., 0., 1., -QM4*ff, -12.*fR5, QM1*fR2, 0., 0., 0., 0., 12.,  1., 1.},
					{ 0., 0., 0., 0., 1.,  QM5*ff, 8.*fR5,  -QM2*fR2, 0., 0., 0., 0., -8., -1. , 1.},
	};
	int k_iter = 0, k, l;
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
	matr[6][10] = -((4.*mu_M+ku_M)*JJP1*sqr(kk_I)-3.*(6.*mu_M+ku_M)*JJP2); 
	matr[6][11] = -((4.*mu_M+ku_M)*JJM1*sqr(kk_I)-3.*(6.*mu_M+ku_M)*JJM2);

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

	mu_MH0 = 0.; 
	nu_H = (1.5*KH*mu_MH0/mu_M-1.)/(3.*KH*mu_MH0/mu_M+1.);
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system(matrix, mu_MH0, nu_H);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH1 = 2.*mu_M/(3.*KH);
	do { 
		mu_MH1 *= 10.; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0.);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H  = (1.5*KH*mu_MH/get_param(NUM_SHEAR)-1.)/(3.*KH*mu_MH/get_param(NUM_SHEAR)+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_MH, nu_H);
		if (optim*sgn0 > 0.) mu_MH0 = mu_MH; else
		if (optim*sgn0 < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);
	mu_MH = (mu_MH0+mu_MH1)*.5; 
	
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6 ) mu_MH = -1.;
	return(mu_M/mu_MH);
}
#undef  Message
