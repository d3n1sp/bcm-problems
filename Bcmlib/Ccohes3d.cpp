#include "stdafx.h"

#include "cgrid_el.h"
#include "ccohes3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

//////////////////////////////////
//...initialization of the blocks;
int  CCohes3D::block_shape_init(Block<double> & B, int id_free)
{
	int k, m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == POLY_BLOCK) {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->add_shape(CreateShapeD(SK3D_ZOOM_SHAPE));
			B.shape->set_shape(get_param(1)*fabs( B.mp[7]));
		}
		else
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShapeD(MP3D_ZOOM_SHAPE));
			B.shape->add_shape(CreateShapeD(SK3D_EXPP_SHAPE, 1));
			B.shape->set_shape(get_param(1)*fabs( B.mp[8]));
		}
		else {
			B.shape->add_shape(CreateShapeD(MP3D_POLY_SHAPE));
			B.shape->set_shape(get_param(1)*fabs( B.mp[7]));
		}

////////////////////////
//...setting parameters;
      B.shape->init1(UnPackInts(get_param(0)), solver.id_norm, 3);
		if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for inclusion!!!
		B.shape->init1(UnPackInts(get_param(0), 1), solver.id_norm, 3);

//////////////////////////////////////////////////////////////////////////////
//...setting acselerator, local system of coordinate and init parametrization;
		B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////
//...setting parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
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

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CCohes3D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2)
{
	if (size_of_param() > NUM_SHEAR+4+NUM_SHIFT) {
      param[NUM_SHEAR-1] = C1;
      param[NUM_SHEAR-1+NUM_SHIFT] = C2;
      param[NUM_SHEAR] = G1;
      param[NUM_SHEAR+NUM_SHIFT] = G2;
      param[NUM_SHEAR+1] = nju1;
      param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
      param[NUM_SHEAR+2] = .25/(nju1-1.);
      param[NUM_SHEAR+2+NUM_SHIFT] = .25/(nju2-1.);;
      param[NUM_SHEAR+3] = sqrt((1.-2.*nju1)*.5/(1.-nju1)*C1/G1);
      param[NUM_SHEAR+3+NUM_SHIFT] = sqrt((1.-2.*nju2)*.5/(1.-nju2)*C2/G2);
      param[NUM_SHEAR+4] = sqrt(C1/G1);
      param[NUM_SHEAR+4+NUM_SHIFT] = sqrt(C2/G2);
	}
}

////////////////////////////////////////////////////////////////////////
//...вычисление эффективных характеристик в градиентной слоистой модели;
double CCohes3D::TakeLayer_kk(int N, double * ff, double * kk, double * ll)
{
	double ** matr = NULL; set_matrix(matr, 4*N, 4*N+1);
	int i1, i2, m1, m2, m;

//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	for ( m2 = 4, i2 = 1, m = m1 = i1 = 0; i1 < N; i1++, i2 = (i1+1)%N, m1 = 4*i1, m2 = 4*i2) {
		matr[m][m1] = ff[i1]; matr[m][m1+1] =  1.0; matr[m][m1+2] = exp(ff[i1]/ll[i1]); matr[m][m1+3] = exp(-ff[i1]/ll[i1]); 
									 matr[m][m2+1] = -1.0; matr[m][m2+2] = -1.0; matr[m][m2+3] = -1.; matr[m][4*N] = i2 ? 0. : 1.; m++;

		matr[m][m1] =  1.0; matr[m][m1+2] = exp(ff[i1]/ll[i1])/ll[i1]; matr[m][m1+3] = -exp(-ff[i1]/ll[i1])/ll[i1]; 
		matr[m][m2] = -1.0; matr[m][m2+2] = -1./ll[i2]; matr[m][m2+3] = 1./ll[i2]; m++;

		matr[m][m1+2] = exp(ff[i1]/ll[i1])*kk[i1]; matr[m][m1+3] = exp(-ff[i1]/ll[i1])*kk[i1]; 
		matr[m][m2+2] = -kk[i2]; matr[m][m2+3] = -kk[i2]; m++;

		if (i2) {
			matr[m][m1] = kk[i1]; matr[m][m2] = -kk[i2]; m++;
		}
	}
	for ( m1 = i1 = 0; i1 < N; i1++, m1 = 4*i1) { //...добавляем условие среднего по периоду;
		matr[m][m1] = 0.5*ff[i1]*ff[i1];	matr[m][m1+1] = ff[i1];
	}

///////////////////////////////////////
//...решаем систему линейных уравнений;
	double sum = 0.;
	int dim_N = 4*N, * ii, i, k, l, k0, l0; ii = (int *)new_struct(4*N*sizeof(int));
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
	for (i = 0; i < N; i++) sum += ff[i]*kk[i]*matr[4*i][dim_N];
err:
	delete_struct(matr); delete_struct(ii);
	return(sum);
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (прямой алгоритм);
double CCohes3D::TakeEshelby_volm_two(double c0)
                                {
	real_E mu_I = real_E(get_param(NUM_SHEAR+NUM_SHIFT)), nu_I = real_E(get_param(NUM_SHEAR+1+NUM_SHIFT)),
			 mu_M = real_E(get_param(NUM_SHEAR)), nu_M = real_E(get_param(NUM_SHEAR+1)),
			 C_I = real_E(get_param(NUM_SHEAR-1+NUM_SHIFT)), ff = real_E(c0), 
			 C_M = real_E(get_param(NUM_SHEAR-1)), AA = real_E(get_param(NUM_GEOMT)),
			 K_I = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K_M = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M);	if (C_I == 0. || C_M == 0.) { //...классическая трехфазная модель;
		return to_double(K_M+ff*(K_I-K_M)/(1.+(1.-ff)*(K_I-K_M)/(K_M+4./3.*mu_M)));
	}
	real_E ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, RR2 = 1./pow(ff, 1./real_E(3.)), 
			 kk_I = sqrt(C_I/ku_I), tt_I = exp(-2.*kk_I),
			 kk_M = sqrt(C_M/ku_M), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), JHP1 =  (kk_M*RR2-1.)*ff*tt_D, 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), JHM1 = -(kk_M*RR2+1.)*ff/tt_D,
			 matr[7][8] = {
					{ 1., -HH1, -1., -1., JJP1, JJM1, 0., 0.},
					{ 0., -HH2,  0.,  3., JJP2, JJM2, 0., 0.},
					{ AA, -ku_I*HH1-AA*(HH2+HH1), 0.,0., ku_M*JJP1, ku_M*JJM1, 0., 0.},
//...новый вариант: антисимметричный, симметрия по второму и третьему индексам;
					{ K_I, (mu_I+2.*ku_I)/3.*HH1, -K_M, 4.*mu_M/3., -(mu_M+2.*ku_M)/3.*JJP1, -(mu_M+2.*ku_M)/3.*JJM1, 0., 0.},
//...прежний вариант: симметрия по первым двум индексам;
					//{ K_I, (4.*mu_I+ku_I)/3.*HH1, -K_M, 4.*mu_M/3., -(4.*mu_M+ku_M)/3.*JJP1, -(4.*mu_M+ku_M)/3.*JJM1, 0., 0.},
					{ 0., 0., 0., 0., JHP1, JHM1,  0., 0.},
					{ 0., 0., K_M, -4.*mu_M/3.*ff, 0., 0., -1., 0.}, //...эффективный модуль с коэффициентом один;
					{ 0., 0., 1., ff, 0., 0., 0., 1.},					 //...перемещения дают единицу в правую часть;
	};

//////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 7, ii[7] = {0, 0, 0, 0, 0, 0, 0}, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		real_E f = 0.;
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
		real_E finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows; 
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matr[k][l0]; matr[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
			}
	}
	return to_double(matr[6][7]);
}


/////////////////////////////////////////////////////////////////////////////////////////////
//...трехфазная (симметричная или некорректная дилатансия) модель для сферического включения;
double CCohes3D::TakeEshelby_volm_sym(double ff)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), 
			 C_M = get_param(NUM_SHEAR-1), AA = get_param(NUM_GEOMT),
			 K_I = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K_M = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M);
	if (C_I == 0. || C_M == 0.) { //...классическая трехфазная модель;
		return to_double(K_M+ff*(K_I-K_M)/(1.+(1.-ff)*(K_I-K_M)/(K_M+4./3.*mu_M)));
	}
	double ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, RR2 = 1./pow(ff, 1./3.), 
			 kk_I = sqrt(C_I/ku_I), tt_I = exp(-2.*kk_I),
			 kk_M = sqrt(C_M/ku_M), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), JHP1 =  (kk_M*RR2-1.)*ff*tt_D, 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), JHM1 = -(kk_M*RR2+1.)*ff/tt_D,
			 matr[7][8] = {
					{ 1., -HH1, -1., -1., JJP1, JJM1, 0., 0.},
					{ 0., -HH2,  0.,  3., JJP2, JJM2, 0., 0.},
					{ AA, -ku_I*HH1-AA*(HH2+HH1), 0.,0., ku_M*JJP1, ku_M*JJM1, 0., 0.},
//...новый вариант: антисимметричный, симметрия по второму и третьему индексам;
					//{ K_I, (mu_I+2.*ku_I)/3.*HH1, -K_M, 4.*mu_M/3., -(mu_M+2.*ku_M)/3.*JJP1, -(mu_M+2.*ku_M)/3.*JJM1, 0., 0.},
//...прежний вариант: симметрия по первым двум индексам;
					{ K_I, (4.*mu_I+ku_I)/3.*HH1, -K_M, 4.*mu_M/3., -(4.*mu_M+ku_M)/3.*JJP1, -(4.*mu_M+ku_M)/3.*JJM1, 0., 0.},
//...новый вариант: антисимметричный, нет симметрии по второму и третьему индексам;
					//{ K_I, (3.*mu_I+ku_I)/3.*HH1, -K_M, 4.*mu_M/3., -(3.*mu_M+ku_M)/3.*JJP1, -(3.*mu_M+ku_M)/3.*JJM1, 0., 0.},
					{ 0., 0., 0., 0., JHP1, JHM1,  0., 0.},
					{ 0., 0., K_M, -4.*mu_M/3.*ff, 0., 0., -1., 0.}, //...эффективный модуль с коэффициентом один;
					{ 0., 0., 1., ff, 0., 0., 0., 1.},					 //...перемещения дают единицу в правую часть;
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
	return(matr[6][7]);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для трехфазной модели (модуль сдвига);
double take_system(double matrix[14][15], double mu_HM, double nu_H)
{
	int dim_N = 14, ii[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;	
	matrix[11][13] *= 2.*(1.-2.*nu_H)/(5.-4.*nu_H);
	matrix[12][12] *= mu_HM;
	matrix[12][13] *= mu_HM*2.*(5.-nu_H)/(5.-4.*nu_H);
	matrix[12][14] *= mu_HM;
	matrix[13][12] *= mu_HM;
	matrix[13][13] *= mu_HM*2.*(1.+nu_H)/(5.-4.*nu_H);
	matrix[13][14] *= mu_HM;
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

/////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби с повышенной точностью;
real_E take_system_E(real_E matrix[14][15], real_E mu_HM, real_E nu_H)
{
	int dim_N = 14, ii[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;	
	matrix[11][13] *= 2.*(1.-2.*nu_H)/(5.-4.*nu_H);
	matrix[12][12] *= mu_HM;
	matrix[12][13] *= mu_HM*2.*(5.-nu_H)/(5.-4.*nu_H);
	matrix[12][14] *= mu_HM;
	matrix[13][12] *= mu_HM;
	matrix[13][13] *= mu_HM*2.*(1.+nu_H)/(5.-4.*nu_H);
	matrix[13][14] *= mu_HM;
///////////////////////////////////////
//...решаем систему линейных уравнений;
	for (i = 0; i < dim_N; i++) {
		real_E f = 0.;
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
		real_E finv = 1./matrix[l0][l0]; matrix[l0][l0] = 1.;
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
double CCohes3D::TakeEshelby_shear_two(double c0, double eps, int max_iter)
{
	real_E mu_I = real_E(get_param(NUM_SHEAR+NUM_SHIFT)), nu_I = real_E(get_param(NUM_SHEAR+1+NUM_SHIFT)),
			 mu_M = real_E(get_param(NUM_SHEAR)), nu_M = real_E(get_param(NUM_SHEAR+1)),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M), ff = real_E(c0),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR2 = 1./pow(ff, 1./real_E(3.)), fR2 = RR2*RR2, fR5 = ff/fR2, 
			 mu_MI = mu_M/mu_I, AA = real_E(get_param(NUM_GEOMT))/mu_I, BB = real_E(get_param(NUM_GEOMT+1))/mu_I,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, sgn0, sgn1, nu_H, mu_HM, mu_HM0, mu_HM1, matrix[14][15],
			 matr[14][15] = {
					{ 1., -2.*QI1, 0., 0., -1., -1., -3., 2.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -1.,		0., 0., -1., -QM3, 2., 1., 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -6.*QI1, 0., 0., -1., 2., 12., 6.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -3.,		0., 0., -1., 2.*QM3, -8., 3., 0., 0., 0., 0., 0., 0., 0.},
					{ AA, -6.*QI1*AA,			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ BB, -3.*BB,				0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 2.*mu_I,  2.*mu_I*QI1, 0., 0., -2.*mu_M,  2.*mu_M*QM4,  24.*mu_M, -2.*mu_M*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 2.*mu_I, -2.*mu_I*QI2, 0., 0., -2.*mu_M, -2.*mu_M*QM5, -16.*mu_M,  2.*mu_M*QM2, 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 1.,      ff,   3.*fR5, -2.*QM1*fR2, 0., 0., 0., 0., -3., -1., 1.}, //...нормировку по радиусу в последних коэффициентах можно убрать!!!
					{ 0., 0., 0., 0., 1.,  QM3*ff,  -2.*fR5,        -fR2, 0., 0., 0., 0.,  2., -1., 1.},
					{ 0., 0., 0., 0., 1., -QM4*ff, -12.*fR5,     QM1*fR2, 0., 0., 0., 0., 12.,  1., 1.},
					{ 0., 0., 0., 0., 1.,  QM5*ff,   8.*fR5,    -QM2*fR2, 0., 0., 0., 0., -8., -1., 1.},
	};
	int k_iter = 0, k = 0, l;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...дописываем когезионную часть в матрице (A0, D0, A*0, D*0, A1, B1, C1, D1, A*1, B*1, C*1, D*1, C2, B2);
	real_E C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), kk_I = sqrt(C_I/ku_I), tt_I = exp(-2.*kk_I), 
			 C_M = get_param(NUM_SHEAR-1), kk_M = sqrt(C_M/ku_M), tt_D = exp(kk_M*(RR2-1.)), 
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), 
			 JHP1 =  (kk_M*RR2-1.)*ff*tt_D, JHP2 = ((sqr(kk_M*RR2)+3.)-3.*kk_M*RR2)*fR5*tt_D, 
			 JHM1 = -(kk_M*RR2+1.)*ff/tt_D, JHM2 = ((sqr(kk_M*RR2)+3.)+3.*kk_M*RR2)*fR5/tt_D,
			 kp_I = sqrt(C_I/mu_I), tp_I = exp(-2.*kp_I), 
			 kp_M = sqrt(C_M/mu_M), tp_D = exp(kp_M*(RR2-1.)),
			 HP1  = ((1.+tp_I)*kp_I-(1.-tp_I))*.5, 
			 HP2  = ((1.-tp_I)*(sqr(kp_I)+3.)-3.*(1.+tp_I)*kp_I)*.5,
			 JPP1 =  (kp_M-1.), JPP2 = ((sqr(kp_M)+3.)-3.*kp_M), 
			 JPM1 = -(kp_M+1.), JPM2 = ((sqr(kp_M)+3.)+3.*kp_M), 
			 JDP1 =  (kp_M*RR2-1.)*ff*tp_D, JDP2 = ((sqr(kp_M*RR2)+3.)-3.*kp_M*RR2)*fR5*tt_D, 
			 JDM1 = -(kp_M*RR2+1.)*ff/tp_D, JDM2 = ((sqr(kp_M*RR2)+3.)+3.*kp_M*RR2)*fR5/tt_D;
////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = -3.*HP2/C_I; 	//...сшивка функций;
	matr[0][3] =  3.*HH2/C_I-HH1/ku_I;
	matr[0][8] = 3.*JPP2/C_M; 
	matr[0][9] = 3.*JPM2/C_M; 
	matr[0][10] = JJP1/ku_M-3.*JJP2/C_M; 
	matr[0][11] = JJM1/ku_M-3.*JJM2/C_M;

	matr[1][2] =  2.*HP2/C_I-HP1/mu_I; 
	matr[1][3] = -2.*HH2/C_I;
	matr[1][8] = JPP1/mu_M-2.*JPP2/C_M; 
	matr[1][9] = JPM1/mu_M-2.*JPM2/C_M; 
	matr[1][10] = 2.*JJP2/C_M; 
	matr[1][11] = 2.*JJM2/C_M;

	//...сшивка нормальных производных;
	matr[2][2] = -3.*(HP1/mu_I-4.*HP2/C_I); 
	matr[2][3] = -(HH2-2.*HH1)/ku_I-12.*HH2/C_I;
	matr[2][8] = 3.*(JPP1/mu_M-4.*JPP2/C_M); 
	matr[2][9] = 3.*(JPM1/mu_M-4.*JPM2/C_M); 
	matr[2][10] = (JJP2-2.*JJP1)/ku_M+12.*JJP2/C_M; 
	matr[2][11] = (JJM2-2.*JJM1)/ku_M+12.*JJM2/C_M;

	matr[3][2] = -(HP2-HP1)/mu_I-8.*HP2/C_I; 
	matr[3][3] = -2.*(HH1/ku_I-4.*HH2/C_I);
	matr[3][8] = (JPP2-JPP1)/mu_M+8.*JPP2/C_M; 
	matr[3][9] = (JPM2-JPM1)/mu_M+8.*JPM2/C_M; 
	matr[3][10] = 2.*(JJP1/ku_M-4.*JJP2/C_M); 
	matr[3][11] = 2.*(JJM1/ku_M-4.*JJM2/C_M);

	//...сшивка моментов;
	matr[4][2] = -ku_I*3.*HP2/C_I-AA*3.*(HP1/mu_I-4.*HP2/C_I); 
	matr[4][3] = -ku_I*(HH1/ku_I-3.*HH2/C_I)-AA*((HH2-2.*HH1)/ku_I+12.*HH2/C_I);
	matr[4][8] = ku_M*3.*JPP2/C_M; 
	matr[4][9] = ku_M*3.*JPM2/C_M;
	matr[4][10] = ku_M*(JJP1/ku_M-3.*JJP2/C_M); 
	matr[4][11] = ku_M*(JJM1/ku_M-3.*JJM2/C_M);

	matr[5][2] = -mu_I*(HP1/mu_I-2.*HP2/C_I)-BB*((HP2-HP1)/mu_I+8.*HP2/C_I); 
	matr[5][3] = -mu_I*2.*HH2/C_I-BB*2.*(HH1/ku_I-4.*HH2/C_I);
	matr[5][8] = mu_M*(JPP1/mu_M-2.*JPP2/C_M); 
	matr[5][9] = mu_M*(JPM1/mu_M-2.*JPM2/C_M);
	matr[5][10] = mu_M*2.*JJP2/C_M; 
	matr[5][11] = mu_M*2.*JJM2/C_M;

	//...новый вариант: антисимметричный, симметрия по второму и третьему индексам;
	matr[6][2] = -9.*((ku_I-mu_I)*HP1/(6.*mu_I)-ku_I*HP2/C_I); 
	matr[6][3] =  ((mu_I+2.*ku_I)*HH1/ku_I-9.*ku_I*HH2/C_I);
	matr[6][8] =  9.*((ku_M-mu_M)*JPP1/(6.*mu_M)-ku_M*JPP2/C_M); 
	matr[6][9] =  9.*((ku_M-mu_M)*JPM1/(6.*mu_M)-ku_M*JPM2/C_M); 
	matr[6][10] = -(mu_M+2.*ku_M)*JJP1/ku_M+9.*ku_M*JJP2/C_M; 
	matr[6][11] = -(mu_M+2.*ku_M)*JJM1/ku_M+9.*ku_M*JJM2/C_M;

	matr[7][2] = (5.*mu_I+ku_I)*HP1/(2.*mu_I)-2.*(4.*mu_I-ku_I)*HP2/C_I; 
	matr[7][3] = (ku_I-mu_I)*HH1/ku_I+2.*(4.*mu_I-ku_I)*HH2/C_I;
	matr[7][8] = -(5.*mu_M-ku_M)*JPP1/(2.*mu_M)+2.*(4.*mu_M-ku_M)*JPP2/C_M; 
	matr[7][9] = -(5.*mu_M-ku_M)*JPM1/(2.*mu_M)+2.*(4.*mu_M-ku_M)*JPM2/C_M; 
	matr[7][10] = -((ku_M-mu_M)*JJP1/ku_M+2.*(4.*mu_M-ku_M)*JJP2/C_M); 
	matr[7][11] = -((ku_M-mu_M)*JJM1/ku_M+2.*(4.*mu_M-ku_M)*JJM2/C_M);

	//...прежний вариант: симметрия по первым двум индексам;
	//matr[6][2] = -3.*(HP1-(6.*mu_I+ku_I)*HP2/C_I); 
	//matr[6][3] =  ((4.*mu_I+ku_I)*HH1/ku_I-3.*(6.*mu_I+ku_I)*HH2/C_I);
	//matr[6][8] =  3.*(JPP1-(6.*mu_M+ku_M)*JPP2/C_M); 
	//matr[6][9] =  3.*(JPM1-(6.*mu_M+ku_M)*JPM2/C_M); 
	//matr[6][10] = -(4.*mu_M+ku_M)*JJP1/ku_M+3.*(6.*mu_M+ku_M)*JJP2/C_M; 
	//matr[6][11] = -(4.*mu_M+ku_M)*JJM1/ku_M+3.*(6.*mu_M+ku_M)*JJM2/C_M;

	//matr[7][2] = 4.*HP1-2.*(10.*mu_I-3.*ku_I)*HP2/C_I; 
	//matr[7][3] = 2.*(ku_I-2.*mu_I)*HH1/ku_I+2.*(10.*mu_I-3.*ku_I)*HH2/C_I;
	//matr[7][8] = -4.*JPP1+2.*(10.*mu_M-3.*ku_M)*JPP2/C_M; 
	//matr[7][9] = -4.*JPM1+2.*(10.*mu_M-3.*ku_M)*JPM2/C_M; 
	//matr[7][10] = -2.*((ku_M-2.*mu_M)*JJP1/ku_M+(10.*mu_M-3.*ku_M)*JJP2/C_M); 
	//matr[7][11] = -2.*((ku_M-2.*mu_M)*JJM1/ku_M+(10.*mu_M-3.*ku_M)*JJM2/C_M);

	//...равенство нулю когезионного поля;
	matr[8][8] = 3.*JDP2/C_M; 
	matr[8][9] = 3.*JDM2/C_M; 
	matr[8][10] = JHP1/ku_M-3.*JHP2/C_M; 
	matr[8][11] = JHM1/ku_M-3.*JHM2/C_M;

	matr[9][8] = JDP1/mu_M-2.*JDP2/C_M; 
	matr[9][9] = JDM1/mu_M-2.*JDM2/C_M; 
	matr[9][10] = 2.*JHP2/C_M; 
	matr[9][11] = 2.*JHM2/C_M;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_HM0 = 0.;
	nu_H = (1.5*KH-mu_HM0*mu_M)/(3.*KH+mu_HM0*mu_M);
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system_E(matrix, mu_HM0, nu_H);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.;

	mu_HM1 = 1.5*KH/mu_M;
	do { 
		mu_HM1 *= 10.; k_iter++; 
		nu_H = (1.5*KH-mu_HM1*mu_M)/(3.*KH+mu_HM1*mu_M);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
		optim = take_system_E(matrix, mu_HM1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_HM = (mu_HM0+mu_HM1)*.5; k_iter++;
		nu_H  = (1.5*KH-mu_HM*mu_M)/(3.*KH+mu_HM*mu_M);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_E(matrix, mu_HM, nu_H);
		if (optim*sgn0 > 0.) mu_HM0 = mu_HM; else
		if (optim*sgn0 < 0.) mu_HM1 = mu_HM; else mu_HM0 = mu_HM1 = mu_HM;
	}
	while(fabs(mu_HM1-mu_HM0) > eps && k_iter < max_iter);
	mu_HM = (mu_HM0+mu_HM1)*.5; 

	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6 ) mu_HM = -mu_HM;
	return to_double(mu_HM*mu_M);
}

/////////////////////////////////////////////////////////////////
//...трехфазная (симметричная) модель для сферического включения;
double CCohes3D::TakeEshelby_shear_sym(double ff, double eps, int max_iter)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 K1 = 2.*mu_I*(1.+nu_I)/(3.-6.*nu_I), K3 = 2.*mu_M*(1.+nu_M)/(3.-6.*nu_M),
			 KH = K3+ff*(K1-K3)/(1.+(1.-ff)*(K1-K3)/(K3+4./3.*mu_M)),
			 RR2 = 1./pow(ff, 1./3.), fR2 = RR2*RR2, fR5 = ff/fR2, 
			 mu_MI = mu_M/mu_I, AA = get_param(NUM_GEOMT)/mu_I, BB = get_param(NUM_GEOMT+1)/mu_I,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,
			 QI1 = 3.*nu_I/(7.-4.*nu_I),
			 QM1 = 3.*nu_M/(7.-4.*nu_M), 
			 QI2 = (7.+2.*nu_I)/(7.-4.*nu_I),
			 QM2 = (7.+2.*nu_M)/(7.-4.*nu_M), 
			 QM3 = 2.*(1.-2.*nu_M)/(5.-4.*nu_M), 
			 QM4 = 2.*(5.-nu_M)/(5.-4.*nu_M),
			 QM5 = 2.*(1.+nu_M)/(5.-4.*nu_M), optim, sgn0, sgn1, nu_H, mu_HM, mu_HM0, mu_HM1, matrix[14][15],
			 matr[14][15] = {
					{ 1., -2.*QI1, 0., 0., -1., -1., -3., 2.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -1.,		0., 0., -1., -QM3, 2., 1., 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -6.*QI1, 0., 0., -1., 2., 12., 6.*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 1., -3.,		0., 0., -1., 2.*QM3, -8., 3., 0., 0., 0., 0., 0., 0., 0.},
					{ AA, -6.*QI1*AA,			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ BB, -3.*BB,				0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 2.*mu_I,  2.*mu_I*QI1, 0., 0., -2.*mu_M,  2.*mu_M*QM4,  24.*mu_M, -2.*mu_M*QM1, 0., 0., 0., 0., 0., 0., 0.},
					{ 2.*mu_I, -2.*mu_I*QI2, 0., 0., -2.*mu_M, -2.*mu_M*QM5, -16.*mu_M,  2.*mu_M*QM2, 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
					{ 0., 0., 0., 0., 1.,      ff,   3.*fR5, -2.*QM1*fR2, 0., 0., 0., 0., -3., -1., 1.}, //...нормировку по радиусу в последних коэффициентах можно убрать!!!
					{ 0., 0., 0., 0., 1.,  QM3*ff,  -2.*fR5,        -fR2, 0., 0., 0., 0.,  2., -1., 1.},
					{ 0., 0., 0., 0., 1., -QM4*ff, -12.*fR5,     QM1*fR2, 0., 0., 0., 0., 12.,  1., 1.},
					{ 0., 0., 0., 0., 1.,  QM5*ff,   8.*fR5,    -QM2*fR2, 0., 0., 0., 0., -8., -1., 1.},
	};
	int k_iter = 0, k = 0, l;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...дописываем когезионную часть в матрице (A0, D0, A*0, D*0, A1, B1, C1, D1, A*1, B*1, C*1, D*1, C2, B2);
	double C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), kk_I = sqrt(C_I/ku_I), tt_I = exp(-2.*kk_I), 
			 C_M = get_param(NUM_SHEAR-1), kk_M = sqrt(C_M/ku_M), tt_D = exp(kk_M*(RR2-1.)), 
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, 
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), 
			 JHP1 =  (kk_M*RR2-1.)*ff*tt_D, JHP2 = ((sqr(kk_M*RR2)+3.)-3.*kk_M*RR2)*fR5*tt_D, 
			 JHM1 = -(kk_M*RR2+1.)*ff/tt_D, JHM2 = ((sqr(kk_M*RR2)+3.)+3.*kk_M*RR2)*fR5/tt_D,
			 kp_I = sqrt(C_I/mu_I), tp_I = exp(-2.*kp_I), 
			 kp_M = sqrt(C_M/mu_M), tp_D = exp(kp_M*(RR2-1.)),
			 HP1  = ((1.+tp_I)*kp_I-(1.-tp_I))*.5, 
			 HP2  = ((1.-tp_I)*(sqr(kp_I)+3.)-3.*(1.+tp_I)*kp_I)*.5,
			 JPP1 =  (kp_M-1.), JPP2 = ((sqr(kp_M)+3.)-3.*kp_M), 
			 JPM1 = -(kp_M+1.), JPM2 = ((sqr(kp_M)+3.)+3.*kp_M), 
			 JDP1 =  (kp_M*RR2-1.)*ff*tp_D, JDP2 = ((sqr(kp_M*RR2)+3.)-3.*kp_M*RR2)*fR5*tt_D, 
			 JDM1 = -(kp_M*RR2+1.)*ff/tp_D, JDM2 = ((sqr(kp_M*RR2)+3.)+3.*kp_M*RR2)*fR5/tt_D;
////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = -3.*HP2/C_I; 	//...сшивка функций;
	matr[0][3] =  3.*HH2/C_I-HH1/ku_I;
	matr[0][8] = 3.*JPP2/C_M; 
	matr[0][9] = 3.*JPM2/C_M; 
	matr[0][10] = JJP1/ku_M-3.*JJP2/C_M; 
	matr[0][11] = JJM1/ku_M-3.*JJM2/C_M;

	matr[1][2] =  2.*HP2/C_I-HP1/mu_I; 
	matr[1][3] = -2.*HH2/C_I;
	matr[1][8] = JPP1/mu_M-2.*JPP2/C_M; 
	matr[1][9] = JPM1/mu_M-2.*JPM2/C_M; 
	matr[1][10] = 2.*JJP2/C_M; 
	matr[1][11] = 2.*JJM2/C_M;

	//...сшивка нормальных производных;
	matr[2][2] = -3.*(HP1/mu_I-4.*HP2/C_I); 
	matr[2][3] = -(HH2-2.*HH1)/ku_I-12.*HH2/C_I;
	matr[2][8] = 3.*(JPP1/mu_M-4.*JPP2/C_M); 
	matr[2][9] = 3.*(JPM1/mu_M-4.*JPM2/C_M); 
	matr[2][10] = (JJP2-2.*JJP1)/ku_M+12.*JJP2/C_M; 
	matr[2][11] = (JJM2-2.*JJM1)/ku_M+12.*JJM2/C_M;

	matr[3][2] = -(HP2-HP1)/mu_I-8.*HP2/C_I; 
	matr[3][3] = -2.*(HH1/ku_I-4.*HH2/C_I);
	matr[3][8] = (JPP2-JPP1)/mu_M+8.*JPP2/C_M; 
	matr[3][9] = (JPM2-JPM1)/mu_M+8.*JPM2/C_M; 
	matr[3][10] = 2.*(JJP1/ku_M-4.*JJP2/C_M); 
	matr[3][11] = 2.*(JJM1/ku_M-4.*JJM2/C_M);

	//...сшивка моментов;
	matr[4][2] = -ku_I*3.*HP2/C_I-AA*3.*(HP1/mu_I-4.*HP2/C_I); 
	matr[4][3] = -ku_I*(HH1/ku_I-3.*HH2/C_I)-AA*((HH2-2.*HH1)/ku_I+12.*HH2/C_I);
	matr[4][8] = ku_M*3.*JPP2/C_M; 
	matr[4][9] = ku_M*3.*JPM2/C_M;
	matr[4][10] = ku_M*(JJP1/ku_M-3.*JJP2/C_M); 
	matr[4][11] = ku_M*(JJM1/ku_M-3.*JJM2/C_M);

	matr[5][2] = -mu_I*(HP1/mu_I-2.*HP2/C_I)-BB*((HP2-HP1)/mu_I+8.*HP2/C_I); 
	matr[5][3] = -mu_I*2.*HH2/C_I-BB*2.*(HH1/ku_I-4.*HH2/C_I);
	matr[5][8] = mu_M*(JPP1/mu_M-2.*JPP2/C_M); 
	matr[5][9] = mu_M*(JPM1/mu_M-2.*JPM2/C_M);
	matr[5][10] = mu_M*2.*JJP2/C_M; 
	matr[5][11] = mu_M*2.*JJM2/C_M;

	//...новый вариант: антисимметричный, симметрия по второму и третьему индексам;
	//matr[6][2] = -9.*((ku_I-mu_I)*HP1/(6.*mu_I)-ku_I*HP2/C_I); 
	//matr[6][3] =  ((mu_I+2.*ku_I)*HH1/ku_I-9.*ku_I*HH2/C_I);
	//matr[6][8] =  9.*((ku_M-mu_M)*JPP1/(6.*mu_M)-ku_M*JPP2/C_M); 
	//matr[6][9] =  9.*((ku_M-mu_M)*JPM1/(6.*mu_M)-ku_M*JPM2/C_M); 
	//matr[6][10] = -(mu_M+2.*ku_M)*JJP1/ku_M+9.*ku_M*JJP2/C_M; 
	//matr[6][11] = -(mu_M+2.*ku_M)*JJM1/ku_M+9.*ku_M*JJM2/C_M;

	//matr[7][2] = (5.*mu_I+ku_I)*HP1/(2.*mu_I)-2.*(4.*mu_I-ku_I)*HP2/C_I; 
	//matr[7][3] = (ku_I-mu_I)*HH1/ku_I+2.*(4.*mu_I-ku_I)*HH2/C_I;
	//matr[7][8] = -(5.*mu_M-ku_M)*JPP1/(2.*mu_M)+2.*(4.*mu_M-ku_M)*JPP2/C_M; 
	//matr[7][9] = -(5.*mu_M-ku_M)*JPM1/(2.*mu_M)+2.*(4.*mu_M-ku_M)*JPM2/C_M; 
	//matr[7][10] = -((ku_M-mu_M)*JJP1/ku_M+2.*(4.*mu_M-ku_M)*JJP2/C_M); 
	//matr[7][11] = -((ku_M-mu_M)*JJM1/ku_M+2.*(4.*mu_M-ku_M)*JJM2/C_M);

	//...прежний вариант: симметрия по первым двум индексам;
	matr[6][2] = -3.*(HP1-(6.*mu_I+ku_I)*HP2/C_I); 
	matr[6][3] =  ((4.*mu_I+ku_I)*HH1/ku_I-3.*(6.*mu_I+ku_I)*HH2/C_I);
	matr[6][8] =  3.*(JPP1-(6.*mu_M+ku_M)*JPP2/C_M); 
	matr[6][9] =  3.*(JPM1-(6.*mu_M+ku_M)*JPM2/C_M); 
	matr[6][10] = -(4.*mu_M+ku_M)*JJP1/ku_M+3.*(6.*mu_M+ku_M)*JJP2/C_M; 
	matr[6][11] = -(4.*mu_M+ku_M)*JJM1/ku_M+3.*(6.*mu_M+ku_M)*JJM2/C_M;

	matr[7][2] = 4.*HP1-2.*(10.*mu_I-3.*ku_I)*HP2/C_I; 
	matr[7][3] = 2.*(ku_I-2.*mu_I)*HH1/ku_I+2.*(10.*mu_I-3.*ku_I)*HH2/C_I;
	matr[7][8] = -4.*JPP1+2.*(10.*mu_M-3.*ku_M)*JPP2/C_M; 
	matr[7][9] = -4.*JPM1+2.*(10.*mu_M-3.*ku_M)*JPM2/C_M; 
	matr[7][10] = -2.*((ku_M-2.*mu_M)*JJP1/ku_M+(10.*mu_M-3.*ku_M)*JJP2/C_M); 
	matr[7][11] = -2.*((ku_M-2.*mu_M)*JJM1/ku_M+(10.*mu_M-3.*ku_M)*JJM2/C_M);

	//...новый вариант: антисимметричный, нет симметрии по второму и третьему индексам;
	//matr[6][2] = -3.*(HP1*.5-(4.*mu_I+ku_I)*HP2/C_I); 
	//matr[6][3] =  ((3.*mu_I+ku_I)*HH1/ku_I-3.*(4.*mu_I+ku_I)*HH2/C_I);
	//matr[6][8] =  3.*(JPP1*.5-(4.*mu_M+ku_M)*JPP2/C_M); 
	//matr[6][9] =  3.*(JPM1*.5-(4.*mu_M+ku_M)*JPM2/C_M); 
	//matr[6][10] = -(3.*mu_M+ku_M)*JJP1/ku_M+3.*(4.*mu_M+ku_M)*JJP2/C_M; 
	//matr[6][11] = -(3.*mu_M+ku_M)*JJM1/ku_M+3.*(4.*mu_M+ku_M)*JJM2/C_M;

	//matr[7][2] = 3.5*HP1-2.*(8.*mu_I-3.*ku_I)*HP2/C_I; 
	//matr[7][3] = (2.*ku_I-3.*mu_I)*HH1/ku_I+2.*(8.*mu_I-3.*ku_I)*HH2/C_I;
	//matr[7][8] = -3.5*JPP1+2.*(8.*mu_M-3.*ku_M)*JPP2/C_M; 
	//matr[7][9] = -3.5*JPM1+2.*(8.*mu_M-3.*ku_M)*JPM2/C_M; 
	//matr[7][10] = -((2.*ku_M-3.*mu_M)*JJP1/ku_M+2.*(8.*mu_M-3.*ku_M)*JJP2/C_M); 
	//matr[7][11] = -((2.*ku_M-3.*mu_M)*JJM1/ku_M+2.*(8.*mu_M-3.*ku_M)*JJM2/C_M);

	//...равенство нулю когезионного поля;
	matr[8][8] = 3.*JDP2/C_M; 
	matr[8][9] = 3.*JDM2/C_M; 
	matr[8][10] = JHP1/ku_M-3.*JHP2/C_M; 
	matr[8][11] = JHM1/ku_M-3.*JHM2/C_M;

	matr[9][8] = JDP1/mu_M-2.*JDP2/C_M; 
	matr[9][9] = JDM1/mu_M-2.*JDM2/C_M; 
	matr[9][10] = 2.*JHP2/C_M; 
	matr[9][11] = 2.*JHM2/C_M;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_HM0 = 0.;
	nu_H = (1.5*KH-mu_HM0*mu_M)/(3.*KH+mu_HM0*mu_M);
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system(matrix, mu_HM0, nu_H);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.;

	mu_HM1 = 1.5*KH/mu_M;
	do { 
		mu_HM1 *= 10.; k_iter++; 
		nu_H = (1.5*KH-mu_HM1*mu_M)/(3.*KH+mu_HM1*mu_M);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_HM1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_HM = (mu_HM0+mu_HM1)*.5; k_iter++;
		nu_H  = (1.5*KH-mu_HM*mu_M)/(3.*KH+mu_HM*mu_M);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system(matrix, mu_HM, nu_H);
		if (optim*sgn0 > 0.) mu_HM0 = mu_HM; else
		if (optim*sgn0 < 0.) mu_HM1 = mu_HM; else mu_HM0 = mu_HM1 = mu_HM;
	}
	while(fabs(mu_HM1-mu_HM0) > eps && k_iter < max_iter);
	mu_HM = (mu_HM0+mu_HM1)*.5; 

	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6 ) mu_HM = -mu_HM;
	return(mu_HM*mu_M);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для трехфазной модели (модуль сдвига);
double take_system_shear(double matrix[14][15], double mu_MH, double nu_H)
{
	int dim_N = 14, ii[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}, i, k, l, k0, l0;
//////////////////////////////////////////
//...заполняем систему линейных уравнений;
	matrix[10][12] *= mu_MH;
	matrix[10][13] *= mu_MH*(1.25-nu_H);
	matrix[10][14] *= mu_MH;
	matrix[11][12] *= mu_MH;
	matrix[11][13] *= mu_MH*(0.5-nu_H);
	matrix[11][14] *= mu_MH;
	matrix[12][13] *= (5.-nu_H)*.5;
	matrix[13][13] *= (1.+nu_H)*.5;
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
double CCohes3D::TakeEshelby_shear(double ff, double nju1, double nju2, double E1, double E2, double l1, double l2)
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
			 QM9 = 1.-2.*nu_M, optim, sgn0, sgn1, nu_H, mu_H, mu_MH, mu_MH0, mu_MH1, matrix[14][15], eps = EE, max_iter = 100, eps0 = 1e-6,
			 matr[14][15] = {
					{ mu_M, -mu_M*QI3*fR2,    0., 0., -mu_I, -mu_I*QM8, -mu_I*2.25, mu_I*QM3*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI1*fR2,    0., 0., -mu_I, -mu_I*QM7,  mu_I*1.5,  mu_I*QM1*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI3*fR2*3., 0., 0., -mu_I,  mu_I*QM6,  mu_I*9.,   mu_I*QM3*fR2*3., 0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI1*fR2*3., 0., 0., -mu_I,  mu_I*QM9, -mu_I*6.,   mu_I*QM1*fR2*3., 0., 0., 0., 0.,    0., 0., 0.},
					{   0.,               0., 0., 0.,    0.,        0.,       0.,                0., 0., 0., 0., 0.,    0., 0., 0.},
					{   0.,               0., 0., 0.,    0.,        0.,       0.,                0., 0., 0., 0., 0.,    0., 0., 0.},
					{   1.,       QI3*fR2*.5, 0., 0.,   -1.,       QM4,       9.,       -QM3*fR2*.5, 0., 0., 0., 0.,    0., 0., 0.},
					{   1.,      -QI2*fR2,    0., 0.,   -1.,      -QM5,      -6.,        QM2*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    0.,        0.,       0.,             0.,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    0.,        0.,       0.,             0.,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    1.,    QM8*c0,  2.25*s0,           -QM3,    0., 0., 0., 0., -2.25, -1., 1.},
					{   0.,            0.,    0., 0.,    1.,    QM7*c0,  -1.5*s0,           -QM1,    0., 0., 0., 0.,   1.5, -1., 1.},
					{   0.,            0.,    0., 0.,    1.,   -QM4*c0,   -9.*s0,            QM3*.5, 0., 0., 0., 0.,    9.,  1., 1.},
					{   0.,            0.,    0., 0.,    1.,    QM5*c0,    6.*s0,           -QM2,    0., 0., 0., 0.,   -6., -1., 1.},
	};
	int k_iter = 0, k, l;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...дописываем когезионную часть в матрице (A0, D0, A*0, C*0, A1, B1, C1, D1, A*1, B*1, C*1, D*1, C2, B2);
	//double l_I = sqr(l1), g_I = (.5-nu_I)/(1.-nu_I), tt_I = exp(-2./l1), tp_I = exp(-2./l1*sqrt(g_I)), 
	//		 l_M = sqr(l2), g_M = (.5-nu_M)/(1.-nu_M), tt_M = exp(-2./l2), tp_M = exp(-2./l2*sqrt(g_M)), 
	//		 tt_D = exp(( 1./RR1-1.)/l2), tp_D = exp(( 1./RR1-1.)/l2*sqrt(g_M)), 
	//		 tt_T = exp((-1./RR1-1.)/l2), tp_T = exp((-1./RR1-1.)/l2*sqrt(g_M)),

	//		 HH1  = ((1.+tt_I)/l1-(1.-tt_I))*.5/c0, 
	//		 HP1  = ((1.+tp_I)/l1*sqrt(g_I)-(1.-tp_I))*.5/c0, 
	//		 HH2  = ((1.-tt_I)*(1./l_I+3.)-3.*(1.+tt_I)/l1)*.5/s0,
	//		 HP2  = ((1.-tp_I)*(1./l_I*g_I+3.)-3.*(1.+tp_I)/l1*sqrt(g_I))*.5/s0, 
	//		 
	//		 JJH1 = ((1.+tt_M)/l2-(1.-tt_M))*.5/c0, 
	//		 JJP1 = ((1.+tp_M)/l2*sqrt(g_M)-(1.-tp_M))*.5/c0, 
	//		 JJH2 = ((1.-tt_M)*(1./l_M+3.)-3.*(1.+tt_M)/l2)*.5/s0,
	//		 JJP2 = ((1.-tp_M)*(1./l_M*g_M+3.)-3.*(1.+tp_M)/l2*sqrt(g_M))*.5/s0,
	//		 
	//		 JJM1 = ((1.-tt_M)/l2-(1.+tt_M))*.5/c0, 
	//		 JJD1 = ((1.-tp_M)/l2*sqrt(g_M)-(1.+tp_M))*.5/c0, 
	//		 JJM2 = ((1.+tt_M)*(1./l_M+3.)-3.*(1.-tt_M)/l2)*.5/s0,
	//		 JJD2 = ((1.+tp_M)*(1./l_M*g_M+3.)-3.*(1.-tp_M)/l2*sqrt(g_M))*.5/s0,

	//		 JHP1 = ((tt_D+tt_T)/(l2*RR1)-(tt_D-tt_T))*.5,
	//		 JDP1 = ((tp_D+tp_T)/(l2*RR1)*sqrt(g_M)-(tp_D-tp_T))*.5,
	//		 JHP2 = ((tt_D-tt_T)*(1./(l_M*fR2)+3.)-3.*(tt_D+tt_T)/(l2*RR1))*.5,
	//		 JDP2 = ((tp_D-tp_T)*(1./(l_M*fR2)*g_M+3.)-3.*(tp_D+tp_T)/(l2*RR1)*sqrt(g_M))*.5,

	//		 JHM1 = ((tt_D-tt_T)/(l2*RR1)-(tt_D+tt_T))*.5,
	//		 JDM1 = ((tp_D-tp_T)/(l2*RR1)*sqrt(g_M)-(tp_D+tp_T))*.5,
	//		 JHM2 = ((tt_D+tt_T)*(1./(l_M*fR2)+3.)-3.*(tt_D-tt_T)/(l2*RR1))*.5,
	//		 JDM2 = ((tp_D+tp_T)*(1./(l_M*fR2)*g_M+3.)-3.*(tp_D-tp_T)/(l2*RR1)*sqrt(g_M))*.5;

///////////////////////////////////////////////////////////
//...вариант упрощенных функций в слое (просто экспоненты);
	double l_I = sqr(l1), g_I = (.5-nu_I)/(1.-nu_I), tt_I = exp(-2./l1), tp_I = exp(-2./l1*sqrt(g_I)), 
			 l_M = sqr(l2), g_M = (.5-nu_M)/(1.-nu_M), tt_M = exp(-2./l2), tp_M = exp(-2./l2*sqrt(g_M)), 
			 tt_D = exp(( 1./RR1-1.)/l2), tp_D = exp(( 1./RR1-1.)/l2*sqrt(g_M)), 
			 tt_T = exp((-1./RR1-1.)/l2), tp_T = exp((-1./RR1-1.)/l2*sqrt(g_M)),

			 HH1  = ((1.+tt_I)/l1-(1.-tt_I))*.5/c0, 
			 HP1  = ((1.+tp_I)/l1*sqrt(g_I)-(1.-tp_I))*.5/c0, 
			 HH2  = ((1.-tt_I)*(1./l_I+3.)-3.*(1.+tt_I)/l1)*.5/s0,
			 HP2  = ((1.-tp_I)*(1./l_I*g_I+3.)-3.*(1.+tp_I)/l1*sqrt(g_I))*.5/s0, 
			 
			 JJH1 = (1./l2-1.)/c0, 
			 JJP1 = (1./l2*sqrt(g_M)-1.)/c0, 
			 JJH2 = ((1./l_M+3.)-3./l2)/s0,
			 JJP2 = ((1./l_M*g_M+3.)-3./l2*sqrt(g_M))*.5/s0,
			 
			 JJM1 = -(1./l2+1.)*tt_M/c0, 
			 JJD1 = -(1./l2*sqrt(g_M)+1.)*tp_M/c0, 
			 JJM2 = ((1./l_M+3.)+3./l2)*tt_M/s0,
			 JJD2 = ((1./l_M*g_M+3.)+3./l2*sqrt(g_M))*tp_M/s0,

			 JHP1 = (1./(l2*RR1)-1.)*tt_D,
			 JDP1 = (1./(l2*RR1)*sqrt(g_M)-1.)*tp_D,
			 JHP2 = ((1./(l_M*fR2)+3.)-3./(l2*RR1))*tt_D,
			 JDP2 = ((1./(l_M*fR2)*g_M+3.)-3./(l2*RR1)*sqrt(g_M))*tp_D,

			 JHM1 = -(1./(l2*RR1)+1.)*tt_T,
			 JDM1 = -(1./(l2*RR1)*sqrt(g_M)+1.)*tp_T,
			 JHM2 = ((1./(l_M*fR2)+3.)+3./(l2*RR1))*tt_T,
			 JDM2 = ((1./(l_M*fR2)*g_M+3.)+3./(l2*RR1)*sqrt(g_M))*tp_T;

////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = mu_M*(g_I*HP1-3.*l_I*HP2); 	//...сшивка функций;
	matr[0][3] = mu_M*3.*l_I*HH2;
	matr[0][8] = -mu_I*(g_M*JJP1-3.*l_M*JJP2); 
	matr[0][9] = -mu_I*(g_M*JJD1-3.*l_M*JJD2); 
	matr[0][10] = -mu_I*3.*l_M*JJH2; 
	matr[0][11] = -mu_I*3.*l_M*JJM2;

	matr[1][2] = mu_M*2.*l_I*HP2; 	
	matr[1][3] = mu_M*(HH1-2.*l_I*HH2);
	matr[1][8] = -mu_I*2.*l_M*JJP2;  
	matr[1][9] = -mu_I*2.*l_M*JJD2;  
	matr[1][10] = -mu_I*(JJH1-2.*l_M*JJH2); 
	matr[1][11] = -mu_I*(JJM1-2.*l_M*JJM2);

	//...сшивка нормальных производных;
	matr[2][2] = mu_M*(g_I*(fR2*HP2-HP1)+12.*l_I*HP2); 
	matr[2][3] = mu_M*3.*(HH1-4.*l_I*HH2);
	matr[2][8] = -mu_I*(g_M*(fR2*JJP2-JJP1)+12.*l_M*JJP2); 
	matr[2][9] = -mu_I*(g_M*(fR2*JJD2-JJD1)+12.*l_M*JJD2); 
	matr[2][10] = -mu_I*3.*(JJH1-4.*l_M*JJH2); 
	matr[2][11] = -mu_I*3.*(JJM1-4.*l_M*JJM2);

	matr[3][2] = mu_M*2.*(g_I*HP1-4.*l_I*HP2); 
	matr[3][3] = mu_M*((fR2*HH2-HH1)+8.*l_I*HH2);
	matr[3][8] = -mu_I*2.*(g_M*JJP1-4.*l_M*JJP2); 
	matr[3][9] = -mu_I*2.*(g_M*JJD1-4.*l_M*JJD2); 
	matr[3][10] = -mu_I*((fR2*JJH2-JJH1)+8.*l_M*JJH2); 
	matr[3][11] = -mu_I*((fR2*JJM2-JJM1)+8.*l_M*JJM2);

	//...сшивка моментов;
	matr[4][2] = g_M*(g_I*HP1-3.*l_I*HP2); 
	matr[4][3] = g_M*3.*l_I*HH2;
	matr[4][8] = -g_I*(g_M*JJP1-3.*l_M*JJP2); 
	matr[4][9] = -g_I*(g_M*JJD1-3.*l_M*JJD2);
	matr[4][10] = -g_I*3.*l_M*JJH2; 
	matr[4][11] = -g_I*3.*l_M*JJM2;

	matr[5][2] = 2.*l_I*HP2; 	 
	matr[5][3] = HH1-2.*l_I*HH2;
	matr[5][8] = -2.*l_M*JJP2; 
	matr[5][9] = -2.*l_M*JJD2;
	matr[5][10] = -JJH1+2.*l_M*JJH2; 
	matr[5][11] = -JJM1+2.*l_M*JJM2;

	//...антисимметричный вариант (симметрия по второму и третьему индексам);
	matr[6][2] = -((2.+g_I)*HP1-9.*l_I/g_I*HP2)*.5; 
	matr[6][3] =  (1./g_I-1.)*.75*HH1-4.5*l_I/g_I*HH2;
	matr[6][8] =  ((2.+g_M)*JJP1-9.*l_M/g_M*JJP2)*.5; 
	matr[6][9] =  ((2.+g_M)*JJD1-9.*l_M/g_M*JJD2)*.5; 
	matr[6][10] = -(1./g_M-1.)*.75*JJH1+4.5*l_M/g_M*JJH2; 
	matr[6][11] = -(1./g_M-1.)*.75*JJM1+4.5*l_M/g_M*JJM2;

	matr[7][2] = -(1.-g_I)*.5*HP1-(4.-1./g_I)*l_I*HP2; 
	matr[7][3] = -(5.+1./g_I)*.25*HH1+(4.-1./g_I)*l_I*HH2;
	matr[7][8] =  (1.-g_M)*.5*JJP1+(4.-1./g_M)*l_M*JJP2; 
	matr[7][9] =  (1.-g_M)*.5*JJD1+(4.-1./g_M)*l_M*JJD2; 
	matr[7][10] = (5.+1./g_M)*.25*JJH1-(4.-1./g_M)*l_M*JJH2; 
	matr[7][11] = (5.+1./g_M)*.25*JJM1-(4.-1./g_M)*l_M*JJM2;

	//...прежний симметричный, но нефизичный вариант (симметрия по первому и второму индексам);
	//matr[6][2] = -((1.+4.*g_I)*HP1-3.*(6.+1./g_I)*l_I*HP2)*.5; 
	//matr[6][3] =  (HH1-(6.+1./g_I)*l_I*HH2)*1.5;
	//matr[6][8] =  ((1.+4.*g_M)*JJP1-3.*(6.+1./g_M)*l_M*JJP2)*.5; 
	//matr[6][9] =  ((1.+4.*g_M)*JJD1-3.*(6.+1./g_M)*l_M*JJD2)*.5; 
	//matr[6][10] = -(JJH1-(6.+1./g_M)*l_M*JJH2)*1.5; 
	//matr[6][11] = -(JJM1-(6.+1./g_M)*l_M*JJM2)*1.5;

	//matr[7][2] = -((1.-2.*g_I)*HP1+(10.-3./g_I)*l_I*HP2); 
	//matr[7][3] = -(2.*HH1-(10.-3./g_I)*l_I*HH2);
	//matr[7][8] =  (1.-2.*g_M)*JJP1+(10.-3./g_M)*l_M*JJP2; 
	//matr[7][9] =  (1.-2.*g_M)*JJD1+(10.-3./g_M)*l_M*JJD2; 
	//matr[7][10] = 2.*JJH1-(10.-3./g_M)*l_M*JJH2; 
	//matr[7][11] = 2.*JJM1-(10.-3./g_M)*l_M*JJM2;

	//...равенство нулю когезионного поля;
	matr[8][8] = g_M*JDP1-3.*l_M*JDP2; 
	matr[8][9] = g_M*JDM1-3.*l_M*JDM2; 
	matr[8][10] = 3.*l_M*JHP2; 
	matr[8][11] = 3.*l_M*JHM2;

	matr[9][8] = 2.*l_M*JDP2; 
	matr[9][9] = 2.*l_M*JDM2; 
	matr[9][10] = JHP1-2.*l_M*JHP2; 
	matr[9][11] = JHM1-2.*l_M*JHM2;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.; nu_H = -1.;
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system_shear(matrix, mu_MH0, nu_H);
	if (optim < 0.) sgn0 = -1.; else sgn0 = 1.;

	mu_MH1 = 2.*mu_M/(3.*KH);
	do { 
		mu_MH1 *= 10.; k_iter++; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_shear(matrix, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H = (1.5*KH*mu_MH/mu_M-1.)/(3.*KH*mu_MH/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_shear(matrix, mu_MH, nu_H);
		if (optim*sgn0 > 0.) mu_MH0 = mu_MH; else
		if (optim*sgn0 < 0.) mu_MH1 = mu_MH; else mu_MH0 = mu_MH1 = mu_MH;
	}
	while(fabs(mu_MH1-mu_MH0) > eps && k_iter < max_iter);

	mu_H = 2.*mu_M/(mu_MH0+mu_MH1);
	if (sgn0*sgn1 > 0. || fabs(optim) > eps0 ) mu_H = -mu_H;
	return(mu_H);
}

/////////////////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (старый итерационный алгоритм);
double CCohes3D::TakeEshelby_shear_old(double ff, double nju1, double nju2, double E1, double E2, double l1, double l2)
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
			 QM9 = 1.-2.*nu_M, optim, sgn0, sgn1, nu_H, mu_H, mu_MH, mu_MH0, mu_MH1, matrix[14][15], eps = EE, max_iter = 100, eps0 = 1e-6,
			 matr[14][15] = {
					{ mu_M, -mu_M*QI3*fR2,    0., 0., -mu_I, -mu_I*QM8, -mu_I*2.25, mu_I*QM3*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI1*fR2,    0., 0., -mu_I, -mu_I*QM7,  mu_I*1.5,  mu_I*QM1*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI3*fR2*3., 0., 0., -mu_I,  mu_I*QM6,  mu_I*9.,   mu_I*QM3*fR2*3., 0., 0., 0., 0.,    0., 0., 0.},
					{ mu_M, -mu_M*QI1*fR2*3., 0., 0., -mu_I,  mu_I*QM9, -mu_I*6.,   mu_I*QM1*fR2*3., 0., 0., 0., 0.,    0., 0., 0.},
					{   0.,               0., 0., 0.,    0.,        0.,       0.,                0., 0., 0., 0., 0.,    0., 0., 0.},
					{   0.,               0., 0., 0.,    0.,        0.,       0.,                0., 0., 0., 0., 0.,    0., 0., 0.},
					{   1.,       QI3*fR2*.5, 0., 0.,   -1.,       QM4,       9.,       -QM3*fR2*.5, 0., 0., 0., 0.,    0., 0., 0.},
					{   1.,      -QI2*fR2,    0., 0.,   -1.,      -QM5,      -6.,        QM2*fR2,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    0.,        0.,       0.,             0.,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    0.,        0.,       0.,             0.,    0., 0., 0., 0.,    0., 0., 0.},
					{   0.,            0.,    0., 0.,    1.,    QM8*c0,  2.25*s0,           -QM3,    0., 0., 0., 0., -2.25, -1., 1.},
					{   0.,            0.,    0., 0.,    1.,    QM7*c0,  -1.5*s0,           -QM1,    0., 0., 0., 0.,   1.5, -1., 1.},
					{   0.,            0.,    0., 0.,    1.,   -QM4*c0,   -9.*s0,            QM3*.5, 0., 0., 0., 0.,    9.,  1., 1.},
					{   0.,            0.,    0., 0.,    1.,    QM5*c0,    6.*s0,           -QM2,    0., 0., 0., 0.,   -6., -1., 1.},
	};
	int k_iter = 0, k, l;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...дописываем когезионную часть в матрице (A0, D0, A*0, C*0, A1, B1, C1, D1, A*1, B*1, C*1, D*1, C2, B2);
	//double l_I = sqr(l1), g_I = (.5-nu_I)/(1.-nu_I), tt_I = exp(-2./l1), tp_I = exp(-2./l1*sqrt(g_I)), 
	//		 l_M = sqr(l2), g_M = (.5-nu_M)/(1.-nu_M), tt_M = exp(-2./l2), tp_M = exp(-2./l2*sqrt(g_M)), 
	//		 tt_D = exp(( 1./RR1-1.)/l2), tp_D = exp(( 1./RR1-1.)/l2*sqrt(g_M)), 
	//		 tt_T = exp((-1./RR1-1.)/l2), tp_T = exp((-1./RR1-1.)/l2*sqrt(g_M)),

	//		 HH1  = ((1.+tt_I)/l1-(1.-tt_I))*.5/c0, 
	//		 HP1  = ((1.+tp_I)/l1*sqrt(g_I)-(1.-tp_I))*.5/c0, 
	//		 HH2  = ((1.-tt_I)*(1./l_I+3.)-3.*(1.+tt_I)/l1)*.5/s0,
	//		 HP2  = ((1.-tp_I)*(1./l_I*g_I+3.)-3.*(1.+tp_I)/l1*sqrt(g_I))*.5/s0, 
	//		 
	//		 JJH1 = ((1.+tt_M)/l2-(1.-tt_M))*.5/c0, 
	//		 JJP1 = ((1.+tp_M)/l2*sqrt(g_M)-(1.-tp_M))*.5/c0, 
	//		 JJH2 = ((1.-tt_M)*(1./l_M+3.)-3.*(1.+tt_M)/l2)*.5/s0,
	//		 JJP2 = ((1.-tp_M)*(1./l_M*g_M+3.)-3.*(1.+tp_M)/l2*sqrt(g_M))*.5/s0,
	//		 
	//		 JJM1 = ((1.-tt_M)/l2-(1.+tt_M))*.5/c0, 
	//		 JJD1 = ((1.-tp_M)/l2*sqrt(g_M)-(1.+tp_M))*.5/c0, 
	//		 JJM2 = ((1.+tt_M)*(1./l_M+3.)-3.*(1.-tt_M)/l2)*.5/s0,
	//		 JJD2 = ((1.+tp_M)*(1./l_M*g_M+3.)-3.*(1.-tp_M)/l2*sqrt(g_M))*.5/s0,

	//		 JHP1 = ((tt_D+tt_T)/(l2*RR1)-(tt_D-tt_T))*.5,
	//		 JDP1 = ((tp_D+tp_T)/(l2*RR1)*sqrt(g_M)-(tp_D-tp_T))*.5,
	//		 JHP2 = ((tt_D-tt_T)*(1./(l_M*fR2)+3.)-3.*(tt_D+tt_T)/(l2*RR1))*.5,
	//		 JDP2 = ((tp_D-tp_T)*(1./(l_M*fR2)*g_M+3.)-3.*(tp_D+tp_T)/(l2*RR1)*sqrt(g_M))*.5,

	//		 JHM1 = ((tt_D-tt_T)/(l2*RR1)-(tt_D+tt_T))*.5,
	//		 JDM1 = ((tp_D-tp_T)/(l2*RR1)*sqrt(g_M)-(tp_D+tp_T))*.5,
	//		 JHM2 = ((tt_D+tt_T)*(1./(l_M*fR2)+3.)-3.*(tt_D-tt_T)/(l2*RR1))*.5,
	//		 JDM2 = ((tp_D+tp_T)*(1./(l_M*fR2)*g_M+3.)-3.*(tp_D-tp_T)/(l2*RR1)*sqrt(g_M))*.5;

///////////////////////////////////////////////////////////
//...вариант упрощенных функций в слое (просто экспоненты);
	double l_I = sqr(l1), g_I = (.5-nu_I)/(1.-nu_I), tt_I = exp(-2./l1), tp_I = exp(-2./l1*sqrt(g_I)), 
			 l_M = sqr(l2), g_M = (.5-nu_M)/(1.-nu_M), tt_M = exp(-2./l2), tp_M = exp(-2./l2*sqrt(g_M)), 
			 tt_D = exp(( 1./RR1-1.)/l2), tp_D = exp(( 1./RR1-1.)/l2*sqrt(g_M)), 
			 tt_T = exp((-1./RR1-1.)/l2), tp_T = exp((-1./RR1-1.)/l2*sqrt(g_M)),

			 HH1  = ((1.+tt_I)/l1-(1.-tt_I))*.5/c0, 
			 HP1  = ((1.+tp_I)/l1*sqrt(g_I)-(1.-tp_I))*.5/c0, 
			 HH2  = ((1.-tt_I)*(1./l_I+3.)-3.*(1.+tt_I)/l1)*.5/s0,
			 HP2  = ((1.-tp_I)*(1./l_I*g_I+3.)-3.*(1.+tp_I)/l1*sqrt(g_I))*.5/s0, 
			 
			 JJH1 = (1./l2-1.)/c0, 
			 JJP1 = (1./l2*sqrt(g_M)-1.)/c0, 
			 JJH2 = ((1./l_M+3.)-3./l2)/s0,
			 JJP2 = ((1./l_M*g_M+3.)-3./l2*sqrt(g_M))*.5/s0,
			 
			 JJM1 = -(1./l2+1.)*tt_M/c0, 
			 JJD1 = -(1./l2*sqrt(g_M)+1.)*tp_M/c0, 
			 JJM2 = ((1./l_M+3.)+3./l2)*tt_M/s0,
			 JJD2 = ((1./l_M*g_M+3.)+3./l2*sqrt(g_M))*tp_M/s0,

			 JHP1 = (1./(l2*RR1)-1.)*tt_D,
			 JDP1 = (1./(l2*RR1)*sqrt(g_M)-1.)*tp_D,
			 JHP2 = ((1./(l_M*fR2)+3.)-3./(l2*RR1))*tt_D,
			 JDP2 = ((1./(l_M*fR2)*g_M+3.)-3./(l2*RR1)*sqrt(g_M))*tp_D,

			 JHM1 = -(1./(l2*RR1)+1.)*tt_T,
			 JDM1 = -(1./(l2*RR1)*sqrt(g_M)+1.)*tp_T,
			 JHM2 = ((1./(l_M*fR2)+3.)+3./(l2*RR1))*tt_T,
			 JDM2 = ((1./(l_M*fR2)*g_M+3.)+3./(l2*RR1)*sqrt(g_M))*tp_T;

////////////////////////////////////////////
//...заполняем строки когезионной матрицы;
	matr[0][2] = mu_M*(g_I*HP1-3.*l_I*HP2); 	//...сшивка функций;
	matr[0][3] = mu_M*3.*l_I*HH2;
	matr[0][8] = -mu_I*(g_M*JJP1-3.*l_M*JJP2); 
	matr[0][9] = -mu_I*(g_M*JJD1-3.*l_M*JJD2); 
	matr[0][10] = -mu_I*3.*l_M*JJH2; 
	matr[0][11] = -mu_I*3.*l_M*JJM2;

	matr[1][2] = mu_M*2.*l_I*HP2; 	
	matr[1][3] = mu_M*(HH1-2.*l_I*HH2);
	matr[1][8] = -mu_I*2.*l_M*JJP2;  
	matr[1][9] = -mu_I*2.*l_M*JJD2;  
	matr[1][10] = -mu_I*(JJH1-2.*l_M*JJH2); 
	matr[1][11] = -mu_I*(JJM1-2.*l_M*JJM2);

	//...сшивка нормальных производных;
	matr[2][2] = mu_M*(g_I*(fR2*HP2-HP1)+12.*l_I*HP2); 
	matr[2][3] = mu_M*3.*(HH1-4.*l_I*HH2);
	matr[2][8] = -mu_I*(g_M*(fR2*JJP2-JJP1)+12.*l_M*JJP2); 
	matr[2][9] = -mu_I*(g_M*(fR2*JJD2-JJD1)+12.*l_M*JJD2); 
	matr[2][10] = -mu_I*3.*(JJH1-4.*l_M*JJH2); 
	matr[2][11] = -mu_I*3.*(JJM1-4.*l_M*JJM2);

	matr[3][2] = mu_M*2.*(g_I*HP1-4.*l_I*HP2); 
	matr[3][3] = mu_M*((fR2*HH2-HH1)+8.*l_I*HH2);
	matr[3][8] = -mu_I*2.*(g_M*JJP1-4.*l_M*JJP2); 
	matr[3][9] = -mu_I*2.*(g_M*JJD1-4.*l_M*JJD2); 
	matr[3][10] = -mu_I*((fR2*JJH2-JJH1)+8.*l_M*JJH2); 
	matr[3][11] = -mu_I*((fR2*JJM2-JJM1)+8.*l_M*JJM2);

	//...сшивка моментов;
	matr[4][2] = g_M*(g_I*HP1-3.*l_I*HP2); 
	matr[4][3] = g_M*3.*l_I*HH2;
	matr[4][8] = -g_I*(g_M*JJP1-3.*l_M*JJP2); 
	matr[4][9] = -g_I*(g_M*JJD1-3.*l_M*JJD2);
	matr[4][10] = -g_I*3.*l_M*JJH2; 
	matr[4][11] = -g_I*3.*l_M*JJM2;

	matr[5][2] = 2.*l_I*HP2; 	 
	matr[5][3] = HH1-2.*l_I*HH2;
	matr[5][8] = -2.*l_M*JJP2; 
	matr[5][9] = -2.*l_M*JJD2;
	matr[5][10] = -JJH1+2.*l_M*JJH2; 
	matr[5][11] = -JJM1+2.*l_M*JJM2;

	//...антисимметричный вариант (симметрия по второму и третьему индексам);
	//matr[6][2] = -((2.+g_I)*HP1-9.*l_I/g_I*HP2)*.5; 
	//matr[6][3] =  (1./g_I-1.)*.75*HH1-4.5*l_I/g_I*HH2;
	//matr[6][8] =  ((2.+g_M)*JJP1-9.*l_M/g_M*JJP2)*.5; 
	//matr[6][9] =  ((2.+g_M)*JJD1-9.*l_M/g_M*JJD2)*.5; 
	//matr[6][10] = -(1./g_M-1.)*.75*JJH1+4.5*l_M/g_M*JJH2; 
	//matr[6][11] = -(1./g_M-1.)*.75*JJM1+4.5*l_M/g_M*JJM2;

	//matr[7][2] = -(1.-g_I)*.5*HP1-(4.-1./g_I)*l_I*HP2; 
	//matr[7][3] = -(5.+1./g_I)*.25*HH1+(4.-1./g_I)*l_I*HH2;
	//matr[7][8] =  (1.-g_M)*.5*JJP1+(4.-1./g_M)*l_M*JJP2; 
	//matr[7][9] =  (1.-g_M)*.5*JJD1+(4.-1./g_M)*l_M*JJD2; 
	//matr[7][10] = (5.+1./g_M)*.25*JJH1-(4.-1./g_M)*l_M*JJH2; 
	//matr[7][11] = (5.+1./g_M)*.25*JJM1-(4.-1./g_M)*l_M*JJM2;

	//...прежний симметричный, но нефизичный вариант (симметрия по первому и второму индексам);
	matr[6][2] = -((1.+4.*g_I)*HP1-3.*(6.+1./g_I)*l_I*HP2)*.5; 
	matr[6][3] =  (HH1-(6.+1./g_I)*l_I*HH2)*1.5;
	matr[6][8] =  ((1.+4.*g_M)*JJP1-3.*(6.+1./g_M)*l_M*JJP2)*.5; 
	matr[6][9] =  ((1.+4.*g_M)*JJD1-3.*(6.+1./g_M)*l_M*JJD2)*.5; 
	matr[6][10] = -(JJH1-(6.+1./g_M)*l_M*JJH2)*1.5; 
	matr[6][11] = -(JJM1-(6.+1./g_M)*l_M*JJM2)*1.5;

	matr[7][2] = -((1.-2.*g_I)*HP1+(10.-3./g_I)*l_I*HP2); 
	matr[7][3] = -(2.*HH1-(10.-3./g_I)*l_I*HH2);
	matr[7][8] =  (1.-2.*g_M)*JJP1+(10.-3./g_M)*l_M*JJP2; 
	matr[7][9] =  (1.-2.*g_M)*JJD1+(10.-3./g_M)*l_M*JJD2; 
	matr[7][10] = 2.*JJH1-(10.-3./g_M)*l_M*JJH2; 
	matr[7][11] = 2.*JJM1-(10.-3./g_M)*l_M*JJM2;

	//...некорректный антисимметричный вариант (нет симметрии по второму и третьему индексам);
	//matr[6][2] = -((1.+3.*g_I)*HP1-3.*(4.+1./g_I)*l_I*HP2)*.5; 
	//matr[6][3] =  (HH1*.5-(4.+1./g_I)*l_I*HH2)*1.5;
	//matr[6][8] =  ((1.+3.*g_M)*JJP1-3.*(4.+1./g_M)*l_M*JJP2)*.5; 
	//matr[6][9] =  ((1.+3.*g_M)*JJD1-3.*(4.+1./g_M)*l_M*JJD2)*.5; 
	//matr[6][10] = -(JJH1*.5-(4.+1./g_M)*l_M*JJH2)*1.5; 
	//matr[6][11] = -(JJM1*.5-(4.+1./g_M)*l_M*JJM2)*1.5;

	//matr[7][2] = -((1.-1.5*g_I)*HP1+(8.-3./g_I)*l_I*HP2); 
	//matr[7][3] = -(1.75*HH1-(8.-3./g_I)*l_I*HH2);
	//matr[7][8] =  (1.-1.5*g_M)*JJP1+(8.-3./g_M)*l_M*JJP2; 
	//matr[7][9] =  (1.-1.5*g_M)*JJD1+(8.-3./g_M)*l_M*JJD2; 
	//matr[7][10] = 1.75*JJH1-(8.-3./g_M)*l_M*JJH2; 
	//matr[7][11] = 1.75*JJM1-(8.-3./g_M)*l_M*JJM2;

	//...равенство нулю когезионного поля;
	matr[8][8] = g_M*JDP1-3.*l_M*JDP2; 
	matr[8][9] = g_M*JDM1-3.*l_M*JDM2; 
	matr[8][10] = 3.*l_M*JHP2; 
	matr[8][11] = 3.*l_M*JHM2;

	matr[9][8] = 2.*l_M*JDP2; 
	matr[9][9] = 2.*l_M*JDM2; 
	matr[9][10] = JHP1-2.*l_M*JHP2; 
	matr[9][11] = JHM1-2.*l_M*JHM2;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_MH0 = 0.; nu_H = -1.;
	for (k = 0; k < 14; k++)
	for (l = 0; l < 15; l++) matrix[k][l] = matr[k][l];
	optim = take_system_shear(matrix, mu_MH0, nu_H);
	if (optim < 0.) sgn0 = -1.; else sgn0 = 1.;

	mu_MH1 = 2.*mu_M/(3.*KH);
	do { 
		mu_MH1 *= 10.; k_iter++; 
		nu_H = (1.5*KH*mu_MH1/mu_M-1.)/(3.*KH*mu_MH1/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_shear(matrix, mu_MH1, nu_H);
	}
	while(optim*sgn0 > 0. && k_iter < max_iter); k_iter = 0;

	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;
	do {
		mu_MH = (mu_MH0+mu_MH1)*.5; k_iter++;
		nu_H = (1.5*KH*mu_MH/mu_M-1.)/(3.*KH*mu_MH/mu_M+1.);
		for (k = 0; k < 14; k++)
		for (l = 0; l < 15; l++)	matrix[k][l] = matr[k][l];
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
