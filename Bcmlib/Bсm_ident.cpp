#include "stdafx.h"

#include "utils.h"
#include "cmapi.h"
#include "clame3d.h"
#include "ccebasic.h"
#include "cgrid_el.h"
#include "ccohes2d.h"
#include "ccohes3d.h"
#include "cblocked.h"

#include "cgrid.h"
#include "ccells.h"
#include "cdraft.h"
#include "kernel.h"
#include "cshapes.h"
#include "csolver.h"

#define ___nMPI_INIT___
#define  Message(Msg)    printf("%s", Msg);  printf("\n")
//----------------------------------------------------------------------------------------
template <typename T> 
void Number_of_Blocks	(void * _pBCM, int & _NBlks);
template <typename T> 
void Blocks_Partitioning(void * _pBCM, int * _Blks);
template <typename T> 
void Blocks_Sparsity		(void * _pBCM, int *& _IA, int *& _JA);
template <typename T> 
void Blocks_Row			(void * _pBCM, int _BlkRow, double * _defc);
template <typename T> 
void Right_Handside		(void * _pBCM, int _BlkRow, double * _refc);
template <typename T> 
void Initial_Guess		(void * _pBCM, int _BlkRow, double * _refc);
template <typename T> 
void Store_Solution		(void * _pBCM, int _BlkRow, double * _refc);
template <typename T> 
void Blocks_SparsitySym	(void * _pBCM, int *& _IA, int *& _JA);
template <typename T> 
void Blocks_RowSym		(void * _pBCM, int _BlkRow, double * _defc);
template <typename T> 
//----------------------------------------------------------------------------------------
#ifdef ___MPI_INIT___
#include "mpi.h"
#include "ExchangeMPI.h"
#include "gsmatrix.h"
extern CMPIComm comm_mpi;
#endif
//---------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////
//...идентификаци€ когезионных параметров (дл€ когезионной или классической модели);
double CIdent(char * name_ini, double R, double A, double * energy, double l1, double l2, double G1, double G2, double nju1, double nju2, double Ad, double Bd, int N0, int id_layer, int id_direct, int id_visual)
{
	char buf[2000];
	int  i, j, k, l, id_reading = 0;
	double X0, Y0, ell_X, ell_Y, rot_Z, aaa_X, aaa_Y, fff;

//////////////////////////////////
//...reading model from data-file;
	CCohes2D * sm = new CCohes2D;
	CGrid_el * nd = new CGrid_el;

	if (sm) {
		sprintf(buf, "Loading model from file '%s'", name_ini);
		Message(" ");
		Message(buf);
		Message("Reading data file ...");

		sm->stru.nodes_in(name_ini);
      sm->bar_condit_in(name_ini);
		sm->LinkUniStruct();
		sm->SetBUniStruct(POLY_BLOCK);

/////////////////////////////////////
//...reading parameters of inclusion;
		if (sm->id_prop && sm->pp_cond)
		for (j = 0; j < sm->id_prop[0]; j++)
		if (sm->id_prop[j*2+2] == BSOURCE_BND) {
			X0		= sm->pp_cond[j*6]; 
			Y0		= sm->pp_cond[j*6+1];
			ell_X = sm->pp_cond[j*6+2];
			ell_Y = sm->pp_cond[j*6+3];
			rot_Z = sm->pp_cond[j*6+4];
			id_reading = 1;
			break;
		}
		Message("Finish!");
	}
	if (id_layer == SPECIAL_STATE && ! id_reading) {//...окружность (дл€ сведени€);
		fff = 4./3.*M_PI*R*R*R/(A*A*A);
		double aaa = (A-R)/(1.-sqrt(fff/M_PI));
		double rad = sqrt(fff/M_PI)*aaa;
		A = aaa;
		R = rad;
//   набор экспериментальных точек:
//   rad		aaa	 aaa/2
//   0.1550 0.8312 0.4156
//   0.1686 0.6726 0.3363
//   0.1794 0.5814 0.2907
//   0.1875 0.5278 0.2639

//   0.1726 0.9422 0.4711
//   0.1890 0.7504 0.3752
//   0.2013 0.6466 0.3233
//   0.2103 0.5888 0.2944

//   0.3370 1.8400 0.9200
//   0.3701 1.4556 0.7278
//   0.3955 1.2466 0.6233
//   0.4084 1.1634 0.5817

//   0.4824 2.6578 1.3289
//   0.5318 2.0820 1.0410
//   0.5701 1.7708 0.8854
//   0.5897 1.6478 0.8239

//   0.7474 4.1176 2.0588
//   0.8216 3.2468 1.6234
//   0.8792 2.7714 1.3857
//   0.9143 2.5486 1.2743
	}

/////////////////////////////////////////////////////////////////////
//...определ€ем размеры €чейки и устанавливаем граничную поверхность;
	double par[6];	sm->SetGeomBounding(par);
	if (! id_reading) {
		X0 = (par[0]+par[1])*.5;
		Y0 = (par[2]+par[3])*.5;
		ell_X = R;
		ell_Y = R;
		rot_Z = 0.;
	}
	if (ell_X != 0. && ell_Y != 0.) {
		CCells * ce = new(CCells);
		ce->cells_new(1, 2, (l = size_of_map(1, CYL_GENUS))+1);
		ce->mp[0] = (CMap)ID_MAP(1, CYL_GENUS);
		ce->mp[1] = X0;
		ce->mp[2] = Y0;
		ce->mp[4] = rot_Z/180.*M_PI;
		ce->mp[7] = ell_X;
		ce->mp[8] = ell_Y;
		ce->mp[l] = (CMap)NULL_CELL;
		sm->bar = new(CCells);
		sm->bar->bar_add(ce);
	}

////////////////////////////////////
//...устанавливаем параметры задачи;
	if (id_layer != NULL_STATE && id_layer != SPECIAL_STATE) {
		fff	= 4./3.*M_PI*R*R*R/(A*A*A);
		if (id_layer == OK_STATE) {
			aaa_X = (A-R)/(1.-fff);
			aaa_Y = fff*aaa_X;
		}
		else {
			aaa_X = A;
			aaa_Y = ell_X;
		}
		sm->TakeLayerModel (aaa_X*.5, par[3]-par[2], aaa_Y*.5, nju1, nju2, G1, G2, l1, l2, sm->get_param(5));
		par[0] = -(par[1] = aaa_X*.5);
		par[2] = -(par[3] = (par[3]-par[2])*.5);
	}
	else {
		sm->set_mpls(PackInts(3, 3)); //...степень мультиполей;
		sm->set_quad(PackInts(8, 4)); //...степень квадратуры;
		sm->set_normaliz(0.92);
		sm->set_lagrange(1e5);
		sm->change_solv(E_PERIODIC_SOLVING);
		if (! l1 || ! l2)					//...material parameters;
		sm->set_fasa_hmg(nju1, nju2, G1, G2); else
		sm->set_fasa_hmg(nju1, nju2, G1, G2, G1/sqr(l1), G2/sqr(l2));
	}
	sm->set_param(3, Ad); //...адгезионные модули Ad и Bd;
	sm->set_param(4, Bd);

////////////////////////////
//...solving of the problem;
//	sm->solver.change(EXTERN_STATE);//...включаем внешний солвер;
//	sm->solver.set_mode(PRINT_MODE | FULLY_MODE/* | REDUCED_MESSAGE | ACCUMULATION*/);
	if (id_direct) sm->solver.change_state(OK_STATE);
	if (id_layer != NULL_STATE && id_layer != SPECIAL_STATE) sm->solver.change_state(OK_STATE);
	else {
		if (sm->computing_kernel(MAPPING_COMPUT) != OK_STATE) {
			Message("Error in sample counting...");
			if (sm) delete sm;
			if (nd) delete nd;

			return(0.);
		}
#ifdef ___MPI_INIT___
		if (sm->solver.changed(EXTERN_STATE)) {
			CSlvParam params;
			params.msglev = 3;
			params.ittype = 2;
			params.sttype = 1;
			params.niter = 500;
			params.eps = 1.0e-7;

			params.tau1 = 1.0e-2;
			params.tau2 = 1.0e-3;
			params.theta = 0.10e0;

			char strbuff[256];
			sprintf (strbuff,"%s%i%s","BsSolver_",comm_mpi.GetMyid(),".dat");

			std::ofstream fout (strbuff);

			sm->shapes_init(NO_STATE); 
			BCM_draft<double> pBCM = {sm, BASIC_COMPUT};
			AbstractParSolver (& pBCM, 
									(void *)&comm_mpi,
									fout, params,
									Number_of_Blocks<double>, Blocks_Partitioning<double>, 
									Blocks_Sparsity<double>, Blocks_Row<double>, 
									Right_Handside<double>, Initial_Guess<double>, Store_Solution<double>);
			sm->shapes_init(OK_STATE); 
		}
#endif
	}

////////////////////////////////////////////////
//..вычисление энергии (или средних напр€жений);
	double EYoung = 0.;
#ifdef ___MPI_INIT___
	if (comm_mpi.GetMyid() == 0) 
#endif
	{
		if (id_layer != NULL_STATE) EYoung = sm->TakeLayer_E1(fff);
		else {
			double K[12], E0, nu, mu, C0, lm; memset(K, 0, 12*sizeof(double));
			sm->GetRigidy(K, SPECIAL2COMPUT);

			C0 = (K[0]-2.*K[2]*K[5]/(K[5]+K[3]))/(K[3]-2.*K[5]*K[5]/(K[5]+K[3]));
			lm = (K[2]-C0*K[5])/(K[5]+K[3]);
			E0 = (C0-lm)*(C0+2.*lm)/(C0+lm);
			nu = lm/(C0+lm);
			mu = K[7]/(K[10]*2.);
			EYoung = E0;

			for (k = 0; k < sm->N; k++) sm->GetEnergyValue(k, energy);
		}
		if (id_visual) {//..visualization;
			nd->zero_grid();

			sm->BlockActivate(NULL_STATE);

			int NX = 100, NY = 100;
			for (i = 0; i <= 2*NX; i++) nd->add_new_point_X(.5*i/NX*(par[1]-par[0])+par[0]);
			for (j = 0; j <= 2*NY; j++) nd->add_new_point_Y(.5*j/NY*(par[3]-par[2])+par[2]);

			nd->hit = (int *)new_struct(nd->N*nd->N1*sizeof(int));

			for (i = 0; i < nd->N;  i++)
			for (j = 0; j < nd->N1; j++) {
				int hit = -1;
				sm->Poly_struc_in2D (hit, nd->X[i], nd->Y[j]);
				sm->StructEllCorrect(hit, nd->X[i], nd->Y[j]);
				nd->hit[i+j*nd->N] = hit;
			}
			char buff[1001], name[2001]; strcpy(buff, name_ini);
			char * pchar = strrchr(buff, '.'); 
			if (pchar != NULL) pchar[0] = 0;
			sprintf(name, "%s_%i_%i", buff, (int)(100*l1), (int)(100*l2));
//			sm->GetSurferFormat(name, nd, ENERGY_VALUE);

			int res = system("del *.grd");
			sm->GetSurferFormat("bb", nd,		ERR_VALUE);
			sm->GetSurferFormat("pp", nd,  DISPL_VALUE);		
			sm->GetSurferFormat("ee", nd, ENERGY_VALUE);
			if (! l1 || ! l2) {
				sm->GetSurferFormat("tt_x", nd, STRESS_X_VALUE);
				sm->GetSurferFormat("tt_y", nd, STRESS_Y_VALUE);
			}
			else {
				sm->GetSurferFormat("tf_x", nd, STRESS_X_VALUE);
				sm->GetSurferFormat("tf_y", nd, STRESS_Y_VALUE);
			}
//			sm->GetSurferFormat("hh_x", nd,				  DILAT_VALUE);
//   		sm->GetSurferFormat("th_x", nd, STRESS_X_COHESION_VALUE);
//			sm->GetSurferFormat("un", nd,				  NORMAL_X_VALUE);
//			sm->GetSurferFormat("uu", nd,			DISPL_CLASSIC_VALUE);
		}
		int id_action = 0;
		if (id_action)	{
			double sum = 0.; //...среднее значение перемещений;
			for (j = nd->N1/2, i = 0; i < nd->N; i++) {
				int hit = nd->hit[i+j*nd->N];
				double X = nd->X[i], 
						 Y = nd->Y[j], 
						 Z = 0., out_F[3];
				sm->GetFuncAllValues (X, Y, Z, out_F, hit, DISPL_VALUE, 0);
				sum += out_F[0];
			}	sum /= nd->N;

			FILE * TST = fopen("Cohes_middle.dat", "w");
			for (j = nd->N1/2, i = 0; i < nd->N; i++) {
				int hit = nd->hit[i+j*nd->N];
				double X = nd->X[i], 
						 Y = nd->Y[j], 
						 Z = 0., out_F[3];
				sm->GetFuncAllValues (X, Y, Z, out_F, hit, DISPL_VALUE, 0);
				fprintf(TST, " %g  %g", X, out_F[0]-sum);

				sm->GetFuncAllValues (X, Y, Z, out_F, hit, ENERGY_VALUE, 0);
				fprintf(TST, " %g", out_F[0]);

				sm->GetFuncAllValues (X, Y, Z, out_F, hit, STRESS_X_VALUE, 0);
				fprintf(TST, " %g", out_F[0]);

				sm->GetFuncAllValues (X, Y, Z, out_F, hit, DILAT_VALUE, 0);
				fprintf(TST, " %g\n", out_F[0]);
			}
			fclose(TST);
		}
	}
	delete sm;
	delete nd;

	return(EYoung);
}

/////////////////////////////////////////////////////////////////////
//...идентификаци€ когезионных параметров дл€ сферического включени€;
double CIdent_sphere3D(double * energy, double rad, double fV, double l1, double l2, double G1, double G2, double nju1, double nju2, int N0, int id_visual)
{
	double AA = pow(4./3.*M_PI/fV, 1./3.)*rad;

////////////////////////////////////////////////
//...образуем модель дл€ сферического включени€;
	CCohes3D * sm = (CCohes3D *)(! l1 || ! l2 ? new CLame3D : new CCohes3D);

	sm->GetSphBoxStruct (AA, AA, AA, rad);
	sm->SetBUniStruct(NULL_BLOCK);

	for (int k = 0; k < sm->N; k++) {
		if (sm->B[k].link[sm->NUM_PHASE] == -1) sm->B[k].type = ZOOM_BLOCK;
		else												 sm->B[k].type = POLY_BLOCK;
	}

////////////////////////////////////
//...устанавливаем параметры задачи;
	sm->set_param(0, PackInts(N0+2, N0));	//...multipoles degree;
	sm->set_param(1, 1.);						//...normalization coeffitient;
	sm->set_param(3, PackInts(20, 0));		//...quadrature degree;
	sm->set_param(sm->size_of_param()-1, 1000000.); //...Lagrange corfficient for LSM;

	if (sm->type() == LAME3D_DRAFT) //...material parameters for interphase layer problem;
	sm->set_fasa_hmg(nju1, nju2, G1, G2, nju1, G1+G2*0.2); else
	sm->set_fasa_hmg(nju1, nju2, G1, G2, G1/sqr(l1), G2/sqr(l2));

///////////////////////////
//...solving of the probem;
//	sm->solver.set_mode(PRINT_MODE);
	if (sm->computing_kernel(MAPPING_COMPUT) != OK_STATE) {
		Message("Error in sample counting...");
		delete sm;
		return(0.);
	}

///////////////////////////////////////////////////////////////////
//..вычисл€ем модуль ёнга и воспроизводим распределение напр€жений;
	double EYoung = 0.;
#ifdef ___MPI_INIT___
	if (comm_mpi.GetMyid() == 0) 
#endif
	{
		double K[24], E0, nu, mu, C0, lm; memset(K, 0, 24*sizeof(double));
		sm->GetRigidy(K, SPECIAL2COMPUT);

		C0 = (K[5]-2.*K[0]*K[6]/(K[6]+K[11]))/(K[11]-2.*K[6]*K[6]/(K[6]+K[11]));
		lm = (K[0]-C0*K[6])/(K[6]+K[11]);
		E0 = (C0-lm)*(C0+2.*lm)/(C0+lm);
		nu = lm/(C0+lm);
		mu = K[16]/(K[22]*2.);
		EYoung = E0;

		if (id_visual) { //..visualization;
			int NX = 50, NY = 50, axis, i, j;
			CGrid_el * nd = new CGrid_el;

			double par[6];
			par[0] = -AA*.5; par[2] = -AA*.5; par[4] = -AA*.5;
			par[1] =  AA*.5; par[3] =  AA*.5; par[5] =  AA*.5;

#ifdef ___LONGITUDINAL_SECTIOM___ 
			for (i = 0; i <= 2*NX; i++) nd->add_new_point_X(.5*i/NX*(par[1]-par[0])+par[0]);
			for (j = 0; j <= 2*NY; j++) nd->add_new_point_Y(.5*j/NY*(par[3]-par[2])+par[2]);

			nd->add_new_point_Z((par[5]+par[4])*.5);
			nd->hit = (int *)new_struct(nd->N*nd->N1*sizeof(int));

			for (i = 0; i < nd->N;  i++)
			for (j = 0; j < nd->N1; j++)
			if (sqr(nd->X[i])+sqr(nd->Y[j])+sqr(nd->Z[0]) < sqr(rad)) 
				nd->hit[i+j*nd->N] = 0; else  
				nd->hit[i+j*nd->N] = 1;	

			axis = AXIS_Z;
#else
			if (1) {
				for (i = 0; i <= 2*NX; i++) nd->add_new_point_X(.5*i/NX*(par[5]-par[4])+par[4]);
				for (j = 0; j <= 2*NY; j++) nd->add_new_point_Y(.5*j/NY*(par[1]-par[0])+par[0]);

				nd->add_new_point_Z((par[3]+par[2])*.5);
				nd->hit = (int *)new_struct(nd->N*nd->N1*sizeof(int));

				for (i = 0; i < nd->N;  i++)
				for (j = 0; j < nd->N1; j++)
				if (sqr(nd->X[i])+sqr(nd->Y[j])+sqr(nd->Z[0]) < sqr(rad)) 
					nd->hit[i+j*nd->N] = 0; else  
					nd->hit[i+j*nd->N] = 1;	

				axis = AXIS_Y;
			}
			else {
				for (i = 0; i <= 2*NX; i++) nd->add_new_point_X(.5*i/NX*(par[3]-par[2])+par[2]);
				for (j = 0; j <= 2*NY; j++) nd->add_new_point_Y(.5*j/NY*(par[5]-par[4])+par[4]);

				nd->add_new_point_Z((par[1]+par[0])*.5);
				nd->hit = (int *)new_struct(nd->N*nd->N1*sizeof(int));

				for (i = 0; i < nd->N;  i++)
				for (j = 0; j < nd->N1; j++)
				if (sqr(nd->X[i])+sqr(nd->Y[j])+sqr(nd->Z[0]) < sqr(rad)) 
					nd->hit[i+j*nd->N] = 0; else  
					nd->hit[i+j*nd->N] = 1;	

				axis = AXIS_X;
			}
#endif
			int res = system("del *.grd");
			sm->GetSurferFormat("pp",  nd,	 DISPL_JUMP_VALUE, 0, axis);
			sm->GetSurferFormat("ttz", nd, STRESS_Z_JUMP_VALUE, 0, axis);
			//sm->GetSurferFormat("pp",  nd,	 DISPL_UNIX_VALUE, AXIS_Z, axis);
			//sm->GetSurferFormat("ttz", nd, STRESS_Z_UNIX_VALUE, AXIS_Z, axis);
		}
	}
	delete sm;
	return(EYoung);
}

/////////////////////////////////////////////////////////////////////
//...одномерна€ модель межфазного сло€ (симметрична€ слоиста€ среда);
double TakeLayerHomogeniz(double L, double H, double l, double k1, double k2, double C1, double C2, double A)
{
	double kk1 = sqrt(C1/k1), t1 = exp(-kk1*(L-l)), tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/k2), t2 = exp(-kk2*l),		tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = k2*al1+k1*al2-A*al1*al2,
			 BB = k2-k1-A*al2,
			 CC = k2-k1+A*al1,
			 QQ = AA*(k2*(L-l)+k1*l+A)-BB*CC;
	return(k1*k2*AA)/QQ;
}
