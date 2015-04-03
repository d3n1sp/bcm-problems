#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "parfbss.h"
#include "parfbss_data.h"
#include "globals.h"
#include "slvparam.h"
#include "smatrix.h"
#include "gsmatrix.h"
#include "tree.h"
#include "mvm.h"

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

static ParFBSS_Array myObjPar;

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_Create (void *_Comm, int _N, int _NZA, int *_ObjId) {

	int i;

	int iobj = -1;

	for (i=0;i<myObjPar.NumberOfObjects;i++) {
		if (myObjPar.IRegistration[i] == 0) {
			iobj = i;
			break;
		};
	};

	if (iobj < 0) return -1;

	myObjPar.IRegistration[iobj] = 1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_Create (_Comm, _N, _NZA, myObjPar.PRegistration[iobj]);

	if (ierr != 0) return ierr;

	*_ObjId = iobj;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_ReCreate (int _ObjId, int _NI, int _NZAI) {

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_ReCreate (myObjPar.PRegistration[_ObjId], _NI, _NZAI);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_LinearEquationsData (int _ObjId, int **_LIST, int **_IA, int **_JA, double **_A, 
											double **_RightHandSide, double **_InitialGuess) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_LinearEquationsData (myObjPar.PRegistration[iobj], *_LIST, *_IA, *_JA, *_A, 
																		*_RightHandSide, *_InitialGuess);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_NormalizationData (int _ObjId, double **_Norm) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_NormalizationData (myObjPar.PRegistration[iobj], *_Norm);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_NonLinearEquationsData (int _ObjId, double **_Poly2) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_NonLinearEquationsData (myObjPar.PRegistration[iobj], *_Poly2);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_SetSolverParameters (int _ObjId, 
											int _NpDistr,
											double _Tau1, double _Tau2, double _Theta, 
											int _NIterMax, double _Eps) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_SetSolverParameters (myObjPar.PRegistration[iobj], 
																		_NpDistr, _Tau1, _Tau2, _Theta, 
																		_NIterMax, _Eps);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_SetMsglevAndFilename (int _ObjId, int _Msglev, char *_FileName) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_SetMsglevAndFilename (myObjPar.PRegistration[iobj], _Msglev, _FileName);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_SetSchurtype (int _ObjId, int _Schurtype, int _Ncycle) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_SetSchurtype (myObjPar.PRegistration[iobj], _Schurtype, _Ncycle);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_SetNcpuext (int _ObjId, int _NcpuExt) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_SetNcpuext (myObjPar.PRegistration[iobj], _NcpuExt);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_SolveSymm (int _ObjId, int _Job, int *_NIterPerf) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_SolveSymm (myObjPar.PRegistration[iobj], _Job, _NIterPerf);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_Solve (int _ObjId, int _Job, int *_NIterPerf) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_Solve (myObjPar.PRegistration[iobj], _Job, _NIterPerf);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_Solution (int _ObjId, double **_Solution) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_Solution (myObjPar.PRegistration[iobj], *_Solution);

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_Distribution (int _ObjId, int *_NIOpt, int **_LISTOpt) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	int nlistloc;
	int *listloc;

	ierr = ParFBSS_Array::ParFBSS_Distribution (myObjPar.PRegistration[iobj], nlistloc, listloc);

	*_NIOpt = nlistloc;
	*_LISTOpt = listloc;

	if (ierr != 0) return ierr;

	return 0;

};

// Author: Kharchenko S.A. 
//========================================================================================
int ParFBSS_Destroy (int _ObjId) {

	int iobj = _ObjId;

	if (iobj < 0 || iobj >= myObjPar.NumberOfObjects) return -1;

	int ierr;

	ierr = ParFBSS_Array::ParFBSS_Destroy (myObjPar.PRegistration[iobj]);

	if (ierr != 0) return ierr;

	myObjPar.IRegistration[iobj] = 0;

	return 0;

};

// Implement interface functions

// Author: Kharchenko S.A. 
// ParFBSS_Array: create an object
//========================================================================================
int ParFBSS_Array::ParFBSS_Create (void *_Comm, int _NI, int _NZAI, CParFbssData *&_PObj) {

	_PObj = new CParFbssData (_Comm, _NI, _NZAI);
	if (!_PObj) return 1;

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: recreate an object
//========================================================================================
int ParFBSS_Array::ParFBSS_ReCreate (void *_PObj, int _NI, int _NZAI) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSMatrixR *pa = pparfbssdata->GetAmtr ();
	CSVector *prhs = pparfbssdata->GetRhs ();
	CSVector *px = pparfbssdata->GetX ();
	CSVector *psol = pparfbssdata->GetSol ();
	CSVector *pnorm = pparfbssdata->GetNorm ();
	CSVector *ppoly2 = pparfbssdata->GetPoly2 ();

	delete pa;
	delete prhs;
	delete px;
	delete psol;
	delete pnorm;
	delete ppoly2;

	pa = new CSMatrixR (_NI, _NZAI);
	prhs = new CSVector (_NI);
	px = new CSVector (_NI);
	psol = new CSVector (_NI);
	pnorm = new CSVector (_NI);
	ppoly2 = new CSVector (_NI);

	int i;
	double *pvect;

	pvect = pnorm->GetVect ();
	for (i=0;i<_NI;i++) pvect[i] = 1.0e0;
	pvect = ppoly2->GetVect ();
	for (i=0;i<_NI;i++) pvect[i] = 0.0e0;

	pparfbssdata->SetAmtr (pa);
	pparfbssdata->SetRhs (prhs);
	pparfbssdata->SetX (px);
	pparfbssdata->SetSol (psol);
	pparfbssdata->SetNorm (pnorm);
	pparfbssdata->SetPoly2 (ppoly2);

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: store data
//========================================================================================
int ParFBSS_Array::ParFBSS_LinearEquationsData (void *_PObj, int *&_LIST, int *&_IA, int *&_JA, double *&_A, 
																double *&_RightHandSide, double *&_InitialGuess) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSMatrixR *pa = pparfbssdata->GetAmtr ();

	_LIST = pa->GetList ();
	_IA = pa->GetIa ();
	_JA = pa->GetJa ();
	_A = pa->GetA ();

	CSVector *prhs = pparfbssdata->GetRhs ();
	CSVector *px = pparfbssdata->GetX ();

	_RightHandSide = prhs->GetVect ();
	_InitialGuess = px->GetVect ();

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: store data
//========================================================================================
int ParFBSS_Array::ParFBSS_NormalizationData (void *_PObj, double *&_Norm) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSVector *pnorm = pparfbssdata->GetNorm ();

	_Norm = pnorm->GetVect ();

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: store data
//========================================================================================
int ParFBSS_Array::ParFBSS_NonLinearEquationsData (void *_PObj, double *&_Poly2) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSVector *ppoly2 = pparfbssdata->GetPoly2 ();

	_Poly2 = ppoly2->GetVect ();

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: get parameters
//========================================================================================
int ParFBSS_Array::ParFBSS_GetParameters (void *_PObj, void *&_params) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSlvParam *pparams = pparfbssdata->GetParams ();

	_params = (void *)pparams;

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: set parameters
//========================================================================================
int ParFBSS_Array::ParFBSS_SetSolverParameters (void *_PObj, 
																int _NpDistr, 
																double _Tau1, double _Tau2, double _Theta, 
																int _NIterMax, double _Eps) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSlvParam *pparams = pparfbssdata->GetParams ();

	pparams->ncpupart = _NpDistr;
	pparams->msglev = 0;
	pparams->collap = 0;

	pparams->ordtype = 1;

	pparams->fcttyp = 2;
	pparams->memory = -0.1e0;
	pparams->nfreer = 5;
	pparams->tau1 = _Tau1;
	pparams->tau2 = _Tau2;
	pparams->theta = _Theta;

	pparams->niter = _NIterMax;

	pparams->eps = _Eps;

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: set parameters
//========================================================================================
int ParFBSS_Array::ParFBSS_SetMsglevAndFilename (void *_PObj, int _Msglev, char *_Filename) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSlvParam *pparams = pparfbssdata->GetParams ();

	pparams->msglev = _Msglev;

	strcpy (pparams->path,_Filename);

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: set parameters
//========================================================================================
int ParFBSS_Array::ParFBSS_SetSchurtype (void *_PObj, int _Schurtype, int _Ncycle) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSlvParam *pparams = pparfbssdata->GetParams ();

	pparams->schurtype = _Schurtype;
	pparams->ncycle = _Ncycle;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: set parameters
//========================================================================================
int ParFBSS_Array::ParFBSS_SetNcpuext (void *_PObj, int _Ncpuext) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSlvParam *pparams = pparfbssdata->GetParams ();

	pparams->ncpuext = _Ncpuext;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: solve
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve (void *_PObj, int _Job, int *_NIterPerf) {

//	const char* funcname = "ParFBSS_Solve";

// Init the data

	int info;

	info = ParFBSS_Solve_Ini (_PObj);
	if (info != 0) return info;

// Switch the cases

	if (_Job == 1) {

		info = ParFBSS_Solve_Distrib (_PObj);
		if (info != 0) return info;

	} else {

		info = ParFBSS_Solve_GetDistrib (_PObj);
		if (info != 0) return info;

	};

// Finilize the computations

	info = ParFBSS_Solve_Transform (_PObj);
	if (info != 0) return info;

	double resfin;

	info = ParFBSS_Solve_LinearIter (_PObj, _NIterPerf, &resfin);
	if (info != 0) return info;

	info = ParFBSS_Solve_ExchangeSolution (_PObj);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: solve
//========================================================================================
int ParFBSS_Array::ParFBSS_SolveSymm (void *_PObj, int _Job, int *_NIterPerf) {

//	const char* funcname = "ParFBSS_SolveSymm";

// Init the data

	int info;

	info = ParFBSS_Solve_Ini (_PObj);
	if (info != 0) return info;

// Switch the cases

	if (_Job == 1) {

		info = ParFBSS_Solve_Distrib (_PObj);
		if (info != 0) return info;

	} else {

		info = ParFBSS_Solve_GetDistrib (_PObj);
		if (info != 0) return info;

	};

// Finilize the computations

	info = ParFBSS_Solve_Transform (_PObj);
	if (info != 0) return info;

	double resfin;

	info = ParFBSS_SolveSymm_LinearIter (_PObj, _NIterPerf, &resfin);
	if (info != 0) return info;

	info = ParFBSS_Solve_ExchangeSolution (_PObj);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: initialize restart solve
//========================================================================================
int ParFBSS_Array::ParFBSS_RestartSolve_Ini (void *_PObj, int _Job) {

//	const char* funcname = "ParFBSS_RestartSolve_Ini";

// Init the data

	int info;

	info = ParFBSS_Solve_Ini (_PObj);
	if (info != 0) return info;

// Switch the cases

	if (_Job == 1) {

		info = ParFBSS_Solve_Distrib (_PObj);
		if (info != 0) return info;

	} else {

		info = ParFBSS_Solve_GetDistrib (_PObj);
		if (info != 0) return info;

	};

// Finalize the computations

	info = ParFBSS_Solve_Transform (_PObj);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: Perform linear iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_RestartSolve_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_RestartSolve_LinearIter";

// Perform linear iterations

	int info;

	info = ParFBSS_Solve_LinearIter (_PObj, _NIterPerf, _ResFin);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: Perform linear iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_RestartSolveSymm_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_RestartSolveSymm_LinearIter";

// Perform linear iterations

	int info;

	info = ParFBSS_SolveSymm_LinearIter (_PObj, _NIterPerf, _ResFin);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: Perform nonlinear iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_RestartSolve_NonLinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_RestartSolve_NonLinearIter";

// Perform nonlinear iterations

	int info;

	info = ParFBSS_Solve_NonLinearIter (_PObj, _NIterPerf, _ResFin);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: finalize restart solve
//========================================================================================
int ParFBSS_Array::ParFBSS_RestartSolve_Fin (void *_PObj) {

//	const char* funcname = "ParFBSS_RestartSolve_Fin";

// Finalize the computations

	int info;

	info = ParFBSS_Solve_ExchangeSolution (_PObj);
	if (info != 0) return info;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: initialize data for solver
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_Ini (void *_PObj) {

//	const char* funcname = "ParFBSS_Solve_ini";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CSMatrixR *pamtr = pparfbssdata->GetAmtr ();
	CSVector *prhs = pparfbssdata->GetRhs ();
	CSVector *px = pparfbssdata->GetX ();
	CSVector *psol = pparfbssdata->GetSol ();
	CSVector *pnorm = pparfbssdata->GetNorm ();
	CSVector *ppoly2 = pparfbssdata->GetPoly2 ();
//	CSlvParam *pparam = pparfbssdata->GetParams ();

// Get nproc and myid

//	int nproc = pcomm->GetNproc ();
//	int myid = pcomm->GetMyid ();

// Check that matrix data are correct

	int *plist, *pia, *pja;
	double *pa;

	int nlistloc = pamtr->GetNlist ();
	int nzjaloc = pamtr->GetNzja ();

	plist = pamtr->GetList ();
	pia = pamtr->GetIa ();
	pja = pamtr->GetJa ();
	pa = pamtr->GetA ();

	int ntot = -1;

	int i;

	for (i=0;i<nlistloc;i++) {
		if (plist[i] > ntot) ntot = plist[i];
	};

	CMPIExchange::ExchangeArrayMPI (*pcomm,
												INTEGERVALUE, MAXIMUM, 1, &ntot, &ntot);

	ntot++;

	pparfbssdata->SetNtot (ntot);

	int j;

	if (pia[0] != 0) return -1;

	for (i=0;i<nlistloc;i++) if (pia[i+1] <= pia[i]) return -1;

	for (i=0;i<nzjaloc;i++) if (pja[i] < 0 || pja[i] >= ntot) return -1;

	pamtr->SetM (ntot);
	pamtr->SetN (ntot);
	pamtr->SetNsupc (ntot);
	pamtr->SetNsupr (ntot);

// Sort column indices

	int nimax = 0;
	int ni;

	for (i=0;i<nlistloc;i++) {
		ni = pia[i+1]-pia[i];
		if (ni > nimax) nimax = ni;
	};

	CIntDouble *idarr;

	idarr = new CIntDouble [nimax];
	if (!idarr) return -1;

	int ibeg, jloc;

	for (i=0;i<nlistloc;i++) {
		ni = pia[i+1]-pia[i];
		ibeg = pia[i];
		for (j=pia[i];j<pia[i+1];j++) {
			jloc = j-ibeg;
			idarr[jloc].intvalue = pja[j];
			idarr[jloc].dvalue = pa[j];
		};
		qsort (idarr,ni,sizeof(CIntDouble),CIntDouble::CompareIntDouble);
		for (j=pia[i];j<pia[i+1];j++) {
			jloc = j-ibeg;
			pja[j] = idarr[jloc].intvalue;
			pa[j] = idarr[jloc].dvalue;
		};
	};

	delete [] idarr;

// Perform local rows sort and store the ordering array

	CIntInt *iiarr;

	iiarr = new CIntInt [nlistloc];
	if (!iiarr) return -1;

	for (i=0;i<nlistloc;i++) {
		iiarr[i].intvalue = plist[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr,nlistloc,sizeof(CIntInt),CIntInt::CompareIntInt);

	int *piorderini = pparfbssdata->GetIorderini ();;
	int *plistini = pparfbssdata->GetListini ();
	int *plistord = pparfbssdata->GetListord ();

	delete [] piorderini;
	delete [] plistini;
	delete [] plistord;

	pparfbssdata->SetNlistini (nlistloc);

	piorderini = new int [nlistloc];
	if (!piorderini) return -1;
	plistini = new int [nlistloc];
	if (!plistini) return -1;
	plistord = new int [nlistloc];
	if (!plistord) return -1;

	pparfbssdata->SetIorderini (piorderini);
	pparfbssdata->SetListini (plistini);
	pparfbssdata->SetListord (plistord);

	for (i=0;i<nlistloc;i++) piorderini[i] = iiarr[i].int2value;
	for (i=0;i<nlistloc;i++) plist[i] = iiarr[i].intvalue;
	for (i=0;i<nlistloc;i++) plistini[i] = iiarr[i].intvalue;

	delete [] iiarr;

	int *ialoc;
	int *jaloc;
	double *aloc;

	ialoc = new int [nlistloc+1];
	if (!ialoc) return -1;
	jaloc = new int [nzjaloc];
	if (!jaloc) return -1;
	aloc = new double [nzjaloc];
	if (!aloc) return -1;

	ialoc[0] = 0;

	int iold;

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		ialoc[i+1] = pia[iold+1]-pia[iold];
	};

	for (i=0;i<nlistloc;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	int ibegold;

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		ibeg = ialoc[i];
		ibegold = pia[iold];
		ni = ialoc[i+1]-ialoc[i];
		for (j=0;j<ni;j++) {
			jaloc[ibeg+j] = pja[ibegold+j];
			aloc[ibeg+j] = pa[ibegold+j];
		};
	};

	for (i=0;i<=nlistloc;i++) pia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) pja[i] = jaloc[i];
	for (i=0;i<nzjaloc;i++) pa[i] = aloc[i];

	delete [] ialoc;
	delete [] jaloc;
	delete [] aloc;

// Reorder rhs, sol and x

	double *vectloc;

	vectloc = new double [nlistloc];
	if (!vectloc) return -1;

	double *pvect = prhs->GetVect ();

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		vectloc[i] = pvect[iold];
	};
	for (i=0;i<nlistloc;i++) pvect[i] = vectloc[i];

	pvect = px->GetVect ();

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		vectloc[i] = pvect[iold];
	};
	for (i=0;i<nlistloc;i++) pvect[i] = vectloc[i];

	pvect = psol->GetVect ();

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		vectloc[i] = pvect[iold];
	};
	for (i=0;i<nlistloc;i++) pvect[i] = vectloc[i];

	pvect = pnorm->GetVect ();

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		vectloc[i] = pvect[iold];
	};
	for (i=0;i<nlistloc;i++) pvect[i] = vectloc[i];

	pvect = ppoly2->GetVect ();

	for (i=0;i<nlistloc;i++) {
		iold = piorderini[i];
		vectloc[i] = pvect[iold];
	};
	for (i=0;i<nlistloc;i++) pvect[i] = vectloc[i];

	delete [] vectloc;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: distribute data
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_Distrib (void *_PObj) {

	const char* funcname = "ParFBSS_Solve_Distrib";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	int ntot = pparfbssdata->GetNtot ();
	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CSMatrixR *pamtr = pparfbssdata->GetAmtr ();
	CTree *ptreedata = pparfbssdata->GetTree ();
	CSlvParam *pparam = pparfbssdata->GetParams ();

// Get nproc and myid

	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Compute partitioning if necessary and new sparsity data

	CTree treeloc;

	CSMatrix mtrfin;

	int nhcells, *hcells, *hcells2cpu;
	int *orderloc;

	CSMatrix mtrdummy;

	int ibegord = 0;
	int iendord = -1;
	int npord = 0;

// Create decomposition partitioning

	int nprocpart = pparam->ncpupart;

	if (nprocpart <= 0) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	if (nprocpart > nproc) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	if (pparam->schurtype > 0) nprocpart = 1;

	int ni = ntot / nprocpart;

	int *blkspart;

	blkspart = new int [nprocpart+1];
	if (!blkspart) return -1;

	blkspart[0] = 0;

	int i;

	for (i=0;i<nprocpart-1;i++) {
		blkspart[i+1] = blkspart[i] + ni;
	};
	blkspart[nprocpart] = ntot;

// Exchange the sparsity to get ready for distribution

	ibegord = 0;
	iendord = -1;
	npord = 0;

	if (myid < nprocpart) {
		if (blkspart[myid+1] > blkspart[myid]) {
			npord = 1;
			ibegord = blkspart[myid];
			iendord = blkspart[myid+1]-1;
		};
	};

	int *ibegordarr;
	int *iendordarr;

	ibegordarr = pparfbssdata->GetIbegpartarr ();
	iendordarr = pparfbssdata->GetIendpartarr ();

	delete [] ibegordarr;
	delete [] iendordarr;

	ibegordarr = new int [1];
	if (!ibegordarr) return -1;
	iendordarr = new int [1];
	if (!iendordarr) return -1;

	ibegordarr[0] = ibegord;
	iendordarr[0] = iendord;

	pparfbssdata->SetNppart (npord);
	pparfbssdata->SetIbegpartarr (ibegordarr);
	pparfbssdata->SetIendpartarr (iendordarr);

	CSMatrix mtrnew;

	mtrnew = ((CSMatrix *)pamtr)->GetSubmatrixViaIntervals (*pcomm, npord, ibegordarr, iendordarr);

// Create local communicator

	int *listcpu;

	listcpu = new int [nprocpart];
	if (!listcpu) return -1;

	for (i=0;i<nprocpart;i++) listcpu[i] = i;

	CMPIComm commloc;

	commloc = pcomm->CreateComm (nprocpart, listcpu);

	delete [] listcpu;

// Distribute via ParMetis

	int nprocdist = 1;
	while (nprocdist*2 <= nproc) nprocdist *= 2;

	int nprocdist_ext = nprocdist;

	int nprocext = pparam->ncpuext;

	if (nprocext > nprocdist_ext) {
		while (nprocdist_ext*2 <= nprocext) nprocdist_ext *= 2;
	};

	pparfbssdata->SetNprocsolve (nprocdist);

	if (myid < nprocpart) {

// Distribute point data

		CTree treeloc_ext;

		mtrnew.OrderNDParMetis (commloc, nprocdist_ext,
										treeloc_ext, orderloc, nhcells, hcells, hcells2cpu);

		if (pparam->schurtype > 0) {

// Redistribute the boundary if necessary

			treeloc_ext.ReassignBlockNumbers (hcells);

			delete [] hcells;
			delete [] hcells2cpu;

			CSMatrix aoloc = mtrnew.OrdMtr (orderloc);

			CSMatrix mtrdummy;

			mtrnew = mtrdummy;

			aoloc.PartBinaryTreeNDSchurBnd_FixedInitialOrdering (treeloc_ext, orderloc, *pparam,
																				nhcells, hcells, hcells2cpu);

		};

		if (nprocdist_ext != nprocdist) {

			if (pparam->schurtype > 0) {
				treeloc.IncludeBinaryTreeSchur (nprocdist, treeloc_ext);
			} else {
				treeloc.IncludeBinaryTree (nprocdist, treeloc_ext);
			};

			delete [] hcells;
			delete [] hcells2cpu;

			nhcells = treeloc.GetNblkstree ();

			hcells = new int [nhcells+1];
			if (!hcells) return -1;
			hcells2cpu = new int [nhcells];
			if (!hcells2cpu) return -1;

			int *pblkstree = treeloc.GetBlkstree ();
			int *pblk2cputree = treeloc.GetBlk2cputree ();

			for (i=0;i<=nhcells;i++) hcells[i] = pblkstree[i];
			for (i=0;i<nhcells;i++) hcells2cpu[i] = pblk2cputree[i];

		} else {
			treeloc = treeloc_ext;
		};

	} else {

		nhcells = 0;

		hcells = new int [nhcells+1];
		if (!hcells) return -1;
		hcells2cpu = new int [nhcells+1];
		if (!hcells2cpu) return -1;

		orderloc = new int [0];
		if (!orderloc) return -1;

	};

// Exchange tree and blocks information

	int isize = 0;
	char *ptree;

	if (myid == 0) {
		treeloc.PackTree (isize,ptree);
	};

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, 1, &isize, &isize);

	if (myid > 0) {
		ptree = new char [isize];
		if (!ptree) MemoryFail (funcname);
		for (i=0;i<isize;i++) ptree[i] = 0;
	};

	CMPIExchange::ExchangeArrayMPI (*pcomm, CHARVALUE, ADD, isize, ptree, ptree);

	treeloc.UnPackTree (isize,ptree);

	delete [] ptree;

	if (myid > 0) {
		nhcells = 0;
	};

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, 1, &nhcells, &nhcells);

	if (myid > 0) {

		delete [] hcells;
		delete [] hcells2cpu;

		hcells = new int [nhcells+1];
		if (!hcells) MemoryFail (funcname);
		hcells2cpu = new int [nhcells];
		if (!hcells2cpu) MemoryFail (funcname);

		for (i=0;i<nhcells+1;i++) hcells[i] = 0;
		for (i=0;i<nhcells;i++) hcells2cpu[i] = 0;

	};

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, nhcells+1, hcells, hcells);
	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, nhcells, hcells2cpu, hcells2cpu);

// Store tree, partitionings and global ordering

	*ptreedata = treeloc;

	int *phcells = pparfbssdata->GetHcells ();
	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

	delete [] phcells;
	delete [] phcells2cpu;

	pparfbssdata->SetNhcells (nhcells);
	pparfbssdata->SetHcells (hcells);
	pparfbssdata->SetHcells2cpu (hcells2cpu);

	int *porderpart = pparfbssdata->GetOrderpart ();

	delete [] porderpart;

	pparfbssdata->SetOrderpart (orderloc);

// Reorder the list

	int nlistini = pparfbssdata->GetNlistini ();
	int *plistini = pparfbssdata->GetListini ();
	int *plistord = pparfbssdata->GetListord ();

	CSMatrix::OrderViaIntervals (*pcomm,
											nlistini, plistini,
											npord, ibegordarr, iendordarr,
											orderloc,
											plistord);

// Create blocks partitioning according to the tree

	int *ibegarr, *iendarr;

	ibegarr = new int [nhcells];
	if (!ibegarr) return -1;
	iendarr = new int [nhcells];
	if (!iendarr) return -1;

	int npartsloc = 0;

	for (i=0;i<nhcells;i++) {
		if (hcells2cpu[i] == myid) {
			if (hcells[i+1] > hcells[i]) {
				ibegarr[npartsloc] = hcells[i];
				iendarr[npartsloc] = hcells[i+1]-1;
				npartsloc++;
			};
		};
	};

// Store a part of global inverse order for own list of rows

	int nlisttree;
	int *invordertree;

	CSMatrix::InvOrderViaIntervals (*pcomm,
												nlistini, plistini, plistord,
												npartsloc, ibegarr, iendarr, 
												nlisttree, invordertree);

	int *listtree;

	listtree = new int [nlisttree];
	if (!listtree) return -1;

	nlisttree = 0;

	int j;

	for (i=0;i<nhcells;i++) {
		if (hcells2cpu[i] == myid) {
			for (j=hcells[i];j<hcells[i+1];j++) {
				listtree[nlisttree] = j;
				nlisttree++;
			};
		};
	};

// Store computed inverse ordering in the structure

	int *plisttree = pparfbssdata->GetListtree();
	int *pinvordertree = pparfbssdata->GetInvordertree();

	delete [] plisttree;
	delete [] pinvordertree;

	pparfbssdata->SetNlisttree(nlisttree);
	pparfbssdata->SetListtree(listtree);
	pparfbssdata->SetInvordertree(invordertree);

// Free work arrays

	delete [] blkspart;
	delete [] ibegarr;
	delete [] iendarr;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: check that data are correct and get ready for solve
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_GetDistrib (void *_PObj) {

//	const char* funcname = "ParFBSS_Solve_GetDistrib";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	int ntot = pparfbssdata->GetNtot ();
	CMPIComm *pcomm = pparfbssdata->GetComm ();

// Get nproc and myid

//	int nproc = pcomm->GetNproc ();
//	int myid = pcomm->GetMyid ();

// Check the size of the linear system

	int nhcells = pparfbssdata->GetNhcells ();
	int *phcells = pparfbssdata->GetHcells ();

	int ntottree = phcells[nhcells];

	if (ntottree != ntot) {
		std::cout << " ParFBSS_Array::ParFBSS_Solve_GetDistrib: incorrect size of the linear system to be solved" << std::endl;
		return -1;
	};

// Recompute list ordering

	int nlistini = pparfbssdata->GetNlistini ();
	int *plistini = pparfbssdata->GetListini ();
	int *plistord = pparfbssdata->GetListord ();
	int nppart = pparfbssdata->GetNppart ();
	int *ibegpartarr = pparfbssdata->GetIbegpartarr ();
	int *iendpartarr = pparfbssdata->GetIendpartarr ();
	int *porderpart = pparfbssdata->GetOrderpart ();
	int nlisttree = pparfbssdata->GetNlisttree ();
	int *plisttree = pparfbssdata->GetListtree ();
	int *pinvordertree = pparfbssdata->GetInvordertree ();

	CSMatrix::OrderViaIntervalsAndInverseOrder (*pcomm,
																nlistini, plistini,
																nlisttree, plisttree, pinvordertree,
																nppart, ibegpartarr, iendpartarr, porderpart,
																plistord);

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: transform the data for the solver
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_Transform (void *_PObj) {

//	const char* funcname = "ParFBSS_Solve_Transform";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CSMatrixR *pamtr = pparfbssdata->GetAmtr ();
	CSVector *prhs = pparfbssdata->GetRhs ();
	CSVector *px = pparfbssdata->GetX ();
	CSVector *psol = pparfbssdata->GetSol ();
	CSVector *pnorm = pparfbssdata->GetNorm ();
	CSVector *ppoly2 = pparfbssdata->GetPoly2 ();
//	CSlvParam *pparam = pparfbssdata->GetParams ();
	CTree *ptree = pparfbssdata->GetTree ();

// Get nproc and myid

//	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Create arrays that describe intervals

	int nhcells = pparfbssdata->GetNhcells ();
	int *phcells = pparfbssdata->GetHcells ();
	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

	int *ibegarr, *iendarr;

	ibegarr = new int [nhcells];
	if (!ibegarr) return -1;
	iendarr = new int [nhcells];
	if (!iendarr) return -1;

	int npartsloc = 0;

	int i;

	for (i=0;i<nhcells;i++) {
		if (phcells2cpu[i] == myid && phcells[i+1] > phcells[i]) {
			ibegarr[npartsloc] = phcells[i];
			iendarr[npartsloc] = phcells[i+1]-1;
			npartsloc++;
		};
	};

// Perform final transformations of the matrix

// Reorder

	int nppart = pparfbssdata->GetNppart ();
	int *ibegpartarr = pparfbssdata->GetIbegpartarr ();
	int *iendpartarr = pparfbssdata->GetIendpartarr ();

	int *plistord = pparfbssdata->GetListord ();
	int nlisttree = pparfbssdata->GetNlisttree ();
	int *plisttree = pparfbssdata->GetListtree ();
	int *pinvordertree = pparfbssdata->GetInvordertree ();

	int *porderpart = pparfbssdata->GetOrderpart ();

	CSMatrixR mtrnew = *pamtr;

	mtrnew.ChangeIndicesViaIntervals (*pcomm,
													plistord,
													nlisttree, plisttree, pinvordertree, 
													nppart, ibegpartarr, iendpartarr, porderpart);

	CSMatrixR mtrdummy;

	*pamtr = mtrdummy;

// Exchange

	CSMatrixR mtrfin;

	mtrfin = mtrnew.GetSubmatrixViaIntervals (*pcomm, npartsloc, ibegarr, iendarr);
/*
   if (true) {
      CSMatrixR amatrtot;
      int np1loc = 0;
      int ibeg1loc = 0;
      int iend1loc = -1;
      if (myid == 0) {
         np1loc = 1;
         ibeg1loc = 0;
         iend1loc = phcells[nhcells]-1;
      };
      amatrtot = mtrfin.GetSubmatrixViaIntervals (*pcomm, np1loc, &ibeg1loc, &iend1loc);
      if (myid == 0) {
         amatrtot.A2Ps (3, "ASlv.ps",0,&i);
      };
   };
*/
	mtrnew = mtrdummy;

// Exchange vector data

	int nlistini = pparfbssdata->GetNlistini ();

	double *pvect;
	double *vectini;

	vectini = new double [nlistini*5];
	if (!vectini) return -1;

	pvect = prhs->GetVect ();
	for (i=0;i<nlistini;i++) vectini[i] = pvect[i];
	pvect = px->GetVect ();
	for (i=0;i<nlistini;i++) vectini[i+nlistini] = pvect[i];
	pvect = psol->GetVect ();
	for (i=0;i<nlistini;i++) vectini[i+nlistini*2] = pvect[i];
	pvect = pnorm->GetVect ();
	for (i=0;i<nlistini;i++) vectini[i+nlistini*3] = pvect[i];
	pvect = ppoly2->GetVect ();
	for (i=0;i<nlistini;i++) vectini[i+nlistini*4] = pvect[i];

	CSVector vectdummy;

	*prhs = vectdummy;
	*px = vectdummy;

	int nlistfin;
	int *listfin;
	double *vectfin;

	CSMatrixR::GetSubvectorsViaIntervals (*pcomm, npartsloc, ibegarr, iendarr,
														5, nlistini, plistord, vectini,
														nlistfin, listfin, vectfin);

	delete [] vectini;

	CSVector *prhsnew;
	CSVector *pxnew;
	CSVector *psolnew;
	CSVector *pnormnew;
	CSVector *ppoly2new;

	prhsnew = new CSVector (nlistfin);
	if (!prhsnew) return -1;
	pxnew = new CSVector (nlistfin);
	if (!pxnew) return -1;
	psolnew = new CSVector (nlistfin);
	if (!psolnew) return -1;
	pnormnew = new CSVector (nlistfin);
	if (!pnormnew) return -1;
	ppoly2new = new CSVector (nlistfin);
	if (!ppoly2new) return -1;

	pvect = prhsnew->GetVect ();
	for (i=0;i<nlistfin;i++) pvect[i] = vectfin[i];
	pvect = pxnew->GetVect ();
	for (i=0;i<nlistfin;i++) pvect[i] = vectfin[nlistfin+i];
	pvect = psolnew->GetVect ();
	for (i=0;i<nlistfin;i++) pvect[i] = vectfin[nlistfin*2+i];
	pvect = pnormnew->GetVect ();
	for (i=0;i<nlistfin;i++) pvect[i] = vectfin[nlistfin*3+i];
	pvect = ppoly2new->GetVect ();
	for (i=0;i<nlistfin;i++) pvect[i] = vectfin[nlistfin*4+i];

	delete [] listfin;
	delete [] vectfin;

// Compute the maximal block sparsity treated by the solver as first order elements

	CSMatrix *paspblk = ptree->GetAblktree ();

	int *pialoc = paspblk->GetIa ();
	int *pjaloc = paspblk->GetJa ();

// Split the matrix into first order and second order parts

	CSMatrixR mtrfin2;

	mtrfin.SplitAccordingBlockSparsity (nhcells, phcells,
														pialoc, pjaloc, mtrfin2);

//	std::cout << " Second order = " << mtrfin2.GetNza () << std::endl;

// Create triangular submatrices

	int ntot = pparfbssdata->GetNtot ();

	CSMatrixR mtral, mtrau;

	mtrfin.CombineLU (*pcomm, nhcells, phcells, phcells2cpu,
							mtral, mtrau);

	mtral.SetNsupc (ntot);
	mtral.SetNsupr (ntot);
	mtrau.SetNsupc (ntot);
	mtrau.SetNsupr (ntot);

	mtrfin = mtrdummy;

// Create CGSMatrixR structures

	CGSMatrixR *pgmtral;
	CGSMatrixR *pgmtrau;

	pgmtral = new CGSMatrixR;
	if (!pgmtral) return -1;
	pgmtrau = new CGSMatrixR;
	if (!pgmtrau) return -1;

	pgmtral->R2GRSubmatrix2IndFastSearch (false, mtral, nhcells, phcells, 
													myid, phcells2cpu);

	mtral = mtrdummy;

	pgmtrau->R2GRSubmatrix2IndFastSearch (false, mtrau, nhcells, phcells, 
													myid, phcells2cpu);

	mtrau = mtrdummy;

// Complete computations with the second order matrix

	CGSMatrixR *pgmtra2;

	pgmtra2 = new CGSMatrixR;
	if (!pgmtra2) return -1;

	pgmtra2->R2GRSubmatrix2IndFastSearch (false, mtrfin2, nhcells, phcells, 
													myid, phcells2cpu);

	mtrfin2 = mtrdummy;

	CSMatrixR *pmtrarr2 = pgmtra2->GetMtrarr ();

	for (i=0;i<nhcells;i++) {
		if (phcells2cpu[i] == myid) {
			pmtrarr2[i].SuperSparsify2Index ();
		};
	};

// Set comm

	ptree->SetCommSchur (*pcomm);

// Store matrix and vector data in the structure

	CGSMatrixR *pmatrtmp;

	pmatrtmp = pparfbssdata->GetGmtral ();
	delete pmatrtmp;
	pparfbssdata->SetGmtral (pgmtral);

	pmatrtmp = pparfbssdata->GetGmtrau ();
	delete pmatrtmp;
	pparfbssdata->SetGmtrau (pgmtrau);

	pmatrtmp = pparfbssdata->GetGmtra2 ();
	delete pmatrtmp;
	pparfbssdata->SetGmtra2 (pgmtra2);

	CSVector *pvecttmp;

	pvecttmp = pparfbssdata->GetRhsnew ();
	delete pvecttmp;
	pparfbssdata->SetRhsnew (prhsnew);

	pvecttmp = pparfbssdata->GetXnew ();
	delete pvecttmp;
	pparfbssdata->SetXnew (pxnew);

	pvecttmp = pparfbssdata->GetSolnew ();
	delete pvecttmp;
	pparfbssdata->SetSolnew (psolnew);

	pvecttmp = pparfbssdata->GetNormnew ();
	delete pvecttmp;
	pparfbssdata->SetNormnew (pnormnew);

	pvecttmp = pparfbssdata->GetPoly2new ();
	delete pvecttmp;
	pparfbssdata->SetPoly2new (ppoly2new);

// Free work arrays

	delete [] ibegarr;
	delete [] iendarr;

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: perform linear solver iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_SolveSymm_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_SolveSymm_LinearIter";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CGSMatrixR *pgmtral = pparfbssdata->GetGmtral ();
	CGSMatrixR *pgmtrau = pparfbssdata->GetGmtrau ();
	CGSMatrixR *pgmtra2 = pparfbssdata->GetGmtra2 ();
	CSVector *prhsnew = pparfbssdata->GetRhsnew ();
	CSVector *pxnew = pparfbssdata->GetXnew ();
	CSVector *psolnew = pparfbssdata->GetSolnew ();
	CSVector *pnormnew = pparfbssdata->GetNormnew ();
//	CSVector *ppoly2new = pparfbssdata->GetPoly2new ();
	CSlvParam *pparam = pparfbssdata->GetParams ();
	CTree *ptree = pparfbssdata->GetTree ();

// Get nproc and myid

	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Create arrays that describe intervals

//	int nhcells = pparfbssdata->GetNhcells ();
//	int *phcells = pparfbssdata->GetHcells ();
//	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

// Solve

	CSlvParam paramloc = *pparam;

	std::ofstream fout;

	if (pparam->msglev > 0) {
		fout.open (pparam->path);
	};

	*psolnew = *pxnew;

	if (paramloc.sttype != -1) {
		if (paramloc.eps < 0) {
			paramloc.sttype = 0;
			paramloc.eps = -paramloc.eps;
		} else {
			paramloc.sttype = 1;
		};
	} else {
		paramloc.sttype = 2;
	};

	paramloc.ittype = 1;

// Create and register solution communicator

	int nprocsolve = 1;

	while (nprocsolve*2 <= nproc) nprocsolve *= 2;

	int nprocsolve_save = pparfbssdata->GetNprocsolve ();

	if (nprocsolve_save != nprocsolve) {
		std::cout << " ParFBSS_Array::ParFBSS_SolveSymm_LinearIter: incorrect number of processors used for linear system solution " << std::endl;
		return -1;
	};

	int *listcpu;

	listcpu = new int [nprocsolve];
	if (!listcpu) return -1;

	int i;

	for (i=0;i<nprocsolve;i++) listcpu[i] = i;

	CMPIComm comm_sol;

	comm_sol = pcomm->CreateComm (nprocsolve, listcpu);

	delete [] listcpu;

// Solve

	CGSMatrixR gmtrl, gmtru;

	if (myid < nprocsolve) {

		ptree->SetCommSchur (comm_sol);

//		psolnew->Solver2Ind (fout, *ptree,
//								*pgmtral, *pgmtrau, *pgmtral, *pgmtrau,
//								true, *pgmtra2,
//								gmtrl, gmtru,
//								*prhsnew, *pxnew, *pnormnew, 
//								paramloc);
		psolnew->SolverSymm2IndSchur (fout, *ptree,
								*pgmtral, *pgmtrau, *pgmtral,
								true, *pgmtra2,
								gmtrl,
								*prhsnew, *pxnew, *pnormnew, 
								paramloc);

		pparam->nitperf = paramloc.nitperf;
		pparam->resfin = paramloc.resfin;

	};

	int niterloc = pparam->nitperf;
	double resfinloc = pparam->resfin;
	if (myid > 0) niterloc = 0;
	if (myid > 0) resfinloc = 0.0e0;

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, 1, &niterloc, &niterloc);
	CMPIExchange::ExchangeArrayMPI (*pcomm, DOUBLEVALUE, ADD, 1, &resfinloc, &resfinloc);

	*_NIterPerf = niterloc;
	*_ResFin = resfinloc;

	*pxnew = *psolnew;

// Close output file

	if (pparam->msglev > 0) {
		fout.close();
	};

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: perform linear solver iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_Solve_LinearIter";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CGSMatrixR *pgmtral = pparfbssdata->GetGmtral ();
	CGSMatrixR *pgmtrau = pparfbssdata->GetGmtrau ();
	CGSMatrixR *pgmtra2 = pparfbssdata->GetGmtra2 ();
	CSVector *prhsnew = pparfbssdata->GetRhsnew ();
	CSVector *pxnew = pparfbssdata->GetXnew ();
	CSVector *psolnew = pparfbssdata->GetSolnew ();
	CSVector *pnormnew = pparfbssdata->GetNormnew ();
//	CSVector *ppoly2new = pparfbssdata->GetPoly2new ();
	CSlvParam *pparam = pparfbssdata->GetParams ();
	CTree *ptree = pparfbssdata->GetTree ();

// Get nproc and myid

	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Create arrays that describe intervals

//	int nhcells = pparfbssdata->GetNhcells ();
//	int *phcells = pparfbssdata->GetHcells ();
//	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

// Solve

	CSlvParam paramloc = *pparam;

	std::ofstream fout;

	if (pparam->msglev > 0) {
		fout.open (pparam->path);
	};

	*psolnew = *pxnew;

	if (paramloc.sttype != -1) {
		if (paramloc.eps < 0) {
			paramloc.sttype = 0;
			paramloc.eps = -paramloc.eps;
		} else {
			paramloc.sttype = 1;
		};
	} else {
		paramloc.sttype = 2;
	};

	paramloc.ittype = 2;
	if (pparam->ittype == 1) paramloc.ittype = 1;

// Create and register solution communicator

	int nprocsolve = 1;

	while (nprocsolve*2 <= nproc) nprocsolve *= 2;

	int nprocsolve_save = pparfbssdata->GetNprocsolve ();

	if (nprocsolve_save != nprocsolve) {
		std::cout << " ParFBSS_Array::ParFBSS_Solve_LinearIter: incorrect number of processors used for linear system solution " << std::endl;
		return -1;
	};

	int *listcpu;

	listcpu = new int [nprocsolve];
	if (!listcpu) return -1;

	int i;

	for (i=0;i<nprocsolve;i++) listcpu[i] = i;

	CMPIComm comm_sol;

	comm_sol = pcomm->CreateComm (nprocsolve, listcpu);

	delete [] listcpu;

// Solve

	CGSMatrixR gmtrl, gmtru;

	if (myid < nprocsolve) {

		ptree->SetCommSchur (comm_sol);

//		psolnew->Solver2Ind (fout, *ptree,
//								*pgmtral, *pgmtrau, *pgmtral, *pgmtrau,
//								true, *pgmtra2,
//								gmtrl, gmtru,
//								*prhsnew, *pxnew, *pnormnew, 
//								paramloc);
		psolnew->Solver2IndSchur (fout, *ptree,
								*pgmtral, *pgmtrau, *pgmtral, *pgmtrau,
								true, *pgmtra2,
								gmtrl, gmtru,
								*prhsnew, *pxnew, *pnormnew, 
								paramloc);

		pparam->nitperf = paramloc.nitperf;
		pparam->resfin = paramloc.resfin;

	};

	int niterloc = pparam->nitperf;
	double resfinloc = pparam->resfin;
	if (myid > 0) niterloc = 0;
	if (myid > 0) resfinloc = 0.0e0;

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, 1, &niterloc, &niterloc);
	CMPIExchange::ExchangeArrayMPI (*pcomm, DOUBLEVALUE, ADD, 1, &resfinloc, &resfinloc);

	*_NIterPerf = niterloc;
	*_ResFin = resfinloc;

	*pxnew = *psolnew;

// Close output file

	if (pparam->msglev > 0) {
		fout.close();
	};

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: perform nonlinear solver iterations
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_NonLinearIter (void *_PObj, int *_NIterPerf, double *_ResFin) {

//	const char* funcname = "ParFBSS_Solve_NonLinearIter";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CGSMatrixR *pgmtral = pparfbssdata->GetGmtral ();
	CGSMatrixR *pgmtrau = pparfbssdata->GetGmtrau ();
	CGSMatrixR *pgmtra2 = pparfbssdata->GetGmtra2 ();
	CSVector *prhsnew = pparfbssdata->GetRhsnew ();
	CSVector *pxnew = pparfbssdata->GetXnew ();
	CSVector *psolnew = pparfbssdata->GetSolnew ();
	CSVector *pnormnew = pparfbssdata->GetNormnew ();
	CSVector *ppoly2new = pparfbssdata->GetPoly2new ();
	CSlvParam *pparam = pparfbssdata->GetParams ();
	CTree *ptree = pparfbssdata->GetTree ();

// Get nproc and myid

	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Create arrays that describe intervals

//	int nhcells = pparfbssdata->GetNhcells ();
//	int *phcells = pparfbssdata->GetHcells ();
//	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

// Solve

	CSlvParam paramloc = *pparam;

	std::ofstream fout;

	if (pparam->msglev > 0) {
		fout.open (pparam->path);
	};

	*psolnew = *pxnew;

	if (paramloc.sttype != -1) {
		if (paramloc.eps < 0) {
			paramloc.sttype = 0;
			paramloc.eps = -paramloc.eps;
		} else {
			paramloc.sttype = 1;
		};
	} else {
		paramloc.sttype = 2;
	};

	paramloc.ittype = 2;

// Create and register solution communicator

	int nprocsolve = 1;

	while (nprocsolve*2 <= nproc) nprocsolve *= 2;

	int nprocsolve_save = pparfbssdata->GetNprocsolve ();

	if (nprocsolve_save != nprocsolve) {
		std::cout << " ParFBSS_Array::ParFBSS_Solve_Fin: incorrect number of processors used for linear system solution " << std::endl;
		return -1;
	};

	int *listcpu;

	listcpu = new int [nprocsolve];
	if (!listcpu) return -1;

	int i;

	for (i=0;i<nprocsolve;i++) listcpu[i] = i;

	CMPIComm comm_sol;

	comm_sol = pcomm->CreateComm (nprocsolve, listcpu);

	delete [] listcpu;

// Solve

	CGSMatrixR gmtrl, gmtru;

	if (myid < nprocsolve) {

		ptree->SetCommSchur (comm_sol);

		int nitersave = pparam->niter;

//		paramloc.msglev = 2;
		paramloc.niter = 100;

		CMvmR mvm (*ptree, *pgmtral, 1);
		CMvmR mvm2 (true, *ptree, *pgmtra2, 1);

		int nblkstree = ptree->GetNblkstree ();
		int *pblkstree = ptree->GetBlkstree ();
		int *pblk2cputree = ptree->GetBlk2cputree ();
		CSMatrix * pablkstr = ptree->GetAblktree ();

		CMvmSchurR mvmschur (&fout,
									nblkstree, pblkstree, pblk2cputree,
									ptree, pablkstr, 
									pgmtral, pgmtrau, 
									pgmtral, pgmtrau, 
									CSMatrixR::AbsMvmBlockColumnAl2Index , CSMatrixR::AbsMvmBlockRowAu2Index,
									CSMatrixR::AbsSolveBlockColumnL2Index, CSMatrixR::AbsSolveBlockRowU2Index,
									&paramloc);

		CParSchurMvmSlv2Index objmvslv (ptree, &mvmschur, true, &mvm2, pgmtra2);

		psolnew->NonlinearSorSchur2Index (fout, *ptree, 
												*pgmtral, *pgmtrau, 
												&objmvslv, CParSchurMvmSlv2Index::MvmA,
												&objmvslv, CParSchurMvmSlv2Index::SolveAL,
												0, 0, 0, 0,
												*prhsnew, *ppoly2new, *pxnew, *pnormnew,
												paramloc);

		*pxnew = *psolnew;

		if (paramloc.resfin > paramloc.eps)
		{

			paramloc.niter = 0;
//			paramloc.niter = 50;

			psolnew->NewtonSolveSchur2Index (fout, *ptree,
									*pgmtral, *pgmtrau, true, *pgmtra2,
									0, 0, 0, 0,
									*prhsnew, *ppoly2new, *pxnew, *pnormnew,
									paramloc);

			*pxnew = *psolnew;

		};

		pparam->niter = nitersave;
		pparam->nitperf = paramloc.nitperf;
		pparam->resfin = paramloc.resfin;

	};

	int niterloc = pparam->nitperf;
	double resfinloc = pparam->resfin;
	if (myid > 0) niterloc = 0;
	if (myid > 0) resfinloc = 0.0e0;

	CMPIExchange::ExchangeArrayMPI (*pcomm, INTEGERVALUE, ADD, 1, &niterloc, &niterloc);
	CMPIExchange::ExchangeArrayMPI (*pcomm, DOUBLEVALUE, ADD, 1, &resfinloc, &resfinloc);

	*_NIterPerf = niterloc;
	*_ResFin = resfinloc;

	*pxnew = *psolnew;

// Close output file

	if (pparam->msglev > 0) {
		fout.close();
	};

	return 0;

};

// Author: Kharchenko S.A.
// ParFBSS_Array: exchange the solution
//========================================================================================
int ParFBSS_Array::ParFBSS_Solve_ExchangeSolution (void *_PObj) {

//	const char* funcname = "ParFBSS_Solve_ExchangeSolution";

// Open the object

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CMPIComm *pcomm = pparfbssdata->GetComm ();
	CSVector *psol = pparfbssdata->GetSol ();
	CGSMatrixR *pgmtral = pparfbssdata->GetGmtral ();
	CGSMatrixR *pgmtrau = pparfbssdata->GetGmtrau ();
	CGSMatrixR *pgmtra2 = pparfbssdata->GetGmtra2 ();
	CSVector *prhsnew = pparfbssdata->GetRhsnew ();
	CSVector *pxnew = pparfbssdata->GetXnew ();
	CSVector *psolnew = pparfbssdata->GetSolnew ();
	CSVector *pnormnew = pparfbssdata->GetNormnew ();
	CSVector *ppoly2new = pparfbssdata->GetPoly2new ();
//	CSlvParam *pparam = pparfbssdata->GetParams ();
//	CTree *ptree = pparfbssdata->GetTree ();

// Get nproc and myid

//	int nproc = pcomm->GetNproc ();
	int myid = pcomm->GetMyid ();

// Create arrays that describe intervals

	int nhcells = pparfbssdata->GetNhcells ();
	int *phcells = pparfbssdata->GetHcells ();
	int *phcells2cpu = pparfbssdata->GetHcells2cpu ();

	int *ibegarr, *iendarr;

	ibegarr = new int [nhcells];
	if (!ibegarr) return -1;
	iendarr = new int [nhcells];
	if (!iendarr) return -1;

	int npartsloc = 0;

	int i;

	for (i=0;i<nhcells;i++) {
		if (phcells2cpu[i] == myid && phcells[i+1] > phcells[i]) {
			ibegarr[npartsloc] = phcells[i];
			iendarr[npartsloc] = phcells[i+1]-1;
			npartsloc++;
		};
	};

// Exchange solution data back

	int nlistini = pparfbssdata->GetNlistini ();
	int *plistord = pparfbssdata->GetListord ();

	double *pvect = psol->GetVect ();
	double *pvectnew = psolnew->GetVect ();

	CSMatrixR::GetSubvectorsViaEntryIntervals (*pcomm,
																npartsloc, ibegarr, iendarr, 
																1, pvectnew,
																nlistini, plistord, pvect);

// Delete matrix/vector data

	delete pgmtral;
	delete pgmtrau;
	delete pgmtra2;
	delete prhsnew;
	delete pxnew;
	delete psolnew;
	delete pnormnew;
	delete ppoly2new;

	pgmtral = new CGSMatrixR;
	if (!pgmtral) return -1;
	pgmtrau = new CGSMatrixR;
	if (!pgmtrau) return -1;
	pgmtra2 = new CGSMatrixR;
	if (!pgmtra2) return -1;

	pparfbssdata->SetGmtral (pgmtral);
	pparfbssdata->SetGmtrau (pgmtrau);
	pparfbssdata->SetGmtra2 (pgmtra2);

	prhsnew = new CSVector;
	if (!prhsnew) return -1;
	pxnew = new CSVector;
	if (!pxnew) return -1;
	psolnew = new CSVector;
	if (!psolnew) return -1;
	pnormnew = new CSVector;
	if (!pnormnew) return -1;
	ppoly2new = new CSVector;
	if (!ppoly2new) return -1;

	pparfbssdata->SetRhsnew (prhsnew);
	pparfbssdata->SetXnew (pxnew);
	pparfbssdata->SetSolnew (psolnew);
	pparfbssdata->SetNormnew (pnormnew);
	pparfbssdata->SetPoly2new (ppoly2new);

// Reorder solution back

	double *vectloc;

	vectloc = new double [nlistini];
	if (!vectloc) return -1;

	pvect = psol->GetVect ();

	int *piorderini = pparfbssdata->GetIorderini ();
	int iold;

	for (i=0;i<nlistini;i++) {
		iold = piorderini[i];
		vectloc[iold] = pvect[i];
	};
	for (i=0;i<nlistini;i++) pvect[i] = vectloc[i];

	delete [] vectloc;

// Free work arrays

	delete [] ibegarr;
	delete [] iendarr;

// Free arrays in the structure that will not be used

	int *plistini = pparfbssdata->GetListini ();

	delete [] piorderini;
	delete [] plistini;
	delete [] plistord;

	nlistini = 0;

	pparfbssdata->SetNlistini (nlistini);

	piorderini = new int [nlistini];
	if (!piorderini) return -1;
	plistini = new int [nlistini];
	if (!plistini) return -1;
	plistord = new int [nlistini];
	if (!plistord) return -1;

	pparfbssdata->SetIorderini (piorderini);
	pparfbssdata->SetListini (plistini);
	pparfbssdata->SetListord (plistord);

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: get solution
//========================================================================================
int ParFBSS_Array::ParFBSS_Solution (void *_PObj, double *&_Solution) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	CSVector *psol = pparfbssdata->GetSol ();

	_Solution = psol->GetVect ();

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: get optimal distribution
//========================================================================================
int ParFBSS_Array::ParFBSS_Distribution (void *_PObj, int &_NIOpt, int *&_LISTOpt) {

	CParFbssData *pparfbssdata = (CParFbssData *) _PObj;

	_NIOpt = pparfbssdata->GetNlisttree ();
	int *pinvordertree = pparfbssdata->GetInvordertree ();

	_LISTOpt = new int [_NIOpt];
	if (!_LISTOpt) return -1;

	int i;

	for (i=0;i<_NIOpt;i++) _LISTOpt[i] = pinvordertree[i];

	qsort (_LISTOpt, _NIOpt, sizeof(int), compint);

	return 0;

};

// Author: Kharchenko S.A. 
// ParFBSS_Array: destroy an object
//========================================================================================
int ParFBSS_Array::ParFBSS_Destroy (CParFbssData *_PObj) {

	delete _PObj;
	return 0;

};

// Author: Kharchenko S.A. 
// Description: Default constructor
// CParFbssData::CParFbssData()
//========================================================================================
CParFbssData::CParFbssData () {

	const char* funcname = "CParFbssData_00";

	comm = new CMPIComm;
	if (!comm) MemoryFail(funcname);

	nprocsolve = 0;
	ntot = 0;
	nlistini = 0;

	nppart = 0;
	ibegpartarr = new int [nppart];
	if (!ibegpartarr) MemoryFail(funcname);
	iendpartarr = new int [nppart];
	if (!iendpartarr) MemoryFail(funcname);

	nhcells = 0;
	nlisttree = 0;

	iorderini = new int [nlistini];
	if (!iorderini) MemoryFail(funcname);
	listini = new int [nlistini];
	if (!listini) MemoryFail(funcname);
	listord = new int [nlistini];
	if (!listord) MemoryFail(funcname);

	hcells = new int [nhcells+1];
	if (!hcells) MemoryFail(funcname);
	hcells2cpu = new int [nhcells];
	if (!hcells2cpu) MemoryFail(funcname);
	orderpart = new int [0];
	if (!orderpart) MemoryFail(funcname);
	listtree = new int [0];
	if (!listtree) MemoryFail(funcname);
	invordertree = new int [0];
	if (!invordertree) MemoryFail(funcname);

	amtr = new CSMatrixR;
	if (!amtr) MemoryFail(funcname);
	rhs = new CSVector;
	if (!rhs) MemoryFail(funcname);
	x = new CSVector;
	if (!x) MemoryFail(funcname);
	norm = new CSVector;
	if (!norm) MemoryFail(funcname);
	poly2 = new CSVector;
	if (!poly2) MemoryFail(funcname);
	sol = new CSVector;
	if (!sol) MemoryFail(funcname);
	tree = new CTree;
	if (!tree) MemoryFail(funcname);
	param = new CSlvParam;
	if (!param) MemoryFail(funcname);
	pgmtral = new CGSMatrixR;
	if (!pgmtral) MemoryFail(funcname);
	pgmtrau = new CGSMatrixR;
	if (!pgmtrau) MemoryFail(funcname);
	pgmtra2 = new CGSMatrixR;
	if (!pgmtra2) MemoryFail(funcname);
	prhsnew = new CSVector;
	if (!prhsnew) MemoryFail(funcname);
	pxnew = new CSVector;
	if (!pxnew) MemoryFail(funcname);
	pnormnew = new CSVector;
	if (!pnormnew) MemoryFail(funcname);
	ppoly2new = new CSVector;
	if (!ppoly2new) MemoryFail(funcname);
	psolnew = new CSVector;
	if (!psolnew) MemoryFail(funcname);

};

// Author: Kharchenko S.A. 
// Description: Constructor
// CParFbssData::CParFbssData()
//========================================================================================
CParFbssData::CParFbssData (void *_comm, int _ni, int _nzai) {

	const char* funcname = "CParFbssData_01";

	nprocsolve = 0;
	ntot = _ni;
	nlistini = _ni;

	nppart = 0;
	ibegpartarr = new int [nppart];
	if (!ibegpartarr) MemoryFail(funcname);
	iendpartarr = new int [nppart];
	if (!iendpartarr) MemoryFail(funcname);

	nhcells = 0;
	nlisttree = 0;

	comm = new CMPIComm;
	if (!comm) MemoryFail(funcname);

	comm->SetCommMPI (_comm);

	iorderini = new int [nlistini];
	if (!iorderini) MemoryFail(funcname);
	listini = new int [nlistini];
	if (!listini) MemoryFail(funcname);
	listord = new int [nlistini];
	if (!listord) MemoryFail(funcname);

	hcells = new int [nhcells+1];
	if (!hcells) MemoryFail(funcname);
	hcells2cpu = new int [nhcells];
	if (!hcells2cpu) MemoryFail(funcname);
	orderpart = new int [0];
	if (!orderpart) MemoryFail(funcname);
	listtree = new int [0];
	if (!listtree) MemoryFail(funcname);
	invordertree = new int [0];
	if (!invordertree) MemoryFail(funcname);

	amtr = new CSMatrixR (_ni,_nzai);
	if (!amtr) MemoryFail(funcname);
	rhs = new CSVector (_ni);
	if (!rhs) MemoryFail(funcname);
	x = new CSVector (_ni);
	if (!x) MemoryFail(funcname);
	norm = new CSVector (_ni);
	if (!norm) MemoryFail(funcname);
	poly2 = new CSVector (_ni);
	if (!poly2) MemoryFail(funcname);
	sol = new CSVector (_ni);
	if (!sol) MemoryFail(funcname);
	tree = new CTree;
	if (!tree) MemoryFail(funcname);
	param = new CSlvParam;
	if (!param) MemoryFail(funcname);
	pgmtral = new CGSMatrixR;
	if (!pgmtral) MemoryFail(funcname);
	pgmtrau = new CGSMatrixR;
	if (!pgmtrau) MemoryFail(funcname);
	pgmtra2 = new CGSMatrixR;
	if (!pgmtra2) MemoryFail(funcname);
	prhsnew = new CSVector;
	if (!prhsnew) MemoryFail(funcname);
	pxnew = new CSVector;
	if (!pxnew) MemoryFail(funcname);
	pnormnew = new CSVector;
	if (!pnormnew) MemoryFail(funcname);
	ppoly2new = new CSVector;
	if (!ppoly2new) MemoryFail(funcname);
	psolnew = new CSVector;
	if (!psolnew) MemoryFail(funcname);

	int i;

	double *pvect;

	pvect = norm->GetVect ();
	for (i=0;i<_ni;i++) pvect[i] = 1.0e0;

	pvect = poly2->GetVect ();
	for (i=0;i<_ni;i++) pvect[i] = 0.0e0;

};

// Author: Kharchenko S.A. 
// Description: Destructor
// CParFbssData::~CParFbssData()
//========================================================================================
CParFbssData::~CParFbssData () {

	nprocsolve = 0;
	ntot = 0;
	nlistini = 0;

	nppart = 0;
	delete [] ibegpartarr;
	delete [] iendpartarr;

	nhcells = 0;

	delete comm;
	delete [] iorderini;
	delete [] listini;
	delete [] listord;
	delete [] hcells;
	delete [] hcells2cpu;
	delete [] orderpart;
	delete [] listtree;
	delete [] invordertree;
	delete amtr;
	delete rhs;
	delete x;
	delete norm;
	delete poly2;
	delete sol;
	delete tree;
	delete param;
	delete pgmtral;
	delete pgmtrau;
	delete pgmtra2;
	delete prhsnew;
	delete pxnew;
	delete pnormnew;
	delete ppoly2new;
	delete psolnew;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
