//------------------------------------------------------------------------------------------------
// File: lanczosc.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "globals.h"
#include "slvparam.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "qrd.h"
#include "corr.h"
#include "tree.h"

using namespace std;


// Used global functions

void IdentityBlockGivens (int _ngiv, // Compute block Givens rotation
									dcmplx *_givarr);
void BlockGivens (int _ngiv, // Compute block Givens rotation
						dcmplx *_giv, dcmplx *_taugiv, 
						dcmplx *_aloc, dcmplx *_givarr);
void ApplyBlockGivens (int _ngiv, // Apply block Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc);
void ApplyBlockGivensH (int _ngiv, // Apply block hermitian transposed Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc);

void SolveTriangular (int _n, int _nrhs, // Solve system with triangular dense matrix
						dcmplx *_rmatr, int _ldrmatr,
						dcmplx *_rhs, int _ldrhs,
						dcmplx *_sol, int _ldsol);

// Author: Kharchenko S.A.
// Solve system with triangular dense matrix by SVD
//========================================================================================
void SolveTriangularBySVD (ofstream &_fout, int _n, int _nrhs, // Solve system with triangular dense matrix by SVD
									dcmplx *_rmatr, int _ldrmatr,
									dcmplx *_rhs, int _ldrhs,
									dcmplx *_sol, int _ldsol) {

	const char *funcname = "SolveTriangularBySVD";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ntot = _n;
	int lworkloc = 10*ntot;
	int info = 0;

	dcmplx *aloc;
	dcmplx *uloc;
	dcmplx *vloc;
	double *svloc;
	dcmplx *workloc;
	double *rworkloc;

	aloc = new dcmplx [ntot*ntot];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [ntot*ntot];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [ntot*ntot];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [ntot];
	if (!svloc) MemoryFail (funcname);
	workloc = new dcmplx [lworkloc];
	if (!workloc) MemoryFail (funcname);
	rworkloc = new double [lworkloc];
	if (!rworkloc) MemoryFail (funcname);

// Compute SVD

	int kii, kjj, irhs;

	for (kii=0;kii<ntot*ntot;kii++) aloc[kii] = czero;

	for (kjj=0;kjj<ntot;kjj++) {
		for (kii=0;kii<=kjj;kii++) {
			aloc[kjj*ntot+kii] = _rmatr[kjj*_ldrmatr+kii];
		};
	};

	zgesvd_ ("A", "A", &ntot, &ntot, 
				aloc, &ntot, svloc, uloc, &ntot, vloc, &ntot, 
				workloc, &lworkloc, rworkloc, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

	csvhyst (cout,  "SvR", ntot, svloc);
	csvhyst (_fout, "SvR", ntot, svloc);

//	OutArr (_fout," SvR = ",ntot,svloc);

// Solve

	dcmplx caux;
	double daux;

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(uloc[kjj*ntot+kii]) * _rhs[irhs*_ldrhs+kii];
			};
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	for (kjj=0;kjj<ntot;kjj++) {
		daux = svloc[kjj];
		daux = 1.0e0 / daux;
		for (irhs=0;irhs<_nrhs;irhs++) {
			caux = aloc[irhs*ntot+kjj] * daux;
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(vloc[kjj*ntot+kii]) * aloc[irhs*ntot+kii];
			};
			_sol[irhs*_ldsol+kjj] = caux;
		};
	};

// Free work arrays

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] workloc;
	delete [] rworkloc;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
//========================================================================================
void CSVectorC::BlockLanczosRight (ofstream &_fout, // Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
									const CTree &_tree, CMvmC &_mvm,
									CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
									CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
									const CSVectorC &_x, const CSVectorC &_b, 
									const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosRight";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx mcone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;

//
// Get the size of the matrix
//
	int n = _x.GetNv ();
//
// Init statistics data
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time0 = clock ();

	ops = 0.0e0;

	nza = 0;

	int ilist;

	for (ilist=0;ilist<_gmtral.nlist;ilist++) {
		int iblk = _gmtral.listb[ilist];
		nza += _gmtral.mtrarr[iblk].GetNzatot();
	};

	nza *= 2;

	nzlu = 0;
	for (ilist=0;ilist<_gmtrl.nlist;ilist++) {
		int iblk = _gmtrl.listb[ilist];
		nzlu += _gmtrl.mtrarr[iblk].GetNzatot();
	};
//
// Init guess to the solution
//
	CSVectorC sol, soli;

	sol = _x;
	soli = _x;

	soli.SetSVect(czero);
//
// Allocate temporary vector data
//
	CSVectorC bloc, r, qwork, u;

	bloc = _x;
	r = _x;
	qwork = _x;
	u = _x;
//
// Allocate all working arrays
//
	dcmplx *taup;
	dcmplx *alpha, *beta;
	dcmplx *alphaarr, *betaarr;
	dcmplx *cc, *ss, *ccbar, *ssbar;
	dcmplx *theta, *rho, *rhobar, *phi, *phibar;
	dcmplx *giv, *giv1, *giv2;
	dcmplx *tau, *cwork;
	double *sv, *rwork;
	double *resi0, *resi1, *resi, *resort0, *resort;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

//	int ldp = n;
	taup = new dcmplx [nrhs];
	if (!taup) MemoryFail (funcname);

	alpha = new dcmplx [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new dcmplx [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new dcmplx [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new dcmplx [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new dcmplx [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new dcmplx [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new dcmplx [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new dcmplx [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new dcmplx [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new dcmplx [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new dcmplx [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new dcmplx [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new dcmplx [nrhs*2];
	if (!tau) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);

	sv = new double [nrhs*2];
	if (!sv) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);

//
// DELETE IT
//
//	dcmplx *uuglob, *vvglob;
//	dcmplx *hglob;
//	int kglob = nrhs * (_param.niter+2);
//	uuglob = new dcmplx [n*kglob];
//	if (!uuglob) MemoryFail (funcname);
//	vvglob = new dcmplx [n*kglob];
//	if (!vvglob) MemoryFail (funcname);
//	hglob = new dcmplx [kglob*kglob];
//	if (!hglob) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structure
//
	CQrdC qrd (treeclp, 2, nrhs);
//
// Compute the array of right hand side norms
//
	int kii, irhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<n;kii++) {
			double auxl = _b.vect[irhs*n+kii].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*n+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs*n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi0, resi);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi0[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi0[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	dcmplx *pr = r.GetVect ();

	qrd.UpdateQrd (treeclp, n, 0, nrhs, pr, n, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						beta, nrhs);
//
// Exchange beta coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, beta, ccbar);
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
	};
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii].real ();
			aux += auxl*auxl;
			auxl = beta[irhs*nrhs+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi1);
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi1[irhs] = sqrt(aux);
		resi[irhs] = resi1[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
		};
	};
//
// Compute initial directions block
//
	int kjj, kkk;

	dcmplx *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
					pr, n, taup, 
					ccbar, nrhs, 
					pu, n);

	ops += 2*n*nrhs_2;
//	for (kii=0;kii<n*nrhs;kii++) uuglob[kii] = pu[kii];
//
// Check multiplications
//
	if (false) {
		int kii, kjj, kkk;
		dcmplx *pqwork = qwork.GetVect ();
		gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, u, qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h A Z = ",nrhs_2,ccbar);
		gmtrdummy.MvmAh (_tree, _mvm, _gmtral, _gmtrau, u, qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h A Z = ",nrhs_2,ccbar);
		_gmtru.SolveU (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h U^{-1} Z = ",nrhs_2,ccbar);
		_gmtru.SolveLConj (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h U^{-1} Z = ",nrhs_2,ccbar);
		_gmtrl.SolveL (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h L^{-1} Z = ",nrhs_2,ccbar);
		_gmtrl.SolveUConj (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h L^{-1} Z = ",nrhs_2,ccbar);
		if (n > 0) throw " Check point";
	};
//
// Multiply by the hermitian transposed preconditioned matrix
//
	CSVectorC v, w, pwork, pwork1;

	v = _b;
	w = _b;
	pwork = _b;
	pwork1 = _b;

	time2 = clock ();

	gmtrdummy.MvmAh (_tree,_mvm,_gmtral,_gmtrau,u,v);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	_gmtru.SolveLConj (_tree,_mvm,v,pwork);
//	pwork = v;

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	time2 = clock ();

	_gmtrl.SolveUConj (_tree,_mvm,pwork,w);
//	w = pwork;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += 2*nzlu*nrhs;
//
// Compute QR decomposition
//
	dcmplx *pw = w.GetVect ();

	qrd.UpdateQrd (treeclp, n, 0, nrhs, pw, n, taup);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	dcmplx *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
					pw, n, taup, 
					ccbar, nrhs, 
					pv, n);

	ops += 2*n*nrhs_2;
//	for (kii=0;kii<n*nrhs;kii++) vvglob[kii] = pv[kii];
//
// Get R part of the QR decomposition
//
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						alpha, nrhs);
//
// Exchange alpha coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, alpha, ccbar);
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
	};
//	OutArr(_fout," Ini alpha ",nrhs_2,alpha);
//
// Hermitian transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[kii*nrhs+irhs].real ();
			aux += auxl*auxl;
			auxl = alpha[kii*nrhs+irhs].imag ();
			aux += auxl*auxl;
		};
		resort0[irhs] = sqrt(aux);
		resort[irhs] = resort0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Init data for the main cycle
//
	w = v;

	for (kii=0;kii<nrhs_2;kii++) phibar[kii] = beta[kii];
	for (kii=0;kii<nrhs_2;kii++) rhobar[kii] = alpha[kii];

	int indxarr = 0;

	for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];

	for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	for (int k=0;k<_param.niter;k++) {
//
// Multiply by the preconditioned matrix
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pv = v.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Orthogonality of V ini = ",nrhs_2,ccbar);
		};
		time2 = clock ();

		_gmtrl.SolveL (_tree,_mvm,v,pwork1);
//		pwork1 = v;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		time2 = clock ();

		_gmtru.SolveU (_tree,_mvm,pwork1,pwork);
//		pwork = pwork1;

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,pwork,qwork);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Compute array of Av norms
//
		dcmplx *pqworkloc = qwork.GetVect ();
		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = 0.0e0;
			for (kii=0;kii<n;kii++) {
				double auxl = pqworkloc[irhs*n+kii].real ();
				aux += auxl*auxl;
				auxl = pqworkloc[irhs*n+kii].imag ();
				aux += auxl*auxl;
			};
			resort[irhs] = aux;
		};

		ops += 2*nrhs_2;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs, resort, resi1);
			for (irhs=0;irhs<nrhs;irhs++) resort[irhs] = resi1[irhs];
		};

		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = resort[irhs];
			resort[irhs] = sqrt(aux);
		};

		if (_tree.myid == _tree.rootcpu) {
			csvhyst (cout,  "AvNorms", nrhs, resort);
			csvhyst (_fout, "AvNorms", nrhs, resort);
		};
//
// Check qwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Current || AV ||^2 = ",nrhs_2,ccbar);
		};
//
// Update U
//
		qwork.DaxpySVectBlk (mcone,alpha,u);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		dcmplx *pqwork = qwork.GetVect ();

		qrd.UpdateQrd (treeclp, n, 0, nrhs, pqwork, n, taup);

		ops += n*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

		qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
						pqwork, n, taup, 
						ccbar, nrhs, 
						pu, n);

		ops += 2*n*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) uuglob[(k+1)*n*nrhs+kii] = pu[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

		qrd.GetRPartQrd (treeclp, 0, 
							0, nrhs, 
							beta, nrhs);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, beta, ccbar);
			for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
		};

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//		OutArr(_fout," Current beta ",nrhs_2,beta);
//
// Multiply by the hermitian transposed preconditioned matrix
//
		time2 = clock ();

		gmtrdummy.MvmAh (_tree,_mvm,_gmtral,_gmtrau,u,qwork);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

		time2 = clock ();

		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h A U for dir = ",nrhs_2,ccbar);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,u,pwork);
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h A U for dir = ",nrhs_2,ccbar);
		};
		_gmtru.SolveLConj (_tree,_mvm,qwork,pwork1);
//		pwork1 = pwork;

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			dcmplx *ppwork1 = pwork1.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Z^h U^{-h} Z for dir = ",nrhs_2,ccbar);
			_gmtru.SolveU (_tree,_mvm,qwork,pwork);
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New Z^h U^{-h} Z for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			_gmtru.SolveU (_tree,_mvm,u,pwork);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,pwork,qwork);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h U^{-h} A^h U for dir = ",nrhs_2,ccbar);
		};
		_gmtrl.SolveUConj (_tree,_mvm,pwork1,pwork);
//		pwork = pwork1;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		ops += 2*nzlu*nrhs;
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			dcmplx *pqwork = qwork.GetVect ();
			dcmplx *ppwork1 = pwork1.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Z^h L^{-h} Z for dir = ",nrhs_2,ccbar);
			_gmtrl.SolveL (_tree,_mvm,pwork1,qwork);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New Z^h L^{-h} Z for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h L^{-h} U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			_gmtrl.SolveL (_tree,_mvm,u,pwork1);
			_gmtru.SolveU (_tree,_mvm,pwork1,qwork);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,qwork,pwork1);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h L^{-h} U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New || AU^{-1} L^{-1} U||^2 for dir = ",nrhs_2,ccbar);
		};
//
// Check pwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," P^h P before update = ",nrhs_2,ccbar);
		};
//
// Compute hermitian transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(beta[kjj*nrhs+kii]);
			};
		};
//		OutArr(_fout," Hermitian transposed beta ",nrhs_2,ccbar);
//
// Update V
//
		pwork.DaxpySVectBlk (mcone,ccbar,v);

		ops += n*nrhs_2;
//
// Check pwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			_fout << " mcone = " << mcone << endl;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h U = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," V^h V = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," first of P^h P = ",nrhs_2,ccbar);
		};
//
// Compute QR decomposition of the result
//
		dcmplx *ppwork = pwork.GetVect ();

		qrd.UpdateQrd (treeclp, n, 0, nrhs, ppwork, n, taup);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

		qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
						ppwork, n, taup, 
						ccbar, nrhs, 
						pv, n);

		ops += 2*n*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) vvglob[(k+1)*n*nrhs+kii] = pv[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

		qrd.GetRPartQrd (treeclp, 0, 
							0, nrhs, 
							alpha, nrhs);

//		OutArr(_fout," Alpha from QR ",nrhs_2,alpha);
		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, alpha, ccbar);
			for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//		OutArr(_fout," Current alpha ",nrhs_2,alpha);
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = czero;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, 2*nrhs, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = czero;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = cone;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, 2*nrhs, tau, 
					giv1, 2*nrhs, giv2, 2*nrhs);

		ops += 16*nrhs_3;

		GetRPartQrd (2*nrhs, 0, 2*nrhs, giv, 2*nrhs, giv1, 2*nrhs);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+kjj]);
				ss[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+nrhs+kjj]);
				ssbar[kjj*nrhs+kii] = mcone * conj(giv2[(kii+nrhs)*nrhs*2+kjj]);
				ccbar[kjj*nrhs+kii] = conj(giv2[(kii+nrhs)*nrhs*2+nrhs+kjj]);
			};
		};
//
// Check block Givens rotation
//
/*
		if (false) {
			for (kii=0;kii<nrhs*2;kii++) {
				for (kjj=0;kjj<nrhs*2;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs*2;kkk++) {
						aux += giv2[kii*nrhs*2+kkk]*giv2[kjj*nrhs*2+kkk];
					};
					giv[kii*nrhs*2+kjj] = aux;
				};
			};
//			OutArr (fout,"Initial orth block giv = ",nrhs_2*4,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += cc[kkk*nrhs+kii]*cc[kkk*nrhs+kjj] + 
								ss[kkk*nrhs+kii]*ss[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Diag 1 orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += ccbar[kkk*nrhs+kii]*ccbar[kkk*nrhs+kjj] + 
								ssbar[kkk*nrhs+kii]*ssbar[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Diag 2 orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += -cc[kkk*nrhs+kii]*ssbar[kkk*nrhs+kjj] + 
								ss[kkk*nrhs+kii]*ccbar[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Offdia orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += cc[kkk*nrhs+kii]*rhobar[kjj*nrhs+kkk] + 
								ss[kkk*nrhs+kii]*beta[kjj*nrhs+kkk];
					};
					giv[kjj*nrhs+kii] = aux-rho[kjj*nrhs+kii];
				};
			};
//			OutArr (fout,"First QR relation = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += ccbar[kkk*nrhs+kii]*beta[kjj*nrhs+kkk] - 
								ssbar[kkk*nrhs+kii]*rhobar[kjj*nrhs+kkk];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Second QR relation = ",nrhs_2,giv);
		};
*/
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ss[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				theta[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\rho}_{i+1} = \bar{c}_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += cc[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				phi[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\phi}_{i+1} = -\tilde{s}_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = mcone*aux;
			};
		};
		for (kii=0;kii<nrhs_2;kii++) phibar[kii] = ccbar[kii];

		ops += 4*nrhs_3;
//
// Compute inverse to the current diagonal matrix
//
// Compute SVD
//
		int info = 0;

		for (kii=0;kii<nrhs_2;kii++) giv[kii] = rho[kii];
		for (kii=0;kii<nrhs_2;kii++) giv2[kii] = rho[kii];

		zgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, sv, giv, &nrhs, giv1, &nrhs, 
					cwork, &lwork, rwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		if (_tree.myid == _tree.rootcpu) {
			csvhyst (cout,  "SvDiaLancz", nrhs, sv);
			csvhyst (_fout, "SvDiaLancz", nrhs, sv);
		};
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) sv[kii] = 1.0e0 / sv[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= sv[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = conj(aux);
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
/*
		if (false) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += rho[kkk*nrhs+kii] * giv2[kjj*nrhs+kkk];
					};
					ccbar[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Inverse check = ",nrhs_2,ccbar);
		};
*/
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * phi[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update solution
//
		soli.DaxpySVectBlk (cone,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * theta[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update w
//
		pwork = v;

		pwork.DaxpySVectBlk (mcone,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			time2 = clock ();

			_gmtrl.SolveL (_tree,_mvm,soli,qwork);

			time3 = clock ();

			timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvl++;

			time2 = clock ();

			_gmtru.SolveU (_tree,_mvm,qwork,pwork1);

			time3 = clock ();

			timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvu++;

			ops += 2*nzlu*nrhs;

			pwork = sol + pwork1;

			time2 = clock ();

			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,pwork,r);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			r = bloc-r;

			ops += n*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = pr[irhs*n+kii].real ();
					aux += auxl*auxl;
					auxl = pr[irhs*n+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = aux;

			};

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															nrhs, resi, resi1);
				for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
			};

			for (irhs=0;irhs<nrhs;irhs++) {
				double aux = resi[irhs];
				resi[irhs] = sqrt(aux);
			};

			if (_tree.myid == _tree.rootcpu) {
				for (irhs=0;irhs<nrhs;irhs++) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				};
			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	time2 = clock ();

	_gmtrl.SolveL (_tree,_mvm,soli,pwork1);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	time2 = clock ();

	_gmtru.SolveU (_tree,_mvm,pwork1,pwork);

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += 2*nzlu*nrhs;

	sol.DaxpySVect (cone,pwork);
//
// Compute SVD of B_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		dcmplx *bmatr;
		double *svloc;
		dcmplx *workloc;
		double *rworkloc;

		bmatr = new dcmplx [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rworkloc = new double [lworkloc];
		if (!rworkloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = czero;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		zgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, rworkloc, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;
		delete [] rworkloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute global nza and ops

	double arr0[2], arr1[2];

	arr0[0] = (double) nza;
	arr0[1] = ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													2, arr0, arr1);
	} else {
		arr1[0] = arr0[0];
		arr1[1] = arr0[1];
	};

//	double dnza = arr1[0];
	ops = arr1[1];

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;


		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};
//	for (kii=0;kii<indxarr*nrhs;kii++) {
//		for (kjj=0;kjj<indxarr*nrhs;kjj++) {
//			dcmplx caux = czero;
//			for (kkk=0;kkk<n;kkk++) {
//				caux += uuglob[kii*n+kkk] * conj(uuglob[kjj*n+kkk]);
//			};
//			hglob[indxarr*nrhs*kii+kjj] = caux;
//		};
//	};
//	OutArr(_fout," Orth of U ",indxarr*indxarr*nrhs_2,hglob);
//	for (kii=0;kii<indxarr*nrhs;kii++) {
//		for (kjj=0;kjj<indxarr*nrhs;kjj++) {
//			dcmplx caux = czero;
//			for (kkk=0;kkk<n;kkk++) {
//				caux += vvglob[kii*n+kkk] * conj(vvglob[kjj*n+kkk]);
//			};
//			hglob[indxarr*nrhs*kii+kjj] = caux;
//		};
//	};
//	OutArr(_fout," Orth of V ",indxarr*indxarr*nrhs_2,hglob);

// Free work arrays

	delete [] taup;
	delete [] alpha;
	delete [] beta;
	delete [] alphaarr;
	delete [] betaarr;
	delete [] cc;
	delete [] ss;
	delete [] ccbar;
	delete [] ssbar;
	delete [] theta;
	delete [] rho;
	delete [] rhobar;
	delete [] phi;
	delete [] phibar;
	delete [] giv;
	delete [] giv1;
	delete [] giv2;
	delete [] tau;
	delete [] cwork;
	delete [] sv;
	delete [] rwork;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;
	delete [] resort0;
	delete [] resort;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
//========================================================================================
void CSVectorC::BlockLanczosRight (ofstream &_fout, // Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
									const CTree &_tree, 
									void *_obja, CMVMC _mvm, CMVMHC _mvmh,
									void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
									const CSVectorC &_x, const CSVectorC &_b, 
									const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosRight";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx mcone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
//
	int n = _x.GetNv ();
	int m = _b.GetNv ();
//
// Init statistics data
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time0 = clock ();

	ops = 0.0e0;

	nza = 0;
	nzlu = 0;
//
// Init guess to the solution
//
	CSVectorC sol, soli;

	sol = _x;
	soli = _x;

	soli.SetSVect(czero);
//
// Allocate temporary vector data
//
	CSVectorC bloc, r, qwork, u;

	bloc = _b;
	r = _b;
	qwork = _b;
	u = _b;
//
// Allocate all working arrays
//
	dcmplx *taup;
	dcmplx *alpha, *beta;
	dcmplx *alphaarr, *betaarr;
	dcmplx *cc, *ss, *ccbar, *ssbar;
	dcmplx *theta, *rho, *rhobar, *phi, *phibar;
	dcmplx *giv, *giv1, *giv2;
	dcmplx *tau, *cwork;
	double *sv, *rwork;
	double *resi0, *resi1, *resi, *resort0, *resort;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

//	int ldp = n;
	taup = new dcmplx [nrhs];
	if (!taup) MemoryFail (funcname);

	alpha = new dcmplx [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new dcmplx [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new dcmplx [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new dcmplx [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new dcmplx [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new dcmplx [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new dcmplx [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new dcmplx [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new dcmplx [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new dcmplx [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new dcmplx [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new dcmplx [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new dcmplx [nrhs*2];
	if (!tau) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);

	sv = new double [nrhs*2];
	if (!sv) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structure
//
	CQrdC qrd (treeclp, 2, nrhs);
//
// Compute the array of right hand side norms
//
	int kii, irhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = _b.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs*m;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi0, resi);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi0[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi0[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute QR decomposition of the residual block
//
	dcmplx *pr = r.GetVect ();

	qrd.UpdateQrd (treeclp, m, 0, nrhs, pr, m, taup);

	ops += m*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						beta, nrhs);
//
// Exchange beta coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, beta, ccbar);
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
	};
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii].real ();
			aux += auxl*auxl;
			auxl = beta[irhs*nrhs+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi1);
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi1[irhs] = sqrt(aux);
		resi[irhs] = resi1[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
		};
	};
//
// Compute initial directions block
//
	int kjj, kkk;

	dcmplx *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrd.MvmQBlk (treeclp, m, nrhs, nrhs,
					pr, m, taup, 
					ccbar, nrhs, 
					pu, m);

	ops += 2*m*nrhs_2;
//	for (kii=0;kii<n*nrhs;kii++) uuglob[kii] = pu[kii];
//
// Check multiplications
//
/*
	if (false) {
		int kii, kjj, kkk;
		dcmplx *pqwork = qwork.GetVect ();
		gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, u, qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h A Z = ",nrhs_2,ccbar);
		gmtrdummy.MvmAh (_tree, _mvm, _gmtral, _gmtrau, u, qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h A Z = ",nrhs_2,ccbar);
		_gmtru.SolveU (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h U^{-1} Z = ",nrhs_2,ccbar);
		_gmtru.SolveLConj (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h U^{-1} Z = ",nrhs_2,ccbar);
		_gmtrl.SolveL (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pqwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," first of Z^h L^{-1} Z = ",nrhs_2,ccbar);
		_gmtrl.SolveUConj (_tree,_mvm,u,qwork);
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
				};
				ccbar[nrhs*kii+kjj] = caux;
			};
		};
		OutArr(_fout," second of Z^h L^{-1} Z = ",nrhs_2,ccbar);
		if (n > 0) throw " Check point";
	};
*/
//
// Multiply by the hermitian transposed preconditioned matrix
//
	CSVectorC v, w, pwork, pwork1;

	v = _x;
	w = _x;
	pwork = _x;
	pwork1 = _x;

	time2 = clock ();

	(_mvmh) (_obja, u, v);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	(_solvel) (_objlu, v, w);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;
//
// Compute QR decomposition
//
	dcmplx *pw = w.GetVect ();

	qrd.UpdateQrd (treeclp, n, 0, nrhs, pw, n, taup);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	dcmplx *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
					pw, n, taup, 
					ccbar, nrhs, 
					pv, n);

	ops += 2*n*nrhs_2;
//	for (kii=0;kii<n*nrhs;kii++) vvglob[kii] = pv[kii];
//
// Get R part of the QR decomposition
//
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						alpha, nrhs);
//
// Exchange alpha coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, alpha, ccbar);
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
	};
//	OutArr(_fout," Ini alpha ",nrhs_2,alpha);
//
// Hermitian transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[kii*nrhs+irhs].real ();
			aux += auxl*auxl;
			auxl = alpha[kii*nrhs+irhs].imag ();
			aux += auxl*auxl;
		};
		resort0[irhs] = sqrt(aux);
		resort[irhs] = resort0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Init data for the main cycle
//
	w = v;

	for (kii=0;kii<nrhs_2;kii++) phibar[kii] = beta[kii];
	for (kii=0;kii<nrhs_2;kii++) rhobar[kii] = alpha[kii];

	int indxarr = 0;

	for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];

	for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	for (int k=0;k<_param.niter;k++) {
//
// Multiply by the preconditioned matrix
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pv = v.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Orthogonality of V ini = ",nrhs_2,ccbar);
		};

		time2 = clock ();

		(_solveu) (_objlu, v, pwork);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		(_mvm) (_obja, pwork, qwork);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + nzlu)*nrhs;
//
// Compute array of Av norms
//
		dcmplx *pqworkloc = qwork.GetVect ();
		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = 0.0e0;
			for (kii=0;kii<m;kii++) {
				double auxl = pqworkloc[irhs*m+kii].real ();
				aux += auxl*auxl;
				auxl = pqworkloc[irhs*m+kii].imag ();
				aux += auxl*auxl;
			};
			resort[irhs] = aux;
		};

		ops += 2*nrhs_2;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs, resort, resi1);
			for (irhs=0;irhs<nrhs;irhs++) resort[irhs] = resi1[irhs];
		};

		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = resort[irhs];
			resort[irhs] = sqrt(aux);
		};

		if (_tree.myid == _tree.rootcpu) {
			if (nrhs > 1) {
				csvhyst (cout,  "AvNorms", nrhs, resort);
				csvhyst (_fout, "AvNorms", nrhs, resort);
			};
		};
//
// Check qwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Current || AV ||^2 = ",nrhs_2,ccbar);
		};
//
// Update U
//
		qwork.DaxpySVectBlk (mcone,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		dcmplx *pqwork = qwork.GetVect ();

		qrd.UpdateQrd (treeclp, m, 0, nrhs, pqwork, m, taup);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

		qrd.MvmQBlk (treeclp, m, nrhs, nrhs,
						pqwork, m, taup, 
						ccbar, nrhs, 
						pu, m);

		ops += 2*m*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) uuglob[(k+1)*n*nrhs+kii] = pu[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

		qrd.GetRPartQrd (treeclp, 0, 
							0, nrhs, 
							beta, nrhs);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, beta, ccbar);
			for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
		};

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//		OutArr(_fout," Current beta ",nrhs_2,beta);
//
// Multiply by the hermitian transposed preconditioned matrix
//
		time2 = clock ();

		(_mvmh) (_obja, u, pwork1);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

/*
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h A U for dir = ",nrhs_2,ccbar);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,u,pwork);
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h A U for dir = ",nrhs_2,ccbar);
		};
*/
		time2 = clock ();

		(_solvel) (_objlu, pwork1, pwork);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		ops += nzlu*nrhs;
/*
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			dcmplx *ppwork1 = pwork1.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Z^h U^{-h} Z for dir = ",nrhs_2,ccbar);
			_gmtru.SolveU (_tree,_mvm,qwork,pwork);
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New Z^h U^{-h} Z for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			_gmtru.SolveU (_tree,_mvm,u,pwork);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,pwork,qwork);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h U^{-h} A^h U for dir = ",nrhs_2,ccbar);
		};
*/

/*
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			dcmplx *pqwork = qwork.GetVect ();
			dcmplx *ppwork1 = pwork1.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Z^h L^{-h} Z for dir = ",nrhs_2,ccbar);
			_gmtrl.SolveL (_tree,_mvm,pwork1,qwork);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New Z^h L^{-h} Z for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h L^{-h} U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			_gmtrl.SolveL (_tree,_mvm,u,pwork1);
			_gmtru.SolveU (_tree,_mvm,pwork1,qwork);
			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,qwork,pwork1);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New U^h L^{-h} U^{-h} A^h U for dir = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork1[kii*n+kkk] * conj(ppwork1[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," New || AU^{-1} L^{-1} U||^2 for dir = ",nrhs_2,ccbar);
		};
*/
//
// Check pwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," P^h P before update = ",nrhs_2,ccbar);
		};
//
// Compute hermitian transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(beta[kjj*nrhs+kii]);
			};
		};
//		OutArr(_fout," Hermitian transposed beta ",nrhs_2,ccbar);
//
// Update V
//
		pwork.DaxpySVectBlk (mcone,ccbar,v);

		ops += n*nrhs_2;
//
// Check pwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			_fout << " mcone = " << mcone << endl;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h U = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," V^h V = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," first of P^h P = ",nrhs_2,ccbar);
		};
//
// Compute QR decomposition of the result
//
		dcmplx *ppwork = pwork.GetVect ();

		qrd.UpdateQrd (treeclp, n, 0, nrhs, ppwork, n, taup);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

		qrd.MvmQBlk (treeclp, n, nrhs, nrhs,
						ppwork, n, taup, 
						ccbar, nrhs, 
						pv, n);

		ops += 2*n*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) vvglob[(k+1)*n*nrhs+kii] = pv[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

		qrd.GetRPartQrd (treeclp, 0, 
							0, nrhs, 
							alpha, nrhs);

//		OutArr(_fout," Alpha from QR ",nrhs_2,alpha);
		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, alpha, ccbar);
			for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//		OutArr(_fout," Current alpha ",nrhs_2,alpha);
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = czero;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, 2*nrhs, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = czero;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = cone;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, 2*nrhs, tau, 
					giv1, 2*nrhs, giv2, 2*nrhs);

		ops += 16*nrhs_3;

		GetRPartQrd (2*nrhs, 0, 2*nrhs, giv, 2*nrhs, giv1, 2*nrhs);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+kjj]);
				ss[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+nrhs+kjj]);
				ssbar[kjj*nrhs+kii] = mcone * conj(giv2[(kii+nrhs)*nrhs*2+kjj]);
				ccbar[kjj*nrhs+kii] = conj(giv2[(kii+nrhs)*nrhs*2+nrhs+kjj]);
			};
		};
//
// Check block Givens rotation
//
/*
		if (false) {
			for (kii=0;kii<nrhs*2;kii++) {
				for (kjj=0;kjj<nrhs*2;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs*2;kkk++) {
						aux += giv2[kii*nrhs*2+kkk]*giv2[kjj*nrhs*2+kkk];
					};
					giv[kii*nrhs*2+kjj] = aux;
				};
			};
//			OutArr (fout,"Initial orth block giv = ",nrhs_2*4,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += cc[kkk*nrhs+kii]*cc[kkk*nrhs+kjj] + 
								ss[kkk*nrhs+kii]*ss[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Diag 1 orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += ccbar[kkk*nrhs+kii]*ccbar[kkk*nrhs+kjj] + 
								ssbar[kkk*nrhs+kii]*ssbar[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Diag 2 orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += -cc[kkk*nrhs+kii]*ssbar[kkk*nrhs+kjj] + 
								ss[kkk*nrhs+kii]*ccbar[kkk*nrhs+kjj];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Offdia orth = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += cc[kkk*nrhs+kii]*rhobar[kjj*nrhs+kkk] + 
								ss[kkk*nrhs+kii]*beta[kjj*nrhs+kkk];
					};
					giv[kjj*nrhs+kii] = aux-rho[kjj*nrhs+kii];
				};
			};
//			OutArr (fout,"First QR relation = ",nrhs_2,giv);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += ccbar[kkk*nrhs+kii]*beta[kjj*nrhs+kkk] - 
								ssbar[kkk*nrhs+kii]*rhobar[kjj*nrhs+kkk];
					};
					giv[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Second QR relation = ",nrhs_2,giv);
		};
*/
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ss[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				theta[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\rho}_{i+1} = \bar{c}_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += cc[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				phi[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\phi}_{i+1} = -\tilde{s}_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = mcone*aux;
			};
		};
		for (kii=0;kii<nrhs_2;kii++) phibar[kii] = ccbar[kii];

		ops += 4*nrhs_3;
//
// Compute inverse to the current diagonal matrix
//
// Compute SVD
//
		int info = 0;

		for (kii=0;kii<nrhs_2;kii++) giv[kii] = rho[kii];
		for (kii=0;kii<nrhs_2;kii++) giv2[kii] = rho[kii];

		zgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, sv, giv, &nrhs, giv1, &nrhs, 
					cwork, &lwork, rwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		if (_tree.myid == _tree.rootcpu) {
			if (nrhs > 1) {
				csvhyst (cout,  "SvDiaLancz", nrhs, sv);
				csvhyst (_fout, "SvDiaLancz", nrhs, sv);
			};
		};
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) sv[kii] = 1.0e0 / sv[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= sv[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = conj(aux);
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
/*
		if (false) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					double aux = 0.0e0;
					for (kkk=0;kkk<nrhs;kkk++) {
						aux += rho[kkk*nrhs+kii] * giv2[kjj*nrhs+kkk];
					};
					ccbar[kjj*nrhs+kii] = aux;
				};
			};
//			OutArr (fout,"Inverse check = ",nrhs_2,ccbar);
		};
*/
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * phi[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update solution
//
		soli.DaxpySVectBlk (cone,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * theta[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update w
//
		pwork = v;

		pwork.DaxpySVectBlk (mcone,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			time2 = clock ();

			(_solveu) (_objlu, soli, pwork1);

			time3 = clock ();

			timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvu++;

			ops += nzlu*nrhs;

			pwork = sol + pwork1;

			time2 = clock ();

			(_mvm) (_obja, pwork, r);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii].real ();
					aux += auxl*auxl;
					auxl = pr[irhs*m+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = aux;

			};

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															nrhs, resi, resi1);
				for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
			};

			for (irhs=0;irhs<nrhs;irhs++) {
				double aux = resi[irhs];
				resi[irhs] = sqrt(aux);
			};

			if (_tree.myid == _tree.rootcpu) {
				for (irhs=0;irhs<nrhs;irhs++) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				};
			};

			time2 = clock ();

			(_mvmh) (_obja, r, pwork1);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			time2 = clock ();

			(_solvel) (_objlu, pwork1, pwork);

			time3 = clock ();

			timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvl++;

			ops += nzlu*nrhs;

			ppwork = pwork.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = ppwork[irhs*n+kii].real ();
					aux += auxl*auxl;
					auxl = ppwork[irhs*n+kii].imag ();
					aux += auxl*auxl;
				};
				resort[irhs] = aux;

			};

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															nrhs, resort, resi1);
				for (irhs=0;irhs<nrhs;irhs++) resort[irhs] = resi1[irhs];
			};

			for (irhs=0;irhs<nrhs;irhs++) {
				double aux = resort[irhs];
				resort[irhs] = sqrt(aux);
			};

			if (_tree.myid == _tree.rootcpu) {
				for (irhs=0;irhs<nrhs;irhs++) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resort || = " << log10(resort[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resort || = " << log10(resort[irhs]) << endl;
				};
			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps || resort[irhs]/resi0[irhs] < _param.eps) iconv++;
//				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	time2 = clock ();

	(_solveu) (_objlu, soli, pwork);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	sol.DaxpySVect (cone,pwork);
//
// Compute SVD of B_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		dcmplx *bmatr;
		double *svloc;
		dcmplx *workloc;
		double *rworkloc;

		bmatr = new dcmplx [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rworkloc = new double [lworkloc];
		if (!rworkloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = czero;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		zgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, rworkloc, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;
		delete [] rworkloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute global nza and ops

	double arr0[2], arr1[2];

	arr0[0] = (double) nza;
	arr0[1] = ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													2, arr0, arr1);
	} else {
		arr1[0] = arr0[0];
		arr1[1] = arr0[1];
	};

//	double dnza = arr1[0];
	ops = arr1[1];

// Output Lanczos statistics

//	perf = 1.0e-6 * ops / tottim;
//	if (nza != 0) ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;


		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] taup;
	delete [] alpha;
	delete [] beta;
	delete [] alphaarr;
	delete [] betaarr;
	delete [] cc;
	delete [] ss;
	delete [] ccbar;
	delete [] ssbar;
	delete [] theta;
	delete [] rho;
	delete [] rhobar;
	delete [] phi;
	delete [] phibar;
	delete [] giv;
	delete [] giv1;
	delete [] giv2;
	delete [] tau;
	delete [] cwork;
	delete [] sv;
	delete [] rwork;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;
	delete [] resort0;
	delete [] resort;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
//========================================================================================
void CSVectorC::BlockLanczosReorth (ofstream &_fout, // Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
									const CTree &_tree, CMvmC &_mvm,
									CGSMatrixCS &_gmtra,
									CGSMatrixCS &_gmtru,
									const CSVectorC &_x, const CSVectorC &_b, 
									const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosReorth";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx mcone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
//
	int n = _x.GetNv ();
	int m = _b.GetNv ();
//
// Init statistics data
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time0 = clock ();

	ops = 0.0e0;

	nza = 0;

	int ilist;

	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		nza += _gmtra.mtrarr[iblk].GetNzatot();
	};

	nzlu = 0;
	for (ilist=0;ilist<_gmtru.nlist;ilist++) {
		int iblk = _gmtru.listb[ilist];
		nzlu += _gmtru.mtrarr[iblk].GetNzatot();
	};
//
// Init guess to the solution
//
	CSVectorC sol, soli;

	sol = _x;
	soli = _x;

	soli.SetSVect(czero);
//
// Allocate temporary vector data
//
	CSVectorC bloc, r, qwork, u;

	bloc = _b;
	r = _b;
	qwork = _b;
	u = _b;
//
// Allocate all working arrays
//
	dcmplx *pp, *tauu, *tauv, *rcolmn;
	dcmplx *alpha, *beta;
	dcmplx *alphaarr, *betaarr;
	dcmplx *cc, *ss, *ccbar, *ssbar;
	dcmplx *theta, *rho, *rhobar, *phi, *phibar;
	dcmplx *giv, *giv1, *giv2;
	dcmplx *tau, *cwork;
	double *sv, *rwork;
	double *resi0, *resi1, *resi, *resort0, *resort;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

//	int ldp = n;
	int ldpp = m;

	pp = new dcmplx [nrhs*2*ldpp];
	if (!pp) MemoryFail (funcname);
	tauu = new dcmplx [nrhs*(_param.niter+2)];
	if (!tauu) MemoryFail (funcname);
	tauv = new dcmplx [nrhs*(_param.niter+2)];
	if (!tauv) MemoryFail (funcname);
	rcolmn = new dcmplx [nrhs_2*(_param.niter+2)];
	if (!rcolmn) MemoryFail (funcname);

	alpha = new dcmplx [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new dcmplx [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new dcmplx [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new dcmplx [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new dcmplx [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new dcmplx [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new dcmplx [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new dcmplx [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new dcmplx [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new dcmplx [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new dcmplx [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new dcmplx [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new dcmplx [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new dcmplx [nrhs*2];
	if (!tau) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);

	sv = new double [nrhs*2];
	if (!sv) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structures
//
	CQrdC qrdu (treeclp, _param.niter+1, nrhs);
	CQrdC qrdv (treeclp, _param.niter+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"LanczosU_",_tree.myid,"_");
		qrdu.SetupFiles (_param.ooctyp, strbuff, m);

		sprintf(strbuff, "%s%s%d%s",_param.path,"LanczosV_",_tree.myid,"_");
		qrdv.SetupFiles (_param.ooctyp, strbuff, n);

	};

//
// Compute the array of right hand side norms
//
	int kii, irhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = _b.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs*m;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi0, resi);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi0[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi0[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	_gmtra.MvmARect (sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute QR decomposition of the residual block
//
	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<m*nrhs;kii++) pp[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<m;kii++) {
			pp[irhs*m+kii] = pr[irhs*m+kii];
		};
	};

	qrdu.UpdateQrd (treeclp, m, 0, nrhs, pp, m, tauu);

	ops += m*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

	qrdu.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						beta, nrhs);
//
// Exchange beta coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, beta, ccbar);
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
	};
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii].real ();
			aux += auxl*auxl;
			auxl = beta[irhs*nrhs+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi1);
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi1[irhs] = sqrt(aux);
		resi[irhs] = resi1[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
		};
	};
//
// Compute initial directions block
//
	int kjj, kkk;

	dcmplx *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrdu.MvmQBlk (treeclp, m, nrhs, nrhs,
					pp, m, tauu, 
					ccbar, nrhs, 
					pu, m);

//	_fout << " U after Qmvm = " << u << endl;
	ops += 2*m*nrhs_2;
//	for (kii=0;kii<n*nrhs;kii++) uuglob[kii] = pu[kii];
//
// Multiply by the hermitian transposed preconditioned matrix
//
	CSVectorC v, w, pwork, pwork1;

	v = _x;
	w = _x;
	pwork = _x;
	pwork1 = _x;

	time2 = clock ();

	_gmtra.MvmAhRect (_tree,u,v);
//	_fout << " V after Ah = " << v << endl;

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	_gmtru.SolveLConj (_tree,_mvm,v,w);
//	pwork = v;
//	_fout << " W after solve = " << w << endl;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;
//
// Compute QR decomposition
//
	dcmplx *pw = w.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) pp[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			pp[irhs*n+kii] = pw[irhs*n+kii];
		};
	};

	qrdv.UpdateQrd (treeclp, n, 0, nrhs, pp, n, tauv);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	dcmplx *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = czero;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = cone;

	qrdv.MvmQBlk (treeclp, n, nrhs, nrhs,
					pp, n, tauv, 
					ccbar, nrhs, 
					pv, n);

	ops += 2*n*nrhs_2;
//	_fout << " V after mvm Qblk = " << v << endl;
//	for (kii=0;kii<n*nrhs;kii++) vvglob[kii] = pv[kii];
//
// Get R part of the QR decomposition
//
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

	qrdv.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						alpha, nrhs);
//
// Exchange alpha coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, alpha, ccbar);
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
	};
//	OutArr(_fout," Ini alpha ",nrhs_2,alpha);
//
// Hermitian transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[kii*nrhs+irhs].real ();
			aux += auxl*auxl;
			auxl = alpha[kii*nrhs+irhs].imag ();
			aux += auxl*auxl;
		};
		resort0[irhs] = sqrt(aux);
		resort[irhs] = resort0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Init data for the main cycle
//
	w = v;

	for (kii=0;kii<nrhs_2;kii++) phibar[kii] = beta[kii];
	for (kii=0;kii<nrhs_2;kii++) rhobar[kii] = alpha[kii];

	int indxarr = 0;

	for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];

	for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	int k, itr = 0;

	for (k=0;k<_param.niter;k++) {

		itr++;

//
// Multiply by the preconditioned matrix
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pv = v.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Orthogonality of V ini = ",nrhs_2,ccbar);
		};
		time2 = clock ();

		_gmtru.SolveU (_tree,_mvm,v,pwork);
//		pwork1 = v;

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		_gmtra.MvmARect (pwork,qwork);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + nzlu)*nrhs;
//
// Compute array of Av norms
//
		dcmplx *pqworkloc = qwork.GetVect ();
		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = 0.0e0;
			for (kii=0;kii<m;kii++) {
				double auxl = pqworkloc[irhs*m+kii].real ();
				aux += auxl*auxl;
				auxl = pqworkloc[irhs*m+kii].imag ();
				aux += auxl*auxl;
			};
			resort[irhs] = aux;
		};

		ops += 2*nrhs_2;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs, resort, resi1);
			for (irhs=0;irhs<nrhs;irhs++) resort[irhs] = resi1[irhs];
		};

		for (irhs=0;irhs<nrhs;irhs++) {
			double aux = resort[irhs];
			resort[irhs] = sqrt(aux);
		};

		if (_tree.myid == _tree.rootcpu) {
			if (nrhs > 1) {
				csvhyst (cout,  "AvNorms", nrhs, resort);
				csvhyst (_fout, "AvNorms", nrhs, resort);
			};
		};
//
// Check qwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *pqwork = qwork.GetVect ();
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pqwork[kii*n+kkk] * conj(pqwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," Current || AV ||^2 = ",nrhs_2,ccbar);
		};
//
// Update U
//
		qwork.DaxpySVectBlk (mcone,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		dcmplx *pqwork = qwork.GetVect ();

		for (kii=0;kii<m*nrhs;kii++) pp[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<m;kii++) {
				pp[irhs*m+kii] = pqwork[irhs*m+kii];
			};
		};

		qrdu.UpdateQrd (treeclp, m, itr*nrhs, itr*nrhs+nrhs, pp, m, tauu);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<(itr+1)*nrhs_2;kii++) rcolmn[kii] = czero;
		for (kii=0;kii<nrhs;kii++) rcolmn[kii*(itr+1)*nrhs+itr*nrhs+kii] = cone;

		qrdu.MvmQBlk (treeclp, m, (itr+1)*nrhs, nrhs,
						pp, m, tauu, 
						rcolmn, (itr+1)*nrhs, 
						pu, m);

		ops += 2*m*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) uuglob[(k+1)*n*nrhs+kii] = pu[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2*(itr+1);kii++) rcolmn[kii] = czero;

		qrdu.GetRPartQrd (treeclp, 0, 
							itr*nrhs, (itr+1)*nrhs, 
							rcolmn, nrhs);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				beta[kjj*nrhs+kii] = rcolmn[kjj*(itr+1)*nrhs+itr*nrhs+kii];
			};
		};

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, beta, ccbar);
			for (kii=0;kii<nrhs_2;kii++) beta[kii] = ccbar[kii];
		};

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//		OutArr(_fout," Current beta ",nrhs_2,beta);
//
// Multiply by the hermitian transposed preconditioned matrix
//
		time2 = clock ();

		_gmtra.MvmAhRect (_tree,u,pwork1);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

		time2 = clock ();

		_gmtru.SolveLConj (_tree,_mvm,pwork1,pwork);
//		pwork1 = pwork;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		ops += nzlu*nrhs;
//
// Compute hermitian transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(beta[kjj*nrhs+kii]);
			};
		};
//		OutArr(_fout," Hermitian transposed beta ",nrhs_2,ccbar);
//
// Update V
//
		pwork.DaxpySVectBlk (mcone,ccbar,v);

		ops += n*nrhs_2;
//
// Check pwork
//
		if (false) {
			int kii, kjj, kkk;
			dcmplx *ppwork = pwork.GetVect ();
			_fout << " mcone = " << mcone << endl;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pu[kii*n+kkk] * conj(pu[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," U^h U = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += pv[kii*n+kkk] * conj(pv[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," V^h V = ",nrhs_2,ccbar);
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<n;kkk++) {
						caux += ppwork[kii*n+kkk] * conj(ppwork[kjj*n+kkk]);
					};
					ccbar[nrhs*kii+kjj] = caux;
				};
			};
			OutArr(_fout," first of P^h P = ",nrhs_2,ccbar);
		};
//
// Compute QR decomposition of the result
//
		dcmplx *ppwork = pwork.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) pp[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				pp[irhs*n+kii] = ppwork[irhs*n+kii];
			};
		};

		qrdv.UpdateQrd (treeclp, n, itr*nrhs, itr*nrhs+nrhs, pp, n, tauv);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<(itr+1)*nrhs_2;kii++) rcolmn[kii] = czero;
		for (kii=0;kii<nrhs;kii++) rcolmn[kii*(itr+1)*nrhs+itr*nrhs+kii] = cone;

		qrdv.MvmQBlk (treeclp, n, (itr+1)*nrhs, nrhs,
						pp, n, tauv, 
						rcolmn, (itr+1)*nrhs, 
						pv, n);

		ops += 2*n*nrhs_2;
//		for (kii=0;kii<n*nrhs;kii++) vvglob[(k+1)*n*nrhs+kii] = pv[kii];
//
// Get R part of the QR decomposition
//
		for (kii=0;kii<nrhs_2*(itr+1);kii++) rcolmn[kii] = czero;

		qrdv.GetRPartQrd (treeclp, 0, 
							itr*nrhs, (itr+1)*nrhs, 
							rcolmn, nrhs);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				alpha[kjj*nrhs+kii] = rcolmn[kjj*(itr+1)*nrhs+itr*nrhs+kii];
			};
		};

//		OutArr(_fout," Alpha from QR ",nrhs_2,alpha);
		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhs_2*2, alpha, ccbar);
			for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = conj(alpha[kjj*nrhs+kii]);
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//		OutArr(_fout," Current alpha ",nrhs_2,alpha);
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = czero;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, 2*nrhs, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = czero;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = cone;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, 2*nrhs, tau, 
					giv1, 2*nrhs, giv2, 2*nrhs);

		ops += 16*nrhs_3;

		GetRPartQrd (2*nrhs, 0, 2*nrhs, giv, 2*nrhs, giv1, 2*nrhs);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+kjj]);
				ss[kjj*nrhs+kii] = conj(giv2[kii*nrhs*2+nrhs+kjj]);
				ssbar[kjj*nrhs+kii] = mcone * conj(giv2[(kii+nrhs)*nrhs*2+kjj]);
				ccbar[kjj*nrhs+kii] = conj(giv2[(kii+nrhs)*nrhs*2+nrhs+kjj]);
			};
		};
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ss[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				theta[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\rho}_{i+1} = \bar{c}_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += cc[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				phi[kjj*nrhs+kii] = aux;
			};
		};
//
// \bar{\phi}_{i+1} = -\tilde{s}_i * \bar{\phi}_{i}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = mcone*aux;
			};
		};
		for (kii=0;kii<nrhs_2;kii++) phibar[kii] = ccbar[kii];

		ops += 4*nrhs_3;
//
// Compute inverse to the current diagonal matrix
//
// Compute SVD
//
		int info = 0;

		for (kii=0;kii<nrhs_2;kii++) giv[kii] = rho[kii];
		for (kii=0;kii<nrhs_2;kii++) giv2[kii] = rho[kii];

		zgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, sv, giv, &nrhs, giv1, &nrhs, 
					cwork, &lwork, rwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		if (_tree.myid == _tree.rootcpu) {
			if (nrhs > 1) {
				csvhyst (cout,  "SvDiaLancz", nrhs, sv);
				csvhyst (_fout, "SvDiaLancz", nrhs, sv);
			};
		};
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) sv[kii] = 1.0e0 / sv[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= sv[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = conj(aux);
			};
		};
		ops += nrhs_3;
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * phi[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update solution
//
		soli.DaxpySVectBlk (cone,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx aux = czero;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += rho[kkk*nrhs+kii] * theta[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = aux;
			};
		};

		ops += nrhs_3;
//
// Update w
//
		pwork = v;

		pwork.DaxpySVectBlk (mcone,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			time2 = clock ();

			_gmtru.SolveU (_tree,_mvm,soli,pwork1);

			time3 = clock ();

			timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvu++;

			ops += nzlu*nrhs;

			pwork = sol + pwork1;

			time2 = clock ();

			_gmtra.MvmARect (pwork,r);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii].real ();
					aux += auxl*auxl;
					auxl = pr[irhs*m+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = aux;

			};

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															nrhs, resi, resi1);
				for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
			};

			for (irhs=0;irhs<nrhs;irhs++) {
				double aux = resi[irhs];
				resi[irhs] = sqrt(aux);
			};

			if (_tree.myid == _tree.rootcpu) {
				for (irhs=0;irhs<nrhs;irhs++) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				};
			};

			time2 = clock ();

			_gmtra.MvmAhRect (_tree,r,pwork1);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			time2 = clock ();

			_gmtru.SolveLConj (_tree,_mvm,pwork1,pwork);
		//	pwork = v;

			time3 = clock ();

			timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvl++;

			ops += nzlu*nrhs;

			ppwork = pwork.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = ppwork[irhs*n+kii].real ();
					aux += auxl*auxl;
					auxl = ppwork[irhs*n+kii].imag ();
					aux += auxl*auxl;
				};
				resort[irhs] = aux;

			};

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															nrhs, resort, resi1);
				for (irhs=0;irhs<nrhs;irhs++) resort[irhs] = resi1[irhs];
			};

			for (irhs=0;irhs<nrhs;irhs++) {
				double aux = resort[irhs];
				resort[irhs] = sqrt(aux);
			};

			if (_tree.myid == _tree.rootcpu) {
				for (irhs=0;irhs<nrhs;irhs++) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resort || = " << log10(resort[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resort || = " << log10(resort[irhs]) << endl;
				};
			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resort[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	time2 = clock ();

	_gmtru.SolveU (_tree,_mvm,soli,pwork);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	sol.DaxpySVect (cone,pwork);
//
// Compute SVD of B_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		dcmplx *bmatr;
		double *svloc;
		dcmplx *workloc;
		double *rworkloc;

		bmatr = new dcmplx [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rworkloc = new double [lworkloc];
		if (!rworkloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = czero;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		zgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, rworkloc, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
			throw " Error inside Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;
		delete [] rworkloc;

	};

//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrdu.CloseFiles ();
		qrdv.CloseFiles ();

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute global nza and ops

	double arr0[2], arr1[2];

	arr0[0] = (double) nza;
	arr0[1] = ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													2, arr0, arr1);
	} else {
		arr1[0] = arr0[0];
		arr1[1] = arr0[1];
	};

//	double dnza = arr1[0];
	ops = arr1[1];

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;


		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] pp;
	delete [] tauu;
	delete [] tauv;
	delete [] rcolmn;
	delete [] alpha;
	delete [] beta;
	delete [] alphaarr;
	delete [] betaarr;
	delete [] cc;
	delete [] ss;
	delete [] ccbar;
	delete [] ssbar;
	delete [] theta;
	delete [] rho;
	delete [] rhobar;
	delete [] phi;
	delete [] phibar;
	delete [] giv;
	delete [] giv1;
	delete [] giv2;
	delete [] tau;
	delete [] cwork;
	delete [] sv;
	delete [] rwork;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;
	delete [] resort0;
	delete [] resort;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
//========================================================================================
void CSVectorC::SOFLanczosNE_old (ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
										const CTree &_tree, 
										void *_obja, CMVMC _mvm, CMVMHC _mvmh,
										void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
										const CSVectorC &_x, const CSVectorC &_b) {

	const char *funcname = "SOFLanczosNE";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx mcone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
//
	int n = _x.GetNv ();
	int m = _b.GetNv ();
//
// Init statistics data
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time0 = clock ();

	ops = 0.0e0;

	nza = 0;
	nzlu = 0;
//
// Init guess to the solution
//
	CSVectorC sol, soli;

	sol = _x;
	soli = _x;

	soli.SetSVect(czero);
//
// Allocate temporary vector data
//
	CSVectorC bloc, r, qwork, u;

	bloc = _b;
	r = _b;
	qwork = _b;
	u = _b;
//
// Allocate all working arrays
//
	dcmplx *ppu, *ppv, *tauu, *tauv, *rmatr, *rcolmn;
	dcmplx *alpha, *beta;
	dcmplx *givarr, *giv, *giv1, *giv2;
	dcmplx *tau, *cwork;
	double *sv, *rwork;
	double *resi0, *resi1, *resi, *resort0, *resort;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ncolmax = _param.nblkmax;
	int nitermax = _param.niter;
	int niterlocal = _param.niterlocal;

	int lwork = 10*nrhs;

//	int ldp = n;
	int ldpu = m;
	int ldpv = n;
	int ldr = nrhs * (ncolmax+2);

	if (_param.ooctyp == 0) {
		ppu = new dcmplx [nrhs*(ncolmax+1)*ldpu];
		if (!ppu) MemoryFail (funcname);
		ppv = new dcmplx [nrhs*(ncolmax+1)*ldpv];
		if (!ppv) MemoryFail (funcname);
	} else {
		ppu = new dcmplx [nrhs*2*ldpu];
		if (!ppu) MemoryFail (funcname);
		ppv = new dcmplx [nrhs*2*ldpv];
		if (!ppv) MemoryFail (funcname);
	};
	tauu = new dcmplx [nrhs*(ncolmax+2)];
	if (!tauu) MemoryFail (funcname);
	tauv = new dcmplx [nrhs*(ncolmax+2)];
	if (!tauv) MemoryFail (funcname);
	rmatr = new dcmplx [ldr*ldr];
	if (!rmatr) MemoryFail (funcname);
	rcolmn = new dcmplx [nrhs*ldr];
	if (!rcolmn) MemoryFail (funcname);

	alpha = new dcmplx [ldr*nrhs];
	if (!alpha) MemoryFail (funcname);
	beta = new dcmplx [ldr*nrhs];
	if (!beta) MemoryFail (funcname);
	givarr = new dcmplx [nrhs_2*4*(ncolmax+2)];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new dcmplx [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new dcmplx [nrhs*2];
	if (!tau) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);

	sv = new double [nrhs*2];
	if (!sv) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);

// Init R matrix and R column

	int kii, igiv;

	for (kii=0;kii<ldr*ldr;kii++) rmatr[kii] = czero;
	for (kii=0;kii<ldr*nrhs;kii++) rcolmn[kii] = czero;
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structures
//
	CQrdC qrdu (treeclp, ncolmax+1, nrhs);
	CQrdC qrdv (treeclp, ncolmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"SOFLancNE_U_",_tree.myid,"_");
		qrdu.SetupFiles (_param.ooctyp, strbuff, m);

		sprintf(strbuff, "%s%s%d%s",_param.path,"SOFLancNE_V_",_tree.myid,"_");
		qrdv.SetupFiles (_param.ooctyp, strbuff, n);

	};

//
// Compute the array of right hand side norms
//
	int irhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = _b.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs*m;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi0, resi);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi0[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi0[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;

	CSVectorC r0 = r;
//
// Compute QR decomposition of the residual block
//
//	_fout << " Residual = " << r << endl;
	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<m*nrhs;kii++) ppu[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<m;kii++) {
			ppu[irhs*m+kii] = pr[irhs*m+kii];
		};
	};

	qrdu.UpdateQrd (treeclp, m, 0, nrhs, ppu, m, tauu);

	ops += m*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;

	qrdu.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						beta, nrhs);
//	OutArr (_fout," beta rhs ",nrhs*nrhs,beta);
//
// Exchange beta coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, beta, alpha);
		for (kii=0;kii<nrhs_2;kii++) beta[kii] = alpha[kii];
	};

// Store beta coefficient in a column

	int kjj;

	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			rcolmn[kjj*ldr+kii] = beta[kjj*nrhs+kii];
		};
	};
//	OutArr (_fout," Initial rcolmn ",nrhs,rcolmn);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii].real ();
			aux += auxl*auxl;
			auxl = beta[irhs*nrhs+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi1);
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi1[irhs] = sqrt(aux);
		resi[irhs] = resi1[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
		};
	};
//
// Compute initial directions block
//
	dcmplx *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;
	for (kii=0;kii<nrhs;kii++) beta[kii*nrhs+kii] = cone;

	qrdu.MvmQBlk (treeclp, m, nrhs, nrhs,
					ppu, m, tauu, 
					beta, nrhs, 
					pu, m);

	ops += 2*m*nrhs_2;
//
// Multiply by the hermitian transposed preconditioned matrix
//
	CSVectorC v, w, pwork, pwork1;

	v = _x;
	w = _x;
	pwork = _x;
	pwork1 = _x;

	time2 = clock ();

	(_mvmh) (_obja, u, v);

//	_fout << " V after Ah = " << v << endl;

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	(_solvel) (_objlu, v, w);
//	pwork = v;
//	_fout << " W after solve = " << w << endl;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;
//	_fout << " Pseudo Residual = " << w << endl;
//
// Compute QR decomposition
//
	dcmplx *pw = w.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) ppv[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			ppv[irhs*n+kii] = pw[irhs*n+kii];
		};
	};

	qrdv.UpdateQrd (treeclp, n, 0, nrhs, ppv, n, tauv);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	dcmplx *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) beta[kii] = czero;
	for (kii=0;kii<nrhs;kii++) beta[kii*nrhs+kii] = cone;

	qrdv.MvmQBlk (treeclp, n, nrhs, nrhs,
					ppv, n, tauv, 
					beta, nrhs, 
					pv, n);

	ops += 2*n*nrhs_2;
//	_fout << " V after mvm Qblk = " << v << endl;
//	for (kii=0;kii<n*nrhs;kii++) vvglob[kii] = pv[kii];
//
// Get R part of the QR decomposition
//
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = czero;

	qrdv.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						alpha, nrhs);
//
// Exchange alpha coefficient
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs_2*2, alpha, beta);
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = beta[kii];
	};
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[kii*nrhs+irhs].real ();
			aux += auxl*auxl;
			auxl = alpha[kii*nrhs+irhs].imag ();
			aux += auxl*auxl;
		};
		resort0[irhs] = sqrt(aux);
		resort[irhs] = resort0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Reso || = " << log10(resort0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Init data for the main cycle
//
	w = v;
//
// Main iterative cycle
//
//	int ichk = _param.ichk;

	int ibs;
	int iterglob = 0, iterlocal;
	int niterlocalmax;
	int iconv;

	int ibegblk = 0;
	int iendblk = -1;
	int inextblk = 0;

	bool svconv = false;
	double svmin_r, svmax_r;

	while (iterglob<nitermax) {

// Implement local iterative cycle

		niterlocalmax = niterlocal;
		if (ncolmax-iendblk < niterlocalmax) niterlocalmax = ncolmax-iendblk;

		iterlocal = 0;
		while (iterlocal<niterlocalmax) {

// Check matrix relations
/*
			if (false) {

				if (nrhs>1) throw " CSVectorC::SOFLanczosNE: the case nrhs > 1 is not implemented yet ";

				dcmplx *yloc = new dcmplx [n*inextblk];
				dcmplx *ayloc = new dcmplx [m*(inextblk+1)];
				dcmplx *wloc = new dcmplx [m*(inextblk+1)];
				dcmplx *rrloc = new dcmplx [(inextblk+1)*(inextblk+1)];

				int i;

				for (i=0;i<inextblk*inextblk;i++) rrloc[i] = czero;
				for (i=0;i<inextblk;i++) rrloc[i*inextblk+i] = cone;

				if (inextblk != 0) {

					for (i=0;i<inextblk;i++) {
						qrdv.MvmQBlk (treeclp, n, inextblk, 1,
											ppv, n, tauv, 
											rrloc+i*inextblk, inextblk, 
											yloc+i*n, n);
					};

				};

				dcmplx *pr0 = r0.GetVect ();
				dcmplx *ppwork = pwork.GetVect ();
				dcmplx *pqwork = qwork.GetVect ();

				int j;

				for (j=0;j<m;j++) {
					ayloc[j] = pr0[j];
				};
				for (i=0;i<inextblk;i++) {
					for (j=0;j<n;j++) {
						ppwork[j] = yloc[i*n+j];
					};
					(_solveu) (_objlu, pwork, pwork1);
					(_mvm) (_obja, pwork1, qwork);
					for (j=0;j<m;j++) {
						ayloc[(i+1)*m+j] = pqwork[j];
					};
				};

				for (i=0;i<(inextblk+1)*(inextblk+1);i++) rrloc[i] = czero;

				qrdu.GetRPartQrd (treeclp, 0, 0, inextblk+1,
										rrloc, inextblk+1);

				OutArr (_fout," rrloc ",(inextblk+1)*(inextblk+1),rrloc);

				for (i=0;i<inextblk+1;i++) {
					qrdu.MvmQBlk (treeclp, m, inextblk+1, 1,
										ppu, m, tauu, 
										rrloc+i*(inextblk+1), inextblk+1, 
										wloc+i*m, m);
				};

				OutArr (_fout," wloc ",(inextblk+1)*m,wloc);

				int ind;
				double diff;
				dcmplx caux;
				for (i=0;i<inextblk+1;i++) {
					diff = 0.0e0;
					for (j=0;j<m;j++) {
						ind = i*m+j;
						caux = wloc[ind]-ayloc[ind];
						diff += caux.x*caux.x + caux.y*caux.y;
					};
					diff = sqrt(diff);
					cout << " Icol in check = " << i << " Diff = " << diff << endl;
				};

				delete [] yloc;
				delete [] ayloc;
				delete [] wloc;
				delete [] rrloc;
			};
*/
// Next iteration

			iterlocal++;
			iterglob++;

			iendblk++;
			inextblk++;
//
// Multiply by the preconditioned matrix
//
			time2 = clock ();

//			_fout << " y dir = " << v << endl;
			(_solveu) (_objlu, v, pwork);

			time3 = clock ();

			timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvu++;

			time2 = clock ();

			(_mvm) (_obja, pwork, qwork);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += (nza + nzlu)*nrhs;
//			_fout << " Ay = " << qwork << endl;
//
// Compute QR decomposition of the result
//
			dcmplx *pqwork = qwork.GetVect ();

			ibs = 0;
			if (_param.ooctyp == 0) ibs = inextblk*nrhs*m;

			for (kii=0;kii<m*nrhs;kii++) ppu[ibs+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<m;kii++) {
					ppu[ibs+irhs*m+kii] = pqwork[irhs*m+kii];
				};
			};

			qrdu.UpdateQrd (treeclp, m, inextblk*nrhs, inextblk*nrhs+nrhs, ppu, m, tauu);

			ops += m*nrhs_2;
//
// Take current block column of H_k
//
			int ibeg = inextblk*nrhs;
			int iend = inextblk*nrhs+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) alpha[kii] = czero;

			qrdu.GetRPartQrd (treeclp, 0, ibeg, iend,
									alpha, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (igiv=0;igiv<iendblk;igiv++) {

				ApplyBlockGivens (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = alpha[irhs*iend+iendblk*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
								giv, tau, giv1, givarr+iendblk*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
									iendblk, givarr+iendblk*nrhs_2*4,
									nrhs, alpha, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<inextblk*nrhs;kii++) {
					rmatr[(iendblk*nrhs+irhs)*ldr+kii] = alpha[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
									iendblk, givarr+iendblk*nrhs_2*4,
									nrhs, rcolmn, ldr, giv);

			ops += 8*nrhs_3;
//
// Compute residual norm
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = rcolmn[irhs*ldr+inextblk*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = rcolmn[irhs*ldr+inextblk*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

			};
//
// Compute current U
//
			pu = u.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) alpha[kii] = czero;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					alpha[kjj*(inextblk+1)*nrhs+inextblk*nrhs+kii] = rcolmn[kjj*ldr+inextblk*nrhs+kii];
				};
			};

			for (igiv=iendblk;igiv>=0;igiv--) {

				ApplyBlockGivensH (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, (inextblk+1)*nrhs, giv);

				ops += 8*nrhs_3;

			};

			qrdu.MvmQBlk (treeclp, m, (inextblk+1)*nrhs, nrhs,
								ppu, m, tauu, 
								alpha, (inextblk+1)*nrhs, 
								pu, m);

			ops += 2*m*(inextblk+1)*nrhs*nrhs;
//
// Multiply by the hermitian transposed preconditioned matrix
//
			time2 = clock ();

			(_mvmh) (_obja, u, pwork1);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			time2 = clock ();

			(_solvel) (_objlu, pwork1, pwork);

			time3 = clock ();

			timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvl++;

			ops += nzlu*nrhs;
//
// Compute QR decomposition of the result
//
			dcmplx *ppwork = pwork.GetVect ();

			ibs = 0;
			if (_param.ooctyp == 0) ibs = inextblk*nrhs*n;

			for (kii=0;kii<n*nrhs;kii++) ppv[ibs+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					ppv[ibs+irhs*n+kii] = ppwork[irhs*n+kii];
				};
			};

			qrdv.UpdateQrd (treeclp, n, inextblk*nrhs, inextblk*nrhs+nrhs, ppv, n, tauv);

			ops += n*(inextblk*nrhs)*nrhs;
//
// Get residual orhogonality
//
			ibeg = inextblk*nrhs;
			iend = inextblk*nrhs+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) alpha[kii] = czero;

			qrdv.GetRPartQrd (treeclp, 0, ibeg, iend,
									alpha, iend);
//
// Compute pseudo residual norm
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<iend;kii++) {
					double auxl = alpha[irhs*iend+kii].real ();
					aux += auxl*auxl;
					auxl = alpha[irhs*iend+kii].imag ();
					aux += auxl*auxl;
				};
				resort[irhs] = sqrt(aux);

			};
//
// Compute current V
//
			pv = v.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) beta[kii] = czero;
			for (kii=0;kii<nrhs;kii++) beta[kii*(inextblk+1)*nrhs+inextblk*nrhs+kii] = cone;

			qrdv.MvmQBlk (treeclp, n, (inextblk+1)*nrhs, nrhs,
							ppv, n, tauv, 
							beta, (inextblk+1)*nrhs, 
							pv, n);

			ops += 2*n*(inextblk*nrhs)*nrhs;
//
// Check convergence
//
			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
//				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
				if (resi[irhs]/resi0[irhs] < _param.eps || (svconv && resort[irhs] / resi0[irhs] < _param.eps)) iconv++;
				if (iterglob%_param.ichk == 0) {
					cout  << " Irh= " << irhs << " Kin= " << nrhs*ibegblk << " Klc= " << nrhs*inextblk << " Itr= " << iterglob << " Log10 ||R|| = " << log10(resi[irhs]) << " Log10 ||ROrt|| = " << log10(resort[irhs]) << endl;
					_fout << " Irh= " << irhs << " Kin= " << nrhs*ibegblk << " Klc= " << nrhs*inextblk << " Itr= " << iterglob << " Log10 ||R|| = " << log10(resi[irhs]) << " Log10 ||ROrt|| = " << log10(resort[irhs]) << endl;
				};
			};

			if (iconv == nrhs) goto exit;

			if (iterglob >= nitermax) goto exit;

		};

// Perform filtering of the matrix relations

//		_fout << " Current residual vector = " << u << endl;
		SOFLanNEFilterRelations_old (_fout, 
											_tree, treeclp, 
											_obja, _mvm, _mvmh,
											_objlu, _solvel, _solveu,
											_param,
											ibegblk, iendblk, inextblk, 
											m, n, qrdu, qrdv, ppu, ppv, tauu, tauv,
											ldr, rmatr, rcolmn, u, w, pwork,
											givarr, giv, giv1,
											svmin_r, svmax_r);

		svconv = false;
		if (svmin_r > _param.xmin*0.9e0 && svmax_r < _param.xmax*1.1e0) svconv = true;

		cout << " svmin_r = " << svmin_r << " svmax_r = " << svmax_r << " svconv = " << svconv << endl;

//
// Compute QR decomposition of the result
//
		pv = v.GetVect ();

		ibs = 0;
		if (_param.ooctyp == 0) ibs = inextblk*nrhs*n;

		for (kii=0;kii<n*nrhs;kii++) ppv[ibs+kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				ppv[ibs+irhs*n+kii] = pv[irhs*n+kii];
			};
		};

		qrdv.UpdateQrd (treeclp, n, inextblk*nrhs, inextblk*nrhs+nrhs, ppv, n, tauv);

		ops += n*(inextblk*nrhs)*nrhs;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) beta[kii] = czero;
		for (kii=0;kii<nrhs;kii++) beta[kii*(inextblk+1)*nrhs+inextblk*nrhs+kii] = cone;

		qrdv.MvmQBlk (treeclp, n, (inextblk+1)*nrhs, nrhs,
						ppv, n, tauv, 
						beta, (inextblk+1)*nrhs, 
						pv, n);

		ops += 2*n*(inextblk*nrhs)*nrhs;
//
// Compute current U
//
		if (false) {

			pu = u.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) alpha[kii] = czero;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					alpha[kjj*(inextblk+1)*nrhs+inextblk*nrhs+kii] = rcolmn[kjj*ldr+inextblk*nrhs+kii];
				};
			};

			for (igiv=iendblk;igiv>=0;igiv--) {

				ApplyBlockGivensH (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, (inextblk+1)*nrhs, giv);

				ops += 8*nrhs_3;

			};

			qrdu.MvmQBlk (treeclp, m, (inextblk+1)*nrhs, nrhs,
								ppu, m, tauu, 
								alpha, (inextblk+1)*nrhs, 
								pu, m);

			_fout << " Recomputed current residual vector = " << u << endl;
//			goto exit;
		};

	};
//
// Update solution vector
//
exit:;
//
// Solve triangular system
//
	int ncolp = (inextblk)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) alpha[kii] = czero;

	SolveTriangularBySVD (_fout, ncolp, nrhs,
									rmatr, ldr,
									rcolmn, ldr,
									alpha, ldeloc);

// Multiply V by solution vector

	pv = v.GetVect ();

	qrdv.MvmQBlk (treeclp, n, inextblk*nrhs, nrhs,
						ppv, n, tauv, 
						alpha, ldeloc, 
						pv, n);

	ops += 2*n*(inextblk)*nrhs*nrhs;

// Update solution

	time2 = clock ();

	(_solveu) (_objlu, v, pwork);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	sol.DaxpySVect (cone,pwork);
//
// Compute the final residual and its norm
//
	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = r.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = r.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		};
	};
// Check orthogonality of the final residual

	time2 = clock ();

	(_mvmh) (_obja, r, v);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	(_solvel) (_objlu, v, w);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<n;kii++) {
			double auxl = w.vect[irhs*n+kii].real ();
			aux += auxl*auxl;
			auxl = w.vect[irhs*n+kii].imag ();
			aux += auxl*auxl;
		};
		resort[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || ResOrt || = " << log10(resort[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || ResOrt || = " << log10(resort[irhs]) << endl;
		};
	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrdu.CloseFiles ();
		qrdv.CloseFiles ();

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute global nza and ops

	double arr0[2], arr1[2];

	arr0[0] = (double) nza;
	arr0[1] = ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													2, arr0, arr1);
	} else {
		arr1[0] = arr0[0];
		arr1[1] = arr0[1];
	};

//	double dnza = arr1[0];
	ops = arr1[1];

// Output Lanczos statistics

//	perf = 1.0e-6 * ops / tottim;
//	ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;


		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] ppu;
	delete [] ppv;
	delete [] tauu;
	delete [] tauv;
	delete [] rmatr;
	delete [] rcolmn;
	delete [] alpha;
	delete [] beta;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] giv2;
	delete [] tau;
	delete [] cwork;
	delete [] sv;
	delete [] rwork;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;
	delete [] resort0;
	delete [] resort;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
//========================================================================================
void CSVectorC::SOFLanNEFilterRelations_old (ofstream &_fout, // Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
														const CTree &_tree, const CTree &_treeclp, 
														void *_obja, CMVMC _mvm, CMVMHC _mvmh,
														void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
														const CSlvParam &_param,
														int &_ibegblk, int &_iendblk, int &_inextblk, 
														int _m, int _n, 
														CQrdC &_qrdu, CQrdC &_qrdv, dcmplx *_ppu, dcmplx *_ppv, 
														dcmplx *_tauu, dcmplx *_tauv,
														int _ldr, dcmplx *_rmatr, dcmplx *_rcolmn,
														CSVectorC &_u, CSVectorC &_v, CSVectorC &_pwork,
														dcmplx *_givarr, dcmplx *_giv, dcmplx *_giv1,
														double &_svmin_r, double &_svmax_r) { 

	const char *funcname = "SOFLanNEFilterRelations";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);

// Take filtering parameters

	double xmin = _param.xmin;
	double xmax = _param.xmax;

// Allocate the memory

	int nrhs_2 = nrhs*nrhs;
	int nloc = (_iendblk-_ibegblk+1)*nrhs;
	int ntot = (_inextblk)*nrhs;
	int lworkloc = 10*nloc;
	int info = 0;

	dcmplx *aloc;
	dcmplx *uloc;
	dcmplx *vloc;
	dcmplx *tauloc;
	double *svloc;
	dcmplx *workloc;
	double *rworkloc;
	dcmplx *solloc;
	dcmplx *uextloc;

	aloc = new dcmplx [nloc*nloc];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [nloc*nloc];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [nloc*nloc];
	if (!vloc) MemoryFail (funcname);
	tauloc = new dcmplx [nloc];
	if (!tauloc) MemoryFail (funcname);
	svloc = new double [nloc];
	if (!svloc) MemoryFail (funcname);
	workloc = new dcmplx [lworkloc];
	if (!workloc) MemoryFail (funcname);
	rworkloc = new double [lworkloc];
	if (!rworkloc) MemoryFail (funcname);
	solloc = new dcmplx [nrhs*nloc];
	if (!solloc) MemoryFail (funcname);
	uextloc = new dcmplx [nloc*ntot];
	if (!uextloc) MemoryFail (funcname);

// Compute new local vectors (singular and solution)

// Take a part of R matrix and right hand side

	int kii, kjj, irhs;
	int ishift = _ibegblk*nrhs;

	for (kii=0;kii<nloc*nloc;kii++) aloc[kii] = czero;

	for (kjj=0;kjj<nloc;kjj++) {
		for (kii=0;kii<=kjj;kii++) {
			aloc[kjj*nloc+kii] = _rmatr[(ishift+kjj)*_ldr+ishift+kii];
		};
	};

	for (kjj=0;kjj<nrhs;kjj++) {
		for (kii=0;kii<nloc;kii++) {
			solloc[kjj*nloc+kii] = _rcolmn[kjj*_ldr+ishift+kii];
		};
	};

// Modify rcolmn and rmatrix back

	int igiv;

	for (igiv=_iendblk;igiv>=_ibegblk;igiv--) {

		ApplyBlockGivensH (nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								nrhs, _rcolmn, _ldr, _giv);

	};

	for (kjj=0;kjj<nrhs;kjj++) {
		for (kii=(_ibegblk+1)*nrhs;kii<(_inextblk+1)*nrhs;kii++) {
			_rcolmn[kjj*_ldr+kii] = czero;
		};
	};

	for (kjj=ishift;kjj<ntot;kjj++) {
		for (kii=0;kii<_ldr;kii++) {
			_rmatr[kjj*_ldr+kii] = czero;
		};
	};

// Compute SVD

	zgesvd_ ("A", "A", &nloc, &nloc, 
				aloc, &nloc, svloc, uloc, &nloc, vloc, &nloc, 
				workloc, &lworkloc, rworkloc, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

	csvhyst (cout,  "SvRSubm", nloc, svloc);
	csvhyst (_fout, "SvRSubm", nloc, svloc);

	_svmax_r = svloc[0];
	_svmin_r = svloc[nloc-1];

//	OutArr (_fout," SvRSubm = ",nloc,svloc);

// Solve

	dcmplx caux;
	double daux;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kjj=0;kjj<nloc;kjj++) {
			caux = czero;
			for (kii=0;kii<nloc;kii++) {
				caux += conj(uloc[kjj*nloc+kii]) * solloc[irhs*nloc+kii];
			};
			aloc[irhs*nloc+kjj] = caux;
		};
	};

	for (kjj=0;kjj<nloc;kjj++) {
		daux = svloc[kjj];
		daux = 1.0e0 / daux;
		for (irhs=0;irhs<nrhs;irhs++) {
			caux = aloc[irhs*nloc+kjj] * daux;
			aloc[irhs*nloc+kjj] = caux;
		};
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kjj=0;kjj<nloc;kjj++) {
			caux = czero;
			for (kii=0;kii<nloc;kii++) {
				caux += conj(vloc[kjj*nloc+kii]) * aloc[irhs*nloc+kii];
			};
			solloc[irhs*nloc+kjj] = caux;
		};
	};

// Filter singular values

// Large:

	int nfiltr = 0;

	for (kii=0;kii<nloc;kii++) {
		if (svloc[kii] > xmax) {
			for (kjj=0;kjj<nloc;kjj++) {
				aloc[nfiltr*nloc+kjj] = conj(vloc[kjj*nloc+kii]);
			};
			nfiltr++;
		};
	};

// Small:

	for (kii=nloc-1;kii>=0;kii--) {
		if (svloc[kii] < xmin) {
			for (kjj=0;kjj<nloc;kjj++) {
				aloc[nfiltr*nloc+kjj] = conj(vloc[kjj*nloc+kii]);
			};
			nfiltr++;
		};
	};

// Additional

	if (nfiltr%nrhs != 0) throw " CSVectorC::SOFLanNEFilterRelations: the case nrhs != 1 not implemented ";

// Add solution data

	if (nfiltr > nloc-nrhs) nfiltr = nloc-nrhs;

	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			aloc[nfiltr*nloc+kjj] = solloc[kii*nloc+kjj];
		};
		nfiltr++;
	};

// Orthonormalize the result

	QrdBlk (nloc, nfiltr, aloc, nloc, tauloc);

	for (kii=0;kii<nfiltr*nfiltr;kii++) vloc[kii] = czero;

	for (kii=0;kii<nfiltr;kii++) vloc[kii*nfiltr+kii] = cone;

	MvmQBlk (nloc, nfiltr, nfiltr,
				aloc, nloc, tauloc, 
				vloc, nfiltr, uloc, nloc);

// Compute new direction vectors

	for (kii=0;kii<nfiltr*ntot;kii++) uextloc[kii] = czero;

	for (kii=0;kii<nfiltr;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			uextloc[kii*ntot+ishift+kjj] = uloc[kii*nloc+kjj];
		};
	};

	int nparts = nfiltr / nrhs;
	int iendblkloc = _ibegblk+nparts-1;
	int inextblkloc = iendblkloc+1;

	dcmplx *vdir;

	vdir = new dcmplx [_m*nfiltr];
	if (!vdir) MemoryFail (funcname);

	int iblk;

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		_qrdv.MvmQBlk (_treeclp, _n, _inextblk*nrhs, nrhs,
						_ppv, _n, _tauv, 
						uextloc+(iblk-_ibegblk)*nrhs*ntot, ntot, 
						vdir+(iblk-_ibegblk)*nrhs*_n, _n);

	};

// Update V part

	int ibs;

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		ibs = 0;
		if (_param.ooctyp == 0) ibs = iblk*nrhs*_n;

		for (kii=0;kii<_n*nrhs;kii++) _ppv[ibs+kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<_n;kii++) {
				_ppv[ibs+irhs*_n+kii] = vdir[((iblk-_ibegblk)*nrhs+irhs)*_n+kii];
			};
		};

		_qrdv.UpdateQrd (_treeclp, _n, iblk*nrhs, iblk*nrhs+nrhs, _ppv, _n, _tauv);

	};
//
// Compute directions explicitely
//
	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		for (kii=0;kii<nrhs*ntot;kii++) uextloc[kii] = czero;

		for (kii=0;kii<nrhs;kii++) {
			uextloc[kii*ntot+ishift+(iblk-_ibegblk)*nrhs+kii] = cone;
		};

		_qrdv.MvmQBlk (_treeclp, _n, inextblkloc*nrhs, nrhs,
						_ppv, _n, _tauv, 
						uextloc, ntot, 
						vdir+(iblk-_ibegblk)*nrhs*_n, _n);

	};

// Multiply by the preconditioned matrix and update QR decomposition

	dcmplx *pv = _v.GetVect ();

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		int ibeg = (iblk+1)*nrhs;
		int iend = (iblk+1)*nrhs+nrhs;

		for (kii=0;kii<nrhs*_n;kii++) {
			pv[kii] = vdir[(iblk-_ibegblk)*nrhs*_n+kii];
		};
//
// Multiply by the preconditioned matrix
//
		(_solveu) (_objlu, _v, _pwork);

		(_mvm) (_obja, _pwork, _u);
//
// Compute QR decomposition of the result
//
		dcmplx *pu = _u.GetVect ();

		ibs = 0;
		if (_param.ooctyp == 0) ibs = ibeg*_m;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<_m;kii++) {
				_ppu[ibs+irhs*_m+kii] = pu[irhs*_m+kii];
			};
		};

		_qrdu.UpdateQrd (_treeclp, _m, ibeg, iend, 
								_ppu, _m, _tauu);
//
// Take current block column of H_k
//
		for (kii=0;kii<iend*nrhs;kii++) uextloc[kii] = czero;

		_qrdu.GetRPartQrd (_treeclp, 0, ibeg, iend,
									uextloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (igiv=0;igiv<iblk;igiv++) {

			ApplyBlockGivens (nrhs, 
									igiv, _givarr+igiv*nrhs_2*4,
									nrhs, uextloc, iend, _giv);

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) _giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				_giv[irhs*nrhs*2+kii] = uextloc[irhs*iend+iblk*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
							_giv, tauloc, _giv1, _givarr+iblk*nrhs_2*4);
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
								iblk, _givarr+iblk*nrhs_2*4,
								nrhs, uextloc, iend, _giv);
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(iblk+1)*nrhs;kii++) {
				_rmatr[(iblk*nrhs+irhs)*_ldr+kii] = uextloc[irhs*iend+kii];
			};
		};

	};

// Update matrix relations (modify rcolmn forward)

	for (igiv=_ibegblk;igiv<=iendblkloc;igiv++) {
		ApplyBlockGivens (nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								nrhs, _rcolmn, _ldr, _giv);
	};

// Reassign indices that control iterative scheme

	_ibegblk = iendblkloc;
	_iendblk = iendblkloc;
	_inextblk = inextblkloc;

// Free work memory

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] tauloc;
	delete [] svloc;
	delete [] workloc;
	delete [] rworkloc;
	delete [] solloc;
	delete [] uextloc;
	delete [] vdir;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
//========================================================================================
void CSVectorC::SOFLanczosNE (ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
										const CTree &_tree, CCorrectorC &_corr,
										void *_obja, CMVMC _mvm, CMVMHC _mvmh,
										void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
										const CSVectorC &_x, const CSVectorC &_b) {

	const char *funcname = "SOFLanczosNE";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx mcone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
//
	int n = _x.GetNv ();
	int m = _b.GetNv ();
//
// Init statistics data
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time0 = clock ();

	ops = 0.0e0;

	nza = 0;
	nzlu = 0;
//
// Init guess to the solution
//
	CSVectorC sol, soli;

	sol = _x;
	soli = _x;

	soli.SetSVect(czero);
//
// Allocate temporary vector data
//
	CSVectorC bloc, r, qwork, u;

	bloc = _b;
	r = _b;
	qwork = _b;
	u = _b;

	CSVectorC v, w, pwork, pwork1;

	v = _x;
	w = _x;
	pwork = _x;
	pwork1 = _x;

//
// Allocate all working arrays
//
	dcmplx *ppu, *ppv, *tauu, *tauv, *rmatr, *rcolmn;
	dcmplx *alpha, *beta;
	dcmplx *givarr, *giv, *giv1, *giv2;
	dcmplx *tau, *cwork;
	double *sv, *rwork;
	double *resi0, *resi1, *resi, *resort0, *resort;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ncolmax = _param.nblkmax;
	int nitermax = _param.niter;
	int niterlocal = _param.niterlocal;

	int lwork = 10*nrhs;

//	int ldp = n;
	int ldpu = m;
	int ldpv = n;
	int ldr = nrhs * (ncolmax+2);

	if (_param.ooctyp == 0) {
		ppu = new dcmplx [nrhs*(ncolmax+1)*ldpu];
		if (!ppu) MemoryFail (funcname);
		ppv = new dcmplx [nrhs*(ncolmax+1)*ldpv];
		if (!ppv) MemoryFail (funcname);
	} else {
		ppu = new dcmplx [nrhs*2*ldpu];
		if (!ppu) MemoryFail (funcname);
		ppv = new dcmplx [nrhs*2*ldpv];
		if (!ppv) MemoryFail (funcname);
	};
	tauu = new dcmplx [nrhs*(ncolmax+2)];
	if (!tauu) MemoryFail (funcname);
	tauv = new dcmplx [nrhs*(ncolmax+2)];
	if (!tauv) MemoryFail (funcname);
	rmatr = new dcmplx [ldr*ldr];
	if (!rmatr) MemoryFail (funcname);
	rcolmn = new dcmplx [nrhs*ldr];
	if (!rcolmn) MemoryFail (funcname);

	alpha = new dcmplx [ldr*nrhs];
	if (!alpha) MemoryFail (funcname);
	beta = new dcmplx [ldr*nrhs];
	if (!beta) MemoryFail (funcname);
	givarr = new dcmplx [nrhs_2*4*(ncolmax+2)];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new dcmplx [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new dcmplx [nrhs*2];
	if (!tau) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);

	sv = new double [nrhs*2];
	if (!sv) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);

// Init R matrix and R column

	int kii, igiv;

	for (kii=0;kii<ldr*ldr;kii++) rmatr[kii] = czero;
	for (kii=0;kii<ldr*nrhs;kii++) rcolmn[kii] = czero;
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structures
//
	CQrdC qrdu (treeclp, ncolmax+1, nrhs);
	CQrdC qrdv (treeclp, ncolmax+1, nrhs);

//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"SOFLancNE_U_",_tree.myid,"_");
		qrdu.SetupFiles (_param.ooctyp, strbuff, m);

		sprintf(strbuff, "%s%s%d%s",_param.path,"SOFLancNE_V_",_tree.myid,"_");
		qrdv.SetupFiles (_param.ooctyp, strbuff, n);

	};

//
// Compute the array of right hand side norms
//
	int irhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = _b.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs*m;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi0, resi);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi0[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi0[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Log10 || Rhs || = " << log10(resi0[irhs]) << endl;
		};
	};

// Init the subspace if necessary

	int jobitrloc = _param.jobitr;
	int jobitrloc_iabs = jobitrloc;
	if (jobitrloc_iabs < 0) jobitrloc_iabs = -jobitrloc_iabs;

	bool svconv = false;

	int ibegblk = 0;

	int ncorrloc = _corr.GetNcorr ();

	if (ncorrloc > 0 && jobitrloc_iabs > 0) {

		int nloc, mloc, ncolloc;
		dcmplx *pycorr, *pwcorr, *pdcorr, *ptaucorr;
		_corr.GetCorr (0, nloc, mloc, ncolloc,
							pycorr, pwcorr, pdcorr, ptaucorr);
		if (nloc != n || mloc != m) throw " CSVectorC::SOFLanczosNE: incompatible subspace sizes ";
		if (ncolloc > ncolmax) throw " CSVectorC::SOFLanczosNE: too many directions stored ";
		int i;
		for (i=0;i<n*ncolloc;i++) ppv[i] = pycorr[i];
		for (i=0;i<ncolloc;i++) tauv[i] = ptaucorr[i];
		ibegblk = ncolloc;
		if (nrhs != 1) throw " CSVectorC::SOFLanczosNE: nrhs > 1 is not implemented ";
		for (i=0;i<ibegblk;i++) {
			qrdv.UpdateQrdRPart (treeclp, n, i, i+1, ppv, n, tauv);
			IdentityBlockGivens (nrhs, givarr+i*nrhs_2*4);
		};
		if (jobitrloc > 0) {
			for (i=0;i<m*ncolloc;i++) ppu[i] = pwcorr[i];
			for (i=0;i<ncolloc;i++) tauu[i] = ptaucorr[ncolloc+i];
			for (i=0;i<ibegblk;i++) {
				qrdu.UpdateQrdRPart (treeclp, m, i, i+1, ppu, m, tauu);
			};
			int j;
			for (j=0;j<ncolloc;j++) {
				for (i=0;i<=j;i++) {
					rmatr[j*ldr+i] = ppu[j*m+i];
				};
			};
		} else {

			dcmplx *pv = v.GetVect ();
			dcmplx *pqwork = qwork.GetVect ();

			int j;
			for (j=0;j<ncolloc;j++) {

				for (kii=0;kii<ncolloc;kii++) beta[kii] = czero;
				beta[j] = cone;

				qrdv.MvmQBlk (treeclp, n, ncolloc, nrhs,
									ppv, n, tauv, 
									beta, ncolloc, 
									pv, n);

				time2 = clock ();

				(_solveu) (_objlu, v, pwork);

				time3 = clock ();

				timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				nslvu++;

				time2 = clock ();

				(_mvm) (_obja, pwork, qwork);

				time3 = clock ();

				timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				nmva++;

				ops += (nza + nzlu)*nrhs;

				int ibs = 0;
				if (_param.ooctyp == 0) ibs = j*nrhs*m;

				for (kii=0;kii<m*nrhs;kii++) ppu[ibs+kii] = czero;

				for (irhs=0;irhs<nrhs;irhs++) {
					for (kii=0;kii<m;kii++) {
						ppu[ibs+irhs*m+kii] = pqwork[irhs*m+kii];
					};
				};

				qrdu.UpdateQrd (treeclp, m, j*nrhs, j*nrhs+nrhs, ppu, m, tauu);

			};

			for (j=0;j<ncolloc;j++) {
				for (i=0;i<=j;i++) {
					rmatr[j*ldr+i] = ppu[j*m+i];
				};
			};

		};
		svconv = true;

	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;

	CSVectorC r0 = r;
//
// Compute QR decomposition of the residual block
//
//	_fout << " Residual = " << r << endl;
	dcmplx *pr = r.GetVect ();

	int ibs = 0;
	if (_param.ooctyp == 0) ibs = ibegblk*nrhs*m;

	for (kii=0;kii<m*nrhs;kii++) ppu[ibs+kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<m;kii++) {
			ppu[ibs+irhs*m+kii] = pr[irhs*m+kii];
		};
	};

	qrdu.UpdateQrd (treeclp, m, ibegblk*nrhs, ibegblk*nrhs+nrhs, ppu, m, tauu);

	ops += m*nrhs_2;
//
// Init the local right hand side
//
	int ibeg = ibegblk*nrhs;
	int iend = ibegblk*nrhs+nrhs;

	for (kii=0;kii<iend*nrhs;kii++) beta[kii] = czero;

	qrdu.GetRPartQrd (treeclp, 0, ibeg, iend,
							beta, iend);

// Store beta coefficient in a column

	int kjj;

	for (kii=0;kii<iend;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			rcolmn[kjj*ldr+kii] = beta[kjj*iend+kii];
		};
	};
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<iend;kii++) {
			double auxl = beta[irhs*iend+kii].real ();
			aux += auxl*auxl;
			auxl = beta[irhs*iend+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi1);
		for (irhs=0;irhs<nrhs;irhs++) resi[irhs] = resi1[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi[irhs];
		resi1[irhs] = sqrt(aux);
		resi[irhs] = resi1[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi1[irhs]) << endl;
		};
	};
//
// Compute initial directions block
//
	dcmplx *pu = u.GetVect ();

	for (kii=0;kii<iend*nrhs;kii++) beta[kii] = czero;
	for (kii=0;kii<nrhs;kii++) beta[kii*iend+ibeg+kii] = cone;

	qrdu.MvmQBlk (treeclp, m, iend, nrhs,
					ppu, m, tauu, 
					beta, iend, 
					pu, m);

	ops += 2*m*nrhs_2;
//	_fout << " Direction = " << u << endl;
//
// Multiply by the hermitian transposed preconditioned matrix
//
	time2 = clock ();

	(_mvmh) (_obja, u, v);

//	_fout << " V after Ah = " << v << endl;

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	(_solvel) (_objlu, v, w);
//	pwork = v;
//	_fout << " W after solve = " << w << endl;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;
//	_fout << " Pseudo Residual = " << w << endl;
//
// Compute QR decomposition
//
	dcmplx *pw = w.GetVect ();

	ibs = 0;
	if (_param.ooctyp == 0) ibs = ibegblk*nrhs*n;

	for (kii=0;kii<n*nrhs;kii++) ppv[ibs+kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			ppv[ibs+irhs*n+kii] = pw[irhs*n+kii];
		};
	};

	qrdv.UpdateQrd (treeclp, n, ibegblk*nrhs, ibegblk*nrhs+nrhs, ppv, n, tauv);

	ops += n*(ibeg*nrhs)*nrhs;
//
// Compute initial directions block
//
	dcmplx *pv = v.GetVect ();

	for (kii=0;kii<iend*nrhs;kii++) beta[kii] = czero;
	for (kii=0;kii<nrhs;kii++) beta[kii*iend+ibeg+kii] = cone;

	qrdv.MvmQBlk (treeclp, n, iend, nrhs,
						ppv, n, tauv, 
						beta, iend, 
						pv, n);

	ops += 2*n*nrhs_2;
//
// Init data for the main cycle
//
	w = v;
//
// Main iterative cycle
//
//	int ichk = _param.ichk;

	int nitrflt = 0;

	int iterglob = 0, iterlocal;
	int niterlocalmax;
	int iconv;

	int iendblk = ibegblk-1;
	int inextblk = iendblk+1;

	double svmin_r, svmax_r;

	while (iterglob<nitermax) {

// Implement local iterative cycle

		niterlocalmax = niterlocal;
		if (ncolmax-iendblk < niterlocalmax) niterlocalmax = ncolmax-iendblk;

		iterlocal = 0;
		while (iterlocal<niterlocalmax) {

// Check matrix relations
/*
			if (false) {

				if (nrhs>1) throw " CSVectorC::SOFLanczosNE: the case nrhs > 1 is not implemented yet ";

				dcmplx *yloc = new dcmplx [n*inextblk];
				dcmplx *ayloc = new dcmplx [m*(inextblk+1)];
				dcmplx *wloc = new dcmplx [m*(inextblk+1)];
				dcmplx *rrloc = new dcmplx [(inextblk+1)*(inextblk+1)];

				int i;

				for (i=0;i<inextblk*inextblk;i++) rrloc[i] = czero;
				for (i=0;i<inextblk;i++) rrloc[i*inextblk+i] = cone;

				if (inextblk != 0) {

					for (i=0;i<inextblk;i++) {
						qrdv.MvmQBlk (treeclp, n, inextblk, 1,
											ppv, n, tauv, 
											rrloc+i*inextblk, inextblk, 
											yloc+i*n, n);
					};

				};

				dcmplx *pr0 = r0.GetVect ();
				dcmplx *ppwork = pwork.GetVect ();
				dcmplx *pqwork = qwork.GetVect ();

				int j;

				for (j=0;j<m;j++) {
					ayloc[j] = pr0[j];
				};
				for (i=0;i<inextblk;i++) {
					for (j=0;j<n;j++) {
						ppwork[j] = yloc[i*n+j];
					};
					(_solveu) (_objlu, pwork, pwork1);
					(_mvm) (_obja, pwork1, qwork);
					for (j=0;j<m;j++) {
						ayloc[(i+1)*m+j] = pqwork[j];
					};
				};

				for (i=0;i<(inextblk+1)*(inextblk+1);i++) rrloc[i] = czero;

				qrdu.GetRPartQrd (treeclp, 0, 0, inextblk+1,
										rrloc, inextblk+1);

				OutArr (_fout," rrloc ",(inextblk+1)*(inextblk+1),rrloc);

				for (i=0;i<inextblk+1;i++) {
					qrdu.MvmQBlk (treeclp, m, inextblk+1, 1,
										ppu, m, tauu, 
										rrloc+i*(inextblk+1), inextblk+1, 
										wloc+i*m, m);
				};

				OutArr (_fout," wloc ",(inextblk+1)*m,wloc);

				int ind;
				double diff;
				dcmplx caux;
				for (i=0;i<inextblk+1;i++) {
					diff = 0.0e0;
					for (j=0;j<m;j++) {
						ind = i*m+j;
						caux = wloc[ind]-ayloc[ind];
						diff += caux.x*caux.x + caux.y*caux.y;
					};
					diff = sqrt(diff);
					cout << " Icol in check = " << i << " Diff = " << diff << endl;
				};

				delete [] yloc;
				delete [] ayloc;
				delete [] wloc;
				delete [] rrloc;
			};
*/
// Next iteration

			iterlocal++;
			iterglob++;

			iendblk++;
			inextblk++;
//
// Multiply by the preconditioned matrix
//
			time2 = clock ();

//			_fout << " y dir = " << v << endl;
			(_solveu) (_objlu, v, pwork);

			time3 = clock ();

			timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvu++;

			time2 = clock ();

			(_mvm) (_obja, pwork, qwork);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += (nza + nzlu)*nrhs;
//			_fout << " Ay = " << qwork << endl;
//
// Compute QR decomposition of the result
//
			dcmplx *pqwork = qwork.GetVect ();

			ibs = 0;
			if (_param.ooctyp == 0) ibs = inextblk*nrhs*m;

			for (kii=0;kii<m*nrhs;kii++) ppu[ibs+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<m;kii++) {
					ppu[ibs+irhs*m+kii] = pqwork[irhs*m+kii];
				};
			};

			qrdu.UpdateQrd (treeclp, m, inextblk*nrhs, inextblk*nrhs+nrhs, ppu, m, tauu);

			ops += m*nrhs_2;
//
// Take current block column of H_k
//
			ibeg = inextblk*nrhs;
			iend = inextblk*nrhs+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) alpha[kii] = czero;

			qrdu.GetRPartQrd (treeclp, 0, ibeg, iend,
									alpha, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (igiv=0;igiv<iendblk;igiv++) {

				ApplyBlockGivens (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = alpha[irhs*iend+iendblk*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
								giv, tau, giv1, givarr+iendblk*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
									iendblk, givarr+iendblk*nrhs_2*4,
									nrhs, alpha, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<inextblk*nrhs;kii++) {
					rmatr[(iendblk*nrhs+irhs)*ldr+kii] = alpha[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
									iendblk, givarr+iendblk*nrhs_2*4,
									nrhs, rcolmn, ldr, giv);

			ops += 8*nrhs_3;
//
// Compute residual norm
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = rcolmn[irhs*ldr+inextblk*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = rcolmn[irhs*ldr+inextblk*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

			};
//
// Compute current U
//
			pu = u.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) alpha[kii] = czero;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					alpha[kjj*(inextblk+1)*nrhs+inextblk*nrhs+kii] = rcolmn[kjj*ldr+inextblk*nrhs+kii];
				};
			};

			for (igiv=iendblk;igiv>=0;igiv--) {

				ApplyBlockGivensH (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, (inextblk+1)*nrhs, giv);

				ops += 8*nrhs_3;

			};

			qrdu.MvmQBlk (treeclp, m, (inextblk+1)*nrhs, nrhs,
								ppu, m, tauu, 
								alpha, (inextblk+1)*nrhs, 
								pu, m);

			ops += 2*m*(inextblk+1)*nrhs*nrhs;
//
// Multiply by the hermitian transposed preconditioned matrix
//
			time2 = clock ();

			(_mvmh) (_obja, u, pwork1);

			time3 = clock ();

			timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nmva++;

			ops += nza*nrhs;

			time2 = clock ();

			(_solvel) (_objlu, pwork1, pwork);

			time3 = clock ();

			timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			nslvl++;

			ops += nzlu*nrhs;
//
// Compute QR decomposition of the result
//
			dcmplx *ppwork = pwork.GetVect ();

			ibs = 0;
			if (_param.ooctyp == 0) ibs = inextblk*nrhs*n;

			for (kii=0;kii<n*nrhs;kii++) ppv[ibs+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					ppv[ibs+irhs*n+kii] = ppwork[irhs*n+kii];
				};
			};

			qrdv.UpdateQrd (treeclp, n, inextblk*nrhs, inextblk*nrhs+nrhs, ppv, n, tauv);

			ops += n*(inextblk*nrhs)*nrhs;
//
// Get residual orhogonality
//
			ibeg = inextblk*nrhs;
			iend = inextblk*nrhs+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) alpha[kii] = czero;

			qrdv.GetRPartQrd (treeclp, 0, ibeg, iend,
									alpha, iend);
//
// Compute pseudo residual norm
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<iend;kii++) {
					double auxl = alpha[irhs*iend+kii].real ();
					aux += auxl*auxl;
					auxl = alpha[irhs*iend+kii].imag ();
					aux += auxl*auxl;
				};
				resort[irhs] = sqrt(aux);

			};
//
// Compute current V
//
			pv = v.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) beta[kii] = czero;
			for (kii=0;kii<nrhs;kii++) beta[kii*(inextblk+1)*nrhs+inextblk*nrhs+kii] = cone;

			qrdv.MvmQBlk (treeclp, n, (inextblk+1)*nrhs, nrhs,
							ppv, n, tauv, 
							beta, (inextblk+1)*nrhs, 
							pv, n);

			ops += 2*n*(inextblk*nrhs)*nrhs;
//
// Check convergence
//
			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
//				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
				if (resi[irhs]/resi0[irhs] < _param.eps || (svconv && resort[irhs] / resi0[irhs] < _param.eps)) iconv++;
				if (iterglob%_param.ichk == 0) {
					cout  << " Irh= " << irhs << " Kin= " << nrhs*ibegblk << " Klc= " << nrhs*inextblk << " Itr= " << iterglob << " Log10 ||R|| = " << log10(resi[irhs]) << " Log10 ||ROrt|| = " << log10(resort[irhs]) << endl;
					_fout << " Irh= " << irhs << " Kin= " << nrhs*ibegblk << " Klc= " << nrhs*inextblk << " Itr= " << iterglob << " Log10 ||R|| = " << log10(resi[irhs]) << " Log10 ||ROrt|| = " << log10(resort[irhs]) << endl;
				};
			};

			if (iconv == nrhs) goto exit;

			if (iterglob >= nitermax) goto exit;

		};

// Perform filtering of the matrix relations

//		_fout << " Current residual vector = " << u << endl;

		if (nitrflt == 0 && jobitrloc == -2) ibegblk = 0;

		SOFLanNEFilterRelations (_fout, 
											_tree, treeclp, 
											_obja, _mvm, _mvmh,
											_objlu, _solvel, _solveu,
											_param,
											ibegblk, iendblk, inextblk, 
											m, n, qrdu, qrdv, ppu, ppv, tauu, tauv,
											ldr, rmatr, rcolmn, 
											r0, u, w, pwork,
											givarr, giv, giv1,
											svmin_r, svmax_r);

		nitrflt++;

		svconv = false;
		if (svmin_r > _param.xmin*0.9e0 && svmax_r < _param.xmax*1.1e0) svconv = true;

//		cout << " svmin_r = " << svmin_r << " svmax_r = " << svmax_r << " svconv = " << svconv << endl;

//
// Compute QR decomposition of the result
//
		pv = v.GetVect ();

		ibs = 0;
		if (_param.ooctyp == 0) ibs = inextblk*nrhs*n;

		for (kii=0;kii<n*nrhs;kii++) ppv[ibs+kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				ppv[ibs+irhs*n+kii] = pv[irhs*n+kii];
			};
		};

		qrdv.UpdateQrd (treeclp, n, inextblk*nrhs, inextblk*nrhs+nrhs, ppv, n, tauv);

		ops += n*(inextblk*nrhs)*nrhs;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) beta[kii] = czero;
		for (kii=0;kii<nrhs;kii++) beta[kii*(inextblk+1)*nrhs+inextblk*nrhs+kii] = cone;

		qrdv.MvmQBlk (treeclp, n, (inextblk+1)*nrhs, nrhs,
						ppv, n, tauv, 
						beta, (inextblk+1)*nrhs, 
						pv, n);

		ops += 2*n*(inextblk*nrhs)*nrhs;
//
// Compute current U
//
		if (false) {

			pu = u.GetVect ();

			for (kii=0;kii<(inextblk+1)*nrhs_2;kii++) alpha[kii] = czero;
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					alpha[kjj*(inextblk+1)*nrhs+inextblk*nrhs+kii] = rcolmn[kjj*ldr+inextblk*nrhs+kii];
				};
			};

			for (igiv=iendblk;igiv>=0;igiv--) {

				ApplyBlockGivensH (nrhs, 
										igiv, givarr+igiv*nrhs_2*4,
										nrhs, alpha, (inextblk+1)*nrhs, giv);

				ops += 8*nrhs_3;

			};

			qrdu.MvmQBlk (treeclp, m, (inextblk+1)*nrhs, nrhs,
								ppu, m, tauu, 
								alpha, (inextblk+1)*nrhs, 
								pu, m);

			_fout << " Recomputed current residual vector = " << u << endl;
//			goto exit;
		};

	};
//
// Update solution vector
//
exit:;
//
// Solve triangular system
//
	int ncolp = (inextblk)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) alpha[kii] = czero;

	SolveTriangularBySVD (_fout, ncolp, nrhs,
									rmatr, ldr,
									rcolmn, ldr,
									alpha, ldeloc);

// Multiply V by solution vector

	pv = v.GetVect ();

	qrdv.MvmQBlk (treeclp, n, inextblk*nrhs, nrhs,
						ppv, n, tauv, 
						alpha, ldeloc, 
						pv, n);

	ops += 2*n*(inextblk)*nrhs*nrhs;

// Update solution

	time2 = clock ();

	(_solveu) (_objlu, v, pwork);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	sol.DaxpySVect (cone,pwork);
//
// Compute the final residual and its norm
//
	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = r.vect[irhs*m+kii].real ();
			aux += auxl*auxl;
			auxl = r.vect[irhs*m+kii].imag ();
			aux += auxl*auxl;
		};
		resi[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		};
	};

// Check orthogonality of the final residual

	time2 = clock ();

	(_mvmh) (_obja, r, v);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	time2 = clock ();

	(_solvel) (_objlu, v, w);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<n;kii++) {
			double auxl = w.vect[irhs*n+kii].real ();
			aux += auxl*auxl;
			auxl = w.vect[irhs*n+kii].imag ();
			aux += auxl*auxl;
		};
		resort[irhs] = sqrt(aux);
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || ResOrt || = " << log10(resort[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || ResOrt || = " << log10(resort[irhs]) << endl;
		};
	};

// Store the subspace if necessary

	ncorrloc = _corr.GetNcorr ();

	if (jobitrloc_iabs > 1) {
		int nsave = ibegblk*nrhs;
		dcmplx *dloc;
		dcmplx *tauvu;
		dloc = new dcmplx [nsave*nsave*4];
		if (!dloc) MemoryFail (funcname);
		tauvu = new dcmplx [nsave*2];
		if (!tauvu) MemoryFail (funcname);
		int i;
		for (i=0;i<nsave;i++) tauvu[i] = tauv[i];
		for (i=0;i<nsave;i++) tauvu[nsave+i] = tauu[i];
		_corr.StoreCorr (0,n,m,nsave,ppv,ppu,dloc,tauvu);
		delete [] dloc;
		delete [] tauvu;
//		_fout << " Computed corr = " << _corr << endl;
	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrdu.CloseFiles ();
		qrdv.CloseFiles ();

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute global nza and ops

	double arr0[2], arr1[2];

	arr0[0] = (double) nza;
	arr0[1] = ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													2, arr0, arr1);
	} else {
		arr1[0] = arr0[0];
		arr1[1] = arr0[1];
	};

//	double dnza = arr1[0];
	ops = arr1[1];

// Output Lanczos statistics

//	perf = 1.0e-6 * ops / tottim;
//	ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;


		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] ppu;
	delete [] ppv;
	delete [] tauu;
	delete [] tauv;
	delete [] rmatr;
	delete [] rcolmn;
	delete [] alpha;
	delete [] beta;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] giv2;
	delete [] tau;
	delete [] cwork;
	delete [] sv;
	delete [] rwork;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;
	delete [] resort0;
	delete [] resort;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
//========================================================================================
void CSVectorC::SOFLanNEFilterRelations (ofstream &_fout, // Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
														const CTree &_tree, const CTree &_treeclp, 
														void *_obja, CMVMC _mvm, CMVMHC _mvmh,
														void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
														const CSlvParam &_param,
														int &_ibegblk, int &_iendblk, int &_inextblk, 
														int _m, int _n, 
														CQrdC &_qrdu, CQrdC &_qrdv, dcmplx *_ppu, dcmplx *_ppv, 
														dcmplx *_tauu, dcmplx *_tauv,
														int _ldr, dcmplx *_rmatr, dcmplx *_rcolmn,
														CSVectorC &_r0, CSVectorC &_u, CSVectorC &_v, CSVectorC &_pwork,
														dcmplx *_givarr, dcmplx *_giv, dcmplx *_giv1,
														double &_svmin_r, double &_svmax_r) { 

	const char *funcname = "SOFLanNEFilterRelations";

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);

// Take filtering parameters

	double xmin = _param.xmin;
	double xmax = _param.xmax;

// Allocate the memory

	int nrhs_2 = nrhs*nrhs;
	int nloc = (_iendblk-_ibegblk+1)*nrhs;
	int ntot = (_inextblk)*nrhs;
	int lworkloc = 10*nloc;
	int info = 0;

	dcmplx *aloc;
	dcmplx *uloc;
	dcmplx *vloc;
	dcmplx *tauloc;
	double *svloc;
	dcmplx *workloc;
	double *rworkloc;
	dcmplx *solloc;
	dcmplx *uextloc;

	aloc = new dcmplx [nloc*nloc];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [nloc*nloc];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [nloc*nloc];
	if (!vloc) MemoryFail (funcname);
	tauloc = new dcmplx [nloc];
	if (!tauloc) MemoryFail (funcname);
	svloc = new double [nloc];
	if (!svloc) MemoryFail (funcname);
	workloc = new dcmplx [lworkloc];
	if (!workloc) MemoryFail (funcname);
	rworkloc = new double [lworkloc];
	if (!rworkloc) MemoryFail (funcname);
	solloc = new dcmplx [nrhs*nloc];
	if (!solloc) MemoryFail (funcname);
	uextloc = new dcmplx [nloc*ntot];
	if (!uextloc) MemoryFail (funcname);

// Compute new local vectors (singular and solution)

// Take a part of R matrix and right hand side

	int kii, kjj, irhs;
	int ishift = _ibegblk*nrhs;

	for (kii=0;kii<nloc*nloc;kii++) aloc[kii] = czero;

	for (kjj=0;kjj<nloc;kjj++) {
		for (kii=0;kii<=kjj;kii++) {
			aloc[kjj*nloc+kii] = _rmatr[(ishift+kjj)*_ldr+ishift+kii];
		};
	};

	for (kjj=0;kjj<nrhs;kjj++) {
		for (kii=0;kii<nloc;kii++) {
			solloc[kjj*nloc+kii] = _rcolmn[kjj*_ldr+ishift+kii];
		};
	};

// Modify rcolmn and rmatrix back

	int igiv;

	for (igiv=_iendblk;igiv>=_ibegblk;igiv--) {

		ApplyBlockGivensH (nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								nrhs, _rcolmn, _ldr, _giv);

	};

	for (kjj=0;kjj<nrhs;kjj++) {
		for (kii=(_ibegblk+1)*nrhs;kii<(_inextblk+1)*nrhs;kii++) {
			_rcolmn[kjj*_ldr+kii] = czero;
		};
	};

	for (kjj=ishift;kjj<ntot;kjj++) {
		for (kii=0;kii<_ldr;kii++) {
			_rmatr[kjj*_ldr+kii] = czero;
		};
	};

// Compute SVD

	zgesvd_ ("A", "A", &nloc, &nloc, 
				aloc, &nloc, svloc, uloc, &nloc, vloc, &nloc, 
				workloc, &lworkloc, rworkloc, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

	csvhyst (cout,  "SvRSubm", nloc, svloc);
	csvhyst (_fout, "SvRSubm", nloc, svloc);

	_svmax_r = svloc[0];
	_svmin_r = svloc[nloc-1];

//	OutArr (_fout," SvRSubm = ",nloc,svloc);

// Solve

	dcmplx caux;
	double daux;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kjj=0;kjj<nloc;kjj++) {
			caux = czero;
			for (kii=0;kii<nloc;kii++) {
				caux += conj(uloc[kjj*nloc+kii]) * solloc[irhs*nloc+kii];
			};
			aloc[irhs*nloc+kjj] = caux;
		};
	};

	for (kjj=0;kjj<nloc;kjj++) {
		daux = svloc[kjj];
		daux = 1.0e0 / daux;
		for (irhs=0;irhs<nrhs;irhs++) {
			caux = aloc[irhs*nloc+kjj] * daux;
			aloc[irhs*nloc+kjj] = caux;
		};
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kjj=0;kjj<nloc;kjj++) {
			caux = czero;
			for (kii=0;kii<nloc;kii++) {
				caux += conj(vloc[kjj*nloc+kii]) * aloc[irhs*nloc+kii];
			};
			solloc[irhs*nloc+kjj] = caux;
		};
	};

// Filter singular values

// Large:

	int nfiltr = 0;

	for (kii=0;kii<nloc;kii++) {
		if (svloc[kii] > xmax) {
			for (kjj=0;kjj<nloc;kjj++) {
				aloc[nfiltr*nloc+kjj] = conj(vloc[kjj*nloc+kii]);
			};
			nfiltr++;
		};
	};

// Small:

	for (kii=nloc-1;kii>=0;kii--) {
		if (svloc[kii] < xmin) {
			for (kjj=0;kjj<nloc;kjj++) {
				aloc[nfiltr*nloc+kjj] = conj(vloc[kjj*nloc+kii]);
			};
			nfiltr++;
		};
	};

// Additional

//	int nfiltr0 = nfiltr;

	if (nfiltr%nrhs != 0) throw " CSVectorC::SOFLanNEFilterRelations: the case nrhs != 1 not implemented ";

// Add solution data

	if (nfiltr > nloc-nrhs) nfiltr = nloc-nrhs;

	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			aloc[nfiltr*nloc+kjj] = solloc[kii*nloc+kjj];
		};
		nfiltr++;
	};

// Orthonormalize the result

	QrdBlk (nloc, nfiltr, aloc, nloc, tauloc);

	for (kii=0;kii<nfiltr*nfiltr;kii++) vloc[kii] = czero;

	for (kii=0;kii<nfiltr;kii++) vloc[kii*nfiltr+kii] = cone;

	MvmQBlk (nloc, nfiltr, nfiltr,
				aloc, nloc, tauloc, 
				vloc, nfiltr, uloc, nloc);

// Compute new direction vectors

	for (kii=0;kii<nfiltr*ntot;kii++) uextloc[kii] = czero;

	for (kii=0;kii<nfiltr;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			uextloc[kii*ntot+ishift+kjj] = uloc[kii*nloc+kjj];
		};
	};

	int nparts = nfiltr / nrhs;
	int iendblkloc = _ibegblk+nparts-1;
	int inextblkloc = iendblkloc+1;

	dcmplx *vdir;

	vdir = new dcmplx [_m*nfiltr];
	if (!vdir) MemoryFail (funcname);

	int iblk;

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		_qrdv.MvmQBlk (_treeclp, _n, _inextblk*nrhs, nrhs,
						_ppv, _n, _tauv, 
						uextloc+(iblk-_ibegblk)*nrhs*ntot, ntot, 
						vdir+(iblk-_ibegblk)*nrhs*_n, _n);

	};

// Update V part

	int ibs;

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		ibs = 0;
		if (_param.ooctyp == 0) ibs = iblk*nrhs*_n;

		for (kii=0;kii<_n*nrhs;kii++) _ppv[ibs+kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<_n;kii++) {
				_ppv[ibs+irhs*_n+kii] = vdir[((iblk-_ibegblk)*nrhs+irhs)*_n+kii];
			};
		};

		_qrdv.UpdateQrd (_treeclp, _n, iblk*nrhs, iblk*nrhs+nrhs, _ppv, _n, _tauv);

	};
//
// Compute directions explicitely
//
	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		for (kii=0;kii<nrhs*ntot;kii++) uextloc[kii] = czero;

		for (kii=0;kii<nrhs;kii++) {
			uextloc[kii*ntot+ishift+(iblk-_ibegblk)*nrhs+kii] = cone;
		};

		_qrdv.MvmQBlk (_treeclp, _n, inextblkloc*nrhs, nrhs,
						_ppv, _n, _tauv, 
						uextloc, ntot, 
						vdir+(iblk-_ibegblk)*nrhs*_n, _n);

	};

// Multiply by the preconditioned matrix and update QR decomposition

	dcmplx *pv = _v.GetVect ();

	for (iblk=_ibegblk;iblk<inextblkloc;iblk++) {

		int ibeg = iblk*nrhs;
		int iend = iblk*nrhs+nrhs;

		for (kii=0;kii<nrhs*_n;kii++) {
			pv[kii] = vdir[(iblk-_ibegblk)*nrhs*_n+kii];
		};
//
// Multiply by the preconditioned matrix
//
		(_solveu) (_objlu, _v, _pwork);

		(_mvm) (_obja, _pwork, _u);
//
// Compute QR decomposition of the result
//
		dcmplx *pu = _u.GetVect ();

		ibs = 0;
		if (_param.ooctyp == 0) ibs = ibeg*_m;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<_m;kii++) {
				_ppu[ibs+irhs*_m+kii] = pu[irhs*_m+kii];
			};
		};

		_qrdu.UpdateQrd (_treeclp, _m, ibeg, iend, 
								_ppu, _m, _tauu);
//
// Take current block column of H_k
//
		for (kii=0;kii<iend*nrhs;kii++) uextloc[kii] = czero;

		_qrdu.GetRPartQrd (_treeclp, 0, ibeg, iend,
									uextloc, iend);
//
// Compute new block Givens rotation
//
		IdentityBlockGivens (nrhs, 
									_givarr+iblk*nrhs_2*4);
//
// Store block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(iblk+1)*nrhs;kii++) {
				_rmatr[(iblk*nrhs+irhs)*_ldr+kii] = uextloc[irhs*iend+kii];
			};
		};

	};

// Add initial residual into the current relations

	int ibeg = inextblkloc*nrhs;
	int iend = inextblkloc*nrhs+nrhs;

	dcmplx *pr0 = _r0.GetVect ();

	ibs = 0;
	if (_param.ooctyp == 0) ibs = ibeg*_m;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<_m;kii++) {
			_ppu[ibs+irhs*_m+kii] = pr0[irhs*_m+kii];
		};
	};

	_qrdu.UpdateQrd (_treeclp, _m, ibeg, iend, 
							_ppu, _m, _tauu);
//
// Take current block column of H_k
//
	for (kii=0;kii<nrhs*_ldr;kii++) _rcolmn[kii] = czero;

	_qrdu.GetRPartQrd (_treeclp, 0, ibeg, iend,
								_rcolmn, _ldr);

// Reassign indices that control iterative scheme

	_ibegblk = iendblkloc;
	_iendblk = iendblkloc;
	_inextblk = inextblkloc;

// Free work memory

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] tauloc;
	delete [] svloc;
	delete [] workloc;
	delete [] rworkloc;
	delete [] solloc;
	delete [] uextloc;
	delete [] vdir;

};

