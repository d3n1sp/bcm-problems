//------------------------------------------------------------------------------------------------
// File: lanczos.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#ifdef __PMPITRACE__
#include "mpi.h"
#include "mpe.h"
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "globals.h"
#include "slvparam.h"
#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "mvm.h"
#include "ExchangeMPI.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned Lanczos iterations
//========================================================================================
void CSVector::Lanczos (ofstream &_fout, // Perform preconditioned Lanczos iterations
						const CSMatrix *_mtra, const CSMatrix *_mtrl, const CSMatrix *_mtru, 
						const CSVector &_x, const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm) { 

//	const char *funcname = "Lanczos";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNza();
	nzlu = _mtrl->GetNza();
//
// Init guess to the solution
//
	CSVector sol(n);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVector bloc(n), r(n);
//
// Compute initial residual and its norm
//
	_mtra->MvmA (sol,r);

	ops += nza;

	bloc = _b;

	r = bloc-r;

	ops += n;

	double resini;
	double errmin;

	errmin = r.MaxProdSVect (_norm);
	resini = r.ScProdSVect (r);

	resini = sqrt(resini);

	ops += n;

	if (_param.msglev > 1) {
		if (_param.sttype == 2) {
			cout  << " Errmin ini = " << errmin << endl;
		} else {
			cout  << " Log10 (resini) = " << log10(resini) << endl;
		};
	};
	if (_param.msglev >= 1) {
		if (_param.sttype == 2) {
			_fout  << " Errmin ini = " << errmin << endl;
		} else {
			_fout  << " Log10 (resini) = " << log10(resini) << endl;
		};
	};
//
// Allocate Additional CSVector data
//
	CSVector goodsol(n), sol1(n), Delta(n), Delta1(n), delta(n), delta1(n), ImGdelta(n);
	CSVector p(n);
//
// Init Lanczos data
//
	goodsol = sol;
	Delta = r;
//
// Multiply by the preconditioner 
//
	_mtrl->SolveL (r,p); 
	_mtru->SolveU (p,delta); 

	ops += 2*nzlu;

//	int numMin    = 0; 
	double gamma  = 0., gamma1, ro = 1.;
	double delDel = 0., delDel1;

	sol1  .SetSVect(0.0e0);
	delta1.SetSVect(0.0e0);
	Delta1.SetSVect(0.0e0);
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	double resi;

	int k;
	int i;

	for (k=1;k<_param.niter;k++) {
//
// Multiply by the coefficient matrix
//
		_mtra->MvmA (delta,r);

		ops += nza;
//
// Multiply by the preconditioner 
//
		_mtrl->SolveL (r,p); 
		_mtru->SolveU (p,ImGdelta); 

		ops += 2*nzlu;
//
// Init gamma and ro
//
		gamma1  = gamma;
		delDel1 = delDel;

		delDel  = delta.ScProdSVect (Delta);

		double dtemp;

		dtemp = ImGdelta.ScProdSVect (Delta);

		gamma = delDel / dtemp;

		if (k>1) ro  = 1.0e0/(1.0e0-(gamma/gamma1)*(delDel/delDel1)/ro);

		ops += 2*n;
//
// Update local solution vector
//
		double rogamma = ro * gamma;
		double unitro  = 1. - ro;

//		p = sol * ro + sol1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = sol.vect[i] * ro + sol1.vect[i] * unitro;
		};

//		sol1 = sol;
		for (i=0;i<n;i++) {
			sol1.vect[i] = sol.vect[i];
		};
//		sol = p + delta * rogamma;
		for (i=0;i<n;i++) {
			sol.vect[i] = p.vect[i] + delta.vect[i] * rogamma;
		};

		ops += 3*n;
//
// Compute current residual
//
		_mtra->MvmA (sol,r);

		ops += nza;

//		r = bloc-r;
		for (i=0;i<n;i++) {
			r.vect[i] = bloc.vect[i]-r.vect[i];
		};

		ops += n;
//
// Check convergence
//
		double errloc = r.MaxProdSVect(_norm);
		resi = r.ScProdSVect (r);
		resi = sqrt(resi);

		ops += n;

		if (errloc < errmin) {
			goodsol = sol;
			errmin = errloc;
		};
		if (_param.msglev > 1) {
			if (_param.sttype == 2) {
				if (k%ichk == 0) cout  << " Iter = " << k << " Errmin = " << errmin << endl;
			} else {
				if (k%ichk == 0) cout  << " Iter = " << k << " Log10 (resi) = " << log10(resi) << endl;
			};
		};
		if (_param.msglev >= 1) {
			if (_param.sttype == 2) {
				if (k%ichk == 0) _fout  << " Iter = " << k << " Errmin = " << errmin << endl;
			} else {
				if (k%ichk == 0) _fout  << " Iter = " << k << " Log10 (resi) = " << log10(resi) << endl;
			};
		};

#ifdef __FlowVisionOut
		if (k%ichk == 0) {
			_param.meth->Err = errmin;
			_param.meth->N = k;
			if (_param.numpar) {
				for (int j=0; j<_param.numpar; j++) {
					_param.methj[j]->Err = _param.meth->Err;
					_param.methj[j]->N   = _param.meth->N;
					ErrorOutput (_param.pregion, _param.methj[j]);
				};
			} else {
				ErrorOutput (_param.pregion, _param.meth);
			};
		};
//		if (k%ichk == 0) cout  << " Errmin Fv = " << errmin << endl;
#endif
		_param.nitperf = k;
		_param.resfin = errmin;

		if (_param.sttype == 2) {
			if (k > _param.niter || errmin <= _param.eps ||
				(fabs (delDel) < 1.e-100 && fabs (delDel1) < 1.e-100)) {
				goto exit;
			};
		} else if (_param.sttype == 0) {
			if (k > _param.niter || resi <= _param.eps ||
				(fabs (delDel) < 1.e-100 && fabs (delDel1) < 1.e-100)) {
				goto exit;
			};
		} else if (_param.sttype == 1) {
			if (k > _param.niter || resi <= _param.eps * resini ||
				(fabs (delDel) < 1.e-100 && fabs (delDel1) < 1.e-100)) {
				goto exit;
			};
		};
//
// Update delta
//
//		p = delta * ro + delta1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = delta.vect[i] * ro + delta1.vect[i] * unitro;
		};
//		delta1 = delta;
		for (i=0;i<n;i++) {
			delta1.vect[i] = delta.vect[i];
		};
//		delta = p - ImGdelta * rogamma;
		for (i=0;i<n;i++) {
			delta.vect[i] = p.vect[i] - ImGdelta.vect[i] * rogamma;
		};

		ops += 3*n;
//
// Multiply by the transposed preconditioner 
//
		_mtru->SolveL (Delta,p);
		_mtrl->SolveU (p,r);

		ops += 2*nzlu;
//
// Multiply by the tranposed coefficient matrix
//
		_mtra->MvmAt (r,ImGdelta);

		ops += nza;
//
// Update Delta
//
//		p = Delta * ro + Delta1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = Delta.vect[i] * ro + Delta1.vect[i] * unitro;
		};
//		Delta1 = Delta;
		for (i=0;i<n;i++) {
			Delta1.vect[i] = Delta.vect[i];
		};
//		Delta = p - ImGdelta * rogamma;
		for (i=0;i<n;i++) {
			Delta.vect[i] = p.vect[i] - ImGdelta.vect[i] * rogamma;
		};

		ops += 3*n;

	};
//
// Return point
//
exit: 

	if (_param.msglev > 1) {
		if (_param.sttype == 2) {
			if (k%ichk > 0) cout  << " Iter = " << k << " Errmin = " << errmin << endl;
		} else {
			if (k%ichk == 0) cout  << " Iter = " << k << " Log10 (resi) = " << log10(resi) << endl;
		};
	};
	if (_param.msglev >= 1) {
		if (_param.sttype == 2) {
			if (k%ichk > 0) _fout  << " Iter = " << k << " Errmin = " << errmin << endl;
		} else {
			if (k%ichk == 0) _fout  << " Iter = " << k << " Log10 (resi) = " << log10(resi) << endl;
		};
	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timeitr = tottim;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	if (_param.msglev > 1) {
		cout  << " Lanczos statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " Lanczos statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Return the solution

	sol = goodsol;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned Lanczos iterations
//========================================================================================
void CSVector::Lanczos (ofstream &_fout, // Perform preconditioned Lanczos iterations
						const CTree &_tree, 
						void *_obja, CMVM _mvm, CMVM _mvmt,
						void *_objl, CSOLVEL _solvel, CSOLVEL _solvelt,
						void *_obju, CSOLVEU _solveu, CSOLVEU _solveut,
						const CSVector &_x, const CSVector &_b, const CSVector &_norm,
						CSlvParam &_param) { 

//	const char *funcname = "Lanczos";

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int lanc_beg;
	static int lanc_end;

	if (i_set == 0) {

		lanc_beg = MPE_Log_get_event_number ();
		lanc_end = MPE_Log_get_event_number ();
		MPE_Describe_state (lanc_beg, lanc_end, "Lanczos", "purple");

		i_set = 1;

	};

	MPE_Log_event (lanc_beg, 0, NULL);
#endif

	int nza, nzlu;
	double ops;
	double perf, tottim;
	clock_t time0;
	clock_t time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the sizes of the matrix
//
//	int ntotal = _gmtral.GetN();
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
	nzlu = 0;

//
// Init guess to the solution
//
	int npartsloc = _x.GetNparts ();

	CSVector sol(n,npartsloc,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVector r, y, w;
	CSVector bloc;

	r = _b;
	y = _b;
	w = _b;
	bloc = _b;
//
// Compute initial residual and its norm
//
	time2 = clock ();

//	_fout << " X bef MvmA ini = " << sol << endl;
	(_mvm) (_obja, sol, r);
//	_fout << " AX aft MvmA ini = " << r << endl;

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;

	double dnormb = 0.0e0;

	double errmin0, errmin, errminabs;

	if (_param.sttype == 2) {
		errmin0 = r.MaxProdSVect (_norm);
	} else {
		dnormb = bloc.ScProdSVect (bloc);
		errmin0 = r.ScProdSVect (r);
	};

	ops += n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		if (_param.sttype == 2) {
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, MAXIMUM,
														1, &errmin0, &errmin);
		} else {
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														1, &errmin0, &errmin);
		};
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													1, &dnormb, &dnormb);
	} else {
		errmin = errmin0;
	};

	errmin0 = errmin;

	if (_param.sttype != 2) {
		errmin0 = sqrt(errmin0);
		errmin = sqrt(errmin);
	};

	dnormb = sqrt(dnormb);

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			if (_param.sttype == 2) {
				cout  << " Errmin ini = " << errmin << endl;
			} else {
				cout  << " Log10 || b || = " << log10(dnormb) << " Log10 || R_0 || = " << log10(errmin0) << endl;
			};
		};
		if (_param.msglev >= 1) {
			if (_param.sttype == 2) {
				_fout  << " Errmin ini = " << errmin << endl;
			} else {
				_fout  << " Log10 || b || = " << log10(dnormb) << " Log10 || R_0 || = " << log10(errmin0) << endl;
			};
		};
	};

	_param.resfin = errmin;
//
// Allocate Additional CSVector data
//
	CSVector goodsol, sol1, Delta;
	CSVector Delta1, delta;
	CSVector delta1, ImGdelta;
	CSVector p;

	goodsol = sol;
	sol1 = sol;
	Delta = sol;
	Delta1 = sol;
	delta = sol;
	delta1 = sol;
	ImGdelta = sol;
	p = sol;
//
// Init Lanczos data
//
	goodsol = sol;
	Delta = r;
//
// Multiply by the preconditioner 
//
	time2 = clock ();

//	_fout << " X bef SlvL ini = " << r << endl;
	(_solvel) (_objl, r, p);
//	_fout << " X aft SlvL p = " << p << endl;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	ops += nzlu*nrhs;

	time2 = clock ();

//	_fout << " X bef SlvU ini = " << p << endl;
	(_solveu) (_obju, p, delta);
//	_fout << " PX aft SlvU ini = " << delta << endl;

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += nzlu*nrhs;

//	int numMin    = 0; 
	double gamma  = 0., gamma1, ro = 1.;
	double delDel = 0., delDel1;

	sol1  .SetSVect(0.0e0);
	delta1.SetSVect(0.0e0);
	Delta1.SetSVect(0.0e0);
//
// Main iterative cycle
//
//	int ichk = _param.ichk;
	int ichk = 1;

	int k, i;
	double aux;

	for (k=1;k<_param.niter;k++) {
//
// Multiply by the coefficient matrix
//
		time2 = clock ();

//		_fout << " X bef MvmA = " << delta << endl;
		(_mvm) (_obja, delta, r);
//		_fout << " AX aft MvmA = " << r << endl;

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;
//
// Multiply by the preconditioner 
//
		time2 = clock ();

//		_fout << " X bef SolveL = " << r << endl;
		(_solvel) (_objl, r, p);
//		_fout << " PX aft SolveL = " << p << endl;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		ops += nzlu*nrhs;

		time2 = clock ();

//		_fout << " X bef SolveU = " << p << endl;
		(_solveu) (_obju, p, ImGdelta);
//		_fout << " PX aft SolveU = " << ImGdelta << endl;

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		ops += nzlu*nrhs;
//
// Init gamma and ro
//
		gamma1  = gamma;
		delDel1 = delDel;

		delDel  = delta.ScProdSVect (Delta);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														1, &delDel, &aux);
			delDel = aux;
		};

		double dtemp;

		dtemp = ImGdelta.ScProdSVect (Delta);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														1, &dtemp, &aux);
			dtemp = aux;
		};

		gamma = delDel / dtemp;

		if (k>1) ro  = 1.0e0/(1.0e0-(gamma/gamma1)*(delDel/delDel1)/ro);

		ops += 2*n;
//
// Update local solution vector
//
		double rogamma = ro * gamma;
		double unitro  = 1. - ro;

//		p = sol * ro + sol1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = sol.vect[i] * ro + sol1.vect[i] * unitro;
		};

//		sol1 = sol;
		for (i=0;i<n;i++) {
			sol1.vect[i] = sol.vect[i];
		};
//		sol = p + delta * rogamma;
		for (i=0;i<n;i++) {
			sol.vect[i] = p.vect[i] + delta.vect[i] * rogamma;
		};

		ops += 3*n;
//
// Compute current residual
//
		time2 = clock ();

		(_mvm) (_obja, sol, r);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

//		r = bloc-r;
		for (i=0;i<n;i++) {
			r.vect[i] = bloc.vect[i]-r.vect[i];
		};

		ops += n;
//
// Check convergence
//
		double errloc;
		if (_param.sttype == 2) {
			errloc = r.MaxProdSVect (_norm);
		} else {
			errloc = r.ScProdSVect (r);
		};

		ops += n;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			if (_param.sttype == 2) {
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, MAXIMUM,
															1, &errloc, &aux);
			} else {
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, ADD,
															1, &errloc, &aux);
			};
			errloc = aux;
		};

		if (_param.sttype != 2) errloc = sqrt(errloc);

		if (errloc < errmin) {
			for (i=0;i<n;i++) {
				goodsol.vect[i] = sol.vect[i];
			};
//			goodsol = sol;
			errmin = errloc;
		};
		if (_tree.myid == 0 && _param.msglev > 1) {
			if (k%ichk == 0) {
				if (_param.sttype == 2) {
					cout  << " Iter = " << k << " Errmin = " << errmin << endl;
				} else {
					cout  << " Iter = " << k << " Log10 || R || = " << log10(errmin) << endl;
				};
			};
		};
		if (_tree.myid == 0 && _param.msglev >= 1) {
			if (k%ichk == 0) {
				if (_param.sttype == 2) {
					_fout  << " Iter = " << k << " Errmin = " << errmin << endl;
				} else {
					_fout  << " Iter = " << k << " Log10 || R || = " << log10(errmin) << endl;
				};
			};
		};

#ifdef __FlowVisionOut
		if (k%ichk == 0) {
			_param.meth->Err = errmin;
			_param.meth->N = k;
			if (_param.numpar) {
				for (int j=0; j<_param.numpar; j++) {
					_param.methj[j]->Err = _param.meth->Err;
					_param.methj[j]->N   = _param.meth->N;
					ErrorOutput (_param.pregion, _param.methj[j]);
				};
			} else {
				ErrorOutput (_param.pregion, _param.meth);
			};
		};
//		if (k%ichk == 0) cout  << " Errmin Fv = " << errmin << endl;
#endif

		errminabs = errmin;

		_param.nitperf = k;
		_param.resfin = errmin;

		if (k > _param.niter || (_param.sttype != 1 && errmin <= _param.eps) || (_param.sttype == 1 && errmin/dnormb <= _param.eps) ||
			(fabs (delDel) < 1.e-100 && fabs (delDel1) < 1.e-100)) {
			goto exit;
		};
//
// Update delta
//
//		p = delta * ro + delta1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = delta.vect[i] * ro + delta1.vect[i] * unitro;
		};
//		delta1 = delta;
		for (i=0;i<n;i++) {
			delta1.vect[i] = delta.vect[i];
		};
//		delta = p - ImGdelta * rogamma;
		for (i=0;i<n;i++) {
			delta.vect[i] = p.vect[i] - ImGdelta.vect[i] * rogamma;
		};

		ops += 3*n;
//
// Multiply by the transposed preconditioner 
//
		time2 = clock ();

		(_solveut) (_obju, Delta, p);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		ops += nzlu*nrhs;

		time2 = clock ();

		(_solvelt) (_objl, p, r);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		ops += nzlu*nrhs;
//
// Multiply by the transposed coefficient matrix
//
		time2 = clock ();

		(_mvmt) (_obja, r, ImGdelta);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;
//
// Update Delta
//
//		p = Delta * ro + Delta1 * unitro;
		for (i=0;i<n;i++) {
			p.vect[i] = Delta.vect[i] * ro + Delta1.vect[i] * unitro;
		};
//		Delta1 = Delta;
		for (i=0;i<n;i++) {
			Delta1.vect[i] = Delta.vect[i];
		};
//		Delta = p - ImGdelta * rogamma;
		for (i=0;i<n;i++) {
			Delta.vect[i] = p.vect[i] - ImGdelta.vect[i] * rogamma;
		};

		ops += 3*n;

	};
//
// Return point
//
exit: 

	if (_tree.myid == 0 && _param.msglev > 1) {
		if (k%ichk > 0) cout  << " Iter = " << k << " Errmin = " << errminabs << endl;
	};
	if (_tree.myid == 0 && _param.msglev >= 1) {
		if (k%ichk > 0) _fout << " Iter = " << k << " Errmin = " << errminabs << endl;
	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timeitr = tottim;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

//	ops = ops / (double) nza;

	if (_tree.myid == 0 && _param.msglev > 1) {
		cout  << " Lanczos statistics: " << endl;
//		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_tree.myid == 0 && _param.msglev >= 1) {
		_fout << " Lanczos statistics: " << endl;
//		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Return the solution

	sol = goodsol;

	*this = sol;

#ifdef __PMPITRACE__
	MPE_Log_event (lanc_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CSVector: Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
//========================================================================================
void CSVector::BlockLanczosRight (ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, // Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
								const CSVector &_x, const CSVector &_b, 
								const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosRight";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrl->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n,1,nrhs), soli(n,1,nrhs);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m,1,nrhs), r(m,1,nrhs), qwork(m,1,nrhs), u(m,1,nrhs);
//
// Allocate all working arrays
//
	double *alpha, *beta;
	double *alphaarr, *betaarr;
	double *cc, *ss, *ccbar, *ssbar;
	double *theta, *rho, *rhobar, *phi, *phibar;
	double *giv, *giv1, *giv2;
	double *tau, *resi0, *resi, *resort0, *resort, *work;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

	alpha = new double [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new double [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new double [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new double [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new double [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new double [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new double [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new double [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new double [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new double [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new double [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new double [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new double [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new double [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new double [nrhs*2];
	if (!tau) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmA (sol,r);

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute QR decomposition of the residual block
//
	double *pr = r.GetVect ();

	QrdBlk (m, nrhs, pr, tau);

	ops += m*nrhs_2;
//
// Compute initial directions block
//
	int kii, kjj, kkk, irhs;

	double *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (m, nrhs, nrhs, pr, tau, 
				ccbar, pu);

	ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (m, nrhs, pr, beta);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n,1,nrhs), w(n,1,nrhs), pwork(n,1,nrhs), pwork1(n,1,nrhs);

	_mtra->CSMatrixRS::MvmAt (u,v);

	_mtru->CSMatrixRS::SolveL (v,pwork);
	_mtrl->CSMatrixRS::SolveU (pwork,w);

	ops += (nza+2*nzlu) * nrhs;
//
// Compute QR decomposition
//
	double *pw = w.GetVect ();

	QrdBlk (n, nrhs, pw, tau);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	double *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (n, nrhs, nrhs, pw, tau, 
				ccbar, pv);

	ops += 2*n*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (n, nrhs, pw, alpha);
//
// Transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[irhs*nrhs+kii];
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
		_mtrl->CSMatrixRS::SolveL (v,pwork1);
		_mtru->CSMatrixRS::SolveU (pwork1,pwork);

		_mtra->CSMatrixRS::MvmA (pwork,qwork);

		ops += (nza + 2*nzlu)*nrhs;
//
// Update U
//
		qwork.DaxpySVectBlk (-1.0e0,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *pqwork = qwork.GetVect ();

		QrdBlk (m, nrhs, pqwork, tau);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (m, nrhs, nrhs, pqwork, tau, 
					ccbar, pu);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (m, nrhs, pqwork, beta);

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//
// Multiply by the transposed preconditioned matrix
//
		_mtra->CSMatrixRS::MvmAt (u,pwork);

		_mtru->CSMatrixRS::SolveL (pwork,pwork1);
		_mtrl->CSMatrixRS::SolveU (pwork1,pwork);

		ops += (nza + 2*nzlu)*nrhs_2;
//
// Compute transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = beta[kjj*nrhs+kii];
			};
		};
//
// Update V
//
		pwork.DaxpySVectBlk (-1.0e0,ccbar,v);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *ppwork = pwork.GetVect ();

		QrdBlk (n, nrhs, ppwork, tau);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (n, nrhs, nrhs, ppwork, tau, 
					ccbar, pv);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (n, nrhs, ppwork, alpha);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = 0.0e0;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = 0.0e0;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = 1.0e0;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, tau, 
					giv1, giv2);

		ops += 16*nrhs_3;

		GetRPartQrd (nrhs*2, nrhs*2, giv, giv1);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = giv2[kii*nrhs*2+kjj];
				ss[kjj*nrhs+kii] = giv2[kii*nrhs*2+nrhs+kjj];
				ssbar[kjj*nrhs+kii] = -giv2[(kii+nrhs)*nrhs*2+kjj];
				ccbar[kjj*nrhs+kii] = giv2[(kii+nrhs)*nrhs*2+nrhs+kjj];
			};
		};
//
// Check block Givens rotation
//
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
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = -aux;
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

		dgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
					work, &lwork, &info);
//		DGESVDDUM (&nrhs, &nrhs, 
//					giv, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
//					work, &lwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
			throw " Error inside Lapack routine DGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		csvhyst (cout,  "SvDiaLancz", nrhs, tau);
		csvhyst (_fout, "SvDiaLancz", nrhs, tau);
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) tau[kii] = 1.0e0 / tau[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= tau[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = aux;
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
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
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
		soli.DaxpySVectBlk (1.0e0,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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

		pwork.DaxpySVectBlk (-1.0e0,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			_mtrl->CSMatrixRS::SolveL (soli,qwork);
			_mtru->CSMatrixRS::SolveU (qwork,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmA (pwork,r);

			ops += (nza+2*nzlu)*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;

			};

			ops += m*nrhs;

			_mtra->CSMatrixRS::MvmAt (r,pwork1);

			_mtru->CSMatrixRS::SolveL (pwork1,qwork);
			_mtrl->CSMatrixRS::SolveU (qwork,pwork);

			ops += (nza+2*nzlu)*nrhs;

			ppwork = pwork.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = ppwork[irhs*n+kii];
					aux += auxl*auxl;
				};
				resort[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;

			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps || resort[irhs]/resort0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	_mtrl->CSMatrixRS::SolveL (soli,pwork1); 
	_mtru->CSMatrixRS::SolveU (pwork1,pwork); 

	ops += 2*nzlu*nrhs;

	sol.DaxpySVect (1.0e0,pwork);
//
// Compute SVD of B_k
//
	if (true) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		double *bmatr;
		double *svloc;
		double *workloc;

		bmatr = new double [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = 0.0e0;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		dgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, &info);

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

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
	delete [] resi0;
	delete [] resi;
	delete [] resort0;
	delete [] resort;
	delete [] work;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform left preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
//========================================================================================
void CSVector::BlockLanczosLeft (ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, // Perform left preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
								const CSVector &_x, const CSVector &_b, 
								const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosLeft";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrl->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n,1,nrhs), soli(n,1,nrhs);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m,1,nrhs), r(m,1,nrhs), qwork(m,1,nrhs), u(m,1,nrhs);
//
// Allocate all working arrays
//
	double *alpha, *beta;
	double *alphaarr, *betaarr;
	double *cc, *ss, *ccbar, *ssbar;
	double *theta, *rho, *rhobar, *phi, *phibar;
	double *giv, *giv1, *giv2;
	double *tau, *resi0, *resi, *resort0, *resort, *work;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

	alpha = new double [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new double [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new double [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new double [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new double [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new double [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new double [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new double [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new double [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new double [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new double [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new double [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new double [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new double [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new double [nrhs*2];
	if (!tau) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmA (sol,r);

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute array of initial residual norms
//
	int kii, kjj, kkk, irhs;

	double *pr = r.GetVect ();

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = pr[irhs*m+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Multiply by the preconditioner
//
	_mtrl->CSMatrixRS::SolveL (r,qwork);
	_mtru->CSMatrixRS::SolveU (qwork,r);

	ops += 2*nzlu*nrhs;
//
// Compute QR decomposition
//
	pr = r.GetVect ();

	QrdBlk (m, nrhs, pr, tau);

	ops += m*nrhs_2;
//
// Compute initial directions block
//
	double *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (m, nrhs, nrhs, pr, tau, 
				ccbar, pu);

	ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (m, nrhs, pr, beta);
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n,1,nrhs), w(n,1,nrhs), pwork(n,1,nrhs), pwork1(n,1,nrhs);

	_mtru->CSMatrixRS::SolveL (u,pwork);
	_mtrl->CSMatrixRS::SolveU (pwork,v);

	_mtra->CSMatrixRS::MvmAt (v,w);

	ops += (nza+2*nzlu) * nrhs;
//
// Compute QR decomposition
//
	double *pw = w.GetVect ();

	QrdBlk (n, nrhs, pw, tau);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	double *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (n, nrhs, nrhs, pw, tau, 
				ccbar, pv);

	ops += 2*n*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (n, nrhs, pw, alpha);
//
// Transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
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
		_mtra->CSMatrixRS::MvmA (v,pwork);

		_mtrl->CSMatrixRS::SolveL (pwork,pwork1);
		_mtru->CSMatrixRS::SolveU (pwork1,qwork);

		ops += (nza + 2*nzlu)*nrhs;
//
// Update U
//
		qwork.DaxpySVectBlk (-1.0e0,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *pqwork = qwork.GetVect ();

		QrdBlk (m, nrhs, pqwork, tau);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (m, nrhs, nrhs, pqwork, tau, 
					ccbar, pu);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (m, nrhs, pqwork, beta);

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//
// Multiply by the transposed preconditioned matrix
//
		_mtru->CSMatrixRS::SolveL (u,pwork);
		_mtrl->CSMatrixRS::SolveU (pwork,pwork1);

		_mtra->CSMatrixRS::MvmAt (pwork1,pwork);

		ops += (nza + 2*nzlu)*nrhs_2;
//
// Compute transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = beta[kjj*nrhs+kii];
			};
		};
//
// Update V
//
		pwork.DaxpySVectBlk (-1.0e0,ccbar,v);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *ppwork = pwork.GetVect ();

		QrdBlk (n, nrhs, ppwork, tau);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (n, nrhs, nrhs, ppwork, tau, 
					ccbar, pv);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (n, nrhs, ppwork, alpha);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = 0.0e0;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = 0.0e0;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = 1.0e0;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, tau, 
					giv1, giv2);

		ops += 16*nrhs_3;

		GetRPartQrd (nrhs*2, nrhs*2, giv, giv1);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = giv2[kii*nrhs*2+kjj];
				ss[kjj*nrhs+kii] = giv2[kii*nrhs*2+nrhs+kjj];
				ssbar[kjj*nrhs+kii] = -giv2[(kii+nrhs)*nrhs*2+kjj];
				ccbar[kjj*nrhs+kii] = giv2[(kii+nrhs)*nrhs*2+nrhs+kjj];
			};
		};
//
// Check block Givens rotation
//
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
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = -aux;
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

		dgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
					work, &lwork, &info);
//		DGESVDDUM (&nrhs, &nrhs, 
//					giv, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
//					work, &lwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
			throw " Error inside Lapack routine DGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		csvhyst (cout,  "SvDiaLancz", nrhs, tau);
		csvhyst (_fout, "SvDiaLancz", nrhs, tau);
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) tau[kii] = 1.0e0 / tau[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= tau[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = aux;
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
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
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
		soli.DaxpySVectBlk (1.0e0,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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

		pwork.DaxpySVectBlk (-1.0e0,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual
//
			pwork = sol + soli;

			_mtra->CSMatrixRS::MvmA (pwork,r);

			ops += (nza)*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;

			};

			ops += m*nrhs;

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

	sol.DaxpySVect (1.0e0,soli);
//
// Compute SVD of B_k
//
	if (true) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		double *bmatr;
		double *svloc;
		double *workloc;

		bmatr = new double [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = 0.0e0;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		dgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, &info);

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

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
	delete [] resi0;
	delete [] resi;
	delete [] resort0;
	delete [] resort;
	delete [] work;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform central preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
//========================================================================================
void CSVector::BlockLanczosCentral (ofstream &_fout, // Perform central preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
									const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, 
									const CSVector &_x, const CSVector &_b, 
									const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosCentral";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrl->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n,1,nrhs), soli(n,1,nrhs);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m,1,nrhs), r(m,1,nrhs), qwork(m,1,nrhs), u(m,1,nrhs);
//
// Allocate all working arrays
//
	double *alpha, *beta;
	double *alphaarr, *betaarr;
	double *cc, *ss, *ccbar, *ssbar;
	double *theta, *rho, *rhobar, *phi, *phibar;
	double *giv, *giv1, *giv2;
	double *tau, *resi0, *resi, *resort0, *resort, *work;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

	alpha = new double [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new double [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new double [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new double [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new double [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new double [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new double [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new double [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new double [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new double [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new double [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new double [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new double [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new double [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new double [nrhs*2];
	if (!tau) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmA (sol,r);

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute array of initial residual norms
//
	int kii, kjj, kkk, irhs;

	double *pr = r.GetVect ();

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<m;kii++) {
			double auxl = pr[irhs*m+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Multiply by the preconditioner
//
	_mtrl->CSMatrixRS::SolveL (r,qwork);

	r = qwork;

	ops += nzlu*nrhs;
//
// Compute QR decomposition
//
	pr = r.GetVect ();

	QrdBlk (m, nrhs, pr, tau);

	ops += m*nrhs_2;
//
// Compute initial directions block
//
	double *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (m, nrhs, nrhs, pr, tau, 
				ccbar, pu);

	ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (m, nrhs, pr, beta);
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n,1,nrhs), w(n,1,nrhs), pwork(n,1,nrhs), pwork1(n,1,nrhs);

	_mtrl->CSMatrixRS::SolveU (u,pwork);

	_mtra->CSMatrixRS::MvmAt (pwork,v);

	_mtru->CSMatrixRS::SolveL (v,w);

	ops += (nza+2*nzlu) * nrhs;
//
// Compute QR decomposition
//
	double *pw = w.GetVect ();

	QrdBlk (n, nrhs, pw, tau);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	double *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (n, nrhs, nrhs, pw, tau, 
				ccbar, pv);

	ops += 2*n*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (n, nrhs, pw, alpha);
//
// Transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
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
		_mtru->CSMatrixRS::SolveU (v,pwork);

		_mtra->CSMatrixRS::MvmA (pwork,pwork1);

		_mtrl->CSMatrixRS::SolveL (pwork1,qwork);

		ops += (nza + 2*nzlu)*nrhs;
//
// Update U
//
		qwork.DaxpySVectBlk (-1.0e0,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *pqwork = qwork.GetVect ();

		QrdBlk (m, nrhs, pqwork, tau);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (m, nrhs, nrhs, pqwork, tau, 
					ccbar, pu);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (m, nrhs, pqwork, beta);

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//
// Multiply by the transposed preconditioned matrix
//
		_mtrl->CSMatrixRS::SolveU (u,pwork);

		_mtra->CSMatrixRS::MvmAt (pwork,pwork1);

		_mtru->CSMatrixRS::SolveL (pwork1,pwork);

		ops += (nza + 2*nzlu)*nrhs_2;
//
// Compute transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = beta[kjj*nrhs+kii];
			};
		};
//
// Update V
//
		pwork.DaxpySVectBlk (-1.0e0,ccbar,v);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *ppwork = pwork.GetVect ();

		QrdBlk (n, nrhs, ppwork, tau);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (n, nrhs, nrhs, ppwork, tau, 
					ccbar, pv);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (n, nrhs, ppwork, alpha);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = 0.0e0;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = 0.0e0;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = 1.0e0;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, tau, 
					giv1, giv2);

		ops += 16*nrhs_3;

		GetRPartQrd (nrhs*2, nrhs*2, giv, giv1);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = giv2[kii*nrhs*2+kjj];
				ss[kjj*nrhs+kii] = giv2[kii*nrhs*2+nrhs+kjj];
				ssbar[kjj*nrhs+kii] = -giv2[(kii+nrhs)*nrhs*2+kjj];
				ccbar[kjj*nrhs+kii] = giv2[(kii+nrhs)*nrhs*2+nrhs+kjj];
			};
		};
//
// Check block Givens rotation
//
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
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = -aux;
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

		dgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
					work, &lwork, &info);
//		DGESVDDUM (&nrhs, &nrhs, 
//					giv, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
//					work, &lwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
			throw " Error inside Lapack routine DGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		csvhyst (cout,  "SvDiaLancz", nrhs, tau);
		csvhyst (_fout, "SvDiaLancz", nrhs, tau);
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) tau[kii] = 1.0e0 / tau[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= tau[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = aux;
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
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
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
		soli.DaxpySVectBlk (1.0e0,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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

		pwork.DaxpySVectBlk (-1.0e0,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual
//
			_mtru->CSMatrixRS::SolveU (soli,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmA (pwork,r);

			ops += (nza+nzlu)*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;

			};

			ops += m*nrhs;

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

	_mtru->CSMatrixRS::SolveU (soli,pwork1);

	ops += (nzlu)*nrhs;

	sol.DaxpySVect (1.0e0,pwork1);
//
// Compute SVD of B_k
//
	if (true) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		double *bmatr;
		double *svloc;
		double *workloc;

		bmatr = new double [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = 0.0e0;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		dgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, &info);

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

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
	delete [] resi0;
	delete [] resi;
	delete [] resort0;
	delete [] resort;
	delete [] work;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned Lanczos iterations with rectangular matrix
//========================================================================================
void CSVector::LanczosLsq (ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrlu, // Perform preconditioned Lanczos iterations with rectangular matrix
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param) { 

//	const char *funcname = "LanczosLsq";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrlu->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n), soli(n);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m), r(m), qwork(m), u(m);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmACol (sol,r);

	ops += nza;

	bloc = _b;

	r = bloc-r;

	ops += m;

	double resi0;

	resi0 = r.ScProdSVect (r);
	resi0 = sqrt(resi0);

	double beta = resi0;

	ops += m;

	cout  << " Initial Log10 || Resi || = " << log10(resi0) << endl;
	_fout << " Initial Log10 || Resi || = " << log10(resi0) << endl;
//
// Compute normalized residual vector
//
	double aux = 1.0e0 / beta;

	u = r;

	u.ScaleSVect (aux);

	ops += m;
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n), w(n), pwork(n), pwork1(n);

	_mtra->CSMatrixRS::MvmAtCol (u,w);

	_mtrlu->CSMatrixRS::SolveL (w,v);

	ops += nza+nzlu;
//
// Compute normalized direction vector
//
	aux = v.ScProdSVect (v);
	aux = sqrt(aux);

	double alpha = aux;

	aux = 1.0e0 / aux;

	v.ScaleSVect (aux);

	ops += n;

	double resort0 = alpha;

	cout  << " Initial Log10 || Reso || = " << log10(resort0) << endl;
	_fout << " Initial Log10 || Reso || = " << log10(resort0) << endl;
//
// Init data for the main cycle
//
	w = v;

	double phibar = beta;
	double rhobar = alpha;

	double resi, resort, rho, c, s, theta, phi;

	resi = resi0;
	resort = resort0;
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	for (int k=0;k<_param.niter;k++) {
//
// Multiply by the preconditioned matrix
//
		_mtrlu->CSMatrixRS::SolveU (v,pwork);

		_mtra->CSMatrixRS::MvmACol (pwork,qwork);

		ops += nza + nzlu;
//
// Update U
//
		qwork.DaxpySVect (-alpha,u);

		ops += m;

		u = qwork;
//
// Normalize U
//
		aux = u.ScProdSVect (u);
		aux = sqrt(aux);

		beta = aux;

		aux = 1.0e0 / aux;

		u.ScaleSVect (aux);

		ops += m;
//
// Multiply by the transposed preconditioned matrix
//
		_mtra->CSMatrixRS::MvmAtCol (u,pwork1);

		_mtrlu->CSMatrixRS::SolveL (pwork1,pwork);

		ops += nza + nzlu;
//
// Update V
//
		pwork.DaxpySVect (-beta,v);

		ops += n;

		v = pwork;
//
// Normalize direction vector
//
		aux = v.ScProdSVect (v);
		aux = sqrt(aux);

		alpha = aux;

		aux = 1.0e0 / aux;

		v.ScaleSVect (aux);

		ops += n;
//
// Compute rotation
//
		rho = rhobar*rhobar + beta*beta;

		rho = sqrt(rho);

		c = rhobar / rho;
		s = beta / rho;
//
// Update scalars
//
		theta = s * alpha;
		rhobar = c * alpha;
		phi = c * phibar;
		phibar = -s * phibar;
//
// Update solution and residual vectors
//
		aux = phi / rho;

		soli.DaxpySVect (aux,w);

		aux = theta / rho;

		pwork = v;

		pwork.DaxpySVect (-aux,w);

		w = pwork;

		ops += 2*n;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual
//
			_mtrlu->CSMatrixRS::SolveU (soli,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmACol (pwork,r);

			ops += nza+nzlu;

			r = bloc-r;

			ops += m;

			resi = r.ScProdSVect (r);
			resi = sqrt(resi);

			ops += m;

			cout  << " Iter = " << k << " Log10 || Resi || = " << log10(resi) << endl;
			_fout << " Iter = " << k << " Log10 || Resi || = " << log10(resi) << endl;

			_mtra->CSMatrixRS::MvmAtCol (r,pwork1);

			_mtrlu->CSMatrixRS::SolveL (pwork1,pwork);

			ops += nza+nzlu;

			aux = pwork.ScProdSVect (pwork);
			aux = sqrt(aux);

			resort = aux;

			cout  << " Iter = " << k << " Log10 || Reso || = " << log10(resort) << endl;
			_fout << " Iter = " << k << " Log10 || Reso || = " << log10(resort) << endl;

		};

		if (resi/resi0 < _param.eps || resort/resort0 < _param.eps) goto exit;

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	_mtrlu->CSMatrixRS::SolveU (soli,pwork); 

	ops += nzlu;

	sol.DaxpySVect (1.0e0,pwork);

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned block Lanczos iterations with rectangular matrix
//========================================================================================
void CSVector::BlockLanczosLsq (ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrlu, // Perform preconditioned block Lanczos iterations with rectangular matrix
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosLsq";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrlu->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n,1,nrhs), soli(n,1,nrhs);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m,1,nrhs), r(m,1,nrhs), qwork(m,1,nrhs), u(m,1,nrhs);
//
// Allocate all working arrays
//
	double *alpha, *beta;
	double *alphaarr, *betaarr;
	double *cc, *ss, *ccbar, *ssbar;
	double *theta, *rho, *rhobar, *phi, *phibar;
	double *giv, *giv1, *giv2;
	double *tau, *resi0, *resi, *resort0, *resort, *work;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

	alpha = new double [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new double [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new double [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new double [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new double [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new double [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new double [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new double [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new double [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new double [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new double [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new double [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new double [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new double [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new double [nrhs*2];
	if (!tau) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmACol (sol,r);

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute QR decomposition of the residual block
//
	double *pr = r.GetVect ();

	QrdBlk (m, nrhs, pr, tau);

	ops += m*nrhs_2;
//
// Compute initial directions block
//
	int kii, kjj, kkk, irhs;

	double *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (m, nrhs, nrhs, pr, tau, 
				ccbar, pu);

	ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (m, nrhs, pr, beta);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n,1,nrhs), w(n,1,nrhs), pwork(n,1,nrhs), pwork1(n,1,nrhs);

	_mtra->CSMatrixRS::MvmAtCol (u,v);

	_mtrlu->CSMatrixRS::SolveL (v,w);

	ops += (nza+nzlu) * nrhs;
//
// Compute QR decomposition
//
	double *pw = w.GetVect ();

	QrdBlk (n, nrhs, pw, tau);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	double *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (n, nrhs, nrhs, pw, tau, 
				ccbar, pv);

	ops += 2*n*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (n, nrhs, pw, alpha);
//
// Transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[irhs*nrhs+kii];
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
		_mtrlu->CSMatrixRS::SolveU (v,pwork);

		_mtra->CSMatrixRS::MvmACol (pwork,qwork);

		ops += (nza + nzlu)*nrhs;
//
// Update U
//
		qwork.DaxpySVectBlk (-1.0e0,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *pqwork = qwork.GetVect ();

		QrdBlk (m, nrhs, pqwork, tau);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (m, nrhs, nrhs, pqwork, tau, 
					ccbar, pu);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (m, nrhs, pqwork, beta);

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//
// Multiply by the transposed preconditioned matrix
//
		_mtra->CSMatrixRS::MvmAtCol (u,pwork1);

		_mtrlu->CSMatrixRS::SolveL (pwork1,pwork);

		ops += (nza + nzlu)*nrhs_2;
//
// Compute transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = beta[kjj*nrhs+kii];
			};
		};
//
// Update V
//
		pwork.DaxpySVectBlk (-1.0e0,ccbar,v);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *ppwork = pwork.GetVect ();

		QrdBlk (n, nrhs, ppwork, tau);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (n, nrhs, nrhs, ppwork, tau, 
					ccbar, pv);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (n, nrhs, ppwork, alpha);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = 0.0e0;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = 0.0e0;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = 1.0e0;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, tau, 
					giv1, giv2);

		ops += 16*nrhs_3;

		GetRPartQrd (nrhs*2, nrhs*2, giv, giv1);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = giv2[kii*nrhs*2+kjj];
				ss[kjj*nrhs+kii] = giv2[kii*nrhs*2+nrhs+kjj];
				ssbar[kjj*nrhs+kii] = -giv2[(kii+nrhs)*nrhs*2+kjj];
				ccbar[kjj*nrhs+kii] = giv2[(kii+nrhs)*nrhs*2+nrhs+kjj];
			};
		};
//
// Check block Givens rotation
//
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
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = -aux;
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

		dgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
					work, &lwork, &info);
//		DGESVDDUM (&nrhs, &nrhs, 
//					giv, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
//					work, &lwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
			throw " Error inside Lapack routine DGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		csvhyst (cout,  "SvDiaLancz", nrhs, tau);
		csvhyst (_fout, "SvDiaLancz", nrhs, tau);
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) tau[kii] = 1.0e0 / tau[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= tau[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = aux;
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
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
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
		soli.DaxpySVectBlk (1.0e0,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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

		pwork.DaxpySVectBlk (-1.0e0,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			_mtrlu->CSMatrixRS::SolveU (soli,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmACol (pwork,r);

			ops += (nza+nzlu)*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;

			};

			ops += m*nrhs;

			_mtra->CSMatrixRS::MvmAtCol (r,pwork1);

			_mtrlu->CSMatrixRS::SolveL (pwork1,pwork);

			ops += (nza+nzlu)*nrhs;

			ppwork = pwork.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = ppwork[irhs*n+kii];
					aux += auxl*auxl;
				};
				resort[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;

			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps || resort[irhs]/resort0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	_mtrlu->CSMatrixRS::SolveU (soli,pwork); 

	ops += nzlu*nrhs;

	sol.DaxpySVect (1.0e0,pwork);
//
// Compute SVD of B_k
//
	if (true) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		double *bmatr;
		double *svloc;
		double *workloc;

		bmatr = new double [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = 0.0e0;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		dgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, &info);

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

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
	delete [] resi0;
	delete [] resi;
	delete [] resort0;
	delete [] resort;
	delete [] work;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned Lanczos iterations with rectangular matrix in out-of-core mode
//========================================================================================
void CSVector::LanczosLsq (ofstream &_fout, // Perform preconditioned Lanczos iterations with rectangular matrix in out-of-core mode
							int _nblks, int *_blks, 
							const CSMatrixRS *_mtra, FILE **_mtrafiles, 
							const CSMatrixRS *_mtrlu, FILE **_mtrlufiles, 
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param) { 

//	const char *funcname = "LanczosLsq";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrlu->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n), soli(n);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m), r(m), qwork(m), u(m);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,sol,r);

	ops += nza;

	bloc = _b;

	r = bloc-r;

	ops += m;

	double resi0;

	resi0 = r.ScProdSVect (r);
	resi0 = sqrt(resi0);

	double beta = resi0;

	ops += m;

	cout  << " Initial Log10 || Resi || = " << log10(resi0) << endl;
	_fout << " Initial Log10 || Resi || = " << log10(resi0) << endl;
//
// Compute normalized residual vector
//
	double aux = 1.0e0 / beta;

	u = r;

	u.ScaleSVect (aux);

	ops += m;
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n), w(n), pwork(n), pwork1(n);

	_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,u,w);

	_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,w,v);

	ops += nza+nzlu;
//
// Compute normalized direction vector
//
	aux = v.ScProdSVect (v);
	aux = sqrt(aux);

	double alpha = aux;

	aux = 1.0e0 / aux;

	v.ScaleSVect (aux);

	ops += n;

	double resort0 = alpha;

	cout  << " Initial Log10 || Reso || = " << log10(resort0) << endl;
	_fout << " Initial Log10 || Reso || = " << log10(resort0) << endl;
//
// Init data for the main cycle
//
	w = v;

	double phibar = beta;
	double rhobar = alpha;

	double resi, resort, rho, c, s, theta, phi;

	resi = resi0;
	resort = resort0;
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	for (int k=0;k<_param.niter;k++) {
//
// Multiply by the preconditioned matrix
//
		_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,v,pwork);

		_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,pwork,qwork);

		ops += nza + nzlu;
//
// Update U
//
		qwork.DaxpySVect (-alpha,u);

		ops += m;

		u = qwork;
//
// Normalize U
//
		aux = u.ScProdSVect (u);
		aux = sqrt(aux);

		beta = aux;

		aux = 1.0e0 / aux;

		u.ScaleSVect (aux);

		ops += m;
//
// Multiply by the transposed preconditioned matrix
//
		_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,u,pwork1);

		_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,pwork1,pwork);

		ops += nza + nzlu;
//
// Update V
//
		pwork.DaxpySVect (-beta,v);

		ops += n;

		v = pwork;
//
// Normalize direction vector
//
		aux = v.ScProdSVect (v);
		aux = sqrt(aux);

		alpha = aux;

		aux = 1.0e0 / aux;

		v.ScaleSVect (aux);

		ops += n;
//
// Compute rotation
//
		rho = rhobar*rhobar + beta*beta;

		rho = sqrt(rho);

		c = rhobar / rho;
		s = beta / rho;
//
// Update scalars
//
		theta = s * alpha;
		rhobar = c * alpha;
		phi = c * phibar;
		phibar = -s * phibar;
//
// Update solution and residual vectors
//
		aux = phi / rho;

		soli.DaxpySVect (aux,w);

		aux = theta / rho;

		pwork = v;

		pwork.DaxpySVect (-aux,w);

		w = pwork;

		ops += 2*n;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual
//
			_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,soli,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,pwork,r);

			ops += nza+nzlu;

			r = bloc-r;

			ops += m;

			resi = r.ScProdSVect (r);
			resi = sqrt(resi);

			ops += m;

			cout  << " Iter = " << k << " Log10 || Resi || = " << log10(resi) << endl;
			_fout << " Iter = " << k << " Log10 || Resi || = " << log10(resi) << endl;

			_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,r,pwork1);

			_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,pwork1,pwork);

			ops += nza+nzlu;

			aux = pwork.ScProdSVect (pwork);
			aux = sqrt(aux);

			resort = aux;

			cout  << " Iter = " << k << " Log10 || Reso || = " << log10(resort) << endl;
			_fout << " Iter = " << k << " Log10 || Reso || = " << log10(resort) << endl;

		};

		if (resi/resi0 < _param.eps || resort/resort0 < _param.eps) goto exit;

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,soli,pwork); 

	ops += nzlu;

	sol.DaxpySVect (1.0e0,pwork);

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform preconditioned block Lanczos iterations with rectangular matrix in the out-of-core mode
//========================================================================================
void CSVector::BlockLanczosLsq (ofstream &_fout, // Perform preconditioned block Lanczos iterations with rectangular matrix in the out-of-core mode
							int _nblks, int *_blks, 
							const CSMatrixRS *_mtra, FILE **_mtrafiles, 
							const CSMatrixRS *_mtrlu, FILE **_mtrlufiles, 
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockLanczosLsq";

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int m = _mtra->GetM();
	int n = _mtra->GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra->GetNzatot();
	nzlu = _mtrlu->GetNzatot();
//
// Init guess to the solution
//
	CSVector sol(n,1,nrhs), soli(n,1,nrhs);

	sol = _x;

	soli.SetSVect(0.0e0);
//
// Allocate temporary vector data
//
	CSVector bloc(m,1,nrhs), r(m,1,nrhs), qwork(m,1,nrhs), u(m,1,nrhs);
//
// Allocate all working arrays
//
	double *alpha, *beta;
	double *alphaarr, *betaarr;
	double *cc, *ss, *ccbar, *ssbar;
	double *theta, *rho, *rhobar, *phi, *phibar;
	double *giv, *giv1, *giv2;
	double *tau, *resi0, *resi, *resort0, *resort, *work;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int lwork = 10*nrhs;

	alpha = new double [nrhs_2];
	if (!alpha) MemoryFail (funcname);
	beta = new double [nrhs_2];
	if (!beta) MemoryFail (funcname);
	alphaarr = new double [nrhs_2*(_param.niter+3)];
	if (!alphaarr) MemoryFail (funcname);
	betaarr = new double [nrhs_2*(_param.niter+3)];
	if (!betaarr) MemoryFail (funcname);
	cc = new double [nrhs_2];
	if (!cc) MemoryFail (funcname);
	ss = new double [nrhs_2];
	if (!ss) MemoryFail (funcname);
	ccbar = new double [nrhs_2];
	if (!ccbar) MemoryFail (funcname);
	ssbar = new double [nrhs_2];
	if (!ssbar) MemoryFail (funcname);
	theta = new double [nrhs_2];
	if (!theta) MemoryFail (funcname);
	rho = new double [nrhs_2];
	if (!rho) MemoryFail (funcname);
	rhobar = new double [nrhs_2];
	if (!rhobar) MemoryFail (funcname);
	phi = new double [nrhs_2];
	if (!phi) MemoryFail (funcname);
	phibar = new double [nrhs_2];
	if (!phibar) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	giv2 = new double [nrhs_2*4];
	if (!giv2) MemoryFail (funcname);

	tau = new double [nrhs*2];
	if (!tau) MemoryFail (funcname);
	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
	resort0 = new double [nrhs];
	if (!resort0) MemoryFail (funcname);
	resort = new double [nrhs];
	if (!resort) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,sol,r);

	ops += nza*nrhs;

	bloc = _b;

	r = bloc-r;

	ops += m*nrhs;
//
// Compute QR decomposition of the residual block
//
	double *pr = r.GetVect ();

	QrdBlk (m, nrhs, pr, tau);

	ops += m*nrhs_2;
//
// Compute initial directions block
//
	int kii, kjj, kkk, irhs;

	double *pu = u.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (m, nrhs, nrhs, pr, tau, 
				ccbar, pu);

	ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (m, nrhs, pr, beta);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = beta[irhs*nrhs+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += nrhs_2;
//
// Multiply by the transposed preconditioned matrix
//
	CSVector v(n,1,nrhs), w(n,1,nrhs), pwork(n,1,nrhs), pwork1(n,1,nrhs);

	_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,u,v);

	_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,v,w);

	ops += (nza+nzlu) * nrhs;
//
// Compute QR decomposition
//
	double *pw = w.GetVect ();

	QrdBlk (n, nrhs, pw, tau);

	ops += n*nrhs_2;
//
// Compute initial directions block
//
	double *pv = v.GetVect ();

	for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
	for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

	MvmQBlk (n, nrhs, nrhs, pw, tau, 
				ccbar, pv);

	ops += 2*n*nrhs_2;
//
// Get R part of the QR decomposition
//
	GetRPartQrd (n, nrhs, pw, alpha);
//
// Transpose alpha
//
	for (kii=0;kii<nrhs;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
		};
	};
	for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];
//
// Compute orthogonality of the initial residual
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = alpha[irhs*nrhs+kii];
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
		_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,v,pwork);

		_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,pwork,qwork);

		ops += (nza + nzlu)*nrhs;
//
// Update U
//
		qwork.DaxpySVectBlk (-1.0e0,alpha,u);

		ops += m*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *pqwork = qwork.GetVect ();

		QrdBlk (m, nrhs, pqwork, tau);

		ops += m*nrhs_2;
//
// Compute current U
//
		pu = u.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (m, nrhs, nrhs, pqwork, tau, 
					ccbar, pu);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (m, nrhs, pqwork, beta);

		indxarr++;

		for (kii=0;kii<nrhs_2;kii++) betaarr[indxarr*nrhs_2+kii] = beta[kii];
//
// Multiply by the transposed preconditioned matrix
//
		_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,u,pwork1);

		_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,pwork1,pwork);

		ops += (nza + nzlu)*nrhs_2;
//
// Compute transposed beta
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = beta[kjj*nrhs+kii];
			};
		};
//
// Update V
//
		pwork.DaxpySVectBlk (-1.0e0,ccbar,v);

		ops += n*nrhs_2;
//
// Compute QR decomposition of the result
//
		double *ppwork = pwork.GetVect ();

		QrdBlk (n, nrhs, ppwork, tau);

		ops += n*nrhs_2;
//
// Compute current V
//
		pv = v.GetVect ();

		for (kii=0;kii<nrhs_2;kii++) ccbar[kii] = 0.0e0;
		for (kii=0;kii<nrhs;kii++) ccbar[kii*nrhs+kii] = 1.0e0;

		MvmQBlk (n, nrhs, nrhs, ppwork, tau, 
					ccbar, pv);

		ops += 2*m*nrhs_2;
//
// Get R part of the QR decomposition
//
		GetRPartQrd (n, nrhs, ppwork, alpha);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				ccbar[kii*nrhs+kjj] = alpha[kjj*nrhs+kii];
			};
		};
		for (kii=0;kii<nrhs_2;kii++) alpha[kii] = ccbar[kii];

		for (kii=0;kii<nrhs_2;kii++) alphaarr[indxarr*nrhs_2+kii] = alpha[kii];
//
// Compute block Givens rotation:
//   _                  _     _            _       _       _
//  |     C        S     |   |  \bar{\rho}  |     |  \rho   |
//  |                    | * |              |  =  |         |
//  | -\bar{S}  \bar{C}  |   |     \beta    |     |    0    |
//  |_                  _|   |_            _|     |_       _|
//
//
		for (kii=0;kii<4*nrhs_2;kii++) giv[kii] = 0.0e0;

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs*2+kii] = rhobar[kjj*nrhs+kii];
				giv[kjj*nrhs*2+nrhs+kii] = beta[kjj*nrhs+kii];
			};
		};

		QrdBlk (2*nrhs, 2*nrhs, giv, tau);

		ops += nrhs*nrhs_2*8;

		for (kii=0;kii<nrhs_2*4;kii++) giv1[kii] = 0.0e0;
		for (kii=0;kii<nrhs*2;kii++) giv1[kii*nrhs*2+kii] = 1.0e0;

		MvmQBlk (2*nrhs, 2*nrhs, 2*nrhs, giv, tau, 
					giv1, giv2);

		ops += 16*nrhs_3;

		GetRPartQrd (nrhs*2, nrhs*2, giv, giv1);

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				rho[kjj*nrhs+kii] = giv1[kjj*nrhs*2+kii];
				cc[kjj*nrhs+kii] = giv2[kii*nrhs*2+kjj];
				ss[kjj*nrhs+kii] = giv2[kii*nrhs*2+nrhs+kjj];
				ssbar[kjj*nrhs+kii] = -giv2[(kii+nrhs)*nrhs*2+kjj];
				ccbar[kjj*nrhs+kii] = giv2[(kii+nrhs)*nrhs*2+nrhs+kjj];
			};
		};
//
// Check block Givens rotation
//
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
//
// Update block coefficients according to the formulae:
//
// \theta_{i+1} = s_i * \alpha_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ccbar[kkk*nrhs+kii] * alpha[kjj*nrhs+kkk];
				};
				rhobar[kjj*nrhs+kii] = aux;
			};
		};
//
// \phi_i = c_i * \bar{\phi}_{i+1}
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += ssbar[kkk*nrhs+kii] * phibar[kjj*nrhs+kkk];
				};
				ccbar[kjj*nrhs+kii] = -aux;
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

		dgesvd_ ("A", "A", &nrhs, &nrhs, 
					rho, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
					work, &lwork, &info);
//		DGESVDDUM (&nrhs, &nrhs, 
//					giv, &nrhs, tau, giv, &nrhs, giv1, &nrhs, 
//					work, &lwork, &info);
		if (info != 0) {
			cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
			throw " Error inside Lapack routine DGESVD";
		};

		ops += 16*nrhs_3;
//
// Output histogram of singular values
//
		csvhyst (cout,  "SvDiaLancz", nrhs, tau);
		csvhyst (_fout, "SvDiaLancz", nrhs, tau);
//
// Inverse
//
		for (kii=0;kii<nrhs;kii++) tau[kii] = 1.0e0 / tau[kii];

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				giv[kjj*nrhs+kii] *= tau[kjj];
			};
		};

		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
				for (kkk=0;kkk<nrhs;kkk++) {
					aux += giv1[kii*nrhs+kkk] * giv[kkk*nrhs+kjj];
				};
				rho[kjj*nrhs+kii] = aux;
			};
		};
		ops += nrhs_3;
//
// Check inverse
//
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
//
// Prepare block coefficient for solution update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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
		soli.DaxpySVectBlk (1.0e0,ccbar,w);

		ops += n*nrhs_2;
//
// Prepare block coefficient for w update
//
		for (kii=0;kii<nrhs;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				double aux = 0.0e0;
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

		pwork.DaxpySVectBlk (-1.0e0,ccbar,w);

		ops += n*nrhs_2;

		w = pwork;
//
// Check convergence
//
		if (k%ichk == 0) {
//
// Compute current residual and its orthogonality
//
			_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,soli,pwork1);

			pwork = sol + pwork1;

			_mtra->CSMatrixRS::MvmACol (_nblks,_blks,_mtrafiles,pwork,r);

			ops += (nza+nzlu)*nrhs;

			r = bloc-r;

			ops += m*nrhs;

			pr = r.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<m;kii++) {
					double auxl = pr[irhs*m+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;

			};

			ops += m*nrhs;

			_mtra->CSMatrixRS::MvmAtCol (_nblks,_blks,_mtrafiles,r,pwork1);

			_mtrlu->CSMatrixRS::SolveL (_nblks,_blks,_mtrlufiles,pwork1,pwork);

			ops += (nza+nzlu)*nrhs;

			ppwork = pwork.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<n;kii++) {
					double auxl = ppwork[irhs*n+kii];
					aux += auxl*auxl;
				};
				resort[irhs] = sqrt(aux);

				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Reso || = " << log10(resort[irhs]) << endl;

			};

			int iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps || resort[irhs]/resort0[irhs] < _param.eps) iconv++;
			};

			if (iconv == nrhs) goto exit;

		};

		if (k > _param.niter) goto exit;

	};
//
// Update solution vector
//
exit:;

	_mtrlu->CSMatrixRS::SolveU (_nblks,_blks,_mtrlufiles,soli,pwork); 

	ops += nzlu*nrhs;

	sol.DaxpySVect (1.0e0,pwork);
//
// Compute SVD of B_k
//
	if (true) {

		int ntot = indxarr * nrhs;
		int mtot = ntot + nrhs;
		int lworkloc = 10*mtot;
		int info = 0;

		double *bmatr;
		double *svloc;
		double *workloc;

		bmatr = new double [mtot*ntot];
		if (!bmatr) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		for (kii=0;kii<mtot*ntot;kii++) bmatr[kii] = 0.0e0;

		for (int isup=0;isup<indxarr;isup++) {
			for (kii=0;kii<nrhs;kii++) {
				for (kjj=0;kjj<nrhs;kjj++) {
					bmatr[(isup*nrhs+kjj)*mtot+isup*nrhs+kii] = alphaarr[isup*nrhs_2+kjj*nrhs+kii];
					bmatr[(isup*nrhs+kjj)*mtot+(isup+1)*nrhs+kii] = betaarr[(isup+1)*nrhs_2+kjj*nrhs+kii];
				};
			};
		};

		dgesvd_ ("N", "N", &mtot, &ntot, 
					bmatr, &mtot, svloc, bmatr, &mtot, bmatr, &mtot, 
					workloc, &lworkloc, &info);

		csvhyst (cout,  "SvCoefLancz", ntot, svloc);
		csvhyst (_fout, "SvCoefLancz", ntot, svloc);

		delete [] bmatr;
		delete [] svloc;
		delete [] workloc;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Lanczos statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " Lanczos statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Lanczos statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

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
	delete [] resi0;
	delete [] resi;
	delete [] resort0;
	delete [] resort;
	delete [] work;

// Return the solution

	*this = sol;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
