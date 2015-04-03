//------------------------------------------------------------------------------------------------
// File: cg.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>

#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "mvm.h"
#include "slvparam.h"
#include "globals.h"


using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A. 
// Description: Perform preconditioned CG iterations
// CSVector::Cg()
//========================================================================================
void CSVector::Cg (ofstream &_fout, // Perform preconditioned CG iterations
					const CSMatrix *_pa, int _istore, CGSMatrix *_pglu, const CSMatrix *_plu, const CSVector &_x, const CSVector &_b, 
					CSlvParam &_param, const CSVector &_norm) { 

	const char *funcname = "Cg";

	int k, it=0, ichk, n;
	double cosmx, rihr0, enrsq = 0., gamold = 0.;
	double gamma0 = 0., enrsq0 = 0., cosphi;
	double alpha, beta, gamma;
	double *oarr;
	double *darr;
	double errmin, errloc;

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
// 
// Allocate coefficients
//
	oarr = new double [_param.niter];
	if (!oarr) MemoryFail (funcname);
	darr = new double [_param.niter];
	if (!darr) MemoryFail (funcname);

//
// Get the size of the matrix
//
	n = _pa->GetN();

// Init statistics data

	time0 = clock ();

	ops = 0.0e0;

	nza = _pa->GetNza();
	nza = 2*nza;

	if (_istore == 0) {
		nzlu = _plu->GetNza();
	} else {
		nzlu = 0;
		int nblksloc = _pglu->GetNblksc ();
		CGSMatrixR *pglur = (CGSMatrixR *)_pglu;
		CSMatrixR *pmtrarr = pglur->GetMtrarr ();
		for (int i=0;i<nblksloc;i++) {
			nzlu += pmtrarr[i].GetNza ();
		};
	};

//
// Init guess to the solution
//
	CSVector sol(n);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVector bloc(n), r(n), p(n), q(n), z(n), goodsol(n);

	goodsol = _x;
//
// Compute initial residual and its norm
//
	_pa->MvmASymm (sol,r);

	ops += nza;

	bloc = _b;

	r = bloc-r;

	ops += n;

	bloc = r;

	if (_param.ittype == 2) {
		errmin = r.MaxProdSVect (_norm);
		ops += n;
	} else {
		errmin = 0.0e0;
	};

	if (_param.msglev > 1) {
		cout  << " Errmin ini = " << errmin << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Errmin ini = " << errmin << endl;
	};
//
// Init CG data
//
	alpha = 1.0e0;
	beta = 0.0e0;
	cosmx = 0.0e0;

	p.SetSVect (0.0e0);

	ichk = _param.ichk;
//
// Main iterative cycle
//
	int i;

	for (k=1;k<_param.niter;k++) {
//
// Multiply by the preconditioner if any
//
		if (_istore == 0) {
			_plu->SolveL (r,q);
			_plu->SolveU (q,z);
		} else {
			_pglu->SolveL (r,q);
			_pglu->SolveU (q,z);
		};

		ops += 2*nzlu;
//
// Compute necessary scalar products
//
		rihr0 = z.ScProdSVect (bloc);
		enrsq = r.ScProdSVect (r);
		gamma = r.ScProdSVect (z);

		ops += 3*n;
//
// Check convergence and output some information
//
		if (k == 1) gamma0 = gamma;
		if (k == 1) enrsq0 = enrsq;
		if (sqrt(gamma) == 0.0e0) goto label;
		cosphi = fabs(rihr0)/(sqrt(gamma)*sqrt(gamma0));
		if (k > 1 && cosphi > cosmx) cosmx = cosphi;
label:	it = k-1;

		if (it%ichk == 0) {
			if (_param.ittype == 0 ) {
				if (_param.msglev > 1) {
					cout  << " It = " << it << " Log10 (Resi/Res0) = " << 0.5*log10(enrsq/enrsq0) << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " It = " << it << " Log10 (Resi/Res0) = " << 0.5*log10(enrsq/enrsq0) << endl;
				};
			} 
			else {
				if (_param.msglev > 1) {
					cout  << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq) << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq) << endl;
				};
			};
		};

		if (_param.ittype == 2) {
			errloc = r.MaxProdSVect(_norm);
			ops += n;
		} else {
			errloc = 0.0e0;
		};

		if (_param.ittype != 2) {
//			goodsol = sol;
			for (i=0;i<n;i++) {
				goodsol.vect[i] = sol.vect[i];
			};
		};
//		if (_param.ittype == 2 && errloc < errmin) {
		if (_param.ittype == 2) {
//			goodsol = sol;
			for (i=0;i<n;i++) {
				goodsol.vect[i] = sol.vect[i];
			};
			errmin = errloc;
		};

		_param.nitperf = k;
		_param.resfin = errmin;

		if (_param.msglev > 1) {
			if (it%ichk == 0) cout  << " Errmin = " << errmin << endl;
		};
		if (_param.msglev >= 1) {
			if (it%ichk == 0) _fout << " Errmin = " << errmin << endl;
		};

#ifdef __FlowVisionOut
///*
		if (it%ichk == 0) {
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
//*/
//		if (it%ichk == 0) cout  << " Errmin Fv = " << errmin << endl;
#endif

		if (it > _param.niter) goto exit;
		if (_param.ittype == 0 && enrsq <= _param.eps*_param.eps*enrsq0) goto exit;
		if (_param.ittype == 1 && enrsq <= _param.eps*_param.eps) goto exit;
		if (_param.ittype == 2 && errmin <= fabs(_param.eps)) goto exit;

		if (k > 1) beta = gamma/gamold;
//
// Update direction vector
//
		p.ScaleSVect (beta);
		p.DaxpySVect (1.0e0,z);
		ops += 2*n;

//
// Update coefficients
//
		oarr[k-1] = beta/alpha*alpha;
		darr[k-1] = beta/alpha;
//
// Multiply by the coefficient matrix
//
		_pa->MvmASymm (p,z);

		ops += nza;

//
// Compute the scalar product
//
		alpha = p.ScProdSVect (z);

		alpha = gamma/alpha;
		darr[k-1] = darr[k-1] + 1.0e0/alpha;
		if (alpha < 0.0e0) {
			if (_param.msglev > 1) {
				cout  << " Matrix is not positive definite " << alpha << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Matrix is not positive definite " << alpha << endl;
			};
			goto exit;
		};

		gamold = gamma;
//
// Update direction vectors
//
		r.DaxpySVect (-alpha,z);
		sol.DaxpySVect (alpha,p);

		ops += 2*n;

	};

exit: 

	if (it%ichk > 0) {
		if (_param.eps > 0.0e0) {
			if (_param.msglev > 1) {
				cout  << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq/enrsq0) << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq/enrsq0) << endl;
			};
		} 
		else {
			if (_param.msglev > 1) {
				cout  << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq) << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " It = " << it << " Log10 (Resi) = " << 0.5*log10(enrsq) << endl;
			};
		};
	};

	delete [] darr;
	delete [] oarr;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timeitr = tottim;

// Output CG statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	if (_param.msglev > 1) {
		cout  << " CG statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " CG statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Return the solution

	sol = goodsol;

	*this = sol;

};

// Author: Kharchenko S.A. 
// Description: Perform preconditioned CG iterations
// CSVector::Cg()
//========================================================================================
void CSVector::Cg (ofstream &_fout, // Perform preconditioned CG iterations
							const CTree &_tree,
							void *_obja, CMVM _mvm,
							void *_objlu, CSOLVELU _solvelu,
							const CSVector &_x, const CSVector &_b, 
							CSlvParam &_param) { 

//	const char *funcname = "Cg";

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

	if (nrhs != 1) {
		throw " CSVectorC::Cg: the case nrhs > 1 is not implemented yet ";
	};
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
	CSVector sol(n);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVector bloc, r, u, au, v, w;

	bloc = _b;
	r = _b;
	u = _b;
	au = _b;
	v = _b;
	w = _b;
//
// Compute initial residual
//
	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	w = _b;

	r = r-w;

	ops += n*nrhs;
//
// Compute residual norm
//
	double resi0, resi;

	double aux = 0.0e0;
	double auxl;
	int i;

	for (i=0;i<n;i++) {
		auxl = r.vect[i];
		aux += auxl*auxl;
	};

	resi0 = aux;

	ops += 2*n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													1, &resi0, &resi);
	} else {
		resi = resi0;
	};

	aux = resi;
	resi0 = sqrt(aux);
	resi = resi0;

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " Initial Log10 || Resi || = " << log10(resi0) << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Initial Log10 || Resi || = " << log10(resi0) << endl;
		};
	};

// Compute initial direction vector

	time2 = clock ();

	(_solvelu) (_objlu, r, u);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	time2 = clock ();

//	_gmtru.SolveU (_tree,_mvm,w,u);

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	double uAu;
	double ru;
	double vAu;
	double acoef, bcoef;

	double daux;
	double darr[2], darr1[2];

	int k, it;

	for (k=1;k<_param.niter;k++) {

		it = k-1;

// Multiply by the coefficient matrix

		time2 = clock ();

		(_mvm) (_obja, u, au);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

// Compute scalar products (Au_i,u_i) and (r_{i-1},u_i)

		uAu = au.ScProdSVect (u);
		ru = r.ScProdSVect (u);

		if (_tree.nproc != 1) {
			darr[0] = uAu;
			darr[1] = ru;
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														2, &darr, &darr1);
			uAu = darr1[0];
			ru = darr1[1];
		};

// Check that matrix is positive definite

		if (uAu < 0.0e0) {
			if (_tree.myid == _tree.rootcpu) {
				if (_param.msglev > 1) {
					cout  << " Matrix is not positive definite " << uAu << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Matrix is not positive definite " << uAu << endl;
				};
			};
			goto exit;
		};

// Update residual vector

		aux = -1.0e0;

		daux = ru * aux;

		acoef = daux / uAu;

		for (i=0;i<n;i++) {
			daux = r.vect[i] + acoef * au.vect[i];
			r.vect[i] = daux;
		};

// Compute norm of the residual

		aux = 0.0e0;
		for (i=0;i<n;i++) {
			auxl = r.vect[i];
			aux += auxl*auxl;
		};

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														1, &aux, &resi);
			aux = resi;
		};

		resi = sqrt(aux);

		if (_tree.myid == _tree.rootcpu) {
			if (it%ichk == 0) {
//				if (_param.ittype == 0 ) {
//					if (_param.msglev > 1) {
//						cout  << " It = " << it << " Resi = " << log10(resi/resi0) << endl;
//					};
//					if (_param.msglev >= 1) {
//						_fout << " It = " << it << " Resi = " << log10(resi/resi0) << endl;
//					};
//				} 
//				else {
					if (_param.msglev > 1) {
						cout  << " It = " << it << " Log10 || Resi || = " << log10(resi) << endl;
					};
					if (_param.msglev >= 1) {
						_fout << " It = " << it << " Log10 || Resi || = " << log10(resi) << endl;
					};
//				};
			};
		};

		if (_param.sttype == 0 && resi <= _param.eps) {
			_param.nitperf = it;
			goto exit;
		};
		if (_param.sttype == 1 && resi <= _param.eps*resi0) {
			_param.nitperf = it;
			goto exit;
		};

// Compute new direction vector

		time2 = clock ();

		(_solvelu) (_objlu, r, v);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		time2 = clock ();

//		_gmtru.SolveU (_tree,_mvm,w,v);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		ops += 2*nzlu*nrhs;

// Compute one more scalar product (v,Au_i)

		vAu = v.ScProdSVect (au);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														1, &vAu, &daux);
			vAu = daux;
		};

// Update direction and solution vectors

//		bcoef = vAu / uAu;
		aux = -1.0e0;
		daux = vAu * aux;
		bcoef = daux / uAu;

		for (i=0;i<n;i++) {
			daux = sol.vect[i] + acoef * u.vect[i];
			sol.vect[i] = daux;
			daux = v.vect[i] + bcoef * u.vect[i];
			u.vect[i] = daux;
		};

	};

exit: 

// Compute the final residual

	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	w = _b;

	r = r-w;

	ops += n*nrhs;
//
// Compute residual norm
//
	aux = 0.0e0;

	for (i=0;i<n;i++) {
		auxl = r.vect[i];
		aux += auxl*auxl;
	};

	ops += 2*n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													1, &aux, &resi);
	} else {
		resi = aux;
	};

	resi = sqrt(resi);

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " Final Log10 || Resi || = " << log10(resi) << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Final Log10 || Resi || = " << log10(resi) << endl;
		};
	};

	_param.resfin = resi;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output CG statistics

	perf = 1.0e-6 * ops / tottim;

//	ops = ops / (double) nza;

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " CG statistics: " << endl;
//			cout  << "     Costs = " << ops << " MvmA flops. " << endl;
//			cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
			cout  << "     Time  = " << tottim << " sec."<< endl;
			cout  << " Mv-Slv functions statistics: " << endl;
			cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
					" MnTimeMvA = " << timemva / (double) nmva << endl;
			cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
					" MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
			cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
					" MnTimeSlU = " << timeslvu / (double)nslvu<< endl;
		};

		if (_param.msglev >= 1) {
			_fout << " CG statistics: " << endl;
//			_fout << "     Costs = " << ops << " MvmA flops. " << endl;
//			_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
			_fout << "     Time  = " << tottim << " sec."<< endl;
			_fout << " Mv-Slv functions statistics: " << endl;
			_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
					" MnTimeMvA = " << timemva / (double) nmva << endl;
			_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
					" MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
			_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
					" MnTimeSlU = " << timeslvu / (double)nslvu<< endl;
		};
	};

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A. 
// Description: Perform preconditioned CG iterations
// CSVectorC::Cg()
//========================================================================================
void CSVectorC::Cg (ofstream &_fout, // Perform preconditioned CG iterations
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtrau, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

//	const char *funcname = "Cg";

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

	if (nrhs != 1) {
		throw " CSVectorC::Cg: the case nrhs > 1 is not implemented yet ";
	};
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

	for (ilist=0;ilist<_gmtrau.nlist;ilist++) {
		int iblk = _gmtrau.listb[ilist];
		nza += _gmtrau.mtrarr[iblk].GetNzatot();
	};

	nza *= 2;

	nzlu = 0;
	for (ilist=0;ilist<_gmtru.nlist;ilist++) {
		int iblk = _gmtru.listb[ilist];
		nzlu += _gmtru.mtrarr[iblk].GetNzatot();
	};

//
// Init guess to the solution
//
	CSVectorC sol(n);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC bloc(n), r(n), u(n), au(n), v(n), w(n);

//
// Compute initial residual
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	gmtrdummy.MvmASymm (_tree, _mvm, _gmtrau, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	w = _b;

	r = r-w;

	ops += n*nrhs;
//
// Compute residual norm
//
	double resi0, resi;

	double aux = 0.0e0;
	double auxl;
	int i;

	for (i=0;i<n;i++) {
		auxl = r.vect[i].x;
		aux += auxl*auxl;
		auxl = r.vect[i].y;
		aux += auxl*auxl;
	};

	resi0 = aux;

	ops += 2*n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													1, &resi0, &resi);
	} else {
		resi = resi0;
	};

	aux = resi;
	resi0 = sqrt(aux);
	resi = resi0;

	if (_tree.myid == _tree.rootcpu) {
		cout  << " Initial Log10 || Resi || = " << log10(resi0) << endl;
		_fout << " Initial Log10 || Resi || = " << log10(resi0) << endl;
	};

// Compute initial direction vector

	time2 = clock ();

	_gmtru.SolveLConj (_tree,_mvm,r,w);

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	time2 = clock ();

	_gmtru.SolveU (_tree,_mvm,w,u);

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Main iterative cycle
//
	int ichk = _param.ichk;

	dcmplx uAu;
	dcmplx ru;
	dcmplx vAu;
	dcmplx acoef, bcoef;

	dcmplx caux;
	dcmplx carr[2], carr1[2];

	int k, it;

	for (k=1;k<_param.niter;k++) {

		it = k-1;

// Multiply by the coefficient matrix

		time2 = clock ();

		gmtrdummy.MvmASymm (_tree,_mvm,_gmtrau,u,au);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += nza*nrhs;

// Compute scalar products (Au_i,u_i) and (r_{i-1},u_i)

		uAu = au.ScProdSVect (u);
		ru = r.ScProdSVect (u);

		if (_tree.nproc != 1) {
			carr[0] = uAu;
			carr[1] = ru;
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														4, &carr, &carr1);
			uAu = carr1[0];
			ru = carr1[1];
		};

// Check that matrix is positive definite

		if (uAu.x < 0.0e0) {
			if (_param.msglev > 1) {
				cout  << " Matrix is not positive definite " << uAu.x << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Matrix is not positive definite " << uAu.x << endl;
			};
			goto exit;
		};

// Update residual vector

		aux = -1.0e0;

		caux = ru * aux;

		acoef = caux / uAu;

		for (i=0;i<n;i++) {
			caux = r.vect[i] + acoef * au.vect[i];
			r.vect[i] = caux;
		};

// Compute norm of the residual

		aux = 0.0e0;
		for (i=0;i<n;i++) {
			auxl = r.vect[i].x;
			aux += auxl*auxl;
			auxl = r.vect[i].y;
			aux += auxl*auxl;
		};
		resi = sqrt(aux);

		if (_tree.myid == _tree.rootcpu) {
			if (it%ichk == 0) {
//				if (_param.ittype == 0 ) {
//					if (_param.msglev > 1) {
//						cout  << " It = " << it << " Resi = " << log10(resi/resi0) << endl;
//					};
//					if (_param.msglev >= 1) {
//						_fout << " It = " << it << " Resi = " << log10(resi/resi0) << endl;
//					};
//				} 
//				else {
					if (_param.msglev > 1) {
						cout  << " It = " << it << " Resi = " << log10(resi) << endl;
					};
					if (_param.msglev >= 1) {
						_fout << " It = " << it << " Resi = " << log10(resi) << endl;
					};
//				};
			};
		};

		if (_param.ittype == 0 && resi <= _param.eps*resi0) goto exit;
		if (_param.ittype == 1 && resi <= _param.eps) goto exit;

// Compute new direction vector

		time2 = clock ();

		_gmtru.SolveLConj (_tree,_mvm,r,w);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		time2 = clock ();

		_gmtru.SolveU (_tree,_mvm,w,v);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		ops += 2*nzlu*nrhs;

// Compute one more scalar product (v,Au_i)

		vAu = v.ScProdSVect (au);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														2, &vAu, &caux);
			vAu = caux;
		};

// Update direction and solution vectors

//		bcoef = vAu / uAu;
		aux = -1.0e0;
		caux = vAu * aux;
		bcoef = caux / uAu;

		for (i=0;i<n;i++) {
			caux = sol.vect[i] + acoef * u.vect[i];
			sol.vect[i] = caux;
			caux = v.vect[i] + bcoef * u.vect[i];
			u.vect[i] = caux;
		};

	};

exit: 

// Compute the final residual

	time2 = clock ();

	gmtrdummy.MvmASymm (_tree, _mvm, _gmtrau, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	w = _b;

	r = r-w;

	ops += n*nrhs;
//
// Compute residual norm
//
	aux = 0.0e0;

	for (i=0;i<n;i++) {
		auxl = r.vect[i].x;
		aux += auxl*auxl;
		auxl = r.vect[i].y;
		aux += auxl*auxl;
	};

	ops += 2*n;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													1, &aux, &resi);
	} else {
		resi = aux;
	};

	resi = sqrt(resi);

	if (_tree.myid == _tree.rootcpu) {
		cout  << " Final Log10 || Resi || = " << log10(resi) << endl;
		_fout << " Final Log10 || Resi || = " << log10(resi) << endl;
	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output CG statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	if (_param.msglev > 1) {
		cout  << " CG statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;
	};

	if (_param.msglev >= 1) {
		_fout << " CG statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;
	};

// Return the solution

	*this = sol;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
