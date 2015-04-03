//------------------------------------------------------------------------------------------------
// File: sor.cpp
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
#include <iomanip>
#include <ctime>

#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "slvparam.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// Compute the value of the polynomial
//========================================================================================
double PolynomialValue (int _type, double _a1, double _a2, double _u) {

	switch (_type) {

		case 1:  return _a1*_u;

		case 2:  return _a2*_u*_u;

		case 3:  return _u*(_a1+_a2*_u);

		default: return 0.;

	};

};

// Author: Kharchenko S.A.
// Compute the derivative of the polynomial
//========================================================================================
double PolynomialDerivative (int _type, double _a1, double _a2, double _u) {

	switch (_type) {

		case 1: {
			return fabs(_a1);
		};

		case 2: {
			double r = _a2*_u;
			r       += r;
			return fabs(r);
		};

		case 3: {
			double r = _a2*_u;
			r       += r;
			r       += _a1;
			return fabs(r);
		};

		default: return 0.;

	};

};

// Author: Kharchenko S.A.
// CSVector: Compute right hand side for nonlinear solver
//========================================================================================
void CSVector::ComputeRhs (int _typePol, const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, // Compute right hand side for nonlinear solver
							const bool *_cmprhs, const CSVector &_sol) { 

//	const char *funcname = "ComputeRhs";

	for (int i=0;i<nv;i++) {
		if (_cmprhs[i]) vect[i] = _poly0.vect[i] + 
							PolynomialValue (_typePol, _poly1.vect[i], _poly2.vect[i], _sol.vect[i]);
	};

};

// Author: Kharchenko S.A.
// CSVector: Perform iterations of the nonlinear SOR algorithm
//========================================================================================
void CSVector::NonlinearSor (const CSMatrixR &_mtra, const CSVector &_x, // Perform iterations of the nonlinear SOR algorithm
								const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, const bool *_wall_E,
								CSlvParam &_param, const CSVector &_norm) {

	const char *funcname = "NonlinearSor";

	int n, nza;

	double ops, perf, tottim;
	clock_t time0, time1;
//
// Get the size of the matrix
//
	n = _mtra.GetN();
//
// Determine polinomial type
//
	double polMax1 = 0., polMax2 = 0.;

	int i;

	for (i=0; i<n; i++) {
		double r1 = fabs(_poly1.vect[i]);
		double r2 = fabs(_poly2.vect[i]);
		if (r1 > polMax1) polMax1 = r1;
		if (r2 > polMax2) polMax2 = r2;
	};

	int typePol = 0;

	if (polMax1 > 1.0e-10) typePol = 1;
	if (polMax2 > 1.0e-10) {
		if (typePol) typePol = 3;
		else         typePol = 2;
	};
	
//
// Allocate and init boolean array
//
	bool *calcResidual, *calcDelta;

	calcResidual = new bool [n];
	if (!calcResidual) MemoryFail (funcname);
	calcDelta = new bool [n];
	if (!calcDelta) MemoryFail (funcname);

	for (i=0;i<n;i++) calcResidual[i] = !_wall_E[i];
	for (i=0;i<n;i++) calcDelta[i]    = !_wall_E[i];
//
// Compute indices of the diagonal entries of a
//
	int *idiag;

	idiag = new int [n];
	if (!idiag) MemoryFail (funcname);

	for (i=0;i<n;i++) idiag[i] = 0;

	for (i=0;i<n;i++) {
		for (int j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
			int jj = _mtra.ja[j];
			if (jj == i) idiag[i] = j;
		};
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = _mtra.GetNzja();
//
// Init guess to the solution
//
	CSVector sol(n);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVector bloc(n), r(n), p(n), q(n), z(n), goodsol(n);
//
// Compute initial residual and its norm
//
//	ofstream fout ("ChkSor.dat");
	_mtra.MvmA (sol,r,calcResidual);
//	fout << " R after mvma = " << r << endl;

//	ops += nza;

	bloc.ComputeRhs (typePol, _poly0, _poly1, _poly2, calcResidual, sol);
//	fout << " bloc after compute rhs = " << bloc << endl;

	for (i=0;i<n;i++) {
		if (calcResidual[i]) {
			r.vect[i] = bloc.vect[i]-r.vect[i];
		};
	};

//	ops += n;

	double errmin;

	errmin = r.MaxProdSVect (_norm, calcResidual);
//	fout << " Norm = " << _norm << endl;
//	fout << " Resi = " << r << endl;
//	fout.close ();

//	ops += n;

	cout << " Errmin ini = " << errmin << endl;
//
// Init SOR data
//
	double invomega = 1.0e0 / _param.omega;
//	double dsclr = invomega-1.0e0;

	int ichk = 5;
//
// Main iterative cycle
//
	CSVector diagl(n);

	int k;

	for (k=1;k<_param.niter;k++) {
//
// Update calcResidual
//
		for (i=0;i<n;i++) calcResidual[i] = calcResidual[i] || calcDelta[i];
//
// Prepare triangular system for solves
//
		for (i=0;i<n;i++) {
			if (calcDelta[i]) {
				double der = PolynomialDerivative (typePol, _poly1.vect[i], _poly2.vect[i], sol.vect[i]);
				double aux = _mtra.a[idiag[i]] * invomega + der;
				aux = 1.0e0 / aux;
				diagl.vect[i] = aux;
			};
		};
//		ops += 2*n;
//
// Solve triangular system with al
//
		_mtra.SolveLA (diagl,r,p,calcDelta);

//		ops += nzal;
//
// Update solution vector
//
		for (i=0;i<n;i++) {
			if (calcDelta[i]) {
				sol.vect[i] += p.vect[i];
				if (sol.vect[i] < _param.permismin) sol.vect[i] = _param.permismin;
			};
		};

//		ops += n;
//
// Compute new right hand side vector
//
		bloc.ComputeRhs (typePol, _poly0, _poly1, _poly2, calcResidual, sol);
//
// Update residual vector
//
		_mtra.MvmA (sol,r,calcResidual);

//		ops += nza;

		for (i=0;i<n;i++) {
			if (calcResidual[i]) {
				r.vect[i] = bloc.vect[i]-r.vect[i];
			};
		};

//		ops += n;
//
// Check convergence
//
		double errloc = 0.0e0;
		double aux;

		for (i=0;i<n;i++) {
			if (calcResidual[i]) {
				aux = fabs(r.vect[i]*_norm.vect[i]);
				if (aux < _param.eps) { 
					calcDelta[i] = false;
				} else {
					calcDelta[i] = true;
				};
				if (aux > errloc) errloc = aux;
			};
		};

//		ops += n;

		if (errloc < errmin) {
			goodsol = sol;
			errmin = errloc;
		};

		if (k%ichk == 0) cout << " Iter = " << k << " Errmin = " << errmin << endl;

		if (k > _param.niter) goto exit;

		if (errmin <= fabs(_param.eps)) {
//
// Perform final check of the solution
//
// Compute initial residual and its norm
//
			for (i=0;i<n;i++) calcResidual[i] = !_wall_E[i];

			_mtra.MvmA (sol,r,calcResidual);

//			ops += nza;

			bloc.ComputeRhs (typePol, _poly0, _poly1, _poly2, calcResidual, sol);

			for (i=0;i<n;i++) {
				if (calcResidual[i]) {
					r.vect[i] = bloc.vect[i]-r.vect[i];
				};
			};

//			ops += n;

			for (i=0;i<n;i++) {
				aux = fabs(r.vect[i]*_norm.vect[i]);
				if (aux < _param.eps || _wall_E[i]) { 
					calcDelta[i] = false;
				} else {
					calcDelta[i] = true;
				};
				if (calcResidual[i] && aux > errloc) errloc = aux;
			};

//			ops += n;

			if (errloc < _param.eps) goto exit;

		};
//
// Modify calcResidual
//
		for (i=0;i<n;i++) calcResidual[i] = false;

		for (i=0;i<n;i++) {
			if (calcDelta[i]) {
				for (int j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
					int kk = _mtra.ja[j];
					if (!calcDelta[kk]) calcResidual[kk] = true;
				};
			};
		};

		for (i=0;i<n;i++) {
			if (_wall_E[i]) calcResidual[i] = false;
		};

	};

exit: 

	if (k%ichk > 0) cout << " Iter = " << k << " Errmin = " << errmin << endl;
//
// Finalize time measurement
//
	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;
//
// Output SOR statistics
//
	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout << " SOR statistics: " << endl;
//	cout << "     Costs = " << ops << " MvmA flops. " << endl;
//	cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	cout << "     Time  = " << tottim << " sec. "<< endl;
//
// Delete working arrays
//
	delete [] calcResidual;
	delete [] calcDelta;
	delete [] idiag;
//
// Return the solution
//
	sol = goodsol;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel iterations of the nonlinear SOR algorithm
//========================================================================================
void CSVector::NonlinearSor (ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau,
										bool _second, CMvmR &_mvm2, CGSMatrixRS &_gmtra2,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param) {

	const char *funcname = "NonlinearSor";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();

// Get the local size of the matrix

	int nloc = _x.GetNv();

// Store initial diagonal data and compute pointers to the diagonal data

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;
	int *spndcl, *spndrl;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);
	_gmtral.GetSups (spndcl, spndrl);

	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	CSMatrixRS *mtralarr, *mtrauarr;

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

	double **dptrldiag, **dptrudiag;
	double *diag;

	dptrldiag = new double * [nloc];
	if (!dptrldiag) MemoryFail (funcname);
	dptrudiag = new double * [nloc];
	if (!dptrudiag) MemoryFail (funcname);
	diag = new double [nloc];
	if (!diag) MemoryFail (funcname);

	int *listl, *ial, *jal, *bsal;
	double *al, *au;

	nloc = 0;

	int blamxloc, iblk, i, ii, j, jj, nsupi, isup, ni, ibs, kii;

	for (ii=0;ii<nlistloc;ii++) {
		iblk = listbloc[ii];
		blamxloc = mtralarr[iblk].GetBlamx ();
		mtralarr[iblk].GetSparsity (listl, ial, jal);
		bsal = mtralarr[iblk].GetBsa ();
		al = mtralarr[iblk].GetA ();
		au = mtrauarr[iblk].GetA ();
		nsupi = bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
		for (i=0;i<nsupi;i++) {
			isup = listl[i];
			ni = spndcl[isup+1]-spndcl[isup];
			for (j=ial[i];j<ial[i+1];j++) {
				jj = jal[j];
				if (blamxloc == 1) {
					ibs = j;
				} else {
					ibs = bsal[j];
				};
				if (isup == jj) {
					for (kii=0;kii<ni;kii++) {
						dptrldiag[nloc] = al+ibs+kii*ni+kii;
						dptrudiag[nloc] = au+ibs+kii*ni+kii;
						nloc++;
					};
				};
			};
		};
	};

	double *dptr;

	for (i=0;i<nloc;i++) {
		dptr = dptrldiag[i];
		diag[i] = *dptr;
	};

// Init guess to the solution

	clock_t time0, time1;

	time0 = clock ();

	CSVector sol(nloc);

	sol = _x;

// Allocate temporary vector data

	CSVector bloc(nloc), r(nloc), p(nloc), q(nloc), z(nloc), goodsol(nloc);

	goodsol = _x;

// Compute initial residual and its norm

	CGSMatrixRS gmtrdummy;

	gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, sol, r);
	if (_second) {
		_gmtra2.MvmA2 (_tree, _mvm2, sol, r);
	};

	double aux;

	if (_itype == 1) {

		int itype = 0;

		for (i=0;i<nloc;i++) {
			aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
			r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
		};

	} else {

		for (i=0;i<nloc;i++) {
			r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
		};

	};

	double errmin;

	errmin = r.MaxProdSVect (_norm);

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, MAXIMUM,
													1, &errmin, &errmin);
	};

	if (_param.msglev > 1 && myid == 0) {
		cout << " Errmin ini = " << errmin << endl;
	};
	if (_param.msglev >= 1 && myid == 0) {
		_fout << " Errmin ini = " << errmin << endl;
	};

// Init SOR data

	double invomega = 1.0e0 / _param.omega;
//	double dsclr = invomega-1.0e0;

	int ichk = 5;

// Main iterative cycle

	int k;

	for (k=1;k<_param.niter;k++) {

// Prepare triangular system for solves

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			double der;
			if (_itype == 1) {

				int itype = 1;

				der = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);

			} else {

				der = 2.0e0 * _poly2.vect[i] * sol.vect[i];

			};
			aux = diag[i] * invomega - der;
			*dptr = 1.0e0 / aux;
		};

// Solve triangular system with al

		_gmtral.SolveL (_tree, _mvm,
								r, p);

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
		};

// Update solution vector

		for (i=0;i<nloc;i++) {
			sol.vect[i] += p.vect[i];
			if (sol.vect[i] < _param.permismin) sol.vect[i] = _param.permismin;
		};

// Compute new right hand side vector

		gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, sol, r);
		if (_second) {
			_gmtra2.MvmA2 (_tree, _mvm2, sol, r);
		};

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
			};

		};

		double errloc;

		errloc = r.MaxProdSVect (_norm);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, MAXIMUM,
														1, &errloc, &errloc);
		};

// Check convergence

		if (errloc < errmin) {
			goodsol = sol;
			errmin = errloc;
		};

		if (k%ichk == 0) {
			if (_param.msglev > 1 && myid == 0) {
				cout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
			if (_param.msglev >= 1 && myid == 0) {
				_fout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
		};

		if (errmin < _param.eps) goto exit;

		if (k > _param.niter) goto exit;

	};

exit: 

	if (k%ichk > 0) {
		if (_param.msglev > 1 && myid == 0) {
			cout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
	};
//
// Finalize time measurement
//
	time1 = clock ();

	double tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output SOR statistics

	if (_param.msglev > 1 && myid == 0) {
		cout << " SOR statistics: " << endl;
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;

// Return the solution and params

	_param.resfin = errmin;
	_param.nitperf = k;

	sol = goodsol;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel iterations of the nonlinear SOR algorithm
//========================================================================================
void CSVector::NonlinearSor2Index (ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
										bool _second, CMvmR &_mvm2, CGSMatrixR &_gmtra2,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param) {

	const char *funcname = "NonlinearSor2Index";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();

// Get the local size of the matrix

	int nloc = _x.GetNv();

// Store initial diagonal data and compute pointers to the diagonal data

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	CSMatrixR *mtralarr, *mtrauarr;

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

	double **dptrldiag, **dptrudiag;
	double *diag;

	dptrldiag = new double * [nloc];
	if (!dptrldiag) MemoryFail (funcname);
	dptrudiag = new double * [nloc];
	if (!dptrudiag) MemoryFail (funcname);
	diag = new double [nloc];
	if (!diag) MemoryFail (funcname);

	int *listl, *list2l, *ial, *jal, *ja2l;
	double *al, *au;

	nloc = 0;

	int i, iblk, iiblk, jjblk, ii, j, jj, nsupi, isup;

	for (ii=0;ii<nlistloc;ii++) {
		iblk = listbloc[ii];
		mtralarr[iblk].GetSparsity (listl, ial, jal);
		list2l = mtralarr[iblk].GetList2 ();
		ja2l = mtralarr[iblk].GetJa2 ();
		al = mtralarr[iblk].GetA ();
		au = mtrauarr[iblk].GetA ();
		nsupi = bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
		for (i=0;i<nsupi;i++) {
			isup = listl[i];
			iiblk = list2l[i];
			for (j=ial[i];j<ial[i+1];j++) {
				jj = jal[j];
				jjblk = ja2l[j];
				if (isup == jj && iiblk == jjblk) {
					dptrldiag[nloc] = al+j;
					dptrudiag[nloc] = au+j;
					nloc++;
				};
			};
		};
	};

	double *dptr;

	for (i=0;i<nloc;i++) {
		dptr = dptrldiag[i];
		diag[i] = *dptr;
	};

// Init guess to the solution

	clock_t time0, time1;

	time0 = clock ();

	CSVector sol(nloc);

	sol = _x;

// Allocate temporary vector data

	CSVector bloc(nloc), r(nloc), p(nloc), q(nloc), z(nloc), goodsol(nloc);

// Compute initial residual and its norm

	CGSMatrixR gmtrdummy;

	gmtrdummy.MvmA2Index (_tree, _mvm, _gmtral, _gmtrau, sol, r);
	if (_second) {
		_gmtra2.MvmA2Index2 (_tree, _mvm2, sol, r);
	};

	double aux;

	if (_itype == 1) {

		int itype = 0;

		for (i=0;i<nloc;i++) {
			aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
			r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
		};

	} else {

		for (i=0;i<nloc;i++) {
			r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
		};

	};

	double errmin;

	errmin = r.MaxProdSVect (_norm);

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, MAXIMUM,
													1, &errmin, &errmin);
	};

	if (_param.msglev > 1 && myid == 0) {
		cout << " Errmin ini = " << errmin << endl;
	};
	if (_param.msglev >= 1 && myid == 0) {
		_fout << " Errmin ini = " << errmin << endl;
	};

// Init SOR data

	double invomega = 1.0e0 / _param.omega;
//	double dsclr = invomega-1.0e0;

	int ichk = 5;

// Main iterative cycle

	int k;

	for (k=1;k<_param.niter;k++) {

// Prepare triangular system for solves

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			double der;
			if (_itype == 1) {

				int itype = 1;

				der = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);

			} else {

				der = 2.0e0 * _poly2.vect[i] * sol.vect[i];

			};
			aux = diag[i] * invomega - der;
			*dptr = 1.0e0 / aux;
		};

// Solve triangular system with al

		_gmtral.SolveL2Index (_tree, _mvm,
								r, p);

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
		};

// Update solution vector

		for (i=0;i<nloc;i++) {
			sol.vect[i] += p.vect[i];
			if (sol.vect[i] < _param.permismin) sol.vect[i] = _param.permismin;
		};

// Compute new right hand side vector

		gmtrdummy.MvmA2Index (_tree, _mvm, _gmtral, _gmtrau, sol, r);
		if (_second) {
			_gmtra2.MvmA2Index2 (_tree, _mvm2, sol, r);
		};

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
			};

		};

		double errloc;

		errloc = r.MaxProdSVect (_norm);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, MAXIMUM,
														1, &errloc, &errloc);
		};

// Check convergence

		if (errloc < errmin) {
			goodsol = sol;
			errmin = errloc;
		};

		if (k%ichk == 0) {
			if (_param.msglev > 1 && myid == 0) {
				cout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
			if (_param.msglev >= 1 && myid == 0) {
				_fout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
		};

		if (errmin < _param.eps) goto exit;

		if (k > _param.niter) goto exit;

	};

exit: 

	if (k%ichk > 0) {
		if (_param.msglev > 1 && myid == 0) {
			cout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
	};
//
// Finalize time measurement
//
	time1 = clock ();

	double tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output SOR statistics

	if (_param.msglev > 1 && myid == 0) {
		cout << " SOR statistics: " << endl;
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;

// Return the solution and params

	_param.resfin = errmin;
	_param.nitperf = k;

	sol = goodsol;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel iterations of the nonlinear SOR algorithm
//========================================================================================
void CSVector::NonlinearSorSchur2Index (ofstream &_fout, const CTree &_tree, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
										void *_obja, CMVM _mvm,
										void *_objal, CSOLVEAL _solveal,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param) {

	const char *funcname = "NonlinearSorSchur2Index";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();

// Get the local size of the matrix

	int nloc = _x.GetNv();

// Store initial diagonal data and compute pointers to the diagonal data

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	CSMatrixR *mtralarr, *mtrauarr;

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

	double **dptrldiag, **dptrudiag;
	double *diag;

	dptrldiag = new double * [nloc];
	if (!dptrldiag) MemoryFail (funcname);
	dptrudiag = new double * [nloc];
	if (!dptrudiag) MemoryFail (funcname);
	diag = new double [nloc];
	if (!diag) MemoryFail (funcname);

	int *listl, *list2l, *ial, *jal, *ja2l;
	double *al, *au;

	nloc = 0;

	int i, iblk, iiblk, jjblk, ii, j, jj, nsupi, isup;

	for (ii=0;ii<nlistloc;ii++) {
		iblk = listbloc[ii];
		mtralarr[iblk].GetSparsity (listl, ial, jal);
		list2l = mtralarr[iblk].GetList2 ();
		ja2l = mtralarr[iblk].GetJa2 ();
		al = mtralarr[iblk].GetA ();
		au = mtrauarr[iblk].GetA ();
		nsupi = bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
		for (i=0;i<nsupi;i++) {
			isup = listl[i];
			iiblk = list2l[i];
			for (j=ial[i];j<ial[i+1];j++) {
				jj = jal[j];
				jjblk = ja2l[j];
				if (isup == jj && iiblk == jjblk) {
					dptrldiag[nloc] = al+j;
					dptrudiag[nloc] = au+j;
					nloc++;
				};
			};
		};
	};

	double *dptr;

	for (i=0;i<nloc;i++) {
		dptr = dptrldiag[i];
		diag[i] = *dptr;
	};

// Init guess to the solution

	clock_t time0, time1;

	time0 = clock ();

	CSVector sol(nloc);

	sol = _x;

// Allocate temporary vector data

	CSVector bloc(nloc), r(nloc), p(nloc), q(nloc), z(nloc), goodsol(nloc);

// Compute initial residual and its norm

	(_mvm) (_obja, sol, r);

	double aux;

	if (_itype == 1) {

		int itype = 0;

		for (i=0;i<nloc;i++) {
			aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
			r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
		};

	} else {

		for (i=0;i<nloc;i++) {
			r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
		};

	};

	double errmin;

	errmin = r.MaxProdSVect (_norm);

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, MAXIMUM,
													1, &errmin, &errmin);
	};

	if (_param.msglev > 1 && myid == 0) {
		cout << " Errmin ini = " << errmin << endl;
	};
	if (_param.msglev >= 1 && myid == 0) {
		_fout << " Errmin ini = " << errmin << endl;
	};

// Init SOR data

	double invomega = 1.0e0 / _param.omega;
//	double dsclr = invomega-1.0e0;

	int ichk = 5;

// Main iterative cycle

	int k;

	for (k=1;k<_param.niter;k++) {

// Prepare triangular system for solves

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			double der;
			if (_itype == 1) {

				int itype = 1;

				der = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);

			} else {

				der = 2.0e0 * _poly2.vect[i] * sol.vect[i];

			};
			aux = diag[i] * invomega - der;
			*dptr = 1.0e0 / aux;
		};

// Solve triangular system with al

		(_solveal) (_objal, r, p);

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
		};

// Update solution vector

		for (i=0;i<nloc;i++) {
			sol.vect[i] += p.vect[i];
			if (sol.vect[i] < _param.permismin) sol.vect[i] = _param.permismin;
		};

// Compute new right hand side vector

		(_mvm) (_obja, sol, r);

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, sol.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i]) - r.vect[i] ;
			};

		};

		double errloc;

		errloc = r.MaxProdSVect (_norm);

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, MAXIMUM,
														1, &errloc, &errloc);
		};

// Check convergence

		if (errloc < errmin) {
			goodsol = sol;
			errmin = errloc;
		};

		if (k%ichk == 0) {
			if (_param.msglev > 1 && myid == 0) {
				cout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
			if (_param.msglev >= 1 && myid == 0) {
				_fout << " Iter = " << k << " Errmin = " << errmin << endl;
			};
		};

		if (errmin < _param.eps) goto exit;

		if (k > _param.niter) goto exit;

	};

exit: 

	if (k%ichk > 0) {
		if (_param.msglev > 1 && myid == 0) {
			cout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Iter = " << k << " Errmin = " << errmin << endl;
		};
	};
//
// Finalize time measurement
//
	time1 = clock ();

	double tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output SOR statistics

	if (_param.msglev > 1 && myid == 0) {
		cout << " SOR statistics: " << endl;
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;

// Return the solution and params

	_param.resfin = errmin;
	_param.nitperf = k;

	sol = goodsol;

	*this = sol;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
