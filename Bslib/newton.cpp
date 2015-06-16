//------------------------------------------------------------------------------------------------
// File: newton.cpp
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
#include "mvm.h"
#include "fct.h"
#include "slvparam.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixR: Perform transformation of the data for nonlinear solver
//========================================================================================
void CSMatrixR::ConvertData (const CSMatrixR &_mtra, const CSVector &_x, // Perform transformation of the data for nonlinear solver
								const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, const bool *_wall_E,
								const CSVector &_norm,
								CSMatrixR &_mtranew, CSVector &_xnew, CSVector &_poly0new, CSVector &_poly2new, CSVector &_normnew) const { 

	const char *funcname = "ConvertData";

// Count the number of nodes, the number of nonzeroes in the matrix and store the order

	int nloc = _mtra.GetN ();

	int i, j, k, jj;

	int *order;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);

	for (i=0;i<nloc;i++) order[i] = -1;

	int nlocnew = 0;
	int nzanew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			for (j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
				jj = _mtra.ja[j];
				if (!_wall_E[jj]) nzanew++;
			};
			order[i] = nlocnew;
			nlocnew++;
		};
	};

// Allocate the matrix and init it

	CSMatrixR temp (nlocnew,nzanew);

	temp.ia[0] = 0;

	nlocnew = 0;
	nzanew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			for (j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
				jj = _mtra.ja[j];
				if (!_wall_E[jj]) {
					k = order[jj];
					temp.ja[nzanew] = k;
					temp.a[nzanew] = _mtra.a[j];
					if (i == jj) {
						temp.a[nzanew] -= _poly1.vect[i];
					};
					nzanew++;
				};
			};
			temp.ia[nlocnew+1] = nzanew;
			nlocnew++;
		};
	};

	for (i=0;i<nlocnew;i++) temp.list[i] = i;

	temp.m = nlocnew;
	temp.n = nlocnew;
	temp.nsupr = nlocnew;
	temp.nsupc = nlocnew;
	temp.nlist = nlocnew;
	temp.nzja = nzanew;
	temp.nza = nzanew;
	temp.nzatot = nzanew;

// Store the result

	_mtranew = temp;

	CSMatrixR mtrdummy;

	temp = mtrdummy;

// Reassign x, poly2 and norm

	CSVector vloc (nlocnew);

	nlocnew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			vloc.vect[nlocnew] = _x.vect[i];
			nlocnew++;
		};
	};

	_xnew = vloc;

	nlocnew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			vloc.vect[nlocnew] = _poly2.vect[i];
			nlocnew++;
		};
	};

	_poly2new = vloc;

	nlocnew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			vloc.vect[nlocnew] = _norm.vect[i];
			nlocnew++;
		};
	};

	_normnew = vloc;

// Compute poly0new

	nlocnew = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			vloc.vect[nlocnew] = _poly0.vect[i];
			for (j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
				jj = _mtra.ja[j];
				if (_wall_E[jj]) {
					vloc.vect[nlocnew] -= _mtra.a[j] * _x.vect[jj];
				};
			};
			nlocnew++;
		};
	};

	_poly0new = vloc;

// Delete working arrays

	delete [] order;

};

// Author: Kharchenko S.A.
// CSVector: Perform Newton iterations for quadratic nonlinearity in the rhs
//========================================================================================
void CSVector::NewtonSolve (ofstream &_fout, CSMatrixR &_mtra, // Perform Newton iterations for quadratic nonlinearity in the rhs
								CSVector &_x, CSVector &_poly0, CSVector &_poly2, 
								CSVector &_norm, CSlvParam &_param) {

	const char *funcname = "NewtonSolve";

	int istep, i, j, jj;

	double tottim;
	clock_t time0, time00, time1, time11;
//
// Init statistics data
//
	time0 = clock ();

// Compute reordering of the matrix

	if (_param.msglev > 0) {
		if (_param.collap != 0) _mtra.A2Ps (_param.collap,"AStrIni.ps",0,&i);
	};

	int nloc = _mtra.GetN ();

	int *order, *iord;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);
	iord = new int [nloc];
	if (!iord) MemoryFail (funcname);

	int ordtyp = _param.ordtype;

	if (ordtyp > 0) {
		_mtra.OrderPrfMtr (ordtyp, order);
	} else {
		for (i=0;i<nloc;i++) order[i] = i;
	};

	for (i=0;i<nloc;i++) iord[order[i]] = i;

// Reorder the matrix in-place

	CSMatrixR mtrao;
	CSMatrixR mtrdummy;

	if (ordtyp > 0) {

		mtrao = _mtra.OrdMtr (order);

		_mtra = mtrao;

		mtrao = mtrdummy;

		if (_param.msglev > 0) {
			if (_param.collap != 0) _mtra.A2Ps (_param.collap,"AStrOrd.ps",0,&i);
		};

	};

// Reorder all other data

	CSVector temp (nloc);

	temp.OrdVect ('D', order, _x);
	_x = temp;

	temp.OrdVect ('D', order, _poly0);
	_poly0 = temp;

	temp.OrdVect ('D', order, _poly2);
	_poly2 = temp;

	temp.OrdVect ('D', order, _norm);
	_norm = temp;

// Store initial diagonal data and compute pointers to the diagonal data

	int *iptrdiag;
	double *diag;

	iptrdiag = new int [nloc];
	if (!iptrdiag) MemoryFail (funcname);
	diag = new double [nloc];
	if (!diag) MemoryFail (funcname);

	for (i=0;i<nloc;i++) {
		for (j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
			jj = _mtra.ja[j];
			if (i == jj) {
				iptrdiag[i] = j;
				diag[i] = _mtra.a[j];
			};
		};
	};

// Set up array that specifies different Newton steps to be checked

	int nstepsmax = 20;

	double *alphaarr;

	alphaarr = new double [nstepsmax];
	if (!alphaarr) MemoryFail (funcname);

	alphaarr[0] = 1.0e0;

	for (i=1;i<nstepsmax;i++) {
		alphaarr[i] = alphaarr[i-1] / 2.0e0;
	};

// Create initial mtrl and mtru

	CSMatrixR mtrl, mtru;

// Main iterative cycle

	CSVector r(nloc), delta(nloc), sol(nloc), solloc(nloc);

	int nitertot = 0;
	double timefcttot = 0.0e0;

	bool conv;

	conv = false;

	int iter = -1;

	time00 = clock ();

	while (!conv) {

		iter++;

		if (_param.msglev > 1) {
			cout  << " Newton step = " << iter << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Newton step = " << iter << endl;
		};

// Compute the nonlinear residual vector

		_mtra.MvmA (_x,r);

		for (i=0;i<nloc;i++) {
			r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
		};

// Modify the diagonal

		for (i=0;i<nloc;i++) {
			j = iptrdiag[i];
			_mtra.a[j] = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
		};

// Compute the preconditioner if necessary

		if (iter%_param.nprec == 0) {

			mtrl = mtrdummy;
			mtru = mtrdummy;

			CSMatrixR mtral, mtrau;

			_mtra.CombineLU (mtral, mtrau);

			if (_param.memory > 0.0e0) {
				mtrdummy.Ilu2 (_fout, _param, 
								mtral, mtrau,
								mtrl, mtru);
			} else if (_param.memory == 0.0e0) {
				mtrdummy.Ilu2Dynamic (_fout, _param, 
										mtral, mtrau,
										mtrl, mtru);
			} else if (_param.memory < 0.0e0) {
				mtrdummy.Ilu2DynamicByBlocks (_fout, _param, 
										mtral, mtrau,
										mtrl, mtru);
			};

			timefcttot += _param.timefct;

			mtral = mtrdummy;
			mtrau = mtrdummy;

		};

// Solve

		if (_param.niter == 0) {

			mtrl.SolveL (r,solloc); 
			mtru.SolveU (solloc,delta);

			nitertot += 1;

		} else {

			sol.SetSVect (0.0e0);

			delta.Lanczos (_fout, (const CSMatrix *) &_mtra, (const CSMatrix *) &mtrl, (const CSMatrix *) &mtru, 
							sol, r, 
							_param, _norm);

			nitertot += _param.nitperf;

		};

// Restore the diagonal of the matrix

		for (i=0;i<nloc;i++) {
			j = iptrdiag[i];
			_mtra.a[j] = diag[i];
		};

// For different shifts compute the norm of the residual, modify the solution and check the convergence

		double value = 0.0e0;

		double valuemin = 0.0e0;

		if (_param.nsteps != 1) {

			solloc = _x;

			for (istep=0;istep<_param.nsteps;istep++) {

				sol = solloc;

				sol.DaxpySVect (alphaarr[istep], delta);

				for (i=0;i<nloc;i++) {
					if (sol.vect[i] < _param.valuemin) sol.vect[i] = _param.valuemin;
				};

				_mtra.MvmA (sol,r);

				for (i=0;i<nloc;i++) {
					r.vect[i] = _poly0.vect[i] - r.vect[i] + _poly2.vect[i] * sol.vect[i] * sol.vect[i];
				};

				value = r.MaxProdSVect(_norm);

				if (_param.msglev > 1) {
					cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};

				if (istep==0 || value < valuemin) {
					valuemin = value;
					_x = sol;
					if (valuemin <= _param.eps2) conv = true;
				};

			};

		} else {

			istep = 0;

			_x.DaxpySVect (alphaarr[istep], delta);

			for (i=0;i<nloc;i++) {
				if (_x.vect[i] < _param.valuemin) _x.vect[i] = _param.valuemin;
			};

			_mtra.MvmA (_x,r);

			for (i=0;i<nloc;i++) {
				r.vect[i] = _poly0.vect[i] - r.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i];
			};

			value = r.MaxProdSVect(_norm);

			if (_param.msglev > 1) {
				cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};

			valuemin = value;
			if (valuemin <= _param.eps2) conv = true;

		};

		_param.nitperf2 = iter;

		if (iter >= _param.niter2) conv = true;

	};

	time11 = clock ();

// Reorder the data back 

	if (ordtyp > 0) {

		mtrao = _mtra.OrdMtr (iord);

		_mtra = mtrao;

		mtrao = mtrdummy;

	};

	temp.OrdVect ('I', order, _x);
	_x = temp;

	temp.OrdVect ('I', order, _poly0);
	_poly0 = temp;

	temp.OrdVect ('I', order, _poly2);
	_poly2 = temp;

	temp.OrdVect ('I', order, _norm);
	_norm = temp;
//
// Finalize time measurement
//
	time1 = clock ();
//
// Output Newton statistics
//
	if (_param.msglev > 1) {
		cout << " Newton statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Newton statistics: " << endl;
	};

	tottim = (double) (time11-time00) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev > 1) {
		cout << "  LocTime  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "  LocTime  = " << tottim << " sec. "<< endl;
	};

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.nitperf = nitertot;
	_param.timefct = timefcttot;
	_param.timeitr = tottim - timefcttot;

	if (_param.msglev > 1) {
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] order;
	delete [] iord;
	delete [] iptrdiag;
	delete [] diag;
	delete [] alphaarr;

// Return solution

	*this = _x;

};

// Author: Kharchenko S.A.
// CSVector: Perform transformation of the data back from the nonlinear solver
//========================================================================================
void CSVector::ConvertDataBack (const CSVector &_x, const CSVector &_sol, const bool *_wall_E) { // Perform transformation of the data back from the nonlinear solver

//	const char *funcname = "ConvertDataBack";

// Count the number of nodes, the number of nonzeroes in the matrix and store the order

	int nloc = _x.GetNv ();

	CSVector solution(nloc);

	int i, j;

	j = 0;

	for (i=0;i<nloc;i++) {
		if (!_wall_E[i]) {
			solution.vect[i] = _sol.vect[j];
			j++;
		} else {
			solution.vect[i] = _x.vect[i];
		};
	};

	*this = solution;

};

// Author: Kharchenko S.A.
// CSVector: Perform Newton iterations for general nonlinearity in the rhs (with one discontinuity)
//========================================================================================
void CSVector::NewtonSolve (ofstream &_fout, CSMatrixR &_mtra, // Perform Newton iterations for general nonlinearity in the rhs (with one discontinuity)
							int _iparam, void **_objs, FUNCVALUE *_funcs,
							CSVector &_x, CSVector &_poly0, 
							CSVector &_norm, CSlvParam &_param) {

	const char *funcname = "NewtonSolve";

	int istep, i, j, jj;

	double tottim;
	clock_t time0, time00, time1, time11;
//
// Init statistics data
//
	time0 = clock ();

// Compute reordering of the matrix

	if (_param.msglev > 0) {
		if (_param.collap != 0) _mtra.A2Ps (_param.collap,"AStrIni.ps",0,&i);
	};

	int nloc = _mtra.GetN ();

	int *order, *iord;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);
	iord = new int [nloc];
	if (!iord) MemoryFail (funcname);

	int ordtyp = _param.ordtype;

	if (ordtyp > 0) {
		_mtra.OrderPrfMtr (ordtyp, order);
	} else {
		for (i=0;i<nloc;i++) order[i] = i;
	};

	for (i=0;i<nloc;i++) iord[order[i]] = i;

// Reorder the matrix in-place

	CSMatrixR mtrao;
	CSMatrixR mtrdummy;

	if (ordtyp > 0) {

		mtrao = _mtra.OrdMtr (order);

		_mtra = mtrao;

		mtrao = mtrdummy;

		if (_param.msglev > 0) {
			if (_param.collap != 0) _mtra.A2Ps (_param.collap,"AStrOrd.ps",0,&i);
		};

	};

// Reorder all other data

	CSVector temp (nloc);

	temp.OrdVect ('D', order, _x);
	_x = temp;

	temp.OrdVect ('D', order, _poly0);
	_poly0 = temp;

	temp.OrdVect ('D', order, _norm);
	_norm = temp;

	void **objsloc;

	FUNCVALUE *funcsloc;

	objsloc = new void * [nloc];
	if (!objsloc) MemoryFail (funcname);
	funcsloc = new FUNCVALUE [nloc];
	if (!funcsloc) MemoryFail (funcname);

	int inew;

	for (i=0;i<nloc;i++) {
		inew = order[i];
		objsloc [inew] = _objs[i];
		funcsloc[inew] = _funcs[i];
	};

// Store initial diagonal data and compute pointers to the diagonal data

	int *iptrdiag;
	double *diag;

	iptrdiag = new int [nloc];
	if (!iptrdiag) MemoryFail (funcname);
	diag = new double [nloc];
	if (!diag) MemoryFail (funcname);

	for (i=0;i<nloc;i++) {
		for (j=_mtra.ia[i];j<_mtra.ia[i+1];j++) {
			jj = _mtra.ja[j];
			if (i == jj) {
				iptrdiag[i] = j;
				diag[i] = _mtra.a[j];
			};
		};
	};

// Set up array that specifies different Newton steps to be checked

	int nstepsmax = 20;

	double *alphaarr;

	alphaarr = new double [nstepsmax];
	if (!alphaarr) MemoryFail (funcname);

	alphaarr[0] = 1.0e0;

	for (i=1;i<nstepsmax;i++) {
		alphaarr[i] = alphaarr[i-1] / 2.0e0;
	};

// Create initial mtrl and mtru

	CSMatrixR mtrl, mtru;

// Main iterative cycle

	CSVector r(nloc), delta(nloc), sol(nloc), solloc(nloc);

	int nitertot = 0;
	double timefcttot = 0.0e0;

	bool conv;

	conv = false;

	int iter = -1;

	time00 = clock ();

	while (!conv) {

		iter++;

		if (_param.msglev > 1) {
			cout  << " Newton step = " << iter << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Newton step = " << iter << endl;
		};

// Compute the nonlinear residual vector

		_mtra.MvmA (_x,r);

		double aux;

		int itype = 0;

		for (i=0;i<nloc;i++) {
			aux = funcsloc[i] (objsloc[i], itype, _iparam, _x.vect[i]);
			r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
		};

// Modify the diagonal

		itype = 1;

		for (i=0;i<nloc;i++) {
			j = iptrdiag[i];
			aux = funcsloc[i] (objsloc[i], itype, _iparam, _x.vect[i]);
			_mtra.a[j] = diag[i] - aux;
		};

// Compute the preconditioner if necessary

		if (iter%_param.nprec == 0) {

			mtrl = mtrdummy;
			mtru = mtrdummy;

			CSMatrixR mtral, mtrau;

			_mtra.CombineLU (mtral, mtrau);

			if (_param.memory > 0.0e0) {
				mtrdummy.Ilu2 (_fout, _param, 
								mtral, mtrau,
								mtrl, mtru);
			} else if (_param.memory == 0.0e0) {
				mtrdummy.Ilu2Dynamic (_fout, _param, 
										mtral, mtrau,
										mtrl, mtru);
			} else if (_param.memory < 0.0e0) {
				mtrdummy.Ilu2DynamicByBlocks (_fout, _param, 
										mtral, mtrau,
										mtrl, mtru);
			};

			timefcttot += _param.timefct;

			mtral = mtrdummy;
			mtrau = mtrdummy;

		};

// Solve

		if (_param.niter == 0) {

			mtrl.SolveL (r,solloc); 
			mtru.SolveU (solloc,delta);

			nitertot += 1;

		} else {

			sol.SetSVect (0.0e0);

			delta.Lanczos (_fout, (const CSMatrix *) &_mtra, (const CSMatrix *) &mtrl, (const CSMatrix *) &mtru, 
							sol, r, 
							_param, _norm);

			nitertot += _param.nitperf;

		};

// Restore the diagonal of the matrix

		for (i=0;i<nloc;i++) {
			j = iptrdiag[i];
			_mtra.a[j] = diag[i];
		};

// For different shifts compute the norm of the residual, modify the solution and check the convergence

		double value = 0.0e0;

		double valuemin = 0.0e0;

		if (_param.nsteps != 1) {

			solloc = _x;

			for (istep=0;istep<_param.nsteps;istep++) {

				sol = solloc;

				sol.DaxpySVect (alphaarr[istep], delta);

				for (i=0;i<nloc;i++) {
					if (sol.vect[i] < _param.valuemin) sol.vect[i] = _param.valuemin;
				};

				_mtra.MvmA (sol,r);

				itype = 0;

				for (i=0;i<nloc;i++) {
					aux = funcsloc[i] (objsloc[i], itype, _iparam, sol.vect[i]);
					r.vect[i] = _poly0.vect[i] - r.vect[i] + aux;
				};

				value = r.MaxProdSVect(_norm);

				if (_param.msglev > 1) {
					cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};

				if (istep==0 || value < valuemin) {
					valuemin = value;
					_x = sol;
					if (valuemin <= _param.eps2) conv = true;
				};

			};

		} else {

			istep = 0;

			_x.DaxpySVect (alphaarr[istep], delta);

			for (i=0;i<nloc;i++) {
				if (_x.vect[i] < _param.valuemin) _x.vect[i] = _param.valuemin;
			};

			_mtra.MvmA (_x,r);

			itype = 0;

			for (i=0;i<nloc;i++) {
				aux = funcsloc[i] (objsloc[i], itype, _iparam, _x.vect[i]);
				r.vect[i] = _poly0.vect[i] - r.vect[i] + aux;
			};

			value = r.MaxProdSVect(_norm);

			if (_param.msglev > 1) {
				cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};

			valuemin = value;
			if (valuemin <= _param.eps2) conv = true;

		};

		_param.nitperf2 = iter;

		if (iter >= _param.niter2) conv = true;

	};

	time11 = clock ();

// Reorder the data back 

	if (ordtyp > 0) {

		mtrao = _mtra.OrdMtr (iord);

		_mtra = mtrao;

		mtrao = mtrdummy;

	};

	temp.OrdVect ('I', order, _x);
	_x = temp;

	temp.OrdVect ('I', order, _poly0);
	_poly0 = temp;

	temp.OrdVect ('I', order, _norm);
	_norm = temp;
//
// Finalize time measurement
//
	time1 = clock ();
//
// Output Newton statistics
//
	if (_param.msglev > 1) {
		cout << " Newton statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Newton statistics: " << endl;
	};

	tottim = (double) (time11-time00) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev > 1) {
		cout << "  LocTime  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "  LocTime  = " << tottim << " sec. "<< endl;
	};

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.nitperf = nitertot;
	_param.timefct = timefcttot;
	_param.timeitr = tottim - timefcttot;

	if (_param.msglev > 1) {
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] order;
	delete [] iord;
	delete [] iptrdiag;
	delete [] diag;
	delete [] alphaarr;
	delete [] objsloc;
	delete [] funcsloc;

// Return solution

	*this = _x;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
//========================================================================================
void CSVector::NewtonSolve (ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau,
										bool _second, CMvmR &_mvm2, CGSMatrixRS &_gmtra2,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param) {

	const char *funcname = "NewtonSolve";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();
//
// Init statistics data
//
	double tottim;
	clock_t time0, time00, time1, time11;

	time0 = clock ();

	int istep, i, j, jj;

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;
	int *spndcl, *spndrl;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);
	_gmtral.GetSups (spndcl, spndrl);

	int nloc = 0;
	int iblk;

	for (i=0;i<nlistloc;i++) {
		iblk = listbloc[i];
		nloc += bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
	};

	CSMatrixRS *mtralarr, *mtrauarr;

	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

// Store initial diagonal data and compute pointers to the diagonal data

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

	int blamxloc, ii, nsupi, isup, ni, ibs, kii;

	for (ii=0;ii<nlistloc;ii++) {
		iblk = listbloc[ii];
		blamxloc = mtralarr[iblk].blamx;
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

// Set up array that specifies different Newton steps to be checked

	int nstepsmax = 20;

	double *alphaarr;

	alphaarr = new double [nstepsmax];
	if (!alphaarr) MemoryFail (funcname);

	alphaarr[0] = 1.0e0;

	for (i=1;i<nstepsmax;i++) {
		alphaarr[i] = alphaarr[i-1] / 2.0e0;
	};

// Create initial mtrl and mtru

	CGSMatrixRS gmtrl, gmtru;

	gmtrl = _gmtral;
	gmtru = _gmtrau;

// Main iterative cycle

	CSVector r, delta, sol, solloc;

	r = _x;
	delta = _x;
	sol = _x;
	solloc = _x;

	int nitertot = 0;
	double timefcttot = 0.0e0;

	bool conv;

	CGSMatrixRS gmtrdummy;

	conv = false;

	int iter = -1;

	time00 = clock ();

	while (!conv) {

		iter++;

		if (_param.msglev > 1 && myid == 0) {
			cout  << " Newton step = " << iter << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Newton step = " << iter << endl;
		};

// Compute the nonlinear residual vector

		gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, _x, r);
		if (_second) {
			_gmtra2.MvmA2 (_tree, _mvm2, _x, r);
		};

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
			};

		};

// Modify the diagonal

		if (_itype == 1) {

			int itype = 1;

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				*dptr = diag[i] - aux;
				dptr = dptrudiag[i];
				*dptr = diag[i] - aux;
			};

		} else {

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
				dptr = dptrudiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
			};

		};

// Compute the preconditioner if necessary

		if (iter%_param.nprec == 0) {

			_gmtral.Ilu2Schur (_fout, _tree, _param,
										_gmtral, _gmtrau,
										gmtrl, gmtru);

			timefcttot += _param.timefct;

		};

// Solve

		if (_param.niter == 0) {

			gmtrl.SolveL (_tree, _mvm,
								r, solloc);
			gmtru.SolveU (_tree, _mvm,
								solloc, delta);

			nitertot += 1;

		} else {

			sol.SetSVect (0.0e0);

			CParMvmA objmva ((CTree *)&_tree, &_mvm, &_gmtral, &_gmtrau,
									_second, &_mvm2, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &_mvm, &gmtrl, &gmtru);

			delta.Lanczos (_fout, _tree,
							&objmva, CParMvmA::MvmA, CParMvmA::MvmAt,
							&objslvlu, CParSolveLU::SolveL, CParSolveLU::SolveLt,
							&objslvlu, CParSolveLU::SolveU, CParSolveLU::SolveUt,
							sol, r, _norm,
							_param);

			nitertot += _param.nitperf;

		};

// Restore the diagonal of the matrix

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
			dptr = dptrudiag[i];
			*dptr = diag[i];
		};

// For different shifts compute the norm of the residual, modify the solution and check the convergence

		double value = 0.0e0;
		double value1;

		double valuemin = 0.0e0;

		if (_param.nsteps != 1) {

			solloc = _x;

			for (istep=0;istep<_param.nsteps;istep++) {

				sol = solloc;

				sol.DaxpySVect (alphaarr[istep], delta);

				for (i=0;i<nloc;i++) {
					if (sol.vect[i] < _param.valuemin) sol.vect[i] = _param.valuemin;
				};

				gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, sol, r);
				if (_second) {
					_gmtra2.MvmA2 (_tree, _mvm2, sol, r);
				};

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

				value = r.MaxProdSVect(_norm);

				if (_tree.nproc != 1) {
					CMPIComm pcomm = _tree.GetComm ();
					CMPIExchange::ExchangeArrayMPI (pcomm,
																DOUBLEVALUE, MAXIMUM,
																1, &value, &value1);
					value = value1;
				};

				if (_param.msglev > 1) {
					cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};

				if (istep==0 || value < valuemin) {
					valuemin = value;
					_x = sol;
					if (valuemin <= _param.eps2) conv = true;
				};

			};

		} else {

			istep = 0;

			_x.DaxpySVect (alphaarr[istep], delta);

			for (i=0;i<nloc;i++) {
				if (_x.vect[i] < _param.valuemin) _x.vect[i] = _param.valuemin;
			};

			gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, _x, r);
			if (_second) {
				_gmtra2.MvmA2 (_tree, _mvm2, _x, r);
			};

			if (_itype == 1) {

				int itype = 0;

				for (i=0;i<nloc;i++) {
					aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
					r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
				};

			} else {

				for (i=0;i<nloc;i++) {
					r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
				};

			};

			value = r.MaxProdSVect(_norm);

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, MAXIMUM,
															1, &value, &value1);
				value = value1;
			};

			if (_param.msglev > 1) {
				cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};

			valuemin = value;

			_param.resfin = valuemin;

			if (valuemin <= _param.eps2) conv = true;

		};

		_param.nitperf2 = iter;

		if (iter >= _param.niter2) conv = true;

	};

	time11 = clock ();
//
// Finalize time measurement
//
	time1 = clock ();
//
// Output Newton statistics
//
	if (_param.msglev > 1) {
		cout << " Newton statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Newton statistics: " << endl;
	};

	tottim = (double) (time11-time00) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev > 1) {
		cout << "  LocTime  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "  LocTime  = " << tottim << " sec. "<< endl;
	};

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.nitperf = nitertot;
	_param.timefct = timefcttot;
	_param.timeitr = tottim - timefcttot;

	if (_param.msglev > 1) {
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;
	delete [] alphaarr;

// Return solution

	*this = _x;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
//========================================================================================
void CSVector::NewtonSolve2Index (ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
										bool _second, CMvmR &_mvm2, CGSMatrixR &_gmtra2,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param) {

	const char *funcname = "NewtonSolve2Index";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();
//
// Init statistics data
//
	double tottim;
	clock_t time0, time00, time1, time11;

	time0 = clock ();

	int istep, i, j, jj;

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;
	int *spndcl, *spndrl;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);
	_gmtral.GetSups (spndcl, spndrl);

	int nloc = 0;
	int iblk;

	for (i=0;i<nlistloc;i++) {
		iblk = listbloc[i];
		nloc += bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
	};

	CSMatrixR *mtralarr, *mtrauarr;

	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

// Store initial diagonal data and compute pointers to the diagonal data

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

	int ii, nsupi, isup;
	int iiblk, jjblk;

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

// Set up array that specifies different Newton steps to be checked

	int nstepsmax = 20;

	double *alphaarr;

	alphaarr = new double [nstepsmax];
	if (!alphaarr) MemoryFail (funcname);

	alphaarr[0] = 1.0e0;

	for (i=1;i<nstepsmax;i++) {
		alphaarr[i] = alphaarr[i-1] / 2.0e0;
	};

// Create initial mtrl and mtru

	CGSMatrixR gmtrl, gmtru;

	gmtrl = _gmtral;
	gmtru = _gmtrau;

// Main iterative cycle

	CSVector r, delta, sol, solloc;

	r = _x;
	delta = _x;
	sol = _x;
	solloc = _x;

	int nitertot = 0;
	double timefcttot = 0.0e0;

	bool conv;

	CGSMatrixR gmtrdummy;

	conv = false;

	int iter = -1;

	time00 = clock ();

	while (!conv) {

		iter++;

		if (_param.msglev > 1 && myid == 0) {
			cout  << " Newton step = " << iter << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Newton step = " << iter << endl;
		};

// Compute the nonlinear residual vector

		gmtrdummy.MvmA2Index (_tree, _mvm, _gmtral, _gmtrau, _x, r);
		if (_second) {
			_gmtra2.MvmA2Index2 (_tree, _mvm2, _x, r);
		};

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
			};

		};

// Modify the diagonal

		if (_itype == 1) {

			int itype = 1;

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				*dptr = diag[i] - aux;
				dptr = dptrudiag[i];
				*dptr = diag[i] - aux;
			};

		} else {

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
				dptr = dptrudiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
			};

		};

// Compute the preconditioner if necessary

		if (iter%_param.nprec == 0) {

			_gmtral.Ilu2SchurDynamic2Index (_fout, _tree, _param,
										_gmtral, _gmtrau,
										gmtrl, gmtru);

			timefcttot += _param.timefct;

		};

// Solve

		if (_param.niter == 0) {

			gmtrl.SolveL2Index (_tree, _mvm,
								r, solloc);
			gmtru.SolveU2Index (_tree, _mvm,
								solloc, delta);

			nitertot += 1;

		} else {

			sol.SetSVect (0.0e0);

			CParMvmA2Index objmva ((CTree *)&_tree, &_mvm, &_gmtral, &_gmtrau,
									_second, &_mvm2, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &_mvm, &gmtrl, &gmtru);

			delta.Lanczos (_fout, _tree,
							&objmva, CParMvmA2Index::MvmA, CParMvmA2Index::MvmAt,
							&objslvlu, CParSolveLU2Index::SolveL, CParSolveLU2Index::SolveLt,
							&objslvlu, CParSolveLU2Index::SolveU, CParSolveLU2Index::SolveUt,
							sol, r, _norm,
							_param);

			nitertot += _param.nitperf;

		};

// Restore the diagonal of the matrix

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
			dptr = dptrudiag[i];
			*dptr = diag[i];
		};

// For different shifts compute the norm of the residual, modify the solution and check the convergence

		double value = 0.0e0;
		double value1;

		double valuemin = 0.0e0;

		if (_param.nsteps != 1) {

			solloc = _x;

			for (istep=0;istep<_param.nsteps;istep++) {

				sol = solloc;

				sol.DaxpySVect (alphaarr[istep], delta);

				for (i=0;i<nloc;i++) {
					if (sol.vect[i] < _param.valuemin) sol.vect[i] = _param.valuemin;
				};

				gmtrdummy.MvmA2Index (_tree, _mvm, _gmtral, _gmtrau, sol, r);
				if (_second) {
					_gmtra2.MvmA2Index2 (_tree, _mvm2, sol, r);
				};

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

				value = r.MaxProdSVect(_norm);

				if (_tree.nproc != 1) {
					CMPIComm pcomm = _tree.GetComm ();
					CMPIExchange::ExchangeArrayMPI (pcomm,
																DOUBLEVALUE, MAXIMUM,
																1, &value, &value1);
					value = value1;
				};

				if (_param.msglev > 1) {
					cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};

				if (istep==0 || value < valuemin) {
					valuemin = value;
					_x = sol;
					if (valuemin <= _param.eps2) conv = true;
				};

			};

		} else {

			istep = 0;

			_x.DaxpySVect (alphaarr[istep], delta);

			for (i=0;i<nloc;i++) {
				if (_x.vect[i] < _param.valuemin) _x.vect[i] = _param.valuemin;
			};

			gmtrdummy.MvmA2Index (_tree, _mvm, _gmtral, _gmtrau, _x, r);
			if (_second) {
				_gmtra2.MvmA2Index2 (_tree, _mvm2, _x, r);
			};

			if (_itype == 1) {

				int itype = 0;

				for (i=0;i<nloc;i++) {
					aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
					r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
				};

			} else {

				for (i=0;i<nloc;i++) {
					r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
				};

			};

			value = r.MaxProdSVect(_norm);

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, MAXIMUM,
															1, &value, &value1);
				value = value1;
			};

			if (_param.msglev > 1) {
				cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};

			valuemin = value;

			_param.resfin = valuemin;

			if (valuemin <= _param.eps2) conv = true;

		};

		_param.nitperf2 = iter;

		if (iter >= _param.niter2) conv = true;

	};

	time11 = clock ();
//
// Finalize time measurement
//
	time1 = clock ();
//
// Output Newton statistics
//
	if (_param.msglev > 1) {
		cout << " Newton statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Newton statistics: " << endl;
	};

	tottim = (double) (time11-time00) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev > 1) {
		cout << "  LocTime  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "  LocTime  = " << tottim << " sec. "<< endl;
	};

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.nitperf = nitertot;
	_param.timefct = timefcttot;
	_param.timeitr = tottim - timefcttot;

	if (_param.msglev > 1) {
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;
	delete [] alphaarr;

// Return solution

	*this = _x;

};

// Author: Kharchenko S.A.
// CSVector: Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
//========================================================================================
void CSVector::NewtonSolveSchur2Index (ofstream &_fout, CTree &_tree, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
													CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
													bool _second, CGSMatrixR &_gmtra2,
													int _itype, 
													int _iparam, void **_objs, FUNCVALUE *_funcs, 
													CSVector &_poly0, CSVector &_poly2, 
													CSVector &_x, CSVector &_norm, 
													CSlvParam &_param) {

	const char *funcname = "NewtonSolveSchur2Index";

	int myid = _tree.GetMyid ();
//	int nproc = _tree.GetNproc ();
//
// Init statistics data
//
	double tottim;
	clock_t time0, time00, time1, time11;

	time0 = clock ();

	int istep, i, j, jj;

	int nlistloc;
	int *listbloc;
	int *bl2ndcloc, *bl2ndrloc;
	int *spndcl, *spndrl;

	nlistloc = _gmtral.GetNlist ();
	_gmtral.GetSparsity (listbloc);
	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);
	_gmtral.GetSups (spndcl, spndrl);

	int nloc = 0;
	int iblk;

	for (i=0;i<nlistloc;i++) {
		iblk = listbloc[i];
		nloc += bl2ndcloc[iblk+1]-bl2ndcloc[iblk];
	};

	CSMatrixR *mtralarr, *mtrauarr;

	_gmtral.GetBlsiz (bl2ndcloc, bl2ndrloc);

	mtralarr = _gmtral.GetMtrarr ();
	mtrauarr = _gmtrau.GetMtrarr ();

// Store initial diagonal data and compute pointers to the diagonal data

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

	int ii, nsupi, isup;
	int iiblk, jjblk;

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

// Set up array that specifies different Newton steps to be checked

	int nstepsmax = 20;

	double *alphaarr;

	alphaarr = new double [nstepsmax];
	if (!alphaarr) MemoryFail (funcname);

	alphaarr[0] = 1.0e0;

	for (i=1;i<nstepsmax;i++) {
		alphaarr[i] = alphaarr[i-1] / 2.0e0;
	};

// Create initial mtrl and mtru

	CGSMatrixR gmtrl, gmtru;

	gmtrl = _gmtral;
	gmtru = _gmtrau;

// Prepare initial data

	CMvmR mvm (_tree, _gmtral, 1);
	CMvmR mvm2 (_second, _tree, _gmtra2, 1);

	int nblkstree = _tree.GetNblkstree ();
	int *pblkstree = _tree.GetBlkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();
	CSMatrix * pablkstr = _tree.GetAblktree ();

	CFctSchurR fctschurr (&_fout, 
									nblkstree, pblkstree, pblk2cputree,
									&_tree, pablkstr, 
									&_gmtral, &_gmtrau, 
									&gmtrl, &gmtru,
									&_param);

	CMvmSchurR mvmschur (&_fout,
								nblkstree, pblkstree, pblk2cputree,
								&_tree, pablkstr, 
								&_gmtral, &_gmtrau, 
								&gmtrl, &gmtru, 
								CSMatrixR::AbsMvmBlockColumnAl2Index , CSMatrixR::AbsMvmBlockRowAu2Index,
								CSMatrixR::AbsSolveBlockColumnL2Index, CSMatrixR::AbsSolveBlockRowU2Index,
								&_param);

	CParSchurMvmSlv2Index objmvslv (&_tree, &mvmschur, _second, &mvm2, &_gmtra2);

	void *pobjmvslv = (void *) (&objmvslv);

// Main iterative cycle

	CSVector r, delta, sol, solloc;

	r = _x;
	delta = _x;
	sol = _x;
	solloc = _x;

	int nitertot = 0;
	double timefcttot = 0.0e0;

	bool conv;

	CGSMatrixR gmtrdummy;

	conv = false;

	int iter = -1;

	time00 = clock ();

	while (!conv) {

		iter++;

		if (_param.msglev > 1 && myid == 0) {
			cout  << " Newton step = " << iter << endl;
		};
		if (_param.msglev >= 1 && myid == 0) {
			_fout << " Newton step = " << iter << endl;
		};

// Compute the nonlinear residual vector

		CParSchurMvmSlv2Index::MvmA (pobjmvslv, _x, r);

		double aux;

		if (_itype == 1) {

			int itype = 0;

			for (i=0;i<nloc;i++) {
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
			};

		} else {

			for (i=0;i<nloc;i++) {
				r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
			};

		};

// Modify the diagonal

		if (_itype == 1) {

			int itype = 1;

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
				*dptr = diag[i] - aux;
				dptr = dptrudiag[i];
				*dptr = diag[i] - aux;
			};

		} else {

			for (i=0;i<nloc;i++) {
				dptr = dptrldiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
				dptr = dptrudiag[i];
				*dptr = diag[i] - 2 * _poly2.vect[i] * _x.vect[i];
			};

		};

// Compute the preconditioner if necessary

		if (iter%_param.nprec == 0) {

			fctschurr.AbstractBlockIlu2 (_tree, *pablkstr, nblkstree, pblk2cputree);

			timefcttot += _param.timefct;

		};

// Solve

		if (_param.niter == 0) {

			CParSchurMvmSlv2Index::SolveLU (pobjmvslv, r, delta);

			nitertot += 1;

		} else {

			sol.SetSVect (0.0e0);

			delta.Lanczos (_fout, _tree, 
						pobjmvslv, CParSchurMvmSlv2Index::MvmA, CParSchurMvmSlv2Index::MvmAt,
						pobjmvslv, CParSchurMvmSlv2Index::SolveL, CParSchurMvmSlv2Index::SolveLt,
						pobjmvslv, CParSchurMvmSlv2Index::SolveU, CParSchurMvmSlv2Index::SolveUt,
						sol, r, _norm, _param);

			nitertot += _param.nitperf;

		};

// Restore the diagonal of the matrix

		for (i=0;i<nloc;i++) {
			dptr = dptrldiag[i];
			*dptr = diag[i];
			dptr = dptrudiag[i];
			*dptr = diag[i];
		};

// For different shifts compute the norm of the residual, modify the solution and check the convergence

		double value = 0.0e0;
		double value1;

		double valuemin = 0.0e0;

		if (_param.nsteps != 1) {

			solloc = _x;

			for (istep=0;istep<_param.nsteps;istep++) {

				sol = solloc;

				sol.DaxpySVect (alphaarr[istep], delta);

				for (i=0;i<nloc;i++) {
					if (sol.vect[i] < _param.valuemin) sol.vect[i] = _param.valuemin;
				};

				CParSchurMvmSlv2Index::MvmA (pobjmvslv, sol, r);

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

				value = r.MaxProdSVect(_norm);

				if (_tree.nproc != 1) {
					CMPIComm pcomm = _tree.GetComm ();
					CMPIExchange::ExchangeArrayMPI (pcomm,
																DOUBLEVALUE, MAXIMUM,
																1, &value, &value1);
					value = value1;
				};

				if (_param.msglev > 1) {
					cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
				};

				if (istep==0 || value < valuemin) {
					valuemin = value;
					_x = sol;
					if (valuemin <= _param.eps2) conv = true;
				};

			};

		} else {

			istep = 0;

			_x.DaxpySVect (alphaarr[istep], delta);

			for (i=0;i<nloc;i++) {
				if (_x.vect[i] < _param.valuemin) _x.vect[i] = _param.valuemin;
			};

			CParSchurMvmSlv2Index::MvmA (pobjmvslv, _x, r);

			if (_itype == 1) {

				int itype = 0;

				for (i=0;i<nloc;i++) {
					aux = _funcs[i] (_objs[i], itype, _iparam, _x.vect[i]);
					r.vect[i] = (_poly0.vect[i] + aux) - r.vect[i] ;
				};

			} else {

				for (i=0;i<nloc;i++) {
					r.vect[i] = (_poly0.vect[i] + _poly2.vect[i] * _x.vect[i] * _x.vect[i]) - r.vect[i] ;
				};

			};

			value = r.MaxProdSVect(_norm);

			if (_tree.nproc != 1) {
				CMPIComm pcomm = _tree.GetComm ();
				CMPIExchange::ExchangeArrayMPI (pcomm,
															DOUBLEVALUE, MAXIMUM,
															1, &value, &value1);
				value = value1;
			};

			if (_param.msglev > 1) {
				cout  << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};
			if (_param.msglev >= 1) {
				_fout << " Iter = " << iter << " Newton shift = " << alphaarr[istep] << " Resi = " << value << endl;
			};

			valuemin = value;

			_param.resfin = valuemin;

			if (valuemin <= _param.eps2) conv = true;

		};

		_param.nitperf2 = iter;

		if (iter >= _param.niter2) conv = true;

	};

	time11 = clock ();
//
// Finalize time measurement
//
	time1 = clock ();
//
// Output Newton statistics
//
	if (_param.msglev > 1) {
		cout << " Newton statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Newton statistics: " << endl;
	};

	tottim = (double) (time11-time00) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev > 1) {
		cout << "  LocTime  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "  LocTime  = " << tottim << " sec. "<< endl;
	};

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.nitperf = nitertot;
	_param.timefct = timefcttot;
	_param.timeitr = tottim - timefcttot;

	if (_param.msglev > 1) {
		cout << "     Time  = " << tottim << " sec. "<< endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec. "<< endl;
	};

// Delete working arrays

	delete [] dptrldiag;
	delete [] dptrudiag;
	delete [] diag;
	delete [] alphaarr;

// Return solution

	*this = _x;

};
