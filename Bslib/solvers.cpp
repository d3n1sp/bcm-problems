//------------------------------------------------------------------------------------------------
// File: solvers.cpp
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
#include <algorithm>

#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "fct.h"
#include "mvm.h"
#include "qrd.h"
#include "globals.h"
#include "slvparam.h"
#include "corr.h"
#ifdef __FV_2004__
#include "../Grid/TimeDebug.h"
#endif

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

class CLsqObjC
{
// Data
	int istate;     // istate is the object state: 0 - not initialized or freed; 1 - initialized
	int nord;       // nord is the size of order array
	int *orderarr;  // orderarr[nord] is the ordering array
	CGSMatrixCS gmtru; // gmtru is the U part of the preconditioner
	CCorrectorC corr;  // corr is the corrector object containing the subspaces information
public:
// Functions
// Constructors and destructor
	CLsqObjC () { // Memory allocation constructor with zero data
		istate = 0;
		nord = 0;
		orderarr = new int [nord];
	};
	~CLsqObjC () { // Destructor
		delete [] orderarr;
	};
// Get/set functions
	int GetIstate () const {return istate;}; // Get istate
	int GetNord () const {return nord;}; // Get nord
	int * GetOrderarr () const {return orderarr;}; // Get orderarr
	CGSMatrixCS * GetGmtru () {return &gmtru;}; // Get pointer to gmtru
	CCorrectorC * GetCorr () {return &corr;}; // Get pointer to corr
	void SetIstate (int _istate) {istate = _istate;}; // Set istate
// Interface functions
	void FreeObject () { // Free the object
		delete [] orderarr;
		istate = 0;
		nord = 0;
		orderarr = new int [nord];
		CGSMatrixCS gmtrdummy;
		gmtru = gmtrdummy;
		CCorrectorC corrdummy;
		corr = corrdummy;
	};
	void AllocOrder (int _nord) { // Allocate the order array
		delete [] orderarr;
		nord = _nord;
		orderarr = new int [nord];
	};
// Input/output
	friend std::ostream &operator<< (std::ostream &_stream, const CLsqObjC &_obj); // Output object
};

// Author: Kharchenko S.A.
// Create object of the LsqObjC class
//========================================================================================
void CreateLsqObjC (void *&_obj) { // Create object of the LsqObjC class

	_obj = new CLsqObjC;

};

// Author: Kharchenko S.A.
// Free data in the LsqObjC object
//========================================================================================
void FreeLsqObjC (void *_obj) { // Free data in the LsqObjC object

	CLsqObjC *pobj = (CLsqObjC *)_obj;

	pobj->FreeObject ();

};

// Author: Kharchenko S.A.
// Delete object of the LsqObjC class
//========================================================================================
void DeleteLsqObjC (void *_obj) { // Delete object of the LsqObjC class

	CLsqObjC *pobj = (CLsqObjC *)_obj;

	delete pobj;

};

// Author: Kharchenko S.A.
// Delete object of the LsqObjC class
//========================================================================================
void SolveLsqSparse (int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
						int _nrhs, dcmplx *_rhs, dcmplx *_sol, int *_iparam, double *_dparam) {

// Setup parameters

	CSlvParam param;

	param.msglev     = _iparam[0];
	param.niterlocal = _iparam[1];
	param.nblkmax    = _iparam[2];
	param.niter      = _iparam[3];

	param.alpha      = _dparam[0];
	param.thrfltmtr  = _dparam[1];
	param.dshift     = _dparam[2];
	param.tau1       = _dparam[3];
	param.tau2       = _dparam[4];
	param.eps        = _dparam[5];
	param.xmin       = _dparam[6];
	param.xmax       = _dparam[7];

	param.ooctyp = 0;
	param.collap = 0;
	param.dshiftim = 0.0e0;
	param.theta    = 0.1e0;

	ofstream ffout ("LsqSlv.log");

	void *pobj;

	CreateLsqObjC (pobj);

	int i;

	for (i=0;i<_nrhs;i++) {
		SolveLsqSparse (ffout, param, pobj,
							_m, _n, _amatr,
							1, _rhs+i*_m, _sol+i*_n);
	};

	DeleteLsqObjC (pobj);

};

// Author: Kharchenko S.A.
// Solve complex LSQ problem
//========================================================================================
void LSQDriver (int _slvtype, double _alpha, // Solve complex LSQ problem
						int _m, int _n, dcmplx *_amatr,
						int _nrhs, dcmplx *_rhs, dcmplx *_sol) {

	const char *funcname = "LSQDriver";

// Switch the cases

	if (_slvtype == 1) {

// Qrd, no scaling

		SolveLsqQrd (_m, _n, _amatr,
							_nrhs, _rhs, _sol);

	} else if (_slvtype == 2) {

// Svd, no scaling

		SolveLsqSvd (_alpha, _m, _n, _amatr,
							_nrhs, _rhs, _sol);

	} else if (_slvtype == 3 || _slvtype == 4) {

// Compute column scaling, scale explicitely, solve and rescale the solution back

		double *sccol, *sccolinv;

		sccol = new double [_n];
		if (!sccol) MemoryFail (funcname);
		sccolinv = new double [_n];
		if (!sccolinv) MemoryFail (funcname);

		int i, j, ind;
		double aux;
		dcmplx caux;

		for (i=0;i<_n;i++) {
			aux = 0.0e0;
			for (j=0;j<_m;j++) {
				ind = i*_m+j;
				aux += (_amatr[ind].x*_amatr[ind].x + _amatr[ind].y*_amatr[ind].y);
			};
			aux = sqrt(aux);
			sccolinv[i] = aux;
			sccol[i] = 1.0e0 / sccolinv[i];
		};

		CHyst hyst;

		hyst.UpdateHystogram (_n, sccolinv);

		cout << " Column Norms: " << hyst << endl;

		for (i=0;i<_m;i++) {
			for (j=0;j<_n;j++) {
				ind = j*_m+i;
				caux = _amatr[ind] * sccol[j];
				_amatr[ind] = caux;
			};
		};

		if (_slvtype == 3) {

			SolveLsqQrd (_m, _n, _amatr,
							_nrhs, _rhs, _sol);

		} else if (_slvtype == 4) {

			SolveLsqSvd (_alpha, _m, _n, _amatr,
								_nrhs, _rhs, _sol);

		};

		for (i=0;i<_n;i++) {
			for (j=0;j<_nrhs;j++) {
				ind = j*_n+i;
				caux = _sol[ind] * sccol[i];
				_sol[ind] = caux;
			};
		};

		delete [] sccol;
		delete [] sccolinv;

	} else if (_slvtype == 5 || _slvtype == 6) {

// Compute rows and column scaling, scale explicitely, solve and rescale the solution back

		double *scrow, *scrowinv;
		double *sccol, *sccolinv;

		scrow = new double [_m];
		if (!scrow) MemoryFail (funcname);
		scrowinv = new double [_m];
		if (!scrowinv) MemoryFail (funcname);
		sccol = new double [_n];
		if (!sccol) MemoryFail (funcname);
		sccolinv = new double [_n];
		if (!sccolinv) MemoryFail (funcname);

		int i, j, ind;
		double aux;
		dcmplx caux;

		for (i=0;i<_m;i++) {
			aux = 0.0e0;
			for (j=0;j<_n;j++) {
				ind = j*_m+i;
				aux += (_amatr[ind].x*_amatr[ind].x + _amatr[ind].y*_amatr[ind].y);
			};
			aux = sqrt(aux);
			scrowinv[i] = aux;
			scrow[i] = 1.0e0 / sqrt(scrowinv[i]);
		};

		CHyst hyst;

		hyst.UpdateHystogram (_m, scrowinv);

		cout << " Rows Norms: " << hyst << endl;

		for (i=0;i<_m;i++) {
			for (j=0;j<_n;j++) {
				ind = j*_m+i;
				caux = _amatr[ind] * scrow[i];
				_amatr[ind] = caux;
			};
		};

		for (i=0;i<_m;i++) {
			for (j=0;j<_nrhs;j++) {
				ind = j*_m+i;
				caux = _rhs[ind] * scrow[i];
				_rhs[ind] = caux;
			};
		};

		for (i=0;i<_n;i++) {
			aux = 0.0e0;
			for (j=0;j<_m;j++) {
				ind = i*_m+j;
				aux += (_amatr[ind].x*_amatr[ind].x + _amatr[ind].y*_amatr[ind].y);
			};
			aux = sqrt(aux);
			sccolinv[i] = aux;
			sccol[i] = 1.0e0 / sccolinv[i];
		};

		CHyst hyst1;

		hyst1.UpdateHystogram (_n, sccolinv);

		cout << " Column Norms: " << hyst1 << endl;

		for (i=0;i<_m;i++) {
			for (j=0;j<_n;j++) {
				ind = j*_m+i;
				caux = _amatr[ind] * sccol[j];
				_amatr[ind] = caux;
			};
		};

		if (_slvtype == 5) {

			SolveLsqQrd (_m, _n, _amatr,
							_nrhs, _rhs, _sol);

		} else if (_slvtype == 6) {

			SolveLsqSvd (_alpha, _m, _n, _amatr,
								_nrhs, _rhs, _sol);

		};

		for (i=0;i<_n;i++) {
			for (j=0;j<_nrhs;j++) {
				ind = j*_n+i;
				caux = _sol[ind] * sccol[i];
				_sol[ind] = caux;
			};
		};

		delete [] scrow;
		delete [] scrowinv;
		delete [] sccol;
		delete [] sccolinv;

	};

};

// Author: Kharchenko S.A.
// Solve complex LSQ problem via QRD
//========================================================================================
void SolveLsqQrd (int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via QRD
						int _nrhs, dcmplx *_rhs, dcmplx *_sol) {

	const char *funcname = "SolveLsqQrd";

// Qr

	dcmplx *tau;

	tau = new dcmplx [_n+1];
	if (!tau) MemoryFail (funcname);

	QrdBlk (_m, _n, _amatr, _m, tau);

	dcmplx *rmatr;
	double *diag;

	rmatr = new dcmplx [_n*_n];
	if (!rmatr) MemoryFail (funcname);
	diag = new double [_n];
	if (!diag) MemoryFail (funcname);

	GetRPartQrd (_m, 0, _n,
						_amatr, _m, rmatr, _n);

	int i, ind;
	double aux;

	for (i=0;i<_n;i++) {
		ind = i*_n+i;
		aux = rmatr[ind].x*rmatr[ind].x + rmatr[ind].y*rmatr[ind].y;
		diag[i] = sqrt(aux);
	};

	CHyst hyst;

	hyst.UpdateHystogram (_n, diag);

	cout << " Diagonal R values: " << hyst << endl;

	dcmplx *solloc;

	solloc = new dcmplx [_m*_nrhs];
	if (!solloc) MemoryFail (funcname);

	MvmQHBlk (_m, _n, _nrhs,
					_amatr, _m, tau, 
					_rhs, _m, solloc, _m);

// Solve triangular system

	int j, k;

	for (k=0;k<_nrhs;k++) {
		for (i=0;i<_n;i++) {
			_sol[k*_n+i] = solloc[k*_m+i];
		};
		for (i=_n-1;i>=0;i--) {
			_sol[k*_n+i] = _sol[k*_n+i] / rmatr[i*_n+i];
			for (j=0;j<i;j++) {
				_sol[k*_n+j] -= _sol[k*_n+i] * rmatr[i*_n+j];
			};
		};
	};

// Free work arrays

	delete [] tau;
	delete [] rmatr;
	delete [] diag;
	delete [] solloc;

};

// Author: Kharchenko S.A.
// Solve complex LSQ problem via SVD
//========================================================================================
void SolveLsqSvd (double _alpha, int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via SVD
						int _nrhs, dcmplx *_rhs, dcmplx *_sol) {

	const char *funcname = "SolveLsqSvd";

// Svd

	dcmplx *tau;

	tau = new dcmplx [_n];
	if (!tau) MemoryFail (funcname);

	QrdBlk (_m, _n, _amatr, _m, tau);

	dcmplx *rmatr;

	rmatr = new dcmplx [_n*_n];
	if (!rmatr) MemoryFail (funcname);

	GetRPartQrd (_m, 0, _n,
						_amatr, _m, rmatr, _n);

	dcmplx *solloc;

	solloc = new dcmplx [_m*_nrhs];
	if (!solloc) MemoryFail (funcname);

	MvmQHBlk (_m, _n, _nrhs,
					_amatr, _m, tau, 
					_rhs, _m, solloc, _m);

	SolveDenseBySVDTikhonov (_alpha, _n, _nrhs,
										rmatr, _n,
										solloc, _m,
										_sol, _n);

// Free work arrays

	delete [] tau;
	delete [] rmatr;
	delete [] solloc;

};

// Author: Kharchenko S.A.
// Compute rectangular SVD
//========================================================================================
void MyDgesvd (int _m, int _n, double *_amatr, // Compute rectangular SVD
					double *_sv, double *_u, double *_v) {

	const char *funcname = "MyDgesvd";

// Svd

	int lwork = 10 * _n;

	double *aloc;
	double *rmatr;
	double *tau;
	double *uloc;
	double *u1loc;
	double *cwork;
	double *rwork;

	aloc = new double [_m*_n];
	if (!aloc) MemoryFail (funcname);
	rmatr = new double [_n*_n];
	if (!rmatr) MemoryFail (funcname);
	tau = new double [_n];
	if (!tau) MemoryFail (funcname);
	uloc = new double [_n*_n];
	if (!uloc) MemoryFail (funcname);
	u1loc = new double [_m*_n];
	if (!u1loc) MemoryFail (funcname);
	cwork = new double [lwork];
	if (!cwork) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

	int i, j;

	for (i=0;i<_m*_n;i++) aloc[i] = _amatr[i];

	QrdBlk (_m, _n, aloc, _m, tau);

	GetRPartQrd (_m, 0, _n,
						aloc, _m, rmatr, _n);

	int info;

	dgesvd_ ("O", "A", &_n, &_n,
				rmatr, &_n, _sv, 
				uloc, &_n, _v, &_n, 
				cwork, &lwork, &info);

	for (i=0;i<_m*_n;i++) u1loc[i] = 0.0e0;

	for (i=0;i<_n;i++) {
		for (j=0;j<_n;j++) {
			u1loc[j*_m+i] = uloc[j*_n+i];
		};
	};

	MvmQBlk (_m, _n, _n,
					aloc, _m, tau, 
					u1loc, _m, _u, _m);

// Free work arrays

	delete [] aloc;
	delete [] rmatr;
	delete [] tau;
	delete [] uloc;
	delete [] u1loc;
	delete [] cwork;
	delete [] rwork;

};

// Author: Kharchenko S.A.
// Compute rectangular SVD
//========================================================================================
void MyZgesvd (int _m, int _n, dcmplx *_amatr, // Compute rectangular SVD
					double *_sv, dcmplx *_u, dcmplx *_v) {

	const char *funcname = "MyZgesvd";

	dcmplx czero (0.0e0,0.0e0);

// Svd

	int lwork = 10 * _n;

	dcmplx *aloc;
	dcmplx *rmatr;
	dcmplx *tau;
	dcmplx *uloc;
	dcmplx *u1loc;
	dcmplx *cwork;
	double *rwork;

	aloc = new dcmplx [_m*_n];
	if (!aloc) MemoryFail (funcname);
	rmatr = new dcmplx [_n*_n];
	if (!rmatr) MemoryFail (funcname);
	tau = new dcmplx [_n];
	if (!tau) MemoryFail (funcname);
	uloc = new dcmplx [_n*_n];
	if (!uloc) MemoryFail (funcname);
	u1loc = new dcmplx [_m*_n];
	if (!u1loc) MemoryFail (funcname);
	cwork = new dcmplx [lwork];
	if (!cwork) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

	int i, j;

	for (i=0;i<_m*_n;i++) aloc[i] = _amatr[i];

	QrdBlk (_m, _n, aloc, _m, tau);

	GetRPartQrd (_m, 0, _n,
						aloc, _m, rmatr, _n);

	int info;

	zgesvd_ ("A", "A", &_n, &_n,
				rmatr, &_n, _sv, 
				uloc, &_n, _v, &_n, 
				cwork, &lwork, rwork, &info);

	for (i=0;i<_m*_n;i++) u1loc[i] = czero;

	for (i=0;i<_n;i++) {
		for (j=0;j<_n;j++) {
			u1loc[j*_m+i] = uloc[j*_n+i];
		};
	};

	MvmQBlk (_m, _n, _n,
					aloc, _m, tau, 
					u1loc, _m, _u, _m);

// Free work arrays

	delete [] aloc;
	delete [] rmatr;
	delete [] tau;
	delete [] uloc;
	delete [] u1loc;
	delete [] cwork;
	delete [] rwork;

};

// Author: Kharchenko S.A.
// Solve complex LSQ problem via sparsification
//========================================================================================
void SolveLsqSparse_old (ofstream &_fout, CSlvParam &_param,
							int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
							int _nrhs, dcmplx *_rhs, dcmplx *_sol) {

	const char *funcname = "SolveLsqSparse";

	if (_m < _n) throw " SolveLsqSparse: error, M < N ! ";
	if (_nrhs > 1) throw " SolveLsqSparse: nrhs > 1 is not implemented yet ";

// Take parameters

	double alpha = _param.alpha;
	double thresh_col = _param.thrfltmtr;

// Compute sparsified matrix

	CSMatrixC ac;

	ac = CSMatrixC::SparsifyRectangularMatrix (alpha, thresh_col,
																_m, _n, _amatr);

	if (_param.msglev > 0) {
		double dnsty = ((double) ac.GetNza ()) / ((double) _m * (double) _n);
		dnsty *= 1.0e2;
		cout  << " Density of the sparse approximation = " << dnsty << " % " << endl;
		_fout << " Density of the sparse approximation = " << dnsty << " % " << endl;
	};

//	_fout << " Ac = " << ac << endl;
/*
	if (false) {
		int *pia, *pja;
		dcmplx *pa;
		pia = ac.GetIa ();
		pja = ac.GetJa ();
		pa = ac.GetA ();
		double diff = 0.0e0;
		dcmplx caux1;
		int i, j, jj;
		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				if (jj < _m) {
					caux1 = _amatr[i*_m+jj]-pa[j];
				} else {
					caux1 = pa[j];
				};
				diff += caux1.x*caux1.x + caux1.y*caux1.y;
			};
		};
		diff = sqrt(diff);
		cout << " Diff in matrix = " << diff << endl;
	};
*/
//	ac.A2Ps ("Spars.ps",0,&_m);

// Transform C into CS

	int nsupmax = _m+_n;
	int *sprndrloc;

	sprndrloc = new int [nsupmax+1];
	if (!sprndrloc) MemoryFail (funcname);

	int i;

	for (i=0;i<=nsupmax;i++) sprndrloc[i] = i;

	CSMatrixCS acs;

	acs = CSMatrixCS::C2CSRectCol (ac, nsupmax, sprndrloc);
/*
	if (false) {
		int *pia, *pja;
		dcmplx *pa;
		pia = acs.GetIa ();
		pja = acs.GetJa ();
		pa = acs.GetA ();
		double diff = 0.0e0;
		dcmplx caux1;
		int i, j, jj;
		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				if (jj < _m) {
					caux1 = _amatr[i*_m+jj]-pa[j];
				} else {
					caux1 = pa[j];
				};
				diff += caux1.x*caux1.x + caux1.y*caux1.y;
			};
		};
		diff = sqrt(diff);
		cout << " Diff in matrixCS = " << diff << endl;
	};
*/
//	_fout << " Acs = " << acs << endl;

	CSMatrixC mtrdummy;

	ac = mtrdummy;

	delete [] sprndrloc;

// Compute AhA

	CSMatrixCS ahacs;

	ahacs = acs.AtAMatrElem (_fout, _param);

	if (_param.msglev > 0) {
		double dnsty = 2.0e0 * ((double) ahacs.GetNza ()) / ((double) _n * (double) _n);
		dnsty *= 1.0e2;
		cout  << " Density of the triangular sparse AhA approximation = " << dnsty << " % " << endl;
		_fout << " Density of the triangular sparse AhA approximation = " << dnsty << " % " << endl;
	};

/*
	if (false) {
		dcmplx *aha;
		aha = new dcmplx [_n*_n];
		int i, j, jj, k, ind1, ind2;
		dcmplx caux1;
		dcmplx czero (0.0e0,0.0e0);
		for (i=0;i<_n;i++) {
			for (j=0;j<_n;j++) {
				caux1 = czero;
				for (k=0;k<_m;k++) {
					ind1 = i*_m+k;
					ind2 = j*_m+k;
					caux1 += conj(_amatr[ind1]) * _amatr[ind2];
				};
				aha[j*_n+i] = caux1;
			};
		};
		int *pia, *pja;
		dcmplx *pa;
		pia = ahacs.GetIa ();
		pja = ahacs.GetJa ();
		pa = ahacs.GetA ();
		double diff = 0.0e0;
		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				caux1 = aha[i*_n+jj]-pa[j];
				diff += caux1.x*caux1.x + caux1.y*caux1.y;
			};
		};
		diff = sqrt(diff);
		cout << " Diff in Aha = " << diff << endl;
	};
*/
	CSMatrixCS mtrdummycs;

	acs = mtrdummycs;

//	_fout << " Ahacs = " << ahacs << endl;
//	ahacs.A2Ps ("AhA.ps",0,&_m);

// Compute norms of the columns

	double *colmnnorms2;

	colmnnorms2 = new double [_n];
	if (!colmnnorms2) MemoryFail (funcname);

	int j, ind;
	double alpha_2 = alpha*alpha;
	double aux;

	for (i=0;i<_n;i++) {
		aux = alpha_2;
		for (j=0;j<_m;j++) {
			ind = i*_m+j;
			aux += _amatr[ind].x*_amatr[ind].x + _amatr[ind].y*_amatr[ind].y;
		};
		colmnnorms2[i] = aux;
	};

// Perform correction of the diagonal values

	int *pia, *pja;
	dcmplx *pa;

	pia = ahacs.GetIa ();
	pja = ahacs.GetJa ();
	pa = ahacs.GetA ();

	int jj;

	for (i=0;i<_n;i++) {
		for (j=pia[i];j<pia[i+1];j++) {
			jj = pja[j];
			if (i == jj) {
//				pa[j].x = colmnnorms2[i];
//				pa[j].y = 0.0e0;
			};
		};
	};

	delete [] colmnnorms2;
/*
	if (false) {
		dcmplx *aha;
		aha = new dcmplx [_n*_n];
		int i, j, jj, k, ind1, ind2;
		dcmplx caux1;
		dcmplx czero (0.0e0,0.0e0);
		for (i=0;i<_n;i++) {
			for (j=0;j<_n;j++) {
				caux1 = czero;
				for (k=0;k<_m;k++) {
					ind1 = i*_m+k;
					ind2 = j*_m+k;
					caux1 += conj(_amatr[ind1]) * _amatr[ind2];
				};
				aha[j*_n+i] = caux1;
			};
		};
		int *pia, *pja;
		dcmplx *pa;
		pia = ahacs.GetIa ();
		pja = ahacs.GetJa ();
		pa = ahacs.GetA ();
		double diff = 0.0e0;
		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				caux1 = aha[i*_n+jj]-pa[j];
				diff += caux1.x*caux1.x + caux1.y*caux1.y;
			};
		};
		diff = sqrt(diff);
		cout << " Diff in Aha after diagonal shift = " << diff << endl;
	};
*/
// Compute new ordering and reorder in-place

	int ordtype = 1;

	int *order;

	order = new int [_n];
	if (!order) MemoryFail (funcname);

	CSMatrixCS ahaocs;

	if (_param.collap > 0) {
		ahacs.A2Ps (_param.collap,"AhASpIni.ps",0,&i);
	};

	ahacs.OrderPrfMtr (ordtype, order);
	ahaocs = ahacs.OrdMtrSymm (order);

	ahacs = mtrdummycs;

	if (_param.collap > 0) {
		ahaocs.A2Ps (_param.collap,"AhASpOrd.ps",0,&i);
	};

//	ahaocs.A2Ps ("AhAOrd.ps",0,&_m);

// Perform incomplete factorization

	int nblksloc = 1;
	int blksloc[2];

	blksloc[0] = 0;
	blksloc[1] = _n;

	CGSMatrixCS gahaocs;

	gahaocs.CS2GCS (ahaocs, nblksloc, blksloc);

	ahaocs = mtrdummycs;

	CGSMatrixCS gmtru;

	CSMatrixCS *pmtrarr = gahaocs.GetMtrarr ();

	int nzaha = pmtrarr->GetNza ();
	dcmplx *paha = pmtrarr->GetA ();

	dcmplx caux;

	for (i=0;i<nzaha;i++) {
		caux = paha[i];
		paha[i] = conj(caux);
	};

	CGSMatrixCS::Ich2 (_fout, _param,
								gahaocs,
								gmtru);

	if (_param.collap > 0) {
		CSMatrixCS *pmtr = gmtru.GetMtrarr ();
		pmtr->A2Ps (_param.collap,"UFct.ps",0,&i);
	};

	if (_param.msglev > 0) {
		CSMatrixCS *pmtr = gmtru.GetMtrarr ();
		double dnsty = 2.0e0 * ((double) pmtr->GetNza ()) / ((double) _n * (double) _n);
		dnsty *= 1.0e2;
		cout  << " Density of the preconditioner = " << dnsty << " % " << endl;
		_fout << " Density of the preconditioner = " << dnsty << " % " << endl;
	};

/*
	if (false) {
		int i, j, jj, k, ind1, ind2;
		dcmplx caux1;
		dcmplx czero (0.0e0,0.0e0);
		dcmplx cone  (1.0e0,0.0e0);
		dcmplx *umatr;
		umatr = new dcmplx [_n*_n];
		for (i=0;i<_n*_n;i++) umatr[i] = czero;
		CSMatrixCS *pmtr = gmtru.GetMtrarr ();
		int *pia, *pja;
		dcmplx *pa;
		pia = pmtr->GetIa ();
		pja = pmtr->GetJa ();
		pa = pmtr->GetA ();
		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				caux1 = pa[j];
				ind = jj*_n+i;
				if (i == jj) {
					umatr[ind] = cone / caux1;
				} else {
					umatr[ind] = caux1;
				};
			};
		};
		dcmplx *aha;
		aha = new dcmplx [_n*_n];
		for (i=0;i<_n;i++) {
			for (j=0;j<_n;j++) {
				caux1 = czero;
				for (k=0;k<_m;k++) {
					ind1 = i*_m+k;
					ind2 = j*_m+k;
					caux1 += conj(_amatr[ind1]) * _amatr[ind2];
				};
				aha[j*_n+i] = caux1;
			};
		};
		dcmplx *uhu;
		uhu = new dcmplx [_n*_n];
		double diff = 0.0e0;
		for (i=0;i<_n;i++) {
			for (j=0;j<_n;j++) {
				caux1 = czero;
				for (k=0;k<_n;k++) {
					ind1 = i*_n+k;
					ind2 = j*_n+k;
					caux1 += conj(umatr[ind1]) * umatr[ind2];
				};
				uhu[j*_n+i] = caux1;
				caux1 = uhu[j*_n+i] - aha[j*_n+i];
				diff += caux1.x*caux1.x + caux1.y*caux1.y;
			};
		};
		diff = sqrt(diff);
		cout << " Diff in Aha after factorization = " << diff << endl;
		OutArr (_fout," Aha via dense ",_n*_n,aha);
		OutArr (_fout," Uhu via dense ",_n*_n,uhu);
	};
*/
	CGSMatrixCS mtrdummygcs;

	gahaocs = mtrdummygcs;

//	_fout << " Prec = " << gmtru << endl;

// Create the tree

	int nchilds = 3;

	double perf2cpu[10];
	double memory2cpu[10];

	for (i=0;i<10;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (1,nchilds,perf2cpu,memory2cpu);

// Create mv and solve objects

	CParMvmALsqDenseC aobj (_m, _n, alpha, _amatr);
	CParSolveIchC luobj (&tree, order, &gmtru);

// Create rhs and initial guess

	dcmplx czero (0.0e0,0.0e0);

	int mn = _m+_n;

	CSVectorC x(_n), sol(_n), rhs(mn);

	dcmplx *px = x.GetVect ();
	dcmplx *prhs = rhs.GetVect ();

	rhs.SetSVect (czero);

	for (i=0;i<_m;i++) prhs[i] = _rhs[i];
	for (i=0;i<_n;i++) px[i] = _sol[i];
/*
	if (false) {
		CSVectorC x1(_n), x2(_n);
		CSVectorC b1(mn), b2(mn);
		CParMvmALsqDenseC::MvmAHC ((void *)&aobj,rhs,x1);
		CParSolveIchC::SolveUH ((void *)&luobj,x1,x2);
		CParSolveIchC::SolveU ((void *)&luobj,x2,x1);
		_fout << " Sol = " << x1 << endl;
	};
*/
// Solve by the iterative scheme

//	sol.BlockLanczosRight (_fout, // Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
//									tree, 
//									(void *)&aobj, CParMvmALsqDenseC::MvmAC, CParMvmALsqDenseC::MvmAHC,
//									(void *)&luobj, CParSolveIchC::SolveUH, CParSolveIchC::SolveU,
//									x, rhs, 
//									_param);

	sol.SOFLanczosNE_old (_fout, _param, // Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
								tree, 
								(void *)&aobj, CParMvmALsqDenseC::MvmAC, CParMvmALsqDenseC::MvmAHC,
								(void *)&luobj, CParSolveIchC::SolveUH, CParSolveIchC::SolveU,
								x, rhs);

// Store the solution

	dcmplx *psol = sol.GetVect ();

	for (i=0;i<_n;i++) _sol[i] = psol[i];

// Free work arrays

	delete [] order;

};

// Author: Kharchenko S.A.
// Solve complex LSQ problem via sparsification
//========================================================================================
void SolveLsqSparse (ofstream &_fout, CSlvParam &_param, void *_obj,
							int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
							int _nrhs, dcmplx *_rhs, dcmplx *_sol) {

	const char *funcname = "SolveLsqSparse";

	if (_m < _n) throw " SolveLsqSparse: error, M < N ! ";
	if (_nrhs > 1) throw " SolveLsqSparse: nrhs > 1 is not implemented yet ";

// Take parameters

	double alpha = _param.alpha;
	double thresh_col = _param.thrfltmtr;

// Get the object

	CLsqObjC *pobj = (CLsqObjC *) _obj;

	int istateloc = pobj->GetIstate ();

// Switch the cases

	if (istateloc == 0) {

// Compute sparsified matrix

		CSMatrixC ac;

		ac = CSMatrixC::SparsifyRectangularMatrix (alpha, thresh_col,
																	_m, _n, _amatr);

		if (_param.msglev > 0) {
			double dnsty = ((double) ac.GetNza ()) / ((double) _m * (double) _n);
			dnsty *= 1.0e2;
			cout  << " Density of the sparse approximation = " << dnsty << " % " << endl;
			_fout << " Density of the sparse approximation = " << dnsty << " % " << endl;
		};

// Transform C into CS

		int nsupmax = _m+_n;
		int *sprndrloc;

		sprndrloc = new int [nsupmax+1];
		if (!sprndrloc) MemoryFail (funcname);

		int i;

		for (i=0;i<=nsupmax;i++) sprndrloc[i] = i;

		CSMatrixCS acs;

		acs = CSMatrixCS::C2CSRectCol (ac, nsupmax, sprndrloc);

		CSMatrixC mtrdummy;

		ac = mtrdummy;

		delete [] sprndrloc;

// Compute AhA

		CSMatrixCS ahacs;

		ahacs = acs.AtAMatrElem (_fout, _param);

		if (_param.msglev > 0) {
			double dnsty = 2.0e0 * ((double) ahacs.GetNza ()) / ((double) _n * (double) _n);
			dnsty *= 1.0e2;
			cout  << " Density of the triangular sparse AhA approximation = " << dnsty << " % " << endl;
			_fout << " Density of the triangular sparse AhA approximation = " << dnsty << " % " << endl;
		};

		CSMatrixCS mtrdummycs;

		acs = mtrdummycs;

//	_fout << " Ahacs = " << ahacs << endl;
//	ahacs.A2Ps ("AhA.ps",0,&_m);

// Compute norms of the columns

		double *colmnnorms2;

		colmnnorms2 = new double [_n];
		if (!colmnnorms2) MemoryFail (funcname);

		int j, ind;
		double alpha_2 = alpha*alpha;
		double aux;

		for (i=0;i<_n;i++) {
			aux = alpha_2;
			for (j=0;j<_m;j++) {
				ind = i*_m+j;
				aux += _amatr[ind].x*_amatr[ind].x + _amatr[ind].y*_amatr[ind].y;
			};
			colmnnorms2[i] = aux;
		};

// Perform correction of the diagonal values

		int *pia, *pja;
		dcmplx *pa;

		pia = ahacs.GetIa ();
		pja = ahacs.GetJa ();
		pa = ahacs.GetA ();

		int jj;

		for (i=0;i<_n;i++) {
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				if (i == jj) {
	//				pa[j].x = colmnnorms2[i];
	//				pa[j].y = 0.0e0;
				};
			};
		};

		delete [] colmnnorms2;

// Compute new ordering and reorder in-place

		pobj->AllocOrder (_n);

		int *order = pobj->GetOrderarr ();

		CSMatrixCS ahaocs;

		if (_param.collap > 0) {
			ahacs.A2Ps (_param.collap,"AhASpIni.ps",0,&i);
		};

		int ordtype = 1;

		ahacs.OrderPrfMtr (ordtype, order);
		ahaocs = ahacs.OrdMtrSymm (order);

		ahacs = mtrdummycs;

		if (_param.collap > 0) {
			ahaocs.A2Ps (_param.collap,"AhASpOrd.ps",0,&i);
		};

// Perform incomplete factorization

		int nblksloc = 1;
		int blksloc[2];

		blksloc[0] = 0;
		blksloc[1] = _n;

		CGSMatrixCS gahaocs;

		gahaocs.CS2GCS (ahaocs, nblksloc, blksloc);

		ahaocs = mtrdummycs;

		CGSMatrixCS *pgmtru = pobj->GetGmtru ();

		CSMatrixCS *pmtrarr = gahaocs.GetMtrarr ();

		int nzaha = pmtrarr->GetNza ();
		dcmplx *paha = pmtrarr->GetA ();

		dcmplx caux;

		for (i=0;i<nzaha;i++) {
			caux = paha[i];
			paha[i] = conj(caux);
		};

		CGSMatrixCS::Ich2 (_fout, _param,
									gahaocs,
									*pgmtru);

		if (_param.collap > 0) {
			CSMatrixCS *pmtr = pgmtru->GetMtrarr ();
			pmtr->A2Ps (_param.collap,"UFct.ps",0,&i);
		};

		if (_param.msglev > 0) {
			CSMatrixCS *pmtr = pgmtru->GetMtrarr ();
			double dnsty = 2.0e0 * ((double) pmtr->GetNza ()) / ((double) _n * (double) _n);
			dnsty *= 1.0e2;
			cout  << " Density of the preconditioner = " << dnsty << " % " << endl;
			_fout << " Density of the preconditioner = " << dnsty << " % " << endl;
		};

		CGSMatrixCS mtrdummygcs;

		gahaocs = mtrdummygcs;

		pobj->SetIstate (1);

	};

// Get data from the object

	int nloc = pobj->GetNord ();
	int *order = pobj->GetOrderarr ();
	CGSMatrixCS *pgmtru = pobj->GetGmtru ();

	if (nloc != _n) throw " SolveLsqSparse: wrong NOrd parameter ";

// Create the tree

	int nchilds = 3;

	double perf2cpu[10];
	double memory2cpu[10];

	int i;

	for (i=0;i<10;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (1,nchilds,perf2cpu,memory2cpu);

// Create mv and solve objects

	CParMvmALsqDenseC aobj (_m, _n, alpha, _amatr);
	CParSolveIchC luobj (&tree, order, pgmtru);

// Create rhs and initial guess

	dcmplx czero (0.0e0,0.0e0);

	int mn = _m+_n;

	CSVectorC x(_n), sol(_n), rhs(mn);

	dcmplx *px = x.GetVect ();
	dcmplx *prhs = rhs.GetVect ();

	rhs.SetSVect (czero);

	for (i=0;i<_m;i++) prhs[i] = _rhs[i];
	for (i=0;i<_n;i++) px[i] = _sol[i];

	CCorrectorC *pcorr = pobj->GetCorr ();

	sol.SOFLanczosNE (_fout, _param, 
									tree, *pcorr,
									(void *)&aobj, CParMvmALsqDenseC::MvmAC, CParMvmALsqDenseC::MvmAHC,
									(void *)&luobj, CParSolveIchC::SolveUH, CParSolveIchC::SolveU,
									x, rhs);

// Store the solution

	dcmplx *psol = sol.GetVect ();

	for (i=0;i<_n;i++) _sol[i] = psol[i];

};

// Author: Kharchenko S.A.
// CSVector: Solve unsymmetric system
//========================================================================================
void CSVector::Solver (ofstream &_fout, CSMatrixR &_mtra, // Solve unsymmetric system
						const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm) {

	const char *funcname = "Solver";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;
//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();

	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [n];
	if (!order) MemoryFail (funcname);
	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
//		_mtra.PartOrdMtr (nparts, blks, order);
		int ordtyp = 1;
		_mtra.OrderPrfMtr (ordtyp, order);
	} else if (_param.ordtype == 0) {
		for (i=0;i<n;i++) order[i] = i;
	} else if (_param.ordtype < 0) {
		int ordtyp = -1;
		_mtra.OrderPrfMtr (ordtyp, order);
	};

	for (i=0;i<n;i++) iord[order[i]] = i;
//
// Reorder the matrix in-place
//
	CSMatrixR mtrao;

	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (order);
		_mtra = mtrao;

	};

	CSMatrixR mtrdummy;

	mtrao = mtrdummy;
	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrOrd.ps",0,&i);
	};
//
// Compute triangular L and U parts of the initial matrix
//
	CSMatrixR mtral, mtrau;

	_mtra.CombineLU (mtral, mtrau);

	_mtra = mtrdummy;
//
// Compute Ilu2
//
	CSMatrixR mtrl, mtru;

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

	_mtra.AddLU (mtral, mtrau);

	mtral = mtrdummy;
	mtrau = mtrdummy;

	if (_param.msglev > 1) {
		if (collap != 0) mtrl.A2Ps (collap,"StrIlu.ps",0,&i);
	};
//
// Assign pointers
//
	CSMatrix *pa, *pl, *pu;

	pa = &_mtra;
	pl = &mtrl;
	pu = &mtru;
//
// Switch solution cases
//
	int nrhsloc = GetNrhs ();
//
// Cycle over the right hand sides
//
	CSVector sol(n), x(n), b(n), norm;

	for (int irhs=0;irhs<nrhsloc;irhs++) {

		if (_param.msglev > 1) {
			cout  << endl;
			cout  << "   Solve Irhs = " << irhs << endl;
			cout  << endl;
		};
		if (_param.msglev >= 1) {
			_fout << endl;
			_fout << "   Solve Irhs = " << irhs << endl;
			_fout << endl;
		};
//
// Init right hand side and initial guess
//
		for (i=0;i<n;i++) {
			b.vect[i] = _b.vect[irhs*n+i];
		};
		for (i=0;i<n;i++) {
			x.vect[i] = vect[irhs*n+i];
		};
		sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
		CSVector vdummy(n);

		vdummy.OrdVect('D', order, x);

		x = vdummy;

		vdummy.OrdVect('D', order, b);

		b = vdummy;

		vdummy.OrdVect('D', order, _norm);

		norm = vdummy;

		double aux=0.0e0;

		aux = b.ScProdSVect (b);

		if (aux > 0.0e0) sol.Lanczos (_fout, pa, pl, pu, x, b, _param, norm);
//
// Reorder solution back
//
		vdummy.OrdVect('I', order, sol);

		x = vdummy;
//
// Store the solution
//
		for (i=0;i<n;i++) {
			vect[irhs*n+i] = x.vect[i];
		};

	};
//
// Reorder the matrix back
//
	mtrl = mtrdummy;
	mtru = mtrdummy;

	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (iord);

		_mtra = mtrao;

	};

	mtrao = mtrdummy;
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_param.msglev > 1) {
		cout  << " Solver statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec." << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Solver statistics: " << endl;
		_fout << "     Time  = " << tottim << " sec." << endl;
	};

};

// Author: Kharchenko S.A.
// CSVector: Solve unsymmetric system
//========================================================================================
void CSVector::Solver (ofstream &_fout, CSMatrixRS &_mtra, // Solve unsymmetric system
						const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm) {

	const char *funcname = "Solver";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;
//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();
	int nsupa = _mtra.GetNsupr();

	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;
	int *sprndcini, *sprndcord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [nsupa];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupa];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupa+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupa+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
		_mtra.PartOrdMtr (nparts, blks, order);
	} else if (_param.ordtype == 0) {
		for (i=0;i<nsupa;i++) order[i] = i;
	} else if (_param.ordtype < 0) {
		_mtra.OrderBrdMtr (nparts, order);
	};

	for (i=0;i<nsupa;i++) iord[order[i]] = i;

	for (i=0;i<=nsupa;i++) sprndcini[i] = _mtra.sprndc[i];
//
// Reorder the matrix in-place
//
	CSMatrixRS mtrao;

	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (order);
		_mtra = mtrao;

	};

	CSMatrixRS mtrdummy;

	mtrao = mtrdummy;
	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrOrd.ps",0,&i);
	};

	for (i=0;i<=nsupa;i++) sprndcord[i] = _mtra.sprndc[i];
//
// Compute triangular L and U parts of the initial matrix
//
	CSMatrixRS mtral, mtrau;

	_mtra.CombineLU (mtral, mtrau);
//
// Compute Ilu2
//
	CSMatrixRS mtrl, mtru;

	if (_param.memory >= 0.0e0) {
		mtrdummy.Ilu2 (_fout, _param,
						mtral, mtrau,
						mtrl, mtru);
	} else {
		mtrdummy.Ilu2DynamicByBlocks (_fout, _param,
						mtral, mtrau,
						mtrl, mtru);
	};

	mtral = mtrdummy;
	mtrau = mtrdummy;

	if (_param.msglev > 1) {
		if (collap != 0) mtrl.A2Ps (collap,"StrIlu.ps",0,&i);
	};
//
// Assign pointers
//
	CSMatrix *pa, *pl, *pu;

	pa = &_mtra;
	pl = &mtrl;
	pu = &mtru;
//
// Switch solution cases
//
	int nrhsloc = GetNrhs ();

	if (_param.ittype == 1) {
//
// Cycle over the right hand sides
//
		CSVector sol(n), x(n), b(n), norm(n);

		for (int irhs=0;irhs<nrhsloc;irhs++) {

			if (_param.msglev > 1) {
				cout  << endl;
				cout  << "   Solve Irhs = " << irhs << endl;
				cout  << endl;
			};
			if (_param.msglev >= 1) {
				_fout << endl;
				_fout << "   Solve Irhs = " << irhs << endl;
				_fout << endl;
			};
//
// Init right hand side and initial guess
//
			for (i=0;i<n;i++) {
				b.vect[i] = _b.vect[irhs*n+i];
			};
			for (i=0;i<n;i++) {
				x.vect[i] = vect[irhs*n+i];
			};
			for (i=0;i<n;i++) {
				norm.vect[i] = _norm.vect[i];
			};
			sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
			x.OrdVect('D', nsupa, sprndcini, sprndcord, order);
			b.OrdVect('D', nsupa, sprndcini, sprndcord, order);
			norm.OrdVect('D', nsupa, sprndcini, sprndcord, order);

			double aux=0.0e0;

			aux = b.ScProdSVect (b);

			if (aux > 0.0e0) sol.Lanczos (_fout, pa, pl, pu, x, b, _param, norm);
//
// Reorder solution back
//
			sol.OrdVect('I', nsupa, sprndcini, sprndcord, order);

			x = sol;
//
// Store the solution
//
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[i];
			};

		};

	} else {
//
// Determine the number of nonzero right hand sides and store the list
//
		int nlistrhs = 0;

		int *listrhs;

		listrhs = new int [nrhsloc];
		if (!listrhs) MemoryFail (funcname);

		int irhs;
		for (irhs=0;irhs<nrhsloc;irhs++) {

			double aux = 0.0e0;

			for (i=0;i<n;i++) {
				double auxl = _b.vect[irhs*n+i];
				aux += auxl*auxl;
			};

			if (aux != 0.0e0) {
				listrhs[nlistrhs] = irhs;
				nlistrhs++;
			};

		};
//
// Solve by block Lanczos
//
		CSVector sol(n,1,nlistrhs), x(n,1,nlistrhs), b(n,1,nlistrhs);
//
// Init right hand side and initial guess
//
		int ilist;
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<n;i++) {
				b.vect[ilist*n+i] = _b.vect[irhs*n+i];
			};
			for (i=0;i<n;i++) {
				x.vect[ilist*n+i] = vect[irhs*n+i];
			};
		};

		sol.SetSVect(0.0e0);
//
// Reorder right hand side and initial guess
//
		x.OrdVect('D', nsupa, sprndcini, sprndcord, order);
		b.OrdVect('D', nsupa, sprndcini, sprndcord, order);

		if (_param.ittype == 2) sol.BlockLanczosCentral (_fout, &_mtra, &mtrl, &mtru, x, b, _param);
		if (_param.ittype == 3) sol.BlockLanczosRight   (_fout, &_mtra, &mtrl, &mtru, x, b, _param);
		if (_param.ittype == 4) sol.BlockLanczosLeft    (_fout, &_mtra, &mtrl, &mtru, x, b, _param);
//
// Reorder solution back
//
		sol.OrdVect('I', nsupa, sprndcini, sprndcord, order);

		x = sol;
//
// Store the solution
//
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[ilist*n+i];
			};
		};

		delete [] listrhs;

	};
//
// Reorder the matrix back
//
	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (iord);

		_mtra = mtrao;

	};

	mtrao = mtrdummy;
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_param.msglev > 1) {
		cout  << " Solver statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec." << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Solver statistics: " << endl;
		_fout << "     Time  = " << tottim << " sec." << endl;
	};

};

// Author: Kharchenko S.A.
// CSVector: Solve unsymmetric system in parallel mode
//========================================================================================
void CSVector::Solver (ofstream &_fout, const CTree &_tree, // Solve unsymmetric system in parallel mode
						CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau, 
						CGSMatrixRS &_gmtralfct, CGSMatrixRS &_gmtraufct,
						bool _second, CGSMatrixRS &_gmtra2,
						CGSMatrixRS &_gmtrl, CGSMatrixRS &_gmtru,
						const CSVector &_b, const CSVector &_x, CSVector &_norm, 
						CSlvParam &_param) {

//	const char *funcname = "Solver";

	double tottim;
	clock_t time0, time1;
//
// Free L and U on entry if necessary
//
	CGSMatrixRS gmtrdummy;

	CSMatrixRS mtrdummy;

	int iblk;
	for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
		_gmtrl.mtrarr[iblk] = mtrdummy;
	};
	for (iblk=0;iblk<_gmtru.nblksr;iblk++) {
		_gmtru.mtrarr[iblk] = mtrdummy;
	};
	_gmtrl = gmtrdummy;
	_gmtru = gmtrdummy;
//
// Init statistics data
//
	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	gmtrdummy.Ilu2Schur (_fout, _tree, _param,
							_gmtralfct, _gmtraufct,
							_gmtrl, _gmtru);

//	_fout << " L in solver " << _gmtrl << endl;
//	_fout << " U in solver " << _gmtru << endl;
//
// Solve by Lanczos
//
	int nrhsloc = GetNrhs ();
//	int npartsloc = GetNparts ();

	SetSVect(0.0e0);

//	_fout << " X = " << _x << endl;
//	_fout << " B = " << _b << endl;
//	_fout << " Norm = " << _norm << endl;

	CMvmR mvm (_tree, _gmtral, nrhsloc);

	if (_param.ittype == 1) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Cg (_fout, 
					_tree,
					&objmva, CParMvmA::MvmA,
					&objslvlu, CParSolveLU::SolveLU,
					_x, _b, _param);

		} else {

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Cg (_fout, 
					_tree,
					&objmva, CParMvmA::MvmA,
					&objslvlu, CParSolveLU::SolveLU,
					_x, _b, _param);

		};
	} else if (_param.ittype == 2) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Lanczos (_fout, _tree, 
						&objmva, CParMvmA::MvmA, CParMvmA::MvmAt,
						&objslvlu, CParSolveLU::SolveL, CParSolveLU::SolveLt,
						&objslvlu, CParSolveLU::SolveU, CParSolveLU::SolveUt,
						_x, _b, _norm, _param);

		} else {

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Lanczos (_fout, _tree, 
						&objmva, CParMvmA::MvmA, CParMvmA::MvmAt,
						&objslvlu, CParSolveLU::SolveL, CParSolveLU::SolveLt,
						&objslvlu, CParSolveLU::SolveU, CParSolveLU::SolveUt,
						_x, _b, _norm, _param);

		};
	} else {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmva, CParMvmA::MvmA,
							&objslvlu, CParSolveLU::SolveLU,
							_x, _b);

		} else {

			CParMvmA objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmva, CParMvmA::MvmA,
							&objslvlu, CParSolveLU::SolveLU,
							_x, _b);

		};
	};

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		if (_param.msglev > 1) {
			cout  << " Solver statistics: " << endl;
			cout  << "     Time  = " << tottim << " sec." << endl;
		};
		if (_param.msglev > 0) {
			_fout << " Solver statistics: " << endl;
			_fout << "     Time  = " << tottim << " sec." << endl;
		};

	};

};

// Author: Kharchenko S.A.
// CSVector: Solve unsymmetric system in parallel mode
//========================================================================================
void CSVector::Solver2Ind (ofstream &_fout, const CTree &_tree, // Solve unsymmetric system in parallel mode
									CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
									CGSMatrixR &_gmtralfct, CGSMatrixR &_gmtraufct,
									bool _second, CGSMatrixR &_gmtra2,
									CGSMatrixR &_gmtrl, CGSMatrixR &_gmtru,
									const CSVector &_b, const CSVector &_x, CSVector &_norm, 
									CSlvParam &_param) {

//	const char *funcname = "Solver2Ind";

	double tottim;
	clock_t time0, time1;
//
// Free L and U on entry if necessary
//
	CGSMatrixR gmtrdummy;

	CSMatrixR mtrdummy;

	int iblk;
	for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
		_gmtrl.mtrarr[iblk] = mtrdummy;
	};
	for (iblk=0;iblk<_gmtru.nblksr;iblk++) {
		_gmtru.mtrarr[iblk] = mtrdummy;
	};
	_gmtrl = gmtrdummy;
	_gmtru = gmtrdummy;
//
// Init statistics data
//
	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	StartTimer ();

//	gmtrdummy.Ilu2Schur2Index (_fout, _tree, _param,
//							_gmtralfct, _gmtraufct,
//							_gmtrl, _gmtru);
#ifdef __FV_2004__
   flowvision::STARTTIMER();
#endif
	gmtrdummy.Ilu2SchurDynamic2Index (_fout, _tree, _param,
							_gmtralfct, _gmtraufct,
							_gmtrl, _gmtru);
#ifdef __FV_2004__
   flowvision::STOPTIMER("          Fct time");
#endif

	double dtime1, dtime2;

	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

//	_fout << " L in solver " << _gmtrl << endl;
//	_fout << " U in solver " << _gmtru << endl;
//
// Solve by Lanczos
//
	int nrhsloc = GetNrhs ();
//	int npartsloc = GetNparts ();

	SetSVect(0.0e0);

//	_fout << " X = " << _x << endl;
//	_fout << " B = " << _b << endl;
//	_fout << " Norm = " << _norm << endl;

	CMvmR mvm (_tree, _gmtral, nrhsloc);

#ifdef __FV_2004__
   flowvision::STARTTIMER();
#endif
	StartTimer ();

	if (_param.ittype == 1) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Cg (_fout, 
					_tree,
					&objmva, CParMvmA2Index::MvmA,
					&objslvlu, CParSolveLU2Index::SolveLU,
					_x, _b, _param);

		} else {

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Cg (_fout, 
					_tree,
					&objmva, CParMvmA2Index::MvmA,
					&objslvlu, CParSolveLU2Index::SolveLU,
					_x, _b, _param);

		};
	} else if (_param.ittype == 2) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Lanczos (_fout, _tree, 
						&objmva, CParMvmA2Index::MvmA, CParMvmA2Index::MvmAt,
						&objslvlu, CParSolveLU2Index::SolveL, CParSolveLU2Index::SolveLt,
						&objslvlu, CParSolveLU2Index::SolveU, CParSolveLU2Index::SolveUt,
						_x, _b, _norm, _param);

		} else {

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			Lanczos (_fout, _tree, 
						&objmva, CParMvmA2Index::MvmA, CParMvmA2Index::MvmAt,
						&objslvlu, CParSolveLU2Index::SolveL, CParSolveLU2Index::SolveLt,
						&objslvlu, CParSolveLU2Index::SolveU, CParSolveLU2Index::SolveUt,
						_x, _b, _norm, _param);

		};
	} else {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm2, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmva, CParMvmA2Index::MvmA,
							&objslvlu, CParSolveLU2Index::SolveLU,
							_x, _b);

		} else {

			CParMvmA2Index objmva ((CTree *)&_tree, &mvm, &_gmtral, &_gmtrau,
									_second, &mvm, &_gmtra2);

			CParSolveLU2Index objslvlu ((CTree *)&_tree, &mvm, &_gmtrl, &_gmtru);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmva, CParMvmA2Index::MvmA,
							&objslvlu, CParSolveLU2Index::SolveLU,
							_x, _b);

		};
	};

#ifdef __FV_2004__
   flowvision::STOPTIMER("          Iter time");
#endif
	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		if (_param.msglev > 1) {
			cout  << " Solver statistics: " << endl;
			cout  << "     Time  = " << tottim << " sec." << endl;
		};
		if (_param.msglev > 0) {
			_fout << " Solver statistics: " << endl;
			_fout << "     Time  = " << tottim << " sec." << endl;
		};

	};

};

// Author: Kharchenko S.A.
// CSVector: Solve symmetric system in parallel mode
//========================================================================================
void CSVector::SolverSymm2IndSchur (ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
									CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
									CGSMatrixR &_gmtralfct,
									bool _second, CGSMatrixR &_gmtra2,
									CGSMatrixR &_gmtrl,
									const CSVector &_b, const CSVector &_x, CSVector &_norm, 
									CSlvParam &_param) {

//	const char *funcname = "SolverSymm2IndSchur";

	double tottim;
	clock_t time0, time1;
//
// Free L and U on entry if necessary
//
	CGSMatrixR gmtrdummy;

	CSMatrixR mtrdummy;

	int iblk;
	for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
		_gmtrl.mtrarr[iblk] = mtrdummy;
	};
	_gmtrl = gmtrdummy;
//
// Init statistics data
//
	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	StartTimer ();
#ifdef __FV_2004__
	flowvision::STARTTIMER();
#endif

	int nblkstree = _tree.GetNblkstree ();
	int *pblkstree = _tree.GetBlkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();
	CSMatrix * pablkstr = _tree.GetAblktree ();

	CFctSymmSchurR fctschurr (&_fout, 
									nblkstree, pblkstree, pblk2cputree,
									&_tree, pablkstr, 
									&_gmtralfct, 
									&_gmtrl,
									&_param);

	fctschurr.AbstractBlockIlu2 (_tree, *pablkstr, nblkstree, pblk2cputree);

//	_fout << " L = " << _gmtrl << endl;
//	_fout << " U = " << _gmtru << endl;

#ifdef __FV_2004__
	flowvision::STOPTIMER("          Fct time");
#endif

	double dtime1, dtime2;

	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

//
// Solve by Lanczos
//
	int nrhsloc = GetNrhs ();

	SetSVect(0.0e0);

	CMvmSchurR mvmschur (&_fout,
									nblkstree, pblkstree, pblk2cputree,
									&_tree, pablkstr, 
									&_gmtral, &_gmtrau, 
									&_gmtrl, &_gmtrl,
									CSMatrixR::AbsMvmBlockColumnAl2Index , CSMatrixR::AbsMvmBlockRowAu2Index,
									CSMatrixR::AbsSolveBlockColumnL2Index, CSMatrixR::AbsSolveBlockRowU2Index,
									&_param);

#ifdef __FV_2004__
	flowvision::STARTTIMER();
#endif
	StartTimer ();

	if (_second) {

		CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

		CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, &mvm2, &_gmtra2);

		Cg (_fout, 
				_tree,
				&objmvslv, CParSchurMvmSlv2Index::MvmA,
				&objmvslv, CParSchurMvmSlv2Index::SolveLU,
				_x, _b, _param);

	} else {

		CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, 0, &_gmtra2);

		Cg (_fout, 
				_tree,
				&objmvslv, CParSchurMvmSlv2Index::MvmA,
				&objmvslv, CParSchurMvmSlv2Index::SolveLU,
				_x, _b, _param);

	};

#ifdef __FV_2004__
	flowvision::STOPTIMER("          Iter time");
#endif
	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		if (_param.msglev > 1) {
			cout  << " Solver statistics: " << endl;
			cout  << "     Time  = " << tottim << " sec." << endl;
		};
		if (_param.msglev > 0) {
			_fout << " Solver statistics: " << endl;
			_fout << "     Time  = " << tottim << " sec." << endl;
		};

	};

};

// Author: Kharchenko S.A.
// CSVector: Solve unsymmetric system in parallel mode
//========================================================================================
void CSVector::Solver2IndSchur (ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
									CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
									CGSMatrixR &_gmtralfct, CGSMatrixR &_gmtraufct,
									bool _second, CGSMatrixR &_gmtra2,
									CGSMatrixR &_gmtrl, CGSMatrixR &_gmtru,
									const CSVector &_b, const CSVector &_x, CSVector &_norm, 
									CSlvParam &_param) {

//	const char *funcname = "Solver2IndSchur";

	double tottim;
	clock_t time0, time1;
//
// Free L and U on entry if necessary
//
	CGSMatrixR gmtrdummy;

	CSMatrixR mtrdummy;

	int iblk;
	for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
		_gmtrl.mtrarr[iblk] = mtrdummy;
	};
	for (iblk=0;iblk<_gmtru.nblksr;iblk++) {
		_gmtru.mtrarr[iblk] = mtrdummy;
	};
	_gmtrl = gmtrdummy;
	_gmtru = gmtrdummy;
//
// Init statistics data
//
	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	StartTimer ();
#ifdef __FV_2004__
	flowvision::STARTTIMER();
#endif

	int nblkstree = _tree.GetNblkstree ();
	int *pblkstree = _tree.GetBlkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();
	CSMatrix * pablkstr = _tree.GetAblktree ();

	CFctSchurR fctschurr (&_fout, 
									nblkstree, pblkstree, pblk2cputree,
									&_tree, pablkstr, 
									&_gmtralfct, &_gmtraufct, 
									&_gmtrl, &_gmtru,
									&_param);

	fctschurr.AbstractBlockIlu2 (_tree, *pablkstr, nblkstree, pblk2cputree);

//	_fout << " L = " << _gmtrl << endl;
//	_fout << " U = " << _gmtru << endl;

#ifdef __FV_2004__
	flowvision::STOPTIMER("          Fct time");
#endif

	double dtime1, dtime2;

	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Incomplete factorization time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

//
// Solve by Lanczos
//
	int nrhsloc = GetNrhs ();

	SetSVect(0.0e0);

	CMvmSchurR mvmschur (&_fout,
									nblkstree, pblkstree, pblk2cputree,
									&_tree, pablkstr, 
									&_gmtral, &_gmtrau, 
									&_gmtrl, &_gmtru,
									CSMatrixR::AbsMvmBlockColumnAl2Index , CSMatrixR::AbsMvmBlockRowAu2Index,
									CSMatrixR::AbsSolveBlockColumnL2Index, CSMatrixR::AbsSolveBlockRowU2Index,
									&_param);

#ifdef __FV_2004__
	flowvision::STARTTIMER();
#endif
	StartTimer ();

	if (_param.ittype == 1) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, &mvm2, &_gmtra2);

			Cg (_fout, 
					_tree,
					&objmvslv, CParSchurMvmSlv2Index::MvmA,
					&objmvslv, CParSchurMvmSlv2Index::SolveLU,
					_x, _b, _param);

		} else {

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, 0, &_gmtra2);

			Cg (_fout, 
					_tree,
					&objmvslv, CParSchurMvmSlv2Index::MvmA,
					&objmvslv, CParSchurMvmSlv2Index::SolveLU,
					_x, _b, _param);

		};
	} else if (_param.ittype == 2) {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, &mvm2, &_gmtra2);

			Lanczos (_fout, _tree, 
						&objmvslv, CParSchurMvmSlv2Index::MvmA, CParSchurMvmSlv2Index::MvmAt,
						&objmvslv, CParSchurMvmSlv2Index::SolveL, CParSchurMvmSlv2Index::SolveLt,
						&objmvslv, CParSchurMvmSlv2Index::SolveU, CParSchurMvmSlv2Index::SolveUt,
						_x, _b, _norm, _param);

		} else {

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, 0, &_gmtra2);

			Lanczos (_fout, _tree, 
						&objmvslv, CParSchurMvmSlv2Index::MvmA, CParSchurMvmSlv2Index::MvmAt,
						&objmvslv, CParSchurMvmSlv2Index::SolveL, CParSchurMvmSlv2Index::SolveLt,
						&objmvslv, CParSchurMvmSlv2Index::SolveU, CParSchurMvmSlv2Index::SolveUt,
						_x, _b, _norm, _param);

		};
	} else {
		if (_second) {

			CMvmR mvm2 (true, _tree, _gmtra2, nrhsloc);

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, &mvm2, &_gmtra2);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmvslv, CParSchurMvmSlv2Index::MvmA,
							&objmvslv, CParSchurMvmSlv2Index::SolveLU,
							_x, _b);

		} else {

			CParSchurMvmSlv2Index objmvslv ((CTree *)&_tree, &mvmschur, _second, 0, &_gmtra2);

			BlockGMRES (_fout, _param, 
							_tree,
							&objmvslv, CParSchurMvmSlv2Index::MvmA,
							&objmvslv, CParSchurMvmSlv2Index::SolveLU,
							_x, _b);

		};
	};

#ifdef __FV_2004__
	flowvision::STOPTIMER("          Iter time");
#endif
	StopTimer (dtime1, dtime2);

	if (_tree.myid == 0 && _param.msglev > 0) {
		cout  << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
		_fout << "    Iterative Solution time: CPU time = " << dtime1 << " Wall time = " << dtime2 << endl;
	};

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		if (_param.msglev > 1) {
			cout  << " Solver statistics: " << endl;
			cout  << "     Time  = " << tottim << " sec." << endl;
		};
		if (_param.msglev > 0) {
			_fout << " Solver statistics: " << endl;
			_fout << "     Time  = " << tottim << " sec." << endl;
		};

	};

};

// Author: Kharchenko S.A.
// CSVectorC: Solve unsymmetric system
//========================================================================================
void CSVectorC::Solver (ofstream &_fout, CSMatrixCS &_mtra, // Solve unsymmetric system
						int _job, int *&_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
						const CSVectorC &_b, 
						const CSlvParam &_param) {

	const char *funcname = "Solver";

	dcmplx czero (0.0e0,0.0e0);

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;

//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();
	int nsupa = _mtra.GetNsupr();

	if (_job == 0) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
	};
//
// Free L and U on entry if necessary
//
	CSMatrixCS mtrdummy;

	if (_job == 0) {
		_mtrl = mtrdummy;
		_mtru = mtrdummy;
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;
	int *sprndcini, *sprndcord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [nsupa];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupa];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupa+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupa+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_job == 0) {
		if (_param.ordtype > 0) {
			_mtra.PartOrdMtr (nparts, blks, order);
		} else if (_param.ordtype == 0) {
			for (i=0;i<nsupa;i++) order[i] = i;
		} else if (_param.ordtype < 0) {
			_mtra.OrderBrdMtr (nparts, order);
		};

		delete [] _order;
		_order = new int [nsupa];
		if (!_order) MemoryFail (funcname);

		for (i=0;i<nsupa;i++) _order[i] = order[i];

	} else {
		for (i=0;i<nsupa;i++) order[i] = _order[i];
	};

	for (i=0;i<nsupa;i++) iord[order[i]] = i;

	for (i=0;i<=nsupa;i++) sprndcini[i] = _mtra.sprndc[i];
//
// Reorder the matrix in-place
//
	CSMatrixCS mtrao;

	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (order);
		_mtra = mtrao;

	};

	mtrao = mtrdummy;

	if (_job == 0) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrOrd.ps",0,&i);
	};

	for (i=0;i<=nsupa;i++) sprndcord[i] = _mtra.sprndc[i];
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	if (_job == 0) {

		CSMatrixCS mtral, mtrau;

		_mtra.CombineLU (mtral, mtrau);
//
// Compute Ilu2
//
		mtrdummy.Ilu2 (_fout, _param,
						mtral, mtrau,
						_mtrl, _mtru);

		if (collap != 0) _mtrl.A2Ps (collap,"StrIlu.ps",0,&i);

		mtral = mtrdummy;
		mtrau = mtrdummy;

	};
//
// Assign pointers
//
	CSMatrix *pa, *pl, *pu;

	pa = &_mtra;
	pl = &_mtrl;
	pu = &_mtru;
//
// Determine the number of nonzero right hand sides and store the list
//
	int nrhsloc = GetNrhs ();

	int nlistrhs = 0;

	int *listrhs;

	listrhs = new int [nrhsloc];
	if (!listrhs) MemoryFail (funcname);

	int irhs;
	for (irhs=0;irhs<nrhsloc;irhs++) {

		double aux = 0.0e0;

		for (i=0;i<n;i++) {
			double auxl = _b.vect[irhs*n+i].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*n+i].imag ();
			aux += auxl*auxl;
		};

		if (aux != 0.0e0) {
			listrhs[nlistrhs] = irhs;
			nlistrhs++;
		};

	};
//
// Solve by block Lanczos
//
	CSVectorC sol(n,1,nlistrhs), x(n,1,nlistrhs), b(n,1,nlistrhs);
//
// Init right hand side and initial guess
//
	int ilist;
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (i=0;i<n;i++) {
			b.vect[ilist*n+i] = _b.vect[irhs*n+i];
		};
		for (i=0;i<n;i++) {
			x.vect[ilist*n+i] = vect[irhs*n+i];
		};
	};

	sol.SetSVect(czero);
//
// Reorder right hand side and initial guess
//
	x.OrdVect('D', nsupa, sprndcini, sprndcord, order);
	b.OrdVect('D', nsupa, sprndcini, sprndcord, order);

	sol.BlockGMRES (_fout, &_mtra, &_mtrl, &_mtru, x, b, _param);
//
// Reorder solution back
//
	sol.OrdVect('I', nsupa, sprndcini, sprndcord, order);

	x = sol;
//
// Store the solution
//
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (i=0;i<n;i++) {
			vect[irhs*n+i] = x.vect[ilist*n+i];
		};
	};

	delete [] listrhs;
//
// Reorder the matrix back
//
	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtr (iord);

		_mtra = mtrao;

	};

	mtrao = mtrdummy;
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	cout  << " Solver statistics: " << endl;
	_fout << " Solver statistics: " << endl;

	cout  << "     Time  = " << tottim << " sec." << endl;
	_fout << "     Time  = " << tottim << " sec." << endl;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve unsymmetric system in the out-of-core mode
//========================================================================================
void CSVectorC::Solver (ofstream &_fout, // Solve unsymmetric system in the out-of-core mode
						int _nblks, int *_blks, FILE **_mtrafiles,
						const CSMatrixCS &_mtra,
						char *_path,
						int _job, int *&_order, 
						FILE **_mtrlfiles, FILE **_mtrufiles,
						CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
						const CSVectorC &_b, 
						const CSlvParam &_param) {

	const char *funcname = "Solver";

	dcmplx czero (0.0e0,0.0e0);

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;
//
// Flush files with the matrix data
//
	FlushDIOFiles (_nblks, _mtrafiles);
//
// Create arrays of files
//
	FILE **mtraofiles, **mtralfiles, **mtraufiles, **mtrqfiles;

	mtraofiles = new FILE * [_nblks];
	if (!mtraofiles) MemoryFail (funcname);
	mtralfiles = new FILE * [_nblks];
	if (!mtralfiles) MemoryFail (funcname);
	mtraufiles = new FILE * [_nblks];
	if (!mtraufiles) MemoryFail (funcname);
	mtrqfiles = new FILE * [_nblks];
	if (!mtrqfiles) MemoryFail (funcname);

	OpenDIOFiles (_nblks, mtraofiles, _path, "AoMatr_", ".bin");
	OpenDIOFiles (_nblks, mtralfiles, _path, "AlMatr_", ".bin");
	OpenDIOFiles (_nblks, mtraufiles, _path, "AuMatr_", ".bin");
	OpenDIOFiles (_nblks, mtrqfiles,  _path, "Gmres_", ".bin");
//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();
	int nsupa = _mtra.GetNsupr();

	if (_job == 0) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
	};
//
// Free L and U on entry if necessary
//
	CSMatrixCS mtrdummy;

	if (_job == 0) {
		_mtrl = mtrdummy;
		_mtru = mtrdummy;
		OpenDIOFiles (_nblks, _mtrlfiles, _path, "LMatr_", ".bin");
		OpenDIOFiles (_nblks, _mtrufiles, _path, "UMatr_", ".bin");
	} else {
		FlushDIOFiles (_nblks, _mtrlfiles);
		FlushDIOFiles (_nblks, _mtrufiles);
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blksord, *order, *iord;
	int *sprndcini, *sprndcord;

	blksord = new int [nparts+2];
	if (!blksord) MemoryFail (funcname);
	order = new int [nsupa];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupa];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupa+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupa+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_job == 0) {
		if (_param.ordtype > 0) {
			_mtra.PartOrdMtr (nparts, blksord, order);
		} else if (_param.ordtype == 0) {
			for (i=0;i<nsupa;i++) order[i] = i;
		} else if (_param.ordtype < 0) {
			_mtra.OrderBrdMtr (nparts, order);
		};

		delete [] _order;
		_order = new int [nsupa];
		if (!_order) MemoryFail (funcname);

		for (i=0;i<nsupa;i++) _order[i] = order[i];

	} else {
		for (i=0;i<nsupa;i++) order[i] = _order[i];
	};

	for (i=0;i<nsupa;i++) iord[order[i]] = i;

	for (i=0;i<=nsupa;i++) sprndcini[i] = _mtra.sprndc[i];
//
// Reorder the matrix in-place
//
	CSMatrixCS mtrao;

	mtrao = _mtra.OrdMtr (_nblks, _blks, _mtrafiles, 
							order,
							mtraofiles);

	if (_job == 0) {
		if (collap != 0) mtrao.A2Ps (collap,"AStrOrd.ps",0,&i);
	};

	for (i=0;i<=nsupa;i++) sprndcord[i] = mtrao.sprndc[i];
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	if (_job == 0) {

		CSMatrixCS mtral, mtrau;

		mtrao.CombineLU (_nblks, _blks, mtraofiles, 
							mtralfiles, mtraufiles, 
							mtral, mtrau);
//
// Compute Ilu2
//
		mtrdummy.Ilu2 (_fout, _param,
						_nblks, _blks,
						mtralfiles, mtraufiles,
						mtral, mtrau,
						_mtrlfiles, _mtrufiles,
						_mtrl, _mtru);

		if (collap != 0) _mtrl.A2Ps (collap,"StrIlu.ps",0,&i);

		mtral = mtrdummy;
		mtrau = mtrdummy;

	};
//
// Assign pointers
//
	CSMatrix *pa, *pl, *pu;

	pa = &mtrao;
	pl = &_mtrl;
	pu = &_mtru;
//
// Determine the number of nonzero right hand sides and store the list
//
	int nrhsloc = GetNrhs ();

	int nlistrhs = 0;

	int *listrhs;

	listrhs = new int [nrhsloc];
	if (!listrhs) MemoryFail (funcname);

	int irhs;
	for (irhs=0;irhs<nrhsloc;irhs++) {

		double aux = 0.0e0;

		for (i=0;i<n;i++) {
			double auxl = _b.vect[irhs*n+i].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*n+i].imag ();
			aux += auxl*auxl;
		};

		if (aux != 0.0e0) {
			listrhs[nlistrhs] = irhs;
			nlistrhs++;
		};

	};
//
// Solve by block Lanczos
//
	CSVectorC sol(n,1,nlistrhs), x(n,1,nlistrhs), b(n,1,nlistrhs);
//
// Init right hand side and initial guess
//
	int ilist;
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (i=0;i<n;i++) {
			b.vect[ilist*n+i] = _b.vect[irhs*n+i];
		};
		for (i=0;i<n;i++) {
			x.vect[ilist*n+i] = vect[irhs*n+i];
		};
	};

	sol.SetSVect(czero);
//
// Reorder right hand side and initial guess
//
	x.OrdVect('D', nsupa, sprndcini, sprndcord, order);
	b.OrdVect('D', nsupa, sprndcini, sprndcord, order);

	sol.BlockGMRES (_fout,
					_nblks, _blks, mtraofiles, &mtrao,
					_mtrlfiles, _mtrufiles,
					&_mtrl, &_mtru, 
					x, b, _param, mtrqfiles);
//
// Reorder solution back
//
	sol.OrdVect('I', nsupa, sprndcini, sprndcord, order);

	x = sol;
//
// Store the solution
//
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (i=0;i<n;i++) {
			vect[irhs*n+i] = x.vect[ilist*n+i];
		};
	};

	delete [] listrhs;
//
// Close work files on disks
//
	CloseDIOFiles (_nblks, mtraofiles);
	CloseDIOFiles (_nblks, mtralfiles);
	CloseDIOFiles (_nblks, mtraufiles);
	CloseDIOFiles (_nblks, mtrqfiles);
//
// Free work arrays
//
	delete [] mtraofiles;
	delete [] mtralfiles;
	delete [] mtraufiles;
	delete [] mtrqfiles;
	delete [] blksord;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	cout  << " Solver statistics: " << endl;
	_fout << " Solver statistics: " << endl;

	cout  << "     Time  = " << tottim << " sec." << endl;
	_fout << "     Time  = " << tottim << " sec." << endl;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve unsymmetric system in parallel mode
//========================================================================================
void CSVectorC::Solver (ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
						CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau,
						int _job, CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
						const CSVectorC &_b, CSVectorC &_ycorr, 
						CSlvParam &_param) {

	const char *funcname = "Solver";

	dcmplx czero (0.0e0,0.0e0);

	double tottim;
	clock_t time0, time1;
//
// Get the sizes of the matrix
//
	int n = _b.nv;
//	int nsupa = _gmtral.GetNsupr();
//
// Free L and U on entry if necessary
//
	CGSMatrixCS gmtrdummy;

	CSMatrixCS mtrdummy;

	if (_job == 0) {
		int iblk;
		for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
			_gmtrl.mtrarr[iblk] = mtrdummy;
		};
		for (iblk=0;iblk<_gmtru.nblksr;iblk++) {
			_gmtru.mtrarr[iblk] = mtrdummy;
		};
		_gmtrl = gmtrdummy;
		_gmtru = gmtrdummy;
	};
//
// Init statistics data
//

	CMPIComm pcomm = _tree.GetComm ();

	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	if (_job == 0) {
//
// Compute Ilu2
//
		if (false) {
			gmtrdummy.Ilu2Schur (_fout, _tree, _param,
									_gmtral, _gmtrau,
									_gmtrl, _gmtru);
		} else {
			int nblkstree = _tree.GetNblkstree ();
			int *pblkstree = _tree.GetBlkstree ();
			int *pblk2cputree = _tree.GetBlk2cputree ();
			CSMatrix * pablkstr = _tree.GetAblktree ();

			CFctSchurC fctschurc (&_fout, 
											nblkstree, pblkstree, pblk2cputree,
											&_tree, pablkstr, 
											&_gmtral, &_gmtrau, 
											&_gmtrl, &_gmtru,
											&_param);

			fctschurc.AbstractBlockIlu2 (_tree, *pablkstr, nblkstree, pblk2cputree);
		};

	};
//
// Determine the number of nonzero right hand sides and store the list
//
///*
	int nrhsloc = GetNrhs ();

	int *listrhs;
	int *list0;
	int *list1;

	listrhs = new int [nrhsloc];
	if (!listrhs) MemoryFail (funcname);
	list0 = new int [nrhsloc];
	if (!list0) MemoryFail (funcname);
	list1 = new int [nrhsloc];
	if (!list1) MemoryFail (funcname);

	int irhs;
	for (irhs=0;irhs<nrhsloc;irhs++) {

		double aux = 0.0e0;

		for (int i=0;i<n;i++) {
			double auxl = _b.vect[irhs*n+i].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*n+i].imag ();
			aux += auxl*auxl;
		};

		if (aux != 0.0e0) {
			list0[irhs] = 1;
		} else {
			list0[irhs] = 0;
		};

	};

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nrhsloc, list0, list1);
	} else {
		for (int i=0;i<nrhsloc;i++) list1[i] = list0[i];
	};

	int nlistrhs = 0;

	for (irhs=0;irhs<nrhsloc;irhs++) {

		if (list1[irhs] != 0) {
			listrhs[nlistrhs] = irhs;
			nlistrhs++;
		};

	};

//
// Solve by block GMRES
//
	int npartsloc = GetNparts ();

	CSVectorC sol(n,npartsloc,nlistrhs), x(n,npartsloc,nlistrhs), b(n,npartsloc,nlistrhs);

	int ilist;
	for (ilist=0;ilist<npartsloc;ilist++) {
		sol.inda[ilist] = inda[ilist];
		x.inda[ilist] = inda[ilist];
		b.inda[ilist] = inda[ilist];
	};
	for (ilist=0;ilist<=npartsloc;ilist++) {
		sol.blksv[ilist] = blksv[ilist];
		x.blksv[ilist] = blksv[ilist];
		b.blksv[ilist] = blksv[ilist];
	};
//
// Init right hand side and initial guess
//
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		int i;
		for (i=0;i<n;i++) {
			b.vect[ilist*n+i] = _b.vect[irhs*n+i];
		};
		for (i=0;i<n;i++) {
			x.vect[ilist*n+i] = vect[irhs*n+i];
		};
	};

	sol.SetSVect(czero);

	if (false) {

		CMvmC mvm (false, _tree, _gmtral, nlistrhs);

		sol.BlockGMRES (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, x, b, _param);

	} else {

		int nblkstree = _tree.GetNblkstree ();
		int *pblkstree = _tree.GetBlkstree ();
		int *pblk2cputree = _tree.GetBlk2cputree ();
		int *pbl2nd = _gmtral.GetBl2ndc ();
		CSMatrix * pablkstr = _tree.GetAblktree ();

		CMvmSchurC mvmschurc (&_fout,
										nblkstree, pblkstree, pbl2nd, pblk2cputree,
										&_tree, pablkstr, 
										((void *)&_gmtral), ((void *)&_gmtrau), 
										((void *)&_gmtrl), ((void *)&_gmtru),
										CParSchurMvmSlvCS::AbsMvmBlockColumnAl, CParSchurMvmSlvCS::AbsMvmBlockRowAu, 
										CParSchurMvmSlvCS::AbsSolveBlockColumnL, CParSchurMvmSlvCS::AbsSolveBlockRowU,
										&_param);

		CParSchurMvmSlvCS objmvmcs (&_tree,
												&_gmtral, &_gmtrau, &_gmtrl, &_gmtru, 
												&mvmschurc);

		mvmschurc.SetSchurmvm (&objmvmcs);

		sol.BlockGMRES (_fout, _param,
								_tree, 
								(void *)(&objmvmcs), CParSchurMvmSlvCS::MvmA,
								(void *)(&objmvmcs), CParSchurMvmSlvCS::SolveLU,
								x, b);

	};
//
// Store the solution
//
	SetSVect(czero);

	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (int i=0;i<n;i++) {
			vect[irhs*n+i] = sol.vect[ilist*n+i];
		};
	};

	delete [] listrhs;
	delete [] list0;
	delete [] list1;

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		cout  << " Solver statistics: " << endl;
		_fout << " Solver statistics: " << endl;

		cout  << "     Time  = " << tottim << " sec." << endl;
		_fout << "     Time  = " << tottim << " sec." << endl;

	};

};

// Author: Kharchenko S.A.
// CSVectorC: Solve unsymmetric system in parallel mode
//========================================================================================
int CSVectorC::Solver (ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
						CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau,
						int _job, CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
						const CSVectorC &_b, CCorrectorC &_corr, 
						const CSlvParam &_param) {

	const char *funcname = "Solver";

	dcmplx czero (0.0e0,0.0e0);

	double tottim;
	clock_t time0, time1;
//
// Get the sizes of the matrix
//
	int n = _b.nv;
//	int nsupa = _gmtral.GetNsupr();
//
// Free L and U on entry if necessary
//
	CGSMatrixCS gmtrdummy;

	CSMatrixCS mtrdummy;

	if (_job == 0) {
		int iblk;
		for (iblk=0;iblk<_gmtrl.nblksr;iblk++) {
			_gmtrl.mtrarr[iblk] = mtrdummy;
		};
		for (iblk=0;iblk<_gmtru.nblksr;iblk++) {
			_gmtru.mtrarr[iblk] = mtrdummy;
		};
		_gmtrl = gmtrdummy;
		_gmtru = gmtrdummy;
	};
//
// Init statistics data
//
	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time0 = clock ();
//
// Compute triangular L and U parts of the initial matrix if necessary
//
	if (_job == 0) {
//
// Compute Ilu2
//
//		gmtrdummy.Ilu2 (_fout, _tree, _param,
//						_gmtral, _gmtrau,
//						_gmtrl, _gmtru);
		gmtrdummy.Ilu2Schur (_fout, _tree, _param,
								_gmtral, _gmtrau,
								_gmtrl, _gmtru);

	};

//
// Perform explicit point scaling of the coefficient matrix
//
	if (_gmtrl.ivect<2) {
		_gmtral.PointScaling (_gmtrl.vector1);
		_gmtrau.PointScaling (_gmtrl.vector1);
	} else {
		_gmtral.BlockScaling (_gmtru.vector1, _gmtrl.vector1);
		_gmtrau.BlockScaling (_gmtrl.vector1, _gmtru.vector1);
	};
//
// Determine the number of nonzero right hand sides and store the list
//
	int nrhsloc = GetNrhs ();

	int *listrhs;
	int *list0;
	int *list1;

	listrhs = new int [nrhsloc];
	if (!listrhs) MemoryFail (funcname);
	list0 = new int [nrhsloc];
	if (!list0) MemoryFail (funcname);
	list1 = new int [nrhsloc];
	if (!list1) MemoryFail (funcname);

	int irhs;
	for (irhs=0;irhs<nrhsloc;irhs++) {

		double aux = 0.0e0;

		for (int i=0;i<n;i++) {
			double auxl = _b.vect[irhs*n+i].real ();
			aux += auxl*auxl;
			auxl = _b.vect[irhs*n+i].imag ();
			aux += auxl*auxl;
		};

		if (aux != 0.0e0) {
			list0[irhs] = 1;
		} else {
			list0[irhs] = 0;
		};

	};

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nrhsloc, list0, list1);
	} else {
		for (int i=0;i<nrhsloc;i++) list1[i] = list0[i];
	};

	int nlistrhs = 0;

	for (irhs=0;irhs<nrhsloc;irhs++) {

		if (list1[irhs] != 0) {
			listrhs[nlistrhs] = irhs;
			nlistrhs++;
		};

	};
//
// Solve by block GMRES
//
	int npartsloc = GetNparts ();

	CSVectorC sol(n,npartsloc,nlistrhs), x(n,npartsloc,nlistrhs), b(n,npartsloc,nlistrhs);

	int ilist;
	for (ilist=0;ilist<npartsloc;ilist++) {
		sol.inda[ilist] = inda[ilist];
		x.inda[ilist] = inda[ilist];
		b.inda[ilist] = inda[ilist];
	};
	for (ilist=0;ilist<=npartsloc;ilist++) {
		sol.blksv[ilist] = blksv[ilist];
		x.blksv[ilist] = blksv[ilist];
		b.blksv[ilist] = blksv[ilist];
	};
//
// Init right hand side and initial guess
//
	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		int i;
		for (i=0;i<n;i++) {
			b.vect[ilist*n+i] = _b.vect[irhs*n+i];
		};
		for (i=0;i<n;i++) {
			x.vect[ilist*n+i] = vect[irhs*n+i];
		};
	};

	sol.SetSVect(czero);

// Perform iterations

	if (_gmtrl.ivect==1) {
		x.ScaleVect ('I', _gmtrl.vector1);
		b.ScaleVect ('D', _gmtrl.vector1);
	} else if (_gmtrl.ivect==2) {
		x.BlockScaleVect ('T', _gmtral, _gmtru.vector2);
		b.BlockScaleVect ('N', _gmtral, _gmtrl.vector1);
	};

	int iconv = 0;

	if (_param.thrfltrhs <= 0.0e0) {

		if (nlistrhs > 0) {

			CMvmC mvm (false, _tree, _gmtral, nlistrhs);

			sol.BlockGMRES (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, x, b, _param);

// Perform explicit point scaling and solve

//			if (_param.ittype == 2) {

//				iconv = sol.BlockGMRESCorr (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, x, b, _corr, _param);

//			} else {

//				sol.BlockLanczosRight (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, x, b, _param);

//			};

		};

	} else {
//
// Create collapsed tree
//
		CTree treeclp;

		treeclp = _tree.CollapseTree ();

// Compute QR decomposition of the reduced rhs

		int nrhsnew = 0;

		CQrdC qrdloc (treeclp, 1, nlistrhs);

		CSVectorC qloc;

		int lwork = 10*nlistrhs;
		int kloc = nlistrhs;

		dcmplx *tauloc, *rloc, *uloc, *vloc, *work;
		double *svloc, *rwork;

		tauloc = new dcmplx [nlistrhs];
		if (!tauloc) MemoryFail (funcname);
		rloc = new dcmplx [nlistrhs*nlistrhs];
		if (!rloc) MemoryFail (funcname);
		uloc = new dcmplx [nlistrhs*nlistrhs];
		if (!uloc) MemoryFail (funcname);
		vloc = new dcmplx [nlistrhs*nlistrhs];
		if (!vloc) MemoryFail (funcname);
		work = new dcmplx [lwork];
		if (!work) MemoryFail (funcname);
		svloc = new double [nlistrhs];
		if (!svloc) MemoryFail (funcname);
		rwork = new double [lwork];
		if (!rwork) MemoryFail (funcname);

		qloc = b;

		if (nlistrhs > 0) {
			qrdloc.UpdateQrd (treeclp, n, 0, nlistrhs, qloc.vect, n, tauloc);
		};

// Process R factor and compute the rotation matrix in parallel mode

		int i, j, irhsnew;
		double aux, aux1, aux2;
		dcmplx caux;

		if (_tree.myid == _tree.rootcpu) {

			if (nlistrhs > 0) {
				qrdloc.GetRPartQrd (treeclp, 0, 0, nlistrhs,
									rloc, nlistrhs);
			};

			for (i=0;i<kloc;i++) {
				aux = 0.0e0;
				for (j=0;j<kloc;j++) {
					aux1 = rloc[i*kloc+j].x;
					aux2 = rloc[i*kloc+j].y;
					aux += aux1 * aux1 + aux2 * aux2;
				};
				if (aux <= 0.0e0) aux = 1.0e0;
				aux = 1.0e0 / sqrt(aux);
				for (j=0;j<kloc;j++) {
					caux = rloc[i*kloc+j] * aux;
					rloc[i*kloc+j] = caux;
				};
			};

			int info=0;

			if (kloc > 0) {
				zgesvd_ ("A", "A", &kloc, &kloc, 
							rloc, &kloc, svloc, uloc, &kloc, vloc, &kloc, 
							work, &lwork, rwork, &info);
			};

			if (info != 0) {
				cout << " Error in the Lapack routine ZGESVD" << endl;
				throw " Error in the Lapack routine ZGESVD";
			};

			csvhyst (cout,  "SvFltRhs", kloc, svloc);
			csvhyst (_fout, "SvFltRhs", kloc, svloc);

			nrhsnew = 0;

			for (i=0;i<kloc;i++) {
				if (svloc[i] > _param.thrfltrhs) nrhsnew++;
			};

			for (i=0;i<nrhsnew;i++) {
				for (j=0;j<kloc;j++) {
					uloc[i*kloc+j] = conj(vloc[j*kloc+i]);
				};
			};

		} else {
			nrhsnew = 0;
			for (i=0;i<kloc*kloc;i++) uloc[i] = czero;
		};

		if (_tree.nproc != 1) {
			int nrhsnewloc;
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &nrhsnew, &nrhsnewloc);
			nrhsnew = nrhsnewloc;
		};

		if (_tree.nproc != 1) {

			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														nrhsnew*kloc*2, uloc, vloc);

			for (i=0;i<nrhsnew*kloc;i++) uloc[i] = vloc[i];

		};

// Rotate rhs and initial guess

		CSVectorC solrot(n,npartsloc,nrhsnew), xrot(n,npartsloc,nrhsnew), brot(n,npartsloc,nrhsnew);

		for (ilist=0;ilist<npartsloc;ilist++) {
			solrot.inda[ilist] = inda[ilist];
			xrot.inda[ilist] = inda[ilist];
			brot.inda[ilist] = inda[ilist];
		};
		for (ilist=0;ilist<=npartsloc;ilist++) {
			solrot.blksv[ilist] = blksv[ilist];
			xrot.blksv[ilist] = blksv[ilist];
			brot.blksv[ilist] = blksv[ilist];
		};

		for (irhsnew=0;irhsnew<nrhsnew;irhsnew++) {
			for (i=0;i<n;i++) {
				caux = czero;
				for (j=0;j<nlistrhs;j++) {
					caux += b.vect[j*n+i] * uloc[irhsnew*kloc+j];
				};
				brot.vect[irhsnew*n+i] = caux;
			};
		};

		for (irhsnew=0;irhsnew<nrhsnew;irhsnew++) {
			for (i=0;i<n;i++) {
				caux = czero;
				for (j=0;j<nlistrhs;j++) {
					caux += sol.vect[j*n+i] * uloc[irhsnew*kloc+j];
				};
				solrot.vect[irhsnew*n+i] = caux;
			};
		};

		for (irhsnew=0;irhsnew<nrhsnew;irhsnew++) {
			for (i=0;i<n;i++) {
				caux = czero;
				for (j=0;j<nlistrhs;j++) {
					caux += x.vect[j*n+i] * uloc[irhsnew*kloc+j];
				};
				xrot.vect[irhsnew*n+i] = caux;
			};
		};

// Solve point scaled reduced problem

		CMvmC mvm (false, _tree, _gmtral, nrhsnew);

		if (nrhsnew != 0) {
			if (_param.nrestarts < 1) {

				if (_param.ittype == 2) {
					iconv = solrot.BlockGMRESCorr (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, xrot, brot, _corr, _param);
				} else {
					solrot.BlockLanczosRight (_fout, _tree, mvm, _gmtral, _gmtrau, _gmtrl, _gmtru, xrot, brot, _param);
				};

			} else {
				int icorrloc = 0;
				int ibegrhs = 0;
				int iendrhs, nrhsrestloc;
				CSlvParam paramloc;
				CCorrectorC corrdummy (_corr.ncorrmax);
				_corr = corrdummy;
				while (ibegrhs < nrhsnew) {
					iendrhs = ibegrhs+_param.nrhsrestart-1;
					if (iendrhs > nrhsnew-1) iendrhs = nrhsnew-1;
					nrhsrestloc = iendrhs-ibegrhs+1;
					CSVectorC solrst(n,npartsloc,nrhsrestloc), xrst(n,npartsloc,nrhsrestloc), brst(n,npartsloc,nrhsrestloc);
					CMvmC mvmloc (false, _tree, _gmtral, nrhsrestloc);
					for (ilist=0;ilist<npartsloc;ilist++) {
						solrst.inda[ilist] = inda[ilist];
						xrst.inda[ilist] = inda[ilist];
						brst.inda[ilist] = inda[ilist];
					};
					for (ilist=0;ilist<=npartsloc;ilist++) {
						solrst.blksv[ilist] = blksv[ilist];
						xrst.blksv[ilist] = blksv[ilist];
						brst.blksv[ilist] = blksv[ilist];
					};
					for (irhsnew=ibegrhs;irhsnew<=iendrhs;irhsnew++) {
						for (i=0;i<n;i++) {
							brst.vect[(irhsnew-ibegrhs)*n+i] = brot.vect[irhsnew*n+i];
							xrst.vect[(irhsnew-ibegrhs)*n+i] = xrot.vect[irhsnew*n+i];
							solrst.vect[(irhsnew-ibegrhs)*n+i] = solrot.vect[irhsnew*n+i];
						};
					};
					for (int irestart=0;irestart<_param.nrestarts;irestart++) {
						paramloc = _param;
						if (icorrloc < _param.ncorrmax) {
							paramloc.jobitr = 2;
						} else {
							paramloc.jobitr = 1;
						};
						int iconvrst = solrst.BlockGMRESCorr (_fout, _tree, mvmloc, _gmtral, _gmtrau, _gmtrl, _gmtru, 
																xrst, brst, _corr, paramloc);
						icorrloc++;
						if (iconvrst>0) break;
						xrst = solrst;
					};
					for (irhsnew=ibegrhs;irhsnew<=iendrhs;irhsnew++) {
						for (i=0;i<n;i++) {
							solrot.vect[irhsnew*n+i] = solrst.vect[(irhsnew-ibegrhs)*n+i];
						};
					};
					ibegrhs = iendrhs+1;
				};
				_corr = corrdummy;
			};
		};

// Rotate guess to the solution back

		for (j=0;j<nlistrhs;j++) {
			for (i=0;i<n;i++) {
				caux = czero;
				for (irhsnew=0;irhsnew<nrhsnew;irhsnew++) {
					caux.PlEq (solrot.vect[irhsnew*n+i] * conj(uloc[irhsnew*kloc+j]));
				};
				sol.vect[j*n+i] = caux;
			};
		};

// Free the memory

		delete [] tauloc;
		delete [] rloc;
		delete [] uloc;
		delete [] vloc;
		delete [] work;
		delete [] svloc;
		delete [] rwork;

	};

// Scale solution vector back

	if (_gmtrl.ivect==1) {
		sol.ScaleVect ('D', _gmtrl.vector1);
	} else if (_gmtrl.ivect==2) {
		sol.BlockScaleVect ('T', _gmtral, _gmtru.vector1);
	};

//
// Store the solution
//
	SetSVect(czero);

	for (ilist=0;ilist<nlistrhs;ilist++) {
		irhs = listrhs[ilist];
		for (int i=0;i<n;i++) {
			vect[irhs*n+i] = sol.vect[ilist*n+i];
		};
	};

	delete [] listrhs;
	delete [] list0;
	delete [] list1;

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_tree.myid == 0) {

		cout  << " Solver statistics: " << endl;
		_fout << " Solver statistics: " << endl;

		cout  << "     Time  = " << tottim << " sec." << endl;
		_fout << "     Time  = " << tottim << " sec." << endl;

	};

	return iconv;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve linear system
//========================================================================================
void CSVectorC::SolveA (ofstream &_fout, // Solve linear system
									CGSMatrixCS &_amatr,
									CSVectorC &_b, 
									CSlvParam &_param) {

	const char *funcname = "SolveA";

	dcmplx czero (0.0e0,0.0e0);

// Compute sparsified matrix

	int nloc = _amatr.GetN ();
	int mloc = _amatr.GetM ();

	if (nloc != mloc) throw " CSVectorC::SolveA: rectangular matrix on entry ";

	int nitermax = 10;
	int ndiag = _param.ndiag;

	CSMatrixCS aloc;

	aloc = _amatr.Sparsify (nitermax, nloc*ndiag);
//	asp.A2Ps (4,"ASp.ps",0,&i);

// Compute optimal ordering of AtA

	int *order;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);

	aloc.OrderPrfMtr (1, order);

// Reorder the matrix

	CSMatrixCS ao;
	int i;

	if (true) {
		ao = aloc.OrdMtr (order);
	} else {
		ao = aloc;
		for (i=0;i<nloc;i++) order[i] = i;
	};
//	aspo.A2Ps (4,"AColOrd.ps",0,&i);

	CSMatrixCS mtrdummycs;

	aloc = mtrdummycs;

// Split matrix into L and U parts

	CSMatrixCS mtral;
	CSMatrixCS mtrau;

	ao.CombineLU (mtral, mtrau);

// Free 

	ao = mtrdummycs;

// Compute incomplete LU factorization

	CSMatrixCS mtrl;
	CSMatrixCS mtru;

	mtrdummycs.Ilu2 (_fout, _param,
							mtral, mtrau,
							mtrl, mtru);

	mtral = mtrdummycs;
	mtrau = mtrdummycs;

// Solve the system

	int nrhsloc = _b.GetNrhs ();

	CSVectorC x(nloc,nrhsloc);

	x.SetSVect (czero);

	CSVectorC sol(nloc,nrhsloc);

//	sol.BlockGMRES (_fout,
//						_amatr, order, mtrl, mtru,
//						x, _b, 
//						_param);

// Create dummy tree

	int nproc = 1;
	int nchilds = 3;

	double perf2cpu[100];
	double memory2cpu[100];

	for (i=0;i<100;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (nproc,nchilds,perf2cpu,memory2cpu);

	CNode *pnodes = tree.GetNode ();

	pnodes[0].SetIndbeg (0);
	pnodes[0].SetIndend (0);
	pnodes[0].SetIndbegtot (0);
	pnodes[0].SetIndendtot (0);

	sol.BlockGMRES (_fout,
						tree, _amatr, order, mtrl, mtru,
						x, _b, 
						_param);

	delete [] order;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve linear system
//========================================================================================
void CSVectorC::SolveA (ofstream &_fout, // Solve linear system
									void *_obj, CGETM _getm, CGETSPARSITYCS _getsparsitycs, CMVMC _mvmc,
									CSVectorC &_b, 
									CSlvParam &_param) {

	const char *funcname = "SolveA";

	dcmplx czero (0.0e0,0.0e0);

// Compute sparsified matrix

	int nloc = (_getm) (_obj);;

	CSMatrixCS aloc;

	aloc = (_getsparsitycs) (_obj, _param);
//	asp.A2Ps (4,"ASp.ps",0,&i);

// Compute optimal ordering of AtA

	int *order;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);

	aloc.OrderPrfMtr (1, order);

// Reorder the matrix

	CSMatrixCS ao;
	int i;

	if (true) {
		ao = aloc.OrdMtr (order);
	} else {
		ao = aloc;
		for (i=0;i<nloc;i++) order[i] = i;
	};
//	aspo.A2Ps (4,"AColOrd.ps",0,&i);

	CSMatrixCS mtrdummycs;

	aloc = mtrdummycs;

// Split matrix into L and U parts

	CSMatrixCS mtral;
	CSMatrixCS mtrau;

	ao.CombineLU (mtral, mtrau);

// Free 

	ao = mtrdummycs;

// Compute incomplete LU factorization

	CSMatrixCS mtrl;
	CSMatrixCS mtru;

	mtrdummycs.Ilu2 (_fout, _param,
							mtral, mtrau,
							mtrl, mtru);

	mtral = mtrdummycs;
	mtrau = mtrdummycs;

// Solve the system

	int nrhsloc = _b.GetNrhs ();

	CSVectorC x(nloc,nrhsloc);

	x.SetSVect (czero);

	CSVectorC sol(nloc,nrhsloc);

//	sol.BlockGMRES (_fout,
//						_amatr, order, mtrl, mtru,
//						x, _b, 
//						_param);

// Create dummy tree

	int nproc = 1;
	int nchilds = 3;

	double perf2cpu[100];
	double memory2cpu[100];

	for (i=0;i<100;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (nproc,nchilds,perf2cpu,memory2cpu);

	CNode *pnodes = tree.GetNode ();

	pnodes[0].SetIndbeg (0);
	pnodes[0].SetIndend (0);
	pnodes[0].SetIndbegtot (0);
	pnodes[0].SetIndendtot (0);

	CMvmCLUOrd mvmcluord (&mtrl,&mtru,order);

	sol.BlockGMRES (_fout, _param,
						tree, 
						_obj, _mvmc,
						&mvmcluord, CMvmCLUOrd::SolveLUOrd,
						x, _b);

	delete [] order;

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve rectangular least squares problem with diagonal shift
//========================================================================================
void CSVectorC::SolveAhA (ofstream &_fout, double _dshift, // Solve rectangular least squares problem with diagonal shift
							CGSMatrixCS &_amatr,
							CSVectorC &_b, 
							CSlvParam &_param) {

	const char *funcname = "SolveAhA";

	dcmplx czero (0.0e0,0.0e0);

// Compute sparsified matrix

	int nloc = _amatr.GetN ();
//	int mloc = _amatr.GetM ();

	int nitermax = 10;
	int ndiag = _param.ndiag;

	CSMatrixCS asp;

	asp = _amatr.Sparsify (nitermax, nloc*ndiag);
//	asp.A2Ps (4,"ASp.ps",0,&i);

// Compute sparsity of AtA

	CSMatrix atasp;

	atasp = asp.AtAMatrSpars ();
//	atasp.A2Ps (4,"AtASp.ps",0,&i);

// Compute optimal ordering of AtA

	int *order;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);

	atasp.OrderPrfMtr (1, order);

	CSMatrix mtrdummy;

	atasp = mtrdummy;

// Reorder columns of the matrix

	CSMatrixCS aspo;
	int i;

	if (true) {
		aspo = asp.OrdMtrRectCol (order);
	} else {
		aspo = asp;
		for (i=0;i<nloc;i++) order[i] = i;
	};
//	aspo.A2Ps (4,"AColOrd.ps",0,&i);

	CSMatrixCS mtrdummycs;

	asp = mtrdummycs;

// Compute ata matrix (its upper triangular part)

	CSMatrixCS ata;

	ata = aspo.AtAMatrElem (_fout, _param);

	aspo = mtrdummycs;

// Conjugate ata matrix

	int nzataloc = ata.GetNza ();
	dcmplx *pata = ata.GetA ();

	dcmplx caux;

	for (i=0;i<nzataloc;i++) {
		caux = pata[i];
		pata[i] = conj(caux);
	};
//	ata.A2Ps (4,"AtASpOrd.ps",0,&i);

// Prepare CGSMatrixCS matrix

	int blksloc[2];

	blksloc[0] = 0;
	blksloc[1] = nloc;

	int nblksloc = 1;

	CGSMatrixCS gata;

	gata.CS2GCS (ata, nblksloc, blksloc);

	ata = mtrdummycs;

// Compute ICH2 factorization of the AtA matrix

	CGSMatrixCS gmtru;

	CGSMatrixCS gmtrdummycs;

	cout << " Before Ich2 " << endl;

	int nproc = 1;
	int nchilds = 3;

	double perf2cpu[100];
	double memory2cpu[100];

	for (i=0;i<100;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (nproc,nchilds,perf2cpu,memory2cpu);

	CNode *pnodes = tree.GetNode ();

	pnodes[0].SetIndbeg (0);
	pnodes[0].SetIndend (0);
	pnodes[0].SetIndbegtot (0);
	pnodes[0].SetIndendtot (0);

	gmtrdummycs.Ich2Schur (_fout, tree, _param,
								gata,
								gmtru);

// Solve

	CSVectorC x(nloc);
	CSVectorC b(nloc);
	CSVectorC bh(nloc);

	dcmplx *pb = b.GetVect ();
	dcmplx *pbh = bh.GetVect ();

	_amatr.MvmAhRect  (_b, bh);

	for (i=0;i<nloc;i++) {
		pb[order[i]] = pbh[i];
	};

	x.SetSVect (czero);

	CSVectorC sol(nloc);

	dcmplx *psol = sol.GetVect ();

	int itype = 1;

	sol.BlockGMRES (_fout,
						itype, 
						gata, 
						_dshift, _amatr, order,
						gmtru,
						x, b, 
						_param);

	int inew;

	for (i=0;i<nloc;i++) {
		inew = order[i];
		pb[i] = psol[inew];
	};

	delete [] order;

	*this = b;

};

// Author: Kharchenko S.A.
// CSVectorC: Solve rectangular least squares problem with diagonal shift
//========================================================================================
void CSVectorC::SolveAhA (ofstream &_fout, double _dshift, // Solve rectangular least squares problem with diagonal shift
							void *_obj, CGETM _getm, CGETN _getn, CGETSPARSITYCS _getsparsitycs, CMVMC _mvmc, CMVMHC _mvmhc,
							CSVectorC &_b, 
							CSlvParam &_param) {

	const char *funcname = "SolveAhA";

	dcmplx czero (0.0e0,0.0e0);

// Compute sparsified matrix

	int nloc = (_getn) (_obj);;
//	int mloc = (_getm) (_obj);;

	CSMatrixCS asp;

	asp = (_getsparsitycs) (_obj, _param);
//	asp.A2Ps (4,"ASp.ps",0,&i);

// Compute sparsity of AtA

	CSMatrix atasp;

	atasp = asp.AtAMatrSpars ();
//	atasp.A2Ps (4,"AtASp.ps",0,&i);

// Compute optimal ordering of AtA

	int *order;

	order = new int [nloc];
	if (!order) MemoryFail (funcname);

	atasp.OrderPrfMtr (1, order);

	CSMatrix mtrdummy;

	atasp = mtrdummy;

// Reorder columns of the matrix

	CSMatrixCS aspo;
	int i;

	if (true) {
		aspo = asp.OrdMtrRectCol (order);
	} else {
		aspo = asp;
		for (i=0;i<nloc;i++) order[i] = i;
	};
//	aspo.A2Ps (4,"AColOrd.ps",0,&i);

	CSMatrixCS mtrdummycs;

	asp = mtrdummycs;

// Compute ata matrix (its upper triangular part)

	CSMatrixCS ata;

	ata = aspo.AtAMatrElem (_fout, _param);

	aspo = mtrdummycs;

// Conjugate ata matrix

	int nzataloc = ata.GetNza ();
	dcmplx *pata = ata.GetA ();

	dcmplx caux;

	for (i=0;i<nzataloc;i++) {
		caux = pata[i];
		pata[i] = conj(caux);
	};
//	ata.A2Ps (4,"AtASpOrd.ps",0,&i);

// Prepare CGSMatrixCS matrix

	int blksloc[2];

	blksloc[0] = 0;
	blksloc[1] = nloc;

	int nblksloc = 1;

	CGSMatrixCS gata;

	gata.CS2GCS (ata, nblksloc, blksloc);

	ata = mtrdummycs;

// Compute ICH2 factorization of the AtA matrix

	CGSMatrixCS gmtru;

	CGSMatrixCS gmtrdummycs;

	cout << " Before Ich2 " << endl;

	int nproc = 1;
	int nchilds = 3;

	double perf2cpu[100];
	double memory2cpu[100];

	for (i=0;i<100;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree tree (nproc,nchilds,perf2cpu,memory2cpu);

	CNode *pnodes = tree.GetNode ();

	pnodes[0].SetIndbeg (0);
	pnodes[0].SetIndend (0);
	pnodes[0].SetIndbegtot (0);
	pnodes[0].SetIndendtot (0);

	gmtrdummycs.Ich2Schur (_fout, tree, _param,
								gata,
								gmtru);

// Solve

	CSVectorC x(nloc);
	CSVectorC b(nloc);
	CSVectorC bh(nloc);

	dcmplx *pb = b.GetVect ();
	dcmplx *pbh = bh.GetVect ();

	(_mvmhc) (_obj, _b, bh);

	for (i=0;i<nloc;i++) {
		pb[order[i]] = pbh[i];
	};

	x.SetSVect (czero);

	CSVectorC sol(nloc);

	dcmplx *psol = sol.GetVect ();

//	int itype = 1;

	sol.BlockGMRES (_fout,
						_obj, _getm, _getn, _mvmc, _mvmhc,_dshift, order,
						gmtru,
						x, b, 
						_param);

	int inew;

	for (i=0;i<nloc;i++) {
		inew = order[i];
		pb[i] = psol[inew];
	};

	delete [] order;

	*this = b;

};

// Author: Kharchenko S.A.
// CSVector: Solve symmetric system
//========================================================================================
void CSVector::SolveSymm (ofstream &_fout, CSMatrixRS &_mtra, // Solve symmetric system
						const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm) {

	const char *funcname = "SolveSymm";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;

//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();
	int nsupa = _mtra.GetNsupr();

	if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;
	int *sprndcini, *sprndcord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [nsupa];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupa];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupa+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupa+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
		_mtra.PartOrdMtr (nparts, blks, order);
	} else if (_param.ordtype == -1) {
		int ordtyp = -1;
		_mtra.OrderPrfMtr (ordtyp, order);
	} else if (_param.ordtype == 0) {
		for (i=0;i<_mtra.nsupr;i++) order[i] = i;
	} else if (_param.ordtype < -1) {
		_mtra.OrderBrdMtr (nparts, order);
	};

	for (i=0;i<_mtra.nsupr;i++) iord[order[i]] = i;

	for (i=0;i<=nsupa;i++) sprndcini[i] = _mtra.sprndc[i];
//
// Reorder the matrix in-place
//
	CSMatrixRS mtrao;

	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtrSymm (order);

		_mtra = mtrao;

	};

	CSMatrixRS mtrrdummy;

	mtrao = mtrrdummy;
	if (collap != 0) _mtra.A2Ps (collap,"AStrOrd.ps",0,&i);

	for (i=0;i<=nsupa;i++) sprndcord[i] = _mtra.sprndc[i];
//
// Compute Ich2
//
	CSMatrixRS ich;

	ich = _mtra.Ich2 (_fout, _param);

	if (collap != 0) ich.A2Ps (collap,"StrIch.ps",0,&i);
//
// Allocate vector data
//
	CSVector x(n), solo(n);
	CSVector bo(n), normo(n);
//
// Solve
//
	bo = _b;
	normo = _norm;
	x = *this;

	if (_param.ordtype != 0) {

		bo.OrdVect('D', nsupa, sprndcini, sprndcord, order);
		normo.OrdVect('D', nsupa, sprndcini, sprndcord, order);
		x.OrdVect('D', nsupa, sprndcini, sprndcord, order);

	};

	int istore = 0;
	CGSMatrixR gich;

	solo.Cg (_fout, &_mtra, istore, &gich, &ich, x, bo, _param, normo);
//
// Reorder the matrix back
//
	if (_param.ordtype != 0) {

		mtrao = _mtra.OrdMtrSymm (iord);

		_mtra = mtrao;

	};

	mtrao = mtrrdummy;
//
// Return the solution
//
	if (_param.ordtype != 0) {
		solo.OrdVect('I', nsupa, sprndcini, sprndcord, order);
	};

	*this = solo;
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	cout  << " Symm solver statistics: " << endl;
	_fout << " Symm solver statistics: " << endl;

	cout  << "     Time  = " << tottim << " sec." << endl;
	_fout << "     Time  = " << tottim << " sec." << endl;

};

// Author: Kharchenko S.A.
// CSVector: Solve symmetric system
//========================================================================================
void CSVector::SolveSymm (ofstream &_fout, CSMatrixR &_mtra, // Solve symmetric system
						const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm) {

	const char *funcname = "SolveSymm";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;

//
// Get the sizes of the matrix
//
	int i = 0;

	int n = _mtra.GetN();

	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrIni.ps",0,&i);
	};
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [n];
	if (!order) MemoryFail (funcname);
	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
//		_mtra.PartOrdMtr (nparts, blks, order);
		int ordtyp = 1;
		_mtra.OrderPrfMtr (ordtyp, order);
	} else if (_param.ordtype == -1) {
		int ordtyp = -1;
		_mtra.OrderPrfMtr (ordtyp, order);
	} else if (_param.ordtype == 0) {
		for (i=0;i<n;i++) order[i] = i;
	} else if (_param.ordtype < -1) {
		_mtra.OrderBrdMtr (nparts, order);
	};

	for (i=0;i<n;i++) iord[order[i]] = i;
//
// Reorder the matrix in-place
//
	CSMatrixR mtrao;

	mtrao = _mtra.OrdMtrSymm (order);

	_mtra = mtrao;

	CSMatrixR mtrrdummy;

	mtrao = mtrrdummy;

	if (_param.msglev > 1) {
		if (collap != 0) _mtra.A2Ps (collap,"AStrOrd.ps",0,&i);
	};
//
// Compute Ich/Ich2
//
	int istore = 1;
//	int istore = 0;

	CSMatrixR ich;
	CGSMatrixR gich;

	if (_param.fcttyp == 1) {
		ich = _mtra.Ich ();
	} else {
		if (_param.memory > 0.0e0) {
			ich = _mtra.Ich2 (_param);
		} else if (_param.memory == 0.0e0) {
			ich = _mtra.Ich2Dynamic (_fout, _param);
		} else if (_param.memory < 0.0e0) {
			_mtra.Ich2DynamicByBlocks (_fout, _param,
						istore, gich, ich);
		};
	};

	if (_param.msglev > 1) {
		if (collap != 0) ich.A2Ps (collap,"StrIch.ps",0,&i);
	};
//
// Allocate vector data
//
	CSVector x(n), solo(n);
	CSVector bo(n), normo(n);
//
// Solve
//
	x.OrdVect ('D', order, *this);
	bo.OrdVect ('D', order, _b);
	normo.OrdVect ('D', order, _norm);

	solo.Cg (_fout, &_mtra, istore, &gich, &ich, x, bo, _param, normo);

	ich = mtrrdummy;

	CGSMatrixR gmtrdummy;

	gich = gmtrdummy;
//
// Reorder the matrix back
//
	mtrao = _mtra.OrdMtrSymm (iord);

	_mtra = mtrao;

	mtrao = mtrrdummy;
//
// Return the solution
//
	OrdVect ('I', order, solo);
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	if (_param.msglev > 1) {
		cout  << " Symm solver statistics: " << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Symm solver statistics: " << endl;
	};

	if (_param.msglev > 1) {
		cout  << "     Time  = " << tottim << " sec." << endl;
	};
	if (_param.msglev >= 1) {
		_fout << "     Time  = " << tottim << " sec." << endl;
	};

};

// Author: Kharchenko S.A.
// CSVector: Solve Lsq problem
//========================================================================================
void CSVector::SolveLsq (ofstream &_fout, CSMatrixRS &_mtra, // Solve Lsq problem
							const CSVector &_b, 
							CSlvParam &_param) {

	const char *funcname = "SolveLsq";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;

//
// Get the sizes of the matrix
//
	int i = 0;

	int m = _mtra.GetM();
	int n = _mtra.GetN();
	int nsupc = _mtra.GetNsupc();
//	int nsupr = _mtra.GetNsupr();

//	_mtra.A2Ps (collap,"AStrIni.ps",0,&i);
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

//	int nza = _mtra.GetNza();
//
// Compute the sparsity of the AtA matrix
//
	CSMatrix atastr;

	atastr = _mtra.AtAMatrSpars ();

	if (collap != 0) atastr.A2Ps (collap,"AtAIni.ps",0,&i);
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;
	int *sprndcini, *sprndcord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [nsupc];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupc];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupc+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupc+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
		atastr.PartOrdMtr (nparts, blks, order);
	} else {
		atastr.OrderBrdMtr (nparts, order);
	};

	for (i=0;i<nsupc;i++) iord[order[i]] = i;

	for (i=0;i<=nsupc;i++) sprndcini[i] = _mtra.sprndc[i];
//// 
// Free AtA sparsity
//
	CSMatrix mtrdummy;

	atastr = mtrdummy;
//
// Reorder the columns of the matrix
//
	CSMatrixRS mtrao;

	mtrao = _mtra.OrdMtrRectCol (order);

	_mtra = mtrao;

//	for (i=0;i<nsupc;i++) order[i] = i;

	CSMatrixRS mtrrsdummy;

	mtrao = mtrrsdummy;

	for (i=0;i<=nsupc;i++) sprndcord[i] = _mtra.sprndc[i];
//
// Compute AtA matrix
//
	CSMatrixRS ata;

	ata = _mtra.AtAMatrElem (_fout);

	if (collap != 0) ata.A2Ps (collap,"AtAOrd.ps",0,&i);
//
// Compute Ich2
//
	CSMatrixRS ich;

	ich = ata.Ich2(_fout, _param);

	if (collap != 0) ich.A2Ps (collap,"StrIch.ps",0,&i);
//
// Free AtA matrix
//
	ata = mtrrsdummy;
//
// Solve
//
	CSMatrixRS *pa, *pich;

	pa = &_mtra;
	pich = &ich;
//
// Switch solution cases
//
	int nrhsloc = GetNrhs ();

	if (_param.ittype == 1) {
//
// Cycle over the right hand sides
//
		CSVector sol(n), x(n), b(m);

		for (int irhs=0;irhs<nrhsloc;irhs++) {

			cout  << endl;
			cout  << "   Solve Irhs = " << irhs << endl;
			cout  << endl;
			_fout << endl;
			_fout << "   Solve Irhs = " << irhs << endl;
			_fout << endl;
//
// Init right hand side and initial guess
//
			for (i=0;i<m;i++) {
				b.vect[i] = _b.vect[irhs*m+i];
			};
			for (i=0;i<n;i++) {
				x.vect[i] = vect[irhs*n+i];
			};
			sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
			x.OrdVect('D', nsupc, sprndcini, sprndcord, order);

			double aux=0.0e0;

			aux = b.ScProdSVect (b);

			if (aux > 0.0e0) sol.LanczosLsq (_fout, pa, pich, x, b, _param);
//
// Reorder solution back
//
			sol.OrdVect('I', nsupc, sprndcini, sprndcord, order);

			x = sol;
//
// Store the solution
//
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[i];
			};

		};

	} else {
//
// Determine the number of nonzero right hand sides and store the list
//
		int nlistrhs = 0;

		int *listrhs;

		listrhs = new int [nrhsloc];
		if (!listrhs) MemoryFail (funcname);

		int irhs;
		for (irhs=0;irhs<nrhsloc;irhs++) {

			double aux = 0.0e0;

			for (i=0;i<m;i++) {
				double auxl = _b.vect[irhs*m+i];
				aux += auxl*auxl;
			};

			if (aux != 0.0e0) {
				listrhs[nlistrhs] = irhs;
				nlistrhs++;
			};

		};
//
// Solve by block Lanczos
//
		CSVector sol(n,1,nlistrhs), x(n,1,nlistrhs), b(m,1,nlistrhs);
//
// Init right hand side and initial guess
//
		int ilist;
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<m;i++) {
				b.vect[ilist*m+i] = _b.vect[irhs*m+i];
			};
			for (i=0;i<n;i++) {
				x.vect[ilist*n+i] = vect[irhs*n+i];
			};
		};

		sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
		x.OrdVect('D', nsupc, sprndcini, sprndcord, order);

		sol.BlockLanczosLsq (_fout, pa, pich, x, b, _param);
//
// Reorder solution back
//
		sol.OrdVect('I', nsupc, sprndcini, sprndcord, order);

		x = sol;
//
// Store the solution
//
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[ilist*n+i];
			};
		};

		delete [] listrhs;

	};
//
// Reorder the matrix back
//
	mtrao = _mtra.OrdMtrRectCol (iord);

	_mtra = mtrao;

	mtrao = mtrrsdummy;
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	cout  << " Lsq Solver statistics: " << endl;
	_fout << " Lsq Solver statistics: " << endl;

	cout  << "     Time  = " << tottim << " sec." << endl;
	_fout << "     Time  = " << tottim << " sec." << endl;

};

// Author: Kharchenko S.A.
// CSVector: Solve Lsq problem in the out-of-core mode
//========================================================================================
void CSVector::SolveLsq (ofstream &_fout, // Solve Lsq problem in the out-of-core mode
							int _nblks, int *_blks, FILE **_mtrfiles,
							const CSMatrixRS &_mtra,
							char *_path,
							const CSVector &_b,
							CSlvParam &_param) {

	const char *funcname = "SolveLsq";

	int collap = _param.collap;

	double ops, tottim;
	clock_t time0, time1;
//
// Flush files with the matrix data
//
	FlushDIOFiles (_nblks, _mtrfiles);
//
// Create arrays of files
//
	FILE **mtraofiles, **mtratafiles, **mtrlufiles;

	mtraofiles = new FILE * [_nblks];
	if (!mtraofiles) MemoryFail (funcname);
	mtratafiles = new FILE * [_nblks];
	if (!mtratafiles) MemoryFail (funcname);
	mtrlufiles = new FILE * [_nblks];
	if (!mtrlufiles) MemoryFail (funcname);

	OpenDIOFiles (_nblks, mtraofiles,  _path, "AoMatr_", ".bin");
	OpenDIOFiles (_nblks, mtratafiles, _path, "AtAMatr_",".bin");
	OpenDIOFiles (_nblks, mtrlufiles,  _path, "LUMatr_", ".bin");
//
// Get the sizes of the matrix
//
	int i = 0;

	int m = _mtra.GetM();
	int n = _mtra.GetN();
	int nsupc = _mtra.GetNsupc();
//	int nsupr = _mtra.GetNsupr();

//	_mtra.A2Ps (collap,"AStrIni.ps",0,&i);
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;
//
// Compute the sparsity of the AtA matrix
//
	CSMatrix atastr;

	atastr = _mtra.AtAMatrSpars ();

	if (collap != 0) atastr.A2Ps (collap,"AtAInOut.ps",0,&i);
//
// Compute new ordering
//
	int nparts = _param.ordtype;
	if (nparts < 0) nparts = -nparts;

	int *blks, *order, *iord;
	int *sprndcini, *sprndcord;

	blks = new int [nparts+2];
	if (!blks) MemoryFail (funcname);
	order = new int [nsupc];
	if (!order) MemoryFail (funcname);
	iord = new int [nsupc];
	if (!iord) MemoryFail (funcname);
	sprndcini = new int [nsupc+1];
	if (!sprndcini) MemoryFail (funcname);
	sprndcord = new int [nsupc+1];
	if (!sprndcord) MemoryFail (funcname);

	if (_param.ordtype > 0) {
		atastr.PartOrdMtr (nparts, blks, order);
	} else {
		atastr.OrderBrdMtr (nparts, order);
	};

	for (i=0;i<nsupc;i++) iord[order[i]] = i;

	for (i=0;i<=nsupc;i++) sprndcini[i] = _mtra.sprndc[i];
// 
// Free AtA sparsity
//
	CSMatrix mtrdummy;

	atastr = mtrdummy;
//
// Reorder the columns of the matrix
//
	CSMatrixRS mtrao;

	mtrao = _mtra.OrdMtrRectCol (_nblks, _blks, _mtrfiles, 
									order, mtraofiles); 

//	for (i=0;i<nsupc;i++) order[i] = i;

	for (i=0;i<=nsupc;i++) sprndcord[i] = mtrao.sprndc[i];
//
// Flush files with the reordered matrix data
//
	FlushDIOFiles (_nblks, mtraofiles);
//
// Compute AtA matrix
//
	CSMatrixRS ata;

	ata = mtrao.AtAMatrElem (_fout, _nblks, _blks, mtraofiles, 
								mtratafiles); 

	if (collap != 0) ata.A2Ps (collap,"AtAOrOut.ps",0,&i);
//
// Flush files with the AtA matrix data
//
	FlushDIOFiles (_nblks, mtratafiles);
//
// Compute Ich2
//
	CSMatrixRS ich;

	ich = ata.Ich2(_fout, 
					_nblks, _blks, 
					mtratafiles, mtrlufiles, 
					_param);

	if (collap != 0) ich.A2Ps (collap,"StIchOut.ps",0,&i);
//
// Flush files with the Ich data
//
	FlushDIOFiles (_nblks, mtrlufiles);
//
// Free AtA matrix
//
	CSMatrixRS mtrrsdummy;

	ata = mtrrsdummy;
//
// Close ata files on disks
//
	CloseDIOFiles (_nblks, mtratafiles);
//
// Solve
//
	CSMatrixRS *pa, *pich;

	pa = &mtrao;
	pich = &ich;
//
// Switch solution cases
//
	int nrhsloc = GetNrhs ();

	if (_param.ittype == 1) {
//
// Cycle over the right hand sides
//
		CSVector sol(n), x(n), b(m);

		for (int irhs=0;irhs<nrhsloc;irhs++) {

			cout  << endl;
			cout  << "   Solve Irhs = " << irhs << endl;
			cout  << endl;
			_fout << endl;
			_fout << "   Solve Irhs = " << irhs << endl;
			_fout << endl;
//
// Init right hand side and initial guess
//
			for (i=0;i<m;i++) {
				b.vect[i] = _b.vect[irhs*m+i];
			};
			for (i=0;i<n;i++) {
				x.vect[i] = vect[irhs*n+i];
			};
			sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
			x.OrdVect('D', nsupc, sprndcini, sprndcord, order);

			double aux=0.0e0;

			aux = b.ScProdSVect (b);

			if (aux > 0.0e0) sol.LanczosLsq (_fout, 
												_nblks, _blks, 
												pa, mtraofiles,
												pich, mtrlufiles,
												x, b, _param);
//
// Reorder solution back
//
			sol.OrdVect('I', nsupc, sprndcini, sprndcord, order);

			x = sol;
//
// Store the solution
//
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[i];
			};

		};

	} else {
//
// Determine the number of nonzero right hand sides and store the list
//
		int nlistrhs = 0;

		int *listrhs;

		listrhs = new int [nrhsloc];
		if (!listrhs) MemoryFail (funcname);

		int irhs;
		for (irhs=0;irhs<nrhsloc;irhs++) {

			double aux = 0.0e0;

			for (i=0;i<m;i++) {
				double auxl = _b.vect[irhs*m+i];
				aux += auxl*auxl;
			};

			if (aux != 0.0e0) {
				listrhs[nlistrhs] = irhs;
				nlistrhs++;
			};

		};
//
// Solve by block Lanczos
//
		CSVector sol(n,1,nlistrhs), x(n,1,nlistrhs), b(m,1,nlistrhs);
//
// Init right hand side and initial guess
//
		int ilist;
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<m;i++) {
				b.vect[ilist*m+i] = _b.vect[irhs*m+i];
			};
			for (i=0;i<n;i++) {
				x.vect[ilist*n+i] = vect[irhs*n+i];
			};
		};

		sol.SetSVect(0.0e0);
//
// Reorder initial guess
//
		x.OrdVect('D', nsupc, sprndcini, sprndcord, order);

		sol.BlockLanczosLsq (_fout, 
								_nblks, _blks, 
								pa, mtraofiles,
								pich, mtrlufiles,
								x, b, _param);
//
// Reorder solution back
//
		sol.OrdVect('I', nsupc, sprndcini, sprndcord, order);

		x = sol;
//
// Store the solution
//
		for (ilist=0;ilist<nlistrhs;ilist++) {
			irhs = listrhs[ilist];
			for (i=0;i<n;i++) {
				vect[irhs*n+i] = x.vect[ilist*n+i];
			};
		};

		delete [] listrhs;

	};
//
// Close ao and ich files on disks
//
	CloseDIOFiles (_nblks, mtraofiles);
	CloseDIOFiles (_nblks, mtrlufiles);
//
// Free work arrays
//
	delete [] blks;
	delete [] order;
	delete [] iord;
	delete [] sprndcini;
	delete [] sprndcord;
	delete [] mtraofiles;
	delete [] mtratafiles;
	delete [] mtrlufiles;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output solver statistics

	cout  << " Lsq Solver statistics: " << endl;
	_fout << " Lsq Solver statistics: " << endl;

	cout  << "     Time  = " << tottim << " sec." << endl;
	_fout << "     Time  = " << tottim << " sec." << endl;

};

// Author: Kharchenko S.A.
// Solve real block symmetric positive system of equations in abstract interface
//========================================================================================
void AbstractSymmSolver (void *_pgener, // Solve real block symmetric positive system of equations in abstract interface
									ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol) {

	const char *funcname = "AbstractSymmSolver";

// Create the matrix

	int nsuploc;

	(_fnsup) (_pgener, nsuploc);

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	(_fsprnds) (_pgener, sprndsloc);

	int nloc = sprndsloc[nsuploc];

	int *ialoc, *jaloc;

	(_fiaja) (_pgener, ialoc, jaloc);

	int nzjaloc = ialoc[nsuploc];

	int i, j, jj, blai, blaj, nzrow;

	int nzaloc = 0;
	int nzrowmax = 0;
	int blamaxloc = 0;

	for (i=0;i<nsuploc;i++) {
		nzrow = 0;
		blai = sprndsloc[i+1]-sprndsloc[i];
		if (blai > blamaxloc) blamaxloc = blai;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
		nzaloc += nzrow;
	};

	CSMatrixRS aloc(nsuploc, nzjaloc, nzaloc);

	aloc.SetM (nloc);
	aloc.SetN (nloc);

	aloc.SetNlist (nsuploc);
	aloc.SetNsupr (nsuploc);
	aloc.SetNsupcBase (nsuploc);

	aloc.SetBlamx (blamaxloc);

	int *psprndr, *psprndc, *pbsa;

	aloc.GetPartition (psprndr, psprndc, pbsa);

	for (i=0;i<=nsuploc;i++) psprndr[i] = sprndsloc[i];
	for (i=0;i<=nsuploc;i++) psprndc[i] = sprndsloc[i];

	int *plist, *pia, *pja;

	aloc.GetSparsity (plist, pia, pja);

	for (i=0;i<nsuploc;i++) plist[i] = i;
	for (i=0;i<=nsuploc;i++) pia[i] = ialoc[i];

	int *bsaloc, *iord;
	CIntInt *iarr;
	double *arow;

	bsaloc = new int [nsuploc];
	if (!bsaloc) MemoryFail (funcname);
	iord = new int [nsuploc];
	if (!iord) MemoryFail (funcname);
	iarr = new CIntInt [nsuploc];
	if (!iarr) MemoryFail (funcname);
	arow = new double [nzrowmax];
	if (!arow) MemoryFail (funcname);

	int nzj, ind, blaij, ibs, k;

	double *pa;
	pa = aloc.GetA();

	nzaloc = 0;

	for (i=0;i<nsuploc;i++) {

		nzj = ialoc[i+1]-ialoc[i];
		blai = sprndsloc[i+1]-sprndsloc[i];

		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			iarr[j-ialoc[i]].intvalue = jj;
			iarr[j-ialoc[i]].int2value = j-ialoc[i];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			bsaloc[jj] = nzrow;
			nzrow += blai*blaj;
		};

// Generate

		(_fblockrow) (_pgener, i, arow);

// Sort data within block row

		qsort (iarr, nzj, sizeof (CIntInt), CIntInt::CompareIntInt);

//		for (j=0;j<nzj;j++) iord[iarr[j].int2value] = j;
		for (j=0;j<nzj;j++) iord[j] = iarr[j].int2value;

		nzrow = 0;

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			ind = iord[j-ialoc[i]];
			jj = jaloc[ialoc[i]+ind];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			blaij = blai*blaj;
			ibs = bsaloc[jj];
			pja[j] = jj;
			for (k=0;k<blaij;k++) {
				pa[nzaloc+nzrow+k] = arow[ibs+k];
			};
			nzrow += blai*blaj;
		};

// Store bsa array

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = pja[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			pbsa[j] = nzaloc;
			nzaloc += blai*blaj;
		};
	};

	aloc.SetNzja (nzjaloc);
	aloc.SetNza (nzaloc);
	aloc.SetNzatot (nzaloc);

// Store rhs and initial guess

	CSVector b(nloc), sol(nloc), x(nloc), norm(nloc);

	double *pvect = b.GetVect ();

	int nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_frhs) (_pgener, i, pvect+nz);
		nz += blai;
	};

	pvect = x.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_finiguess) (_pgener, i, pvect+nz);
		nz += blai;
	};

	norm.SetSVect (1.0e0);

	sol = x;

// Solve the system by the solver

	sol.SolveSymm (_fout, aloc,
						b, 
						_param, norm);

// Store the solution via interface function

	pvect = sol.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_fstoresol) (_pgener, i, pvect+nz);
		nz += blai;
	};

// Free work arrays

	delete [] sprndsloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;
	delete [] iord;
	delete [] iarr;
	delete [] arow;

};

// Author: Kharchenko S.A.
// Solve real block system of equations in abstract interface
//========================================================================================
void AbstractSolver (void *_pgener, // Solve real block system of equations in abstract interface
									ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol) {

	const char *funcname = "AbstractSolver";

// Create the matrix

	int nsuploc;

	(_fnsup) (_pgener, nsuploc);

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	(_fsprnds) (_pgener, sprndsloc);

	int nloc = sprndsloc[nsuploc];

	int *ialoc, *jaloc;

	(_fiaja) (_pgener, ialoc, jaloc);

	int nzjaloc = ialoc[nsuploc];

	int i, j, jj, blai, blaj, nzrow;

	int nzaloc = 0;
	int nzrowmax = 0;
	int blamaxloc = 0;

	for (i=0;i<nsuploc;i++) {
		nzrow = 0;
		blai = sprndsloc[i+1]-sprndsloc[i];
		if (blai > blamaxloc) blamaxloc = blai;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
		nzaloc += nzrow;
	};

	CSMatrixRS aloc(nsuploc, nzjaloc, nzaloc);

	aloc.SetM (nloc);
	aloc.SetN (nloc);

	aloc.SetNlist (nsuploc);
	aloc.SetNsupr (nsuploc);
	aloc.SetNsupcBase (nsuploc);

	aloc.SetBlamx (blamaxloc);

	int *psprndr, *psprndc, *pbsa;

	aloc.GetPartition (psprndr, psprndc, pbsa);

	for (i=0;i<=nsuploc;i++) psprndr[i] = sprndsloc[i];
	for (i=0;i<=nsuploc;i++) psprndc[i] = sprndsloc[i];

	int *plist, *pia, *pja;

	aloc.GetSparsity (plist, pia, pja);

	for (i=0;i<nsuploc;i++) plist[i] = i;
	for (i=0;i<=nsuploc;i++) pia[i] = ialoc[i];

	int *bsaloc, *iord;
	CIntInt *iarr;
	double *arow;

	bsaloc = new int [nsuploc];
	if (!bsaloc) MemoryFail (funcname);
	iord = new int [nsuploc];
	if (!iord) MemoryFail (funcname);
	iarr = new CIntInt [nsuploc];
	if (!iarr) MemoryFail (funcname);
	arow = new double [nzrowmax];
	if (!arow) MemoryFail (funcname);

	int nzj, ind, blaij, ibs, k;

	double *pa;
	pa = aloc.GetA();

	nzaloc = 0;

	for (i=0;i<nsuploc;i++) {

		nzj = ialoc[i+1]-ialoc[i];
		blai = sprndsloc[i+1]-sprndsloc[i];

		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			iarr[j-ialoc[i]].intvalue = jj;
			iarr[j-ialoc[i]].int2value = j-ialoc[i];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			bsaloc[jj] = nzrow;
			nzrow += blai*blaj;
		};

// Generate

		(_fblockrow) (_pgener, i, arow);

// Sort data within block row

		qsort (iarr, nzj, sizeof (CIntInt), CIntInt::CompareIntInt);

//		for (j=0;j<nzj;j++) iord[iarr[j].int2value] = j;
		for (j=0;j<nzj;j++) iord[j] = iarr[j].int2value;

		nzrow = 0;

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			ind = iord[j-ialoc[i]];
			jj = jaloc[ialoc[i]+ind];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			blaij = blai*blaj;
			ibs = bsaloc[jj];
			pja[j] = jj;
			for (k=0;k<blaij;k++) {
				pa[nzaloc+nzrow+k] = arow[ibs+k];
			};
			nzrow += blai*blaj;
		};

// Store bsa array

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = pja[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			pbsa[j] = nzaloc;
			nzaloc += blai*blaj;
		};
	};

	aloc.SetNzja (nzjaloc);
	aloc.SetNza (nzaloc);
	aloc.SetNzatot (nzaloc);

// Store rhs and initial guess

	CSVector b(nloc), sol(nloc), x(nloc), norm(nloc);

	double *pvect = b.GetVect ();

	int nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_frhs) (_pgener, i, pvect+nz);
		nz += blai;
	};

	pvect = x.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_finiguess) (_pgener, i, pvect+nz);
		nz += blai;
	};

	norm.SetSVect (1.0e0);

	sol = x;

// Solve the system by the solver

	sol.Solver (_fout, aloc,
						b, 
						_param, norm);

// Store the solution via interface function

	pvect = sol.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_fstoresol) (_pgener, i, pvect+nz);
		nz += blai;
	};

// Free work arrays

	delete [] sprndsloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;
	delete [] iord;
	delete [] iarr;
	delete [] arow;

};

// Author: Kharchenko S.A.
// Solve real block system of equations in parallel mode in abstract interface
//========================================================================================
void AbstractParSolver (void *_pgener, // Solve real block system of equations in parallel mode in abstract interface
									void *_comm, ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol) {

	const char *funcname = "AbstractParSolver";

	int nproc = ((CMPIComm *)_comm)->GetNproc ();
	int myid = ((CMPIComm *)_comm)->GetMyid ();

// Create the sparsity

	int nsuploc;

	(_fnsup) (_pgener, nsuploc);

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	(_fsprnds) (_pgener, sprndsloc);

	int ntot = sprndsloc[nsuploc];

	int *ialoc, *jaloc;

	(_fiaja) (_pgener, ialoc, jaloc);

	int nzjaloc = ialoc[nsuploc];

// Pack the sparsity into the matrix

	CSMatrix asploc (nsuploc,nzjaloc);

	int *plistsp = asploc.GetList ();
	int *piasp = asploc.GetIa ();
	int *pjasp = asploc.GetJa ();

	int i;

	for (i=0;i<nsuploc;i++) plistsp[i] = i;
	for (i=0;i<=nsuploc;i++) piasp[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) pjasp[i] = jaloc[i];

	asploc.SetM (nsuploc);
	asploc.SetN (nsuploc);
	asploc.SetNsupc(nsuploc);
	asploc.SetNsupr (nsuploc);
	asploc.SetNlist (nsuploc);

	CSMatrix asymm = asploc.SymmMtr ();

//	_fout << " Asymm = " << asymm << endl;

	CSMatrix mtrdummy;

	asploc = mtrdummy;

// Create the binary tree, new ordering and partitioning

	int nprocloc = nproc;
	int nchilds = 2;

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nprocloc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocloc+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

	CTree tree (nprocloc, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	tree.SortNodes ();

	tree.CpuListsForNodes ();

	int *order;

	order = new int [nsuploc];
	if (!order) MemoryFail (funcname);

	int nblks, *blks, *blk2cpu;

	asymm.PartBinaryTreeND (tree, order, _param,
									nblks, blks, blk2cpu);

//	for (i=0;i<nsuploc;i++) order[i] = i;
//	for (i=0;i<nsuploc;i++) order[i] = nsuploc-1-i;

//	OutArr (_fout," Order = ",nsuploc,order);
	tree.SetComm (*((CMPIComm *)_comm));

// Create inverse order

	int *iorder;

	iorder = new int [nsuploc];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nsuploc;i++) iorder[order[i]] = i;

// Create ordered sprnds

	int *sprndso;

	sprndso = new int [nsuploc+1];
	if (!sprndso) MemoryFail (funcname);

	int iold;

	sprndso[0] = 0;
	for (i=0;i<nsuploc;i++) {
		iold = iorder[i];
		sprndso[i+1] = sprndso[i]+(sprndsloc[iold+1]-sprndsloc[iold]);
	};

//	OutArr (_fout," Spndso = ",nsuploc+1,sprndso);

// Create the corresponding part of the matrix

	int nlistloc = 0;

	for (i=0;i<nblks;i++) {
		if (blk2cpu[i] == myid) nlistloc += (blks[i+1]-blks[i]);
	};

	OutArr (_fout," blks = ",nblks+1,blks);
	OutArr (_fout," blk2cpu = ",nblks,blk2cpu);

	int *listloc;

	listloc = new int [nlistloc];
	if (!listloc) MemoryFail (funcname);

	nlistloc = 0;

	int j;

	for (i=0;i<nblks;i++) {
		if (blk2cpu[i] == myid) {
			for (j=blks[i];j<blks[i+1];j++) {
				listloc[nlistloc] = iorder[j];
				nlistloc++;
			};
		};
	};

	_fout << " nlistloc = " << nlistloc << endl;

	int jj, blai, blaj, nzrow, ilist;

	nzjaloc = 0;
	int nzaloc = 0;
	int nzrowmax = 0;
	int blamaxloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		nzrow = 0;
		blai = sprndsloc[i+1]-sprndsloc[i];
		if (blai > blamaxloc) blamaxloc = blai;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
		nzjaloc += (ialoc[i+1]-ialoc[i]);
		nzaloc += nzrow;
	};

	CSMatrixRS aloc(nlistloc, nzjaloc, nzaloc);
	cout << " nlistloc = " << nlistloc << endl;

	aloc.SetM (ntot);
	aloc.SetN (ntot);

	aloc.SetNlist (nlistloc);
	aloc.SetNsupr (nsuploc);
	aloc.SetNsupcBase (nsuploc);

	aloc.SetBlamx (blamaxloc);

	int *psprndr, *psprndc, *pbsa;

	aloc.GetPartition (psprndr, psprndc, pbsa);

	delete [] psprndr;
	delete [] psprndc;

	psprndr = new int [nsuploc+1];
	if (!psprndr) MemoryFail (funcname);
	psprndc = new int [nsuploc+1];
	if (!psprndc) MemoryFail (funcname);

	for (i=0;i<=nsuploc;i++) psprndr[i] = sprndso[i];
	for (i=0;i<=nsuploc;i++) psprndc[i] = sprndso[i];

	aloc.SetSprndr (psprndr);
	aloc.SetSprndc (psprndc);

	int *plist, *pia, *pja;

	aloc.GetSparsity (plist, pia, pja);

	for (i=0;i<nlistloc;i++) plist[i] = order[listloc[i]];

	int *bsaloc, *iord;
	CIntInt *iarr;
	double *arow;

	bsaloc = new int [nsuploc];
	if (!bsaloc) MemoryFail (funcname);
	iord = new int [nsuploc];
	if (!iord) MemoryFail (funcname);
	iarr = new CIntInt [nsuploc];
	if (!iarr) MemoryFail (funcname);
	arow = new double [nzrowmax];
	if (!arow) MemoryFail (funcname);

	int nzj, ind, blaij, ibs, k;

	double *pa;
	pa = aloc.GetA();

	nzjaloc = 0;
	nzaloc = 0;

	pia[0] = 0;

	int jjnew;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = plist[ilist];
		iold = iorder[i];
		nzj = ialoc[iold+1]-ialoc[iold];
		blai = sprndsloc[iold+1]-sprndsloc[iold];

		nzrow = 0;
		for (j=ialoc[iold];j<ialoc[iold+1];j++) {
			jj = jaloc[j];
			jjnew = order[jj];
			iarr[j-ialoc[iold]].intvalue = jjnew;
			iarr[j-ialoc[iold]].int2value = j-ialoc[iold];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			bsaloc[jjnew] = nzrow;
			nzrow += blai*blaj;
		};

// Generate

		(_fblockrow) (_pgener, iold, arow);

// Sort data within block row

		qsort (iarr, nzj, sizeof (CIntInt), CIntInt::CompareIntInt);

		for (j=0;j<nzj;j++) iord[j] = iarr[j].int2value;

		nzrow = 0;

		for (j=ialoc[iold];j<ialoc[iold+1];j++) {
			ind = iord[j-ialoc[iold]];
			jj = jaloc[ialoc[iold]+ind];
			jjnew = order[jj];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			blaij = blai*blaj;
			ibs = bsaloc[jjnew];
			pja[nzjaloc] = jjnew;
			pbsa[nzjaloc] = nzaloc+nzrow;
			nzjaloc++;
			for (k=0;k<blaij;k++) {
				pa[nzaloc+nzrow+k] = arow[ibs+k];
			};
			nzrow += blai*blaj;
		};
		pia[ilist+1] = nzjaloc;
		nzaloc += nzrow;

	};

	aloc.SetNzja (nzjaloc);
	aloc.SetNza (nzaloc);
	aloc.SetNzatot (nzaloc);

//	_fout << " Aloc = " << aloc << endl;

// Store rhs and initial guess

	int nloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		nloc += (sprndsloc[i+1]-sprndsloc[i]);
	};

	CSVector bloc(nloc), solloc(nloc), xloc(nloc), normloc(nloc);

	double *pvect = bloc.GetVect ();

	int nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_frhs) (_pgener, i, pvect+nz);
		nz += blai;
	};

	pvect = xloc.GetVect ();

	nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_finiguess) (_pgener, i, pvect+nz);
		nz += blai;
	};

	normloc.SetSVect (1.0e0);

	solloc = xloc;

// Split matrix into L and U parts

	CSMatrixRS mtralloc, mtrauloc;

	aloc.CombineLU (*((CMPIComm *)_comm), nblks, blks, blk2cpu,
									mtralloc, mtrauloc);
//	_fout << " mtralloc = " << mtralloc << endl;
//	_fout << " mtrauloc = " << mtrauloc << endl;

	CSMatrixRS mtrdummyrs;

	aloc = mtrdummyrs;

// Create global structures

	CGSMatrixRS gmtral, gmtrau;

	gmtral.RS2GRS (false, mtralloc, nblks, blks,
							myid, blk2cpu);
//	_fout << " gmtral = " << gmtral << endl;

	gmtrau.RS2GRS (false, mtrauloc, nblks, blks,
							myid, blk2cpu);
//	_fout << " gmtrau = " << gmtrau << endl;

	mtralloc = mtrdummyrs;
	mtrauloc = mtrdummyrs;

// Solve the system by the solver

	CGSMatrixRS gmtrl, gmtru;

//	_fout << " bloc = " << bloc << endl;

	solloc.Solver (_fout, tree, // Solve unsymmetric system in parallel mode
						gmtral, gmtrau,
						gmtral, gmtrau,
						false, gmtral,
						gmtrl, gmtru,
						bloc, xloc, normloc, 
						_param);

// Exchange the solution

	CSVector sol(ntot);

	sol.SetSVect (0.0e0);

	double *psol = sol.GetVect ();

	pvect = solloc.GetVect ();

	nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		ibs = sprndsloc[i];
		for (j=0;j<blai;j++) psol[ibs+j] = pvect[nz+j];
		nz += blai;
	};

	CMPIExchange::ExchangeArrayMPI (*((CMPIComm *)_comm),
												DOUBLEVALUE, ADD, ntot, psol, psol);

// Store the solution via calls to an unterface funtion

	for (i=0;i<nsuploc;i++) {
		iold = iorder[i];
		blai = sprndsloc[iold+1]-sprndsloc[iold];
		ibs = sprndsloc[iold];
		(_fstoresol) (_pgener, iold, psol+ibs);
	};

// Free work arrays

	delete [] sprndsloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] order;
	delete [] blks;
	delete [] blk2cpu;
	delete [] iorder;
	delete [] sprndso;
	delete [] listloc;
	delete [] bsaloc;
	delete [] iord;
	delete [] iarr;
	delete [] arow;

};

// Author: Kharchenko S.A.
// Solve real block system of equations in parallel mode in abstract interface
//========================================================================================
void AbstractParSolver (void *_pgener, // Solve real block system of equations in parallel mode in abstract interface
									void *_comm, ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROWC _fblockrow, 
									CRHSC _frhs, CINIGUESSC _finiguess, CSTORESOLC _fstoresol) {

	const char *funcname = "AbstractParSolverC";

	dcmplx czero (0.0e0,0.0e0);

	int nproc = ((CMPIComm *)_comm)->GetNproc ();
	int myid = ((CMPIComm *)_comm)->GetMyid ();

// Create the sparsity

	int nsuploc;

	(_fnsup) (_pgener, nsuploc);

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	(_fsprnds) (_pgener, sprndsloc);

	int ntot = sprndsloc[nsuploc];

	int *ialoc, *jaloc;

	(_fiaja) (_pgener, ialoc, jaloc);

	int nzjaloc = ialoc[nsuploc];

// Pack the sparsity into the matrix

	CSMatrix asploc (nsuploc,nzjaloc);

	int *plistsp = asploc.GetList ();
	int *piasp = asploc.GetIa ();
	int *pjasp = asploc.GetJa ();

	int i;

	for (i=0;i<nsuploc;i++) plistsp[i] = i;
	for (i=0;i<=nsuploc;i++) piasp[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) pjasp[i] = jaloc[i];

	asploc.SetM (nsuploc);
	asploc.SetN (nsuploc);
	asploc.SetNsupc(nsuploc);
	asploc.SetNsupr (nsuploc);
	asploc.SetNlist (nsuploc);

	CSMatrix asymm = asploc.SymmMtr ();

//	_fout << " Asymm = " << asymm << endl;

	CSMatrix mtrdummy;

	asploc = mtrdummy;

// Create the binary tree's

	int nprocext = _param.ncpuext;
	if (nprocext < nproc) nprocext = nproc;

	int nprocloc = nprocext;
	int nchilds = 2;

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nprocloc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocloc+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

	CTree treeext (nprocloc, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	treeext.SortNodes ();

	treeext.CpuListsForNodes ();

	nprocloc = nproc;
	nchilds = 2;

	memory2cpu = new double [nprocloc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocloc+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

	CTree tree (nprocloc, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	tree.SortNodes ();

	tree.CpuListsForNodes ();

// Create the set of communicators

	CMPIComm commloc;

	for (i=0;i<nproc;i++) {
		CMPIComm comm;
		comm = ((CMPIComm *)_comm)->CreateComm (1, &i);
		if (i==myid) commloc = comm;
	};

// Compute new ordering

	int *order;

	order = new int [nsuploc];
	if (!order) MemoryFail (funcname);

	int nblks, *blks, *blk2cpu;

	if (false) {

// Store new ordering and partitioning

		asymm.PartBinaryTreeND (treeext, order, _param,
										nblks, blks, blk2cpu);

		CSMatrix *pamtr = treeext.GetAblktree ();

		treeext.IncludeBinaryTreeSchur (tree, nblks, blks, *pamtr);

		tree.Block2CpuSchur (blk2cpu);

	} else {

// Find node to set of nodes correspondense

		int nprocpart = 1;

		int nnodesloc = tree.GetNnodes ();
		int nlevloc = tree.GetNlev ();
		CNode *pnodesloc = tree.GetNode ();
		int nnodeslocext = treeext.GetNnodes ();
		CNode *pnodeslocext = treeext.GetNode ();

		int *nd2ndext;
		int *jnd2ndext;

		nd2ndext = new int [nnodesloc+1];
		if (!nd2ndext) MemoryFail (funcname);
		jnd2ndext = new int [nnodeslocext];
		if (!jnd2ndext) MemoryFail (funcname);

		for (i=0;i<=nnodesloc;i++) nd2ndext[i] = 0;

		int iploc = nnodeslocext-1;

		int nratio = nprocext/nproc;
		int nndlow = nratio*2-1;

		int ilev, icount;

		for (i=nnodesloc-1;i>=0;i--) {
			ilev = pnodesloc[i].GetNodelv ();
			if (ilev < nlevloc-1) {
				if (pnodeslocext[iploc].GetNodelv () == ilev) {
					nd2ndext[i+1]++;
					iploc--;
				} else {
					throw " ParFBSS_Array::ParFBSS_Solve_Distrib: node not found ";
				};
			} else {
				icount = 0;
				while (iploc >= 0) {
					if (pnodeslocext[iploc].GetNodelv () >= ilev && icount < nndlow) {
						nd2ndext[i+1]++;
						iploc--;
						icount++;
					} else {
						break;
					};
				};
			};
		};

		for (i=0;i<nnodesloc;i++) nd2ndext[i+1] = nd2ndext[i]+nd2ndext[i+1];

		int *iptr;

		iptr = new int [nnodesloc];
		if (!iptr) MemoryFail (funcname);

		for (i=0;i<nnodesloc;i++) iptr[i] = nd2ndext[i];

		iploc = nnodeslocext-1;

		int k, jj;

		for (i=nnodesloc-1;i>=0;i--) {
			ilev = pnodesloc[i].GetNodelv ();
			if (ilev < nlevloc-1) {
				if (pnodeslocext[iploc].GetNodelv () == ilev) {
					k = iptr[i];
					jnd2ndext[k] = iploc;
					iptr[i]++;
					iploc--;
				} else {
					throw " ParFBSS_Array::ParFBSS_Solve_Distrib: node not found ";
				};
			} else {
				icount = 0;
				while (iploc >= 0) {
					if (pnodeslocext[iploc].GetNodelv () >= ilev && icount < nndlow) {
						k = iptr[i];
						jnd2ndext[k] = iploc;
						iptr[i]++;
						iploc--;
						icount++;
					} else {
						break;
					};
				};
			};
		};

		for (i=0;i<nnodesloc;i++) {
			std::sort (jnd2ndext+nd2ndext[i], jnd2ndext+nd2ndext[i+1]);
		};

		delete [] iptr;

// Distribute point data

		int nsuploc = asymm.GetNlist ();

		int *sprndsloc;
		int *sp2cpuloc;
		int *orderloc;
		int *orderlocnew;
		int *ordersp;

		sprndsloc = new int [nsuploc+1];
		if (!sprndsloc) MemoryFail (funcname);
		sp2cpuloc = new int [nsuploc];
		if (!sp2cpuloc) MemoryFail (funcname);
		orderloc = new int [nsuploc];
		if (!orderloc) MemoryFail (funcname);
		orderlocnew = new int [nsuploc];
		if (!orderlocnew) MemoryFail (funcname);
		ordersp = new int [nsuploc];
		if (!ordersp) MemoryFail (funcname);

		for (i=0;i<=nsuploc;i++) sprndsloc[i] = i;

		asymm.SplitBlockMatrixRecursive (commloc, nprocext,
														nsuploc, sprndsloc, 
														sp2cpuloc,
														orderloc, ordersp);

		delete [] sprndsloc;
		delete [] ordersp;

		for (i=0;i<nsuploc;i++) order[i] = orderloc[i];

		nblks = nprocext;

		blks = new int [nblks+1];
		if (!blks) MemoryFail (funcname);
		blk2cpu = new int [nblks];
		if (!blk2cpu) MemoryFail (funcname);

		for (i=0;i<=nblks;i++) blks[i] = 0;

		int icpu;

		for (i=0;i<nsuploc;i++) {
			icpu = sp2cpuloc[i];
			blks[icpu+1]++;
		};

		CMPIExchange::ExchangeArrayMPI (commloc, INTEGERVALUE, ADD, nblks+1, blks, blks);

		for (i=0;i<nblks;i++) blks[i+1] = blks[i]+blks[i+1];

		for (i=0;i<nblks;i++) blk2cpu[i] = i;

		delete [] sp2cpuloc;

// Reorder the matrix in-place

		asymm.ChangeIndices (commloc, orderloc);

// Exchange matrix data according to the blocks partitioning

		int myidloc = 0;

		int niblk = nprocext / nprocpart;

		int iblkbeg = myidloc*niblk;
		int iblkend = (myidloc+1)*niblk-1;

		int ibegloc = blks[iblkbeg];
		int iendloc = blks[iblkend+1]-1;

		int nploc = 0;
		if (iendloc >= ibegloc) nploc = 1;

		CSMatrix mtrnew1;

		mtrnew1 = asymm.GetSubmatrixViaIntervals (commloc,
																nploc, &ibegloc, &iendloc);

		asymm = mtrnew1;
		mtrnew1 = mtrdummy;

		if (_param.decomptype > 0) {

// Create local arrays

			int nhcellsloc = nblks;
			int *hcellsloc = blks;
			int *hcells2cpuloc = blk2cpu;

// Call the new ordering and partitioning routine

			int *blk2cpu1;

			blk2cpu1 = new int [nprocext];
			if (!blk2cpu1) MemoryFail (funcname);

			int j;

			for (i=0;i<nprocpart;i++) {
				for (j=i*niblk;j<(i+1)*niblk;j++) {
					blk2cpu1[j] = i;
				};
			};

			int nlist1loc = asymm.GetNlist ();

			int *order1loc;

			order1loc = new int [nlist1loc];
			if (!order1loc) MemoryFail (funcname);

			int *blk2cpunew;

			int nhcells_ext;
			int *blknew2ext;

			asymm.DecompSchurBndRecursive (commloc, _param.decomptype-1, _param.ncycle, nprocext,
														nhcellsloc, hcellsloc, blk2cpu1, hcells2cpuloc,
														order1loc,
														nhcells_ext, blknew2ext, blks, blk2cpunew, blk2cpu);

			nblks = blknew2ext[nhcells_ext];

// Store collapsed blocks partitioning in the tree

			int *blkstree;

			blkstree = new int [nnodesloc+1];
			if (!blkstree) MemoryFail (funcname);

			for (i=0;i<=nnodesloc;i++) blkstree[i] = 0;

			for (i=0;i<nnodesloc;i++) {
				for (j=nd2ndext[i];j<nd2ndext[i+1];j++) {
					jj = jnd2ndext[j];
					blkstree[i+1] += blknew2ext[jj+1]-blknew2ext[jj];
				};
			};

			for (i=0;i<nnodesloc;i++) blkstree[i+1] = blkstree[i]+blkstree[i+1];

			for (i=0;i<nnodesloc;i++) {
				int indbegloc = blkstree[i];
				int indendloc = blkstree[i+1]-1;
				pnodesloc[i].SetIndbeg (indbegloc);
				pnodesloc[i].SetIndend (indendloc);
				pnodesloc[i].SetIndbegtot (indbegloc);
				pnodesloc[i].SetIndendtot (indendloc);
			};

			delete [] blkstree;

// Reorder the matrix once again

			asymm.ChangeIndices (commloc, order1loc);

			int *ibeg1arr;
			int *iend1arr;

			ibeg1arr = new int [nblks];
			if (!ibeg1arr) MemoryFail (funcname);
			iend1arr = new int [nblks];
			if (!iend1arr) MemoryFail (funcname);

			int np1 = 0;

			for (i=0;i<nblks;i++) {
				if (blk2cpunew[i] == myidloc && blks[i] != blks[i+1]) {
					ibeg1arr[np1] = blks[i];
					iend1arr[np1] = blks[i+1]-1;
					np1++;
				};
			};

			delete [] blk2cpunew;

			mtrnew1 = asymm.GetSubmatrixViaIntervals (commloc,
																		np1, ibeg1arr, iend1arr);

			asymm = mtrnew1;
			mtrnew1 = mtrdummy;

// Reassign the ordering array

// Compute the direct ordering for the local partitioning

			int nlistord = 0;

			for (i=0;i<nhcellsloc;i++) {
				if (blk2cpu1[i] == myidloc) {
					nlistord += (hcellsloc[i+1]-hcellsloc[i]);
				};
			};

			int *listord;

			listord = new int [nlistord];
			if (!listord) MemoryFail (funcname);

			nlistord = 0;

			for (i=0;i<nhcellsloc;i++) {
				if (blk2cpu1[i] == myidloc) {
					for (j=hcellsloc[i];j<hcellsloc[i+1];j++) {
						listord[nlistord] = j;
						nlistord++;
					};
				};
			};

			int nlistpart;
			int *orderpart;

			CSMatrix::DirOrderViaIntervals (commloc,
														nlistord, listord, order1loc,
														np1, ibeg1arr, iend1arr, 
														nlistpart, orderpart);

			CSMatrix::OrderViaIntervals (commloc,
													nsuploc, orderloc,
													np1, ibeg1arr, iend1arr, 
													orderpart,
													order);

			delete [] hcellsloc;
			delete [] hcells2cpuloc;
			delete [] blk2cpu1;

			delete [] ibeg1arr;
			delete [] iend1arr;

			delete [] orderpart;
			delete [] order1loc;

			delete [] listord;
			delete [] blknew2ext;

		};

		delete [] orderlocnew;
		delete [] orderloc;

// Compute block sparsity of a matrix

		CSMatrix ablkstr;

		if (_param.decomptype <= 0) {
			ablkstr = tree.MaximalBlockSparsityForSortedTree (nblks);
		} else {
			ablkstr = asymm.ComputeBlockSparsity (commloc, nblks, blks);
		};

		if (myid == 0) ablkstr.A2Ps ("ABlkStr.ps",0,&i);
		if (_param.collap > 0) {
			if (myid == 0) asymm.A2Ps (_param.collap,"AMatrStr.ps",nblks,blks);
		};

		CSMatrix *pablkstr = tree.GetAblktree ();

		*pablkstr = ablkstr;

// Condense blocks rows distribution

		int niblkloc = nprocext / nproc;

		int *cpuext2cpu;

		cpuext2cpu = new int [nprocext];
		if (!cpuext2cpu) MemoryFail (funcname);

		int j;

		for (i=0;i<nproc;i++) {
			for (j=0;j<niblkloc;j++) {
				cpuext2cpu[i*niblkloc+j] = i;
			};
		};

		int icpunew;

		for (i=0;i<nblks;i++) {
			icpu = blk2cpu[i];
			icpunew = cpuext2cpu[icpu];
			blk2cpu[i] = icpunew;
		};

		delete [] cpuext2cpu;

		tree.StorePartitioning (nblks, blks, blk2cpu);

//		std::ofstream fffout ("ChkDecomp.dat");
//		OutArr (fffout," Hcells ",nhcells+1,hcells);
//		OutArr (fffout," Hcells2cpu ",nhcells,hcells2cpu);
//		fffout << " Tree = " << treeloc << std::endl;

//		if (myid > 0) delete [] hcells;
//		if (myid > 0) delete [] hcells2cpu;

		delete [] nd2ndext;
		delete [] jnd2ndext;

	};

// Set communicator

	tree.SetComm (*((CMPIComm *)_comm));

// Create inverse order

	int *iorder;

	iorder = new int [nsuploc];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nsuploc;i++) iorder[order[i]] = i;

// Create ordered sprnds

	int *sprndso;

	sprndso = new int [nsuploc+1];
	if (!sprndso) MemoryFail (funcname);

	int iold;

	sprndso[0] = 0;
	for (i=0;i<nsuploc;i++) {
		iold = iorder[i];
		sprndso[i+1] = sprndso[i]+(sprndsloc[iold+1]-sprndsloc[iold]);
	};

//	OutArr (_fout," Spndso = ",nsuploc+1,sprndso);

// Create the corresponding part of the matrix

	int nlistloc = 0;

	for (i=0;i<nblks;i++) {
		if (blk2cpu[i] == myid) nlistloc += (blks[i+1]-blks[i]);
	};

	OutArr (_fout," blks = ",nblks+1,blks);
	OutArr (_fout," blk2cpu = ",nblks,blk2cpu);

	int *listloc;

	listloc = new int [nlistloc];
	if (!listloc) MemoryFail (funcname);

	nlistloc = 0;

	int j;

	for (i=0;i<nblks;i++) {
		if (blk2cpu[i] == myid) {
			for (j=blks[i];j<blks[i+1];j++) {
				listloc[nlistloc] = iorder[j];
				nlistloc++;
			};
		};
	};

	_fout << " nlistloc = " << nlistloc << endl;

	int jj, blai, blaj, nzrow, ilist;

	nzjaloc = 0;
	int nzaloc = 0;
	int nzrowmax = 0;
	int blamaxloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		nzrow = 0;
		blai = sprndsloc[i+1]-sprndsloc[i];
		if (blai > blamaxloc) blamaxloc = blai;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
		nzjaloc += (ialoc[i+1]-ialoc[i]);
		nzaloc += nzrow;
	};

	CSMatrixCS aloc(nlistloc, nzjaloc, nzaloc);
	cout << " nlistloc = " << nlistloc << endl;

	aloc.SetM (ntot);
	aloc.SetN (ntot);

	aloc.SetNlist (nlistloc);
	aloc.SetNsupr (nsuploc);
	aloc.SetNsupcBase (nsuploc);

	aloc.SetBlamx (blamaxloc);

	int *psprndr, *psprndc, *pbsa;

	aloc.GetPartition (psprndr, psprndc, pbsa);

	delete [] psprndr;
	delete [] psprndc;

	psprndr = new int [nsuploc+1];
	if (!psprndr) MemoryFail (funcname);
	psprndc = new int [nsuploc+1];
	if (!psprndc) MemoryFail (funcname);

	for (i=0;i<=nsuploc;i++) psprndr[i] = sprndso[i];
	for (i=0;i<=nsuploc;i++) psprndc[i] = sprndso[i];

	aloc.SetSprndr (psprndr);
	aloc.SetSprndc (psprndc);

	int *plist, *pia, *pja;

	aloc.GetSparsity (plist, pia, pja);

	for (i=0;i<nlistloc;i++) plist[i] = order[listloc[i]];

	int *bsaloc, *iord;
	CIntInt *iarr;
	dcmplx *arow;

	bsaloc = new int [nsuploc];
	if (!bsaloc) MemoryFail (funcname);
	iord = new int [nsuploc];
	if (!iord) MemoryFail (funcname);
	iarr = new CIntInt [nsuploc];
	if (!iarr) MemoryFail (funcname);
	arow = new dcmplx [nzrowmax];
	if (!arow) MemoryFail (funcname);

	int nzj, ind, blaij, ibs, k;

	dcmplx *pa;
	pa = aloc.GetA();

	nzjaloc = 0;
	nzaloc = 0;

	pia[0] = 0;

	int jjnew;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = plist[ilist];
		iold = iorder[i];
		nzj = ialoc[iold+1]-ialoc[iold];
		blai = sprndsloc[iold+1]-sprndsloc[iold];

		nzrow = 0;
		for (j=ialoc[iold];j<ialoc[iold+1];j++) {
			jj = jaloc[j];
			jjnew = order[jj];
			iarr[j-ialoc[iold]].intvalue = jjnew;
			iarr[j-ialoc[iold]].int2value = j-ialoc[iold];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			bsaloc[jjnew] = nzrow;
			nzrow += blai*blaj;
		};

// Generate

		(_fblockrow) (_pgener, iold, arow);

// Sort data within block row

		qsort (iarr, nzj, sizeof (CIntInt), CIntInt::CompareIntInt);

		for (j=0;j<nzj;j++) iord[j] = iarr[j].int2value;

		nzrow = 0;

		for (j=ialoc[iold];j<ialoc[iold+1];j++) {
			ind = iord[j-ialoc[iold]];
			jj = jaloc[ialoc[iold]+ind];
			jjnew = order[jj];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			blaij = blai*blaj;
			ibs = bsaloc[jjnew];
			pja[nzjaloc] = jjnew;
			pbsa[nzjaloc] = nzaloc+nzrow;
			nzjaloc++;
			for (k=0;k<blaij;k++) {
				pa[nzaloc+nzrow+k] = arow[ibs+k];
			};
			nzrow += blai*blaj;
		};
		pia[ilist+1] = nzjaloc;
		nzaloc += nzrow;

	};

	aloc.SetNzja (nzjaloc);
	aloc.SetNza (nzaloc);
	aloc.SetNzatot (nzaloc);

//	_fout << " Aloc = " << aloc << endl;

// Store rhs and initial guess

	int nloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		nloc += (sprndsloc[i+1]-sprndsloc[i]);
	};

	CSVectorC bloc(nloc), solloc(nloc), xloc(nloc);
	CSVector normloc(nloc);

	dcmplx *pvect = bloc.GetVect ();

	int nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_frhs) (_pgener, i, pvect+nz);
		nz += blai;
	};

	pvect = xloc.GetVect ();

	nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_finiguess) (_pgener, i, pvect+nz);
		nz += blai;
	};

	normloc.SetSVect (1.0e0);

	solloc = xloc;

// Split matrix into L and U parts

	CSMatrixCS mtralloc, mtrauloc;

	aloc.CombineLU (*((CMPIComm *)_comm), nblks, blks, blk2cpu,
									mtralloc, mtrauloc);
//	_fout << " mtralloc = " << mtralloc << endl;
//	_fout << " mtrauloc = " << mtrauloc << endl;

	CSMatrixCS mtrdummycs;

	aloc = mtrdummycs;

// Create global structures

	CGSMatrixCS gmtral, gmtrau;

	gmtral.CS2GCS (false, mtralloc, nblks, blks,
							myid, blk2cpu);
//	_fout << " gmtral = " << gmtral << endl;

	gmtrau.CS2GCS (false, mtrauloc, nblks, blks,
							myid, blk2cpu);
//	_fout << " gmtrau = " << gmtrau << endl;

	mtralloc = mtrdummycs;
	mtrauloc = mtrdummycs;

// Solve the system by the solver

	CGSMatrixCS gmtrl, gmtru;

//	_fout << " bloc = " << bloc << endl;

	int job = 0;

	solloc.Solver (_fout, tree, // Solve unsymmetric system in parallel mode
						gmtral, gmtrau,
						job, gmtrl, gmtru,
						bloc, xloc,
						_param);

// Exchange the solution

	CSVectorC sol(ntot);

	sol.SetSVect (czero);

	dcmplx *psol = sol.GetVect ();

	pvect = solloc.GetVect ();

	nz = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		i = listloc[ilist];
		blai = sprndsloc[i+1]-sprndsloc[i];
		ibs = sprndsloc[i];
		for (j=0;j<blai;j++) psol[ibs+j] = pvect[nz+j];
		nz += blai;
	};

	CMPIExchange::ExchangeArrayMPI (*((CMPIComm *)_comm),
												DOUBLEVALUE, ADD, ntot*2, psol, psol);

// Store the solution via calls to an unterface funtion

	for (i=0;i<nsuploc;i++) {
		iold = iorder[i];
		blai = sprndsloc[iold+1]-sprndsloc[iold];
		ibs = sprndsloc[iold];
		(_fstoresol) (_pgener, iold, psol+ibs);
	};

// Free work arrays

	delete [] sprndsloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] order;
	delete [] blks;
	delete [] blk2cpu;
	delete [] iorder;
	delete [] sprndso;
	delete [] listloc;
	delete [] bsaloc;
	delete [] iord;
	delete [] iarr;
	delete [] arow;

};

// Author: Kharchenko S.A.
// Solve complex block system of equations in abstract interface
//========================================================================================
void AbstractSolver (void *_pgener, // Solve complex block system of equations in abstract interface
									ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROWC _fblockrow, 
									CRHSC _frhs, CINIGUESSC _finiguess, CSTORESOLC _fstoresol) {

	const char *funcname = "AbstractSolver";

	dcmplx czero (0.0e0,0.0e0);

// Create the matrix

	int nsuploc;

	(_fnsup) (_pgener, nsuploc);

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	(_fsprnds) (_pgener, sprndsloc);

	int nloc = sprndsloc[nsuploc];

	int *ialoc, *jaloc;

	(_fiaja) (_pgener, ialoc, jaloc);

	int nzjaloc = ialoc[nsuploc];

	int i, j, jj, blai, blaj, nzrow;

	int nzaloc = 0;
	int nzrowmax = 0;
	int blamaxloc = 0;

	for (i=0;i<nsuploc;i++) {
		nzrow = 0;
		blai = sprndsloc[i+1]-sprndsloc[i];
		if (blai > blamaxloc) blamaxloc = blai;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
		nzaloc += nzrow;
	};

	CSMatrixCS aloc(nsuploc, nzjaloc, nzaloc);

	aloc.SetM (nloc);
	aloc.SetN (nloc);

	aloc.SetNlist (nsuploc);
	aloc.SetNsupr (nsuploc);
	aloc.SetNsupcBase (nsuploc);

	aloc.SetBlamx (blamaxloc);

	int *psprndr, *psprndc, *pbsa;

	aloc.GetPartition (psprndr, psprndc, pbsa);

	for (i=0;i<=nsuploc;i++) psprndr[i] = sprndsloc[i];
	for (i=0;i<=nsuploc;i++) psprndc[i] = sprndsloc[i];

	int *plist, *pia, *pja;

	aloc.GetSparsity (plist, pia, pja);

	for (i=0;i<nsuploc;i++) plist[i] = i;
	for (i=0;i<=nsuploc;i++) pia[i] = ialoc[i];

	int *bsaloc, *iord;
	CIntInt *iarr;
	dcmplx *arow;

	bsaloc = new int [nsuploc];
	if (!bsaloc) MemoryFail (funcname);
	iord = new int [nsuploc];
	if (!iord) MemoryFail (funcname);
	iarr = new CIntInt [nsuploc];
	if (!iarr) MemoryFail (funcname);
	arow = new dcmplx [nzrowmax];
	if (!arow) MemoryFail (funcname);

	int nzj, ind, blaij, ibs, k;

	dcmplx *pa;
	pa = aloc.GetA();

	nzaloc = 0;

	for (i=0;i<nsuploc;i++) {

		nzj = ialoc[i+1]-ialoc[i];
		blai = sprndsloc[i+1]-sprndsloc[i];

		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			iarr[j-ialoc[i]].intvalue = jj;
			iarr[j-ialoc[i]].int2value = j-ialoc[i];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			bsaloc[jj] = nzrow;
			nzrow += blai*blaj;
		};

// Generate

		(_fblockrow) (_pgener, i, arow);

// Sort data within block row

		qsort (iarr, nzj, sizeof (CIntInt), CIntInt::CompareIntInt);

//		for (j=0;j<nzj;j++) iord[iarr[j].int2value] = j;
		for (j=0;j<nzj;j++) iord[j] = iarr[j].int2value;

		nzrow = 0;

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			ind = iord[j-ialoc[i]];
			jj = jaloc[ialoc[i]+ind];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			blaij = blai*blaj;
			ibs = bsaloc[jj];
			pja[j] = jj;
			for (k=0;k<blaij;k++) {
				pa[nzaloc+nzrow+k] = arow[ibs+k];
			};
			nzrow += blai*blaj;
		};

// Store bsa array

		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = pja[j];
			blaj = sprndsloc[jj+1]-sprndsloc[jj];
			pbsa[j] = nzaloc;
			nzaloc += blai*blaj;
		};
	};

	aloc.SetNzja (nzjaloc);
	aloc.SetNza (nzaloc);
	aloc.SetNzatot (nzaloc);

// Store rhs and initial guess

	CSVectorC b(nloc), sol(nloc), x(nloc);
	CSVector norm(nloc);

	dcmplx *pvect = b.GetVect ();

	int nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_frhs) (_pgener, i, pvect+nz);
		nz += blai;
	};

	pvect = x.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_finiguess) (_pgener, i, pvect+nz);
		nz += blai;
	};

	norm.SetSVect (1.0e0);

	sol = x;

// Solve the system by the solver

	CSMatrixCS mtrl, mtru;
	int *order1loc;

	order1loc = new int [1];

	sol.Solver (_fout, aloc,
					0, order1loc, mtrl, mtru,
						b, 
						_param);

	delete [] order1loc;

// Store the solution via interface function

	pvect = sol.GetVect ();

	nz = 0;

	for (i=0;i<nsuploc;i++) {
		blai = sprndsloc[i+1]-sprndsloc[i];
		(_fstoresol) (_pgener, i, pvect+nz);
		nz += blai;
	};

// Free work arrays

	delete [] sprndsloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;
	delete [] iord;
	delete [] iarr;
	delete [] arow;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
