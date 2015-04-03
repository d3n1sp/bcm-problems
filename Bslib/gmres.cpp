//------------------------------------------------------------------------------------------------
// File: gmres.cpp
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

#include "tree.h"
#include "slvparam.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "qrd.h"
#include "mvm.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// Description: Compute block Givens rotation
// BlockGivens()
//========================================================================================
void BlockGivens (int _ngiv, // Compute block Givens rotation
					double *_giv, double *_taugiv, 
					double *_aloc, double *_givarr) {

//	const char *funcname = "BlockGivens";

	double dzero = 0.0e0;
	double done = 1.0e0;

	int nloc = _ngiv*2;

// Compute QR decomposition

	QrdBlk (nloc, nloc, _giv, nloc, _taugiv);

// Compute Q factor explicitely

	int kii;

	for (kii=0;kii<nloc*nloc;kii++) _aloc[kii] = dzero;
	for (kii=0;kii<nloc;kii++) _aloc[kii*nloc+kii] = done;

	MvmQBlk (nloc, nloc, nloc,
				_giv, nloc, _taugiv,
				_aloc, nloc, _givarr, nloc);

};

// Author: Kharchenko S.A.
// Description: Apply block Givens rotation
// ApplyBlockGivens()
//========================================================================================
void ApplyBlockGivens (int _ngiv, // Apply block Givens rotation
						int _indgiv, double *_giv, 
						int _nrhs, double *_rhs, int _ldrhs,
						double *_aloc) {

//	const char *funcname = "ApplyBlockGivens";

	double dzero = 0.0e0;

	int nloc = _ngiv*2;

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs*nloc;kii++) _aloc[kii] = dzero;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			double daux = dzero;
			for (kkk=0;kkk<nloc;kkk++) {
//				caux += conj(_giv[kjj*nloc+kkk]) * _rhs[kii*_ldrhs+_indgiv*_ngiv+kkk];
				daux += (_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * _giv[kjj*nloc+kkk]);
//				caux.PlEq (_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * conj(_giv[kjj*nloc+kkk]));
			};
			_aloc[kii*nloc+kjj] = daux;
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			_rhs[kii*_ldrhs+_indgiv*_ngiv+kjj] = _aloc[kii*nloc+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Apply block hermitian transposed Givens rotation
// ApplyBlockGivensH()
//========================================================================================
void ApplyBlockGivensH (int _ngiv, // Apply block hermitian transposed Givens rotation
						int _indgiv, double *_giv, 
						int _nrhs, double *_rhs, int _ldrhs,
						double *_aloc) {

//	const char *funcname = "ApplyBlockGivensH";

	double dzero = 0.0e0;

	int nloc = _ngiv*2;

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs*nloc;kii++) _aloc[kii] = dzero;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			double daux = dzero;
			for (kkk=0;kkk<nloc;kkk++) {
//				caux += conj(_giv[kjj*nloc+kkk]) * _rhs[kii*_ldrhs+_indgiv*_ngiv+kkk];
				daux +=_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * _giv[kkk*nloc+kjj];
//				caux.PlEq (_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * _giv[kkk*nloc+kjj]);
			};
			_aloc[kii*nloc+kjj] = daux;
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			_rhs[kii*_ldrhs+_indgiv*_ngiv+kjj] = _aloc[kii*nloc+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Solve system with triangular dense matrix
// SolveTriangular()
//========================================================================================
void SolveTriangular (int _n, int _nrhs, // Solve system with triangular dense matrix
						double *_rmatr, int _ldrmatr,
						double *_rhs, int _ldrhs,
						double *_sol, int _ldsol) {

//	const char *funcname = "SolveTriangular";

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<_n;kjj++) {
			_sol[kii*_ldsol+kjj] = _rhs[kii*_ldrhs+kjj];
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=_n-1;kjj>=0;kjj--) {
			double daux = _sol[kii*_ldsol+kjj];
			double daux1 = daux / _rmatr[kjj*_ldrmatr+kjj];
			_sol[kii*_ldsol+kjj] = daux1;
			for (kkk=0;kkk<kjj;kkk++) {
				daux = _rmatr[kjj*_ldrmatr+kkk] * daux1;
				_sol[kii*_ldsol+kkk] -= daux;
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Compute identity block Givens rotation
// IdentityBlockGivens()
//========================================================================================
void IdentityBlockGivens (int _ngiv, // Compute block Givens rotation
									dcmplx *_givarr) {

//	const char *funcname = "BlockGivens";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nloc = _ngiv*2;

	int i;

	for (i=0;i<nloc*nloc;i++) _givarr[i] = czero;
	for (i=0;i<nloc;i++) _givarr[i*nloc+i] = cone;

};

// Author: Kharchenko S.A.
// Description: Compute block Givens rotation
// BlockGivens()
//========================================================================================
void BlockGivens (int _ngiv, // Compute block Givens rotation
					dcmplx *_giv, dcmplx *_taugiv, 
					dcmplx *_aloc, dcmplx *_givarr) {

//	const char *funcname = "BlockGivens";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nloc = _ngiv*2;

// Compute QR decomposition

	QrdBlk (nloc, nloc, _giv, nloc, _taugiv);

// Compute Q factor explicitely

	int kii;

	for (kii=0;kii<nloc*nloc;kii++) _aloc[kii] = czero;
	for (kii=0;kii<nloc;kii++) _aloc[kii*nloc+kii] = cone;

	MvmQBlk (nloc, nloc, nloc,
				_giv, nloc, _taugiv,
				_aloc, nloc, _givarr, nloc);

};

// Author: Kharchenko S.A.
// Description: Apply block Givens rotation
// ApplyBlockGivens()
//========================================================================================
void ApplyBlockGivens (int _ngiv, // Apply block Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc) {

//	const char *funcname = "ApplyBlockGivens";

	dcmplx czero (0.0e0,0.0e0);

	int nloc = _ngiv*2;

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs*nloc;kii++) _aloc[kii] = czero;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<nloc;kkk++) {
//				caux += conj(_giv[kjj*nloc+kkk]) * _rhs[kii*_ldrhs+_indgiv*_ngiv+kkk];
				caux.PlEq (_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * conj(_giv[kjj*nloc+kkk]));
			};
			_aloc[kii*nloc+kjj] = caux;
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			_rhs[kii*_ldrhs+_indgiv*_ngiv+kjj] = _aloc[kii*nloc+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Apply block hermitian transposed Givens rotation
// ApplyBlockGivensH()
//========================================================================================
void ApplyBlockGivensH (int _ngiv, // Apply block hermitian transposed Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc) {

//	const char *funcname = "ApplyBlockGivensH";

	dcmplx czero (0.0e0,0.0e0);

	int nloc = _ngiv*2;

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs*nloc;kii++) _aloc[kii] = czero;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<nloc;kkk++) {
//				caux += conj(_giv[kjj*nloc+kkk]) * _rhs[kii*_ldrhs+_indgiv*_ngiv+kkk];
				caux.PlEq (_rhs[kii*_ldrhs+_indgiv*_ngiv+kkk] * _giv[kkk*nloc+kjj]);
			};
			_aloc[kii*nloc+kjj] = caux;
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			_rhs[kii*_ldrhs+_indgiv*_ngiv+kjj] = _aloc[kii*nloc+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Solve system with triangular dense matrix
// SolveTriangular()
//========================================================================================
void SolveTriangular (int _n, int _nrhs, // Solve system with triangular dense matrix
						dcmplx *_rmatr, int _ldrmatr,
						dcmplx *_rhs, int _ldrhs,
						dcmplx *_sol, int _ldsol) {

//	const char *funcname = "SolveTriangular";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk;

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=0;kjj<_n;kjj++) {
			_sol[kii*_ldsol+kjj] = _rhs[kii*_ldrhs+kjj];
		};
	};

	for (kii=0;kii<_nrhs;kii++) {
		for (kjj=_n-1;kjj>=0;kjj--) {
			dcmplx caux = _sol[kii*_ldsol+kjj];
			dcmplx caux1 = caux / _rmatr[kjj*_ldrmatr+kjj];
			_sol[kii*_ldsol+kjj] = caux1;
			for (kkk=0;kkk<kjj;kkk++) {
				caux = _rmatr[kjj*_ldrmatr+kkk] * caux1;
				_sol[kii*_ldsol+kkk] -= caux;
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
// CSVector::BlockGMRES()
//========================================================================================
void CSVector::BlockGMRES (ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
										const CTree &_tree, 
										void *_obja, CMVM _mvm,
										void *_objlu, CSOLVELU _solvelu,
										const CSVector &_x, const CSVector &_b) {

	const char *funcname = "BlockGMRES";

	double dzero = 0.0e0;
	double done = 1.0e0;

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
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
// Determine extended length
//
	int nitmax = _param.niter;

	int nlocext = n;
	if ((nitmax+1)*nrhs > nlocext) nlocext = (nitmax+1)*nrhs;
//
// Allocate temporary vector data
//
	CSVector r(n,npartsloc,nrhs), y(n,npartsloc,nrhs), w(n,npartsloc,nrhs);
//
// Allocate all working arrays
//
	double *p, *yp, *taup;
	double *givarr, *giv, *giv1, *taugiv;
	double *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.ooctyp == 0) {
		p = new double [ldp*(nitmax+1)*nrhs];
		if (!p) MemoryFail (funcname);
	} else {
		p = new double [2*ldp*nrhs];
		if (!p) MemoryFail (funcname);
	};

	yp = new double [ldp*nrhs];
	if (!yp) MemoryFail (funcname);
	taup = new double [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new double [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new double [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new double [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new double [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new double [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new double [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new double [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structure
//
	CQrd qrd (treeclp, nitmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.ooctyp, strbuff, ldp);

	};
//
// Compute initial residual and its norm
//
	CGSMatrixRS gmtrdummy;

	time2 = clock ();

	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	double *pr = r.GetVect ();

	for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = dzero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			p[irhs*nlocext+kii] = pr[irhs*n+kii];
		};
	};

	qrd.UpdateQrd (treeclp, nlocext, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = dzero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii];
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs_2;

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
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi0[irhs] > resimax) resimax = resi0[irhs];
		};
		if (_param.msglev > 1) {
			cout  << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
		};
	};
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = dzero;
//
// Main iterative cycle
//
	int ncolp, iconv;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = dzero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = done;
		};

		double *py = y.GetVect ();

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		time2 = clock ();

		(_solvelu) (_objlu, y, w);
		y = w;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;
		nslvu++;

		time2 = clock ();

//		_gmtra.MvmARect (y,w);
		(_mvm) (_obja, y, w);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		double *pw = w.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = dzero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = dzero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;
//
// Local computations on the root node
//
		if (_tree.myid == _tree.rootcpu) {
//
// Take current block column of H_k
//
			int ibeg=ncolp;
			int iend=ncolp+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = dzero;

			qrd.GetRPartQrd (treeclp, 0, ibeg, iend,
								eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (int igiv=0;igiv<k;igiv++) {

				ApplyBlockGivens (nrhs, 
									igiv, givarr+igiv*nrhs_2*4,
									nrhs, eloc, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = dzero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
							giv, taugiv, giv1, givarr+k*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<(k+1)*nrhs;kii++) {
					rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, bloc, ldbloc, giv);

			ops += 8*nrhs_3;
//
// Check convergence
//
			double resimax = 0.0e0;

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii];
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);
				if (resi[irhs] > resimax) resimax = resi[irhs];

			};

			if (iter%_param.ichk == 0) {

				if (_param.msglev > 1) {
					cout  << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
				};
				if (_param.msglev >= 1) {
					_fout << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
				};
			};

			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (_param.sttype == 0) {
					if (resi[irhs] < _param.eps) iconv++;
				} else {
					if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
				};
			};

		} else {
			iconv = 0;
		};

// Add convergence control numbers on all processors

		int iconvtot = 0;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &iconv, &iconvtot);
		} else {
			iconvtot = iconv;
		};

// Check convergence and exit

		if (iconvtot == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	if (true && _tree.myid == 0) {
		double *hystarr;
		hystarr = new double [ncolp];
		if (!hystarr) MemoryFail (funcname);
		double daux;
		for (kii=0;kii<ncolp;kii++) {
			daux = rmatr[kii*ldrmatr+kii];
			hystarr[kii] = sqrt(daux*daux);
		};
		csvhyst (cout,  "DiagRGmres", ncolp, hystarr);
		csvhyst (_fout, "DiagRGmres", ncolp, hystarr);
		delete [] hystarr;
	};

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = dzero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	double *py = y.GetVect ();

	qrd.MvmQBlk (treeclp, nlocext, ncolp+nrhs, nrhs,
					p, ldp, taup, 
					eloc, ldeloc, 
					yp, nlocext);

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			py[irhs*n+kii] = yp[irhs*nlocext+kii];
		};
	};

	ops += 2*nlocext*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	time2 = clock ();

	(_solvelu) (_objlu, y, w);
	y = w;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;
	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (done,y);

	ops += n*nrhs;

	time2 = clock ();

//	_gmtra.MvmARect (sol,r);
	(_mvm) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;

	pr = r.GetVect ();

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<n;kii++) {
			double auxl = pr[irhs*n+kii];
			aux += auxl*auxl;
		};
		resi[irhs] = aux;
	};

	ops += 2*n*nrhs;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi0);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi0[irhs] = resi[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi0[irhs];
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs] > resimax) resimax = resi[irhs];
		};
		if (_param.msglev > 1) {
			cout  << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
		};
		if (_param.msglev >= 1) {
			_fout << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
		};
	};
//
// Compute SVD of R_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc;
		double *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new double [lworkloc];
		if (!workloc) MemoryFail (funcname);

		dgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine DGESVD" << endl;
			throw " Error in the Lapack routine DGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;

	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrd.CloseFiles ();

	};
//
// Finalize time measurement
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

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

	double dnza = arr1[0];
	ops = arr1[1];

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	if (dnza > 0.0e0) {
		ops = ops / dnza;
	} else {
		ops = 0.0e0;
	};

	if (_tree.myid == _tree.rootcpu) {

		cout  << " GMRES statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;

		_fout << " GMRES statistics: " << endl;
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

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, const CSMatrixCS *_mtra, const CSMatrixCS *_mtrl, const CSMatrixCS *_mtru, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

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

	nza = _mtra->GetNzatot();
	nzlu = _mtrl->GetNzatot();
//
// Init guess to the solution
//
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	p = new dcmplx [ldp*(nitmax+1)*nrhs];
	if (!p) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixCS::MvmA (sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) p[kii] = pr[kii];

	UpdateQrd (n, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, nrhs, 
					p, ldp, bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, ncolp, nrhs,
					p, ldp, taup, 
					eloc, ncolp, 
					py, n);

		ops += 2*n*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		_mtrl->CSMatrixCS::SolveL (y,w);
		_mtru->CSMatrixCS::SolveU (w,y);

		_mtra->CSMatrixCS::MvmA (y,w);

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) p[n*ncolp+kii] = pw[kii];

		UpdateQrd (n, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, ibeg, iend,
						p, ldp, eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, ncolp+nrhs, nrhs,
				p, ldp, taup, 
				eloc, ldeloc, 
				py, n);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_mtrl->CSMatrixCS::SolveL (y,w); 
	_mtru->CSMatrixCS::SolveU (w,y); 

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	_mtra->CSMatrixCS::MvmA (sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] p;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix in the out-of-core mode
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix in the out-of-core mode
							int _nblks, int *_blks, FILE **_mtrafiles, const CSMatrixCS *_mtra, 
							FILE **_mtrlfiles, FILE **_mtrufiles, const CSMatrixCS *_mtrl, const CSMatrixCS *_mtru, 
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param, FILE **_qfiles) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

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

	nza = _mtra->GetNzatot();
	nzlu = _mtrl->GetNzatot();
//
// Init guess to the solution
//
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	int ncolbuff = nrhs*2;

	dcmplx *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;
	dcmplx *qbuff;
	int *bl2file, *bsblk;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);

	qbuff = new dcmplx [ldp*ncolbuff];
	if (!qbuff) MemoryFail (funcname);

	bl2file = new int [nitmax+2+_nblks];
	if (!bl2file) MemoryFail (funcname);
	bsblk = new int [nitmax+2+_nblks];
	if (!bsblk) MemoryFail (funcname);
//
// Init bl2file and bsblk
//
	int ni = (nitmax+1+_nblks-1)/_nblks;

	if (ni<1) ni = 1;

	int nz;

	int ip=0;

	for (int iblk=0;iblk<_nblks-1;iblk++) {
		nz = 0;
		for (int j=0;j<ni;j++) {
			bl2file[ip] = iblk;
			bsblk[ip] = j*ldp*nrhs;
			ip++;
		};
	};
	nz = 0;
	while (ip<nitmax+1) {
		bl2file[ip] = _nblks-1;
		bsblk[ip] = nz;
		ip++;
		nz += ldp*nrhs;
	};
//
// Compute initial residual and its norm
//
	_mtra->CSMatrixCS::MvmA (_nblks, _blks, _mtrafiles, sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	int ifile = bl2file[0];
	int ibs = bsblk[0];

	FPut(_qfiles[ifile],n*nrhs,pr,ibs);

	UpdateQrd (n, 0, nrhs, 
				_qfiles, bl2file, bsblk, ldp, taup, 
				qbuff);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, 1, nrhs, 
					_qfiles, bl2file, bsblk, ldp, bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, k+1, nrhs, nrhs,
					_qfiles, bl2file, bsblk, ldp, taup,
					eloc, ncolp,
					py, n,
					qbuff);

		ops += 2*n*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		_mtrl->CSMatrixCS::SolveL (_nblks, _blks, _mtrlfiles, y,w);
		_mtru->CSMatrixCS::SolveU (_nblks, _blks, _mtrufiles, w,y);

		_mtra->CSMatrixCS::MvmA (_nblks, _blks, _mtrafiles, y,w);

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		ifile = bl2file[k+1];
		ibs = bsblk[k+1];

		FPut (_qfiles[ifile],n*nrhs,pw,ibs);

		UpdateQrd (n, k+1, nrhs, 
					_qfiles, bl2file, bsblk, ldp, taup, 
					qbuff);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
//		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, k+1, k+2, nrhs, 
						_qfiles, bl2file, bsblk, ldp, eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, iter+2, nrhs, nrhs,
				_qfiles, bl2file, bsblk, ldp, taup,
				eloc, ldeloc, 
				py, n,
				qbuff);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_mtrl->CSMatrixCS::SolveL (_nblks, _blks, _mtrlfiles, y,w); 
	_mtru->CSMatrixCS::SolveU (_nblks, _blks, _mtrufiles, w,y); 

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	_mtra->CSMatrixCS::MvmA (_nblks, _blks, _mtrafiles, sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;
	delete [] qbuff;
	delete [] bl2file;
	delete [] bsblk;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							const CGSMatrixCS &_gmtral, const CGSMatrixCS &_gmtrau, 
							const CGSMatrixCS &_gmtrl, const CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int n = _gmtral.GetN();
//
// Init statistics data
//
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
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	p = new dcmplx [ldp*(nitmax+1)*nrhs];
	if (!p) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	gmtrdummy.MvmA (_gmtral, _gmtrau, sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) p[kii] = pr[kii];

	UpdateQrd (n, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, nrhs, 
					p, ldp, bloc, ldbloc);
//	OutArr(_fout," bloc on return from getr ",ldbloc*nrhs,bloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, ncolp, nrhs,
					p, ldp, taup, 
					eloc, ncolp, 
					py, n);

		ops += 2*n*ncolp*nrhs;
//		_fout << " y on return from MvmQBlk" << y << endl;
//
// Multiply by the preconditioned matrix
//
		_gmtrl.SolveL (y,w);
//		_fout << " w on return from SolveL" << w << endl;
		_gmtru.SolveU (w,y);
//		_fout << " y on return from SolveU" << y << endl;

		gmtrdummy.MvmA (_gmtral,_gmtrau,y,w);
//		_fout << " w on return from MvmA" << w << endl;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) p[n*ncolp+kii] = pw[kii];

		UpdateQrd (n, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, ibeg, iend,
						p, ldp, eloc, iend);
//		OutArr(_fout," R aft update = ",iend*nrhs,eloc);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, ncolp+nrhs, nrhs,
				p, ldp, taup, 
				eloc, ldeloc, 
				py, n);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_gmtrl.SolveL (y,w); 
	_gmtru.SolveU (w,y); 

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	gmtrdummy.MvmA (_gmtral,_gmtrau,sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] p;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
										CGSMatrixCS &_gmtra, 
										int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
										const CSVectorC &_x, const CSVectorC &_b, 
										const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int n = _gmtra.GetN();
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = 0;

	int ilist;

	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		nza += _gmtra.mtrarr[iblk].GetNzatot();
	};

	nzlu = 0;
	nzlu += _mtrl.GetNzatot();
	nzlu += _mtru.GetNzatot();

//
// Init guess to the solution
//
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	p = new dcmplx [ldp*(nitmax+1)*nrhs];
	if (!p) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	_gmtra.MvmARect (sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) p[kii] = pr[kii];

	UpdateQrd (n, 0, nrhs, p, ldp, taup);
//	OutArr(_fout," Ploc on return from Update qrd ",ldp*nrhs,p);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, nrhs, 
					p, ldp, bloc, ldbloc);
//	OutArr(_fout," bloc on return from getr ",ldbloc*nrhs,bloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, ncolp, nrhs,
					p, ldp, taup, 
					eloc, ncolp, 
					py, n);

		ops += 2*n*ncolp*nrhs;
//		_fout << " y on return from MvmQBlk" << y << endl;
//
// Multiply by the preconditioned matrix
//
		_mtrl.SolveLUOrd (_order, _mtrl, _mtru, y, w);

		y = w;

		_gmtra.MvmARect (y,w);

//		_fout << " w on return from MvmA" << w << endl;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) p[n*ncolp+kii] = pw[kii];

		UpdateQrd (n, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, ibeg, iend,
						p, ldp, eloc, iend);
//		OutArr(_fout," R aft update = ",iend*nrhs,eloc);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, ncolp+nrhs, nrhs,
				p, ldp, taup, 
				eloc, ldeloc, 
				py, n);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_mtrl.SolveLUOrd (_order, _mtrl, _mtru, y, w);

	y = w;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	_gmtra.MvmARect (sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] p;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
							const CTree &_tree, 
							CGSMatrixCS &_gmtra, 
							int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
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

	int ilist;

	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		nza += _gmtra.mtrarr[iblk].GetNzatot();
	};

	nzlu = 2*_mtrl.GetNzatot();
//
// Init guess to the solution
//
	int npartsloc = _x.GetNparts ();

	CSVectorC sol(n,npartsloc,nrhs);

	sol = _x;
//
// Determine extended length
//
	int nitmax = _param.niter;

	int nlocext = n;
	if ((nitmax+1)*nrhs > nlocext) nlocext = (nitmax+1)*nrhs;
//
// Allocate temporary vector data
//
	CSVectorC r(n,npartsloc,nrhs), y(n,npartsloc,nrhs), w(n,npartsloc,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *yp, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.ooctyp == 0) {
		p = new dcmplx [ldp*(nitmax+1)*nrhs];
		if (!p) MemoryFail (funcname);
	} else {
		p = new dcmplx [2*ldp*nrhs];
		if (!p) MemoryFail (funcname);
	};

	yp = new dcmplx [ldp*nrhs];
	if (!yp) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structure
//
	CQrdC qrd (treeclp, nitmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.ooctyp, strbuff, ldp);

	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	_gmtra.MvmARect (sol,r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			p[irhs*nlocext+kii] = pr[irhs*n+kii];
		};
	};

	qrd.UpdateQrd (treeclp, nlocext, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs_2;

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
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi0[irhs] > resimax) resimax = resi0[irhs];
		};
		cout  << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
		_fout << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
	};
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp, iconv;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		time2 = clock ();

		_mtrl.SolveLUOrd (_order, _mtrl, _mtru, y, w);
		y = w;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;
		nslvu++;

		time2 = clock ();

		_gmtra.MvmARect (y,w);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;
//
// Local computations on the root node
//
		if (_tree.myid == _tree.rootcpu) {
//
// Take current block column of H_k
//
			int ibeg=ncolp;
			int iend=ncolp+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

			qrd.GetRPartQrd (treeclp, 0, ibeg, iend,
								eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (int igiv=0;igiv<k;igiv++) {

				ApplyBlockGivens (nrhs, 
									igiv, givarr+igiv*nrhs_2*4,
									nrhs, eloc, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
							giv, taugiv, giv1, givarr+k*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<(k+1)*nrhs;kii++) {
					rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, bloc, ldbloc, giv);

			ops += 8*nrhs_3;
//
// Check convergence
//
			double resimax = 0.0e0;

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);
				if (resi[irhs] > resimax) resimax = resi[irhs];

			};

			if (iter%_param.ichk == 0) {
				cout  << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
				_fout << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
			};

			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

		} else {
			iconv = 0;
		};

// Add convergence control numbers on all processors

		int iconvtot = 0;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &iconv, &iconvtot);
		} else {
			iconvtot = iconv;
		};

// Check convergence and exit

		if (iconvtot == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	if (true) {
		double *hystarr;
		hystarr = new double [ncolp];
		if (!hystarr) MemoryFail (funcname);
		dcmplx caux;
		for (kii=0;kii<ncolp;kii++) {
			caux = rmatr[kii*ldrmatr+kii];
			hystarr[kii] = sqrt(caux.x*caux.x+caux.y*caux.y);
		};
		csvhyst (cout,  "DiagRGmres", ncolp, hystarr);
		csvhyst (_fout, "DiagRGmres", ncolp, hystarr);
		delete [] hystarr;
	};

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	qrd.MvmQBlk (treeclp, nlocext, ncolp+nrhs, nrhs,
					p, ldp, taup, 
					eloc, ldeloc, 
					yp, nlocext);

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			py[irhs*n+kii] = yp[irhs*nlocext+kii];
		};
	};

	ops += 2*nlocext*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	time2 = clock ();

	_mtrl.SolveLUOrd (_order, _mtrl, _mtru, y, w);
	y = w;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;
	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	time2 = clock ();

	_gmtra.MvmARect (sol,r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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

	ops += 2*n*nrhs;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi0);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi0[irhs] = resi[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi0[irhs];
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs] > resimax) resimax = resi[irhs];
		};
		cout  << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
		_fout << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
	};
//
// Compute SVD of R_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrd.CloseFiles ();

	};
//
// Finalize time measurement
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

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

	double dnza = arr1[0];
	ops = arr1[1];

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " GMRES statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;

		_fout << " GMRES statistics: " << endl;
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

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
										const CTree &_tree, 
										void *_obja, CMVMC _mvmc,
										void *_objlu, CSOLVELUC _solvelu,
										const CSVectorC &_x, const CSVectorC &_b) {

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
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

	CSVectorC sol(n,npartsloc,nrhs);

	sol = _x;
//
// Determine extended length
//
	int nitmax = _param.niter;

	int nlocext = n;
	if ((nitmax+1)*nrhs > nlocext) nlocext = (nitmax+1)*nrhs;
//
// Allocate temporary vector data
//
	CSVectorC r(n,npartsloc,nrhs), y(n,npartsloc,nrhs), w(n,npartsloc,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *yp, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.ooctyp == 0) {
		p = new dcmplx [ldp*(nitmax+1)*nrhs];
		if (!p) MemoryFail (funcname);
	} else {
		p = new dcmplx [2*ldp*nrhs];
		if (!p) MemoryFail (funcname);
	};

	yp = new dcmplx [ldp*nrhs];
	if (!yp) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

//	treeclp = _tree.CollapseTree ();
	treeclp = _tree;
//
// Create QRD structure
//
	CQrdC qrd (treeclp, nitmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.ooctyp, strbuff, ldp);

	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	(_mvmc) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

//	_fout << " Residual vector = " << r << endl;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			p[irhs*nlocext+kii] = pr[irhs*n+kii];
		};
	};

	qrd.UpdateQrd (treeclp, nlocext, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs_2;

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
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi0[irhs] > resimax) resimax = resi0[irhs];
		};
		cout  << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
		_fout << " Initial Log10 || ResiBlk || = " << log10(resimax) << endl;
	};
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp, iconv;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;

//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		time2 = clock ();

//		_fout << " Before solve by precond y = " << y << endl;

		(_solvelu) (_objlu, y, w);
		y = w;

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;
		nslvu++;

		time2 = clock ();

//		_fout << " Y bef mvm = " << y << endl;

//		_gmtra.MvmARect (y,w);
		(_mvmc) (_obja, y, w);
//		if (true) throw " After mvm ";

//		_fout << " W aft mvm = " << w << endl;

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;
//
// Local computations on the root node
//
		if (_tree.myid == _tree.rootcpu) {
//
// Take current block column of H_k
//
			int ibeg=ncolp;
			int iend=ncolp+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

			qrd.GetRPartQrd (treeclp, 0, ibeg, iend,
								eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (int igiv=0;igiv<k;igiv++) {

				ApplyBlockGivens (nrhs, 
									igiv, givarr+igiv*nrhs_2*4,
									nrhs, eloc, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
							giv, taugiv, giv1, givarr+k*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<(k+1)*nrhs;kii++) {
					rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, bloc, ldbloc, giv);

			ops += 8*nrhs_3;
//
// Check convergence
//
			double resimax = 0.0e0;

			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);
				if (resi[irhs] > resimax) resimax = resi[irhs];

			};

			if (iter%_param.ichk == 0) {
				cout  << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
				_fout << " Iter = " << k << " Log10 || Resiblk || = " << log10(resimax) << endl;
			};

			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

		} else {
			iconv = 0;
		};

// Add convergence control numbers on all processors

		int iconvtot = 0;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &iconv, &iconvtot);
		} else {
			iconvtot = iconv;
		};

// Check convergence and exit

		if (iconvtot == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	if (true && _tree.myid == 0) {
		double *hystarr;
		hystarr = new double [ncolp];
		if (!hystarr) MemoryFail (funcname);
		dcmplx caux;
		for (kii=0;kii<ncolp;kii++) {
			caux = rmatr[kii*ldrmatr+kii];
			hystarr[kii] = sqrt(caux.x*caux.x+caux.y*caux.y);
		};
		csvhyst (cout,  "DiagRGmres", ncolp, hystarr);
		csvhyst (_fout, "DiagRGmres", ncolp, hystarr);
		delete [] hystarr;
	};

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	qrd.MvmQBlk (treeclp, nlocext, ncolp+nrhs, nrhs,
					p, ldp, taup, 
					eloc, ldeloc, 
					yp, nlocext);

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			py[irhs*n+kii] = yp[irhs*nlocext+kii];
		};
	};

	ops += 2*nlocext*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	time2 = clock ();

	(_solvelu) (_objlu, y, w);
	y = w;

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;
	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	time2 = clock ();

//	_gmtra.MvmARect (sol,r);
	(_mvmc) (_obja, sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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

	ops += 2*n*nrhs;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi0);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi0[irhs] = resi[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi0[irhs];
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		double resimax = 0.0e0;
		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs] > resimax) resimax = resi[irhs];
		};
		cout  << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
		_fout << " Final Log10 || Resiblk || = " << log10(resimax) << endl;
	};
//
// Compute SVD of R_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrd.CloseFiles ();

	};
//
// Finalize time measurement
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

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

	double dnza = arr1[0];
	ops = arr1[1];

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	if (dnza > 0.0e0) {
		ops = ops / dnza;
	} else {
		ops = 0.0e0;
	};

	if (_tree.myid == _tree.rootcpu) {

		cout  << " GMRES statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;

		_fout << " GMRES statistics: " << endl;
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

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							int _ivar, 
							CGSMatrixCS &_gmtrau, 
							double _dshift, CGSMatrixCS &_gmtrarect, int *_order, 
							CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int n = _gmtrau.GetN();
//
// Init statistics data
//
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
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	p = new dcmplx [ldp*(nitmax+1)*nrhs];
	if (!p) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	if (_ivar == 0) {
		_gmtrau.MvmASymm (sol,r);
	} else {
		_gmtrarect.MvmAhAOrdShiftSymm (_dshift, _order, sol,r);
	};

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) p[kii] = pr[kii];

	UpdateQrd (n, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, nrhs, 
					p, ldp, bloc, ldbloc);
//	OutArr(_fout," bloc on return from getr ",ldbloc*nrhs,bloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, ncolp, nrhs,
					p, ldp, taup, 
					eloc, ncolp, 
					py, n);

		ops += 2*n*ncolp*nrhs;
//		_fout << " y on return from MvmQBlk" << y << endl;
//
// Multiply by the preconditioned matrix
//
		_gmtru.SolveLConj (y,w);
//		_fout << " w on return from SolveL" << w << endl;
		_gmtru.SolveU (w,y);
//		_fout << " y on return from SolveU" << y << endl;

		if (_ivar == 0) {
			_gmtrau.MvmASymm (y,w);
		} else {
			_gmtrarect.MvmAhAOrdShiftSymm (_dshift, _order, y,w);
		};

//		_fout << " w on return from MvmA" << w << endl;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) p[n*ncolp+kii] = pw[kii];

		UpdateQrd (n, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, ibeg, iend,
						p, ldp, eloc, iend);
//		OutArr(_fout," R aft update = ",iend*nrhs,eloc);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, ncolp+nrhs, nrhs,
				p, ldp, taup, 
				eloc, ldeloc, 
				py, n);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_gmtru.SolveLConj (y,w); 
	_gmtru.SolveU (w,y); 

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	if (_ivar == 0) {
		_gmtrau.MvmASymm (sol,r);
	} else {
		_gmtrarect.MvmAhAOrdShiftSymm (_dshift, _order, sol,r);
	};

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nza;

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] p;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
										void *_obj, CGETM _getm, CGETN _getn, CMVMC _mvmc, CMVMHC _mvmhc,
										double _dshift, int *_order, 
										CGSMatrixCS &_gmtru,
										const CSVectorC &_x, const CSVectorC &_b, 
										const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;

//
// Get the size of the matrix
//
	int n = (_getn) (_obj);
//
// Init statistics data
//
	time0 = clock ();

	ops = 0.0e0;

	nza = 0;

	nzlu = 0;

	int ilist;
	for (ilist=0;ilist<_gmtru.nlist;ilist++) {
		int iblk = _gmtru.listb[ilist];
		nzlu += _gmtru.mtrarr[iblk].GetNzatot();
	};
//
// Init guess to the solution
//
	CSVectorC sol(n,1,nrhs);

	sol = _x;
//
// Allocate temporary vector data
//
	CSVectorC r(n,1,nrhs), y(n,1,nrhs), w(n,1,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int nitmax = _param.niter;

	int ldp = n;
	p = new dcmplx [ldp*(nitmax+1)*nrhs];
	if (!p) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	MvmAhAOrdShiftSymm (_obj, _getm, _getn, _mvmc, _mvmhc,
								_dshift, _order, sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<n*nrhs;kii++) p[kii] = pr[kii];

	UpdateQrd (n, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	GetRPartQrd (n, 0, nrhs, 
					p, ldp, bloc, ldbloc);
//	OutArr(_fout," bloc on return from getr ",ldbloc*nrhs,bloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
		cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
	};

	ops += 2*nrhs_2;
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		MvmQBlk (n, ncolp, nrhs,
					p, ldp, taup, 
					eloc, ncolp, 
					py, n);

		ops += 2*n*ncolp*nrhs;
//		_fout << " y on return from MvmQBlk" << y << endl;
//
// Multiply by the preconditioned matrix
//
		_gmtru.SolveLConj (y,w);
//		_fout << " w on return from SolveL" << w << endl;
		_gmtru.SolveU (w,y);
//		_fout << " y on return from SolveU" << y << endl;

		MvmAhAOrdShiftSymm (_obj, _getm, _getn, _mvmc, _mvmhc,
									_dshift, _order, y,w);

//		_fout << " w on return from MvmA" << w << endl;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		for (kii=0;kii<n*nrhs;kii++) p[n*ncolp+kii] = pw[kii];

		UpdateQrd (n, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*n*ncolp*nrhs + n*nrhs_2;
//
// Take current block column of H_k
//
		int ibeg=ncolp;
		int iend=ncolp+nrhs;

		for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

		GetRPartQrd (n, ibeg, iend,
						p, ldp, eloc, iend);
//		OutArr(_fout," R aft update = ",iend*nrhs,eloc);
//
// Apply previous block Givens rotations to the current block column of H_k
//
		for (int igiv=0;igiv<k;igiv++) {

			ApplyBlockGivens (nrhs, 
								igiv, givarr+igiv*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;

		};
//
// Compute new block Givens rotation
//
		for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<nrhs*2;kii++) {
				giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
			};
		};

		BlockGivens (nrhs, 
						giv, taugiv, giv1, givarr+k*nrhs_2*4);

		ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, eloc, iend, giv);

		ops += 8*nrhs_3;
//
// Store transformed block column in R
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<(k+1)*nrhs;kii++) {
				rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
			};
		};
//
// Apply new block Givens rotation to the reduced block rhs
//
		ApplyBlockGivens (nrhs, 
							k, givarr+k*nrhs_2*4,
							nrhs, bloc, ldbloc, giv);

		ops += 8*nrhs_3;
//
// Check convergence
//
		for (irhs=0;irhs<nrhs;irhs++) {

			double aux = 0.0e0;

			for (kii=0;kii<nrhs;kii++) {
				double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
				aux += auxl*auxl;
				auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
				aux += auxl*auxl;
			};
			resi[irhs] = sqrt(aux);

			if (iter%_param.ichk == 0) {
				cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
			};

		};

		int iconv = 0;

		for (irhs=0;irhs<nrhs;irhs++) {
			if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
		};

		if (iconv == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	MvmQBlk (n, ncolp+nrhs, nrhs,
				p, ldp, taup, 
				eloc, ldeloc, 
				py, n);

	ops += 2*n*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	_gmtru.SolveLConj (y,w); 
	_gmtru.SolveU (w,y); 

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	MvmAhAOrdShiftSymm (_obj, _getm, _getn, _mvmc, _mvmhc,
								_dshift, _order, sol,r);

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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
		resi[irhs] = sqrt(aux);
		cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
	};

	ops += 2*n*nrhs;
//
// Compute SVD of R_k
//
	if (true) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	if (nza != 0) {
		ops = ops / (double) nza;
	} else {
		ops = 0.0e0;
	};

	cout  << " GMRES statistics: " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " GMRES statistics: " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Free work arrays

	delete [] p;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix in parallel mode
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
							CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
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
	int npartsloc = _x.GetNparts ();

	CSVectorC sol(n,npartsloc,nrhs);

	sol = _x;
//
// Determine extended length
//
	int nitmax = _param.niter;

	int nlocext = n;
	if ((nitmax+1)*nrhs > nlocext) nlocext = (nitmax+1)*nrhs;
//
// Allocate temporary vector data
//
	CSVectorC r(n,npartsloc,nrhs), y(n,npartsloc,nrhs), w(n,npartsloc,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *yp, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.ooctyp == 0) {
		p = new dcmplx [ldp*(nitmax+1)*nrhs];
		if (!p) MemoryFail (funcname);
	} else {
		p = new dcmplx [2*ldp*nrhs];
		if (!p) MemoryFail (funcname);
	};

	yp = new dcmplx [ldp*nrhs];
	if (!yp) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

//	treeclp = _tree.CollapseTree ();
	treeclp = _tree;
//
// Create QRD structure
//
	CQrdC qrd (treeclp, nitmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.ooctyp, strbuff, ldp);

	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	gmtrdummy.MvmA (_tree, _mvm, _gmtral, _gmtrau, sol,r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			p[irhs*nlocext+kii] = pr[irhs*n+kii];
		};
	};

	qrd.UpdateQrd (treeclp, nlocext, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs_2;

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
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp, iconv;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		time2 = clock ();

		_gmtrl.SolveL (_tree,_mvm,y,w);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		time2 = clock ();

		_gmtru.SolveU (_tree,_mvm,w,y);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,y,w);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;
//
// Local computations on the root node
//
		if (_tree.myid == _tree.rootcpu) {
//
// Take current block column of H_k
//
			int ibeg=ncolp;
			int iend=ncolp+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

			qrd.GetRPartQrd (treeclp, 0, ibeg, iend,
								eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (int igiv=0;igiv<k;igiv++) {

				ApplyBlockGivens (nrhs, 
									igiv, givarr+igiv*nrhs_2*4,
									nrhs, eloc, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
							giv, taugiv, giv1, givarr+k*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<(k+1)*nrhs;kii++) {
					rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, bloc, ldbloc, giv);

			ops += 8*nrhs_3;
//
// Check convergence
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				if (iter%_param.ichk == 0) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				};

			};

			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

		} else {
			iconv = 0;
		};

// Add convergence control numbers on all processors

		int iconvtot = 0;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &iconv, &iconvtot);
		} else {
			iconvtot = iconv;
		};

// Check convergence and exit

		if (iconvtot == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	qrd.MvmQBlk (treeclp, nlocext, ncolp+nrhs, nrhs,
					p, ldp, taup, 
					eloc, ldeloc, 
					yp, nlocext);

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			py[irhs*n+kii] = yp[irhs*nlocext+kii];
		};
	};

	ops += 2*nlocext*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	time2 = clock ();

	_gmtrl.SolveL (_tree,_mvm,y,w); 

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	time2 = clock ();

	_gmtru.SolveU (_tree,_mvm,w,y); 

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	time2 = clock ();

	gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,sol,r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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

	ops += 2*n*nrhs;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi0);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi0[irhs] = resi[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi0[irhs];
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		};
	};
//
// Compute SVD of R_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrd.CloseFiles ();

	};
//
// Finalize time measurement
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

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

	double dnza = arr1[0];
	ops = arr1[1];

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " GMRES statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;

		_fout << " GMRES statistics: " << endl;
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

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// Description: Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix in parallel mode
// CSVectorC::BlockGMRES()
//========================================================================================
void CSVectorC::BlockGMRES (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtrau, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param) { 

	const char *funcname = "BlockGMRES";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmva=0, nslvl=0, nslvu=0;
//
// Get the size of the matrix
//
//	int ntotal = _gmtrau.GetN();
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
	int npartsloc = _x.GetNparts ();

	CSVectorC sol(n,npartsloc,nrhs);

	sol = _x;
//
// Determine extended length
//
	int nitmax = _param.niter;

	int nlocext = n;
	if ((nitmax+1)*nrhs > nlocext) nlocext = (nitmax+1)*nrhs;
//
// Allocate temporary vector data
//
	CSVectorC r(n,npartsloc,nrhs), y(n,npartsloc,nrhs), w(n,npartsloc,nrhs);
//
// Allocate all working arrays
//
	dcmplx *p, *yp, *taup;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.ooctyp == 0) {
		p = new dcmplx [ldp*(nitmax+1)*nrhs];
		if (!p) MemoryFail (funcname);
	} else {
		p = new dcmplx [2*ldp*nrhs];
		if (!p) MemoryFail (funcname);
	};

	yp = new dcmplx [ldp*nrhs];
	if (!yp) MemoryFail (funcname);
	taup = new dcmplx [(nitmax+1)*nrhs];
	if (!taup) MemoryFail (funcname);

	givarr = new dcmplx [nitmax*nrhs_2*4];
	if (!givarr) MemoryFail (funcname);
	giv = new dcmplx [nrhs_2*4];
	if (!giv) MemoryFail (funcname);
	giv1 = new dcmplx [nrhs_2*4];
	if (!giv1) MemoryFail (funcname);
	taugiv = new dcmplx [nrhs*2];
	if (!taugiv) MemoryFail (funcname);

	int ldbloc = (nitmax+1)*nrhs;
	bloc = new dcmplx [ldbloc*nrhs];
	if (!bloc) MemoryFail (funcname);
	eloc = new dcmplx [ldbloc*nrhs];
	if (!eloc) MemoryFail (funcname);
	int ldrmatr = nitmax*nrhs;
	rmatr = new dcmplx [ldrmatr*ldrmatr];
	if (!rmatr) MemoryFail (funcname);

	resi0 = new double [nrhs];
	if (!resi0) MemoryFail (funcname);
	resi = new double [nrhs];
	if (!resi) MemoryFail (funcname);
//
// Create collapsed tree
//
	CTree treeclp;

	treeclp = _tree.CollapseTree ();
//
// Create QRD structure
//
	CQrdC qrd (treeclp, nitmax+1, nrhs);
//
// Setup files if necessary
//
	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.ooctyp, strbuff, ldp);

	};
//
// Compute initial residual and its norm
//
	CGSMatrixCS gmtrdummy;

	time2 = clock ();

	gmtrdummy.MvmASymm (_tree, _mvm, _gmtrau, sol, r);
//	_gmtrau.MvmASymm (sol, r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

	ops += n*nrhs;
//
// Compute QR decomposition of the residual block
//
	int kii, irhs;

	dcmplx *pr = r.GetVect ();

	for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = czero;

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			p[irhs*nlocext+kii] = pr[irhs*n+kii];
		};
	};

	qrd.UpdateQrd (treeclp, nlocext, 0, nrhs, p, ldp, taup);

	ops += n*nrhs_2;
//
// Init the local right hand side
//
	for (kii=0;kii<ldbloc*nrhs;kii++) bloc[kii] = czero;

	qrd.GetRPartQrd (treeclp, 0, 
						0, nrhs, 
						bloc, ldbloc);
//
// Compute array of initial residual norms
//
	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = 0.0e0;
		for (kii=0;kii<nrhs;kii++) {
			double auxl = bloc[irhs*ldbloc+kii].real ();
			aux += auxl*auxl;
			auxl = bloc[irhs*ldbloc+kii].imag ();
			aux += auxl*auxl;
		};
		resi0[irhs] = aux;
	};

	ops += 2*nrhs_2;

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
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Initial Log10 || Resi || = " << log10(resi0[irhs]) << endl;
		};
	};
//
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int ncolp, iconv;
	int iter = 0;

	for (int k=0;k<nitmax;k++) {

		iter = k;

		ncolp = (k+1)*nrhs;
//
// Compute current initial directions block
//
		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[irhs*ncolp+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+k*nrhs+irhs] = cone;
		};

		dcmplx *py = y.GetVect ();

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;
//
// Multiply by the preconditioned matrix
//
		time2 = clock ();

		_gmtru.SolveLConj (_tree,_mvm,y,w);

		time3 = clock ();

		timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvl++;

		time2 = clock ();

		_gmtru.SolveU (_tree,_mvm,w,y);

		time3 = clock ();

		timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nslvu++;

		time2 = clock ();

		gmtrdummy.MvmASymm (_tree,_mvm,_gmtrau,y,w);
//		_gmtrau.MvmASymm (y,w);

		time3 = clock ();

		timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmva++;

		ops += (nza + 2*nzlu)*nrhs;
//
// Update QR decomposition for the result
//
		dcmplx *pw = w.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;
//
// Local computations on the root node
//
		if (_tree.myid == _tree.rootcpu) {
//
// Take current block column of H_k
//
			int ibeg=ncolp;
			int iend=ncolp+nrhs;

			for (kii=0;kii<iend*nrhs;kii++) eloc[kii] = czero;

			qrd.GetRPartQrd (treeclp, 0, ibeg, iend,
								eloc, iend);
//
// Apply previous block Givens rotations to the current block column of H_k
//
			for (int igiv=0;igiv<k;igiv++) {

				ApplyBlockGivens (nrhs, 
									igiv, givarr+igiv*nrhs_2*4,
									nrhs, eloc, iend, giv);

				ops += 8*nrhs_3;

			};
//
// Compute new block Givens rotation
//
			for (kii=0;kii<nrhs_2*4;kii++) giv[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<nrhs*2;kii++) {
					giv[irhs*nrhs*2+kii] = eloc[irhs*iend+k*nrhs+kii];
				};
			};

			BlockGivens (nrhs, 
							giv, taugiv, giv1, givarr+k*nrhs_2*4);

			ops += 8*nrhs_3;
//
// Apply new block Givens rotation to the transformed block column of H_k
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, eloc, iend, giv);

			ops += 8*nrhs_3;
//
// Store transformed block column in R
//
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<(k+1)*nrhs;kii++) {
					rmatr[(k*nrhs+irhs)*ldrmatr+kii] = eloc[irhs*iend+kii];
				};
			};
//
// Apply new block Givens rotation to the reduced block rhs
//
			ApplyBlockGivens (nrhs, 
								k, givarr+k*nrhs_2*4,
								nrhs, bloc, ldbloc, giv);

			ops += 8*nrhs_3;
//
// Check convergence
//
			for (irhs=0;irhs<nrhs;irhs++) {

				double aux = 0.0e0;

				for (kii=0;kii<nrhs;kii++) {
					double auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].real ();
					aux += auxl*auxl;
					auxl = bloc[irhs*ldbloc+(k+1)*nrhs+kii].imag ();
					aux += auxl*auxl;
				};
				resi[irhs] = sqrt(aux);

				if (iter%_param.ichk == 0) {
					cout  << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
					_fout << " Irhs = " << irhs << " Iter = " << k << " Log10 || Resi || = " << log10(resi[irhs]) << endl;
				};

			};

			iconv = 0;

			for (irhs=0;irhs<nrhs;irhs++) {
				if (resi[irhs]/resi0[irhs] < _param.eps) iconv++;
			};

		} else {
			iconv = 0;
		};

// Add convergence control numbers on all processors

		int iconvtot = 0;

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														INTEGERVALUE, ADD,
														1, &iconv, &iconvtot);
		} else {
			iconvtot = iconv;
		};

// Check convergence and exit

		if (iconvtot == nrhs) goto exitgmr;

		if (k > _param.niter) goto exitgmr;

	};
//
// Update solution
//
exitgmr:;
//
// Solve triangular system
//
	ncolp = (iter+1)*nrhs;

	int ldeloc = ncolp + nrhs;

	for (kii=0;kii<nrhs*ldeloc;kii++) eloc[kii] = czero;

	SolveTriangular (ncolp, nrhs,
						rmatr, ldrmatr,
						bloc, ldbloc,
						eloc, ldeloc);
//
// Compute directions block
//
	dcmplx *py = y.GetVect ();

	qrd.MvmQBlk (treeclp, nlocext, ncolp+nrhs, nrhs,
					p, ldp, taup, 
					eloc, ldeloc, 
					yp, nlocext);

	for (irhs=0;irhs<nrhs;irhs++) {
		for (kii=0;kii<n;kii++) {
			py[irhs*n+kii] = yp[irhs*nlocext+kii];
		};
	};

	ops += 2*nlocext*(ncolp+nrhs)*nrhs;
//
// Update solution block
//
	time2 = clock ();

	_gmtru.SolveLConj (_tree,_mvm,y,w); 

	time3 = clock ();

	timeslvl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvl++;

	time2 = clock ();

	_gmtru.SolveU (_tree,_mvm,w,y); 

	time3 = clock ();

	timeslvu += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nslvu++;

	ops += 2*nzlu*nrhs;
//
// Compute the final residual
//
	sol.DaxpySVect (cone,y);

	ops += n*nrhs;

	time2 = clock ();

	gmtrdummy.MvmASymm (_tree,_mvm,_gmtrau,sol,r);
//	_gmtrau.MvmASymm (sol,r);

	time3 = clock ();

	timemva += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmva++;

	ops += nza*nrhs;

	y = _b;

	r = y-r;

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

	ops += 2*n*nrhs;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi, resi0);
	} else {
		for (irhs=0;irhs<nrhs;irhs++) resi0[irhs] = resi[irhs];
	};

	for (irhs=0;irhs<nrhs;irhs++) {
		double aux = resi0[irhs];
		resi0[irhs] = sqrt(aux);
		resi[irhs] = resi0[irhs];
	};

	if (_tree.myid == _tree.rootcpu) {
		for (irhs=0;irhs<nrhs;irhs++) {
			cout  << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
			_fout << " Irhs = " << irhs << " Final Log10 || Resi || = " << log10(resi[irhs]) << endl;
		};
	};
//
// Compute SVD of R_k
//
	if (true && _tree.myid == _tree.rootcpu) {

		int ntot = ncolp;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rmatr, &ldrmatr, svloc, rmatr, &ldrmatr, rmatr, &ldrmatr, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};
//
// Close files if necessary
//
	if (_param.ooctyp > 0) {

		qrd.CloseFiles ();

	};
//
// Finalize time measurement
//
	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

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

	double dnza = arr1[0];
	ops = arr1[1];

// Output GMRES statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnza;

	if (_tree.myid == _tree.rootcpu) {

		cout  << " GMRES statistics: " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		cout  << " Mv-Slv functions statistics: " << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< endl;

		_fout << " GMRES statistics: " << endl;
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

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi;

// Return the solution

	*this = sol;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
