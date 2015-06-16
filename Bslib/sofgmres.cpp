//------------------------------------------------------------------------------------------------
// File: sofgmres.cpp
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
#include "corr.h"
#include "globals.h"

using namespace std;

// Interface routines 

void BlockGivens (int _ngiv, // Compute block Givens rotation
					dcmplx *_giv, dcmplx *_taugiv, 
					dcmplx *_aloc, dcmplx *_givarr);

void ApplyBlockGivens (int _ngiv, // Apply block Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc);

void ApplyBlockGivensH (int _ngiv, // Apply Hermitian transposed block Givens rotation
						int _indgiv, dcmplx *_giv, 
						int _nrhs, dcmplx *_rhs, int _ldrhs,
						dcmplx *_aloc);

void SolveTriangular (int _n, int _nrhs, // Solve system with triangular dense matrix
						dcmplx *_rmatr, int _ldrmatr,
						dcmplx *_rhs, int _ldrhs,
						dcmplx *_sol, int _ldsol);

// Author: Kharchenko S.A.
// CSVectorC: Compute corrector subspace directions
//========================================================================================
void CSVectorC::BlockCorrector (ofstream &_fout, // Compute corrector subspace directions
								const CTree &_tree, 
								int _n, int _nlocext, int _nrhs, int _iter, 
								CQrdC &_qrd, dcmplx *_p, int _ldp, dcmplx *_taup,
								int _ldrmatr, dcmplx *_rmatr, dcmplx *_givarr,
								CSVectorC &_ycorr,
								const CSlvParam &_param, double &_ops) { 

	const char *funcname = "BlockCorrector";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);
	dcmplx chalph (0.5e0,0.0e0);

// Allocate work data

	int kloc = _iter * _nrhs;
	int kextloc = (_iter+1) * _nrhs;

	dcmplx *tplus, *pp;
	dcmplx *aloc, *uloc, *vloc, *work;
	double *svloc, *rwork;
	dcmplx *ytemp;

	int lwork = 10 * kextloc;

	tplus = new dcmplx [kextloc*kextloc];
	if (!tplus) MemoryFail (funcname);
	pp = new dcmplx [kextloc*kextloc];
	if (!pp) MemoryFail (funcname);
	aloc = new dcmplx [kextloc*kextloc];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [kextloc*kextloc];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [kextloc*kextloc];
	if (!vloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);

	svloc = new double [kextloc];
	if (!svloc) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

	int nrhs_2 = _nrhs * _nrhs;

	int kii, kjj, kkk, info, irhs, igiv;

	int ndir;

	int ndirtot;

	CSVectorC vectdummy;

// Compute matrix T_+

	if (_tree.myid == _tree.rootcpu) {

		for (kii=0;kii<kextloc*kextloc;kii++) aloc[kii] = czero;

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kextloc+kii] = conj(_rmatr[kii*_ldrmatr+kjj]);
			};
		};

		for (igiv=0;igiv<_iter;igiv++) {

			ApplyBlockGivens (_nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								kloc, aloc, kextloc, uloc);

		};

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				tplus[kjj*kloc+kii] = (aloc[kjj*kextloc+kii] + conj (aloc[kii*kextloc+kjj])) * chalph;
			};
		};

// Compute singular value and spectral decompositions
// and determine the block to be multiplied

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kloc+kii] = _rmatr[kjj*_ldrmatr+kii];
			};
		};

		zgesvd_ ("A", "A", &kloc, &kloc, 
					aloc, &kloc, svloc, uloc, &kloc, vloc, &kloc, 
					work, &lwork, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		_ops += 8 * kloc * kloc * kloc;

		csvhyst (cout,  "CorrSvR", kloc, svloc);
		csvhyst (_fout, "CorrSvR", kloc, svloc);

		ndir = 0;

		for (kii=0;kii<kloc;kii++) {
			if (svloc[kii] >= _param.xmax) {
				for (kjj=0;kjj<kloc;kjj++) {
					pp[ndir*kloc+kjj] = conj(vloc[kjj*kloc+kii]);
				};
				ndir++;
			};
		};

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kloc+kii] = tplus[kjj*kloc+kii];
			};
		};

		zheev_ ("V", "U", &kloc, 
					aloc, &kloc, svloc,
					work, &lwork, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZHEEV" << endl;
			throw " Error in the Lapack routine ZHEEV";
		};

		_ops += 8 * kloc * kloc * kloc;

		csvhyst (cout,  "CorrEigT+", kloc, svloc);
		csvhyst (_fout, "CorrEigT+", kloc, svloc);

		kkk = 0;

		for (kii=0;kii<kloc;kii++) {
			if (svloc[kii] <= _param.xmin) {
				kkk++;
			};
		};

		ndirtot = kkk+ndir;

		if (ndirtot > kloc) kkk = kloc - ndir;

		ndirtot = kkk+ndir;

		kii = (kloc - ndirtot) % _nrhs;

		kkk += kii;

		ndirtot = kkk+ndir;

		if (ndirtot % _nrhs != 0) {
			throw " Wrong number of direction vectors";
		};

		for (kii=0;kii<kkk;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				pp[ndir*kloc+kjj] = aloc[kii*kloc+kjj];
			};
			ndir++;
		};

// Orthonormalize directions

		if (ndir != 0) {

			zgesvd_ ("O", "A", &kloc, &ndir,
						pp, &kloc, svloc, pp, &kloc, vloc, &kloc, 
						work, &lwork, rwork, &info);

			if (info != 0) {
				cout << " Error in the Lapack routine ZGESVD" << endl;
				throw " Error in the Lapack routine ZGESVD";
			};

			_ops += 8 * kloc * kloc * kloc;

		};

	} else {

		ndir = 0;

	};

// Exchange directions data among processors

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													1, &ndir, &ndirtot);
		ndir = ndirtot;
	};

	if (ndir == 0) goto ExitCorr;

	if (_tree.myid != _tree.rootcpu) {
		for (kii=0;kii<ndir*kloc;kii++) pp[kii] = czero;
	};

	if (_tree.nproc != 1) {

		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													ndir*kloc*2, pp, aloc);

		for (kii=0;kii<ndir*kloc;kii++) pp[kii] = aloc[kii];

	};

// Multiply and locally store

	{

		ytemp = new dcmplx [_nlocext*ndir];
		if (!ytemp) MemoryFail (funcname);

		for (int icol=0;icol<ndir/_nrhs;icol++) {

			_qrd.MvmQBlk (_tree, _nlocext, kloc, _nrhs,
							_p, _ldp, _taup, 
							pp+kloc*icol*_nrhs, kloc,
							ytemp+_nlocext*icol*_nrhs, _nlocext);

		};

		_ops += 2*_nlocext*kloc*ndir;

// Store the result in ycorr along with the previous directions

		CSVectorC ywork(_n, 1, ndir+_ycorr.nrhs);

		kkk = 0;

		for (irhs=0;irhs<_ycorr.nrhs;irhs++) {
			for (kii=0;kii<_n;kii++) {
				ywork.vect[kkk*_n+kii] = _ycorr.vect[irhs*_n+kii];
			};
			kkk++;
		};

		for (irhs=0;irhs<ndir;irhs++) {
			for (kii=0;kii<_n;kii++) {
				ywork.vect[kkk*_n+kii] = ytemp[irhs*_nlocext+kii];
			};
			kkk++;
		};

		delete [] ytemp;

		_ycorr = ywork;

		ywork = vectdummy;

	};

// Free work data

ExitCorr:;

	delete [] tplus;
	delete [] pp;
	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;

};

// Author: Kharchenko S.A.
// CSVectorC: Compute corrector
//========================================================================================
void CSVectorC::BlockCorrector (ofstream &_fout, // Compute corrector
								const CTree &_tree, 
								int _n, int _nlocext, int _nrhs, int _iter, 
								CQrdC &_qrd, dcmplx *_p, int _ldp, dcmplx *_taup,
								int _ldrmatr, dcmplx *_rmatr, dcmplx *_givarr,
								CCorrectorC &_corr,
								const CSlvParam &_param, double &_ops) { 

	const char *funcname = "BlockCorrector";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);
	dcmplx chalph (0.5e0,0.0e0);
	dcmplx cmone (-1.0e0,0.0e0);

// Allocate work data

	int kloc = _iter * _nrhs;
	int kextloc = (_iter+1) * _nrhs;

	dcmplx *tplus, *qloc, *tauloc, *pp;
	dcmplx *aloc, *rloc, *uloc, *vloc, *work;
	double *svloc, *rwork;
	dcmplx *ytemp;

	int lwork = 10 * kextloc;

	tplus = new dcmplx [kextloc*kextloc];
	if (!tplus) MemoryFail (funcname);
	qloc = new dcmplx [kextloc*kextloc];
	if (!qloc) MemoryFail (funcname);
	tauloc = new dcmplx [kextloc];
	if (!tauloc) MemoryFail (funcname);
	pp = new dcmplx [kextloc*kextloc];
	if (!pp) MemoryFail (funcname);
	aloc = new dcmplx [kextloc*kextloc];
	if (!aloc) MemoryFail (funcname);
	rloc = new dcmplx [kextloc*kextloc];
	if (!rloc) MemoryFail (funcname);
	uloc = new dcmplx [kextloc*kextloc];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [kextloc*kextloc];
	if (!vloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);

	svloc = new double [kextloc];
	if (!svloc) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

	int nrhs_2 = _nrhs * _nrhs;

	int kii, kjj, kkk, info, irhs, igiv;

	int ndir;

	int ndirtot;

// Compute matrix T_+

	if (_tree.myid == _tree.rootcpu) {

		for (kii=0;kii<kextloc*kextloc;kii++) aloc[kii] = czero;

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kextloc+kii] = conj(_rmatr[kii*_ldrmatr+kjj]);
			};
		};

		for (igiv=0;igiv<_iter;igiv++) {

			ApplyBlockGivens (_nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								kloc, aloc, kextloc, uloc);

		};

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				tplus[kjj*kloc+kii] = (aloc[kjj*kextloc+kii] + conj (aloc[kii*kextloc+kjj])) * chalph;
			};
		};

// Compute singular value and spectral decompositions
// and determine the block to be multiplied

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kloc+kii] = _rmatr[kjj*_ldrmatr+kii];
			};
		};

		zgesvd_ ("A", "A", &kloc, &kloc, 
					aloc, &kloc, svloc, uloc, &kloc, vloc, &kloc, 
					work, &lwork, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		_ops += 8 * kloc * kloc * kloc;

		csvhyst (cout,  "CorrSvR", kloc, svloc);
		csvhyst (_fout, "CorrSvR", kloc, svloc);

		ndir = 0;

		for (kii=0;kii<kloc;kii++) {
			if (svloc[kii] >= _param.xmax) {
				for (kjj=0;kjj<kloc;kjj++) {
					pp[ndir*kloc+kjj] = conj(vloc[kjj*kloc+kii]);
				};
				ndir++;
			};
		};

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				aloc[kjj*kloc+kii] = tplus[kjj*kloc+kii];
			};
		};

		zheev_ ("V", "U", &kloc, 
					aloc, &kloc, svloc,
					work, &lwork, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZHEEV" << endl;
			throw " Error in the Lapack routine ZHEEV";
		};

		_ops += 8 * kloc * kloc * kloc;

		csvhyst (cout,  "CorrEigT+", kloc, svloc);
		csvhyst (_fout, "CorrEigT+", kloc, svloc);

		kkk = 0;

		for (kii=0;kii<kloc;kii++) {
			if (svloc[kii] <= _param.xmin) {
				kkk++;
			};
		};

		ndirtot = kkk+ndir;

		if (ndirtot > kloc) kkk = kloc - ndir;

		ndirtot = kkk+ndir;

		kii = (kloc - ndirtot) % _nrhs;

		kkk += kii;

		ndirtot = kkk+ndir;

		if (ndirtot % _nrhs != 0) {
			throw " Wrong number of direction vectors";
		};

		for (kii=0;kii<kkk;kii++) {
			for (kjj=0;kjj<kloc;kjj++) {
				pp[ndir*kloc+kjj] = aloc[kii*kloc+kjj];
			};
			ndir++;
		};

// Orthonormalize directions

		if (ndir != 0) {

			zgesvd_ ("O", "A", &kloc, &ndir,
						pp, &kloc, svloc, pp, &kloc, vloc, &kloc, 
						work, &lwork, rwork, &info);

			if (info != 0) {
				cout << " Error in the Lapack routine ZGESVD" << endl;
				throw " Error in the Lapack routine ZGESVD";
			};

			_ops += 8 * kloc * ndir * ndir;

		};

	} else {

		ndir = 0;

	};

// Exchange directions data among processors

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													1, &ndir, &ndirtot);
		ndir = ndirtot;
	};

	if (ndir == 0) goto ExitCorr;

	if (_tree.myid != _tree.rootcpu) {
		for (kii=0;kii<ndir*kloc;kii++) pp[kii] = czero;
	};

	if (_tree.nproc != 1) {

		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													ndir*kloc*2, pp, aloc);

		for (kii=0;kii<ndir*kloc;kii++) pp[kii] = aloc[kii];

	};

// Multiply and locally store

	{

// Allocate memory for corrector data and register

		dcmplx *ycorr, *wcorr, *dcorr;

		ycorr = new dcmplx [_n*ndir];
		if (!ycorr) MemoryFail (funcname);
		wcorr = new dcmplx [_n*ndir];
		if (!wcorr) MemoryFail (funcname);
		dcorr = new dcmplx [4*ndir*ndir];
		if (!dcorr) MemoryFail (funcname);

		if (_corr.ncorr >= _corr.ncorrmax-1) {
			throw " Error in the number of correctors";
		};

		_corr.ldcorrarr[_corr.ncorr] = _n;
		_corr.corrsizes[_corr.ncorr] = ndir;

		_corr.yarr[_corr.ncorr] = ycorr;
		_corr.warr[_corr.ncorr] = wcorr;
		_corr.darr[_corr.ncorr] = dcorr;

		_corr.ncorr++;

// Compute Y part

		ytemp = new dcmplx [_nlocext*ndir];
		if (!ytemp) MemoryFail (funcname);

		int icol;

		for (icol=0;icol<ndir/_nrhs;icol++) {

			_qrd.MvmQBlk (_tree, _nlocext, kloc, _nrhs,
							_p, _ldp, _taup, 
							pp+kloc*icol*_nrhs, kloc,
							ytemp+_nlocext*icol*_nrhs, _nlocext);

		};

		_ops += 2*_nlocext*kloc*ndir;

// Store the result in ycorr

		kkk = 0;

		for (irhs=0;irhs<ndir;irhs++) {
			for (kii=0;kii<_n;kii++) {
				ycorr[kkk*_n+kii] = ytemp[irhs*_nlocext+kii];
			};
			kkk++;
		};

// Transform the coefficients

		for (kii=0;kii<kloc;kii++) {
			for (kjj=0;kjj<ndir;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<kloc;kkk++) {
					caux += _rmatr[kkk*_ldrmatr+kii] * pp[kjj*kloc+kkk];
				};
				qloc[kjj*kloc+kii] = caux;
			};
		};

		QrdBlk (kloc, ndir, qloc, kloc, tauloc);

		for (kii=0;kii<ndir*ndir;kii++) aloc[kii] = czero;

		for (kii=0;kii<ndir;kii++) aloc[kii*ndir+kii] = cone;

		for (kii=0;kii<ndir*kextloc;kii++) pp[kii] = czero;

		MvmQBlk (kloc, ndir, ndir, 
					qloc, kloc, tauloc,
					aloc, ndir, pp, kextloc);

		for (igiv=_iter-1;igiv>=0;igiv--) {

			ApplyBlockGivensH (_nrhs, 
								igiv, _givarr+igiv*nrhs_2*4,
								ndir, pp, kextloc, uloc);

		};

// Compute W part

		for (icol=0;icol<ndir/_nrhs;icol++) {

			_qrd.MvmQBlk (_tree, _nlocext, kextloc, _nrhs,
							_p, _ldp, _taup, 
							pp+kextloc*icol*_nrhs, kextloc,
							ytemp+_nlocext*icol*_nrhs, _nlocext);

		};

		_ops += 2*_nlocext*kextloc*ndir;

// Store the result in wcorr

		kkk = 0;

		for (irhs=0;irhs<ndir;irhs++) {
			for (kii=0;kii<_n;kii++) {
				wcorr[kkk*_n+kii] = ytemp[irhs*_nlocext+kii];
			};
			kkk++;
		};

		delete [] ytemp;

// Compute D part in parallel

// Y^H * W

		for (kii=0;kii<ndir*ndir;kii++) tplus[kii] = czero;

		for (kii=0;kii<ndir;kii++) {
			for (kjj=0;kjj<ndir;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<_n;kkk++) {
					caux += wcorr[kjj*_n+kkk] * conj(ycorr[kii*_n+kkk]);
				};
				tplus[kjj*ndir+kii] = caux;
			};
		};

		_ops += _n * ndir * ndir;

		if (_tree.nproc != 1) {

			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														ndir*ndir*2, tplus, aloc);

			for (kii=0;kii<ndir*ndir;kii++) tplus[kii] = aloc[kii];

		};
//
// Check singular values of Y^H * W if necessary
//
		if (true && _tree.myid == _tree.rootcpu && ndir != 0) {

			int ntot = ndir;

			for (kii=0;kii<ndir*ndir;kii++) aloc[kii] = tplus[kii];

			zgesvd_ ("N", "N", &ntot, &ntot, 
						aloc, &ntot, svloc, aloc, &ntot, aloc, &ntot, 
						work, &lwork, rwork, &info);

			if (info != 0) {
				cout << " Error in the Lapack routine ZGESVD" << endl;
				throw " Error in the Lapack routine ZGESVD";
			};

			csvhyst (cout,  "SvYtW", ntot, svloc);
			csvhyst (_fout, "SvYtW", ntot, svloc);

		};
//
// R^{-1}
//
		if (_tree.myid == _tree.rootcpu && ndir != 0) {

// Compute singular value decomposition

			GetRPartQrd (kloc, 0, ndir,
							qloc, kloc, aloc, ndir);

			int ntot = ndir;
//			int lworkloc = 10*ntot;
			int info = 0;

			for (kii=0;kii<ntot*ntot;kii++) dcorr[kii] = aloc[kii];

			zgesvd_ ("A", "A", &ntot, &ntot,
						aloc, &ntot, svloc, 
						uloc, &ntot, vloc, &ntot, 
						work, &lwork, rwork, &info);

			if (info != 0) throw " Error in the Lapack routine ZGESVD";

			_ops += 8 * ntot * ntot * ntot;

			csvhyst (cout,  "SvRcorr", ntot, svloc);
			csvhyst (_fout, "SvRcorr", ntot, svloc);

// Modify singular values

			for (kii=0;kii<ntot;kii++) {
				svloc[kii] = 1.0e0/svloc[kii];
			};

// Compute inverse matrix

			for (kii=0;kii<ntot;kii++) {
				for (kjj=0;kjj<ntot;kjj++) {
					uloc[kjj*ntot+kii] *= svloc[kjj];
				};
			};

			for (kii=0;kii<ntot;kii++) {
				for (kjj=0;kjj<ntot;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<ntot;kkk++) {
						caux += vloc[kii*ntot+kkk] * uloc[kkk*ntot+kjj];
					};
					aloc[kjj*ntot+kii] = conj(caux);
				};
			};

			_ops += ntot * ntot * ntot;

// Check

			if (false) {
				for (kii=0;kii<ntot;kii++) {
					for (kjj=0;kjj<ntot;kjj++) {
						dcmplx caux = czero;
						for (kkk=0;kkk<ntot;kkk++) {
							caux += aloc[kkk*ntot+kii] * dcorr[kjj*ntot+kkk];
						};
						uloc[kjj*ntot+kii] = caux;
					};
				};
				for (kii=0;kii<ntot;kii++) uloc[kii*ntot+kii] -= cone;
				double fnorm = 0.0e0;
				for (kii=0;kii<ntot*ntot;kii++) {
					double auxl = uloc[kii].real ();
					fnorm += auxl*auxl;
					auxl = uloc[kii].imag ();
					fnorm += auxl*auxl;
				};
				cout << " Check of Rinv = " << sqrt(fnorm) << endl;
			};

		} else {
			for (kii=0;kii<ndir*ndir;kii++) aloc[kii] = czero;
		};

		if (_tree.nproc != 1) {

			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														ndir*ndir*2, aloc, dcorr);

			for (kii=0;kii<ndir*ndir;kii++) aloc[kii] = dcorr[kii];

		};
//
// Compute Dcorr matrix
//
		for (kii=0;kii<4*ndir*ndir;kii++) dcorr[kii] = czero;

		for (kii=0;kii<2*ndir;kii++) dcorr[kii*2*ndir+kii] = cmone;

		for (kii=0;kii<ndir;kii++) {
			for (kjj=0;kjj<ndir;kjj++) {
				dcorr[(ndir+kjj)*2*ndir+kii] = aloc[kjj*ndir+kii] + tplus[kjj*ndir+kii];
			};
		};

	};

// Free work data

ExitCorr:;

	delete [] tplus;
	delete [] qloc;
	delete [] tauloc;
	delete [] pp;
	delete [] aloc;
	delete [] rloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix in parallel mode
//========================================================================================
void CSVectorC::BlockGMRESCorr (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix
								const CTree &_tree, CMvmC &_mvm,
								CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
								CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
								const CSVectorC &_x, const CSVectorC &_b, CSVectorC &_ycorr,
								const CSlvParam &_param) { 

	const char *funcname = "BlockGMRESCorr";

// jobitr == 0 - zero out corrector, do not construct the new one
// jobitr == 1 - use available corrector
// jobitr == 2 - use available corrector and add the new one

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);
	dcmplx cmone (-1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timecorr=0.0e0, timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmvcorr=0, nmva=0, nslvl=0, nslvu=0;
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

	if (_param.jobitr == 0) {
		CSVectorC vectdummy;
		_ycorr = vectdummy;
	};
//
// Determine the size of corrector and check
//
	if (_ycorr.nrhs%nrhs != 0) {
		throw " Wrong number of initial vectors on entry to GMRESCorr";
	};

	if (_ycorr.nrhs != 0 && _ycorr.nv != n) {
		throw " Wrong vector dimension on entry to GMRESCorr";
	};

	int ndircorr = _ycorr.nrhs/nrhs;

	int ncorr = ndircorr * nrhs;
	int ldcorr = ncorr * 2;
//
// Init guess to the solution
//
	int npartsloc = _x.GetNparts ();

	CSVectorC sol;

	sol = _x;
//
// Determine extended length
//
	int nitmax = _param.niter;
	if (ndircorr > nitmax) nitmax = ndircorr;

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
	dcmplx *ytw, *rinv, *dcorr, *coef1, *coef2;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi1, *resi;

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

	ytw = new dcmplx [ldcorr*ldcorr];
	if (!ytw) MemoryFail (funcname);
	rinv = new dcmplx [ldcorr*ldcorr];
	if (!rinv) MemoryFail (funcname);
	dcorr = new dcmplx [ldcorr*ldcorr];
	if (!dcorr) MemoryFail (funcname);
	coef1 = new dcmplx [ldcorr*ldcorr];
	if (!coef1) MemoryFail (funcname);
	coef2 = new dcmplx [ldcorr*ldcorr];
	if (!coef2) MemoryFail (funcname);

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
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
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
// Compute corrector matrix
//
	CGSMatrixCS gmtrdummy;

	CSVectorC wcorr(n,npartsloc,nrhs*ndircorr);

	int icol, irhs, kii, kjj, kkk;

	int ncoly=0, ncolp=0;
//
// Check orthogonality before
//
	if (false) {
		for (kii=0;kii<ldcorr*ldcorr;kii++) ytw[kii] = czero;
		for (kii=0;kii<ncorr;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += _ycorr.vect[kii*n+kkk] * conj(_ycorr.vect[kjj*n+kkk]);
				};
				ytw[kjj*ncorr+kii] = caux;
			};
		};
		for (kii=0;kii<ncorr;kii++) ytw[kii*ncorr+kii] -= cone;
		double fnorm = 0.0e0;
		for (kii=0;kii<ncorr*ncorr;kii++) {
			double auxl = ytw[kii].real ();
			fnorm += auxl*auxl;
			auxl = ytw[kii].imag ();
			fnorm += auxl*auxl;
		};
		cout << " Orthogonality of Y = " << sqrt(fnorm) << endl;
	};
//
// Compute Y QR decomposition for the current block
//
	for (icol=0;icol<ndircorr;icol++) {

		dcmplx *pycorr = _ycorr.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<nlocext*nrhs;kii++) p[ncoly*nlocext+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[ncoly*nlocext+irhs*nlocext+kii] = pycorr[ncoly*n+irhs*n+kii];
				};
			};

		} else {

			for (kii=0;kii<nlocext*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pycorr[ncoly*n+irhs*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncoly, ncoly+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncoly*nrhs + nlocext*nrhs_2;
//
// Compute the block to be multiplied
//
		ncoly += nrhs;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncoly;kii++) {
				eloc[ncoly*irhs+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncoly+icol*nrhs+irhs] = cone;
		};

		qrd.MvmQBlk (treeclp, nlocext, ncoly, nrhs,
						p, ldp, taup, 
						eloc, ncoly, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				pycorr[(ncoly-nrhs)*n+irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncoly*nrhs;

	};
//
// Check orthogonality after
//
	if (false) {
		_ycorr.BlockDot(_ycorr,ncorr,ytw);
		for (kii=0;kii<ncorr;kii++) ytw[kii*ncorr+kii] -= cone;
		double fnorm = 0.0e0;
		for (kii=0;kii<ncorr*ncorr;kii++) {
			double auxl = ytw[kii].real ();
			fnorm += auxl*auxl;
			auxl = ytw[kii].imag ();
			fnorm += auxl*auxl;
		};
		cout << " Orthogonality of Y after = " << sqrt(fnorm) << endl;
	};
//
// Multiply by the preconditioned matrix
//
	for (icol=0;icol<ndircorr;icol++) {

		dcmplx *pycorr = _ycorr.GetVect ();
		dcmplx *pwcorr = wcorr.GetVect ();
		dcmplx *py = y.GetVect ();
		dcmplx *pw = w.GetVect ();

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				py[irhs*n+kii] = pycorr[icol*nrhs*n+irhs*n+kii];
			};
		};

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

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				pwcorr[icol*nrhs*n+irhs*n+kii] = pw[irhs*n+kii];
			};
		};

	};
//
// Compute QR decomposition of the result
//
	for (icol=0;icol<ndircorr;icol++) {

		ncolp = icol*nrhs;

		dcmplx *pwcorr = wcorr.GetVect ();

		if (_param.ooctyp == 0) {

			for (kii=0;kii<ldp*nrhs;kii++) p[ldp*ncolp+kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[(irhs+ncolp)*nlocext+kii] = pwcorr[(irhs+ncolp)*n+kii];
				};
			};

		} else {

			for (kii=0;kii<ldp*nrhs;kii++) p[kii] = czero;

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					p[irhs*nlocext+kii] = pwcorr[(irhs+ncolp)*n+kii];
				};
			};

		};

		qrd.UpdateQrd (treeclp, nlocext, ncolp, ncolp+nrhs, p, ldp, taup);

		ops += 2*nlocext*ncolp*nrhs + nlocext*nrhs_2;

// Compute and store W explicitely

		ncolp += nrhs;

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<ncolp;kii++) {
				eloc[ncolp*irhs+kii] = czero;
			};
		};
		for (irhs=0;irhs<nrhs;irhs++) {
			eloc[irhs*ncolp+icol*nrhs+irhs] = cone;
		};

		qrd.MvmQBlk (treeclp, nlocext, ncolp, nrhs,
						p, ldp, taup, 
						eloc, ncolp, 
						yp, nlocext);

		for (irhs=0;irhs<nrhs;irhs++) {
			for (kii=0;kii<n;kii++) {
				pwcorr[(ncolp-nrhs)*n+irhs*n+kii] = yp[irhs*nlocext+kii];
			};
		};

		ops += 2*nlocext*ncolp*nrhs;

	};
//
// Check orthogonality of W
//
	if (false) {
		for (kii=0;kii<ldcorr*ldcorr;kii++) ytw[kii] = czero;
		for (kii=0;kii<ncorr;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += wcorr.vect[kii*n+kkk] * conj(wcorr.vect[kjj*n+kkk]);
				};
				ytw[kjj*ncorr+kii] = caux;
			};
		};
		for (kii=0;kii<ncorr;kii++) ytw[kii*ncorr+kii] -= cone;
		double fnorm = 0.0e0;
		for (kii=0;kii<ncorr*ncorr;kii++) {
			double auxl = ytw[kii].real ();
			fnorm += auxl*auxl;
			auxl = ytw[kii].imag ();
			fnorm += auxl*auxl;
		};
		cout << " Orthogonality of W after = " << sqrt(fnorm) << endl;
	};
//
// Compute corrector matrix in parallel
//
// Y^H * W
//
	dcmplx *pycorr = _ycorr.GetVect ();
	dcmplx *pwcorr = wcorr.GetVect ();

	for (kii=0;kii<ncorr*ncorr;kii++) ytw[kii] = czero;

	for (kii=0;kii<ncorr;kii++) {
		for (kjj=0;kjj<ncorr;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux += pwcorr[kjj*n+kkk] * conj(pycorr[kii*n+kkk]);
			};
			ytw[kjj*ncorr+kii] = caux;
		};
	};

	ops += n * ncorr * ncorr;

	if (treeclp.nproc != 1) {

		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													ncorr*ncorr*2, ytw, rinv);

		for (kii=0;kii<ncorr*ncorr;kii++) ytw[kii] = rinv[kii];

	};
//
// Check singular values of Y^H * W if necessary
//
	if (true && _tree.myid == _tree.rootcpu && ncorr != 0) {

		int ntot = ncorr;
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

		for (kii=0;kii<ncorr*ncorr;kii++) rinv[kii] = ytw[kii];

		zgesvd_ ("N", "N", &ntot, &ntot, 
					rinv, &ntot, svloc, rinv, &ntot, rinv, &ntot, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvYtW", ntot, svloc);
		csvhyst (_fout, "SvYtW", ntot, svloc);

		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};
//
// R^{-1}
//
	qrd.GetRPartQrd (treeclp, 0, 0, ncorr,
						rinv, ncorr);

	if (treeclp.myid == treeclp.rootcpu && ncorr != 0) {

// Compute singular value decomposition

		int ntot = ncorr;
		int lworkloc = 10*ntot;
		int info = 0;

		double *svloc, *rwork;
		dcmplx *workloc;

		dcmplx *uloc, *vloc;

		uloc = new dcmplx [ntot*ntot];
		if (!uloc) MemoryFail (funcname);
		vloc = new dcmplx [ntot*ntot];
		if (!vloc) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		for (kii=0;kii<ntot*ntot;kii++) dcorr[kii] = rinv[kii];

		zgesvd_ ("A", "A", &ntot, &ntot,
					rinv, &ntot, svloc, 
					uloc, &ntot, vloc, &ntot, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) throw " Error in the Lapack routine ZGESVD";

		ops += 8 * ntot * ntot * ntot;

		csvhyst (cout,  "SvRcorr", ntot, svloc);
		csvhyst (_fout, "SvRcorr", ntot, svloc);

// Modify singular values

		for (kii=0;kii<ntot;kii++) {
			svloc[kii] = 1.0e0/svloc[kii];
		};

// Compute inverse matrix

		for (kii=0;kii<ntot;kii++) {
			for (kjj=0;kjj<ntot;kjj++) {
				uloc[kjj*ntot+kii] *= svloc[kjj];
			};
		};

		for (kii=0;kii<ntot*ntot;kii++) rinv[kii] = czero;

		for (kii=0;kii<ntot;kii++) {
			for (kjj=0;kjj<ntot;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<ntot;kkk++) {
					caux += vloc[kii*ntot+kkk] * uloc[kkk*ntot+kjj];
				};
				rinv[kjj*ntot+kii] = conj(caux);
			};
		};

		ops += ntot * ntot * ntot;

		delete [] uloc;
		delete [] vloc;
		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

// Check

		if (false) {
			for (kii=0;kii<ntot;kii++) {
				for (kjj=0;kjj<ntot;kjj++) {
					dcmplx caux = czero;
					for (kkk=0;kkk<ntot;kkk++) {
						caux += rinv[kkk*ntot+kii] * dcorr[kjj*ntot+kkk];
					};
					rmatr[kjj*ntot+kii] = caux;
				};
			};
			for (kii=0;kii<ncorr;kii++) rmatr[kii*ncorr+kii] -= cone;
			double fnorm = 0.0e0;
			for (kii=0;kii<ncorr*ncorr;kii++) {
				double auxl = rmatr[kii].real ();
				fnorm += auxl*auxl;
				auxl = rmatr[kii].imag ();
				fnorm += auxl*auxl;
			};
			cout << " Check of Rinv = " << sqrt(fnorm) << endl;
		};

	} else {
		for (kii=0;kii<ncorr*ncorr;kii++) rinv[kii] = czero;
	};

	if (treeclp.nproc != 1) {

		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													ncorr*ncorr*2, rinv, dcorr);

		for (kii=0;kii<ncorr*ncorr;kii++) rinv[kii] = dcorr[kii];

	};
//
// Compute Dcorr matrix
//
	for (kii=0;kii<ldcorr*ldcorr;kii++) dcorr[kii] = czero;

	for (kii=0;kii<ldcorr;kii++) dcorr[kii*ldcorr+kii] = cmone;

	for (kii=0;kii<ncorr;kii++) {
		for (kjj=0;kjj<ncorr;kjj++) {
			dcorr[(ncorr+kjj)*ldcorr+kii] = rinv[kjj*ncorr+kii] + ytw[kjj*ncorr+kii];
		};
	};
	_fout << " Y part in Vect" << _ycorr << endl;
	_fout << " W part in Vect" << wcorr << endl;
	OutArr(_fout," D part in Vect",ldcorr*ldcorr,dcorr);
//
// Check matrix relations
//
	if (false) {
//
// Check that A(LU)^{-1} Y = W R
//
		CSVectorC zcheck(n,npartsloc,nrhs*ndircorr);

		for (icol=0;icol<ndircorr;icol++) {

			dcmplx *pycorr = _ycorr.GetVect ();
			dcmplx *pzchk  = zcheck.GetVect ();
			dcmplx *py = y.GetVect ();
			dcmplx *pw = w.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					py[irhs*n+kii] = pycorr[icol*nrhs*n+irhs*n+kii];
				};
			};

			_gmtrl.SolveL (_tree,_mvm,y,w);

			_gmtru.SolveU (_tree,_mvm,w,y);

			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,y,w);

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					pzchk[icol*nrhs*n+irhs*n+kii] = pw[irhs*n+kii];
				};
			};

		};

		qrd.GetRPartQrd (treeclp, 0, 0, ncorr,
						ytw, ncorr);

/*
		double fnorm = 0.0e0;
		for (kii=0;kii<n;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<ncorr;kkk++) {
					caux += wcorr.vect[kkk*n+kii] * ytw[kjj*ncorr+kkk];
				};
				caux -= zcheck.vect[kjj*n+kii];
				double auxl = caux.real ();
				fnorm += auxl*auxl;
				auxl = caux.imag ();
				fnorm += auxl*auxl;
			};
		}; */
		for (kii=0;kii<ncorr*ncorr;kii++) ytw[kii] = -ytw[kii];

		zcheck.BlockDaxpy (wcorr, ncorr, ytw);

		double fnorm = 0.0e0;
		for (kii=0;kii<n;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				caux = zcheck.vect[kjj*n+kii];
				double auxl = caux.real ();
				fnorm += auxl*auxl;
				auxl = caux.imag ();
				fnorm += auxl*auxl;
			};
		};
		cout << " || AY-WR || = " << sqrt(fnorm) << endl;
//
// Check that K W = Y R^{-1} directly
//
		zcheck = wcorr;

		zcheck.MvCorr (treeclp, _ycorr, wcorr,
						ldcorr, dcorr, 
						ldcorr, coef1, coef2);

		fnorm = 0.0e0;
		for (kii=0;kii<n;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<ncorr;kkk++) {
					caux += _ycorr.vect[kkk*n+kii] * rinv[kjj*ncorr+kkk];
				};
				caux -= zcheck.vect[kjj*n+kii];
				double auxl = caux.real ();
				fnorm += auxl*auxl;
				auxl = caux.imag ();
				fnorm += auxl*auxl;
			};
		};
		cout << " || KW-YR^{-1} direct || = " << sqrt(fnorm) << endl;

//
// Check that K W = Y R^{-1}
//
		for (icol=0;icol<ndircorr;icol++) {

			dcmplx *pwcorr = wcorr.GetVect ();
			dcmplx *pzchk  = zcheck.GetVect ();
			dcmplx *py = y.GetVect ();
//			dcmplx *pw = w.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					py[irhs*n+kii] = pwcorr[icol*nrhs*n+irhs*n+kii];
				};
			};

			y.MvCorr (treeclp, _ycorr, wcorr,
						ldcorr, dcorr, 
						ldcorr, coef1, coef2);

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					pzchk[icol*nrhs*n+irhs*n+kii] = py[irhs*n+kii];
				};
			};

		};

		fnorm = 0.0e0;
		for (kii=0;kii<n;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<ncorr;kkk++) {
					caux += _ycorr.vect[kkk*n+kii] * rinv[kjj*ncorr+kkk];
				};
				caux -= zcheck.vect[kjj*n+kii];
				double auxl = caux.real ();
				fnorm += auxl*auxl;
				auxl = caux.imag ();
				fnorm += auxl*auxl;
			};
		};
		cout << " || KW-YR^{-1} || = " << sqrt(fnorm) << endl;
//
// Check that W^H A(LU)^{-1} K W = I
//
		for (icol=0;icol<ndircorr;icol++) {

			dcmplx *pwcorr = wcorr.GetVect ();
			dcmplx *pzchk  = zcheck.GetVect ();
			dcmplx *py = y.GetVect ();
			dcmplx *pw = w.GetVect ();

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					py[irhs*n+kii] = pwcorr[icol*nrhs*n+irhs*n+kii];
				};
			};

			y.MvCorr (treeclp, _ycorr, wcorr,
						ldcorr, dcorr, 
						ldcorr, coef1, coef2);

			_gmtrl.SolveL (_tree,_mvm,y,w);

			_gmtru.SolveU (_tree,_mvm,w,y);

			gmtrdummy.MvmA (_tree,_mvm,_gmtral,_gmtrau,y,w);

			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<n;kii++) {
					pzchk[icol*nrhs*n+irhs*n+kii] = pw[irhs*n+kii];
				};
			};

		};

		for (kii=0;kii<ncorr;kii++) {
			for (kjj=0;kjj<ncorr;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<n;kkk++) {
					caux += zcheck.vect[kjj*n+kkk] * conj(wcorr.vect[kii*n+kkk]);
				};
				ytw[kjj*ncorr+kii] = caux;
			};
		};

		for (kii=0;kii<ncorr;kii++) ytw[kii*ncorr+kii] -= cone;

		fnorm = 0.0e0;

		for (kii=0;kii<ncorr*ncorr;kii++) {
			double auxl = ytw[kii].real ();
			fnorm += auxl*auxl;
			auxl = ytw[kii].imag ();
			fnorm += auxl*auxl;
		};

		cout << " || W^H A(LU)^{-1} K W - I || = " << sqrt(fnorm) << endl;

	};
//
// Compute the array of right hand side norms
//
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
		resi1[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi1, resi);
	} else {
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
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int iconv;
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
// Multiply by the corrected preconditioned matrix
//
		time2 = clock ();

		y.MvCorr (treeclp, _ycorr, wcorr,
					ldcorr, dcorr, 
					ldcorr, coef1, coef2);

		ops += n*ncorr*4*nrhs;

		time3 = clock ();

		timecorr += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmvcorr++;

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

	y.MvCorr (treeclp, _ycorr, wcorr,
				ldcorr, dcorr, 
				ldcorr, coef1, coef2);

	ops += n*ncorr*4*nrhs;

	time3 = clock ();

	timecorr += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmvcorr++;

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

		dcmplx *aloc;
		double *svloc, *rwork;
		dcmplx *workloc;

		aloc = new dcmplx [ntot*ntot];
		if (!aloc) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		for (kii=0;kii<ntot;kii++) {
			for (kjj=0;kjj<ntot;kjj++) {
				aloc[kjj*ntot+kii] = rmatr[kjj*ldrmatr+kii];
			};
		};

		zgesvd_ ("N", "N", &ntot, &ntot, 
					aloc, &ntot, svloc, aloc, &ntot, aloc, &ntot, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] aloc;
		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

//
// Compute/update corrector subspace if necessary
//
	if (_param.jobitr == 2) {

		CSVectorC vectdummy;

		vectdummy.BlockCorrector (_fout, treeclp, 
									n, nlocext, nrhs, iter+1, 
									qrd, p, ldp, taup,
									ldrmatr, rmatr, givarr,
									_ycorr,
									_param, ops);

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
		cout  << "     NMCor = " << nmvcorr  << " TimeMCor = " << timecorr  << " sec." <<
				 " MnTimeMCr = " << timecorr / (double) nmvcorr << " sec." << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

		_fout << " GMRES statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMCor = " << nmvcorr  << " TimeMCor = " << timecorr  << " sec." <<
				 " MnTimeMCr = " << timecorr / (double) nmvcorr << " sec." << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] ytw;
	delete [] rinv;
	delete [] dcorr;
	delete [] coef1;
	delete [] coef2;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;

// Return the solution

	*this = sol;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix in parallel mode
//========================================================================================
int CSVectorC::BlockGMRESCorr (ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix
								const CTree &_tree, CMvmC &_mvm,
								CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
								CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
								const CSVectorC &_x, const CSVectorC &_b, 
								CCorrectorC &_corr,
								const CSlvParam &_param) { 

	const char *funcname = "BlockGMRESCorr";

// jobitr == 0 - zero out corrector, do not construct the new one
// jobitr == 1 - use available corrector
// jobitr == 2 - use available corrector and add the new one

	int iconvglob = 0;

	dcmplx czero ( 0.0e0,0.0e0);
	dcmplx cone  ( 1.0e0,0.0e0);

	int nza, nzlu;
	double ops, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timecorr=0.0e0, timemva=0.0e0, timeslvl=0.0e0, timeslvu=0.0e0;
	int nmvcorr=0, nmva=0, nslvl=0, nslvu=0;
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

	if (_param.jobitr == 0) {
		CCorrectorC corrdummy (_corr.ncorrmax);
		_corr = corrdummy;
	};
//
// Determine the maximal size of corrector data, the total size of correctors and check
//
	int ncorrtot=0;
	int ncorrmax=0;

	if (_corr.ncorr != 0 && _corr.ldcorrarr[0] != n) {
		throw " Wrong vector dimension on entry to GMRESCorr";
	};

	for (ilist=0;ilist<_corr.ncorr;ilist++) {
		if (_corr.corrsizes[ilist] > ncorrmax) ncorrmax = _corr.corrsizes[ilist];
		ncorrtot += _corr.corrsizes[ilist];
		if (_corr.corrsizes[ilist]%nrhs != 0) {
			throw " Wrong number of corrector vectors on entry to GMRESCorr";
		};
	};

	int ldcorr = ncorrmax * 2;
//
// Init guess to the solution
//
	int npartsloc = _x.GetNparts ();

	CSVectorC sol;

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
	dcmplx *coef1, *coef2;
	dcmplx *givarr, *giv, *giv1, *taugiv;
	dcmplx *bloc, *eloc, *rmatr;
	double *resi0, *resi1, *resi;

	int nrhs_2 = nrhs*nrhs;
	int nrhs_3 = nrhs*nrhs_2;

	int ldp = nlocext;

	if (_param.oocitr == 0) {
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

	coef1 = new dcmplx [ldcorr*ldcorr];
	if (!coef1) MemoryFail (funcname);
	coef2 = new dcmplx [ldcorr*ldcorr];
	if (!coef2) MemoryFail (funcname);

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
	resi1 = new double [nrhs];
	if (!resi1) MemoryFail (funcname);
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
	if (_param.oocitr > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"Gmres_",_tree.myid,"_");

		qrd.SetupFiles (_param.oocitr, strbuff, ldp);

	};
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
		resi1[irhs] = aux;
	};

	ops += 2*nrhs_2;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													nrhs, resi1, resi);
	} else {
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
// Init R matrix
//
	for (kii=0;kii<ldrmatr*ldrmatr;kii++) rmatr[kii] = czero;
//
// Main iterative cycle
//
	int iconv;
	int iter = 0;
	int ncolp, kjj;

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
// Multiply by the corrected preconditioned matrix
//
		time2 = clock ();

		y.MvCorr (treeclp, _corr, 
					ldcorr, coef1, coef2);

		ops += n*ncorrtot*4*nrhs;

		time3 = clock ();

		timecorr += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		nmvcorr++;

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

		if (_param.oocitr == 0) {

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

		if (iconvtot == nrhs) {
			iconvglob = 1;
			goto exitgmr;
		};

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

	y.MvCorr (treeclp, _corr,
				ldcorr, coef1, coef2);

	ops += n*ncorrtot*4*nrhs;

	time3 = clock ();

	timecorr += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

	nmvcorr++;

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

		dcmplx *aloc;
		double *svloc, *rwork;
		dcmplx *workloc;

		aloc = new dcmplx [ntot*ntot];
		if (!aloc) MemoryFail (funcname);
		svloc = new double [ntot];
		if (!svloc) MemoryFail (funcname);
		workloc = new dcmplx [lworkloc];
		if (!workloc) MemoryFail (funcname);
		rwork = new double [lworkloc];
		if (!rwork) MemoryFail (funcname);

		for (kii=0;kii<ntot;kii++) {
			for (kjj=0;kjj<ntot;kjj++) {
				aloc[kjj*ntot+kii] = rmatr[kjj*ldrmatr+kii];
			};
		};

		zgesvd_ ("N", "N", &ntot, &ntot, 
					aloc, &ntot, svloc, aloc, &ntot, aloc, &ntot, 
					workloc, &lworkloc, rwork, &info);

		if (info != 0) {
			cout << " Error in the Lapack routine ZGESVD" << endl;
			throw " Error in the Lapack routine ZGESVD";
		};

		csvhyst (cout,  "SvGmres", ntot, svloc);
		csvhyst (_fout, "SvGmres", ntot, svloc);

		delete [] aloc;
		delete [] svloc;
		delete [] workloc;
		delete [] rwork;

	};

//
// Compute/update corrector subspace if necessary
//
	if (_param.jobitr == 2) {

		CSVectorC vectdummy;

		vectdummy.BlockCorrector (_fout, treeclp, 
									n, nlocext, nrhs, iter+1, 
									qrd, p, ldp, taup,
									ldrmatr, rmatr, givarr,
									_corr,
									_param, ops);

	};
//
// Close files if necessary
//
	if (_param.oocitr > 0) {

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
		cout  << "     NMCor = " << nmvcorr  << " TimeMCor = " << timecorr  << " sec." <<
				 " MnTimeMCr = " << timecorr / (double) nmvcorr << " sec." << endl;
		cout  << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		cout  << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		cout  << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

		_fout << " GMRES statistics: " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		_fout << " Mv-Slv functions statistics: " << endl;
		_fout << "     NMCor = " << nmvcorr  << " TimeMCor = " << timecorr  << " sec." <<
				 " MnTimeMCr = " << timecorr / (double) nmvcorr << " sec." << endl;
		_fout << "     NMvmA = " << nmva  << " TimeMvmA = " << timemva  << " sec." <<
				 " MnTimeMvA = " << timemva / (double) nmva << " sec." << endl;
		_fout << "     NSlvL = " << nslvl << " TimeSlvL = " << timeslvl << " sec." <<
				 " MnTimeSlL = " << timeslvl / (double)nslvl<< " sec." << endl;
		_fout << "     NSlvU = " << nslvu << " TimeSlvU = " << timeslvu << " sec." <<
				 " MnTimeSlU = " << timeslvu / (double)nslvu<< " sec." << endl;

	};

// Free work arrays

	delete [] p;
	delete [] yp;
	delete [] taup;
	delete [] coef1;
	delete [] coef2;
	delete [] givarr;
	delete [] giv;
	delete [] giv1;
	delete [] taugiv;
	delete [] bloc;
	delete [] eloc;
	delete [] rmatr;
	delete [] resi0;
	delete [] resi1;
	delete [] resi;

// Return the solution

	*this = sol;

	return iconvglob;

};
