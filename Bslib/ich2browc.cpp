//------------------------------------------------------------------------------------------------
// File: ich2browc.cpp
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

#include "tree.h"
#include "fct.h"
#include "fctdiag.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "slvparam.h"
#include "mvm.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CFctDiagC: Init scaling for the current block row
//========================================================================================
void CFctDiagC::Ich2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init scaling for the current block row
									const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau,
									CFctC &_fct) {

	const char *funcname = "Ich2ScaleBlkRow";

	dcmplx cone (1.0e0,0.0e0);
	dcmplx czero (0.0e0,0.0e0);

// Determine the size of the arrays to be stored

	int isup, isupbeg, isupend, nlistl, blai, kii;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	int nz=0;
	int nlocblk=0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		nz += blai*blai;
		nlocblk += blai;
	};

// Allocate and register new data

	dcmplx *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc;
	double *dpivloc;

	sclptloc = new dcmplx [nlocblk];
	if (!sclptloc) MemoryFail (funcname);
	scllloc = new dcmplx [nz];
	if (!scllloc) MemoryFail (funcname);
	scluloc = new dcmplx [nz];
	if (!scluloc) MemoryFail (funcname);
	sclinvlloc = new dcmplx [nz];
	if (!sclinvlloc) MemoryFail (funcname);
	sclinvuloc = new dcmplx [nz];
	if (!sclinvuloc) MemoryFail (funcname);
	diaggloc = new dcmplx [nz];
	if (!diaggloc) MemoryFail (funcname);
	dpivloc = new double [nlocblk];
	if (!dpivloc) MemoryFail (funcname);

	delete [] sclpt  [_iblkrow];
	delete [] scll   [_iblkrow];
	delete [] sclu   [_iblkrow];
	delete [] sclinvl[_iblkrow];
	delete [] sclinvu[_iblkrow];
	delete [] diagg  [_iblkrow];
	delete [] dpiv   [_iblkrow];

	diasize[_iblkrow] = nz;
	sclpt  [_iblkrow] = sclptloc;
	scll   [_iblkrow] = scllloc;
	sclu   [_iblkrow] = scluloc;
	sclinvl[_iblkrow] = sclinvlloc;
	sclinvu[_iblkrow] = sclinvuloc;
	diagg  [_iblkrow] = diaggloc;
	dpiv   [_iblkrow] = dpivloc;

// Init point scaling data for the current block row if necessary

	int j, kjj, kiib;
	int ibs, ibsv, blai_2;

	double *peig;

	dcmplx *paloc, *pa, *pscl;

	dcmplx caux;
	double daux;

	if (_sclpttype == 1) {

		nz = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {

			blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
			blai_2 = blai*blai;

// Search diagonal supernode

			j = _mtrau.ia[isup-isupbeg];
			ibs = _mtrau.bsa[j];

			ibsv = _gmtra.sprndr[isup]-_gmtra.sprndr[isupbeg];

			pa = _mtrau.a+ibs;

			for (kii=0;kii<blai;kii++) {
				caux = pa[kii*blai+kii];
				daux = abs(caux);
				dpivloc[nz+kii] = daux;
				daux = sqrt(daux);
				daux = 1.0e0 / daux;
				sclptloc[nz+kii].x = daux;
				sclptloc[nz+kii].y = 0.0e0;
			};

			nz += blai;

		};

	} else {
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = cone;
		for (kii=0;kii<nlocblk;kii++) dpivloc[kii] = 1.0e0;
	};

// Update hystogram

	hyst2.UpdateHystogram (nlocblk,dpivloc);

// Init scaling data for the current block row

	int npt=0;

	nz = 0;
//	ofstream fout ("CheckScl.dat");

	for (isup=isupbeg;isup<=isupend;isup++) {

//		fout << " isup = " << isup << endl;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		blai_2 = blai*blai;

// Assign diagonal supernode

		j = _mtrau.ia[isup-isupbeg];
		ibs = _mtrau.bsa[j];

		ibsv = _gmtra.sprndr[isup]-_gmtra.sprndr[isupbeg];

		paloc = _fct.aloc;
		pa = _mtrau.a+ibs;
//		for (kii=0;kii<blai_2;kii++) _aloc[kii] = a[ibs+kii];
		for (kii=0;kii<blai_2;kii++) *paloc++ = *pa++;

// Perform point scaling

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
				caux = sclptloc[npt+kii] * sclptloc[npt+kjj];
				_fct.aloc[kjj*blai+kii] *= caux;
			};
		};

// Compute its eigenvalue decomposition

		int blal = blai;
		int info = 0;
//		OutArr (fout," A bef zheev ",blal*blal,_fct.aloc);

//		int MyZheev (char *, char *, int *, dcmplx *, 
//				int *, double *, dcmplx *, int *, double *, int *);
//		MyZheev ("V", "U", &blal, 
//					_fct.aloc, &blal, _fct.eig, 
//					_fct.work, &_fct.lwork, _fct.rwork, &info);
		zheev_ ("V", "U", &blal, 
					_fct.aloc, &blal, _fct.eig, 
					_fct.work, &_fct.lwork, _fct.rwork, &info);

		if (info != 0) throw " Error in the Lapack routine ZHEEV";

		_fct.ops += 8 * blai_2*blai;

//		OutArr (fout," Eig aft zheev ",blal,_fct.eig);
//		OutArr (fout," A aft zheev ",blal*blal,_fct.aloc);

// Modify the pivots

//		fout << " Sclmin = " << _sclmin << endl;
		for (kii=0;kii<blai;kii++) {
			if (_fct.eig[kii] < _sclmin) _fct.eig[kii] = _sclmin;
		};
//		OutArr (fout," Modified Eig aft zheev ",blal,_fct.eig);

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			daux = _fct.eig[kii];
			daux = sqrt(daux);
			_fct.eig[kii] = daux;
			dpivloc[ibsv+kii] = _fct.eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			peig = _fct.eig;
			paloc = _fct.aloc+kii;
			pscl = sclinvlloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvl[ibsd+kjj*blai+kii] = _uloc[kjj*blai+kii] * _eig[kjj];
				*pscl = *paloc * *peig++;
				paloc += blai;
				pscl += blai;
			};
		};

		for (kii=0;kii<blai_2;kii++) sclinvuloc[nz+kii] = conj(sclinvlloc[nz+kii]);

		for (kii=0;kii<blai;kii++) {
			daux = _fct.eig[kii];
			daux = 1.0e0 / daux;
			_fct.eig[kii] = daux;
		};

		peig = _fct.eig;
		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			paloc = _fct.aloc+kiib;
			pscl = scllloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scll[ibsd+kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//				*pscl = conj(*paloc++) * *peig;
				dcmplx temp = conj(*paloc++);
				*pscl = temp * *peig;
				pscl += blai;
			};
			peig++;
			kiib += blai;
		};

		for (kii=0;kii<blai_2;kii++) scluloc[nz+kii] = conj(scllloc[nz+kii]);

		_fct.ops += 4*blai_2 + 2*blai;

		nz += blai*blai;
		npt += blai;

	};

// Init diagg

	for (kii=0;kii<nz;kii++) diaggloc[kii] = czero;

// Update hystogram

	hyst.UpdateHystogram (nlocblk,dpivloc);

// Store data to the disk if necessary

	if (nfiles > 0) {

		int ibsloc = ibsfile[ifile];

		blk2file[5*_iblkrow] = ifile;
		blk2bs  [5*_iblkrow] = ibsloc;

		FPut (files[ifile], nz, scllloc, ibsloc);

		ibsloc += nz;

		blk2file[5*_iblkrow+1] = ifile;
		blk2bs  [5*_iblkrow+1] = ibsloc;

		FPut (files[ifile], nz, scluloc, ibsloc);

		ibsloc += nz;

		blk2file[5*_iblkrow+2] = ifile;
		blk2bs  [5*_iblkrow+2] = ibsloc;

		FPut (files[ifile], nz, sclinvlloc, ibsloc);

		ibsloc += nz;

		blk2file[5*_iblkrow+3] = ifile;
		blk2bs  [5*_iblkrow+3] = ibsloc;

		FPut (files[ifile], nz, sclinvuloc, ibsloc);

		ibsloc += nz;

		blk2file[5*_iblkrow+4] = ifile;
		blk2bs  [5*_iblkrow+4] = ibsloc;

		FPut (files[ifile], nz, diaggloc, ibsloc);

		ibsloc += nz;

		ibsfile[ifile] = ibsloc;

		ifile++;

		if (ifile >= nfiles) ifile = 0;

		delete [] scllloc;
		delete [] scluloc;
		delete [] sclinvlloc;
		delete [] sclinvuloc;
		delete [] diaggloc;

		scllloc = new dcmplx [0];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new dcmplx [0];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new dcmplx [0];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new dcmplx [0];
		if (!sclinvuloc) MemoryFail (funcname);
		diaggloc = new dcmplx [0];
		if (!diaggloc) MemoryFail (funcname);

		scll   [_iblkrow] = scllloc;
		sclu   [_iblkrow] = scluloc;
		sclinvl[_iblkrow] = sclinvlloc;
		sclinvu[_iblkrow] = sclinvuloc;
		diagg  [_iblkrow] = diaggloc;

	};

};

/*
//========================================================================================
// Author: Kharchenko S.A.
// CFctDiagC: Store point scaling data in the vector array
void CFctDiagC::StorePointScaling (CGSMatrixCS &_mtrl) const { // Store point scaling data in the vector array

	const char *funcname = "StorePointScaling";

// Determine the size of the arrays to be stored

	int nblkloc = 0;
	int nloc = 0;

	int i, ni, kii;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			nloc += (_mtrl.bl2ndc[i+1]-_mtrl.bl2ndc[i]);
			nblkloc++;
		};
	};

// Allocate temporary SVectorC

	CSVectorC temp (nloc, nblkloc, 1);

// Write scaling into SVectorC

	nblkloc = 0;
	nloc = 0;

	temp.blksv[0] = 0;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			ni = _mtrl.bl2ndc[i+1]-_mtrl.bl2ndc[i];
			for (kii=0;kii<ni;kii++) temp.vect[nloc+kii] = sclpt[i][kii];
			temp.blksv[nblkloc+1] = temp.blksv[nblkloc] + ni;
			temp.inda[nblkloc] = _mtrl.bl2ndc[i];
			nloc += ni;
			nblkloc++;
		};
	};

// Store

	_mtrl.vector1 = temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Perform point scaling of the block row
//========================================================================================
void CSMatrixCS::PointScaling (int *_sprndr, int *_ibsscl, dcmplx *_scl) { // Perform point scaling of the block row

	const char *funcname = "PointScaling";

// Scan all local supernodes

	int ilist, isup, jsup, j, blai, blaj, ibs, ibsv, jbsv, kii, kjj;
	dcmplx caux;

	for (ilist=0;ilist<nlist;ilist++) {
		isup = list[ilist];
		blai = _sprndr[isup+1]-_sprndr[isup];
		ibsv = _ibsscl[isup];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jsup = ja[j];
			blaj = _sprndr[jsup+1]-_sprndr[jsup];
			jbsv = _ibsscl[jsup];
			if (ibsv < 0 || jbsv < 0) {
				throw " PointScaling: required supernode is not found in scaling array";
			};
			ibs = bsa[j];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					caux = _scl[ibsv+kii] * _scl[jbsv+kjj];
					a[ibs+kjj*blai+kii] *= caux;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Perform supernode scaling of the block row
//========================================================================================
void CSMatrixCS::BlockScaling (int *_sprndr, int *_ibsscl, dcmplx *_scll, dcmplx *_sclu) { // Perform supernode scaling of the block row

	const char *funcname = "BlockScaling";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work array

	dcmplx *aloc;
	dcmplx *vloc;

	aloc = new dcmplx [blamx*blamx];
	if (!aloc) MemoryFail (funcname);
	vloc = new dcmplx [blamx*blamx];
	if (!vloc) MemoryFail (funcname);

// Scan all local supernodes

	int ilist, isup, jsup, j, blai, blaj, ibs, ibsv, jbsv, kii, kjj;
	int kjjb, kkk, blaij;
	dcmplx aux;
	dcmplx *paloc, *pvloc, *pscloc;

	for (ilist=0;ilist<nlist;ilist++) {
		isup = list[ilist];
		blai = _sprndr[isup+1]-_sprndr[isup];
		ibsv = _ibsscl[isup];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jsup = ja[j];
			blaj = _sprndr[jsup+1]-_sprndr[jsup];
			jbsv = _ibsscl[jsup];
			if (ibsv < 0 || jbsv < 0) {
				throw " BlockScaling: required supernode is not found in scaling array";
			};
			ibs = bsa[j];
			blaij = blai*blaj;
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					paloc = a+ibs+kjjb;
					pscloc = _scll+ibsv+kii;
					aux = czero;
					for (kkk=0;kkk<blai;kkk++) {
//						aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//						aux += *paloc++ * *pscloc;
						aux.PlEq (*paloc++ * *pscloc);
						pscloc += blai;
					};
					aloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					aux = czero;
					paloc = aloc+kii;
					pscloc = _sclu+jbsv+kjj;
					for (kkk=0;kkk<blaj;kkk++) {
//						aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//						aux += *paloc * *pscloc;
						aux.PlEq (*paloc * *pscloc);
						pscloc += blaj;
						paloc += blai;
					};
					vloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			paloc = a+ibs;
			pvloc = vloc;
			for (kii=0;kii<blaij;kii++) *paloc++ = *pvloc++;
		};
	};

// Free work array

	delete [] aloc;
	delete [] vloc;

};

*/
// Author: Kharchenko S.A.
// CFctC: Init current block row
//========================================================================================
void CFctC::Ich2InitBlkRow (int _iblkrow, // Init current block row
							double _dshift, double _dshiftim,
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, CSMatrixCS &_mtru) {

//	const char *funcname = "Ich2InitBlkRow";

	dcmplx czero (0.0e0,0.0e0);

// Check the data on entry

	if (_iblkrow >= _gmtra.nblksc || _iblkrow < 0) {
		cout << " Wrong block row number " << endl;
		throw " Wrong block row number ";
	};

	int ilist, isup, isupbeg, isupend, nlistl;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	if (nlistl != _mtrau.nlist) {
		cout << " Wrong number of supernode rows in the block row" << endl;
		throw " Wrong number of supernode rows in the block row";
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = _mtrau.list[ilist];
		if (isup < isupbeg || isup > isupend) {
			cout << " Wrong supernode row in the block row" << endl;
			throw " Wrong supernode row in the block row";
		};
	};

// Copy the data

	_mtru = _mtrau;

	_mtru.nsupr = _mtru.nlist;
	_mtru.nsupc = _mtru.nlist;

// Scale 

	int ibs, ibssc, jbssc, blai, blaj, blaij, blaj_2;
	int kii, kjj, kkk, kjjb, ibsv, jbsv;

	dcmplx *paloc, *pscloc, *pvloc;
	dcmplx aux;

//	ofstream fout ("CheckIni.dat");
	for (ilist=0;ilist<nlistl;ilist++) {

		isup = _mtru.list[ilist];
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];

		for (int j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {

			int jsup = _mtru.ja[j];

			ibs = _mtru.bsa[j];
			ibssc = bsdia[isup];
			jbssc = bsdia[jsup];
			blaj = _gmtra.sprndr[jsup+1]-_gmtra.sprndr[jsup];
			ibsv = adrfmask[isup];
			jbsv = adrfmask[jsup];

			blaij = blai * blaj;
			blaj_2 = blaj*blaj;
//			OutArr (fout," A bef scale ",blaij,_mtru.a+ibs);

// Local point scaling (U part)

			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					aux = _mtru.a[ibs+kjj*blai+kii] * sclpt[ibsv+kii];
					vloc[kjj*blai+kii] = aux * sclpt[jbsv+kjj];
				};
			};
//			OutArr (fout," A aft point scale ",blaij,_mtru.a+ibs);

// Local supernode scaling (U part)

			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					paloc = vloc+kjjb;
					pscloc = scll+ibssc+kii;
					aux = czero;
					for (kkk=0;kkk<blai;kkk++) {
//						aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//						aux += *paloc++ * *pscloc;
						aux.PlEq (*paloc++ * *pscloc);
						pscloc += blai;
					};
					aloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					aux = czero;
					paloc = aloc+kii;
					pscloc = sclu+jbssc+kjj;
					for (kkk=0;kkk<blaj;kkk++) {
//						aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//						aux += *paloc * *pscloc;
						aux.PlEq (*paloc * *pscloc);
						pscloc += blaj;
						paloc += blai;
					};
					vloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			paloc = _mtru.a+ibs;
			pvloc = vloc;
			for (kii=0;kii<blaij;kii++) *paloc++ = *pvloc++;
//			OutArr (fout," A aft block scale ",blaij,_mtru.a+ibs);

			ops += blai*(blaij+blaj_2);

// Shift diagonal data if necessary

			paloc = _mtru.a+ibs;

			if (isup == jsup) {
				for (kii=0;kii<blai;kii++) paloc[kii*blai+kii].x += _dshift;
				for (kii=0;kii<blai;kii++) paloc[kii*blai+kii].y += _dshiftim;
			};

		};
	};

};

// Author: Kharchenko S.A.
// CFctC: Init current supernode row
//========================================================================================
void CFctC::Ich2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau) {

//	const char *funcname = "Ich2InitRow";

	dcmplx czero (0.0e0,0.0e0);

	int kii, ibs;
	int blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2;
	dcmplx aux;

	dcmplx *paloc;
	dcmplx *pfmask;
	dcmplx *pvloc;

	icycle++;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	blai_2 = blai*blai;
	ibsd = bsdia[_irow];
	ibsv = adrfmask[_irow];

	nlist = 0;

	lstloc[nlist] = _irow;
	imask[_irow] = icycle;

	int j = _mtrau.ia[_ilist];
	ibs = _mtrau.bsa[j];
	pfmask = fmasku+ibsv*blai;
	pvloc = diagg+ibsd;
	paloc = _mtrau.a+ibs;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	if (_type == 'F') {
		for (kii=0;kii<blai_2;kii++) *pfmask++ = *paloc++ + *pvloc++;
	} else {
		for (kii=0;kii<blai_2;kii++) *pfmask++ = *paloc++;
	};

	nlist++;

	for (j=_mtrau.ia[_ilist]+1;j<_mtrau.ia[_ilist+1];j++) {
		int jj = _mtrau.ja[j];
		lstloc[nlist] = jj;
		imask[jj] = icycle;

// Local copy

		ibs = _mtrau.bsa[j];
		blaj = _gmtra.sprndr[jj+1]-_gmtra.sprndr[jj];
		jbsv = adrfmask[jj];
		blaij = blai * blaj;

		pfmask = fmasku+jbsv*blai;
		paloc = _mtrau.a+ibs;
		for (kii=0;kii<blaij;kii++) *pfmask++ = *paloc++;

		nlist++;

	};

};

// Author: Kharchenko S.A.
// CFctC: Init mask and update arrays
//========================================================================================
void CFctC::Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ich2InitMask";

	dcmplx czero (0.0e0,0.0e0);

	int kii;
	int blai, blaj, jbsv, blaij;
	dcmplx *pfmask;

	int *addrgloc;
	dcmplx *guloc;
	int k, kloc;
	int *jgloc;
	int jcolmn;

	nupd = 0;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];

	for (k=iv[_ilstprv];k<ig[_ilstprv+1];k++) {

		kloc = k-ig[_ilstprv];
		
		jgloc = jg[_ilstprv];
		jcolmn = jgloc[kloc];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (imask[jcolmn] != icycle) {
			imask[jcolmn] = icycle;
			lstloc[nlist] = jcolmn;
			blaj = _gmtra.sprndr[jcolmn+1]-_gmtra.sprndr[jcolmn];
			jbsv = adrfmask[jcolmn];
			blaij = blai*blaj;
			pfmask = fmasku+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = czero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		addrgloc = addrg[_ilstprv];
		guloc = gu[_ilstprv];

		adrupdu[nupd] = guloc+addrgloc[kloc];
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctC: Update current row
//========================================================================================
void CFctC::Ich2UpdateRow (int _irow, int _irwprv, // Update current row
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv, blaik;

	dcmplx *pfmloc;
	dcmplx *plloc;
	dcmplx *pgloc;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	blak = _gmtra.sprndr[_irwprv+1]-_gmtra.sprndr[_irwprv];
	blaik = blai*blak;

	for (int k=0;k<nupd;k++) {
		int icol = lstupd[k];
		blaj = _gmtra.sprndr[icol+1]-_gmtra.sprndr[icol];
		jbsv = adrfmask[icol];
		ibs = jbsv*blai;

		dcmplx *elemg = adrupdu[k];
		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = fmasku+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = uelem+kikb;
				pgloc = elemg+kjkb;
				dcmplx *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
//					*pfmloc -= *plloc * *pgloc;
					(*pfmloc).MnEq (*plloc * *pgloc);
				};
				pfmloc += blai;
				kjkb += blak;
			};
			kikb += blak;
		};

		ops += blaik * blaj;

	};

};

// Author: Kharchenko S.A.
// CFctC: Perform filtering of the current row
//========================================================================================
void CFctC::Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ich2FiltrRow";

	int nlist2 = 0;
	int kii, kjj, blai, blaj, ibsv, jbsv, jbsd, blaij;
	dcmplx *pfmasku, *pfmaskuu, *pdiag;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	ibsv = adrfmask[_irow];

	for (int i=1;i<nlist;i++) {

		int j = lstloc[i];

		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		jbsv = adrfmask[j];
		jbsd = bsdia[j];
		blaij = blai*blaj;

		double aux = 0.0e0;

		pfmasku = fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = abs(*pfmasku++);
			if (auxloc > aux) aux = auxloc;
		};

		if (aux >= _tau2) {
			lstupd[nlist2] = j;
			nlist2++;
		} else {
			pfmaskuu = fmasku+ibsv*blai;
			for (kii=0;kii<blai;kii++) {
				pfmasku = fmasku+jbsv*blai+kii;
				pdiag = diagg+jbsd;
				for (kjj=0;kjj<blaj;kjj++) {
//					double auxloc = _fmasku[jbsv*blai+kjj*blai+kii];
					double auxloc = abs(*pfmasku);
					auxloc *= _theta;
//					_fmasku[ibsv*blai+kii*blai+kii] += auxloc;
					*pfmaskuu += auxloc;
//					_diagg[jbsd+kjj*blaj+kjj] += auxloc;
					*pdiag += auxloc;
					if (kjj<blaj-1) {
						pfmasku += blai;
						pdiag += blaj+1;
					};
				};
				if (kii < blai-1) pfmaskuu += blai+1;
			};
		};

	};

	nlist = nlist2;

};

// CFctC: Compute diagonal pivot and scale current row
// Author: Kharchenko S.A.
//========================================================================================
void CFctC::Ich2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
								const CGSMatrix &_gmtra) {

//	const char *funcname = "Ich2PivScaleRow";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, ibs;
	int blai, blaj, ibsv, jbsv;
	dcmplx caux;

// Copy diagonal supernode

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	ibsv = adrfmask[_irow];
	ibs = ibsv * blai;

	for (kii=0;kii<blai*blai;kii++) aloc[kii] = fmasku[ibs+kii];

// Compute spectral decomposition of the diagonal supernode

	int blal = blai;
	int info = 0;

//	int MyZheev (char *, char *, int *, dcmplx *, 
//			int *, double *, dcmplx *, int *, double *, int *);
//	MyZheev ("V", "U", &blal,
//				aloc, &blal, eig, 
//				work, &lwork, rwork, &info);
	zheev_ ("V", "U", &blal,
				aloc, &blal, eig, 
				work, &lwork, rwork, &info);

	if (info != 0) throw " Error in the Lapack routine ZHEEV";

	ops += 8 * blai * blai * blai;

//	ofstream fout ("CheckPivIch.dat");
//	fout << " irow = " << _irow << endl;
//	OutArr(fout," Eig ",blal,eig);
//	cout << " irow = " << _irow << endl;
//	OutArr(cout," Eig ",blal,eig);

// Check the eigenvalues

	for (kii=0;kii<blai;kii++) {
		if (eig[kii] < 0.0e0) {
			throw " Ich2PivScaleRow: indefinite supernode ";
		};
	};

// Modify eigenvalues if necessary

	for (kii=0;kii<blai;kii++) {
		eig[kii] = sqrt(eig[kii]);
		if (eig[kii] < _pivmin) {
			eig[kii] = _pivmin;
			nmodsv++;
		};
		dpiv[ibsv+kii] = eig[kii];
		eig[kii] = 1.0e0/eig[kii];
	};

// Compute diagonal supernodes

	int kiib;
	dcmplx *pdeloc;
	dcmplx *pfmloc;
	dcmplx *palloc;

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = deleml+kii;
		palloc = aloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_deleml[kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//			*pdeloc = eig[kii] * conj(*palloc++);
			dcmplx temp = conj(*palloc++);
			*pdeloc = eig[kii] * temp;
			pdeloc += blai;
		};
		kiib += blai;
	};

	for (kii=0;kii<blai*blai;kii++) delemu[kii] = conj(deleml[kii]);

	ops += blai * blai;

// Perform scaling of the data in the list

	int kjjb, blaij;

	for (int i=0;i<nlist;i++) {

		int j = lstupd[i];
		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		jbsv = adrfmask[j];
		blaij = blai*blaj;

		ibs = jbsv*blai;

		for (kii=0;kii<blai;kii++) {
			palloc = aloc+kii;
			for (kjj=0;kjj<blaj;kjj++) {
				kjjb = kjj*blai;
				pfmloc = fmasku+ibs+kjjb;
				pdeloc = deleml+kii;
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _deleml[kkk*blai+kii] * _fmasku[ibs+kjj*blai+kkk];
//					caux += *pdeloc * *pfmloc++;
					caux.PlEq (*pdeloc * *pfmloc++);
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = caux;
				palloc += blai;
			};
		};

		pfmloc = fmasku+ibs;
		palloc = aloc;
//		for (kii=0;kii<blaij;kii++) _fmasku[ibs+kii] = _aloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmloc++ = *palloc++;

		ops += blai * blaij;

	};

};

// Author: Kharchenko S.A.
// CFctC: Store current row
//========================================================================================
void CFctC::Ich2StoreRow (int _ilist, int _irow, double _tau1, // Store current row
							int &_nzrtot, 
							const CGSMatrix &_gmtra) {

	const char *funcname = "Ich2StoreRow";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	dcmplx *pg, *pfmask, *pd;

// Count the number of elements to be stored

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];

	int nz = blai*blai;

	int i;
	for (i=0;i<nlist;i++) {
		int j = lstupd[i];
		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		nz += blai*blaj;
	};

// Allocate the memory

	int *jgloc;
	int *adrgloc;
	dcmplx *guloc;

	jgloc = new int [nlist+1];
	if (!jgloc) MemoryFail (funcname);
	adrgloc = new int [nlist+1];
	if (!adrgloc) MemoryFail (funcname);
	guloc = new dcmplx [nz];
	if (!guloc) MemoryFail (funcname);

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gu[_ilist] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

// Store current diagonal supernode

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	jgloc[nzjgloc] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = delemu+kiib;
		pg = guloc+nzloc+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gu [_nznew+kjj*blai+kii] = _delemu[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
	adrgloc[nzjgloc] = nzloc;
	nzjgloc++;
	nzloc += blai*blai;

	for (i=0;i<nlist;i++) {
		int j = lstupd[i];
		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		jbsv = adrfmask[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = abs(*pfmask++);
			if (auxloc > aux) aux = auxloc;
		};
		adrgloc[nzjgloc] = nzloc;
		if (aux >= _tau1) {
			jgloc[nzjgloc] = j;
		} else {
			jgloc[nzjgloc] = -j;
			_nzrtot += blaij;
		};
		pg = guloc+nzloc;
		pfmask = fmasku+jbsv*blai;
//		for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
		for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

		nzjgloc++;
		nzloc += blaij;
	};

	ig[_ilist+1] = ig[_ilist]+nzjgloc;

};

// Author: Kharchenko S.A.
// CFctC: Store the resulting U objects
//========================================================================================
void CFctC::Ich2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
								const CGSMatrix &_gmtra, 
								CSMatrixCS &_mtru, CSMatrixCS &_mtru2) {

//	const char *funcname = "Ich2StoreScaleU";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	dcmplx caux;

	dcmplx *pscloc;
	dcmplx *pgloc;
	dcmplx *paloc;

// Count the numbers of elements in the matrices to be stored

	int isupbeg, isupend, nlistl;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	int *jgloc;

	int ilist, isup, j, jloc, jj, jjloc, ibs;

	int nzjlu=0, nzlu=0, nzjlu2=0, nzlu2=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jj >= 0) {
				nzjlu++;
				nzlu += blai*blaj;
			};
		};
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jjloc > isupend) {
				nzjlu2++;
				nzlu2 += blai*blaj;
			};
		};
	};

// Allocate and init first order matrices

	CSMatrixCS mtrdummy;

	_mtru = mtrdummy;

	CSMatrixCS temp (nlistl,nzjlu,nzlu);

	temp.m = _gmtra.m;
	temp.n = _gmtra.n;

	temp.blamx = blamx;

	temp.sprndr[0] = 0;
	temp.sprndc[0] = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		temp.sprndr[ilist+1] = 0;
		temp.sprndc[ilist+1] = 0;
	};

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jj >= 0) {
				temp.ja[nzjlu] = jjloc;
				temp.bsa[nzjlu] = nzlu;
				nzjlu++;
				nzlu += blai*blaj;
			};
		};
		temp.list[ilist] = isup;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nzja = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	dcmplx *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jj >= 0) {
				ibs = adrgloc[jloc];
				for (kii=0;kii<blai*blaj;kii++) _mtru.a[nzlu+kii] = guloc[ibs+kii];
				nzjlu++;
				nzlu += blai*blaj;
			};
		};
	};

// Scale

	if (_param.sclpttype < 2) {

		for (ilist=0;ilist<nlistl;ilist++) {
			isup = isupbeg+ilist;
			blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
			ibssc = bsdia[isup];
			for (int j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {
				jbs = _mtru.bsa[j];
				int jj = _mtru.ja[j];
				blaj = _gmtra.sprndr[jj+1]-_gmtra.sprndr[jj];
				blaij = blai*blaj;
				jjbssc = bsdia[jj];
				if (jj == isup) {

					kiib = 0;
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = aloc+kii;
						for (kjj=0;kjj<blai;kjj++) {
							pgloc = _mtru.a+jbs+kjjb;
							pscloc = sclu+ibssc+kiib;
							caux = czero;
							for (kkk=0;kkk<blai;kkk++) {
//								caux += *pscloc * *pgloc;
								caux.PlEq (*pscloc * *pgloc);
								pgloc++;
								pscloc++;
							};
							*paloc = caux;
							paloc += blai;
							kjjb += blai;
						};
						kiib += blai;
					};
					paloc = aloc;
					pgloc = _mtru.a+jbs;
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;

	 				ops += 2*blaij*blai;

				} else {

					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = aloc+kii;
						for (kjj=0;kjj<blaj;kjj++) {
							caux = czero;
							pgloc = _mtru.a+jbs+kii;
							pscloc = sclinvu+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								caux += *pgloc * *pscloc;
								caux.PlEq (*pgloc * *pscloc);
								pgloc += blai;
								pscloc += blaj;
							};
							*paloc = caux;
							paloc += blai;
							kjjb += blai;
						};
					};
					paloc = aloc;
					pgloc = _mtru.a+jbs;
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;

					ops += blaij*blaj;
				};
			};
		};

	};

// Allocate and init second order matrices

	_mtru2 = mtrdummy;

	CSMatrixCS temp2 (nlistl,nzjlu2,nzlu2);

	temp2.m = _gmtra.m;
	temp2.n = _gmtra.n;

	temp2.blamx = blamx;

	temp2.sprndr[0] = 0;
	temp2.sprndc[0] = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		temp2.sprndr[ilist+1] = 0;
		temp2.sprndc[ilist+1] = 0;
	};

	temp2.ia[0] = 0;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jjloc > isupend) {
				temp2.ja[nzjlu2] = jj;
				temp2.bsa[nzjlu2] = nzlu2;
				nzjlu2++;
				nzlu2 += blai*blaj;
			};
		};
		temp2.list[ilist] = isup;
		temp2.ia[ilist+1] = nzjlu2;
	};

	int nlistn=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = temp2.list[ilist];
		if (temp2.ia[ilist+1] > temp2.ia[ilist]) {
			temp2.list[nlistn] = isup;
			temp2.ia[nlistn+1] = temp2.ia[ilist+1];
			nlistn++;
		};
	};

	temp2.nlist = nlistn;
	temp2.nzja = nzjlu2;
	temp2.nza = nzlu2;
	temp2.nzatot = nzlu2;

	_mtru2 = temp2;

	temp2 = mtrdummy;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jjloc > isupend) {
				ibs = adrgloc[jloc];
				for (kii=0;kii<blai*blaj;kii++) _mtru2.a[nzlu2+kii] = guloc[ibs+kii];
				nzjlu2++;
				nzlu2 += blai*blaj;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctC: Perform factorization of the current block row
//========================================================================================
void CFctC::Ich2FctBlkRow (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, CSMatrixCS &_mtru,
							CSMatrixCS &_mtru2) {

//	const char *funcname = "Ich2FctBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int nzrtot;
	int irwprv, ilstprv, icolmn, j, jloc, ibs;
	int isupbeg, isupend, nlistl, blai, blaj, kii;
	bool elmisr;
	int *jgloc;
	int *adrgloc;
	dcmplx *guloc;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	nzrtot = 0;

	for (ilist=0;ilist<=_mtrau.nlist;ilist++) ibegm[ilist] = -1;

	ig[0] = 0;

	for (ilist=0;ilist<_mtrau.nlist;ilist++) {

		irow = _mtrau.list[ilist];

		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];

// Init current scaled row

		Ich2InitRow ('F', ilist, irow, _gmtra, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			jgloc = jg[ilstprv];
			adrgloc = addrg[ilstprv];
			guloc = gu[ilstprv];

			j = iv[ilstprv];

			jloc = j-ig[ilstprv];

			icolmn = jgloc[jloc];
			ibs = adrgloc[jloc];

			irwprv = _mtrau.list[ilstprv];
			int blailoc = _gmtra.sprndr[irwprv+1]-_gmtra.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _gmtra.sprndr[icolmnl+1]-_gmtra.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = conj(guloc[ibs+kii]);

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask (elmisr, irow, ilstprv, _gmtra);

// Perform updates

			Ich2UpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						_gmtra);

// Compute and store the pivot value and scale current supernode row

		Ich2PivScaleRow (irow, _param.pivmin,
							_gmtra);

// Sort the list of indices

		if (nlist != 0) qsort (lstupd, nlist, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist != 0) {
			iv[ilist] = ig[ilist]+1;
			madj[ilist] = ibegm[lstupd[0]-isupbeg];
			ibegm[lstupd[0]-isupbeg] = ilist;
		};

// Store the data in the global arrays

		Ich2StoreRow (ilist, irow, _param.tau1, 
						nzrtot,
						_gmtra);

	};

// Store the resulting U objects

	Ich2StoreScaleU (_iblkrow, _param, 
						_gmtra,
						_mtru, _mtru2);

	nzr += nzrtot;

// Store the hystogram of values

	int nloc = 0;

	for (ilist=0;ilist<_mtrau.nlist;ilist++) {
		irow = _mtrau.list[ilist];
		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];
		nloc += blai;
	};

	double *dpivarr;

	dpivarr = new double [nloc];

	nloc = 0;

	for (ilist=0;ilist<_mtrau.nlist;ilist++) {
		irow = _mtrau.list[ilist];
		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];
		int ibsv = adrfmask[irow];;
		for (kii=0;kii<blai;kii++) dpivarr[nloc+kii] = dpiv[ibsv+kii];
		nloc += blai;
	};

	hyst.UpdateHystogram (nloc, dpivarr);

	delete [] dpivarr;

// Free local factorization arrays 

	for (ilist=0;ilist<nlistl;ilist++) {
		delete [] jg[ilist];
		delete [] addrg[ilist];
//		delete [] gl[ilist];
		delete [] gu[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctC: Init mask and update arrays
//========================================================================================
void CFctC::Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtru2) {

//	const char *funcname = "Ich2InitMask";

	dcmplx czero (0.0e0,0.0e0);

	int kii;
	int blai, blaj, jbsv, blaij, ibs;
	dcmplx *pfmask;

	int k;
	int jcolmn;

	nupd = 0;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];

	for (k=iv[_ilstprv];k<_mtru2.ia[_ilstprv+1];k++) {

		jcolmn = _mtru2.ja[k];
		ibs = _mtru2.bsa[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (imask[jcolmn] != icycle) {
			imask[jcolmn] = icycle;
			lstloc[nlist] = jcolmn;
			blaj = _gmtra.sprndr[jcolmn+1]-_gmtra.sprndr[jcolmn];
			jbsv = adrfmask[jcolmn];
			blaij = blai*blaj;
			pfmask = fmasku+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = czero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		ibs = _mtru2.bsa[k];

		adrupdu[nupd] = _mtru2.a+ibs;
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctC: Store current row
//========================================================================================
void CFctC::Ich2StoreRow (int _ilist, int _irow, // Store current row
							const CGSMatrix &_gmtra) {

	const char *funcname = "Ich2StoreRow";

	int kii;
	int blai, blaj, jbsv, blaij;
	dcmplx *pg, *pfmask;

// Count the number of elements to be stored

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];

	int nz = 0;

	int i;
	for (i=0;i<nlist;i++) {
		int j = lstloc[i];
		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		nz += blai*blaj;
	};

// Allocate the memory

	int *jgloc;
	int *adrgloc;
	dcmplx *guloc;

	jgloc = new int [nlist+1];
	if (!jgloc) MemoryFail (funcname);
	adrgloc = new int [nlist+1];
	if (!adrgloc) MemoryFail (funcname);
	guloc = new dcmplx [nz];
	if (!guloc) MemoryFail (funcname);

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gu[_ilist] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

	for (i=0;i<nlist;i++) {
		int j = lstloc[i];
		blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
		jbsv = adrfmask[j];
		blaij = blai*blaj;
		adrgloc[nzjgloc] = nzloc;
		jgloc[nzjgloc] = j;
		pg = guloc+nzloc;
		pfmask = fmasku+jbsv*blai;
//		for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
		for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

		nzjgloc++;
		nzloc += blaij;
	};

	ig[_ilist+1] = ig[_ilist]+nzjgloc;

};

// Author: Kharchenko S.A.
// CFctC: Store the resulting U objects
//========================================================================================
void CFctC::Ich2StoreU (int _iblkrow, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixCS &_mtru) {

//	const char *funcname = "Ich2StoreU";

	dcmplx czero (0.0e0,0.0e0);

	int kii;
	int blai, blaj;
	dcmplx caux;

// Count the numbers of elements in the matrices to be stored

	int isupbeg, isupend, nlistl;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	int *jgloc;

	int ilist, isup, j, jloc, jj, jjloc, ibs;

	int nzjlu=0, nzlu=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			nzjlu++;
			nzlu += blai*blaj;
		};
	};

// Allocate and init matrices

	CSMatrixCS mtrdummy;

	_mtru = mtrdummy;

	CSMatrixCS temp (nlistl,nzjlu,nzlu);

	temp.m = _gmtra.m;
	temp.n = _gmtra.n;
	temp.nsupr = nlistl;
	temp.nsupc = nlistl;

	for (ilist=0;ilist<=nlistl;ilist++) temp.sprndr[ilist] = 0;
	for (ilist=0;ilist<=nlistl;ilist++) temp.sprndc[ilist] = 0;

	temp.blamx = blamx;

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			temp.ja[nzjlu] = jjloc;
			temp.bsa[nzjlu] = nzlu;
			nzjlu++;
			nzlu += blai*blaj;
		};
		temp.list[ilist] = isup;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nzja = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	dcmplx *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			ibs = adrgloc[jloc];
			for (kii=0;kii<blai*blaj;kii++) _mtru.a[nzlu+kii] = guloc[ibs+kii];
			nzjlu++;
			nzlu += blai*blaj;
		};
	};

};

// Author: Kharchenko S.A.
// CFctC: Perform updates of the current block row
//========================================================================================
void CFctC::Ich2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
								const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, const CSMatrixCS &_mtru2, 
								CSMatrixCS &_mtru) {

//	const char *funcname = "Ich2UpdateBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int irwprv, ilstprv, icolmn, j, ibs;
	int isupbeg, isupend, nlistl, blai, blaj, kii;
	bool elmisr;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

// Init arrays that support transposed structures search

	Ilu2InitTransp (isupbeg, isupend, _mtru2);

	ig[0] = 0;

	for (ilist=0;ilist<_mtrau.nlist;ilist++) {

		irow = _mtrau.list[ilist];

		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];

// Init current row

		Ich2InitRow ('U', ilist, irow, _gmtra, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			j = iv[ilstprv];

			icolmn = _mtru2.ja[j];
			ibs = _mtru2.bsa[j];

			irwprv = _mtru2.list[ilstprv];
			int blailoc = _gmtra.sprndr[irwprv+1]-_gmtra.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _gmtra.sprndr[icolmnl+1]-_gmtra.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = conj(_mtru2.a[ibs+kii]);

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask (elmisr, irow, ilstprv, _gmtra, _mtru2);

// Perform updates

			Ich2UpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend, _mtru2);

// Sort the list of indices

		if (nlist != 0) qsort (lstloc, nlist, sizeof(int), compint);

// Store the data in the global arrays

		Ich2StoreRow (ilist, irow, _gmtra);

	};

// Store the resulting U objects

	Ich2StoreU (_iblkrow, 
					_gmtra, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<nlistl;ilist++) {
		delete [] jg[ilist];
		delete [] addrg[ilist];
		delete [] gl[ilist];
		delete [] gu[ilist];
	};

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
