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
#include "svector.h"
#include "corr.h"
#include "slvparam.h"
#include "hyst.h"
#include "mvm.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A. 
// Description: Perform in-place filtering of the system of equations
// CGSMatrixCS::FilterLinearSystem()
//========================================================================================
void CGSMatrixCS::FilterLinearSystem (// Perform in-place filtering of the system of equations
													CGSMatrixCS &_mtra, CSVectorC &_rhs,
													int *_sprndsflt) {

	const char *funcname = "FilterLinearSystem";

	dcmplx czero (0.0e0,0.0e0);

// Main cycle over the block columns

	int i, isupbeg, isupend, nsupl;
	int isuploc, isup, j, jsup, blaiold, blainew, blajold, blajnew;
	int ibs, ibsnew, iblk, kii, kjj, ibsvold, ibsvnew;

	for (i=0;i<_mtra.nblksc;i++) {

// Read current block

		_mtra.ReadBlock (i);

// Cycle over supernodes

		iblk = i;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

		ibsnew = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blaiold = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			blainew = _sprndsflt[isup+1]-_sprndsflt[isup];
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				ibs = _mtra.mtrarr[iblk].bsa[j];
				blajold = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				blajnew = _sprndsflt[jsup+1]-_sprndsflt[jsup];
				for (kjj=0;kjj<blainew;kjj++) {
					for (kii=0;kii<blajnew;kii++) {
						_mtra.mtrarr[iblk].a[ibsnew+kjj*blajnew+kii] = _mtra.mtrarr[iblk].a[ibs+kjj*blajold+kii];
					};
				};
				_mtra.mtrarr[iblk].bsa[j] = ibsnew;
				ibsnew += blainew * blajnew;
			};
		};
		_mtra.mtrarr[iblk].nza = ibsnew;
		_mtra.mtrarr[iblk].nzatot = ibsnew;

// Free current block

		_mtra.ReWriteBlock (i);
		_mtra.FreeBlock (i);

	};

// Modify rhs

	for (isup=0;isup<_mtra.nsupc;isup++) {
		ibsvold = _mtra.sprndc[isup];
		ibsvnew = _sprndsflt[isup];
		blainew = _sprndsflt[isup+1]-_sprndsflt[isup];
		for (kii=0;kii<blainew;kii++) _rhs.vect[ibsvnew+kii] = _rhs.vect[ibsvold+kii];
	};

	_rhs.nv = _sprndsflt[_mtra.nsupc];

// Store new sprnds partitioning

	for (i=0;i<=_mtra.nsupc;i++) _mtra.sprndc[i] = _sprndsflt[i];
	for (i=0;i<=_mtra.nsupc;i++) _mtra.sprndr[i] = _sprndsflt[i];

// Recompute blocks sizes

	if (_rhs.nparts != _mtra.nblksc) {
		throw " CGSMatrixCS::FilterLinearSystem: wrong number of blocks ";
	};

	for (i=0;i<_mtra.nblksc;i++) {
		isupbeg = _mtra.blksc[i];
		isupend = _mtra.blksc[i+1]-1;
		_mtra.bl2ndr[i+1] = _mtra.sprndc[isupend+1];
		_mtra.bl2ndc[i+1] = _mtra.sprndc[isupend+1];
		_rhs.inda[i] = _mtra.sprndc[isupbeg];
		_rhs.blksv[i+1] = _mtra.sprndc[isupend+1];
	};

};

// Author: Kharchenko S.A. 
// Description: Perform in-place symmetric block scaling of the square matrix stored by block columns
// CGSMatrixCS::ExplicitSymmetricBlockScaling()
//========================================================================================
void CGSMatrixCS::ExplicitSymmetricBlockScaling (ofstream &_fout, const CSlvParam &_param, // Perform in-place symmetric block scaling of the square matrix stored by block columns
																	CGSMatrixCS &_mtra, CSVectorC &_rhs,
																	int *_sprndsflt,
																	CCorrectorC &_scla) {

	const char *funcname = "ExplicitSymmetricBlockScaling";

	dcmplx czero (0.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtra.nsupc];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;

	for (i=0;i<_mtra.nblksc;i++) {
		for (j=_mtra.blksc[i];j<_mtra.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Init the mask

	int icycle = -1;

	int *imask;

	imask = new int [_mtra.nsupc];
	if (!imask) MemoryFail (funcname);

	for (i=0;i<_mtra.nsupc;i++) imask[i] = icycle;

// Estimate the memory

	int ilist, iblk, isupbeg, isupend, nsupl, blai, blaj;
	int isuploc, isup, jsup, nzal;

	int nzscl = 0;
	int nzacolmax = 0;
	int blamax = 0;

	int *pldcorrarr = _scla.GetLdcorrarr ();
	int *pcorrsizes = _scla.GetCorrsizes ();
	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		nzscl = 0;

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			if (blai > blamax) blamax = blai;
			nzscl += blai*blai;
			nzal = 0;
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				if (blaj > blamax) blamax = blaj;
				nzal += blai*blaj;
			};
			if (nzal > nzacolmax) nzacolmax = nzal;
		};

// Store the memory

		pldcorrarr[ilist] = nzscl;
		pcorrsizes[ilist] = 1;

	};

// Init the hystogram

	CHyst hyst;

// Compute the scaling

	double threshflt = _param.thrfltmtr;

	CGSMatrixCS gmtrdummy;

	for (i=0;i<=_mtra.nsupc;i++) _sprndsflt[i] = 0;

	for (i=0;i<_mtra.nblksc;i++) {

// Read current block

		_mtra.ReadBlock (i);

// Allocate the scaling 

		_scla.AllocateCorr (i);

// Compute scaling

		gmtrdummy.ComputeSymmetricBlockScaling (i, threshflt, hyst,
																_mtra, 
																_sprndsflt, _scla);

// Write the scaling to the disk

		_scla.WriteCorr (i);
		_scla.FreeCorr (i);

// Free current block

		_mtra.FreeBlock (i);

	};

	for (i=0;i<_mtra.nsupc;i++) _sprndsflt[i+1] = _sprndsflt[i]+_sprndsflt[i+1];

// Output hystogram

	hyst.MakeGlobalHystogram ();

	_fout << " ExplSymmBlockScaling: " << hyst;
	cout  << " ExplSymmBlockScaling: " << hyst;

	int nfiltr = _mtra.sprndc[_mtra.nsupc]-_sprndsflt[_mtra.nsupc];

	_fout << " ThreshFlt = " << threshflt << " Nfiltr = " << nfiltr << endl;
	cout  << " ThreshFlt = " << threshflt << " Nfiltr = " << nfiltr << endl;

// Perform scaling block column by block column

	int *listscl, *ibsscl;

	listscl = new int [_mtra.nblksc];
	if (!listscl) MemoryFail (funcname);
	ibsscl = new int [_mtra.nsupc+1];
	if (!ibsscl) MemoryFail (funcname);

	ibsscl[0] = 0;

	for (isup=0;isup<_mtra.nsupc;isup++) {
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		ibsscl[isup+1] = ibsscl[isup] + blai*blai;
	};

	int jblk, nlistblk;;

	for (i=0;i<_mtra.nblksc;i++) {

		cout  << " Scale Iblk = " << i << endl;

// Read current block

		_mtra.ReadBlock (i);

// Compute the list

		icycle++;

		iblk = i;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

		nlistblk = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				jblk = sp2blk[jsup];
				if (imask[jblk] != icycle) {
					listscl[nlistblk] = jblk;
					nlistblk++;
					imask[jblk] = icycle;
				};
			};
		};

// Read from the disk the list of scaling blocks if necessary

		for (j=0;j<nlistblk;j++) {
			jblk = listscl[j];
			_scla.ReadCorr (jblk);
		};

// Scale

		gmtrdummy.SymmetricScaleBlockColumn (i, _mtra, sp2blk, _rhs, ibsscl, _scla);

// Free the list of blocks

		for (j=0;j<nlistblk;j++) {
			jblk = listscl[j];
			_scla.FreeCorr (jblk);
		};

// Write current block

		_mtra.ReWriteBlock (i);
		_mtra.FreeBlock (i);

	};

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout << " Scaling block rows statistics: " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Scaling block rows statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work arrays

	delete [] sp2blk;
	delete [] imask;
	delete [] listscl;
	delete [] ibsscl;

};

// Author: Kharchenko S.A. 
// Description: Compute scaling of the current block
// CGSMatrixCS::ComputeSymmetricBlockScaling()
//========================================================================================
void CGSMatrixCS::ComputeSymmetricBlockScaling (int _iblk, double _threshflt, // Compute scaling of the current block
																CHyst &_hyst,
																CGSMatrixCS &_mtra, 
																int *_sprndsflt, CCorrectorC &_scla) {

	const char *funcname = "ComputeSymmetricBlockScaling";

	dcmplx czero (0.0e0,0.0e0);

// Determine the maximal size

	int i, j, isup, jsup, blai, blaj;

	int nlocmax = 0;

	for (i=0;i<_mtra.mtrarr[_iblk].nlist;i++) {
		isup = _mtra.mtrarr[_iblk].list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtra.mtrarr[_iblk].ia[i];j<_mtra.mtrarr[_iblk].ia[i+1];j++) {
			jsup = _mtra.mtrarr[_iblk].ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;
	int nsupl = isupend-isupbeg+1;

// Allocate work arrays

	dcmplx *aloc, *uloc, *vloc, *work;
	double *svloc, *rwork;
	int lwork = nlocmax * 10;
	int info;

	aloc = new dcmplx [nlocmax*nlocmax];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [nlocmax*nlocmax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [nlocmax*nlocmax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [nlocmax];
	if (!svloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Main cycle over the supernode columns

	int ibs;
	int kii, kjj;
	int isuploc;
	double daux;
	dcmplx caux;

	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	int nzscl = 0;

	for (isuploc=0;isuploc<nsupl;isuploc++) {
		isup = isuploc+isupbeg;

		blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];

// Scan current block data

		int icheck = 0;

		for (j=_mtra.mtrarr[_iblk].ia[isuploc];j<_mtra.mtrarr[_iblk].ia[isuploc+1];j++) {

			jsup = _mtra.mtrarr[_iblk].ja[j];

			if (jsup == isup) {
				icheck = 1;

				ibs = _mtra.mtrarr[_iblk].bsa[j];

				for (kii=0;kii<blai*blai;kii++) aloc[kii] = conj(_mtra.mtrarr[_iblk].a[ibs+kii]);

// Compute symmetric eigendecomposition

				zheev_ ("V", "U", &blai,
							aloc, &blai, svloc, 
							work, &lwork, rwork, &info);

				if (info != 0) throw " Error in the Lapack routine ZHEEV";

				_hyst.UpdateHystogram (blai,svloc);

				int niflt = 0;

				for (kii=0;kii<blai;kii++) {
					if (svloc[kii] > _threshflt) niflt++;
				};

				_sprndsflt[isup+1] = niflt;

// Compute and store the scaling and its inverse

				for (kii=0;kii<blai;kii++) {
					daux = svloc[kii];
					daux = sqrt(daux);
					for (kjj=0;kjj<blai;kjj++) {
						uloc[kii*blai+kjj] = conj(aloc[kjj*blai+kii]) * daux;
					};
				};

				for (kii=0;kii<blai*blai;kii++) pyarr[_iblk][nzscl+kii] = uloc[kii];

				for (kii=0;kii<blai;kii++) {
					daux = svloc[kii];
					daux = sqrt(daux);
					daux = 1.0e0 / daux;
					for (kjj=0;kjj<blai;kjj++) {
						caux = aloc[kii*blai+kjj] * daux;
						uloc[kii*blai+kjj] = caux;
					};
				};

				for (kii=0;kii<blai*blai;kii++) pwarr[_iblk][nzscl+kii] = uloc[kii];

			};

		};

		if (icheck == 0) throw " CGSMatrixCS::ComputeSymmetricBlockScaling: index not found ";

		nzscl += blai*blai;

	};

// Free work arrays

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;

};

// Author: Kharchenko S.A. 
// Description: Perform symmetric scaling of the current block
// CGSMatrixCS::SymmetricScaleBlockColumn()
//========================================================================================
void CGSMatrixCS::SymmetricScaleBlockColumn (int _iblk, // Perform symmetric scaling of the current block
															CGSMatrixCS &_mtra, int *_sp2blk, CSVectorC &_rhs,
															int *_ibsscl, CCorrectorC &_scla) {

	const char *funcname = "SymmetricScaleBlockColumn";

	dcmplx czero (0.0e0,0.0e0);

// Determine the maximal size

	int i, j, isup, jsup, blai, blaj;

	int nlocmax = 0;

	for (i=0;i<_mtra.mtrarr[_iblk].nlist;i++) {
		isup = _mtra.mtrarr[_iblk].list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtra.mtrarr[_iblk].ia[i];j<_mtra.mtrarr[_iblk].ia[i+1];j++) {
			jsup = _mtra.mtrarr[_iblk].ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;
	int nsupl = isupend-isupbeg+1;

// Allocate work arrays

	dcmplx *aloc, *uloc;

	aloc = new dcmplx [nlocmax*nlocmax];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [nlocmax*nlocmax];
	if (!uloc) MemoryFail (funcname);

// Main cycle over the supernode columns

	int ibs, ibsv, ibssc, jbssc, jblk;
	int kii, kjj, kkk, isuploc, jsupbeg;
	dcmplx caux;

	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	int nzscl = 0;

	for (isuploc=0;isuploc<nsupl;isuploc++) {

		isup = isuploc+isupbeg;

		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

// Scan current column

		ibssc = _ibsscl[isup]-_ibsscl[isupbeg];

		for (j=_mtra.mtrarr[_iblk].ia[isuploc];j<_mtra.mtrarr[_iblk].ia[isuploc+1];j++) {

			jsup = _mtra.mtrarr[_iblk].ja[j];
			ibs = _mtra.mtrarr[_iblk].bsa[j];

			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];

			jblk = _sp2blk[jsup];
			jsupbeg = _mtra.blksc[jblk];

			jbssc = _ibsscl[jsup]-_ibsscl[jsupbeg];

// Scale

			for (kii=0;kii<blaj;kii++) {
				for (kjj=0;kjj<blai;kjj++) {
					caux = czero;
					for (kkk=0;kkk<blai;kkk++) {
						caux += _mtra.mtrarr[_iblk].a[ibs+kii*blai+kkk] * pwarr[_iblk][ibssc+kjj*blai+kkk];
					};
					aloc[kjj*blaj+kii] = caux;
				};
			};

			for (kii=0;kii<blaj;kii++) {
				for (kjj=0;kjj<blai;kjj++) {
					caux = czero;
					for (kkk=0;kkk<blaj;kkk++) {
						caux += conj(pwarr[jblk][jbssc+kii*blaj+kkk]) * aloc[kjj*blaj+kkk];
					};
					uloc[kii*blai+kjj] = caux;
				};
			};

			for (kii=0;kii<blai*blaj;kii++) _mtra.mtrarr[_iblk].a[ibs+kii] = uloc[kii];

		};

// Scale right hand side

		ibsv = _mtra.sprndr[isup];

		for (kii=0;kii<blai;kii++) {
			caux = czero;
			for (kjj=0;kjj<blai;kjj++) {
				caux += conj(pwarr[_iblk][ibssc+kii*blai+kjj]) * _rhs.vect[ibsv+kjj];
			};
			aloc[kii] = caux;
		};

		for (kii=0;kii<blai;kii++) _rhs.vect[ibsv+kii] = aloc[kii];

	};

// Free work arrays

	delete [] aloc;
	delete [] uloc;

};

// Author: Kharchenko S.A. 
// Description: Perform in-place scaling of the block rows of the square matrix stored by block columns
// CGSMatrixCS::ScaleBlockRows()
//========================================================================================
void CGSMatrixCS::ScaleBlockRows (ofstream &_fout, const CSlvParam &_param, // Perform in-place scaling of the block rows of the square matrix stored by block columns
												CSMatrix &_mtrasp, CGSMatrixCS &_mtra, CSVectorC &_rhs) {

	const char *funcname = "ScaleBlockRows";

	dcmplx czero (0.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtra.nsupc];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;

	for (i=0;i<_mtra.nblksc;i++) {
		for (j=_mtra.blksc[i];j<_mtra.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Allocate index array

	int *indarr;

	indarr = new int [_mtrasp.nzja];
	if (!indarr) MemoryFail (funcname);

// Compute transposed matrix

	CSMatrix at;

	at = _mtrasp.TranspMtrRect (indarr);

// Init the mask

	int icycle = -1;

	int *imask;
	int *jamask;

	imask = new int [_mtra.nsupc];
	if (!imask) MemoryFail (funcname);
	jamask = new int [_mtrasp.nzja];
	if (!jamask) MemoryFail (funcname);

	for (i=0;i<_mtra.nsupc;i++) imask[i] = icycle;
	for (i=0;i<_mtrasp.nzja;i++) jamask[i] = 1;

// Perform scaling block column by block column

	CGSMatrixCS gmtrdummy;

	CSMatrixCS mtraoff;

	int *ja2isup, *ja2blk, *ja2ind, *ja2bs;

	for (i=0;i<_mtra.nblksc;i++) {

		cout  << " Scale Iblk = " << i << endl;

// Read current block

		_mtra.ReadBlock (i);

// Compute indices

		gmtrdummy.IndicesForScaling (i,
												_mtrasp, _mtra, 
												jamask, sp2blk, icycle, imask,
												at, indarr,
												mtraoff, 
												ja2isup, ja2blk, ja2ind, ja2bs);

// Read off-block data

		gmtrdummy.ReadBlockRows (_mtra, mtraoff, 
											ja2isup, ja2blk, ja2ind, ja2bs,
											icycle, imask);

// Filter block rows

		gmtrdummy.ScaleBlockRows (i, // Filter some block rows
											_mtrasp, at, indarr,
											_mtra, jamask, sp2blk,
											mtraoff, _rhs,
											ja2isup, ja2blk, ja2ind, ja2bs);

// Write current block

		_mtra.ReWriteBlock (i);
		_mtra.FreeBlock (i);

// Write off-block data

		gmtrdummy.WriteBlockRows (_mtra, mtraoff, 
											ja2isup, ja2blk, ja2ind, ja2bs,
											icycle, imask);

// Free work arrays

		delete [] ja2isup;
		delete [] ja2blk;
		delete [] ja2ind;
		delete [] ja2bs;

	};

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout << " Scaling block rows statistics: " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Scaling block rows statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work arrays

	delete [] sp2blk;
	delete [] indarr;
	delete [] imask;
	delete [] jamask;

};

// Author: Kharchenko S.A. 
// Description: Perform in-place filtering of the block rows of the rectangular matrix
// CGSMatrixCS::FilterBlockRows()
//========================================================================================
void CGSMatrixCS::FilterBlockRows (ofstream &_fout, const CSlvParam &_param, // Perform in-place filtering of the block rows of the rectangular matrix
												CSMatrix &_mtrasp, CGSMatrixCS &_mtra, CSVectorC &_rhs,
												int *_iamask, int *_jamask) {

	const char *funcname = "FilterBlockRows";

	dcmplx czero (0.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtra.nsupc];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;

	for (i=0;i<_mtra.nblksc;i++) {
		for (j=_mtra.blksc[i];j<_mtra.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Allocate index array

	int *indarr;

	indarr = new int [_mtrasp.nzja];
	if (!indarr) MemoryFail (funcname);

// Compute transposed matrix

	CSMatrix at;

	at = _mtrasp.TranspMtrRect (indarr);

// Init the mask

	int icycle = -1;

	int *imask;

	imask = new int [_mtra.nsupc];
	if (!imask) MemoryFail (funcname);

	for (i=0;i<_mtra.nsupc;i++) imask[i] = icycle;

// Perform filtering block by block

	CGSMatrixCS gmtrdummy;

	CSMatrixCS mtraoff;

	int *ja2isup, *ja2blk, *ja2ind, *ja2bs;

	for (i=0;i<_mtra.nblksc;i++) {

		cout  << " Filter Iblk = " << i << endl;

// Read current block

		_mtra.ReadBlock (i);

// Compute indices

		gmtrdummy.IndicesForFilteringOrScaling (i,
													_mtrasp, _mtra, 
													_jamask, sp2blk, icycle, imask,
													at, indarr,
													mtraoff, 
													ja2isup, ja2blk, ja2ind, ja2bs);

// Read off-block data

		gmtrdummy.ReadBlockRows (_mtra, mtraoff, 
											ja2isup, ja2blk, ja2ind, ja2bs,
											icycle, imask);

// Filter block rows

		gmtrdummy.FilterBlockRows (i, // Filter some block rows
											_mtrasp, at, indarr,
											_mtra, _jamask, sp2blk,
											mtraoff, _rhs,
											ja2isup, ja2blk, ja2ind, ja2bs);

// Write current block

		_mtra.ReWriteBlock (i);
		_mtra.FreeBlock (i);

// Write off-block data

		gmtrdummy.WriteBlockRows (_mtra, mtraoff, 
											ja2isup, ja2blk, ja2ind, ja2bs,
											icycle, imask);

// Free work arrays

		delete [] ja2isup;
		delete [] ja2blk;
		delete [] ja2ind;
		delete [] ja2bs;

	};

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout << " Filtering block rows statistics: " << endl;
		_fout << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Filtering block rows statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work arrays

	delete [] sp2blk;
	delete [] indarr;
	delete [] imask;

};

// Author: Kharchenko S.A. 
// Description: Compute indices of the nonlocal data for block rows filtering or scaling
// CGSMatrixCS::IndicesForFilteringOrScaling()
//========================================================================================
void CGSMatrixCS::IndicesForFilteringOrScaling (int _iblk, // Compute indices of the nonlocal data for block rows filtering or scaling
													CSMatrix &_mtrasp, CGSMatrixCS &_mtra, 
													int *_jamask, int *_sp2blk,
													int &_icycle, int *_imask,
													CSMatrix &_mtrasptran, int *_indarr,
													CSMatrixCS &_mtraoff, 
													int *&_ja2isup, int *&_ja2blk, int *&_ja2ind, int *&_ja2bs) {

	const char *funcname = "IndicesForFilteringOrScaling";

// Initial scan of the current block row and estimate memory for the matrix

	int isup, j, jsup, k, ksup, kblk;
	int blai, blaj, blak;

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;

	int nzjaloc = 0;
	int nzaloc = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		for (j=_mtrasp.ia[isup];j<_mtrasp.ia[isup+1];j++) {
			jsup = _mtrasp.ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (_jamask[j] > 0) {
				for (k=_mtrasptran.ia[jsup];k<_mtrasptran.ia[jsup+1];k++) {
					ksup = _mtrasptran.ja[k];
					blak = _mtra.sprndc[ksup+1]-_mtra.sprndc[ksup];
					kblk = _sp2blk[ksup];
					if (kblk != _iblk) {
						nzjaloc++;
						nzaloc += blaj*blak;
					};
				};
			};
		};
	};

// Allocate the memory for int arrays

	_ja2isup = new int [nzjaloc];
	if (!_ja2isup) MemoryFail (funcname);
	_ja2blk = new int [nzjaloc];
	if (!_ja2blk) MemoryFail (funcname);
	_ja2ind = new int [nzjaloc];
	if (!_ja2ind) MemoryFail (funcname);
	_ja2bs = new int [nzjaloc];
	if (!_ja2bs) MemoryFail (funcname);

// Fill int arrays

	nzjaloc = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		for (j=_mtrasp.ia[isup];j<_mtrasp.ia[isup+1];j++) {
			jsup = _mtrasp.ja[j];
			if (_jamask[j] > 0) {
				for (k=_mtrasptran.ia[jsup];k<_mtrasptran.ia[jsup+1];k++) {
					ksup = _mtrasptran.ja[k];
					kblk = _sp2blk[ksup];
					if (kblk != _iblk) {
						_ja2blk[nzjaloc] = kblk;
						_ja2ind[nzjaloc] = ksup;
						_ja2bs [nzjaloc] = jsup;
						nzjaloc++;
					};
				};
			};
		};
	};

// Create and sort the list of nonzero supernode columns in the matrix

	int *listloc, *imaskloc;

	int i;

	int nlistloc = 0;

	listloc = new int [_mtrasp.nsupc];
	if (!listloc) MemoryFail (funcname);
	imaskloc = new int [_mtrasp.nsupc];
	if (!imaskloc) MemoryFail (funcname);

	_icycle++;

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		if (_imask[ksup] != _icycle) {
			listloc[nlistloc] = ksup;
			nlistloc++;
			_imask[ksup] = _icycle;
		};
	};

	if (nlistloc != 0) qsort (listloc, nlistloc, sizeof(int), compint);

	for (i=0;i<nlistloc;i++) {
		isup = listloc[i];
		imaskloc[isup] = i;
	};

// Create the matrix except the values

	CSMatrixCS aloc (nlistloc, nzjaloc, nzaloc);

	for (i=0;i< nlistloc;i++) aloc.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) aloc.sprndr[i] = 0;
	for (i=0;i<=nlistloc;i++) aloc.sprndc[i] = 0;

	for (i=0;i<=nlistloc;i++) aloc.ia[i] = 0;

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		j = imaskloc[ksup];
		aloc.ia[j+1]++;
	};

	for (i=0;i<nlistloc;i++) aloc.ia[i+1] = aloc.ia[i] + aloc.ia[i+1];

	for (i=0;i<nlistloc;i++) listloc[i] = aloc.ia[i];

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		jsup = _ja2bs[i];
		j = imaskloc[ksup];
		k = listloc[j];
		aloc.ja[k] = jsup;
		listloc[j]++;
	};

	int ibs = 0;

	for (i=0;i<nlistloc;i++) {
		isup = aloc.list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jsup = aloc.ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			aloc.bsa[j] = ibs;
			ibs += blai*blaj;
		};
	};

	aloc.nsupc  = nlistloc;
	aloc.nsupr  = nlistloc;
	aloc.nlist  = nlistloc;
	aloc.nzja   = nzjaloc;
	aloc.nza    = nzaloc;
	aloc.nzatot = nzaloc;

	_mtraoff = aloc;

// Fill int arrays

	int iloc, i1beg, i2beg, i1end, i2end, i1ind, i2ind, jsup1, jsup2;

	for (i=0;i<nlistloc;i++) {
		isup = aloc.list[i];
		kblk = _sp2blk[isup];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			_ja2isup[j] = isup;
			_ja2blk[j] = kblk;
		};
		isupbeg = _mtra.blksc[kblk];
		iloc = isup-isupbeg;
		i1beg = aloc.ia[i];
		i1end = aloc.ia[i+1]-1;
		i2beg = _mtra.mtrarr[kblk].ia[iloc];
		i2end = _mtra.mtrarr[kblk].ia[iloc+1]-1;
		i1ind = i1beg;
		i2ind = i2beg;
		while (i1ind <= i1end && i2ind <= i2end) {
			jsup1 = aloc.ja[i1ind];
			jsup2 = _mtra.mtrarr[kblk].ja[i2ind];
			if (jsup1 == jsup2) {
				_ja2ind[i1ind] = i2ind;
				_ja2bs [i1ind] = _mtra.mtrarr[kblk].bsa[i2ind];
				i1ind++;
				i2ind++;
			} else if (jsup1 < jsup2) {
				i1ind++;
			} else if (jsup1 > jsup2) {
				i2ind++;
			};
		};
	};

// Free work arrays

	delete [] listloc;
	delete [] imaskloc;

};

// Author: Kharchenko S.A. 
// Description: Compute indices of the nonlocal data for block rows scaling
// CGSMatrixCS::IndicesForScaling()
//========================================================================================
void CGSMatrixCS::IndicesForScaling (int _iblk, // Compute indices of the nonlocal data for block rows scaling
													CSMatrix &_mtrasp, CGSMatrixCS &_mtra, 
													int *_jamask, int *_sp2blk,
													int &_icycle, int *_imask,
													CSMatrix &_mtrasptran, int *_indarr,
													CSMatrixCS &_mtraoff, 
													int *&_ja2isup, int *&_ja2blk, int *&_ja2ind, int *&_ja2bs) {

	const char *funcname = "IndicesForScaling";

// Initial scan of the current block row and estimate memory for the matrix

	int isup, j, jsup, k, ksup, kblk;
	int blai, blaj;

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;

	int nzjaloc = 0;
	int nzaloc = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
		for (j=_mtrasptran.ia[isup];j<_mtrasptran.ia[isup+1];j++) {
			jsup = _mtrasptran.ja[j];
			blaj = _mtra.sprndc[jsup+1]-_mtra.sprndc[jsup];
			kblk = _sp2blk[jsup];
			if (kblk != _iblk) {
				nzjaloc++;
				nzaloc += blaj*blai;
			};
		};
	};

// Allocate the memory for int arrays

	_ja2isup = new int [nzjaloc];
	if (!_ja2isup) MemoryFail (funcname);
	_ja2blk = new int [nzjaloc];
	if (!_ja2blk) MemoryFail (funcname);
	_ja2ind = new int [nzjaloc];
	if (!_ja2ind) MemoryFail (funcname);
	_ja2bs = new int [nzjaloc];
	if (!_ja2bs) MemoryFail (funcname);

// Fill int arrays

	nzjaloc = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
		for (j=_mtrasptran.ia[isup];j<_mtrasptran.ia[isup+1];j++) {
			jsup = _mtrasptran.ja[j];
			blaj = _mtra.sprndc[jsup+1]-_mtra.sprndc[jsup];
			kblk = _sp2blk[jsup];
			if (kblk != _iblk) {
				_ja2blk[nzjaloc] = kblk;
				_ja2ind[nzjaloc] = jsup;
				_ja2bs [nzjaloc] = isup;
				nzjaloc++;
			};
		};
	};

// Create and sort the list of nonzero supernode columns in the matrix

	int *listloc, *imaskloc;

	int i;

	int nlistloc = 0;

	listloc = new int [_mtrasp.nsupc];
	if (!listloc) MemoryFail (funcname);
	imaskloc = new int [_mtrasp.nsupc];
	if (!imaskloc) MemoryFail (funcname);

	_icycle++;

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		if (_imask[ksup] != _icycle) {
			listloc[nlistloc] = ksup;
			nlistloc++;
			_imask[ksup] = _icycle;
		};
	};

	if (nlistloc != 0) qsort (listloc, nlistloc, sizeof(int), compint);

	for (i=0;i<nlistloc;i++) {
		isup = listloc[i];
		imaskloc[isup] = i;
	};

// Create the matrix except the values

	CSMatrixCS aloc (nlistloc, nzjaloc, nzaloc);

	for (i=0;i< nlistloc;i++) aloc.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) aloc.sprndr[i] = 0;
	for (i=0;i<=nlistloc;i++) aloc.sprndc[i] = 0;

	for (i=0;i<=nlistloc;i++) aloc.ia[i] = 0;

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		j = imaskloc[ksup];
		aloc.ia[j+1]++;
	};

	for (i=0;i<nlistloc;i++) aloc.ia[i+1] = aloc.ia[i] + aloc.ia[i+1];

	for (i=0;i<nlistloc;i++) listloc[i] = aloc.ia[i];

	for (i=0;i<nzjaloc;i++) {
		ksup = _ja2ind[i];
		jsup = _ja2bs[i];
		j = imaskloc[ksup];
		k = listloc[j];
		aloc.ja[k] = jsup;
		listloc[j]++;
	};

	int ibs = 0;

	for (i=0;i<nlistloc;i++) {
		isup = aloc.list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jsup = aloc.ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			aloc.bsa[j] = ibs;
			ibs += blai*blaj;
		};
	};

	aloc.nsupc  = nlistloc;
	aloc.nsupr  = nlistloc;
	aloc.nlist  = nlistloc;
	aloc.nzja   = nzjaloc;
	aloc.nza    = nzaloc;
	aloc.nzatot = nzaloc;

	_mtraoff = aloc;

// Fill int arrays

	int iloc, i1beg, i2beg, i1end, i2end, i1ind, i2ind, jsup1, jsup2;

	for (i=0;i<nlistloc;i++) {
		isup = aloc.list[i];
		kblk = _sp2blk[isup];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			_ja2isup[j] = isup;
			_ja2blk[j] = kblk;
		};
		isupbeg = _mtra.blksc[kblk];
		iloc = isup-isupbeg;
		i1beg = aloc.ia[i];
		i1end = aloc.ia[i+1]-1;
		i2beg = _mtra.mtrarr[kblk].ia[iloc];
		i2end = _mtra.mtrarr[kblk].ia[iloc+1]-1;
		i1ind = i1beg;
		i2ind = i2beg;
		while (i1ind <= i1end && i2ind <= i2end) {
			jsup1 = aloc.ja[i1ind];
			jsup2 = _mtra.mtrarr[kblk].ja[i2ind];
			if (jsup1 == jsup2) {
				_ja2ind[i1ind] = i2ind;
				_ja2bs [i1ind] = _mtra.mtrarr[kblk].bsa[i2ind];
				i1ind++;
				i2ind++;
			} else if (jsup1 < jsup2) {
				i1ind++;
			} else if (jsup1 > jsup2) {
				i2ind++;
			};
		};
	};

// Free work arrays

	delete [] listloc;
	delete [] imaskloc;

};

// Author: Kharchenko S.A. 
// Description: Read some block rows
// CGSMatrixCS::ReadBlockRows()
//========================================================================================
void CGSMatrixCS::ReadBlockRows ( // Read some block rows
											CGSMatrixCS &_mtra, 
											CSMatrixCS &_mtraoff, 
											int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs,
											int &_icycle, int *_imask) {

	const char *funcname = "ReadBlockRows";

// Count the number and list of blocks to be read

	int *listloc;

	int i, j, jblk;

	int nblkloc = 0;

	listloc = new int [_mtra.nblksc];
	if (!listloc) MemoryFail (funcname);

	_icycle++;

	for (j=0;j<_mtraoff.nzja;j++) {
		jblk = _ja2blk[j];
		if (_imask[jblk] != _icycle) {
			listloc[nblkloc] = jblk;
			nblkloc++;
			_imask[jblk] = _icycle;
		};
	};

	if (nblkloc != 0) qsort (listloc, nblkloc, sizeof(int), compint);

// Read the data block by block and fill the matrix

	int kblk, isup, jsup, blai, blaj, k, ibs, jbs;

	int nzjaloc = 0;

	for (i=0;i<nblkloc;i++) {
		jblk = listloc[i];
		_mtra.ReadBlock (jblk);
		if (nzjaloc < _mtraoff.nzja) {
			kblk = _ja2blk[nzjaloc];
		} else {
			kblk = -1;
		};
		while (nzjaloc < _mtraoff.nzja && kblk == jblk) {
			isup = _ja2isup[nzjaloc];
			jsup = _mtraoff.ja[nzjaloc];
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			ibs = _mtraoff.bsa[nzjaloc];
			jbs = _ja2bs[nzjaloc];
			for (k=0;k<blai*blaj;k++) {
				_mtraoff.a[ibs+k] = _mtra.mtrarr[jblk].a[jbs+k];
			};
			nzjaloc++;
			if (nzjaloc < _mtraoff.nzja) {
				kblk = _ja2blk[nzjaloc];
			} else {
				kblk = -1;
			};
		};
		_mtra.FreeBlock (jblk);
	};

// Free work arrays

	delete [] listloc;

};

// Author: Kharchenko S.A. 
// Description: Write some block rows
// CGSMatrixCS::WriteBlockRows()
//========================================================================================
void CGSMatrixCS::WriteBlockRows ( // Write some block rows
											CGSMatrixCS &_mtra, 
											CSMatrixCS &_mtraoff, 
											int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs,
											int &_icycle, int *_imask) {

	const char *funcname = "WriteBlockRows";

// Count the number and list of blocks to be written

	int *listloc;

	int i, j, jblk;

	int nblkloc = 0;

	listloc = new int [_mtra.nblksc];
	if (!listloc) MemoryFail (funcname);

	_icycle++;

	for (j=0;j<_mtraoff.nzja;j++) {
		jblk = _ja2blk[j];
		if (_imask[jblk] != _icycle) {
			listloc[nblkloc] = jblk;
			nblkloc++;
			_imask[jblk] = _icycle;
		};
	};

	if (nblkloc != 0) qsort (listloc, nblkloc, sizeof(int), compint);

// Write the data block by block and fill the matrix

	int kblk, isup, jsup, blai, blaj, k, ibs, jbs;

	int nzjaloc = 0;

	for (i=0;i<nblkloc;i++) {
		jblk = listloc[i];
		_mtra.ReadBlock (jblk);
		if (nzjaloc < _mtraoff.nzja) {
			kblk = _ja2blk[nzjaloc];
		} else {
			kblk = -1;
		};
		while (nzjaloc < _mtraoff.nzja && kblk == jblk) {
			isup = _ja2isup[nzjaloc];
			jsup = _mtraoff.ja[nzjaloc];
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			ibs = _mtraoff.bsa[nzjaloc];
			jbs = _ja2bs[nzjaloc];
			for (k=0;k<blai*blaj;k++) {
				_mtra.mtrarr[jblk].a[jbs+k] = _mtraoff.a[ibs+k];
			};
			nzjaloc++;
			if (nzjaloc < _mtraoff.nzja) {
				kblk = _ja2blk[nzjaloc];
			} else {
				kblk = -1;
			};
		};
		_mtra.ReWriteBlock (jblk);
		_mtra.FreeBlock (jblk);
	};

// Free work arrays

	delete [] listloc;

};

// Author: Kharchenko S.A. 
// Description: Filter some block rows
// CGSMatrixCS::FilterBlockRows()
//========================================================================================
void CGSMatrixCS::FilterBlockRows (int _iblk, // Filter some block rows
												CSMatrix &_mtrasp, CSMatrix &_mtrasptran, int *_indarr,
												CGSMatrixCS &_mtra, int *_jamask, int *_sp2blk,
												CSMatrixCS &_mtraoff, CSVectorC &_rhs,
												int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs) {

	const char *funcname = "FilterBlockRows";

	dcmplx czero (0.0e0,0.0e0);

// Determine the maximal size

	int i, j, isup, jsup, blai, blaj;

	int nlocmax = 0;

	for (i=0;i<_mtra.mtrarr[_iblk].nlist;i++) {
		isup = _mtra.mtrarr[_iblk].list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtra.mtrarr[_iblk].ia[i];j<_mtra.mtrarr[_iblk].ia[i+1];j++) {
			jsup = _mtra.mtrarr[_iblk].ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};
	for (i=0;i<_mtraoff.nlist;i++) {
		isup = _mtraoff.list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtraoff.ia[i];j<_mtraoff.ia[i+1];j++) {
			jsup = _mtraoff.ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};

// Allocate work arrays

	dcmplx *aloc, *bloc, *cloc, *uloc, *vloc, *work;
	double *svloc, *sclloc, *rwork;
	int lwork = nlocmax * 10;
	int info, k;

	aloc = new dcmplx [nlocmax*nlocmax];
	if (!aloc) MemoryFail (funcname);
	bloc = new dcmplx [nlocmax*nlocmax];
	if (!bloc) MemoryFail (funcname);
	cloc = new dcmplx [nlocmax*nlocmax];
	if (!cloc) MemoryFail (funcname);
	uloc = new dcmplx [nlocmax*nlocmax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [nlocmax*nlocmax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [nlocmax];
	if (!svloc) MemoryFail (funcname);
	sclloc = new double [nlocmax];
	if (!sclloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Create mask array

	int *imaskloc;

	imaskloc = new int [_mtra.nsupc];
	if (!imaskloc) MemoryFail (funcname);

	for (i=0;i<_mtraoff.nlist;i++) {
		isup = _mtraoff.list[i];
		imaskloc[isup] = i;
	};

// Main cycle over the supernode columns

	int isuploc, niloc, njloc, ibs, ibsv;
	int kii, kjj, jbs, j1, j2, ind1, ind2, ksup, kksup, jloc, kblk;
	int iiloc, jjloc, blak, kkk;
	double aux;
	dcmplx caux, caux1;
	bool bpart, boff;

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;

	int ishift = _mtrasp.ia[isupbeg];

	for (isup=isupbeg;isup<=isupend;isup++) {
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		isuploc = isup-isupbeg;
		for (j=_mtra.mtrarr[_iblk].ia[isuploc];j<_mtra.mtrarr[_iblk].ia[isuploc+1];j++) {
			jsup = _mtra.mtrarr[_iblk].ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			niloc = blaj;
			ibs = _mtra.mtrarr[_iblk].bsa[j];
			if (_jamask[j+ishift] > 0) {

// Read current blocks data

				if (blai != blaj) {
					throw " CGSMatrixCS::FilterBlockRows: wrong block sizes ";
				};

				for (kii=0;kii<niloc;kii++) {
					for (kjj=0;kjj<niloc;kjj++) {
						aloc[kjj*niloc+kii] = _mtra.mtrarr[_iblk].a[ibs+kii*niloc+kjj];
					};
				};

				if (_mtrasptran.ia[jsup+1]-_mtrasptran.ia[jsup] == 1) {
					njloc = -1;
					jbs = -1;
					bpart = false;
					boff = false;
				} else {
					bpart = true;
					if (_mtrasptran.ia[jsup+1]-_mtrasptran.ia[jsup] != 2) {
						throw " CGSMatrixCS::FilterBlockRows: wrong number of nonzero blocks in block row ";
					};
					j1 = _mtrasptran.ia[jsup];
					j2 = _mtrasptran.ia[jsup]+1;
					ind1 = _mtrasptran.ja[j1];
					ind2 = _mtrasptran.ja[j2];
					if (ind1 == isup) {
						ksup = ind2;
						jloc = j2;
					} else {
						ksup = ind1;
						jloc = j1;
					};
					kblk = _sp2blk[ksup];
					if (kblk != _iblk) {
						boff = true;
						iiloc = imaskloc[ksup];
						jbs = -1;
						for (k=_mtraoff.ia[iiloc];k<_mtraoff.ia[iiloc+1];k++) {
							kksup = _mtraoff.ja[k];
							if (kksup == jsup) {
								jbs = _mtraoff.bsa[k];
							};
						};
						if (jbs == -1) {
							throw " CGSMatrixCS::FilterBlockRows: wrong base address ";
						};
					} else {
						boff = false;
						jjloc = _indarr[jloc];
						jbs = _mtra.mtrarr[_iblk].bsa[jjloc-ishift];
					};
					blak = _mtra.sprndc[ksup+1]-_mtra.sprndc[ksup];
					njloc = blak;
					if (boff) {
						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								bloc[kjj*niloc+kii] = _mtraoff.a[jbs+kii*njloc+kjj];
							};
						};
					} else {
						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								bloc[kjj*niloc+kii] = _mtra.mtrarr[_iblk].a[jbs+kii*njloc+kjj];
							};
						};
					};
				};

// Filter current block row

// Via SVD

				if (false) {

// Compute SVD

					zgesvd_ ("A", "A", &niloc, &niloc,
								aloc, &niloc, svloc, 
								uloc, &niloc, vloc, &niloc, 
								work, &lwork, rwork, &info);

					if (info != 0) throw " Error in the Lapack routine ZGESVD";

// Compute projector size and scaling array

					double thresh = 1.0e-20;

					double svmax = svloc[0];

					int nsv = 0;

					for (kii=0;kii<niloc;kii++) {
						if (svloc[kii] >= svmax * thresh) nsv++;
	//					aux = svloc[kii] / svmax;
						aux = svloc[kii];
	//					sclloc[kii] = sqrt(aux);
						sclloc[kii] = aux;
					};

// Multiply from the left and zero out remaining part

// A part

					for (kii=0;kii<nsv;kii++) {
						aux = svloc[kii] * sclloc[kii];
						for (kjj=0;kjj<niloc;kjj++) {
							caux = vloc[kjj*niloc+kii];
							aloc[kjj*niloc+kii] = caux * aux;
						};
					};
					for (kii=nsv;kii<niloc;kii++) {
						for (kjj=0;kjj<niloc;kjj++) {
							aloc[kjj*niloc+kii] = czero;
						};
					};

// Modify U

					for (kii=0;kii<nsv;kii++) {
						aux = sclloc[kii];
						for (kjj=0;kjj<niloc;kjj++) {
							caux = uloc[kii*niloc+kjj];
							uloc[kii*niloc+kjj] = caux * aux;
						};
					};

// Modify B part if necessary

					if (bpart) {

						for (kii=0;kii<nsv;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								caux = czero;
								for (kkk=0;kkk<niloc;kkk++) {
									caux1 = conj(uloc[kii*niloc+kkk]) * bloc[kjj*niloc+kkk];
									caux += caux1;
								};
								cloc[kjj*niloc+kii] = caux;
							};
						};

						for (kii=nsv;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								cloc[kjj*niloc+kii] = czero;
							};
						};

						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								bloc[kjj*niloc+kii] = cloc[kjj*niloc+kii];
							};
						};

					};

// Rhs part

					ibsv = _mtra.sprndr[jsup];

					for (kii=0;kii<nsv;kii++) {
						caux = czero;
						for (kkk=0;kkk<niloc;kkk++) {
							caux1 = conj(uloc[kii*niloc+kkk]) * _rhs.vect[ibsv+kkk];
							caux += caux1;
						};
						cloc[kii] = caux;
					};

					for (kii=nsv;kii<niloc;kii++) {
						cloc[kii] = czero;
					};

					for (kii=0;kii<niloc;kii++) {
						_rhs.vect[ibsv+kii] = cloc[kii];
					};

				} else {

// Amatr

					int iind, jind;

					for (kii=0;kii<niloc*niloc;kii++) {
						uloc[kii] = conj(aloc[kii]);
					};

					for (kii=0;kii<niloc;kii++) {
						for (kjj=0;kjj<niloc;kjj++) {
							caux = czero;
							iind = kii*niloc;
							jind = kjj*niloc;
							for (kkk=0;kkk<niloc;kkk++) {
//								caux1 = uloc[kii*niloc+kkk] * aloc[kjj*niloc+kkk];
								caux1 = uloc[iind++] * aloc[jind++];
								caux += caux1;
							};
							cloc[kjj*niloc+kii] = caux;
						};
					};

					for (kii=0;kii<niloc*niloc;kii++) aloc[kii] = cloc[kii];

// Modify B part if necessary

					if (bpart) {

						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								caux = czero;
								iind = kii*niloc;
								jind = kjj*niloc;
								for (kkk=0;kkk<niloc;kkk++) {
//									caux1 = uloc[kii*niloc+kkk] * bloc[kjj*niloc+kkk];
									caux1 = uloc[iind++] * bloc[jind++];
									caux += caux1;
								};
								cloc[kjj*niloc+kii] = caux;
							};
						};

						for (kii=0;kii<niloc*njloc;kii++) bloc[kii] = cloc[kii];

					};

// Modify rhs

					ibsv = _mtra.sprndr[jsup];

					for (kii=0;kii<niloc;kii++) {
						caux = czero;
						for (kkk=0;kkk<niloc;kkk++) {
							caux1 = uloc[kii*niloc+kkk] * _rhs.vect[ibsv+kkk];
							caux += caux1;
						};
						cloc[kii] = caux;
					};

					for (kii=0;kii<niloc;kii++) {
						_rhs.vect[ibsv+kii] = cloc[kii];
					};

				};

// Write current blocks data

				for (kii=0;kii<niloc;kii++) {
					for (kjj=0;kjj<niloc;kjj++) {
						_mtra.mtrarr[_iblk].a[ibs+kii*niloc+kjj] = aloc[kjj*niloc+kii];
					};
				};

				if (bpart) {
					if (boff) {
						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								_mtraoff.a[jbs+kii*njloc+kjj] = bloc[kjj*niloc+kii];
							};
						};
					} else {
						for (kii=0;kii<niloc;kii++) {
							for (kjj=0;kjj<njloc;kjj++) {
								_mtra.mtrarr[_iblk].a[jbs+kii*njloc+kjj] = bloc[kjj*niloc+kii];
							};
						};
					};
				};

			};
		};
	};

// Free work arrays

	delete [] aloc;
	delete [] bloc;
	delete [] cloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] sclloc;
	delete [] work;
	delete [] rwork;
	delete [] imaskloc;

};

// Author: Kharchenko S.A. 
// Description: Scale some block rows
// CGSMatrixCS::ScaleBlockRows()
//========================================================================================
void CGSMatrixCS::ScaleBlockRows (int _iblk, // Scale some block rows
												CSMatrix &_mtrasp, CSMatrix &_mtrasptran, int *_indarr,
												CGSMatrixCS &_mtra, int *_jamask, int *_sp2blk,
												CSMatrixCS &_mtraoff, CSVectorC &_rhs,
												int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs) {

	const char *funcname = "ScaleBlockRows";

	dcmplx czero (0.0e0,0.0e0);

// Determine the maximal size

	int i, j, isup, jsup, blai, blaj, niloc;

	int nlocmax = 0;

	for (i=0;i<_mtra.mtrarr[_iblk].nlist;i++) {
		isup = _mtra.mtrarr[_iblk].list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtra.mtrarr[_iblk].ia[i];j<_mtra.mtrarr[_iblk].ia[i+1];j++) {
			jsup = _mtra.mtrarr[_iblk].ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};
	for (i=0;i<_mtraoff.nlist;i++) {
		isup = _mtraoff.list[i];
		blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
		if (blai > nlocmax) nlocmax = blai;
		for (j=_mtraoff.ia[i];j<_mtraoff.ia[i+1];j++) {
			jsup = _mtraoff.ja[j];
			blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
			if (blaj > nlocmax) nlocmax = blaj;
		};
	};

	int isupbeg = _mtra.blksc[_iblk];
	int isupend = _mtra.blksc[_iblk+1]-1;

	int nilocmax = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {
		niloc = 0;
		for (j=_mtrasptran.ia[isup];j<_mtrasptran.ia[isup+1];j++) {
			jsup = _mtrasptran.ja[j];
			blaj = _mtra.sprndc[jsup+1]-_mtra.sprndc[jsup];
			niloc += blaj;
		};
		if (nilocmax < niloc) nilocmax = niloc;
	};

// Allocate work arrays

	dcmplx *aloc, *qloc, *tauloc, *cloc, *uloc, *vloc, *work;
	double *svloc, *rwork;
	int lwork = nlocmax * 10;
	int info;

	aloc = new dcmplx [nlocmax*nilocmax];
	if (!aloc) MemoryFail (funcname);
	qloc = new dcmplx [nlocmax*nilocmax];
	if (!qloc) MemoryFail (funcname);
	tauloc = new dcmplx [nlocmax];
	if (!tauloc) MemoryFail (funcname);
	cloc = new dcmplx [nlocmax*nlocmax];
	if (!cloc) MemoryFail (funcname);
	uloc = new dcmplx [nlocmax*nlocmax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [nlocmax*nlocmax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [nlocmax];
	if (!svloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Create mask array

	int *imaskloc;

	imaskloc = new int [_mtra.nsupc];
	if (!imaskloc) MemoryFail (funcname);

	for (i=0;i<_mtraoff.nlist;i++) {
		isup = _mtraoff.list[i];
		imaskloc[isup] = i;
	};

// Main cycle over the supernode rows

	int ibs, jloc;
	int kii, kjj;
	int kkk, ishift, jsuploc, isup1;
	double daux;
	dcmplx caux;

	for (isup=isupbeg;isup<=isupend;isup++) {

		blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
		ishift = _mtra.sprndr[isup];

// Read current blocks data

		niloc = 0;

		for (j=_mtrasptran.ia[isup];j<_mtrasptran.ia[isup+1];j++) {
			jsup = _mtrasptran.ja[j];
			blaj = _mtra.sprndc[jsup+1]-_mtra.sprndc[jsup];

			if (jsup >= isupbeg && jsup <= isupend) {

				jsuploc = jsup-isupbeg;

				ibs = -1;

				for (jloc=_mtra.mtrarr[_iblk].ia[jsuploc];jloc<_mtra.mtrarr[_iblk].ia[jsuploc+1];jloc++) {
					isup1 = _mtra.mtrarr[_iblk].ja[jloc];
					if (isup1 == isup) {
						ibs = _mtra.mtrarr[_iblk].bsa[jloc];
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								aloc[niloc*blai+kjj*blai+kii] = _mtra.mtrarr[_iblk].a[ibs+kii*blaj+kjj];
							};
						};
					};
				};

				if (ibs < 0) throw " CGSMatrixCS::ScaleBlockRows: index not found ";

			} else {

				jsuploc = imaskloc[jsup];

				ibs = -1;

				for (jloc=_mtraoff.ia[jsuploc];jloc<_mtraoff.ia[jsuploc+1];jloc++) {
					isup1 = _mtraoff.ja[jloc];
					if (isup1 == isup) {
						ibs = _mtraoff.bsa[jloc];
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								aloc[niloc*blai+kjj*blai+kii] = _mtraoff.a[ibs+kii*blaj+kjj];
							};
						};
					};
				};

				if (ibs < 0) throw " CGSMatrixCS::ScaleBlockRows: index not found ";

			};

			niloc += blaj;

		};

// Hermitian transpose the data

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<niloc;kjj++) {
				qloc[kii*niloc+kjj] = conj (aloc[kjj*blai+kii]);
			};
		};

// Compute QR decomposition

		QrdBlk (niloc, blai, qloc, niloc, tauloc);

// Take R matrix

		GetRPartQrd (niloc, 0, blai,
							qloc, niloc, cloc, blai);

// Compute SVD

		zgesvd_ ("A", "A", &blai, &blai,
					cloc, &blai, svloc, 
					uloc, &blai, vloc, &blai, 
					work, &lwork, rwork, &info);

		if (info != 0) throw " Error in the Lapack routine ZGESVD";

//		hyst?;

// Compute scaling matrix

		for (kii=0;kii<blai;kii++) {
			daux = 1.0e0 / sqrt(svloc[kii]);
			for (kjj=0;kjj<blai;kjj++) {
				caux = uloc[kii*blai+kjj];
				uloc[kii*blai+kjj] = caux * daux;
			};
		};

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
					caux += uloc[kkk*blai+kii] * vloc[kjj*blai+kkk];
				};
				cloc[kjj*blai+kii] = caux;
			};
		};

// Scale the matrix data

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<niloc;kjj++) {
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
					caux += cloc[kkk*blai+kii] * aloc[kjj*blai+kkk];
				};
				qloc[kjj*blai+kii] = caux;
			};
		};

// Scale rhs data

		for (kii=0;kii<blai;kii++) {
			caux = czero;
			for (kjj=0;kjj<blai;kjj++) {
				caux += cloc[kjj*blai+kii] * _rhs.vect[ishift+kjj];
			};
			tauloc[kii] = caux;
		};

		for (kii=0;kii<blai;kii++) _rhs.vect[ishift+kii] = tauloc[kii];

// Rewrite the data back

		niloc = 0;

		for (j=_mtrasptran.ia[isup];j<_mtrasptran.ia[isup+1];j++) {
			jsup = _mtrasptran.ja[j];
			blaj = _mtra.sprndc[jsup+1]-_mtra.sprndc[jsup];

			if (jsup >= isupbeg && jsup <= isupend) {

				jsuploc = jsup-isupbeg;

				ibs = -1;

				for (jloc=_mtra.mtrarr[_iblk].ia[jsuploc];jloc<_mtra.mtrarr[_iblk].ia[jsuploc+1];jloc++) {
					isup1 = _mtra.mtrarr[_iblk].ja[jloc];
					if (isup1 == isup) {
						ibs = _mtra.mtrarr[_iblk].bsa[jloc];
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								_mtra.mtrarr[_iblk].a[ibs+kii*blaj+kjj] = qloc[niloc*blai+kjj*blai+kii];
							};
						};
					};
				};

				if (ibs < 0) throw " CGSMatrixCS::ScaleBlockRows: index not found ";

			} else {

				jsuploc = imaskloc[jsup];

				ibs = -1;

				for (jloc=_mtraoff.ia[jsuploc];jloc<_mtraoff.ia[jsuploc+1];jloc++) {
					isup1 = _mtraoff.ja[jloc];
					if (isup1 == isup) {
						ibs = _mtraoff.bsa[jloc];
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								_mtraoff.a[ibs+kii*blaj+kjj] = qloc[niloc*blai+kjj*blai+kii];
							};
						};
					};
				};

				if (ibs < 0) throw " CGSMatrixCS::ScaleBlockRows: index not found ";

			};

			niloc += blaj;

		};

	};

// Free work arrays

	delete [] aloc;
	delete [] qloc;
	delete [] tauloc;
	delete [] cloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;
	delete [] imaskloc;

};

// Author: Kharchenko S.A. 
// Description: Compute point scaling of the rectangular matrix A and scale explicitely
// CGSMatrixCS::ExplicitSymmetricPointScaling()
//========================================================================================
void CGSMatrixCS::ExplicitSymmetricPointScaling (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute point scaling of the rectangular matrix A and scale explicitely
														CGSMatrixCS &_mtra, CCorrectorC &_scla) {

	const char *funcname = "ExplicitSymmetricPointScaling";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Estimate the memory

	int ilist, iblk, isupbeg, isupend, nsupl, blai, blaj;
	int isuploc, isup, jsup, nzal, j;

	int nzscl = 0;
	int nzacolmax = 0;
	int blamax = 0;

	int *pldcorrarr = _scla.GetLdcorrarr ();
	int *pcorrsizes = _scla.GetCorrsizes ();
	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		nzscl = 0;

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			if (blai > blamax) blamax = blai;
//			nzscl += blai*blai;
			nzscl += blai;
			nzal = 0;
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				if (blaj > blamax) blamax = blaj;
				nzal += blai*blaj;
			};
			if (nzal > nzacolmax) nzacolmax = nzal;
		};

// Store the memory

		pldcorrarr[ilist] = nzscl;
		pcorrsizes[ilist] = 1;

	};

// Allocate the scaling data

	dcmplx *aloc, *qloc, *tauloc, *uloc, *vloc, *work;
	double *svloc, *rwork;

	int lwork = 10 * blamax;

	aloc = new dcmplx [nzacolmax];
	if (!aloc) MemoryFail (funcname);
	qloc = new dcmplx [nzacolmax];
	if (!qloc) MemoryFail (funcname);
	tauloc = new dcmplx [blamax];
	if (!tauloc) MemoryFail (funcname);
	uloc = new dcmplx [blamax*blamax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [blamax*blamax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [blamax];
	if (!svloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Init the hystogram

	CHyst hyst;

// Compute the scaling block column by block column

	int kii, kjj, ibs;
	double daux, dauxinv;
	dcmplx caux, caux1, caux2;

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

//		cout << " Myid = " << myid << " Scale Iblk = " << iblk << endl;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Read the block column from the disk if necessary

		_mtra.ReadBlock (iblk);

// Allocate the scaling 

		_scla.AllocateCorr (ilist);

// Cycle over the supernodes

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

// Compute scaling

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					ibs = _mtra.mtrarr[iblk].bsa[j];
					if (isup == jsup) {
						for (kii=0;kii<blai;kii++) {
							caux = _mtra.mtrarr[iblk].a[ibs+kii*blai+kii];
							daux = caux.x * caux.x + caux.y * caux.y;
							daux = sqrt(daux);
							svloc[kii] = daux;
							daux = sqrt(daux);
							dauxinv = 1.0e0 / daux;
							caux.x = daux;
							caux.y = 0.0e0;
							pyarr[ilist][nzscl+kii] = caux;
							caux.x = dauxinv;
							caux.y = 0.0e0;
							pwarr[ilist][nzscl+kii] = caux;
						};
					};
//				};
			};

			hyst.UpdateHystogram (blai,svloc);

			nzscl += blai;

		};

// Free the block column

		_mtra.FreeBlock (iblk);

// Write the scaling to the disk

		_scla.WriteCorr (ilist);
		_scla.FreeCorr (ilist);

	};

// Exchange the scaling among processors

	if (nproc != 1) throw " CGSMatrixCS::ExplicitSymmetricPointScaling: the case nproc > 1 is not implemented yet";

// Read scaling data into main memory

	pyarr = _scla.GetYarr ();
	pwarr = _scla.GetWarr ();

	dcmplx *sclarr;

	int nloc = _mtra.GetN ();

	sclarr = new dcmplx [nloc];
	if (!sclarr) MemoryFail (funcname);

	int niloc, ibeg, iishift, jjshift;

	for (iblk=0;iblk<_mtra.nblksc;iblk++) {

		_scla.ReadCorr (iblk);

		niloc = _mtra.bl2ndc[iblk+1]-_mtra.bl2ndc[iblk];
		ibeg = _mtra.bl2ndc[iblk];

		for (kii=0;kii<niloc;kii++) {
			sclarr[ibeg+kii] = pwarr[iblk][kii];
		};

		_scla.FreeCorr (iblk);

	};

// Scale explicitely

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

//		cout << " Myid = " << myid << " Scale Iblk = " << iblk << endl;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Read the block column from the disk if necessary

		_mtra.ReadBlock (iblk);

// Cycle over the supernodes

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;
			iishift = _mtra.sprndc[isup];

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

// Scale block column

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				jjshift = _mtra.sprndr[jsup];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				ibs = _mtra.mtrarr[iblk].bsa[j];
				for (kii=0;kii<blai;kii++) {
					for (kjj=0;kjj<blaj;kjj++) {
						caux1 = sclarr[iishift+kii] * sclarr[jjshift+kjj];
						caux = _mtra.mtrarr[iblk].a[ibs+kjj*blai+kii] * caux1;
						_mtra.mtrarr[iblk].a[ibs+kjj*blai+kii] = caux;
					};
				};
			};

		};

// Rewrite the block column to the disk if necessary

		_mtra.ReWriteBlock (iblk);
		_mtra.FreeBlock (iblk);

	};

// Output scaling hystogram

	CMPIComm pcomm = _tree.GetComm ();

	hyst.MakeGlobalHystogram (pcomm);

	_fout << " ExplSymmPointScaling: " << hyst;
	cout << " ExplSymmPointScaling: " << hyst;

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout  << " Point Scaling generation statistics: " << endl;
		_fout  << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Point Scaling generation statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work memory

	delete [] aloc;
	delete [] qloc;
	delete [] tauloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;
	delete [] sclarr;

};

// Author: Kharchenko S.A. 
// Description: Compute point scaling of the rectangular matrix A and scale explicitely
// CGSMatrixCS::ExplicitPointScaling()
//========================================================================================
void CGSMatrixCS::ExplicitPointScaling (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute point scaling of the rectangular matrix A and scale explicitely
														CGSMatrixCS &_mtra, int *_iamask, int *_jamask,
														CCorrectorC &_scla) {

	const char *funcname = "ExplicitPointScaling";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

	int myid = _tree.GetMyid ();

// Estimate the memory

	int ilist, iblk, isupbeg, isupend, nsupl, blai, blaj;
	int isuploc, isup, jsup, nzal, j;

	int nzscl = 0;
	int nzacolmax = 0;
	int blamax = 0;

	int *pldcorrarr = _scla.GetLdcorrarr ();
	int *pcorrsizes = _scla.GetCorrsizes ();
	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		nzscl = 0;

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			if (blai > blamax) blamax = blai;
//			nzscl += blai*blai;
			nzscl += blai;
			nzal = 0;
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				if (blaj > blamax) blamax = blaj;
				nzal += blai*blaj;
			};
			if (nzal > nzacolmax) nzacolmax = nzal;
		};

// Store the memory

		pldcorrarr[ilist] = nzscl;
		pcorrsizes[ilist] = 1;

	};

// Allocate the scaling data

	dcmplx *aloc, *qloc, *tauloc, *uloc, *vloc, *work;
	double *svloc, *rwork;

	int lwork = 10 * blamax;

	aloc = new dcmplx [nzacolmax];
	if (!aloc) MemoryFail (funcname);
	qloc = new dcmplx [nzacolmax];
	if (!qloc) MemoryFail (funcname);
	tauloc = new dcmplx [blamax];
	if (!tauloc) MemoryFail (funcname);
	uloc = new dcmplx [blamax*blamax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [blamax*blamax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [blamax];
	if (!svloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Init the hystogram

	CHyst hyst;

// Compute the scaling block column by block column and scale explicitely

	int kii, kjj, ni, nj, njloc, ibs, ishift, jold;
	double daux, dauxloc, dauxinv;
	dcmplx caux, caux1, caux2;

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

//		cout << " Myid = " << myid << " Scale Iblk = " << iblk << endl;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		ishift = _iamask[isupbeg];

		nsupl = isupend-isupbeg+1;

// Read the block column from the disk if necessary

		_mtra.ReadBlock (iblk);

// Allocate the scaling 

		_scla.AllocateCorr (ilist);

// Cycle over the supernodes

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

// Count the number of elements in the block row

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

			ni = blai;
			nj = 0;

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					nj += blaj;
//				};
			};

// Take block column

			njloc = 0;

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					ibs = _mtra.mtrarr[iblk].bsa[j];
					for (kii=0;kii<blai;kii++) {
						for (kjj=0;kjj<blaj;kjj++) {
							aloc[kii*nj+njloc+kjj] = _mtra.mtrarr[iblk].a[ibs+kjj*ni+kii];
						};
					};
					njloc += blaj;
//				};
			};

// Compute scaling and scale

			for (kii=0;kii<ni;kii++) {
				daux = 0.0e0;
				for (kjj=0;kjj<nj;kjj++) {
					dauxloc = aloc[kii*nj+kjj].x;
					daux += dauxloc*dauxloc;
					dauxloc = aloc[kii*nj+kjj].y;
					daux += dauxloc*dauxloc;
				};
				daux = sqrt(daux);
				svloc[kii] = daux;
				dauxinv = 1.0e0 / daux;
				caux.x = daux;
				caux.y = 0.0e0;
				pyarr[ilist][nzscl+kii] = caux;
				caux.x = dauxinv;
				caux.y = 0.0e0;
				pwarr[ilist][nzscl+kii] = caux;
				for (kjj=0;kjj<nj;kjj++) {
					caux = aloc[kii*nj+kjj] * dauxinv;
					aloc[kii*nj+kjj] = caux;
				};
			};

			hyst.UpdateHystogram (ni,svloc);

			njloc = 0;

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					ibs = _mtra.mtrarr[iblk].bsa[j];
					for (kii=0;kii<blai;kii++) {
						for (kjj=0;kjj<blaj;kjj++) {
							_mtra.mtrarr[iblk].a[ibs+kjj*ni+kii] = aloc[kii*nj+njloc+kjj];
						};
					};
					njloc += blaj;
//				} else {
//					jsup = _mtra.mtrarr[iblk].ja[j];
//					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
//					ibs = _mtra.mtrarr[iblk].bsa[j];
//					for (kii=0;kii<blai;kii++) {
//						for (kjj=0;kjj<blaj;kjj++) {
//							caux = _mtra.mtrarr[iblk].a[ibs+kjj*ni+kii] * pwarr[ilist][nzscl+kii];
//							_mtra.mtrarr[iblk].a[ibs+kjj*ni+kii] = caux;
//						};
//					};
//				};
			};

			nzscl += blai;

		};

// Rewrite the block column to the disk if necessary

		_mtra.ReWriteBlock (iblk);
		_mtra.FreeBlock (iblk);

// Write the scaling to the disk

		_scla.WriteCorr (ilist);
		_scla.FreeCorr (ilist);

	};

// Output scaling hystogram

	CMPIComm pcomm = _tree.GetComm ();

	hyst.MakeGlobalHystogram (pcomm);

	_fout << " ExplicitPointScaling: " << hyst;
	cout << " ExplicitPointScaling: " << hyst;

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout  << " Point Scaling generation statistics: " << endl;
		_fout  << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Point Scaling generation statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work memory

	delete [] aloc;
	delete [] qloc;
	delete [] tauloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;

};

// Author: Kharchenko S.A. 
// Description: Compute block scaling of the rectangular matrix A and scale explicitely
// CGSMatrixCS::ExplicitBlockScaling()
//========================================================================================
void CGSMatrixCS::ExplicitBlockScaling (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute block scaling of the rectangular matrix A and scale explicitely
														CGSMatrixCS &_mtra, int *_iamask, int *_jamask,
														CCorrectorC &_scla) {

	const char *funcname = "ExplicitBlockScaling";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

	int myid = _tree.GetMyid ();

// Estimate the memory

	int ilist, iblk, isupbeg, isupend, nsupl, blai, blaj;
	int isuploc, isup, jsup, nzal, j;

	int nzscl = 0;
	int nzacolmax = 0;
	int blamax = 0;

	int *pldcorrarr = _scla.GetLdcorrarr ();
	int *pcorrsizes = _scla.GetCorrsizes ();
	dcmplx **pyarr = _scla.GetYarr ();
	dcmplx **pwarr = _scla.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		nzscl = 0;

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			if (blai > blamax) blamax = blai;
			nzscl += blai*blai;
			nzal = 0;
			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jsup = _mtra.mtrarr[iblk].ja[j];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				if (blaj > blamax) blamax = blaj;
				nzal += blai*blaj;
			};
			if (nzal > nzacolmax) nzacolmax = nzal;
		};

// Store the memory

		pldcorrarr[ilist] = nzscl;
		pcorrsizes[ilist] = 1;

	};

// Allocate the scaling data

	dcmplx *aloc, *qloc, *tauloc, *uloc, *vloc, *work;
	double *svloc, *rwork;

	int lwork = 10 * blamax;

	aloc = new dcmplx [nzacolmax];
	if (!aloc) MemoryFail (funcname);
	qloc = new dcmplx [nzacolmax];
	if (!qloc) MemoryFail (funcname);
	tauloc = new dcmplx [blamax];
	if (!tauloc) MemoryFail (funcname);
	uloc = new dcmplx [blamax*blamax];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [blamax*blamax];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [blamax];
	if (!svloc) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

// Init the hystogram

	CHyst hyst;

// Compute the scaling block column by block column and scale explicitely

	bool bmult = false;

	int kii, kjj, ni, nj, njloc, ibs, ishift, jold;
	int info, kkk;
	double aux;
	dcmplx caux, caux1, caux2;

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		cout << " Myid = " << myid << " Scale Iblk = " << iblk << endl;

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		ishift = _iamask[isupbeg];

		nsupl = isupend-isupbeg+1;

// Read the block column from the disk if necessary

		_mtra.ReadBlock (iblk);

// Allocate the scaling 

		_scla.AllocateCorr (ilist);

// Cycle over the supernodes

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

// Count the number of elements in the block row

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

			ni = blai;
			nj = 0;

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					nj += blaj;
//				};
			};

// Take block column

			njloc = 0;

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
//				if (_jamask[jold] == 1) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					ibs = _mtra.mtrarr[iblk].bsa[j];
					for (kii=0;kii<blai;kii++) {
						for (kjj=0;kjj<blaj;kjj++) {
							aloc[kii*nj+njloc+kjj] = _mtra.mtrarr[iblk].a[ibs+kjj*ni+kii];
						};
					};
					njloc += blaj;
//				};
			};

// Compute QR decomposition

			QrdBlk (nj, ni, aloc, nj, tauloc);

// Take R matrix

			GetRPartQrd (nj, 0, ni,
								aloc, nj, uloc, ni);

			for (j=0;j<blai*blai;j++) pyarr[ilist][nzscl+j] = uloc[j];

// Compute Q explicitely

			if (!bmult) {

				for (j=0;j<ni*ni;j++) uloc[j] = czero;

				for (j=0;j<ni;j++) uloc[j*ni+j] = cone;

				MvmQBlk (nj, ni, ni,
							aloc, nj, tauloc, 
							uloc, ni, qloc, nj);

// Store Q matrix

				njloc = 0;

				for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
					jold = j + ishift;
//					if (_jamask[jold] == 1) {
						jsup = _mtra.mtrarr[iblk].ja[j];
						blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
						ibs = _mtra.mtrarr[iblk].bsa[j];
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								_mtra.mtrarr[iblk].a[ibs+kjj*ni+kii] = qloc[kii*nj+njloc+kjj];
							};
						};
						njloc += blaj;
//					};
				};

			};

// Next block

			nzscl += blai*blai;

		};

// Cycle over the supernodes

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {
			isup = isuploc+isupbeg;
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

// Init the matrix

			for (kii=0;kii<blai*blai;kii++) aloc[kii] = pyarr[ilist][nzscl+kii];

// Compute SVD

			zgesvd_ ("A", "A", &blai, &blai,
						aloc, &blai, svloc, 
						uloc, &blai, vloc, &blai, 
						work, &lwork, rwork, &info);

			if (info != 0) throw " Error in the Lapack routine ZGESVD";

// Update hystogram

			hyst.UpdateHystogram (blai,svloc);

// Compute inverse matrix

			for (kii=0;kii<blai;kii++) {
				aux = 1.0e0 / svloc[kii];
				for (kjj=0;kjj<blai;kjj++) {
					caux = uloc[kii*blai+kjj];
					uloc[kii*blai+kjj] = caux * aux;
				};
			};

			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blai;kjj++) {
					caux = czero;
					for (kkk=0;kkk<blai;kkk++) {
						caux1 = vloc[kii*blai+kkk] * uloc[kkk*blai+kjj];
						caux2 = caux + conj(caux1);
						caux = caux2;
					};
					aloc[kjj*blai+kii] = caux;
				};
			};

// Store inverse matrix

			for (kii=0;kii<blai*blai;kii++) pwarr[ilist][nzscl+kii] = aloc[kii];

			nzscl += blai*blai;

		};

// Scale remaining (not scaled yet) blocks of the block column
/*
		dcmplx *pa;
		dcmplx *pw;

		int iind, jind;

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

			for (j=_mtra.mtrarr[iblk].ia[isuploc];j<_mtra.mtrarr[iblk].ia[isuploc+1];j++) {
				jold = j + ishift;
				if (_jamask[jold] == 0 || bmult) {
					jsup = _mtra.mtrarr[iblk].ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					ibs = _mtra.mtrarr[iblk].bsa[j];
					pa = _mtra.mtrarr[iblk].a+ibs;
					pw = pwarr[ilist]+nzscl;
					for (kii=0;kii<blai;kii++) {
						for (kjj=0;kjj<blaj;kjj++) {
							caux = czero;
							iind = kii*blaj;
							jind = kjj*blai;
							for (kkk=0;kkk<blai;kkk++) {
//								caux1 = _mtra.mtrarr[iblk].a[ibs+kii*blaj+kkk] * pwarr[ilist][nzscl+kjj*blai+kkk];
								caux1 = pa[iind++] * pw[jind++];
								caux += caux1;
							};
							aloc[kii*blaj+kjj] = caux;
						};
					};
					for (kii=0;kii<blai*blaj;kii++) _mtra.mtrarr[iblk].a[ibs+kii] = aloc[kii];
				};
			};

			nzscl += blai*blai;

		};
*/
// Rewrite the block column to the disk if necessary

		_mtra.ReWriteBlock (iblk);
		_mtra.FreeBlock (iblk);

// Write the scaling to the disk

		_scla.WriteCorr (ilist);
		_scla.FreeCorr (ilist);

	};

// Output scaling hystogram

	CMPIComm pcomm = _tree.GetComm ();

	hyst.MakeGlobalHystogram (pcomm);

	_fout << " ExplicitScaling: " << hyst;
	cout << " ExplicitScaling: " << hyst;

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout  << " Scaling generation statistics: " << endl;
		_fout  << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " Scaling generation statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work memory

	delete [] aloc;
	delete [] qloc;
	delete [] tauloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] work;
	delete [] rwork;

};

// Author: Kharchenko S.A. 
// Description: Perform explicit point scaling of the solution
// CGSMatrixCS::VectorPointScale()
//========================================================================================
void CGSMatrixCS::VectorPointScale (int _ivar, // Perform explicit point scaling of the solution
													CGSMatrixCS &_mtra,
													CCorrectorC &_scl, CSVectorC &_sol) {

	const char *funcname = "VectorPointScale";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ilist, iblk, isupbeg, isupend, nsupl;
	int isuploc, isup, blai, kii;

	int blamax = 0;

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		if (_mtra.mtrarr[iblk].blamx > blamax) blamax = _mtra.mtrarr[iblk].blamx;

	};

	dcmplx *vectloc;

	vectloc = new dcmplx [blamax];
	if (!vectloc) MemoryFail (funcname);

// Main cycle over the supernodes

	if (_sol.nrhs != 1) {
		throw " CGSMatrixCS::VectorBlockScale: nrhs > 1 is not implemented yet ";
	};

	dcmplx caux, caux1;

	int nzscl;
	int nloc = 0;

	dcmplx **pyarr = _scl.GetYarr ();
	dcmplx **pwarr = _scl.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

		_scl.ReadCorr (ilist);

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

			if (_ivar == 0) {

				for (kii=0;kii<blai;kii++) {
					caux = pyarr[ilist][nzscl+kii] * _sol.vect[nloc+kii];
					vectloc[kii] = caux;
				};

			} else {

				for (kii=0;kii<blai;kii++) {
					caux = pwarr[ilist][nzscl+kii] * _sol.vect[nloc+kii];
					vectloc[kii] = caux;
				};

			};

			for (kii=0;kii<blai;kii++) _sol.vect[nloc+kii] = vectloc[kii];

			nzscl += blai;
			nloc += blai;

		};

		_scl.FreeCorr (ilist);

	};

// Free work arrays

	delete [] vectloc;

};

// Author: Kharchenko S.A. 
// Description: Perform explicit block scaling of the solution
// CGSMatrixCS::VectorBlockScale()
//========================================================================================
void CGSMatrixCS::VectorBlockScale (int _ivar,  // Perform explicit block scaling of the solution
													CGSMatrixCS &_mtra,
													CCorrectorC &_scl, CSVectorC &_sol) {

	const char *funcname = "VectorBlockScale";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ilist, iblk, isupbeg, isupend, nsupl;
	int isuploc, isup, blai, kii, kjj;

	int blamax = 0;

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		if (_mtra.mtrarr[iblk].blamx > blamax) blamax = _mtra.mtrarr[iblk].blamx;

	};

	dcmplx *vectloc;

	vectloc = new dcmplx [blamax];
	if (!vectloc) MemoryFail (funcname);

// Main cycle over the supernodes

	if (_sol.nrhs != 1) {
		throw " CGSMatrixCS::VectorBlockScale: nrhs > 1 is not implemented yet ";
	};

	dcmplx caux, caux1;

	int nzscl;
	int nloc = 0;

	dcmplx **pyarr = _scl.GetYarr ();
	dcmplx **pwarr = _scl.GetWarr ();

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

		_scl.ReadCorr (ilist);

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];

			if (_ivar == 0) {

				for (kii=0;kii<blai;kii++) {
					caux = czero;
					for (kjj=0;kjj<blai;kjj++) {
						caux1 = caux + pyarr[ilist][nzscl+kjj*blai+kii] * _sol.vect[nloc+kjj];
						caux = caux1;
					};
					vectloc[kii] = caux;
				};

			} else {

				for (kii=0;kii<blai;kii++) {
					caux = czero;
					for (kjj=0;kjj<blai;kjj++) {
						caux1 = caux + pwarr[ilist][nzscl+kjj*blai+kii] * _sol.vect[nloc+kjj];
						caux = caux1;
					};
					vectloc[kii] = caux;
				};

			};

			for (kii=0;kii<blai;kii++) _sol.vect[nloc+kii] = vectloc[kii];

			nzscl += blai * blai;
			nloc += blai;

		};

		_scl.FreeCorr (ilist);

	};

// Free work arrays

	delete [] vectloc;

};

// Author: Kharchenko S.A. 
// Description: Perform explicit block scaling of the solution with extension
// CGSMatrixCS::ExtendedVectorBlockScale()
//========================================================================================
void CGSMatrixCS::ExtendedVectorBlockScale (int _ivar,  // Perform explicit block scaling of the solution with extension
													CGSMatrixCS &_mtraflt,
													int *_sprndsini, CCorrectorC &_scl, 
													CSVectorC &_sol, CSVectorC &_solext) {

	const char *funcname = "ExtendedVectorBlockScale";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ilist, iblk, isupbeg, isupend, nsupl;
	int isuploc, isup, blaiflt, blaiini, kii, kjj;

	int blamax = 0;

	for (ilist=0;ilist<_mtraflt.nlist;ilist++) {

		iblk = _mtraflt.listb[ilist];

		if (_mtraflt.mtrarr[iblk].blamx > blamax) blamax = _mtraflt.mtrarr[iblk].blamx;

	};

	dcmplx *vectloc;

	vectloc = new dcmplx [blamax];
	if (!vectloc) MemoryFail (funcname);

// Main cycle over the supernodes

	if (_sol.nrhs != 1) {
		throw " CGSMatrixCS::ExtendedVectorBlockScale: nrhs > 1 is not implemented yet ";
	};

	dcmplx caux, caux1;

	int nzscl;
	int nloc = 0;
	int nlocext = 0;

	dcmplx **pyarr = _scl.GetYarr ();
	dcmplx **pwarr = _scl.GetWarr ();

	for (ilist=0;ilist<_mtraflt.nlist;ilist++) {

		iblk = _mtraflt.listb[ilist];

		isupbeg = _mtraflt.blksc[iblk];
		isupend = _mtraflt.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

		_scl.ReadCorr (ilist);

		nzscl = 0;

		for (isuploc=0;isuploc<nsupl;isuploc++) {

			isup = isuploc+isupbeg;

			blaiflt = _mtraflt.sprndc[isup+1]-_mtraflt.sprndc[isup];
			blaiini = _sprndsini[isup+1]-_sprndsini[isup];

			if (_ivar == 0) {

				for (kii=0;kii<blaiini;kii++) {
					caux = czero;
					for (kjj=0;kjj<blaiflt;kjj++) {
						caux1 = caux + pyarr[ilist][nzscl+kjj*blaiini+kii] * _sol.vect[nloc+kjj];
						caux = caux1;
					};
					vectloc[kii] = caux;
				};

			} else {

				for (kii=0;kii<blaiini;kii++) {
					caux = czero;
					for (kjj=0;kjj<blaiflt;kjj++) {
						caux1 = caux + pwarr[ilist][nzscl+kjj*blaiini+kii] * _sol.vect[nloc+kjj];
						caux = caux1;
					};
					vectloc[kii] = caux;
				};

			};

			for (kii=0;kii<blaiini;kii++) _solext.vect[nlocext+kii] = vectloc[kii];

			nzscl += blaiini * blaiini;
			nloc += blaiflt;
			nlocext += blaiini;

		};

		_scl.FreeCorr (ilist);

	};

// Free work arrays

	delete [] vectloc;

};

// Author: Kharchenko S.A. 
// Description: Compute complex A^hA matrix for rectangular A in serial mode
// CGSMatrixCS::AtAMatr()
//========================================================================================
void CGSMatrixCS::AtAMatr (ofstream &_fout, const CSlvParam &_param, // Compute complex A^hA matrix for rectangular A in serial mode
									CGSMatrixCS &_mtra,
									CSMatrix &_mtratasp,
									CGSMatrixCS &_mtratau) { 

	const char *funcname = "AtAMatr";

	dcmplx czero (0.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtra.nsupc];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;
	for (i=0;i<_mtra.nblksc;i++) {
		for (j=_mtra.blksc[i];j<_mtra.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Create the block sparsity

	CSMatrix ablkloc;

	ablkloc = _mtratasp.BlockSparsity (_mtra.nblksc,sp2blk);

// Init atau matrix

	_mtratau = _mtra;

// Modify block rows data

	_mtratau.m = _mtratau.n;
	for (i=0;i<=_mtratau.nblksc;i++) {
		_mtratau.blksr[i] = _mtratau.blksc[i];
		_mtratau.bl2ndr[i] = _mtratau.bl2ndc[i];
	};
	_mtratau.nblksr = _mtratau.nblksc;
	for (i=0;i<=_mtratau.nsupc;i++) {
		_mtratau.sprndr[i] = _mtratau.sprndc[i];
	};
	_mtratau.nsupr = _mtratau.nsupc;

// Init out of core data if necessary

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"AtAMatr_",0,"_");
		_mtratau.SetupFiles (_param.ooctyp, strbuff);

	};

// Store sparsity of the block rows

	int iblk, isupbeg, isupend, nsupl, nzjal, nzal, isup, blai, jsup;
	int blaj, ni;

	CSMatrixCS mtrdummy;

	int blamxloc = 0;

	for (i=0;i<_mtratau.nsupc;i++) {
		ni = _mtratau.sprndc[i+1]-_mtratau.sprndc[i];
		if (ni > blamxloc) blamxloc = ni;
	};

	for (iblk=0;iblk<_mtratau.nblksc;iblk++) {

		isupbeg = _mtratau.blksc[iblk];
		isupend = _mtratau.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		nzjal = 0;
		nzal = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtratau.sprndc[isup+1]-_mtratau.sprndc[isup];
			for (j=_mtratasp.ia[isup];j<_mtratasp.ia[isup+1];j++) {
				jsup = _mtratasp.ja[j];
				if (jsup >= isup) {
					blaj = _mtratau.sprndc[jsup+1]-_mtratau.sprndc[jsup];
					nzjal++;
					nzal += blai*blaj;
				};
			};
		};

// Allocate the data

		CSMatrixCS temp (nsupl, nzjal, nzal);

// Init the matrix

		temp.m = _mtratau.m;
		temp.n = _mtratau.n;
		temp.nsupr = nsupl;
		temp.nsupc = nsupl;
		temp.nlist = nsupl;
		temp.nzja = nzjal;
		temp.nza = nzal;
		temp.nzatot = nzal;
		temp.blamx = blamxloc;

		temp.sprndr[0] = 0;
		temp.sprndc[0] = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			temp.list[isup-isupbeg] = isup;
			temp.sprndr[isup-isupbeg+1] = 0;
			temp.sprndc[isup-isupbeg+1] = 0;
		};

		temp.ia[0] = 0;

		nzjal = 0;
		nzal = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtratau.sprndc[isup+1]-_mtratau.sprndc[isup];
			for (j=_mtratasp.ia[isup];j<_mtratasp.ia[isup+1];j++) {
				jsup = _mtratasp.ja[j];
				if (jsup >= isup) {
					blaj = _mtratau.sprndc[jsup+1]-_mtratau.sprndc[jsup];
					temp.ja[nzjal] = jsup;
					temp.bsa[nzjal] = nzal;
					nzjal++;
					nzal += blai*blaj;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

		};

// Store the matrix

		_mtratau.mtrarr[iblk] = temp;

		_mtratau.FreeBlock (iblk);

		temp = mtrdummy;

	};

// Compute and store block rows of atau

	int jblk;

	for (iblk=0;iblk<_mtratau.nblksc;iblk++) {

// Allocate data

		_mtratau.AllocateBlock (iblk);

// Init data by zeroes

		for (i=0;i<_mtratau.mtrarr[iblk].nzatot;i++) _mtratau.mtrarr[iblk].a[i] = czero;

// Read current block of A

		_mtra.ReadBlock (iblk);

// Cycle over the other blocks of A

		for (j=ablkloc.ia[iblk];j<ablkloc.ia[iblk+1];j++) {

			jblk = ablkloc.ja[j];

// Read current block if necessary

			if (jblk != iblk) {
				_mtra.ReadBlock (jblk);
			};

// Update the elements

			_mtratau.UpdateAtAMatrixBlock (iblk, jblk, _mtra, sp2blk);

// Free current block if necessary

			if (jblk != iblk) {
				_mtra.FreeBlock (jblk);
			};

		};

// Free current block of A

		_mtra.FreeBlock (iblk);

// Write current block of AtA

		_mtratau.WriteBlock (iblk);
		_mtratau.FreeBlock (iblk);

	};

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout  << " AtA generation statistics: " << endl;
		_fout  << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " AtA generation statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work arrays

	delete [] sp2blk;

};

// Author: Kharchenko S.A. 
// Description: Compute complex A^hA matrix for rectangular A in parallel mode
// CGSMatrixCS::AtAMatr()
//========================================================================================
void CGSMatrixCS::AtAMatr (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute complex A^hA matrix for rectangular A in parallel mode
									CSMatrix &_mtrasp,
									CGSMatrixCS &_mtra,
									CSMatrix &_mtratasp,
									CGSMatrixCS &_mtratau) { 

	const char *funcname = "AtAMatr";

	dcmplx czero (0.0e0,0.0e0);

// Init time measurement

	clock_t time0, time1;
	double tottim;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtra.nsupc];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;
	for (i=0;i<_mtra.nblksc;i++) {
		for (j=_mtra.blksc[i];j<_mtra.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Open classes

	int nnodestree = _tree.GetNnodes ();
	CNode *pnode = _tree.GetNode ();

	int nproc = _tree.GetNproc ();
	int myid = _tree.GetMyid ();

// Create the block sparsity

	CSMatrix ablkloc;

	ablkloc = _mtratasp.BlockSparsity (_mtra.nblksc,sp2blk);

// Init atau matrix

	_mtratau = _mtra;

// Modify block rows data

	_mtratau.m = _mtratau.n;
	for (i=0;i<=_mtratau.nblksc;i++) {
		_mtratau.blksr[i] = _mtratau.blksc[i];
		_mtratau.bl2ndr[i] = _mtratau.bl2ndc[i];
	};
	_mtratau.nblksr = _mtratau.nblksc;
	for (i=0;i<=_mtratau.nsupc;i++) {
		_mtratau.sprndr[i] = _mtratau.sprndc[i];
	};
	_mtratau.nsupr = _mtratau.nsupc;

// Init out of core data if necessary

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"AtAMatr_",myid,"_");
		_mtratau.SetupFiles (_param.ooctyp, strbuff);

	};

// Store sparsity of the block rows

	int ilist, iblk, isupbeg, isupend, nsupl, nzjal, nzal, isup, blai, jsup;
	int blaj, ni;

	CSMatrixCS mtrdummy;

	int blamxloc = 0;

	for (i=0;i<_mtratau.nsupc;i++) {
		ni = _mtratau.sprndc[i+1]-_mtratau.sprndc[i];
		if (ni > blamxloc) blamxloc = ni;
	};

	for (ilist=0;ilist<_mtratau.nlist;ilist++) {

		iblk = _mtratau.listb[ilist];

		isupbeg = _mtratau.blksc[iblk];
		isupend = _mtratau.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		nzjal = 0;
		nzal = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtratau.sprndc[isup+1]-_mtratau.sprndc[isup];
			for (j=_mtratasp.ia[isup];j<_mtratasp.ia[isup+1];j++) {
				jsup = _mtratasp.ja[j];
				if (jsup >= isup) {
					blaj = _mtratau.sprndc[jsup+1]-_mtratau.sprndc[jsup];
					nzjal++;
					nzal += blai*blaj;
				};
			};
		};

// Allocate the data

		CSMatrixCS temp (nsupl, nzjal, nzal);

// Init the matrix

		temp.m = _mtratau.m;
		temp.n = _mtratau.n;
		temp.nsupr = nsupl;
		temp.nsupc = nsupl;
		temp.nlist = nsupl;
		temp.nzja = nzjal;
		temp.nza = nzal;
		temp.nzatot = nzal;
		temp.blamx = blamxloc;

		temp.sprndr[0] = 0;
		temp.sprndc[0] = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			temp.list[isup-isupbeg] = isup;
			temp.sprndr[isup-isupbeg+1] = 0;
			temp.sprndc[isup-isupbeg+1] = 0;
		};

		temp.ia[0] = 0;

		nzjal = 0;
		nzal = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtratau.sprndc[isup+1]-_mtratau.sprndc[isup];
			for (j=_mtratasp.ia[isup];j<_mtratasp.ia[isup+1];j++) {
				jsup = _mtratasp.ja[j];
				if (jsup >= isup) {
					blaj = _mtratau.sprndc[jsup+1]-_mtratau.sprndc[jsup];
					temp.ja[nzjal] = jsup;
					temp.bsa[nzjal] = nzal;
					nzjal++;
					nzal += blai*blaj;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

		};

// Store the matrix

		_mtratau.mtrarr[iblk] = temp;

		_mtratau.FreeBlock (iblk);

		temp = mtrdummy;

	};

// Compute distribution of blocks over cpu's

	int *blk2cpu;

	blk2cpu = new int [_mtratau.nblksc];
	if (!blk2cpu) MemoryFail (funcname);

	for (i=0;i<_mtratau.nblksc;i++) blk2cpu[i] = -1;

	int inode, iproc, indbeg, indend;

	for (inode=0;inode<nnodestree;inode++) {
		iproc = pnode[inode].GetNodecpu ();
		indbeg = pnode[inode].GetIndbeg ();
		indend = pnode[inode].GetIndend ();
		for (j=indbeg;j<=indend;j++) {
			blk2cpu[j] = iproc;
		};
	};

// For each block row determine which its blocks are necessary on the other processors

	int *imaskja;

	int nzjaloc = _mtrasp.GetNzja ();

	imaskja = new int [nzjaloc];
	if (!imaskja) MemoryFail (funcname);

	for (i=0;i<nzjaloc;i++) imaskja[i] = -1;

	int jblk, jproc;

	for (iblk=0;iblk<_mtratau.nblksc;iblk++) {

		iproc = blk2cpu[iblk];

// Cycle over the other blocks of A

		if (iproc != myid) {

			for (j=ablkloc.ia[iblk];j<ablkloc.ia[iblk+1];j++) {

				jblk = ablkloc.ja[j];

				jproc = blk2cpu[jblk];

	// Update the structural elements if necessary

				if (jproc == myid) {

					_mtrasp.UpdateAtAMatrixStruct (_mtratasp, 
																iblk, jblk, _mtratau.blksc, sp2blk, imaskja);

				};

			};

		};

	};

// Create local structures to be send

	CGSMatrixCS mtrasndrcv;

	mtrasndrcv = _mtra;

	blamxloc = 0;

	for (i=0;i<_mtra.nsupc;i++) {
		ni = _mtra.sprndc[i+1]-_mtra.sprndc[i];
		if (ni > blamxloc) blamxloc = ni;
	};
	for (i=0;i<_mtra.nsupr;i++) {
		ni = _mtra.sprndr[i+1]-_mtra.sprndr[i];
		if (ni > blamxloc) blamxloc = ni;
	};

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		isupbeg = _mtra.blksc[iblk];
		isupend = _mtra.blksc[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

		nzjal = 0;
		nzal = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			for (j=_mtrasp.ia[isup];j<_mtrasp.ia[isup+1];j++) {
				jsup = _mtrasp.ja[j];
				if (imaskja[j] >= 0) {
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					nzjal++;
					nzal += blai*blaj;
				};
			};
		};

// Allocate the data

		CSMatrixCS temp (nsupl, nzjal, nzal);

// Read current block of A

		_mtra.ReadBlock (iblk);

// Init the matrix

		temp.m = _mtra.m;
		temp.n = _mtra.n;
		temp.nsupr = nsupl;
		temp.nsupc = nsupl;
		temp.nlist = nsupl;
		temp.nzja = nzjal;
		temp.nza = nzal;
		temp.nzatot = nzal;
		temp.blamx = blamxloc;

		temp.sprndr[0] = 0;
		temp.sprndc[0] = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			temp.list[isup-isupbeg] = isup;
			temp.sprndr[isup-isupbeg+1] = 0;
			temp.sprndc[isup-isupbeg+1] = 0;
		};

		temp.ia[0] = 0;

		nzjal = 0;
		nzal = 0;
		int nzaini = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {
			blai = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			for (j=_mtrasp.ia[isup];j<_mtrasp.ia[isup+1];j++) {
				jsup = _mtrasp.ja[j];
				blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
				if (imaskja[j] >= 0) {
					temp.ja[nzjal] = jsup;
					temp.bsa[nzjal] = nzal;
					nzjal++;
					for (i=0;i<blai*blaj;i++) temp.a[nzal+i] = _mtra.mtrarr[iblk].a[nzaini+i];
					nzal += blai*blaj;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
				nzaini += blai*blaj;
			};

		};

// Free current block of A

		_mtra.FreeBlock (iblk);

// Store the matrix

		mtrasndrcv.mtrarr[iblk] = temp;

		temp = mtrdummy;

	};

// Prepare send and receive lists

	int *imaskblk;
	int *listblk;
	int *iablk;
	int *jablk;

	imaskblk = new int [_mtra.nblksc];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [_mtra.nblksc];
	if (!listblk) MemoryFail (funcname);
	iablk = new int [nproc+1];
	if (!iablk) MemoryFail (funcname);
	jablk = new int [nproc*_mtra.nblksc];
	if (!jablk) MemoryFail (funcname);

	for (i=0;i<_mtra.nblksc;i++) imaskblk[i] = -1;

	iablk[0] = 0;

	int icycle = -1;

	int iiproc, nlistloc;

	int nz = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		icycle++;
		nlistloc = 0;
		if (iproc != myid) {
			for (iblk=0;iblk<_mtratau.nblksc;iblk++) {
				iiproc = blk2cpu[iblk];

// Cycle over the other blocks of A

				if (iiproc == iproc) {
					for (j=ablkloc.ia[iblk];j<ablkloc.ia[iblk+1];j++) {
						jblk = ablkloc.ja[j];
						jproc = blk2cpu[jblk];
						if (jproc == myid) {
							if (imaskblk[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imaskblk[jblk] = icycle;
							};
						};
					};
				};

			};

// Sort blocks and store

			if (nlistloc != 0) qsort (listblk, nlistloc, sizeof(int), compint);

			for (j=0;j<nlistloc;j++) jablk[nz+j] = listblk[j];

			nz += nlistloc;

		};
		iablk[iproc+1] = nz;
	};

// Prepare send data for each cpu (combine submatrices)

	CSMatrixCS *sndarr;
	CSMatrixCS *rcvarr;

	sndarr = new CSMatrixCS [nproc];
	if (!sndarr) MemoryFail (funcname);
	rcvarr = new CSMatrixCS [nproc];
	if (!rcvarr) MemoryFail (funcname);

	int ibeg;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			nlistloc = iablk[iproc+1]-iablk[iproc];
			ibeg = iablk[iproc];
			if (nlistloc > 0) {
				_mtra.CombineMatrices (nlistloc, jablk+ibeg,
												mtrasndrcv.mtrarr,
												sndarr[iproc]);
			};
		};
	};

// Free prepared initial blocks

	for (ilist=0;ilist<_mtra.nlist;ilist++) {

		iblk = _mtra.listb[ilist];

		mtrasndrcv.mtrarr[iblk] = mtrdummy;

	};

// Prepare receive lists

	iablk[0] = 0;

	nz = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		icycle++;
		nlistloc = 0;
		if (iproc != myid) {
			for (iblk=0;iblk<_mtratau.nblksc;iblk++) {
				iiproc = blk2cpu[iblk];

// Cycle over the other blocks of A

				if (iiproc == myid) {
					for (j=ablkloc.ia[iblk];j<ablkloc.ia[iblk+1];j++) {
						jblk = ablkloc.ja[j];
						jproc = blk2cpu[jblk];
						if (jproc == iproc) {
							if (imaskblk[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imaskblk[jblk] = icycle;
							};
						};
					};
				};

			};

// Sort blocks and store

			if (nlistloc != 0) qsort (listblk, nlistloc, sizeof(int), compint);

			for (j=0;j<nlistloc;j++) jablk[nz+j] = listblk[j];

			nz += nlistloc;

		};
		iablk[iproc+1] = nz;
	};

// Perform sends and receives

	for (iproc=0;iproc<myid;iproc++) {
		mtrasndrcv.ReceiveMatrix (_tree.comm, iproc, myid, rcvarr[iproc]);
	};
	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			mtrasndrcv.SendMatrix (_tree.comm, iproc, iproc, sndarr[iproc]);
		};
	};
	for (iproc=myid+1;iproc<nproc;iproc++) {
		mtrasndrcv.ReceiveMatrix (_tree.comm, iproc, myid, rcvarr[iproc]);
	};

// Unpack received data into mtra arrays

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			nlistloc = iablk[iproc+1]-iablk[iproc];
			ibeg = iablk[iproc];
			if (nlistloc > 0) {
				_mtra.SplitMatrix (nlistloc, jablk+ibeg,
											rcvarr[iproc],
											_mtra.mtrarr);
			};
			for (j=0;j<nlistloc;j++) {
				iblk = jablk[ibeg+j];
				_mtra.WriteBlock (iblk);
				_mtra.FreeBlock (iblk);
			};
		};
	};

// Compute and store block rows of atau

	for (iblk=0;iblk<_mtratau.nblksc;iblk++) {

		iproc = blk2cpu[iblk];

		if (iproc == myid) {

// Allocate data

			_mtratau.AllocateBlock (iblk);

// Init data by zeroes

			for (i=0;i<_mtratau.mtrarr[iblk].nzatot;i++) _mtratau.mtrarr[iblk].a[i] = czero;

// Read current block of A

			_mtra.ReadBlock (iblk);

// Cycle over the other blocks of A

			for (j=ablkloc.ia[iblk];j<ablkloc.ia[iblk+1];j++) {

				jblk = ablkloc.ja[j];

// Read current block if necessary

				if (jblk != iblk) {
					_mtra.ReadBlock (jblk);
				};

// Update the elements

				_mtratau.UpdateAtAMatrixBlock (iblk, jblk, _mtra, sp2blk);

// Free current block if necessary

				if (jblk != iblk) {
					_mtra.FreeBlock (jblk);
				};

			};

// Free current block of A

			_mtra.FreeBlock (iblk);

// Write current block of AtA

			_mtratau.WriteBlock (iblk);
			_mtratau.FreeBlock (iblk);

		};

	};

// Free unpacked received data

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			for (j=iablk[iproc];j<iablk[iproc+1];j++) {
				iblk = jablk[j];
				_mtra.mtrarr[iblk] = mtrdummy;
			};
		};
	};

// Output statistics

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	if (_param.msglev >= 1) {
		_fout  << " AtA generation statistics: " << endl;
		_fout  << "     Time  = " << tottim << " sec. " << endl;
	};
	if (_param.msglev > 1) {
		cout  << " AtA generation statistics: " << endl;
		cout  << "     Time  = " << tottim << " sec. " << endl;
	};

// Free work arrays

	delete [] sp2blk;
	delete [] blk2cpu;
	delete [] imaskja;
	delete [] imaskblk;
	delete [] listblk;
	delete [] iablk;
	delete [] jablk;
	delete [] sndarr;
	delete [] rcvarr;

};

// Author: Kharchenko S.A. 
// Description: Update mask of the current jblk block elements of the complex A^hA matrix for rectangular A
// CSMatrix::UpdateAtAMatrixStruct()
//========================================================================================
int CSMatrix::UpdateAtAMatrixStruct (CSMatrix &_mtratasp, // Update mask of the current jblk block elements of the complex A^hA matrix for rectangular A
													int _iblk, int _jblk, 
													int *_blksc, int *_sp2blk, int *_imaskja) {

	const char *funcname = "UpdateAtAMatrixStruct";

// Cycle over the supernodes

	int isup, jsup, jblk;
	int ibegi, iendi, ibegj, iendj, j;
	int indi, indj, ksupi, ksupj;

	int nadd = 0;

	for (isup=_blksc[_iblk];isup<_blksc[_iblk+1];isup++) {
		for (j=_mtratasp.ia[isup];j<_mtratasp.ia[isup+1];j++) {
			jsup = _mtratasp.ja[j];
			jblk = _sp2blk[jsup];

			if (jblk == _jblk) {

// Scan both block columns

				ibegi = ia[isup];
				iendi = ia[isup+1]-1;
				ibegj = ia[jsup];
				iendj = ia[jsup+1]-1;

				indi = ibegi;
				indj = ibegj;

				while (indi <= iendi && indj <= iendj) {
					ksupi = ja[indi];
					ksupj = ja[indj];
					if (ksupi == ksupj) {
						nadd++;
						_imaskja[indj]++;
						indi++;
						indj++;
					} else if (ksupi < ksupj) {
						indi++;
					} else if (ksupi > ksupj) {
						indj++;
					};
				};

			};

		};
	};

	return nadd;

};

// Author: Kharchenko S.A. 
// Description: Compute current block of the complex A^hA matrix for rectangular A
// CGSMatrixCS::UpdateAtAMatrixBlock()
//========================================================================================
void CGSMatrixCS::UpdateAtAMatrixBlock (int _iblk, int _jblk, // Compute current block of the complex A^hA matrix for rectangular A
														CGSMatrixCS &_mtra, int *_sp2blk) {

	const char *funcname = "UpdateAtAMatrixBlock";

// Cycle over the supernodes

	int ilist, isup, ni, nj, nk;
	int kii, kjj, kkk, ibs, ibsi, ibsj, jsup;
	int ibegi, iendi, ibegj, iendj, j, jblkloc, jlist;
	int indi, indj, ksupi, ksupj;
	dcmplx caux;

	for (ilist=0;ilist<mtrarr[_iblk].nlist;ilist++) {
		isup = mtrarr[_iblk].list[ilist];
		ni = sprndc[isup+1]-sprndc[isup];
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jsup = mtrarr[_iblk].ja[j];
			ibs = mtrarr[_iblk].bsa[j];
			nj = sprndc[jsup+1]-sprndc[jsup];
			jblkloc = _sp2blk[jsup];
			if (jblkloc == _jblk) {

// Find block number

				jlist = jsup - blksc[_jblk];

// Scan both block columns

				ibegi = _mtra.mtrarr[_iblk].ia[ilist];
				iendi = _mtra.mtrarr[_iblk].ia[ilist+1]-1;
				ibegj = _mtra.mtrarr[_jblk].ia[jlist];
				iendj = _mtra.mtrarr[_jblk].ia[jlist+1]-1;

				indi = ibegi;
				indj = ibegj;

				while (indi <= iendi && indj <= iendj) {
					ksupi = _mtra.mtrarr[_iblk].ja[indi];
					ksupj = _mtra.mtrarr[_jblk].ja[indj];
					if (ksupi == ksupj) {
						ibsi = _mtra.mtrarr[_iblk].bsa[indi];
						ibsj = _mtra.mtrarr[_jblk].bsa[indj];
						nk = _mtra.sprndr[ksupi+1]-_mtra.sprndr[ksupi];
						for (kii=0;kii<ni;kii++) {
							for (kjj=0;kjj<nj;kjj++) {
								caux = mtrarr[_iblk].a[ibs+kjj*ni+kii];
								for (kkk=0;kkk<nk;kkk++) {
									caux += conj(_mtra.mtrarr[_iblk].a[ibsi+kkk*ni+kii]) * 
														_mtra.mtrarr[_jblk].a[ibsj+kkk*nj+kjj];
								};
								mtrarr[_iblk].a[ibs+kjj*ni+kii] = caux;
							};
						};
						indi++;
						indj++;
					} else if (ksupi < ksupj) {
						indi++;
					} else if (ksupi > ksupj) {
						indj++;
					};
				};

			};
		};
	};

};

// Author: Kharchenko S.A. 
// Description: Compute complex A^hA matrix in parallel mode
// CGSMatrixCS::AtAMatr()
//========================================================================================
void CGSMatrixCS::AtAMatr (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute complex A^hA matrix in parallel mode
							CGSMatrixCS &_mtral, CGSMatrixCS &_mtrau,
							CGSMatrixCS &_mtratau) { 

	const char *funcname = "AtAMatr";

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

// Init time measurement

	clock_t time0;

	time0 = clock ();

// Allocate the memory

	int *sp2blk;

	sp2blk = new int [_mtral.nsupr];
	if (!sp2blk) MemoryFail (funcname);

// Compute sp2blk array

	int i, j;
	for (i=0;i<_mtral.nblksc;i++) {
		for (j=_mtral.blksc[i];j<_mtral.blksc[i+1];j++) {
			sp2blk[j] = i;
		};
	};

// Compute global sparsity of local data

	CSMatrix alloc, auloc;

	alloc = _mtral.GlobalSparsity ();
	auloc = _mtrau.GlobalSparsity ();

// Transpose au sparsity

	CSMatrix autloc;

	autloc = auloc.TranspMtrMask ();

	CSMatrix mtrdummy;

	auloc = mtrdummy;

// Create global sparsity by supernode columns

	CSMatrix aloc;

	aloc = mtrdummy.CombineStructures (alloc, autloc);

	alloc = mtrdummy;
	autloc = mtrdummy;

// Create the block sparsity

	CSMatrix ablkloc;

	ablkloc = aloc.BlockSparsity (_mtral.nblksc,sp2blk);

// Compute locally Ai, Bi, Ci and D and store them by block columns

// Init global submatrices

	CGSMatrixCS gmtrai, gmtrbi, gmtrci, gmtrdloc;

	gmtrai = _mtral;
	gmtrbi = _mtral;
	gmtrci = _mtral;
	gmtrdloc = _mtral;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"AiMatr_",_tree.myid,"_");
		gmtrai.SetupFiles (_param.ooctyp, strbuff);
		sprintf(strbuff, "%s%s%d%s",_param.path,"BiMatr_",_tree.myid,"_");
		gmtrbi.SetupFiles (_param.ooctyp, strbuff);
		sprintf(strbuff, "%s%s%d%s",_param.path,"CiMatr_",_tree.myid,"_");
		gmtrci.SetupFiles (_param.ooctyp, strbuff);
		sprintf(strbuff, "%s%s%d%s",_param.path,"DLocMatr_",_tree.myid,"_");
		gmtrdloc.SetupFiles (_param.ooctyp, strbuff);

	};

// Store separately Ai, Bi, Ci and local D data by block columns

	CGSMatrixCS gmtrdummy;

//	gmtrdummy.StoreABCDBlockColumns (_tree, _mtral, _mtrau, 
//										aloc, ablkloc,
//										gmtrai, gmtrbi, gmtrci, gmtrdloc);

// Exchange D submatrices from different CPU's

//	????

// Compute locally Ci^hCi

// Compute on root processor D^hD

// Collect the results on root cpu

// Send partial results to the corresponding cpu's

// Finalize computation of ata

// Close files

	if (_param.ooctyp > 0) {

		gmtrai.CloseFiles();
		gmtrbi.CloseFiles();
		gmtrci.CloseFiles();
		gmtrdloc.CloseFiles();

	};

// Output statistics

// Free work arrays

	delete [] sp2blk;

};

