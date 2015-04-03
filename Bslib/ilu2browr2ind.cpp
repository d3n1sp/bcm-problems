//------------------------------------------------------------------------------------------------
// File: ilu2browr2ind.cpp
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

#include "globals.h"
#include "tree.h"
#include "fct.h"
#include "fctdiag.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "slvparam.h"
#include "mvm.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CFctDiagR: Init scaling for the current block row
//========================================================================================
void CFctDiagR::Ilu2ScaleBlkRow2Index (int _iblkrow, int _sclpttype, double _sclmin, // Init scaling for the current block row
										const CGSMatrix &_gmtra, const CSMatrixR &_mtrau,
										CFctR &_fct) {

	const char *funcname = "Ilu2ScaleBlkRow2Index";

	double done = 1.0e0;
	double dzero = 0.0e0;

// Determine the size of the arrays to be stored

	int indbeg, indend, nlistl;

	indbeg = _gmtra.blksr[_iblkrow];
	indend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

	int nz = nlistl;
	int nlocblk = nlistl;

// Allocate and register new data

	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc;
	double *dpivloc;

	sclptloc = new double [nlocblk];
	if (!sclptloc) MemoryFail (funcname);
	scllloc = new double [nz];
	if (!scllloc) MemoryFail (funcname);
	scluloc = new double [nz];
	if (!scluloc) MemoryFail (funcname);
	sclinvlloc = new double [nz];
	if (!sclinvlloc) MemoryFail (funcname);
	sclinvuloc = new double [nz];
	if (!sclinvuloc) MemoryFail (funcname);
	diaggloc = new double [nz];
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

	int ilist, j, kii;

	double *peig;

	double *paloc, *pscl;

	double daux;

	if (_sclpttype == 1) {

		nz = 0;

		for (ilist=0;ilist<nlistl;ilist++) {

// Search diagonal supernode

			j = _mtrau.ia[ilist];
			daux = _mtrau.a[j];
			if (daux < 0.0e0) daux = -daux;
			dpivloc[ilist] = daux;
			daux = sqrt(daux);
			daux = done / daux;
			sclptloc[ilist] = daux;

		};

	} else {
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = done;
		for (kii=0;kii<nlocblk;kii++) dpivloc[kii] = done;
	};

// Update hystogram

	hyst2.UpdateHystogram (nlocblk,dpivloc);

// Init scaling data for the current block row

	for (ilist=0;ilist<nlistl;ilist++) {

// Assign diagonal supernode

		j = _mtrau.ia[ilist];

		paloc = _fct.aloc;
		*paloc = _mtrau.a[j] * sclptloc[ilist] * sclptloc[ilist];

		if (ilist == 1354) {
			cout << " Wait " << endl;
		};

// Compute its singular value decomposition

		if (_fct.aloc[0] >= 0) {
			_fct.eig[0] = _fct.aloc[0];
			_fct.uloc[0] = 1.0e0;
			_fct.vloc[0] = 1.0e0;
		} else {
			_fct.eig[0] = -_fct.aloc[0];
			_fct.uloc[0] = -1.0e0;
			_fct.vloc[0] = 1.0e0;
		};

		_fct.ops += 8;

// Modify the pivots

		if (_fct.eig[0] < _sclmin) _fct.eig[0] = _sclmin;

// Store scaling factors

		daux = _fct.eig[0];
		daux = sqrt(daux);
		_fct.eig[0] = daux;
		dpivloc[ilist] = _fct.eig[0];

		peig = _fct.eig;
		paloc = _fct.uloc;
		pscl = sclinvlloc+ilist;
		*pscl = *paloc * *peig;

		peig = _fct.eig;
		paloc = _fct.vloc;
		pscl = sclinvuloc+ilist;
		*pscl = *paloc * *peig++;

		daux = _fct.eig[0];
		daux = 1.0e0 / daux;
		_fct.eig[0] = daux;

		peig = _fct.eig;
		paloc = _fct.uloc;
		pscl = scllloc+ilist;
		*pscl = *paloc * *peig;

		peig = _fct.eig;
		paloc = _fct.vloc;
		pscl = scluloc+ilist;
		*pscl = *paloc * *peig;

		_fct.ops += 6;

	};

// Init diagg

	for (kii=0;kii<nlistl;kii++) diaggloc[kii] = dzero;

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

		scllloc = new double [0];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new double [0];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new double [0];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new double [0];
		if (!sclinvuloc) MemoryFail (funcname);
		diaggloc = new double [0];
		if (!diaggloc) MemoryFail (funcname);

		scll   [_iblkrow] = scllloc;
		sclu   [_iblkrow] = scluloc;
		sclinvl[_iblkrow] = sclinvlloc;
		sclinvu[_iblkrow] = sclinvuloc;
		diagg  [_iblkrow] = diaggloc;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init current block row
//========================================================================================
void CFctR::Ilu2InitBlkRow2Index (int _iblkrow, // Init current block row
									double _dshift,
									const CSMatrixR &_mtral, const CSMatrixR &_mtrau, 
									CSMatrixR &_mtrl, CSMatrixR &_mtru) {

//	const char *funcname = "Ilu2InitBlkRow2Index";

// Check the data on entry

	int ilist, nlistl, irow, iblk;

	nlistl = _mtral.nlist;

	if (nlistl != _mtral.nlist2 || nlistl != _mtrau.nlist || nlistl != _mtrau.nlist2) {
		cout << " CFctR::Ilu2InitBlkRow2Index: Wrong number of supernode rows in the block row" << endl;
		assert (false); throw " CFctR::Ilu2InitBlkRow2Index: Wrong number of supernode rows in the block row";
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = _mtral.list[ilist];
		iblk = _mtral.list2[ilist];
		if (irow != ilist || iblk != _iblkrow) {
			cout << " CFctR::Ilu2InitBlkRow2Index: Wrong supernode row in the block row" << endl;
			assert (false); throw " CFctR::Ilu2InitBlkRow2Index: Wrong supernode row in the block row";
		};
		irow = _mtrau.list[ilist];
		iblk = _mtrau.list2[ilist];
		if (irow != ilist || iblk != _iblkrow) {
			cout << " CFctR::Ilu2InitBlkRow2Index: Wrong supernode row in the block row" << endl;
			assert (false); throw " CFctR::Ilu2InitBlkRow2Index: Wrong supernode row in the block row";
		};
	};

// Copy the data

	_mtrl = _mtral;
	_mtru = _mtrau;

	_mtrl.nsupr = _mtrl.nlist;
	_mtrl.nsupc = _mtrl.nlist;
	_mtru.nsupr = _mtru.nlist;
	_mtru.nsupc = _mtru.nlist;

	_mtrl.AllocateJa3 (_mtral.nzja);
	_mtru.AllocateJa3 (_mtrau.nzja);

	_mtrl.SetNzja3 (_mtrau.nzja);
	_mtru.SetNzja3 (_mtrau.nzja);

// Scale 

	int j, jcol, jblk, ibsd, jbsd;
	double aux;

	for (ilist=0;ilist<nlistl;ilist++) {

		irow = _mtrl.list[ilist];
		iblk = _mtrl.list2[ilist];

		for (j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {

			jcol = _mtrl.ja[j];
			jblk = _mtrl.ja2[j];

			ibsd = ibsdia[iblk]+irow;
			jbsd = ibsdia[jblk]+jcol;

// Assign Ja3

			_mtrl.ja3[j] = 1;
			_mtru.ja3[j] = 1;

// Local point scaling (U part)

			aux = _mtru.a[j] * sclpt[ibsd] * sclpt[jbsd] * scll[ibsd] * sclu[jbsd];
			_mtru.a[j] = aux;

// Local point scaling (L part)

			aux = _mtrl.a[j] * sclpt[ibsd] * sclpt[jbsd] * sclu[ibsd] * scll[jbsd];
			_mtrl.a[j] = aux;

// Shift diagonal data if necessary

			if (irow == jcol && iblk == jblk) _mtrl.a[j] += _dshift;
			if (irow == jcol && iblk == jblk) _mtru.a[j] += _dshift;

		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Init current block row
//========================================================================================
void CFctR::Ich2InitBlkRow2Index (int _iblkrow, // Init current block row
									double _dshift,
									const CSMatrixR &_mtral,
									CSMatrixR &_mtrl) {

//	const char *funcname = "Ilu2InitBlkRow2Index";

// Check the data on entry

	int ilist, nlistl, irow, iblk;

	nlistl = _mtral.nlist;

	if (nlistl != _mtral.nlist2) {
		cout << " CFctR::Ich2InitBlkRow2Index: Wrong number of supernode rows in the block row" << endl;
		assert (false); throw " CFctR::Ich2InitBlkRow2Index: Wrong number of supernode rows in the block row";
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = _mtral.list[ilist];
		iblk = _mtral.list2[ilist];
		if (irow != ilist || iblk != _iblkrow) {
			cout << " CFctR::Ich2InitBlkRow2Index: Wrong supernode row in the block row" << endl;
			assert (false); throw " CFctR::Ich2InitBlkRow2Index: Wrong supernode row in the block row";
		};
	};

// Copy the data

	_mtrl = _mtral;

	_mtrl.nsupr = _mtrl.nlist;
	_mtrl.nsupc = _mtrl.nlist;

	_mtrl.AllocateJa3 (_mtral.nzja);

	_mtrl.SetNzja3 (_mtral.nzja);

// Scale 

	int j, jcol, jblk, ibsd, jbsd;
	double aux;

	for (ilist=0;ilist<nlistl;ilist++) {

		irow = _mtrl.list[ilist];
		iblk = _mtrl.list2[ilist];

		for (j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {

			jcol = _mtrl.ja[j];
			jblk = _mtrl.ja2[j];

			ibsd = ibsdia[iblk]+irow;
			jbsd = ibsdia[jblk]+jcol;

// Assign Ja3

			_mtrl.ja3[j] = 1;

// Local point scaling (L part)

			aux = _mtrl.a[j] * sclpt[ibsd] * sclpt[jbsd] * sclu[ibsd] * scll[jbsd];
			_mtrl.a[j] = aux;

// Shift diagonal data if necessary

			if (irow == jcol && iblk == jblk) _mtrl.a[j] += _dshift;

		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Init current supernode row
//========================================================================================
void CFctR::Ilu2InitRow2Index (char _type, int _iblk, int _irow, // Init current supernode row
								const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau) {

//	const char *funcname = "Ilu2InitRow2Index";

	double dzero = 0.0e0;

	int j, ibsd, ibsi, ibsf, jbsi, jbsf, jbsd;

	double *paloc;
	double *pfmask;
	double *pvloc;

	icycle++;

	ibsd = ibsdia  [_iblk]+_irow;
	ibsi = ibsimask[_iblk]+_irow;
	ibsf = ibsfmask[_iblk]+_irow;

	nlist = 0;
	lstloc [nlist] = _irow;
	lstloc2[nlist] = _iblk;
	imask[ibsi] = icycle;

	pfmask = fmaskl+ibsf;
	pvloc = diagg+ibsd;
	*pfmask++ = dzero;

	j = _mtrau.ia[_irow];
	pfmask = fmasku+ibsf;
	pvloc = diagg+ibsd;
	paloc = _mtrau.a+j;
	if (_type == 'F') {
		*pfmask = *paloc + *pvloc;
	} else {
		*pfmask = *paloc ;
	};
	imasklev[ibsi] = _mtral.ja3[j];

	nlist++;

	for (j=_mtral.ia[_irow]+1;j<_mtral.ia[_irow+1];j++) {
		int jj = _mtral.ja[j];
		int jj2 = _mtral.ja2[j];
		lstloc[nlist] = jj;
		lstloc2[nlist] = jj2;

		jbsi = ibsimask[jj2]+jj;
		jbsf = ibsfmask[jj2]+jj;
		jbsd = ibsdia  [jj2]+jj;

		imask[jbsi] = icycle;

// Local copy

		imasklev[jbsi] = _mtral.ja3[j];

		fmasku[jbsf] = _mtrau.a[j];
		fmaskl[jbsf] = _mtral.a[j];

		nlist++;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init current supernode row
//========================================================================================
void CFctR::Ich2InitRow2Index (char _type, int _iblk, int _irow, // Init current supernode row
											const CGSMatrix &_gmtra, const CSMatrixR &_mtral) {

//	const char *funcname = "Ich2InitRow2Index";

	double dzero = 0.0e0;

	int j, ibsd, ibsi, ibsf, jbsi, jbsf, jbsd;

	double *paloc;
	double *pfmask;
	double *pvloc;

	icycle++;

	ibsd = ibsdia  [_iblk]+_irow;
	ibsi = ibsimask[_iblk]+_irow;
	ibsf = ibsfmask[_iblk]+_irow;

	nlist = 0;
	lstloc [nlist] = _irow;
	lstloc2[nlist] = _iblk;
	imask[ibsi] = icycle;

	pfmask = fmasku+ibsf;
	pvloc = diagg+ibsd;
	*pfmask++ = dzero;

	j = _mtral.ia[_irow];
	pfmask = fmaskl+ibsf;
	pvloc = diagg+ibsd;
	paloc = _mtral.a+j;
	if (_type == 'F') {
		*pfmask = *paloc + *pvloc;
	} else {
		*pfmask = *paloc ;
	};
	imasklev[ibsi] = _mtral.ja3[j];

	nlist++;

	for (j=_mtral.ia[_irow]+1;j<_mtral.ia[_irow+1];j++) {
		int jj = _mtral.ja[j];
		int jj2 = _mtral.ja2[j];
		lstloc[nlist] = jj;
		lstloc2[nlist] = jj2;

		jbsi = ibsimask[jj2]+jj;
		jbsf = ibsfmask[jj2]+jj;
		jbsd = ibsdia  [jj2]+jj;

		imask[jbsi] = icycle;

// Local copy

		imasklev[jbsi] = _mtral.ja3[j];

		fmaskl[jbsf] = _mtral.a[j];

		nlist++;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2InitMask2Index (bool _elmisr, int _ilstprv) { // Init mask and update arrays

//	const char *funcname = "Ilu2InitMask2Index";

	double dzero = 0.0e0;

	int jbsi, jbsf;

	double *glloc;
	double *guloc;
	int k, kloc;
	int *jgloc, *jg2loc, *jg3loc;
	int jcolmn, jcolblk, jlev, jlevnew;

	nupd = 0;

	for (k=iv[_ilstprv];k<ig[_ilstprv+1];k++) {

		kloc = k-ig[_ilstprv];

		jgloc = jg[_ilstprv];
		jcolmn = jgloc[kloc];
		jg2loc = jg2[_ilstprv];
		jcolblk = jg2loc[kloc];
		jg3loc = jg3[_ilstprv];
		jlev = jg3loc[kloc];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn-1;

		jbsi = ibsimask[jcolblk]+jcolmn;
		jbsf = ibsfmask[jcolblk]+jcolmn;

		jlevnew = jlev + elemlev;

		if (imask[jbsi] != icycle) {
			imask[jbsi] = icycle;
			imasklev[jbsi] = jlevnew;
			lstloc[nlist] = jcolmn;
			lstloc2[nlist] = jcolblk;
			fmaskl[jbsf] = dzero;
			fmasku[jbsf] = dzero;
			nlist++;
		} else {
			if (jlevnew < imasklev[jbsi]) imasklev[jbsi] = jlevnew;
		};
		lstupd[nupd] = jcolmn;
		lstupd2[nupd] = jcolblk;

//		glloc = gl[_ilstprv];
//		guloc = gu[_ilstprv];
//		adrupdl[nupd] = glloc+kloc;
//		adrupdu[nupd] = guloc+kloc;

		glloc = gl[_ilstprv]+kloc;
		guloc = gu[_ilstprv]+kloc;

		fmasku[jbsf] -= *lelem * *guloc;
		fmaskl[jbsf] -= *uelem * *glloc;

		nupd++;

label1:;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ich2InitMask2Index (bool _elmisr, int _ilstprv) { // Init mask and update arrays

//	const char *funcname = "Ilu2InitMask2Index";

	double dzero = 0.0e0;

	int jbsi, jbsf;

	double *glloc;
	int k, kloc;
	int *jgloc, *jg2loc, *jg3loc;
	int jcolmn, jcolblk, jlev, jlevnew;

	nupd = 0;

	for (k=iv[_ilstprv];k<ig[_ilstprv+1];k++) {

		kloc = k-ig[_ilstprv];

		jgloc = jg[_ilstprv];
		jcolmn = jgloc[kloc];
		jg2loc = jg2[_ilstprv];
		jcolblk = jg2loc[kloc];
		jg3loc = jg3[_ilstprv];
		jlev = jg3loc[kloc];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn-1;

		jbsi = ibsimask[jcolblk]+jcolmn;
		jbsf = ibsfmask[jcolblk]+jcolmn;

		jlevnew = jlev + elemlev;

		if (imask[jbsi] != icycle) {
			imask[jbsi] = icycle;
			imasklev[jbsi] = jlevnew;
			lstloc[nlist] = jcolmn;
			lstloc2[nlist] = jcolblk;
			fmaskl[jbsf] = dzero;
			nlist++;
		} else {
			if (jlevnew < imasklev[jbsi]) imasklev[jbsi] = jlevnew;
		};
		lstupd[nupd] = jcolmn;
		lstupd2[nupd] = jcolblk;

//		glloc = gl[_ilstprv];
//		guloc = gu[_ilstprv];
//		adrupdl[nupd] = glloc+kloc;
//		adrupdu[nupd] = guloc+kloc;

		glloc = gl[_ilstprv]+kloc;

		fmaskl[jbsf] -= *uelem * *glloc;

		nupd++;

label1:;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Update current row
//========================================================================================
void CFctR::Ilu2UpdateRow2Index () { // Update current row

//	const char *funcname = "Ilu2UpdateRow2Index";
/*
	int jbsf;

//	double *elemg; 
	int k, icol, icolblk;

	for (k=0;k<nupd;k++) {

		icol = lstupd[k];
		icolblk = lstupd2[k];

		jbsf = ibsfmask[icolblk]+icol;

//		elemg = adrupdu[k];
		fmasku[jbsf] -= *lelem * *adrupdu[k];

//		elemg = adrupdl[k];
		fmaskl[jbsf] -= *uelem * *adrupdl[k];

	};

	ops += 2*nupd;
*/
};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTransp2Index (int _iblk, int _ilist) { // Update arrays that support transposed structures search

//	const char *funcname = "Ilu2UpdateTransp2Index";

	int ilstprv = ibegm[_ilist];

	while (ilstprv != -1) {

		int ilstprv1 = madj[ilstprv];

		if (iv[ilstprv] >= ig[ilstprv+1]-1) {

			madj[ilstprv] = -1;

		} else {

			int j = iv[ilstprv]+1;
			int *jgloc = jg[ilstprv];
			int *jg2loc = jg2[ilstprv];
			int jloc = j-ig[ilstprv];
			int icolmn = jgloc[jloc];
			int icolblk = jg2loc[jloc];
			if (icolmn < 0) icolmn = -icolmn-1;
			if (icolblk != _iblk) {
				madj[ilstprv] = -1;
			} else {
				int icolmnl = icolmn;
				madj[ilstprv] = ibegm[icolmnl];
				ibegm[icolmnl] = ilstprv;
			};
		};

		iv[ilstprv]++;

		ilstprv = ilstprv1;

	};

};

// Author: Kharchenko S.A.
// CFctR: Perform filtering of the current row
//========================================================================================
void CFctR::Ilu2FiltrRow2Index (int _iblk, int _irow, // Perform filtering of the current row
								int _fcttyp, double _theta, double _tau2, int _ndiagfct2, int _nlevfct2, 
								const CSMatrix &_mtraini) {

//	const char *funcname = "Ilu2FiltrRow2Index";

	int nlist2 = 0;
	int ibsf, jbsf, jbsd, jbsi;
	int i, i1;
	int ibeg1, iend1, ibeg2, iend2, ip1, ip2, jj1, jj2, jjblk1, jjblk2, jlev;
	double aux, auxloc, auxloc1;

	ibsf = ibsfmask[_iblk]+_irow;

	if (_fcttyp == 0) {

		if (nlist > 0) qsort (lstloc+1, nlist-1, sizeof(int), compint);

		ibeg1 = 1;
		iend1 = nlist-1;
		ibeg2 = _mtraini.ia[_irow]+1;
		iend2 = _mtraini.ia[_irow+1]-1;
		ip1 = ibeg1;
		ip2 = ibeg2;
		while (ip1 <= iend1 && ip2 <= iend2) {
			jj1 = lstloc[ip1];
			jjblk1 = lstloc2[ip1];
			jj2 = _mtraini.ja[ip2];
			jjblk2 = _mtraini.ja2[ip2];
			if (jj1 == jj2 && jjblk1 == jjblk2) {
				lstupd[nlist2] = jj1;
				lstupd2[nlist2] = jjblk1;
				nlist2++;
				ip1++;
				ip2++;
			} else if (jjblk1 < jjblk2 || (jjblk1 == jjblk2 && jj1 < jj2)) {

				jbsd = ibsdia[jjblk1]+jj1;
				jbsf = ibsfmask[jjblk1]+jj1;

				aux = 0.0e0;

				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;
				auxloc = fmasku[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				auxloc = aux * _theta;
				fmasku[ibsf] += auxloc;
				diagg[jbsd] += auxloc;

				ip1++;
			} else {
				ip2++;
			};
		};

	} else {

		for (i=1;i<nlist;i++) {

			jj1 = lstloc[i];
			jjblk1 = lstloc2[i];

			jbsd = ibsdia[jjblk1]+jj1;
			jbsi = ibsimask[jjblk1]+jj1;
			jbsf = ibsfmask[jjblk1]+jj1;

			jlev = imasklev[jbsi];

			aux = 0.0e0;

			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;

			if (aux >= _tau2 && !(_fcttyp == 5 && jlev > _nlevfct2)) {
				lstupd[nlist2] = jj1;
				lstupd2[nlist2] = jjblk1;
				nlist2++;
			} else {
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc1 = fmasku[jbsf];
				if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
				if (auxloc1 > auxloc) auxloc = auxloc1;
				auxloc *= _theta;
				fmasku[ibsf] += auxloc;
				diagg[jbsd] += auxloc;
			};

		};

	};

	if ((_fcttyp == 3 || _fcttyp == 4) && nlist2 > _ndiagfct2) {
		int nlist3 = 0;
		int nlistsrt = 0;
		for (i=0;i<nlist2;i++) {
			jj1 = lstupd[i];
			jjblk1 = lstupd2[i];
			jbsd = ibsdia[jjblk1]+jj1;
			jbsf = ibsfmask[jjblk1]+jj1;
			if (jj1 == _irow && jjblk1 == _iblk) {
				lstloc[nlist3] = jj1;
				lstloc2[nlist3] = jjblk1;
				nlist3++;
			} else {
				aux = 0.0e0;
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;
				auxloc = fmasku[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				diarr[nlistsrt].dvalue = aux;
				diarr[nlistsrt].intvalue = i;
				nlistsrt++;
			};
		};

		if (nlistsrt > _ndiagfct2) {

			if (nlistsrt != 0) qsort (diarr, nlistsrt, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

			for (i=0;i<nlistsrt-_ndiagfct2;i++) {
				i1 = diarr[i].intvalue;
				jj1 = lstupd[i1];
				jjblk1 = lstupd2[i1];
				jbsd = ibsdia[jjblk1]+jj1;
				jbsf = ibsfmask[jjblk1]+jj1;
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc1 = fmasku[jbsf];
				if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
				if (auxloc1 > auxloc) auxloc = auxloc1;
				auxloc *= _theta;
				fmasku[ibsf] += auxloc;
				diagg[jbsd] += auxloc;
			};
			for (i=nlistsrt-_ndiagfct2;i<nlistsrt;i++) {
				i1 = diarr[i].intvalue;
				jj1 = lstupd[i1];
				jjblk1 = lstupd2[i1];
				lstloc[nlist3] = jj1;
				lstloc2[nlist3] = jjblk1;
				nlist3++;
			};

			nlist2 = nlist3;
			for (i=0;i<nlist2;i++) lstupd[i] = lstloc[i];
			for (i=0;i<nlist2;i++) lstupd2[i] = lstloc2[i];

		};
	};

	nlist = nlist2;

};

// Author: Kharchenko S.A.
// CFctR: Perform filtering of the current row
//========================================================================================
void CFctR::Ich2FiltrRow2Index (int _iblk, int _irow, // Perform filtering of the current row
								int _fcttyp, double _theta, double _tau2, int _ndiagfct2, int _nlevfct2, 
								const CSMatrix &_mtraini) {

//	const char *funcname = "Ich2FiltrRow2Index";

	int nlist2 = 0;
	int ibsf, jbsf, jbsd, jbsi;
	int i, i1;
	int ibeg1, iend1, ibeg2, iend2, ip1, ip2, jj1, jj2, jjblk1, jjblk2, jlev;
	double aux, auxloc;

	ibsf = ibsfmask[_iblk]+_irow;

	if (_fcttyp == 0) {

		if (nlist > 0) qsort (lstloc+1, nlist-1, sizeof(int), compint);

		ibeg1 = 1;
		iend1 = nlist-1;
		ibeg2 = _mtraini.ia[_irow]+1;
		iend2 = _mtraini.ia[_irow+1]-1;
		ip1 = ibeg1;
		ip2 = ibeg2;
		while (ip1 <= iend1 && ip2 <= iend2) {
			jj1 = lstloc[ip1];
			jjblk1 = lstloc2[ip1];
			jj2 = _mtraini.ja[ip2];
			jjblk2 = _mtraini.ja2[ip2];
			if (jj1 == jj2 && jjblk1 == jjblk2) {
				lstupd[nlist2] = jj1;
				lstupd2[nlist2] = jjblk1;
				nlist2++;
				ip1++;
				ip2++;
			} else if (jjblk1 < jjblk2 || (jjblk1 == jjblk2 && jj1 < jj2)) {

				jbsd = ibsdia[jjblk1]+jj1;
				jbsf = ibsfmask[jjblk1]+jj1;

				aux = 0.0e0;

				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				auxloc = aux * _theta;
				fmaskl[ibsf] += auxloc;
				diagg[jbsd] += auxloc;

				ip1++;
			} else {
				ip2++;
			};
		};

	} else {

		for (i=1;i<nlist;i++) {

			jj1 = lstloc[i];
			jjblk1 = lstloc2[i];

			jbsd = ibsdia[jjblk1]+jj1;
			jbsi = ibsimask[jjblk1]+jj1;
			jbsf = ibsfmask[jjblk1]+jj1;

			jlev = imasklev[jbsi];

			aux = 0.0e0;

			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;

			if (aux >= _tau2 && !(_fcttyp == 5 && jlev > _nlevfct2)) {
				lstupd[nlist2] = jj1;
				lstupd2[nlist2] = jjblk1;
				nlist2++;
			} else {
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc *= _theta;
				fmaskl[ibsf] += auxloc;
				diagg[jbsd] += auxloc;
			};

		};

	};

	if ((_fcttyp == 3 || _fcttyp == 4) && nlist2 > _ndiagfct2) {
		int nlist3 = 0;
		int nlistsrt = 0;
		for (i=0;i<nlist2;i++) {
			jj1 = lstupd[i];
			jjblk1 = lstupd2[i];
			jbsd = ibsdia[jjblk1]+jj1;
			jbsf = ibsfmask[jjblk1]+jj1;
			if (jj1 == _irow && jjblk1 == _iblk) {
				lstloc[nlist3] = jj1;
				lstloc2[nlist3] = jjblk1;
				nlist3++;
			} else {
				aux = 0.0e0;
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				diarr[nlistsrt].dvalue = aux;
				diarr[nlistsrt].intvalue = i;
				nlistsrt++;
			};
		};

		if (nlistsrt > _ndiagfct2) {

			if (nlistsrt != 0) qsort (diarr, nlistsrt, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

			for (i=0;i<nlistsrt-_ndiagfct2;i++) {
				i1 = diarr[i].intvalue;
				jj1 = lstupd[i1];
				jjblk1 = lstupd2[i1];
				jbsd = ibsdia[jjblk1]+jj1;
				jbsf = ibsfmask[jjblk1]+jj1;
				auxloc = fmaskl[jbsf];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc *= _theta;
				fmaskl[ibsf] += auxloc;
				diagg[jbsd] += auxloc;
			};
			for (i=nlistsrt-_ndiagfct2;i<nlistsrt;i++) {
				i1 = diarr[i].intvalue;
				jj1 = lstupd[i1];
				jjblk1 = lstupd2[i1];
				lstloc[nlist3] = jj1;
				lstloc2[nlist3] = jjblk1;
				nlist3++;
			};

			nlist2 = nlist3;
			for (i=0;i<nlist2;i++) lstupd[i] = lstloc[i];
			for (i=0;i<nlist2;i++) lstupd2[i] = lstloc2[i];

		};
	};

	nlist = nlist2;

};

// Author: Kharchenko S.A.
// CFctR: Compute diagonal pivot and scale current row
//========================================================================================
void CFctR::Ilu2PivScaleRow2Index (int _iblk, int _irow, double _pivmin) { // Compute diagonal pivot and scale current row

//	const char *funcname = "Ilu2PivScaleRow2Index";

	int ibsf, jbsf;

// Copy diagonal supernode

	ibsf = ibsfmask[_iblk] + _irow;

	aloc[0] = fmasku[ibsf];

// Compute singular value decomposition of the diagonal supernode

	if (aloc[0] >= 0) {
		eig[0] = aloc[0];
		uloc[0] = 1.0e0;
		vloc[0] = 1.0e0;
	} else {
		eig[0] = -aloc[0];
		uloc[0] = -1.0e0;
		vloc[0] = 1.0e0;
	};

	ops += 1;

// Modify singular values if necessary

	eig[0] = sqrt(eig[0]);
	if (eig[0] < _pivmin) {
		eig[0] = _pivmin;
		nmodsv++;
	};
	dpiv[ibsf] = eig[0];
	eig[0] = 1.0e0/eig[0];

// Compute diagonal supernodes

	double *pdeloc;
	double *palloc;

	pdeloc = deleml;
	palloc = uloc;
	double temp = *palloc++;
	*pdeloc = eig[0] * temp;

	pdeloc = delemu;
	palloc = vloc;
	temp = *palloc;
	*pdeloc = eig[0] * temp;

	ops += 2;

// Perform scaling of the data in the list

	int i, j, jblk;

	for (i=0;i<nlist;i++) {

		j = lstupd[i];
		jblk = lstupd2[i];

		jbsf = ibsfmask[jblk] + j;

		fmaskl[jbsf] = *delemu * fmaskl[jbsf];

		fmasku[jbsf] = *deleml * fmasku[jbsf];

	};

	ops += 2*nlist;

};

// Author: Kharchenko S.A.
// CFctR: Compute diagonal pivot and scale current row
//========================================================================================
void CFctR::Ich2PivScaleRow2Index (int _iblk, int _irow, double _pivmin) { // Compute diagonal pivot and scale current row

//	const char *funcname = "Ich2PivScaleRow2Index";

	int ibsf, jbsf;

// Copy diagonal supernode

	ibsf = ibsfmask[_iblk] + _irow;

	aloc[0] = fmaskl[ibsf];

// Compute singular value decomposition of the diagonal supernode

	if (aloc[0] >= 0) {
		eig[0] = aloc[0];
		uloc[0] = 1.0e0;
		vloc[0] = 1.0e0;
	} else {
		eig[0] = -aloc[0];
		uloc[0] = -1.0e0;
		vloc[0] = 1.0e0;
	};

	ops += 1;

// Modify singular values if necessary

	eig[0] = sqrt(eig[0]);
	if (eig[0] < _pivmin) {
		eig[0] = _pivmin;
		nmodsv++;
	};
	dpiv[ibsf] = eig[0];
	eig[0] = 1.0e0/eig[0];

// Compute diagonal supernodes

	double *pdeloc;
	double *palloc;

	pdeloc = deleml;
	palloc = uloc;
	double temp = *palloc++;
	*pdeloc = eig[0] * temp;

	pdeloc = delemu;
	palloc = vloc;
	temp = *palloc;
	*pdeloc = eig[0] * temp;

	ops += 2;

// Perform scaling of the data in the list

	int i, j, jblk;

	for (i=0;i<nlist;i++) {

		j = lstupd[i];
		jblk = lstupd2[i];

		jbsf = ibsfmask[jblk] + j;

		fmaskl[jbsf] = *delemu * fmaskl[jbsf];

	};

	ops += 2*nlist;

};

// Author: Kharchenko S.A.
// CFctR: Store current row
//========================================================================================
void CFctR::Ilu2StoreRow2Index (int _iblk, int _irow, // Store current row
									int _fcttyp, double _tau1, int _ndiagfct, int _nlevfct, 
									int &_nzrtot) {

	const char *funcname = "Ilu2StoreRow2Index";

// Allocate the memory

	int i, j, jbsi, jlev;

	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		jg2loc = new int [nzalloc];
		if (!jg2loc) MemoryFail (funcname);
		jg3loc = new int [nzalloc];
		if (!jg3loc) MemoryFail (funcname);
		adrgloc = new int [0];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [nzalloc];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		jg2blk[iblkalloc]  = jg2loc;
		jg3blk[iblkalloc]  = jg3loc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
//		cout << " Iblkalloc = " << iblkalloc << " Iblk = " << _iblk << " Irow = " << _irow << endl;
//		cout << " Control: left = " << jgloc[-1] << " right = " << jgloc[nzalloc] << endl;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		jg2loc  = jg2blk[iblkalloc]+nzused;
		jg3loc  = jg3blk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc];
		glloc   = glblk[iblkalloc]+nzused;
		guloc   = gublk[iblkalloc]+nzused;
	};

	jg[_irow] = jgloc;
	jg2[_irow] = jg2loc;
	jg3[_irow] = jg3loc;
	addrg[_irow] = adrgloc;
	gl[_irow] = glloc;
	gu[_irow] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

// Store current diagonal supernode

	double *pd, *pg;

	jgloc[nzjgloc] = _irow;
	jg2loc[nzjgloc] = _iblk;
	jbsi = ibsimask[_iblk] + _irow;
	jg3loc[nzjgloc] = imasklev[jbsi];
	pd = deleml;
	pg = glloc+nzloc;
	*pg = *pd;
	pd = delemu;
	pg = guloc+nzloc;
	*pg = *pd;
//	adrgloc[nzjgloc] = nzloc;
	nzjgloc++;
	nzloc += 1;

	int jblk, jbsf;
	double aux, auxloc;

	if ((_fcttyp != 3 && _fcttyp != 4) || nlist < _ndiagfct) {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jlev = imasklev[jbsi];
			jbsf = ibsfmask[jblk] + j;
			aux = 0.0e0;
			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			if (_fcttyp == 0 || (aux >= _tau1 && !(_fcttyp == 5 && jlev > _nlevfct))) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j-1;
				_nzrtot += 1;
			};
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			glloc[nzloc] = fmaskl[jbsf];
			guloc[nzloc] = fmasku[jbsf];

			nzjgloc++;
			nzloc += 1;
		};

	} else {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jbsf = ibsfmask[jblk] + j;
			aux = 0.0e0;
			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			diarr[i].dvalue = aux;
			diarr[i].intvalue = i;
		};

		if (nlist != 0) qsort (diarr, nlist, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

		for (i=0;i<nlist;i++) lstloc[diarr[i].intvalue] = i;

		int i1, iord;

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jlev = imasklev[jbsi];
			jbsf = ibsfmask[jblk] + j;
			iord = lstloc[i];
			i1 = nlist-iord-1;
			if (i1 < _ndiagfct) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j-1;
				_nzrtot += 1;
			};
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			glloc[nzloc] = fmaskl[jbsf];
			guloc[nzloc] = fmasku[jbsf];

			nzjgloc++;
			nzloc += 1;

		};

	};

	ig[_irow+1] = ig[_irow]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store current row
//========================================================================================
void CFctR::Ich2StoreRow2Index (int _iblk, int _irow, // Store current row
									int _fcttyp, double _tau1, int _ndiagfct, int _nlevfct, 
									int &_nzrtot) {

	const char *funcname = "Ich2StoreRow2Index";

// Allocate the memory

	int i, j, jbsi, jlev;

	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		jg2loc = new int [nzalloc];
		if (!jg2loc) MemoryFail (funcname);
		jg3loc = new int [nzalloc];
		if (!jg3loc) MemoryFail (funcname);
		adrgloc = new int [0];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [0];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		jg2blk[iblkalloc]  = jg2loc;
		jg3blk[iblkalloc]  = jg3loc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
//		cout << " Iblkalloc = " << iblkalloc << " Iblk = " << _iblk << " Irow = " << _irow << endl;
//		cout << " Control: left = " << jgloc[-1] << " right = " << jgloc[nzalloc] << endl;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		jg2loc  = jg2blk[iblkalloc]+nzused;
		jg3loc  = jg3blk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc];
		glloc   = glblk[iblkalloc]+nzused;
	};

	jg[_irow] = jgloc;
	jg2[_irow] = jg2loc;
	jg3[_irow] = jg3loc;
	addrg[_irow] = adrgloc;
	gl[_irow] = glloc;

	int nzjgloc = 0;
	int nzloc = 0;

// Store current diagonal supernode

	double *pd, *pg;

	jgloc[nzjgloc] = _irow;
	jg2loc[nzjgloc] = _iblk;
	jbsi = ibsimask[_iblk] + _irow;
	jg3loc[nzjgloc] = imasklev[jbsi];
	pd = deleml;
	pg = glloc+nzloc;
	*pg = *pd;
	pd = delemu;
	pg = glloc+nzloc;
	*pg = *pd;
//	adrgloc[nzjgloc] = nzloc;
	nzjgloc++;
	nzloc += 1;

	int jblk, jbsf;
	double aux, auxloc;

	if ((_fcttyp != 3 && _fcttyp != 4) || nlist < _ndiagfct) {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jlev = imasklev[jbsi];
			jbsf = ibsfmask[jblk] + j;
			aux = 0.0e0;
			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			if (_fcttyp == 0 || (aux >= _tau1 && !(_fcttyp == 5 && jlev > _nlevfct))) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j-1;
				_nzrtot += 1;
			};
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			glloc[nzloc] = fmaskl[jbsf];

			nzjgloc++;
			nzloc += 1;
		};

	} else {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jbsf = ibsfmask[jblk] + j;
			aux = 0.0e0;
			auxloc = fmaskl[jbsf];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			diarr[i].dvalue = aux;
			diarr[i].intvalue = i;
		};

		if (nlist != 0) qsort (diarr, nlist, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

		for (i=0;i<nlist;i++) lstloc[diarr[i].intvalue] = i;

		int i1, iord;

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jblk = lstupd2[i];
			jbsi = ibsimask[jblk] + j;
			jlev = imasklev[jbsi];
			jbsf = ibsfmask[jblk] + j;
			iord = lstloc[i];
			i1 = nlist-iord-1;
			if (i1 < _ndiagfct) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j-1;
				_nzrtot += 1;
			};
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			glloc[nzloc] = fmaskl[jbsf];

			nzjgloc++;
			nzloc += 1;

		};

	};

	ig[_irow+1] = ig[_irow]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2StoreScaleU2Index (int _nlistl, int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
									CSMatrixR &_mtrl, CSMatrixR &_mtru,
									CSMatrixR &_mtrl2, CSMatrixR &_mtru2) {

//	const char *funcname = "Ilu2StoreScaleU2Index";

// Count the numbers of elements in the matrices to be stored

	int nlistl;

	nlistl = _nlistl;

	int *jgloc;
	int *jg2loc;

	int ilist, irow, j, jloc, jj, jjloc, ibeg, jjblk;

	int nzjlu=0, nzlu=0, nzjlu2=0, nzlu2=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			if (jj >= 0) {
				nzjlu++;
				nzlu++;
			};
			if (jjblk > _iblkrow) {
				nzjlu2++;
				nzlu2++;
			};
		};
	};

// Allocate and init first order matrices

	CSMatrixR mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixR temp (nlistl,nlistl,nzjlu,nzjlu,nzlu);

	temp.m = nlistl;
	temp.n = nlistl;

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			if (jj >= 0) {
				temp.ja[nzjlu] = jjloc;
				temp.ja2[nzjlu] = jjblk;
				nzjlu++;
				nzlu++;
			};
		};
		temp.list[ilist] = irow;
		temp.list2[ilist] = _iblkrow;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nlist = nlistl;
	temp.nlist2 = nlistl;
	temp.nzja = nzjlu;
	temp.nzja2 = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	double *glloc;
	double *guloc;

	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		jgloc = jg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			if (jj >= 0) {
				_mtrl.a[nzlu] = glloc[jloc];
				_mtru.a[nzlu] = guloc[jloc];
				nzlu++;
			};
		};
	};

// Scale

	if (_param.sclpttype < 2) {

		int ibsd, jbsd;

		for (ilist=0;ilist<nlistl;ilist++) {
			irow = ilist;
			ibsd = ibsdia[_iblkrow]+irow;
			for (int j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {
				jj = _mtrl.ja[j];
				jjblk = _mtrl.ja2[j];
				jbsd = ibsdia[jjblk]+jj;
				if (jj == irow && jjblk == _iblkrow) {

					_mtrl.a[j] = _mtrl.a[j] * scll[ibsd];
					_mtru.a[j] = _mtru.a[j] * sclu[ibsd];

	 				ops += 2;

				} else {

					_mtrl.a[j] = _mtrl.a[j] * sclinvl[jbsd];
					_mtru.a[j] = _mtru.a[j] * sclinvu[jbsd];

					ops += 2;
				};
			};
		};

	};

// Allocate and init second order matrices

	_mtrl2 = mtrdummy;
	_mtru2 = mtrdummy;

	CSMatrixR temp2 (nlistl,nlistl,nzjlu2,nzjlu2,nzlu2);

	temp2.m = nlistl;
	temp2.n = nlistl;

	temp2.ia[0] = 0;

	temp2.AllocateJa3 (nzjlu2);

	nzjlu2 = 0;
	nzlu2 = 0;

	int *jg3loc;
	int jlev;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		jg3loc = jg3[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jlev = jg3loc[jloc];
			if (jjblk > _iblkrow) {
				temp2.ja[nzjlu2] = jj;
				temp2.ja2[nzjlu2] = jjblk;
				temp2.ja3[nzjlu2] = jlev;
				nzjlu2++;
				nzlu2++;
			};
		};
		temp2.list[ilist] = irow;
		temp2.list2[ilist] = _iblkrow;
		temp2.ia[ilist+1] = nzjlu2;
	};

	int nlistn=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = temp2.list[ilist];
		if (temp2.ia[ilist+1] > temp2.ia[ilist]) {
			temp2.list[nlistn] = irow;
			temp2.list2[nlistn] = _iblkrow;
			temp2.ia[nlistn+1] = temp2.ia[ilist+1];
			nlistn++;
		};
	};

	temp2.nlist = nlistn;
	temp2.nlist2 = nlistn;
	temp2.nzja = nzjlu2;
	temp2.nzja2 = nzjlu2;
	temp2.nzja3 = nzjlu2;
	temp2.nza = nzlu2;
	temp2.nzatot = nzlu2;

	_mtrl2 = temp2;
	_mtru2 = temp2;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		jg2loc = jg2[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jjblk = jg2loc[jloc];
			if (jjblk > _iblkrow) {
				_mtrl2.a[nzlu2] = glloc[jloc];
				_mtru2.a[nzlu2] = guloc[jloc];
				nzjlu2++;
				nzlu2++;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ich2StoreScaleU2Index (int _nlistl, int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
									CSMatrixR &_mtrl,
									CSMatrixR &_mtrl2) {

//	const char *funcname = "Ich2StoreScaleU2Index";

// Count the numbers of elements in the matrices to be stored

	int nlistl;

	nlistl = _nlistl;

	int *jgloc;
	int *jg2loc;

	int ilist, irow, j, jloc, jj, jjloc, ibeg, jjblk;

	int nzjlu=0, nzlu=0, nzjlu2=0, nzlu2=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			if (jj >= 0) {
				nzjlu++;
				nzlu++;
			};
			if (jjblk > _iblkrow) {
				nzjlu2++;
				nzlu2++;
			};
		};
	};

// Allocate and init first order matrices

	CSMatrixR mtrdummy;

	_mtrl = mtrdummy;

	CSMatrixR temp (nlistl,nlistl,nzjlu,nzjlu,nzlu);

	temp.m = nlistl;
	temp.n = nlistl;

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			if (jj >= 0) {
				temp.ja[nzjlu] = jjloc;
				temp.ja2[nzjlu] = jjblk;
				nzjlu++;
				nzlu++;
			};
		};
		temp.list[ilist] = irow;
		temp.list2[ilist] = _iblkrow;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nlist = nlistl;
	temp.nlist2 = nlistl;
	temp.nzja = nzjlu;
	temp.nzja2 = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtrl = temp;

	temp = mtrdummy;

	double *glloc;

	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		jgloc = jg[ilist];
		glloc = gl[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			if (jj >= 0) {
				_mtrl.a[nzlu] = glloc[jloc];
				nzlu++;
			};
		};
	};

// Scale

	if (_param.sclpttype < 2) {

		int ibsd, jbsd;

		for (ilist=0;ilist<nlistl;ilist++) {
			irow = ilist;
			ibsd = ibsdia[_iblkrow]+irow;
			for (int j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {
				jj = _mtrl.ja[j];
				jjblk = _mtrl.ja2[j];
				jbsd = ibsdia[jjblk]+jj;
				if (jj == irow && jjblk == _iblkrow) {

					_mtrl.a[j] = _mtrl.a[j] * scll[ibsd];

	 				ops += 2;

				} else {

					_mtrl.a[j] = _mtrl.a[j] * sclinvl[jbsd];

					ops += 2;
				};
			};
		};

	};

// Allocate and init second order matrices

	_mtrl2 = mtrdummy;

	CSMatrixR temp2 (nlistl,nlistl,nzjlu2,nzjlu2,nzlu2);

	temp2.m = nlistl;
	temp2.n = nlistl;

	temp2.ia[0] = 0;

	temp2.AllocateJa3 (nzjlu2);

	nzjlu2 = 0;
	nzlu2 = 0;

	int *jg3loc;
	int jlev;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		jg3loc = jg3[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jlev = jg3loc[jloc];
			if (jjblk > _iblkrow) {
				temp2.ja[nzjlu2] = jj;
				temp2.ja2[nzjlu2] = jjblk;
				temp2.ja3[nzjlu2] = jlev;
				nzjlu2++;
				nzlu2++;
			};
		};
		temp2.list[ilist] = irow;
		temp2.list2[ilist] = _iblkrow;
		temp2.ia[ilist+1] = nzjlu2;
	};

	int nlistn=0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = temp2.list[ilist];
		if (temp2.ia[ilist+1] > temp2.ia[ilist]) {
			temp2.list[nlistn] = irow;
			temp2.list2[nlistn] = _iblkrow;
			temp2.ia[nlistn+1] = temp2.ia[ilist+1];
			nlistn++;
		};
	};

	temp2.nlist = nlistn;
	temp2.nlist2 = nlistn;
	temp2.nzja = nzjlu2;
	temp2.nzja2 = nzjlu2;
	temp2.nzja3 = nzjlu2;
	temp2.nza = nzlu2;
	temp2.nzatot = nzlu2;

	_mtrl2 = temp2;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		jg2loc = jg2[ilist];
		glloc = gl[ilist];
		ibeg = ig[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ibeg;
			jjblk = jg2loc[jloc];
			if (jjblk > _iblkrow) {
				_mtrl2.a[nzlu2] = glloc[jloc];
				nzjlu2++;
				nzlu2++;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform factorization of the current block row
//========================================================================================
void CFctR::Ilu2FctBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau, const CSMatrixR &_mtraini,
									CSMatrixR &_mtrl, CSMatrixR &_mtru,
									CSMatrixR &_mtrl2, CSMatrixR &_mtru2) {

//	const char *funcname = "Ilu2FctBlkRow2Index";

// Main cycle over the rows

	int ilist, irow;
	int nzrtot;
	int irwprv, ilstprv, icolmn, icolblk, j, jloc, ibs, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;
	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	double *glloc, *guloc;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

	nzrtot = 0;

	for (ilist=0;ilist<=_mtral.nlist;ilist++) ibegm[ilist] = -1;

	ig[0] = 0;

	int iblk;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];
		iblk = _mtral.list2[ilist];

// Init current scaled row

		Ilu2InitRow2Index ('F', iblk, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			jgloc = jg[ilstprv];
			jg2loc = jg2[ilstprv];
			jg3loc = jg3[ilstprv];
			glloc = gl[ilstprv];
			guloc = gu[ilstprv];

			j = iv[ilstprv];

			jloc = j-ig[ilstprv];

			icolmn = jgloc[jloc];
			icolblk = jg2loc[jloc];

			ibs = jloc;

			irwprv = _mtral.list[ilstprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl-1;

			elemlev = jg3loc[jloc];
			lelem[0] = glloc[ibs];
			uelem[0] = guloc[ibs];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2InitMask2Index (elmisr, ilstprv);

// Perform updates

			Ilu2UpdateRow2Index ();

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp2Index (_iblkrow, ilist);

// Perform filtering of the row with diagonal modification if necessary 

		Ilu2FiltrRow2Index (iblk, irow, 
									_param.fcttyp, _param.theta, _param.tau2, _param.ndiagfct2, _param.nlevfct2,
									_mtraini);

// Compute and store the pivot value and scale current supernode row

		Ilu2PivScaleRow2Index (iblk, irow, _param.pivmin);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstupd[i];
			jjblk = lstupd2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstupd[i] = ind2arr[i].indy;
			lstupd2[i] = ind2arr[i].indx;
		};

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist != 0 && lstupd2[0] == _iblkrow) {
			iv[ilist] = ig[ilist]+1;
			madj[ilist] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = ilist;
		};

// Store the data in the global arrays

		Ilu2StoreRow2Index (iblk, irow, 
								_param.fcttyp, _param.tau1, _param.ndiagfct, _param.nlevfct, 
								nzrtot);

	};

// Store the resulting U objects

	Ilu2StoreScaleU2Index (nlistl, _iblkrow, _param, 
							_mtrl, _mtru,
							_mtrl2, _mtru2);

	nzr += nzrtot;

// Store the hystogram of values

	int nloc = nlistl;

	double *dpivarr;

	dpivarr = new double [nloc];

	nloc = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = _mtral.list[ilist];
		int ibsd = ibsdia[_iblkrow]+irow;
		dpivarr[nloc] = dpiv[ibsd];
		nloc++;
	};

	hyst.UpdateHystogram (nloc, dpivarr);

	delete [] dpivarr;

// Free local factorization arrays 

	for (ilist=0;ilist<=iblkalloc;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform factorization of the current block row
//========================================================================================
void CFctR::Ich2FctBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtraini,
									CSMatrixR &_mtrl,
									CSMatrixR &_mtrl2) {

//	const char *funcname = "Ich2FctBlkRow2Index";

// Main cycle over the rows

	int ilist, irow;
	int nzrtot;
	int irwprv, ilstprv, icolmn, icolblk, j, jloc, ibs, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;
	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	double *glloc, *guloc;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

	nzrtot = 0;

	for (ilist=0;ilist<=_mtral.nlist;ilist++) ibegm[ilist] = -1;

	ig[0] = 0;

	int iblk;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];
		iblk = _mtral.list2[ilist];

// Init current scaled row

		Ich2InitRow2Index ('F', iblk, irow, _gmtra, _mtral);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			jgloc = jg[ilstprv];
			jg2loc = jg2[ilstprv];
			jg3loc = jg3[ilstprv];
			glloc = gl[ilstprv];
			guloc = gu[ilstprv];

			j = iv[ilstprv];

			jloc = j-ig[ilstprv];

			icolmn = jgloc[jloc];
			icolblk = jg2loc[jloc];

			ibs = jloc;

			irwprv = _mtral.list[ilstprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl-1;

			elemlev = jg3loc[jloc];
			lelem[0] = glloc[ibs];
			uelem[0] = lelem[0];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask2Index (elmisr, ilstprv);

// Perform updates

			Ilu2UpdateRow2Index ();

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp2Index (_iblkrow, ilist);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow2Index (iblk, irow, 
									_param.fcttyp, _param.theta, _param.tau2, _param.ndiagfct2, _param.nlevfct2,
									_mtraini);

// Compute and store the pivot value and scale current supernode row

		Ich2PivScaleRow2Index (iblk, irow, _param.pivmin);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstupd[i];
			jjblk = lstupd2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstupd[i] = ind2arr[i].indy;
			lstupd2[i] = ind2arr[i].indx;
		};

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist != 0 && lstupd2[0] == _iblkrow) {
			iv[ilist] = ig[ilist]+1;
			madj[ilist] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = ilist;
		};

// Store the data in the global arrays

		Ich2StoreRow2Index (iblk, irow, 
								_param.fcttyp, _param.tau1, _param.ndiagfct, _param.nlevfct, 
								nzrtot);

	};

// Store the resulting U objects

	Ich2StoreScaleU2Index (nlistl, _iblkrow, _param, 
							_mtrl,
							_mtrl2);

	nzr += nzrtot;

// Store the hystogram of values

	int nloc = nlistl;

	double *dpivarr;

	dpivarr = new double [nloc];

	nloc = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = _mtral.list[ilist];
		int ibsd = ibsdia[_iblkrow]+irow;
		dpivarr[nloc] = dpiv[ibsd];
		nloc++;
	};

	hyst.UpdateHystogram (nloc, dpivarr);

	delete [] dpivarr;

// Free local factorization arrays 

	for (ilist=0;ilist<=iblkalloc;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Init arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2InitTransp2Index (int _iblk, int _nlist, // Init arrays that support transposed structures search
							const CSMatrixR &_mtrl2) { 

//	const char *funcname = "Ilu2InitTransp2Index";

	int ilist, irow, j, jj, jjblk, nlistl, ilistn;

	for (ilist=0;ilist<_mtrl2.nlist;ilist++) {
		irow = _mtrl2.list[ilist];
		iv[ilist] = _mtrl2.ia[ilist+1];
		for (j=_mtrl2.ia[ilist];j<_mtrl2.ia[ilist+1];j++) {
			jj = _mtrl2.ja[j];
			jjblk = _mtrl2.ja2[j];
			if (jj < 0) jj = -jj-1;
			if (jjblk == _iblk) {
				iv[ilist] = j;
				break;
			};
		};
	};

	nlistl = _nlist;

	for (ilist=0;ilist<nlistl;ilist++) {
		ibegm[ilist] = -1;
	};

	for (ilist=0;ilist<_mtrl2.nlist;ilist++) {
		madj[ilist] = -1;
	};

	for (ilist=0;ilist<_mtrl2.nlist;ilist++) {
		j = iv[ilist];
		if (j<_mtrl2.ia[ilist+1]) {
			jj = _mtrl2.ja[j];
			jjblk = _mtrl2.ja2[j];
			if (jj < 0) jj = -jj-1;
			if (jjblk == _iblk) {
				ilistn = jj;
				madj[ilist] = ibegm[ilistn];
				ibegm[ilistn] = ilist;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, // Init mask and update arrays
								const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2) {

//	const char *funcname = "Ilu2InitMask2Index";

	double dzero = 0.0e0;

	int jbsi, jbsf;

	int k;
	int jcolmn, jcolblk, ilev, jlev, jlevnew;
	double *glloc, *guloc;

	nupd = 0;

	k = iv[_ibs+_ilstprv];

	ilev = _mtrl2.ja3[k];

	for (k=iv[_ibs+_ilstprv];k<_mtrl2.ia[_ilstprv+1];k++) {

		jcolmn = _mtrl2.ja[k];
		jcolblk = _mtrl2.ja2[k];
		jlev = _mtrl2.ja3[k];

		jlevnew = ilev + jlev;

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn-1;

		jbsi = ibsimask[jcolblk]+jcolmn;
		jbsf = ibsfmask[jcolblk]+jcolmn;

		if (imask[jbsi] != icycle) {
			imask[jbsi] = icycle;
			imasklev[jbsi] = jlevnew;
			lstloc[nlist] = jcolmn;
			lstloc2[nlist] = jcolblk;
			fmaskl[jbsf] = dzero;
			fmasku[jbsf] = dzero;
			nlist++;
		} else {
			if (jlevnew < imasklev[jbsi]) imasklev[jbsi] = jlevnew;
		};

		lstupd[nupd] = jcolmn;
		lstupd2[nupd] = jcolblk;

//		adrupdl[nupd] = _mtrl2.a+k;
//		adrupdu[nupd] = _mtru2.a+k;

		glloc = _mtrl2.a+k;
		guloc = _mtru2.a+k;

		fmasku[jbsf] -= *lelem * *guloc;
		fmaskl[jbsf] -= *uelem * *glloc;

		nupd++;

label1:;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTransp2Index (int _iblk, int _ilist, // Update arrays that support transposed structures search
								const CSMatrixR &_mtrl2) {

//	const char *funcname = "Ilu2UpdateTransp2Index";

	int ilstprv = ibegm[_ilist];

	while (ilstprv != -1) {

		int ilstprv1 = madj[ilstprv];

		if (iv[ilstprv] >= _mtrl2.ia[ilstprv+1]-1) {

			madj[ilstprv] = -1;

		} else {

			int j = iv[ilstprv]+1;
			int icolmn = _mtrl2.ja[j];
			int icolblk = _mtrl2.ja2[j];
			if (icolmn < 0) icolmn = -icolmn-1;
			if (icolblk > _iblk) {
				madj[ilstprv] = -1;
			} else {
				int icolmnl=icolmn;
				madj[ilstprv] = ibegm[icolmnl];
				ibegm[icolmnl] = ilstprv;
			};
		};

		iv[ilstprv]++;

		ilstprv = ilstprv1;

	};

};

// Author: Kharchenko S.A.
// CFctR: Store current row
//========================================================================================
void CFctR::Ilu2StoreRow2Index (int _irow, int _nlevschur) { // Store current row

	const char *funcname = "Ilu2StoreRow2Index";

	double *pg, *pfmask;

// Allocate the memory

	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		jg2loc = new int [nzalloc];
		if (!jg2loc) MemoryFail (funcname);
		jg3loc = new int [nzalloc];
		if (!jg3loc) MemoryFail (funcname);
		adrgloc = new int [0];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [nzalloc];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		jg2blk[iblkalloc]   = jg2loc;
		jg3blk[iblkalloc]   = jg3loc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		jg2loc   = jg2blk[iblkalloc]+nzused;
		jg3loc   = jg3blk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc];
		glloc   = glblk[iblkalloc]+nzused;
		guloc   = gublk[iblkalloc]+nzused;
	};

	jg[_irow] = jgloc;
	jg2[_irow] = jg2loc;
	jg3[_irow] = jg3loc;
	addrg[_irow] = adrgloc;
	gl[_irow] = glloc;
	gu[_irow] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

	int i, j, jblk, jbsf, jbsi, jlev;

	for (i=0;i<nlist;i++) {
		j = lstloc[i];
		jblk = lstloc2[i];
		jbsi = ibsimask[jblk]+j;
		jbsf = ibsfmask[jblk]+j;
		jlev = imasklev[jbsi];
		if (_nlevschur < 0 || jlev <= _nlevschur) {
			jgloc[nzjgloc] = j;
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			pg = glloc+nzloc;
			pfmask = fmaskl+jbsf;
			*pg = *pfmask;
			pg = guloc+nzloc;
			pfmask = fmasku+jbsf;
			*pg = *pfmask;

			nzjgloc++;
			nzloc++;
		};
	};

	ig[_irow+1] = ig[_irow]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store current row
//========================================================================================
void CFctR::Ich2StoreRow2Index (int _irow, int _nlevschur) { // Store current row

	const char *funcname = "Ich2StoreRow2Index";

	double *pg, *pfmask;

// Allocate the memory

	int *jgloc;
	int *jg2loc;
	int *jg3loc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		jg2loc = new int [nzalloc];
		if (!jg2loc) MemoryFail (funcname);
		jg3loc = new int [nzalloc];
		if (!jg3loc) MemoryFail (funcname);
		adrgloc = new int [0];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [0];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		jg2blk[iblkalloc]   = jg2loc;
		jg3blk[iblkalloc]   = jg3loc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		jg2loc   = jg2blk[iblkalloc]+nzused;
		jg3loc   = jg3blk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc];
		glloc   = glblk[iblkalloc]+nzused;
	};

	jg[_irow] = jgloc;
	jg2[_irow] = jg2loc;
	jg3[_irow] = jg3loc;
	addrg[_irow] = adrgloc;
	gl[_irow] = glloc;
	gu[_irow] = glloc;

	int nzjgloc = 0;
	int nzloc = 0;

	int i, j, jblk, jbsf, jbsi, jlev;

	for (i=0;i<nlist;i++) {
		j = lstloc[i];
		jblk = lstloc2[i];
		jbsi = ibsimask[jblk]+j;
		jbsf = ibsfmask[jblk]+j;
		jlev = imasklev[jbsi];
		if (_nlevschur < 0 || jlev <= _nlevschur) {
			jgloc[nzjgloc] = j;
			jg2loc[nzjgloc] = jblk;
			jg3loc[nzjgloc] = jlev;
			pg = glloc+nzloc;
			pfmask = fmaskl+jbsf;
			*pg = *pfmask;

			nzjgloc++;
			nzloc++;
		};
	};

	ig[_irow+1] = ig[_irow]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2StoreU2Index (int _nlistl, int _iblkrow, // Store the resulting U objects
										CSMatrixR &_mtrl, CSMatrixR &_mtru) {

//	const char *funcname = "Ilu2StoreU2Index";

// Count the numbers of elements in the matrices to be stored

	int nlistl;

	nlistl = _nlistl;

	int *jgloc;
	int *jg2loc;
	int *jg3loc;

	int ilist, irow, j, jloc, jj, jjblk, jjloc, jlev;

	int nzjlu = ig[nlistl];
	int nzlu = ig[nlistl];

// Allocate and init matrices

	CSMatrixR mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixR temp (nlistl,nlistl,nzjlu,nzjlu,nzlu);

	temp.AllocateJa3 (nzjlu);

	temp.m = nlistl;
	temp.n = nlistl;
	temp.nsupr = nlistl;
	temp.nsupc = nlistl;

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		jg3loc = jg3[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jlev = jg3loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			temp.ja[nzjlu] = jjloc;
			temp.ja2[nzjlu] = jjblk;
			temp.ja3[nzjlu] = jlev;
			nzjlu++;
			nzlu++;
		};
		temp.list[ilist] = irow;
		temp.list2[ilist] = _iblkrow;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nlist = nlistl;
	temp.nlist2 = nlistl;
	temp.nzja = nzjlu;
	temp.nzja2 = nzjlu;
	temp.nzja3 = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	double *glloc;
	double *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			_mtrl.a[nzlu] = glloc[jloc];
			_mtru.a[nzlu] = guloc[jloc];
			nzjlu++;
			nzlu++;
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ich2StoreU2Index (int _nlistl, int _iblkrow, // Store the resulting U objects
										CSMatrixR &_mtrl) {

//	const char *funcname = "Ich2StoreU2Index";

// Count the numbers of elements in the matrices to be stored

	int nlistl;

	nlistl = _nlistl;

	int *jgloc;
	int *jg2loc;
	int *jg3loc;

	int ilist, irow, j, jloc, jj, jjblk, jjloc, jlev;

	int nzjlu = ig[nlistl];
	int nzlu = ig[nlistl];

// Allocate and init matrices

	CSMatrixR mtrdummy;

	_mtrl = mtrdummy;

	CSMatrixR temp (nlistl,nlistl,nzjlu,nzjlu,nzlu);

	temp.AllocateJa3 (nzjlu);

	temp.m = nlistl;
	temp.n = nlistl;
	temp.nsupr = nlistl;
	temp.nsupc = nlistl;

	temp.ia[0] = 0;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		jg3loc = jg3[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jlev = jg3loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			temp.ja[nzjlu] = jjloc;
			temp.ja2[nzjlu] = jjblk;
			temp.ja3[nzjlu] = jlev;
			nzjlu++;
			nzlu++;
		};
		temp.list[ilist] = irow;
		temp.list2[ilist] = _iblkrow;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nlist = nlistl;
	temp.nlist2 = nlistl;
	temp.nzja = nzjlu;
	temp.nzja2 = nzjlu;
	temp.nzja3 = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	_mtrl = temp;

	double *glloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		irow = ilist;
		jgloc = jg[ilist];
		jg2loc = jg2[ilist];
		glloc = gl[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjblk = jg2loc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			_mtrl.a[nzlu] = glloc[jloc];
			nzjlu++;
			nzlu++;
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row
//========================================================================================
void CFctR::Ilu2UpdateBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
									const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2,
									CSMatrixR &_mtrl, CSMatrixR &_mtru) {

//	const char *funcname = "Ilu2UpdateBlkRow2Index";

// Main cycle over the rows

	int ilist, irow;
	int ilstprv, icolmn, j, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

// Init arrays that support transposed structures search

	Ilu2InitTransp2Index (_iblkrow, nlistl, _mtrl2);

	ig[0] = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

// Init current row

		Ilu2InitRow2Index ('U', _iblkrow, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			j = iv[ilstprv];

			icolmn = _mtrl2.ja[j];

			lelem[0] = _mtrl2.a[j];
			uelem[0] = _mtru2.a[j];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2InitMask2Index (elmisr, 0, ilstprv, _mtrl2, _mtru2);

// Perform updates

			Ilu2UpdateRow2Index ();

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp2Index (_iblkrow, ilist, _mtrl2);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstloc[i];
			jjblk = lstloc2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstloc[i] = ind2arr[i].indy;
			lstloc2[i] = ind2arr[i].indx;
		};

// Store the data in the global arrays

		Ilu2StoreRow2Index (irow, _param.nlevschur);

	};

// Store the resulting U objects

	Ilu2StoreU2Index (nlistl, _iblkrow, 
					_mtrl, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Init arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2InitTransp2Index (int _iblk, int _nlist, // Init arrays that support transposed structures search
												int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr) { 

//	const char *funcname = "Ilu2InitTransp2Index";

	int ilist, irow, j, jj, jjblk, jblkupd, ibs, nlistl, ilistn, ilistupd;

	for (ilistupd=0;ilistupd<_nlistupd;ilistupd++) {

		jblkupd = _listupd[ilistupd];

		ibs = ibsblk[ilistupd];

		for (ilist=0;ilist<_mtrl2arr[jblkupd].nlist;ilist++) {
			irow = _mtrl2arr[jblkupd].list[ilist];
			iv[ibs+ilist] = _mtrl2arr[jblkupd].ia[ilist+1];
			for (j=_mtrl2arr[jblkupd].ia[ilist];j<_mtrl2arr[jblkupd].ia[ilist+1];j++) {
				jj = _mtrl2arr[jblkupd].ja[j];
				jjblk = _mtrl2arr[jblkupd].ja2[j];
				if (jj < 0) jj = -jj-1;
				if (jjblk == _iblk) {
					iv[ibs+ilist] = j;
					break;
				};
			};
		};

	};

	nlistl = _nlist;

	for (ilist=0;ilist<nlistl;ilist++) {
		ibegm[ilist] = -1;
		ibegm2[ilist] = -1;
	};

	int nlistgl = ibsblk[_nlistupd];

	for (ilist=0;ilist<nlistgl;ilist++) {
		madj[ilist] = -1;
		madj2[ilist] = -1;
	};

	for (ilistupd=0;ilistupd<_nlistupd;ilistupd++) {

		jblkupd = _listupd[ilistupd];

		ibs = ibsblk[ilistupd];

		for (ilist=0;ilist<_mtrl2arr[jblkupd].nlist;ilist++) {
			j = iv[ibs+ilist];
			if (j<_mtrl2arr[jblkupd].ia[ilist+1]) {
				jj = _mtrl2arr[jblkupd].ja[j];
				jjblk = _mtrl2arr[jblkupd].ja2[j];
				if (jj < 0) jj = -jj-1;
				if (jjblk == _iblk) {
					ilistn = jj;
					madj[ibs+ilist] = ibegm[ilistn];
					madj2[ibs+ilist] = ibegm2[ilistn];
					ibegm[ilistn] = ilist;
					ibegm2[ilistn] = ilistupd;
				};
			};
		};

	};

};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTransp2Index (int _iblk, int _ilist, // Update arrays that support transposed structures search
												int *_listupd, CSMatrixR *_mtrl2arr) {

//	const char *funcname = "Ilu2UpdateTransp2Index";

	int ilstprv = ibegm[_ilist];
	int iblkupdprv = ibegm2[_ilist];

	while (ilstprv != -1) {

		int jblkupd = _listupd[iblkupdprv];
		int ibs = ibsblk[iblkupdprv];

		int ilstprv1 = madj[ibs+ilstprv];
		int iblkupdprv1 = madj2[ibs+ilstprv];

		if (iv[ibs+ilstprv] >= _mtrl2arr[jblkupd].ia[ilstprv+1]-1) {

			madj[ibs+ilstprv] = -1;
			madj2[ibs+ilstprv] = -1;

		} else {

			int j = iv[ibs+ilstprv]+1;
			int icolmn = _mtrl2arr[jblkupd].ja[j];
			int icolblk = _mtrl2arr[jblkupd].ja2[j];
			if (icolmn < 0) icolmn = -icolmn-1;
			if (icolblk > _iblk) {
				madj[ibs+ilstprv] = -1;
				madj2[ibs+ilstprv] = -1;
			} else {
				int icolmnl=icolmn;
				madj[ibs+ilstprv] = ibegm[icolmnl];
				madj2[ibs+ilstprv] = ibegm2[icolmnl];
				ibegm[icolmnl] = ilstprv;
				ibegm2[icolmnl] = iblkupdprv;
			};
		};

		iv[ibs+ilstprv]++;

		ilstprv = ilstprv1;
		iblkupdprv = iblkupdprv1;

	};

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row by the set of previous block rows
//========================================================================================
void CFctR::Ilu2UpdateBlkRowBySetOfBlocks2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
									int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr, CSMatrixR *_mtru2arr,
									CSMatrixR &_mtrl, CSMatrixR &_mtru) {

	const char *funcname = "Ilu2UpdateBlkRowBySetOfBlocks2Index";

// Main cycle over the rows

	int ilist, irow;
	int ilstprv, icolmn, j, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

// Reallocate the lists sizes and init ibsblk array

	delete [] madj;
	delete [] madj2;
	delete [] ibegm;
	delete [] ibegm2;
	delete [] ibsblk;

	ibsblk = new int [_nlistupd+1];
	if (!ibsblk) MemoryFail (funcname);

	ibsblk[0] = 0;

	int jblk;

	for (ilist=0;ilist<_nlistupd;ilist++) {
		jblk = _listupd[ilist];
		ibsblk[ilist+1] = ibsblk[ilist] + _mtrl2arr[jblk].nlist;
	};

	int nlistgl = ibsblk[_nlistupd];
	if (nlistgl < nlistl) nlistgl = nlistl;

	ibegm = new int [nlistl+1];
	if (!ibegm) MemoryFail (funcname);
	ibegm2 = new int [nlistl+1];
	if (!ibegm2) MemoryFail (funcname);
	madj = new int [nlistgl+1];
	if (!madj) MemoryFail (funcname);
	madj2 = new int [nlistgl+1];
	if (!madj2) MemoryFail (funcname);

// Init arrays that support transposed structures search

	Ilu2InitTransp2Index (_iblkrow, nlistl, _nlistupd, _listupd, _mtrl2arr);

	ig[0] = 0;

	int iblkupdprv, ibs, jblkupd;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

// Init current row

		Ilu2InitRow2Index ('U', _iblkrow, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];
		iblkupdprv = ibegm2[ilist];

		while (ilstprv != -1) {

			jblkupd = _listupd[iblkupdprv];
			ibs = ibsblk[iblkupdprv];

			j = iv[ibs+ilstprv];

			icolmn = _mtrl2arr[jblkupd].ja[j];

			lelem[0] = _mtrl2arr[jblkupd].a[j];
			uelem[0] = _mtru2arr[jblkupd].a[j];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2InitMask2Index (elmisr, ibs, ilstprv, _mtrl2arr[jblkupd], _mtru2arr[jblkupd]);

// Perform updates

			Ilu2UpdateRow2Index ();

			iblkupdprv = madj2[ibs+ilstprv];
			ilstprv = madj[ibs+ilstprv];

		};

// Update arrays iv, madj, madj2, ibegm and ibegm2 for the next column

		Ilu2UpdateTransp2Index (_iblkrow, ilist, _listupd, _mtrl2arr);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstloc[i];
			jjblk = lstloc2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstloc[i] = ind2arr[i].indy;
			lstloc2[i] = ind2arr[i].indx;
		};

// Store the data in the global arrays

		Ilu2StoreRow2Index (irow, _param.nlevfct);

	};

// Store the resulting U objects

	Ilu2StoreU2Index (nlistl, _iblkrow, 
					_mtrl, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask, update arrays and perform computations
//========================================================================================
void CFctR::Ilu2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, int _irowprv, // Init mask, update arrays and perform computations
											const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2) {

//	const char *funcname = "Ilu2InitMask2Index";

	double dzero = 0.0e0;

	int jbsi, jbsf;

	int k;
	int jcolmn, jcolblk, ilev, jlev, jlevnew;
	double *glloc, *guloc;

	nupd = 0;

	k = ivgl[_ibs+_irowprv];

	ilev = _mtrl2.ja3[k];

	for (k=ivgl[_ibs+_irowprv];k<_mtrl2.ia[_ilstprv+1];k++) {

		jcolmn = _mtrl2.ja[k];
		jcolblk = _mtrl2.ja2[k];
		jlev = _mtrl2.ja3[k];

		jlevnew = ilev + jlev;

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn-1;

		jbsi = ibsimask[jcolblk]+jcolmn;
		jbsf = ibsfmask[jcolblk]+jcolmn;

		if (imask[jbsi] != icycle) {
			imask[jbsi] = icycle;
			imasklev[jbsi] = jlevnew;
			lstloc[nlist] = jcolmn;
			lstloc2[nlist] = jcolblk;
			fmaskl[jbsf] = dzero;
			fmasku[jbsf] = dzero;
			nlist++;
		} else {
			if (jlevnew < imasklev[jbsi]) imasklev[jbsi] = jlevnew;
		};

		lstupd[nupd] = jcolmn;
		lstupd2[nupd] = jcolblk;

//		adrupdl[nupd] = _mtrl2.a+k;
//		adrupdu[nupd] = _mtru2.a+k;

		glloc = _mtrl2.a+k;
		guloc = _mtru2.a+k;

		fmasku[jbsf] -= *lelem * *guloc;
		fmaskl[jbsf] -= *uelem * *glloc;

		nupd++;

label1:;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Init mask, update arrays and perform computations
//========================================================================================
void CFctR::Ich2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, int _irowprv, // Init mask, update arrays and perform computations
											const CSMatrixR &_mtrl2) {

//	const char *funcname = "Ich2InitMask2Index";

	double dzero = 0.0e0;

	int jbsi, jbsf;

	int k;
	int jcolmn, jcolblk, ilev, jlev, jlevnew;
	double *glloc;

	nupd = 0;

	k = ivgl[_ibs+_irowprv];

	ilev = _mtrl2.ja3[k];

	for (k=ivgl[_ibs+_irowprv];k<_mtrl2.ia[_ilstprv+1];k++) {

		jcolmn = _mtrl2.ja[k];
		jcolblk = _mtrl2.ja2[k];
		jlev = _mtrl2.ja3[k];

		jlevnew = ilev + jlev;

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn-1;

		jbsi = ibsimask[jcolblk]+jcolmn;
		jbsf = ibsfmask[jcolblk]+jcolmn;

		if (imask[jbsi] != icycle) {
			imask[jbsi] = icycle;
			imasklev[jbsi] = jlevnew;
			lstloc[nlist] = jcolmn;
			lstloc2[nlist] = jcolblk;
			fmaskl[jbsf] = dzero;
//			fmasku[jbsf] = dzero;
			nlist++;
		} else {
			if (jlevnew < imasklev[jbsi]) imasklev[jbsi] = jlevnew;
		};

		lstupd[nupd] = jcolmn;
		lstupd2[nupd] = jcolblk;

//		adrupdl[nupd] = _mtrl2.a+k;
//		adrupdu[nupd] = _mtru2.a+k;

		glloc = _mtrl2.a+k;

//		fmasku[jbsf] -= *lelem * *guloc;
		fmaskl[jbsf] -= *uelem * *glloc;

		nupd++;

label1:;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTranspGlobal2Index (int _iblk, int _irow, // Update arrays that support transposed structures search
														int *_listupd, CSMatrixR *_mtrl2arr) {

//	const char *funcname = "Ilu2UpdateTranspGlobal2Index";

	int ibs = ibsblkgl[_iblk];
	int iirow = ibs+_irow;

	int irowprv = ibegmgl[iirow];
	int iblkprv = ibegm2gl[iirow];

	while (irowprv != -1) {

		int jbs = ibsblkgl[iblkprv];
		int jjrow = jbs+irowprv;

		int irowprv1 = madjgl[jjrow];
		int iblkprv1 = madj2gl[jjrow];

		int ilstprv = irw2lst[jjrow];

		if (ivgl[jjrow] >= _mtrl2arr[iblkprv].ia[ilstprv+1]-1) {

			madjgl[jjrow] = -1;
			madj2gl[jjrow] = -1;

		} else {

			int j = ivgl[jjrow]+1;
			int icolmn = _mtrl2arr[iblkprv].ja[j];
			int icolblk = _mtrl2arr[iblkprv].ja2[j];
			if (icolmn < 0) icolmn = -icolmn-1;

			int ibs1 = ibsblkgl[icolblk];
			int iirow1 = ibs1+icolmn;

			madjgl[jjrow] = ibegmgl[iirow1];
			madj2gl[jjrow] = ibegm2gl[iirow1];
			ibegmgl[iirow1] = irowprv;
			ibegm2gl[iirow1] = iblkprv;

		};

		ivgl[jjrow]++;

		irowprv = irowprv1;
		iblkprv = iblkprv1;

	};

	ibegmgl[iirow] = -1;
	ibegm2gl[iirow] = -1;

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row by the set of previous block rows with the global transpose search
//========================================================================================
void CFctR::Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows with the global transpose search
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
									int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr, CSMatrixR *_mtru2arr,
									CSMatrixR &_mtrl, CSMatrixR &_mtru) {

//	const char *funcname = "Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index";

// Main cycle over the rows

	int ilist, irow;
	int ilstprv, j, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

// Main cycle

	ig[0] = 0;

	int iblkprv, ibs, iirow, irowprv, jcol, jblk, jbs, jjrow, jlist;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];
		ibs = ibsblkgl[_iblkrow];
		iirow = ibs + irow;

// Init current row

		Ilu2InitRow2Index ('U', _iblkrow, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		irowprv = ibegmgl[iirow];
		iblkprv = ibegm2gl[iirow];

		while (irowprv != -1) {

			jbs = ibsblkgl[iblkprv];
			jjrow = jbs+irowprv;

			j = ivgl[jjrow];

			jlist = irw2lst[jjrow];

			if (j >= _mtrl2arr[iblkprv].ia[jlist+1]) {
				throw " CFctR::Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index: Wrong array ivgl ";
			};

			jcol = _mtrl2arr[iblkprv].ja[j];
			jblk = _mtrl2arr[iblkprv].ja2[j];

			elmisr = jcol < 0;

			if (jcol < 0) jcol = -jcol-1;

			if (jcol != irow || jblk != _iblkrow) {
				throw " CFctR::Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index: Wrong array ivgl ";
			};

			lelem[0] = _mtrl2arr[iblkprv].a[j];
			uelem[0] = _mtru2arr[iblkprv].a[j];

// Init and register if necessary the mask and perform updates

			ilstprv = irw2lst[jjrow];

			Ilu2InitMask2Index (elmisr, jbs, ilstprv, irowprv, _mtrl2arr[iblkprv], _mtru2arr[iblkprv]);

			irowprv = madjgl[jjrow];
			iblkprv = madj2gl[jjrow];

		};

// Update arrays iv, madj, madj2, ibegm and ibegm2 for the next column

		Ilu2UpdateTranspGlobal2Index (_iblkrow, irow, _listupd, _mtrl2arr);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstloc[i];
			jjblk = lstloc2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstloc[i] = ind2arr[i].indy;
			lstloc2[i] = ind2arr[i].indx;
		};

// Store the data in the global arrays

		Ilu2StoreRow2Index (irow, _param.nlevschur);

	};

// Store the resulting U objects

	Ilu2StoreU2Index (nlistl, _iblkrow, 
					_mtrl, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row by the set of previous block rows with the global transpose search
//========================================================================================
void CFctR::Ich2UpdateBlkRowBySetOfBlocksGlobal2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows with the global transpose search
									const CGSMatrix &_gmtra, const CSMatrixR &_mtral, 
									int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr,
									CSMatrixR &_mtrl) {

//	const char *funcname = "Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index";

// Main cycle over the rows

	int ilist, irow;
	int ilstprv, j, i, jj, jjblk;
	int indbeg, indend, nlistl;
	bool elmisr;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	indbeg = _gmtra.blksc[_iblkrow];
	indend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = indend-indbeg+1;

// Main cycle

	ig[0] = 0;

	int iblkprv, ibs, iirow, irowprv, jcol, jblk, jbs, jjrow, jlist;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];
		ibs = ibsblkgl[_iblkrow];
		iirow = ibs + irow;

// Init current row

		Ich2InitRow2Index ('U', _iblkrow, irow, _gmtra, _mtral);

// Update current row by the previous ones

		irowprv = ibegmgl[iirow];
		iblkprv = ibegm2gl[iirow];

		while (irowprv != -1) {

			jbs = ibsblkgl[iblkprv];
			jjrow = jbs+irowprv;

			j = ivgl[jjrow];

			jlist = irw2lst[jjrow];

			if (j >= _mtrl2arr[iblkprv].ia[jlist+1]) {
				throw " CFctR::Ich2UpdateBlkRowBySetOfBlocksGlobal2Index: Wrong array ivgl ";
			};

			jcol = _mtrl2arr[iblkprv].ja[j];
			jblk = _mtrl2arr[iblkprv].ja2[j];

			elmisr = jcol < 0;

			if (jcol < 0) jcol = -jcol-1;

			if (jcol != irow || jblk != _iblkrow) {
				throw " CFctR::Ich2UpdateBlkRowBySetOfBlocksGlobal2Index: Wrong array ivgl ";
			};

			lelem[0] = _mtrl2arr[iblkprv].a[j];
			uelem[0] = lelem[0];

// Init and register if necessary the mask and perform updates

			ilstprv = irw2lst[jjrow];

//			Ilu2InitMask2Index (elmisr, jbs, ilstprv, irowprv, _mtrl2arr[iblkprv], _mtru2arr[iblkprv]);
			Ich2InitMask2Index (elmisr, jbs, ilstprv, irowprv, _mtrl2arr[iblkprv]);

			irowprv = madjgl[jjrow];
			iblkprv = madj2gl[jjrow];

		};

// Update arrays iv, madj, madj2, ibegm and ibegm2 for the next column

		Ilu2UpdateTranspGlobal2Index (_iblkrow, irow, _listupd, _mtrl2arr);

// Sort the list of indices

		for (i=0;i<nlist;i++) {
			jj = lstloc[i];
			jjblk = lstloc2[i];
			ind2arr[i].indx = jjblk;
			ind2arr[i].indy = jj;
			ind2arr[i].intvalue = i;
		};

		if (nlist != 0) qsort (ind2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

		for (i=0;i<nlist;i++) {
			lstloc[i] = ind2arr[i].indy;
			lstloc2[i] = ind2arr[i].indx;
		};

// Store the data in the global arrays

		Ich2StoreRow2Index (irow, _param.nlevschur);

	};

// Store the resulting U objects

	Ich2StoreU2Index (nlistl, _iblkrow, 
					_mtrl);

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] jg2blk[ilist];
		delete [] jg3blk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Add two block rows/columns assuming that the list of indices coincide
//========================================================================================
CSMatrixR CFctR::AddBlocks2Index ( // Add two block rows/columns assuming that the list of indices coincide
									CSMatrixR &_mtru, CSMatrixR &_mtru2) {

	const char *funcname = "AddBlocks2Index";

//	double dzero = 0.0e0;

// Check that the lists of indices coincide

	if (_mtru.nlist != _mtru2.nlist) {
		assert (false); throw " CFctR::AddBlocks2Index: The lists of indices do not coincide";
	};
	if (_mtru.nlist2 != _mtru2.nlist2) {
		assert (false); throw " CFctR::AddBlocks2Index: The lists of indices do not coincide";
	};
	if (_mtru.nlist != _mtru.nlist2) {
		assert (false); throw " CFctR::AddBlocks2Index: The lists of indices do not coincide";
	};
	int i;
	for (i=0;i<_mtru.nlist;i++) {
		if (_mtru.list[i] != _mtru2.list[i]) {
			assert (false); throw " CFctR::AddBlocks2Index: The lists of indices do not coincide";
		};
		if (_mtru.list2[i] != _mtru2.list2[i]) {
			assert (false); throw " CFctR::AddBlocks2Index: The lists of indices do not coincide";
		};
	};

// Count the total number of nonzero supernodes

	int *ialoc;

	ialoc = new int [_mtru.nlist+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	int nzjaloc = 0;

	int ilist, j, jj, jjblk, jbsi, jjloc;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		icyclegl++;
		for (j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {
			jj = _mtru.ja[j];
			jjblk = _mtru.ja2[j];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			jbsi = ibsblkgl[jjblk]+jjloc;
			if (imaskgl[jbsi] != icyclegl) {
				imaskgl[jbsi] = icyclegl;
				nzjaloc++;
			};
		};
		for (j=_mtru2.ia[ilist];j<_mtru2.ia[ilist+1];j++) {
			jj = _mtru2.ja[j];
			jjblk = _mtru2.ja2[j];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			jbsi = ibsblkgl[jjblk]+jjloc;
			if (imaskgl[jbsi] != icyclegl) {
				imaskgl[jbsi] = icyclegl;
				nzjaloc++;
			};
		};
		ialoc[ilist+1] = nzjaloc;
	};

// Allocate and init array jaloc

	int *jaloc, *jaloc2;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);
	jaloc2 = new int [nzjaloc];
	if (!jaloc2) MemoryFail (funcname);

	nzjaloc = 0;

	int iloc;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		icyclegl++;
		for (j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {
			jj = _mtru.ja[j];
			jjblk = _mtru.ja2[j];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			jbsi = ibsblkgl[jjblk]+jjloc;
			if (imaskgl[jbsi] != icyclegl) {
				imaskgl[jbsi] = icyclegl;
				jaloc[nzjaloc] = jjloc;
				jaloc2[nzjaloc] = jjblk;
				nzjaloc++;
			};
		};
		for (j=_mtru2.ia[ilist];j<_mtru2.ia[ilist+1];j++) {
			jj = _mtru2.ja[j];
			jjblk = _mtru2.ja2[j];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc-1;
			jbsi = ibsblkgl[jjblk]+jjloc;
			if (imaskgl[jbsi] != icyclegl) {
				imaskgl[jbsi] = icyclegl;
				jaloc[nzjaloc] = jjloc;
				jaloc2[nzjaloc] = jjblk;
				nzjaloc++;
			};
		};
	};

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		int nz = ialoc[ilist+1]-ialoc[ilist];
		int ibs = ialoc[ilist];
		for (i=ialoc[ilist];i<ialoc[ilist+1];i++) {
			iloc = i-ibs;
			jj = jaloc[i];
			jjblk = jaloc2[i];
			ind2arrgl[iloc].indx = jjblk;
			ind2arrgl[iloc].indy = jj;
			ind2arrgl[iloc].intvalue = iloc;
		};
		if (nz != 0) qsort (ind2arrgl, nz, sizeof(CInd2Int), CInd2Int::CompareInd2Int);
		for (i=ialoc[ilist];i<ialoc[ilist+1];i++) {
			iloc = i-ibs;
			jaloc[i] = ind2arrgl[iloc].indy;
			jaloc2[i] = ind2arrgl[iloc].indx;
		};
	};

// Store computed data

	CSMatrixR temp (_mtru.nlist,_mtru.nlist,nzjaloc,nzjaloc,nzjaloc);

	temp.AllocateJa3 (nzjaloc);

	for (i=0;i< _mtru.nlist;i++) temp.list[i] = _mtru.list[i];
	for (i=0;i< _mtru.nlist;i++) temp.list2[i] = _mtru.list2[i];
	for (i=0;i<=_mtru.nlist;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];
	for (i=0;i<nzjaloc;i++) temp.ja2[i] = jaloc2[i];

// Free work memory

	delete [] ialoc;
	delete [] jaloc;
	delete [] jaloc2;

// Compute array a

	int ibeg1, iend1, ibeg2, iend2, ind1, ind2;
	int jsup1, jblk1, jlev1, jsup2, jblk2, jlev2, jlev;

	nzjaloc = 0;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		ibeg1 = _mtru.ia[ilist];
		iend1 = _mtru.ia[ilist+1]-1;
		ibeg2 = _mtru2.ia[ilist];
		iend2 = _mtru2.ia[ilist+1]-1;
		ind1 = ibeg1;
		ind2 = ibeg2;
		while (ind1 <= iend1 || ind2 <= iend2) {
			if (ind1 <= iend1 && ind2 <= iend2) {
				jsup1 = _mtru.ja[ind1];
				if (jsup1 < 0) jsup1 = -jsup1-1;
				jblk1 = _mtru.ja2[ind1];
				jlev1 = _mtru.ja3[ind1];
				jsup2 = _mtru2.ja[ind2];
				if (jsup2 < 0) jsup2 = -jsup2-1;
				jblk2 = _mtru2.ja2[ind2];
				jlev2 = _mtru2.ja3[ind2];
				if (jsup1 == jsup2 && jblk1 == jblk2) {
					jlev = jlev1;
					if (jlev > jlev2) jlev = jlev2;
					temp.ja3[nzjaloc] = jlev;
					temp.a[nzjaloc] = _mtru.a[ind1]+_mtru2.a[ind2];
					nzjaloc++;
					ind1++;
					ind2++;
				} else if (jblk1 < jblk2 || (jblk1 == jblk2 && jsup1 < jsup2)) {
					temp.ja3[nzjaloc] = jlev1;
					temp.a[nzjaloc] = _mtru.a[ind1];
					nzjaloc++;
					ind1++;
				} else if (jblk1 > jblk2 || (jblk1 == jblk2 && jsup1 > jsup2)) {
					temp.ja3[nzjaloc] = jlev2;
					temp.a[nzjaloc] = _mtru2.a[ind2];
					nzjaloc++;
					ind2++;
				};
			} else if (ind1 <= iend1) {
				temp.ja3[nzjaloc] = _mtru.ja3[ind1];
				temp.a[nzjaloc] = _mtru.a[ind1];
				nzjaloc++;
				ind1++;
			} else if (ind2 <= iend2) {
				temp.ja3[nzjaloc] = _mtru2.ja3[ind2];
				temp.a[nzjaloc] = _mtru2.a[ind2];
				nzjaloc++;
				ind2++;
			};
		};
	};

// Fill control data

	temp.m      = _mtru.m;
	temp.n      = _mtru.n;
	temp.nsupr  = _mtru.m;
	temp.nsupc  = _mtru.m;
	temp.nlist  = _mtru.nlist;
	temp.nlist2 = _mtru.nlist;
	temp.nzja   = nzjaloc;
	temp.nzja2  = nzjaloc;
	temp.nzja3  = nzjaloc;
	temp.nza    = nzjaloc;
	temp.nzatot = nzjaloc;

	return temp;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Add and replace set of block rows/columns
//========================================================================================
void CGSMatrixR::AddReplaceSetOfBlocks2Index (CFctR &_fct, int _nlistblk, int *_listblk, // Add and replace set of block rows/columns
																CSMatrixR *_mtrlarr, CSMatrixR *_mtrl2arr, CSMatrixR *_mtruarr, CSMatrixR *_mtru2arr) {

	const char *funcname = "AddReplaceSetOfBlocks2Index";

//	double dzero = 0.0e0;

// Create the list of blocks involved

	int nblksloc = nblksr;

	int *imaskblk;
	int *listblk;

	imaskblk = new int [nblksloc];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [nblksloc];
	if (!listblk) MemoryFail (funcname);

	int i, iblk, iiblk, ilistblk;

	for (i=0;i<nblksloc;i++) imaskblk[i] = -1;

	int icycle = -1;

	icycle++;

	int nlistblk = 0;

	CSMatrixR *pmtrlblk, *pmtrl2blk, *pmtrublk, *pmtru2blk;

	for (ilistblk=0;ilistblk<_nlistblk;ilistblk++) {
		iiblk = _listblk[ilistblk];
		pmtrublk = _mtruarr+iiblk;
		pmtru2blk = _mtru2arr+iiblk;
		for (i=0;i<pmtrublk->nlist2;i++) {
			iblk = pmtrublk->list2[i];
			if (imaskblk[iblk] != icycle) {
				listblk[nlistblk] = iblk;
				nlistblk++;
				imaskblk[iblk] = icycle;
			};
		};
		for (i=0;i<pmtrublk->nzja2;i++) {
			iblk = pmtrublk->ja2[i];
			if (imaskblk[iblk] != icycle) {
				listblk[nlistblk] = iblk;
				nlistblk++;
				imaskblk[iblk] = icycle;
			};
		};
		for (i=0;i<pmtru2blk->nlist2;i++) {
			iblk = pmtru2blk->list2[i];
			if (imaskblk[iblk] != icycle) {
				listblk[nlistblk] = iblk;
				nlistblk++;
				imaskblk[iblk] = icycle;
			};
		};
		for (i=0;i<pmtru2blk->nzja2;i++) {
			iblk = pmtru2blk->ja2[i];
			if (imaskblk[iblk] != icycle) {
				listblk[nlistblk] = iblk;
				nlistblk++;
				imaskblk[iblk] = icycle;
			};
		};
	};

	qsort (listblk, nlistblk, sizeof(int), compint);

// Main cycle over the blocks

	int jj, ni;

	for (ilistblk=0;ilistblk<_nlistblk;ilistblk++) {

		iiblk = _listblk[ilistblk];
		ni = blksr[iiblk+1] - blksr[iiblk];

		pmtrlblk = _mtrlarr+iiblk;
		pmtrl2blk = _mtrl2arr+iiblk;
		pmtrublk = _mtruarr+iiblk;
		pmtru2blk = _mtru2arr+iiblk;

		if (pmtrublk->nlist != ni) {
			assert (false); throw " CGSMatrixR::AddReplaceSetOfBlocks2Index: Not all rows in the list";
		};

// Check inclusion of the lists

		_fct.icyclegl++;

		int ibs = _fct.ibsblkgl[iiblk];

		int *imaskia;

		imaskia = new int [pmtrublk->nlist];
		if (!imaskia) MemoryFail (funcname);

		for (i=0;i<pmtrublk->nlist;i++) imaskia[i] = -1;

		for (i=0;i<pmtrublk->nlist;i++) {
			jj = pmtrublk->list[i];
			_fct.imaskgl[ibs+jj] = _fct.icyclegl;
		};
		for (i=0;i<pmtru2blk->nlist;i++) {
			jj = pmtru2blk->list[i];
			if (_fct.imaskgl[ibs+jj] != _fct.icyclegl) {
				assert (false); throw " CGSMatrixR::AddReplaceSetOfBlocks2Index: Assumed inclusion of lists is not correct";
			};
			imaskia[jj] = i;
		};

// Count the total number of nonzero supernodes

		int *ialoc;

		ialoc = new int [pmtrublk->nlist+1];
		if (!ialoc) MemoryFail (funcname);

		ialoc[0] = 0;

		int nzjaloc = 0;

		int ilist, ilist1, j, jj, jjblk, jbsi, jjloc;

		for (ilist=0;ilist<pmtrublk->nlist;ilist++) {
			_fct.icyclegl++;
			for (j=pmtrublk->ia[ilist];j<pmtrublk->ia[ilist+1];j++) {
				jj = pmtrublk->ja[j];
				jjblk = pmtrublk->ja2[j];
				jjloc = jj;
				if (jjloc < 0) jjloc = -jjloc-1;
				jbsi = _fct.ibsblkgl[jjblk]+jjloc;
				if (_fct.imaskgl[jbsi] != _fct.icyclegl) {
					_fct.imaskgl[jbsi] = _fct.icyclegl;
					nzjaloc++;
				};
			};
			if (imaskia[ilist] >= 0) {
				ilist1 = imaskia[ilist];
				for (j=pmtru2blk->ia[ilist1];j<pmtru2blk->ia[ilist1+1];j++) {
					jj = pmtru2blk->ja[j];
					jjblk = pmtru2blk->ja2[j];
					jjloc = jj;
					if (jjloc < 0) jjloc = -jjloc-1;
					jbsi = _fct.ibsblkgl[jjblk]+jjloc;
					if (_fct.imaskgl[jbsi] != _fct.icyclegl) {
						_fct.imaskgl[jbsi] = _fct.icyclegl;
						nzjaloc++;
					};
				};
			};
			ialoc[ilist+1] = nzjaloc;
		};

// Allocate and init array jaloc

		int *jaloc, *jaloc2;

		jaloc = new int [nzjaloc];
		if (!jaloc) MemoryFail (funcname);
		jaloc2 = new int [nzjaloc];
		if (!jaloc2) MemoryFail (funcname);

		nzjaloc = 0;

		int iloc;

		for (ilist=0;ilist<pmtrublk->nlist;ilist++) {
			_fct.icyclegl++;
			for (j=pmtrublk->ia[ilist];j<pmtrublk->ia[ilist+1];j++) {
				jj = pmtrublk->ja[j];
				jjblk = pmtrublk->ja2[j];
				jjloc = jj;
				if (jjloc < 0) jjloc = -jjloc-1;
				jbsi = _fct.ibsblkgl[jjblk]+jjloc;
				if (_fct.imaskgl[jbsi] != _fct.icyclegl) {
					_fct.imaskgl[jbsi] = _fct.icyclegl;
					jaloc[nzjaloc] = jjloc;
					jaloc2[nzjaloc] = jjblk;
					nzjaloc++;
				};
			};
			if (imaskia[ilist] >= 0) {
				ilist1 = imaskia[ilist];
				for (j=pmtru2blk->ia[ilist1];j<pmtru2blk->ia[ilist1+1];j++) {
					jj = pmtru2blk->ja[j];
					jjblk = pmtru2blk->ja2[j];
					jjloc = jj;
					if (jjloc < 0) jjloc = -jjloc-1;
					jbsi = _fct.ibsblkgl[jjblk]+jjloc;
					if (_fct.imaskgl[jbsi] != _fct.icyclegl) {
						_fct.imaskgl[jbsi] = _fct.icyclegl;
						jaloc[nzjaloc] = jjloc;
						jaloc2[nzjaloc] = jjblk;
						nzjaloc++;
					};
				};
			};
		};

		for (ilist=0;ilist<_mtruarr[iiblk].nlist;ilist++) {
			int nz = ialoc[ilist+1]-ialoc[ilist];
			int ibs = ialoc[ilist];
			for (i=ialoc[ilist];i<ialoc[ilist+1];i++) {
				iloc = i-ibs;
				jj = jaloc[i];
				jjblk = jaloc2[i];
				_fct.ind2arrgl[iloc].indx = jjblk;
				_fct.ind2arrgl[iloc].indy = jj;
				_fct.ind2arrgl[iloc].intvalue = iloc;
			};
			if (nz != 0) qsort (_fct.ind2arrgl, nz, sizeof(CInd2Int), CInd2Int::CompareInd2Int);
			for (i=ialoc[ilist];i<ialoc[ilist+1];i++) {
				iloc = i-ibs;
				jaloc[i] = _fct.ind2arrgl[iloc].indy;
				jaloc2[i] = _fct.ind2arrgl[iloc].indx;
			};
		};

// Store computed data

		CSMatrixR temp (pmtrublk->nlist,pmtrublk->nlist,nzjaloc,nzjaloc,nzjaloc);

		temp.AllocateJa3 (nzjaloc);

		for (i=0;i< pmtrublk->nlist;i++) temp.list[i] = pmtrublk->list[i];
		for (i=0;i< pmtrublk->nlist;i++) temp.list2[i] = pmtrublk->list2[i];
		for (i=0;i<=pmtrublk->nlist;i++) temp.ia[i] = ialoc[i];
		for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];
		for (i=0;i<nzjaloc;i++) temp.ja2[i] = jaloc2[i];

// Fill control data

		temp.m      = _mtruarr[iiblk].m;
		temp.n      = _mtruarr[iiblk].n;
		temp.nsupr  = _mtruarr[iiblk].m;
		temp.nsupc  = _mtruarr[iiblk].m;
		temp.nlist  = _mtruarr[iiblk].nlist;
		temp.nlist2 = _mtruarr[iiblk].nlist;
		temp.nzja   = nzjaloc;
		temp.nzja2  = nzjaloc;
		temp.nzja3  = nzjaloc;
		temp.nza    = nzjaloc;
		temp.nzatot = nzjaloc;

// Free work memory

		delete [] ialoc;
		delete [] jaloc;
		delete [] jaloc2;

// Compute array a

		double *alloc;

		alloc = new double [nzjaloc];
		if (!alloc) MemoryFail (funcname);

		int ibeg1, iend1, ibeg2, iend2, ind1, ind2;
		int jsup1, jblk1, jlev1, jsup2, jblk2, jlev2, jlev;

		nzjaloc = 0;

		for (ilist=0;ilist<pmtrublk->nlist;ilist++) {
			ibeg1 = pmtrublk->ia[ilist];
			iend1 = pmtrublk->ia[ilist+1]-1;
			ilist1 = imaskia[ilist];
			if (ilist1 >= 0) {
				ibeg2 = _mtru2arr[iiblk].ia[ilist1];
				iend2 = _mtru2arr[iiblk].ia[ilist1+1]-1;
			} else {
				ibeg2 = 1;
				iend2 = 0;
			};
			ind1 = ibeg1;
			ind2 = ibeg2;
			while (ind1 <= iend1 || ind2 <= iend2) {
				if (ind1 <= iend1 && ind2 <= iend2) {
					jsup1 = pmtrublk->ja[ind1];
					if (jsup1 < 0) jsup1 = -jsup1-1;
					jblk1 = pmtrublk->ja2[ind1];
					jlev1 = pmtrublk->ja3[ind1];
					jsup2 = pmtru2blk->ja[ind2];
					if (jsup2 < 0) jsup2 = -jsup2-1;
					jblk2 = pmtru2blk->ja2[ind2];
					jlev2 = pmtru2blk->ja3[ind2];
					if (jsup1 == jsup2 && jblk1 == jblk2) {
						jlev = jlev1;
						if (jlev > jlev2) jlev = jlev2;
						temp.ja3[nzjaloc] = jlev;
						temp.a[nzjaloc] = pmtrublk->a[ind1]+pmtru2blk->a[ind2];
						alloc[nzjaloc] = pmtrlblk->a[ind1]+pmtrl2blk->a[ind2];
						nzjaloc++;
						ind1++;
						ind2++;
					} else if (jblk1 < jblk2 || (jblk1 == jblk2 && jsup1 < jsup2)) {
						temp.ja3[nzjaloc] = jlev1;
						temp.a[nzjaloc] = pmtrublk->a[ind1];
						alloc[nzjaloc] = pmtrlblk->a[ind1];
						nzjaloc++;
						ind1++;
					} else if (jblk1 > jblk2 || (jblk1 == jblk2 && jsup1 > jsup2)) {
						temp.ja3[nzjaloc] = jlev2;
						temp.a[nzjaloc] = pmtru2blk->a[ind2];
						alloc[nzjaloc] = pmtrl2blk->a[ind2];
						nzjaloc++;
						ind2++;
					};
				} else if (ind1 <= iend1) {
					temp.ja3[nzjaloc] = pmtrublk->ja3[ind1];
					temp.a[nzjaloc] = pmtrublk->a[ind1];
					alloc[nzjaloc] = pmtrlblk->a[ind1];
					nzjaloc++;
					ind1++;
				} else if (ind2 <= iend2) {
					temp.ja3[nzjaloc] = pmtru2blk->ja3[ind2];
					temp.a[nzjaloc] = pmtru2blk->a[ind2];
					alloc[nzjaloc] = pmtrl2blk->a[ind2];
					nzjaloc++;
					ind2++;
				};
			};
		};

// Free mask array

		delete [] imaskia;

		_mtruarr[iiblk] = temp;
		_mtrlarr[iiblk] = temp;

		for (i=0;i<nzjaloc;i++) _mtrlarr[iiblk].a[i] = alloc[i];

		delete [] alloc;

	};

// Free work arrays

	delete [] imaskblk;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Filter list of Schur blocks according to the block indices
//========================================================================================
void CGSMatrixR::FilterListOfSchurBlocks2Index (int _nlistblk, int *_listblk, // Filter list of Schur blocks according to the block indices
																int *_iabextflt, int *_jabextflt, double _theta,
																CFctR &_fct) {

	const char *funcname = "FilterListOfSchurBlocks2Index";

// Create block mask array

	int nblksloc = nblksr;

	int *imaskblk;

	imaskblk = new int [nblksloc];
	if (!imaskblk) MemoryFail (funcname);

	int icycle = -1;

	int i;

	for (i=0;i<nblksloc;i++) imaskblk[i] = icycle;

// Main cycle over the blocks

	int ilistblk, iiblk, j, jjblk;

	for (ilistblk=0;ilistblk<_nlistblk;ilistblk++) {

		iiblk = _listblk[ilistblk];

// Modify block mask

		icycle++;

		for (j=_iabextflt[iiblk];j<_iabextflt[iiblk+1];j++) {
			jjblk = _jabextflt[j];
			imaskblk[jjblk] = icycle;
		};

// Create new list of supernode rows

		int *listcount;

		int nlistloc = mtrarr[iiblk].nlist;

		listcount = new int [nlistloc];
		if (!listcount) MemoryFail (funcname);

		int ilist, nz;

		for (ilist=0;ilist<nlistloc;ilist++) {
			nz = 0;
			for (j=mtrarr[iiblk].ia[ilist];j<mtrarr[iiblk].ia[ilist+1];j++) {
				jjblk = mtrarr[iiblk].ja2[j];
				if (imaskblk[jjblk] == icycle) nz++;
			};
			listcount[ilist] = nz;
		};

		int nlistnew = 0;

		for (ilist=0;ilist<nlistloc;ilist++) {
			if (listcount[ilist]>0) nlistnew++;
		};

		int *listnew;
		int *list2new;

		listnew = new int [nlistnew];
		if (!listnew) MemoryFail (funcname);
		list2new = new int [nlistnew];
		if (!list2new) MemoryFail (funcname);

		nlistnew = 0;

		for (ilist=0;ilist<nlistloc;ilist++) {
			if (listcount[ilist]>0) {
				listnew[nlistnew] = mtrarr[iiblk].list[ilist];
				list2new[nlistnew] = mtrarr[iiblk].list2[ilist];
				nlistnew++;
			};
		};

// Create ia, ja, ja2 and ja3 arrays

		int *ialoc;

		ialoc = new int [nlistnew+1];
		if (!ialoc) MemoryFail (funcname);

		ialoc[0] = 0;

		nlistnew = 0;

		for (ilist=0;ilist<nlistloc;ilist++) {
			if (listcount[ilist]>0) {
				ialoc[nlistnew+1] = ialoc[nlistnew] + listcount[ilist];
				nlistnew++;
			};
		};

		int nzjanew = ialoc[nlistnew];

		int *jaloc;
		int *ja2loc;
		int *ja3loc;

		jaloc = new int [nzjanew];
		if (!jaloc) MemoryFail (funcname);
		ja2loc = new int [nzjanew];
		if (!ja2loc) MemoryFail (funcname);
		ja3loc = new int [nzjanew];
		if (!ja3loc) MemoryFail (funcname);

		int nzjaloc = 0;

		for (ilist=0;ilist<nlistloc;ilist++) {
			for (j=mtrarr[iiblk].ia[ilist];j<mtrarr[iiblk].ia[ilist+1];j++) {
				jjblk = mtrarr[iiblk].ja2[j];
				if (imaskblk[jjblk] == icycle) {
					jaloc[nzjaloc] = mtrarr[iiblk].ja[j];
					ja2loc[nzjaloc] = jjblk;
					ja3loc[nzjaloc] = mtrarr[iiblk].ja3[j];
					nzjaloc++;
				};
			};
		};

// Create a and modify the corresponding diagg array

		double *aloc;

		aloc = new double [nzjaloc];
		if (!aloc) MemoryFail (funcname);

		nzjaloc = 0;

		int irow, iblk, jcol, jjblk, ibsd, jbsd;
		double auxloc;

		for (ilist=0;ilist<nlistloc;ilist++) {
			irow = mtrarr[iiblk].list[ilist];
			iblk = mtrarr[iiblk].list2[ilist];
			for (j=mtrarr[iiblk].ia[ilist];j<mtrarr[iiblk].ia[ilist+1];j++) {
				jcol = mtrarr[iiblk].ja[j];
				jjblk = mtrarr[iiblk].ja2[j];
				if (imaskblk[jjblk] == icycle) {
					aloc[nzjaloc] = mtrarr[iiblk].a[j];
					nzjaloc++;
				} else {
					auxloc = mtrarr[iiblk].a[j];
					if (auxloc < 0.0e0) auxloc = -auxloc;
					auxloc *= _theta;
					ibsd = _fct.ibsdia[iblk] +irow;
					jbsd = _fct.ibsdia[jjblk]+jcol;
					_fct.diagg[ibsd] += auxloc;
					_fct.diagg[jbsd] += auxloc;
				};
			};
		};

// Free previous data of the block

		delete [] mtrarr[iiblk].list;
		delete [] mtrarr[iiblk].list2;
		delete [] mtrarr[iiblk].ia;
		delete [] mtrarr[iiblk].ja;
		delete [] mtrarr[iiblk].ja2;
		delete [] mtrarr[iiblk].ja3;
		delete [] mtrarr[iiblk].a;

// Register new data

		mtrarr[iiblk].list   = listnew;
		mtrarr[iiblk].list2  = list2new;
		mtrarr[iiblk].ia     = ialoc;
		mtrarr[iiblk].ja     = jaloc;
		mtrarr[iiblk].ja2    = ja2loc;
		mtrarr[iiblk].ja3    = ja3loc;
		mtrarr[iiblk].a      = aloc;

		mtrarr[iiblk].nsupc  = nlistnew;
		mtrarr[iiblk].nsupr  = nlistnew;
		mtrarr[iiblk].nlist  = nlistnew;
		mtrarr[iiblk].nlist2 = nlistnew;
		mtrarr[iiblk].nzja   = nzjaloc;
		mtrarr[iiblk].nzja2  = nzjaloc;
		mtrarr[iiblk].nzja3  = nzjaloc;
		mtrarr[iiblk].nza    = nzjaloc;
		mtrarr[iiblk].nzatot = nzjaloc;

// Free work memory

		delete [] listcount;

	};

	delete [] imaskblk;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Create diagonal update matrix
//========================================================================================
CSMatrixR CGSMatrixR::CreateDiagonalUpdate2Index (int _nlistblk, int *_listblk, // Create diagonal update matrix
																	CFctDiagR &_fctdiag) {

	const char *funcname = "CreateDiagonalUpdate2Index";

	double dzero = 0.0e0;

// Count maximal size and allocate the mask

	int nloc = 0;

	int i, iblk;

	for (i=0;i<_nlistblk;i++) {
		iblk = _listblk[i];
		nloc += blksr[iblk+1]-blksr[iblk];
	};

	int *imask;

	imask = new int [nloc];
	if (!imask) MemoryFail (funcname);

	for (i=0;i<nloc;i++) imask[i] = 1;

// Fill the mask and count the number of elements

	int nflt = 0;

	nloc = 0;

	int j, jloc;

	for (i=0;i<_nlistblk;i++) {
		iblk = _listblk[i];
		for (j=blksr[iblk];j<blksr[iblk+1];j++) {
			jloc = j-blksr[iblk];
			if (_fctdiag.diagg[iblk][jloc] == dzero) {
				imask[nloc] = 0;
			} else {
				nflt++;
			};
			nloc++;
		};
	};

// Create new matrix

	CSMatrixR temp (nflt, nflt, nflt, nflt, nflt);

	nflt = 0;
	nloc = 0;

	temp.ia[0] = 0;

	for (i=0;i<_nlistblk;i++) {
		iblk = _listblk[i];
		for (j=blksr[iblk];j<blksr[iblk+1];j++) {
			jloc = j-blksr[iblk];
			if (imask[nloc] > 0) {
				temp.list[nflt] = jloc;
				temp.list2[nflt] = iblk;
				temp.ia[nflt+1] = nflt+1;
				temp.ja[nflt] = jloc;
				temp.ja2[nflt] = iblk;
				temp.a[nflt] = _fctdiag.diagg[iblk][jloc];
				nflt++;
			};
			_fctdiag.diagg[iblk][jloc] = dzero;
			nloc++;
		};
	};

	temp.nsupc  = nflt;
	temp.nsupr  = nflt;
	temp.nlist  = nflt;
	temp.nlist2 = nflt;
	temp.nzja   = nflt;
	temp.nzja2  = nflt;
	temp.nzja3  = 0;
	temp.nza    = nflt;
	temp.nzatot = nflt;

// Free work arrays

	delete [] imask;

	return temp;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Update diagonal data
//========================================================================================
void CGSMatrixR::UpdateDiagonal2Index (CSMatrixR &_mtrd, // Update diagonal data
													CFctDiagR &_fctdiag) {

//	const char *funcname = "UpdateDiagonal2Index";

// Update diagonal values

	int i, irow, iblk;

	for (i=0;i<_mtrd.nlist;i++) {
		irow = _mtrd.list[i];
		iblk = _mtrd.list2[i];
		_fctdiag.diagg[iblk][irow] += _mtrd.a[i];
	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Load imask data
//========================================================================================
void CFctDiagR::Ilu2LoadIMask (// Load imask data
								int _nlstblk, int *_lstblk, int *_iab, int *_jab,
								const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2LoadIMask";

// Add lists

	int *imaskblk, *listblk;

	int icycleloc = -1;

	imaskblk = new int [_gmtra.nblksr];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [_gmtra.nblksr];
	if (!listblk) MemoryFail (funcname);

	int i;

	for (i=0;i<_gmtra.nblksr;i++) {
		imaskblk[i] = icycleloc;
	};

	int nlistloc = 0;

	icycleloc++;

	int ilist, iblk, j, jblk;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		iblk = _lstblk[ilist];
		for (j=_iab[iblk];j<_iab[iblk+1];j++) {
			jblk = _jab[j];
			if (imaskblk[jblk] != icycleloc) {
				listblk[nlistloc] = jblk;
				imaskblk[jblk] = icycleloc;
				nlistloc++;
			};
		};
	};

	qsort (listblk, nlistloc, sizeof(int), compint);

// Assign reference array

	int jndbeg, jndend;

	int nzdia = 0;
	for (ilist=0;ilist<nlistloc;ilist++) {
		jblk = listblk[ilist];
		jndbeg = _gmtra.blksc[jblk];
		jndend = _gmtra.blksc[jblk+1]-1;
		int nloc = jndend-jndbeg+1;
		_fct.ibsimask[jblk] = nzdia;
		nzdia += nloc;
	};

// Allocate integer arrays

	delete [] _fct.imask;
	delete [] _fct.imasklev;
	delete [] _fct.ind2arr;

	_fct.imask = new int [nzdia];
	if (!_fct.imask) MemoryFail (funcname);
	_fct.imasklev = new int [nzdia];
	if (!_fct.imasklev) MemoryFail (funcname);
	_fct.ind2arr = new CInd2Int [nzdia];
	if (!_fct.ind2arr) MemoryFail (funcname);

// Free work arrays

	delete [] imaskblk;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CFctDiagR: Unload imask data
//========================================================================================
void CFctDiagR::Ilu2UnLoadIMask (CFctR &_fct) {// Unload imask data

	const char *funcname = "Ilu2UnLoadIMask";

	delete [] _fct.imask;
	delete [] _fct.imasklev;
	delete [] _fct.ind2arr;

	_fct.imask = new int [0];
	if (!_fct.imask) MemoryFail (funcname);
	_fct.imasklev = new int [0];
	if (!_fct.imasklev) MemoryFail (funcname);
	_fct.ind2arr = new CInd2Int [0];
	if (!_fct.ind2arr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CFctDiagR: Load diagonal data for the computations with the block row
//========================================================================================
void CFctDiagR::Ilu2LoadDiagonalData2Index (int _iblkrow, // Load diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk, int _nlstblkcol, int *_lstblkcol,
										const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2LoadDiagonalData2Index";

// Delete previously stored in _fct arrays with undefined length

	delete [] _fct.dpiv;
	delete [] _fct.sclpt;
	delete [] _fct.scll;
	delete [] _fct.sclu;
	delete [] _fct.sclinvl;
	delete [] _fct.sclinvu;
	delete [] _fct.diagg;
	delete [] _fct.fmaskl;
	delete [] _fct.fmasku;

// Compute the memory required to store all required arrays

	int nzscl=0;
	int nzdpiv=0;
	int nzfmask=0;

	int ilist, jblk, jndbeg, jndend, nj;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		jndbeg = _gmtra.blksc[jblk];
		jndend = _gmtra.blksc[jblk+1]-1;
		nj = jndend-jndbeg+1;
		nzscl += nj;
		nzdpiv += nj;
	};

	nzfmask = nzdpiv;

// Allocate all arrays

	double *dpivloc;
	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc, *fmasklloc, *fmaskuloc;

	dpivloc = new double [nzdpiv];
	if (!dpivloc) MemoryFail (funcname);
	sclptloc = new double [nzdpiv];
	if (!sclptloc) MemoryFail (funcname);
	scllloc = new double [nzscl];
	if (!scllloc) MemoryFail (funcname);
	scluloc = new double [nzscl];
	if (!scluloc) MemoryFail (funcname);
	sclinvlloc = new double [nzscl];
	if (!sclinvlloc) MemoryFail (funcname);
	sclinvuloc = new double [nzscl];
	if (!sclinvuloc) MemoryFail (funcname);
	diaggloc = new double [nzscl];
	if (!diaggloc) MemoryFail (funcname);
	fmasklloc = new double [nzfmask];
	if (!fmasklloc) MemoryFail (funcname);
	fmaskuloc = new double [nzfmask];
	if (!fmaskuloc) MemoryFail (funcname);

	_fct.dpiv    = dpivloc;
	_fct.sclpt   = sclptloc;
	_fct.scll    = scllloc;
	_fct.sclu    = scluloc;
	_fct.sclinvl = sclinvlloc;
	_fct.sclinvu = sclinvuloc;
	_fct.diagg   = diaggloc;
	_fct.fmaskl  = fmasklloc;
	_fct.fmasku  = fmaskuloc;

// Load necessary data

	int kii;
	double *darr;

	nzscl = 0;
	nzdpiv = 0;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		jndbeg = _gmtra.blksc[jblk];
		jndend = _gmtra.blksc[jblk+1]-1;
		int nzloc = diasize[jblk];
		int nloc = jndend-jndbeg+1;
		darr = sclpt[jblk];
		for (kii=0;kii<nloc;kii++) sclptloc[nzdpiv+kii] = darr[kii];
		if (nfiles == 0) {

			darr = scll[jblk];
			for (kii=0;kii<nzloc;kii++) scllloc[nzscl+kii] = darr[kii];
			darr = sclu[jblk];
			for (kii=0;kii<nzloc;kii++) scluloc[nzscl+kii] = darr[kii];
			darr = sclinvl[jblk];
			for (kii=0;kii<nzloc;kii++) sclinvlloc[nzscl+kii] = darr[kii];
			darr = sclinvu[jblk];
			for (kii=0;kii<nzloc;kii++) sclinvuloc[nzscl+kii] = darr[kii];
			darr = diagg[jblk];
			for (kii=0;kii<nzloc;kii++) diaggloc[nzscl+kii] = darr[kii];

		} else {

			int ifileloc  = blk2file[jblk*5];
			int ibsloc    = blk2bs  [jblk*5];

			FGet (files[ifileloc],nzloc,scllloc+nzscl,ibsloc);

			ifileloc  = blk2file[jblk*5+1];
			ibsloc    = blk2bs  [jblk*5+1];

			FGet (files[ifileloc],nzloc,scluloc+nzscl,ibsloc);

			ifileloc  = blk2file[jblk*5+2];
			ibsloc    = blk2bs  [jblk*5+2];

			FGet (files[ifileloc],nzloc,sclinvlloc+nzscl,ibsloc);

			ifileloc  = blk2file[jblk*5+3];
			ibsloc    = blk2bs  [jblk*5+3];

			FGet (files[ifileloc],nzloc,sclinvuloc+nzscl,ibsloc);

			ifileloc  = blk2file[jblk*5+4];
			ibsloc    = blk2bs  [jblk*5+4];

			FGet (files[ifileloc],nzloc,diaggloc+nzscl,ibsloc);

		};

		nzscl += diasize[jblk];
		nzdpiv += nloc;

	};

// Assign reference arrays

	int nzdia = 0;
	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		jndbeg = _gmtra.blksc[jblk];
		jndend = _gmtra.blksc[jblk+1]-1;
		int nloc = jndend-jndbeg+1;
		_fct.ibsdia[jblk] = nzdia;
		_fct.ibsimask[jblk] = nzdia;
		_fct.ibsfmask[jblk] = nzdia;
		nzdia += nloc;
	};

// Allocate integer arrays

	int nrowstot = 0;

	for (ilist=0;ilist<_nlstblkcol;ilist++) {
		jblk = _lstblkcol[ilist];
		jndbeg = _gmtra.blksc[jblk];
		jndend = _gmtra.blksc[jblk+1]-1;
		int nloc = jndend-jndbeg+1;
		nrowstot += nloc;
	};

	int nmax = nrowstot;
	if (nmax < nzdpiv) nmax = nzdpiv;

	delete [] _fct.imask;
	delete [] _fct.imasklev;
	delete [] _fct.lstloc;
	delete [] _fct.lstloc2;
	delete [] _fct.iv;
	delete [] _fct.madj;
	delete [] _fct.ibegm;
	delete [] _fct.lstupd;
	delete [] _fct.lstupd2;
	delete [] _fct.ig;
	delete [] _fct.jg;
	delete [] _fct.jg2;
	delete [] _fct.jg3;
	delete [] _fct.addrg;
	delete [] _fct.jgblk;
	delete [] _fct.jg2blk;
	delete [] _fct.jg3blk;
	delete [] _fct.adrgblk;
	delete [] _fct.diarr;
	delete [] _fct.ind2arr;

	_fct.imask = new int [nmax+1];
	if (!_fct.imask) MemoryFail (funcname);
	_fct.imasklev = new int [nmax+1];
	if (!_fct.imasklev) MemoryFail (funcname);
	_fct.lstloc = new int [nmax+1];
	if (!_fct.lstloc) MemoryFail (funcname);
	_fct.lstloc2 = new int [nmax+1];
	if (!_fct.lstloc2) MemoryFail (funcname);
	_fct.iv = new int [nmax+1];
	if (!_fct.iv) MemoryFail (funcname);
	_fct.madj = new int [nmax+1];
	if (!_fct.madj) MemoryFail (funcname);
	_fct.ibegm = new int [nmax+1];
	if (!_fct.ibegm) MemoryFail (funcname);
	_fct.lstupd = new int [nmax+1];
	if (!_fct.lstupd) MemoryFail (funcname);
	_fct.lstupd2 = new int [nmax+1];
	if (!_fct.lstupd2) MemoryFail (funcname);
	_fct.ig = new int [nmax+1];
	if (!_fct.ig) MemoryFail (funcname);
	_fct.jg = new int * [nmax+1];
	if (!_fct.jg) MemoryFail (funcname);
	_fct.jg2 = new int * [nmax+1];
	if (!_fct.jg2) MemoryFail (funcname);
	_fct.jg3 = new int * [nmax+1];
	if (!_fct.jg3) MemoryFail (funcname);
	_fct.addrg = new int * [nmax+1];
	if (!_fct.addrg) MemoryFail (funcname);
	_fct.jgblk = new int * [nmax+1];
	if (!_fct.jgblk) MemoryFail (funcname);
	_fct.jg2blk = new int * [nmax+1];
	if (!_fct.jg2blk) MemoryFail (funcname);
	_fct.jg3blk = new int * [nmax+1];
	if (!_fct.jg3blk) MemoryFail (funcname);
	_fct.adrgblk = new int * [nmax+1];
	if (!_fct.adrgblk) MemoryFail (funcname);
	_fct.diarr = new CDoubleInt [nmax+1];
	if (!_fct.diarr) MemoryFail (funcname);
	_fct.ind2arr = new CInd2Int [nmax+1];
	if (!_fct.ind2arr) MemoryFail (funcname);

// Allocate pointers to float arrays

	delete [] _fct.gl;
	delete [] _fct.gu;
	delete [] _fct.glblk;
	delete [] _fct.gublk;
	delete [] _fct.adrupdl;
	delete [] _fct.adrupdu;

	_fct.gl = new double * [nmax+1];
	if (!_fct.gl) MemoryFail (funcname);
	_fct.gu = new double * [nmax+1];
	if (!_fct.gu) MemoryFail (funcname);
	_fct.glblk = new double * [nmax+1];
	if (!_fct.glblk) MemoryFail (funcname);
	_fct.gublk = new double * [nmax+1];
	if (!_fct.gublk) MemoryFail (funcname);
	_fct.adrupdl = new double * [nmax+1];
	if (!_fct.adrupdl) MemoryFail (funcname);
	_fct.adrupdu = new double * [nmax+1];
	if (!_fct.adrupdu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CFctDiagR: Unload diagonal data for the computations with the block row
//========================================================================================
void CFctDiagR::Ilu2UnloadDiagonalData2Index (int _iblkrow, // Unload diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk,
										const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2UnloadDiagonalData2Index";

// Store diagg data

	int ilist, jblk, kii;
	double *carr;

	int nzscl = 0;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		int nzloc = diasize[jblk];

		if (nfiles == 0) {
			carr = diagg[jblk];
			for (kii=0;kii<nzloc;kii++) carr[kii] = _fct.diagg[nzscl+kii];
		} else {

			int ifileloc  = blk2file[jblk*5+4];
			int ibsloc    = blk2bs  [jblk*5+4];

			FPut (files[ifileloc],nzloc,_fct.diagg+nzscl,ibsloc);

		};


		nzscl += diasize[jblk];

	};

// Delete previously stored in _fct arrays with undefined length

	delete [] _fct.dpiv;
	delete [] _fct.sclpt;
	delete [] _fct.scll;
	delete [] _fct.sclu;
	delete [] _fct.sclinvl;
	delete [] _fct.sclinvu;
	delete [] _fct.diagg;
	delete [] _fct.fmaskl;
	delete [] _fct.fmasku;

// Allocate all arrays with dummy length

	double *dpivloc;
	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc, *fmasklloc, *fmaskuloc;

	dpivloc = new double [1];
	if (!dpivloc) MemoryFail (funcname);
	sclptloc = new double [1];
	if (!sclptloc) MemoryFail (funcname);
	scllloc = new double [1];
	if (!scllloc) MemoryFail (funcname);
	scluloc = new double [1];
	if (!scluloc) MemoryFail (funcname);
	sclinvlloc = new double [1];
	if (!sclinvlloc) MemoryFail (funcname);
	sclinvuloc = new double [1];
	if (!sclinvuloc) MemoryFail (funcname);
	diaggloc = new double [1];
	if (!diaggloc) MemoryFail (funcname);
	fmasklloc = new double [1];
	if (!fmasklloc) MemoryFail (funcname);
	fmaskuloc = new double [1];
	if (!fmaskuloc) MemoryFail (funcname);

	_fct.dpiv    = dpivloc;
	_fct.sclpt   = sclptloc;
	_fct.scll    = scllloc;
	_fct.sclu    = scluloc;
	_fct.sclinvl = sclinvlloc;
	_fct.sclinvu = sclinvuloc;
	_fct.diagg   = diaggloc;
	_fct.fmaskl  = fmasklloc;
	_fct.fmasku  = fmaskuloc;

	int nmax = 1;

	delete [] _fct.imask;
	delete [] _fct.imasklev;
	delete [] _fct.lstloc;
	delete [] _fct.lstloc2;
	delete [] _fct.iv;
	delete [] _fct.madj;
	delete [] _fct.madj2;
	delete [] _fct.ibegm;
	delete [] _fct.ibegm2;
	delete [] _fct.ibsblk;
	delete [] _fct.lstupd;
	delete [] _fct.lstupd2;
	delete [] _fct.ig;
	delete [] _fct.jg;
	delete [] _fct.jg2;
	delete [] _fct.jg3;
	delete [] _fct.addrg;
	delete [] _fct.jgblk;
	delete [] _fct.jg2blk;
	delete [] _fct.jg3blk;
	delete [] _fct.adrgblk;
	delete [] _fct.diarr;
	delete [] _fct.ind2arr;

	_fct.imask = new int [nmax+1];
	if (!_fct.imask) MemoryFail (funcname);
	_fct.imasklev = new int [nmax+1];
	if (!_fct.imasklev) MemoryFail (funcname);
	_fct.lstloc = new int [nmax+1];
	if (!_fct.lstloc) MemoryFail (funcname);
	_fct.lstloc2 = new int [nmax+1];
	if (!_fct.lstloc2) MemoryFail (funcname);
	_fct.iv = new int [nmax+1];
	if (!_fct.iv) MemoryFail (funcname);
	_fct.madj = new int [nmax+1];
	if (!_fct.madj) MemoryFail (funcname);
	_fct.madj2 = new int [nmax+1];
	if (!_fct.madj2) MemoryFail (funcname);
	_fct.ibegm = new int [nmax+1];
	if (!_fct.ibegm) MemoryFail (funcname);
	_fct.ibegm2 = new int [nmax+1];
	if (!_fct.ibegm2) MemoryFail (funcname);
	_fct.ibsblk = new int [nmax+1];
	if (!_fct.ibsblk) MemoryFail (funcname);
	_fct.lstupd = new int [nmax+1];
	if (!_fct.lstupd) MemoryFail (funcname);
	_fct.lstupd2 = new int [nmax+1];
	if (!_fct.lstupd2) MemoryFail (funcname);
	_fct.ig = new int [nmax+1];
	if (!_fct.ig) MemoryFail (funcname);
	_fct.jg = new int * [nmax+1];
	if (!_fct.jg) MemoryFail (funcname);
	_fct.jg2 = new int * [nmax+1];
	if (!_fct.jg2) MemoryFail (funcname);
	_fct.jg3 = new int * [nmax+1];
	if (!_fct.jg3) MemoryFail (funcname);
	_fct.addrg = new int * [nmax+1];
	if (!_fct.addrg) MemoryFail (funcname);
	_fct.jgblk = new int * [nmax+1];
	if (!_fct.jgblk) MemoryFail (funcname);
	_fct.jg2blk = new int * [nmax+1];
	if (!_fct.jg2blk) MemoryFail (funcname);
	_fct.jg3blk = new int * [nmax+1];
	if (!_fct.jg3blk) MemoryFail (funcname);
	_fct.adrgblk = new int * [nmax+1];
	if (!_fct.adrgblk) MemoryFail (funcname);
	_fct.diarr = new CDoubleInt [nmax+1];
	if (!_fct.diarr) MemoryFail (funcname);
	_fct.ind2arr = new CInd2Int [nmax+1];
	if (!_fct.ind2arr) MemoryFail (funcname);

// Allocate pointers to float arrays

	delete [] _fct.gl;
	delete [] _fct.gu;
	delete [] _fct.glblk;
	delete [] _fct.gublk;
	delete [] _fct.adrupdl;
	delete [] _fct.adrupdu;

	_fct.gl = new double * [nmax+1];
	if (!_fct.gl) MemoryFail (funcname);
	_fct.gu = new double * [nmax+1];
	if (!_fct.gu) MemoryFail (funcname);
	_fct.glblk = new double * [nmax+1];
	if (!_fct.glblk) MemoryFail (funcname);
	_fct.gublk = new double * [nmax+1];
	if (!_fct.gublk) MemoryFail (funcname);
	_fct.adrupdl = new double * [nmax+1];
	if (!_fct.adrupdl) MemoryFail (funcname);
	_fct.adrupdu = new double * [nmax+1];
	if (!_fct.adrupdu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CFctR: Init zero data diagonal block
//========================================================================================
CSMatrixR CFctR::CreateBlock2Index (int _iblk, const CGSMatrix &_gmtra) { // Init zero data diagonal block

//	const char *funcname = "CreateBlock2Index";

	double dzero = 0.0e0;

// Count the total number of nonzero elements

	int nlistloc = _gmtra.blksr[_iblk+1]-_gmtra.blksr[_iblk];

// Store the data

	CSMatrixR temp (nlistloc,nlistloc,nlistloc,nlistloc,nlistloc);

	temp.AllocateJa3 (nlistloc);

	int i;
	for (i=0;i<nlistloc;i++) temp.list[i] = i;
	for (i=0;i<nlistloc;i++) temp.list2[i] = _iblk;
	for (i=0;i<=nlistloc;i++) temp.ia[i] = i;
	for (i=0;i<nlistloc;i++) temp.ja[i] = i;
	for (i=0;i<nlistloc;i++) temp.ja2[i] = _iblk;
	for (i=0;i<nlistloc;i++) temp.ja3[i] = 1;
	for (i=0;i<nlistloc;i++) temp.a[i] = dzero;

// Fill control data

	temp.m      = _gmtra.m;
	temp.n      = _gmtra.n;
	temp.nsupr  = nlistloc;
	temp.nsupc  = nlistloc;
	temp.nlist  = nlistloc;
	temp.nlist2 = nlistloc;
	temp.nzja   = nlistloc;
	temp.nzja2  = nlistloc;
	temp.nzja3  = nlistloc;
	temp.nza    = nlistloc;
	temp.nzatot = nlistloc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Filter block row according to the list of indices
//========================================================================================
void CSMatrixR::FilterBlockRow2Index (bool _is_3index, int _nlistblk, int *_listblk, // Filter block row according to the list of indices
													int &_icycleblk, int *_imaskblk) {

	const char *funcname = "FilterBlockRow2Index";

// Init the block mask

	int i, jjblk;

	_icycleblk++;
	for (i=0;i<_nlistblk;i++) {
		jjblk = _listblk[i];
		_imaskblk[jjblk] = _icycleblk;
	};

// Scan block row and count the new number of elements

	int *listcount;

	listcount = new int [nlist];
	if (!listcount) MemoryFail (funcname);

	int nzflt = 0;
	int nz = 0;

	int ilist, nzloc, j, jj2, jj3, jj;

	for (ilist=0;ilist<nlist;ilist++) {
		nzloc = nz;
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj2 = ja2[j];
			if (_imaskblk[jj2] == _icycleblk) {
				nz++;
			} else {
				nzflt++;
			};
		};
		listcount[ilist] = nz-nzloc;
	};

	if (nzflt == 0) {
		delete [] listcount;
		return;
	};

// Create ia, ja and ja2 arrays

	int *ialoc;

	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		ialoc[ilist+1] = ialoc[ilist] + listcount[ilist];
	};

	int nzjanew = ialoc[nlist];

	int *jaloc;
	int *ja2loc;
	int *ja3loc;
	double *aloc;

	jaloc = new int [nzjanew];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjanew];
	if (!ja2loc) MemoryFail (funcname);
	if (_is_3index) {
		ja3loc = new int [nzjanew];
		if (!ja3loc) MemoryFail (funcname);
	} else {
		ja3loc = new int [0];
		if (!ja3loc) MemoryFail (funcname);
	};
	aloc = new double [nzjanew];
	if (!aloc) MemoryFail (funcname);

	int nzjaloc = 0;
	jj3 = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			jj2 = ja2[j];
			if (_is_3index) {
				jj3 = ja3[j];
			};
			if (_imaskblk[jj2] == _icycleblk) {
				jaloc[nzjaloc] = jj;
				ja2loc[nzjaloc] = jj2;
				if (_is_3index) {
					ja3loc[nzjaloc] = jj3;
				};
				aloc[nzjaloc] = a[j];
				nzjaloc++;
			};
		};
	};

// Free previous data of the block

	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;
	delete [] a;

// Register new data

	ia     = ialoc;
	ja     = jaloc;
	ja2    = ja2loc;
	ja3    = ja3loc;
	a      = aloc;

	nzja   = nzjaloc;
	nzja2  = nzjaloc;
	if (_is_3index) {
		nzja3  = nzjaloc;
	} else {
		nzja3  = 0;
	};
	nza    = nzjaloc;
	nzatot = nzjaloc;

// Free work memory

	delete [] listcount;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Filter Schur block row/column according to the indices
//========================================================================================
void CGSMatrixR::FilterSchurBlock2Index (int _iblk, // Filter Schur block row/column according to the indices
											int _iblkbeg, int _iblkend, int _iblkbegrest, double _theta,
											CFctR &_fct) {

	const char *funcname = "FilterSchurBlock2Index";

// Create new list of supernode rows

	int *listcount;

	int nlistloc = mtrarr[_iblk].nlist;

	listcount = new int [nlistloc];
	if (!listcount) MemoryFail (funcname);

	int ilist, nz, j, jj, jj2;

	for (ilist=0;ilist<nlistloc;ilist++) {
		nz = 0;
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj2 = mtrarr[_iblk].ja2[j];
			if ((jj2>=_iblkbeg && jj2 <=_iblkend) || jj2>=_iblkbegrest) nz++;
		};
		listcount[ilist] = nz;
	};

	int nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) nlistnew++;
	};

	int *listnew;
	int *list2new;

	listnew = new int [nlistnew];
	if (!listnew) MemoryFail (funcname);
	list2new = new int [nlistnew];
	if (!list2new) MemoryFail (funcname);

	nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) {
			listnew[nlistnew] = mtrarr[_iblk].list[ilist];
			list2new[nlistnew] = mtrarr[_iblk].list2[ilist];
			nlistnew++;
		};
	};

// Create ia, ja and ja2 arrays

	int *ialoc;

	ialoc = new int [nlistnew+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) {
			ialoc[nlistnew+1] = ialoc[nlistnew] + listcount[ilist];
			nlistnew++;
		};
	};

	int nzjanew = ialoc[nlistnew];

	int *jaloc;
	int *ja2loc;
	int *ja3loc;

	jaloc = new int [nzjanew];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjanew];
	if (!ja2loc) MemoryFail (funcname);
	ja3loc = new int [nzjanew];
	if (!ja3loc) MemoryFail (funcname);

	int nzjaloc = 0;

	int jlev;

	for (ilist=0;ilist<nlistloc;ilist++) {
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj = mtrarr[_iblk].ja[j];
			jj2 = mtrarr[_iblk].ja2[j];
			jlev = mtrarr[_iblk].ja3[j];
			if ((jj2>=_iblkbeg && jj2 <=_iblkend) || jj2>=_iblkbegrest) {
				jaloc[nzjaloc] = jj;
				ja2loc[nzjaloc] = jj2;
				ja3loc[nzjaloc] = jlev;
				nzjaloc++;
			};
		};
	};

// Create a and modify the corresponding diagg array

	double *aloc;

	aloc = new double [nzjaloc];
	if (!aloc) MemoryFail (funcname);

	int nzaloc = 0;

	int irow, iblk;

	for (ilist=0;ilist<nlistloc;ilist++) {
		irow = mtrarr[_iblk].list[ilist];
		iblk = mtrarr[_iblk].list2[ilist];
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj = mtrarr[_iblk].ja[j];
			jj2 = mtrarr[_iblk].ja2[j];
			if ((jj2>=_iblkbeg && jj2 <=_iblkend) || jj2>=_iblkbegrest) {
				aloc[nzaloc] = mtrarr[_iblk].a[j];
				nzaloc++;
			} else {
				int ibsd = _fct.ibsdia[iblk]+irow;
				int jbsd = _fct.ibsdia[jj2]+jj;
				double auxloc = mtrarr[_iblk].a[j];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc *= _theta;
				_fct.diagg[ibsd] += auxloc;
				_fct.diagg[jbsd] += auxloc;
			};
		};
	};

// Free previous data of the block

	delete [] mtrarr[_iblk].list;
	delete [] mtrarr[_iblk].list2;
	delete [] mtrarr[_iblk].ia;
	delete [] mtrarr[_iblk].ja;
	delete [] mtrarr[_iblk].ja2;
	delete [] mtrarr[_iblk].ja3;
	delete [] mtrarr[_iblk].a;

// Register new data

	mtrarr[_iblk].list   = listnew;
	mtrarr[_iblk].list2  = list2new;
	mtrarr[_iblk].ia     = ialoc;
	mtrarr[_iblk].ja     = jaloc;
	mtrarr[_iblk].ja2    = ja2loc;
	mtrarr[_iblk].ja3    = ja3loc;
	mtrarr[_iblk].a      = aloc;

	mtrarr[_iblk].nsupc  = nlistnew;
	mtrarr[_iblk].nsupr  = nlistnew;
	mtrarr[_iblk].nlist  = nlistnew;
	mtrarr[_iblk].nlist2 = nlistnew;
	mtrarr[_iblk].nzja   = nzjaloc;
	mtrarr[_iblk].nzja2  = nzjaloc;
	mtrarr[_iblk].nzja3  = nzjaloc;
	mtrarr[_iblk].nza    = nzaloc;
	mtrarr[_iblk].nzatot = nzaloc;

// Free work memory

	delete [] listcount;

};

// Author: Kharchenko S.A.
// CFct: Init arrays that support global update factorization routines
//========================================================================================
void CFct::AllocateGlobalArrays2Index (int _nblks, int *_blks, int _nlistschur, int *_listschur) { // Init arrays that support global update factorization routines

	const char *funcname = "AllocateGlobalArrays2Index";

// Determine the size of the arrays to be stored

	delete [] ibsblkgl;

	ibsblkgl = new int [_nblks+1];
	if (!ibsblkgl) MemoryFail (funcname);

	int i;

	for (i=0;i<_nblks+1;i++) {
		ibsblkgl[i] = -1;
	};

	int nloc = 0;

	int iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		ibsblkgl[iblk] = nloc;
		nloc += _blks[iblk+1]-_blks[iblk];
	};

// Free and reallocate all arrays

	delete [] irw2lst;
	delete [] ivgl;
	delete [] madjgl;
	delete [] madj2gl;
	delete [] ibegmgl;
	delete [] ibegm2gl;
	delete [] imaskgl;
	delete [] ind2arrgl;

	irw2lst = new int [nloc+1];
	if (!irw2lst) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		irw2lst[i] = -1;
	};
	ivgl = new int [nloc+1];
	if (!ivgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ivgl[i] = -1;
	};
	madjgl = new int [nloc+1];
	if (!madjgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		madjgl[i] = -1;
	};
	madj2gl = new int [nloc+1];
	if (!madj2gl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		madj2gl[i] = -1;
	};
	ibegmgl = new int [nloc+1];
	if (!ibegmgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ibegmgl[i] = -1;
	};
	ibegm2gl = new int [nloc+1];
	if (!ibegm2gl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ibegm2gl[i] = -1;
	};
	icyclegl = -1;
	imaskgl = new int [nloc+1];
	if (!imaskgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		imaskgl[i] = -1;
	};
	ind2arrgl = new CInd2Int [nloc+1];
	if (!ind2arrgl) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CFct: Free arrays that support global update factorization routines
//========================================================================================
void CFct::FreeGlobalArrays2Index () { // Free arrays that support global update factorization routines

	const char *funcname = "FreeGlobalArrays2Index";

// Free and reallocate all arrays

	int nloc = 0;

	delete [] irw2lst;
	delete [] ivgl;
	delete [] madjgl;
	delete [] madj2gl;
	delete [] ibegmgl;
	delete [] ibegm2gl;
	delete [] ibsblkgl;
	delete [] imaskgl;
	delete [] ind2arrgl;

	int i;

	irw2lst = new int [nloc+1];
	if (!irw2lst) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		irw2lst[i] = -1;
	};
	ivgl = new int [nloc+1];
	if (!ivgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ivgl[i] = -1;
	};
	madjgl = new int [nloc+1];
	if (!madjgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		madjgl[i] = -1;
	};
	madj2gl = new int [nloc+1];
	if (!madj2gl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		madj2gl[i] = -1;
	};
	ibegmgl = new int [nloc+1];
	if (!ibegmgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ibegmgl[i] = -1;
	};
	ibegm2gl = new int [nloc+1];
	if (!ibegm2gl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ibegm2gl[i] = -1;
	};
	ibsblkgl = new int [nloc+1];
	if (!ibsblkgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		ibsblkgl[i] = -1;
	};
	icyclegl = -1;
	imaskgl = new int [nloc+1];
	if (!imaskgl) MemoryFail (funcname);
	for (i=0;i<nloc+1;i++) {
		imaskgl[i] = -1;
	};
	ind2arrgl = new CInd2Int [nloc+1];
	if (!ind2arrgl) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CFct: Init global arrays by block data
//========================================================================================
void CFct::InitGlobalArraysByBlockData2Index (CSMatrix &_mtra) { // Init global arrays by block data

//	const char *funcname = "InitGlobalArraysByBlockData2Index";

// Scan all rows

	int nlistloc = _mtra.GetNlist ();
	int *plist = _mtra.GetList ();
	int *plist2 = _mtra.GetList2 ();
	int *pia = _mtra.GetIa ();
	int *pja = _mtra.GetJa ();
	int *pja2 = _mtra.GetJa2 ();

	int i, irow, iblk, ibs, iirow, j, jcol, jblk, jbs, jjcol;

	for (i=0;i<nlistloc;i++) {
		irow = plist[i];
		iblk = plist2[i];
		ibs = ibsblkgl[iblk];
		iirow = ibs+irow;
		irw2lst[iirow] = i;
		j = pia[i];
		ivgl[iirow] = j;
		if (j < pia[i+1]) {
			jcol = pja[j];
			if (jcol < 0) jcol = -jcol-1;
			jblk = pja2[j];
			jbs = ibsblkgl[jblk];
			jjcol = jbs+jcol;
			madjgl[iirow] = ibegmgl[jjcol];
			madj2gl[iirow] = ibegm2gl[jjcol];
			ibegmgl[jjcol] = irow;
			ibegm2gl[jjcol] = iblk;
		};
	};

};

// Author: Kharchenko S.A.
// CFct: ReInit global arrays by block data
//========================================================================================
void CFct::ReInitGlobalArraysByBlockData2Index (CSMatrix &_mtra) { // ReInit global arrays by block data

//	const char *funcname = "ReInitGlobalArraysByBlockData2Index";

// Scan all rows

	int nlistloc = _mtra.GetNlist ();
	int *plist = _mtra.GetList ();
	int *plist2 = _mtra.GetList2 ();
	int *pia = _mtra.GetIa ();

	int i, irow, iblk, ibs, iirow, j;

	for (i=0;i<nlistloc;i++) {
		irow = plist[i];
		iblk = plist2[i];
		ibs = ibsblkgl[iblk];
		iirow = ibs+irow;
		irw2lst[iirow] = i;
		j = pia[i];
		ivgl[iirow] = j;
	};

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif// 
