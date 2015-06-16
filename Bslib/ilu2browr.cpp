//------------------------------------------------------------------------------------------------
// File: ilu2browr.cpp
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

// Author: Kharchenko S.A.
// CFctDiagR: Find supernodes with small singular values
//========================================================================================
void CFctDiagR::FindSups (int _iblkrow, int _nthresh, double *_thresh, // Find supernodes with small singular values
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau,
							int *_nlistarr, int **_listarr,
							CFctR &_fct) {

	const char *funcname = "FindSups";

//	double dzero = 0.0e0;

// Determine the size of the arrays to be stored

	int isup, isupbeg, isupend, nlistl, blai, kii;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	int *listtot;

	listtot = new int [nlistl*_nthresh];
	if (!listtot) MemoryFail (funcname);

// Init scaling data for the current block row

	int j;
	int ibs, ibsv, blai_2;

	double *paloc, *pa;

	for (j=0;j<_nthresh;j++) _nlistarr[j] = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {

		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		blai_2 = blai*blai;

// Assign diagonal supernode

		j = _mtrau.ia[isup-isupbeg];
		ibs = _mtrau.bsa[j];

		ibsv = _gmtra.sprndr[isup]-_gmtra.sprndr[isupbeg];

		paloc = _fct.aloc;
		pa = _mtrau.a+ibs;
		for (kii=0;kii<blai_2;kii++) *paloc++ = *pa++;

// Compute its singular value decomposition

		int blal = blai;
		int info = 0;

		dgesvd_ ("A", "A", &blal, &blal,
					_fct.aloc, &blal, _fct.eig, 
					_fct.uloc, &blal, _fct.vloc, &blal, 
					_fct.work, &_fct.lwork, &info);

		if (info != 0) throw " Error in the Lapack routine DGESVD";

		_fct.ops += 8 * blai_2*blai;

		for (kii=0;kii<blai;kii++) _fct.eig[kii] = sqrt(_fct.eig[kii]);

// Check threshold data and update the lists

		for (int ithresh=0;ithresh<_nthresh;ithresh++) {
			int icount=0;
			for (kii=0;kii<blai;kii++) {
				if (_fct.eig[kii] < _thresh[ithresh]) icount++;
			};
			if (icount > 0) {
				j = _nlistarr[ithresh];
				listtot[ithresh*nlistl+j] = isup;
				_nlistarr[ithresh]++;
			};
		};

	};

// Copy data to the return arrays

	int *listloc;

	for (int ithresh=0;ithresh<_nthresh;ithresh++) {
		j = _nlistarr[ithresh];
		listloc = new int [j];
		if (!listloc) MemoryFail (funcname);
		for (int kii=0;kii<j;kii++) {
			listloc[kii] = listtot[ithresh*nlistl+kii];
		};
		_listarr[ithresh] = listloc;
	};

// Free work memory

	delete [] listtot;

};

// Author: Kharchenko S.A.
// CFctDiagR: Init scaling for the current block row
//========================================================================================
void CFctDiagR::Ilu2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init scaling for the current block row
									const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau,
									CFctR &_fct) {

	const char *funcname = "Ilu2ScaleBlkRow";

	double done = 1.0e0;
	double dzero = 0.0e0;

// Determine the size of the arrays to be stored

	int blamxloc = _mtrau.blamx;

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

	int j, kjj, kiib;
	int ibs, ibsv, blai_2;

	double *peig;

	double *paloc, *pa, *pscl;

	double caux;
	double daux;

	if (_sclpttype == 1) {

		nz = 0;

		for (isup=isupbeg;isup<=isupend;isup++) {

			blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
			blai_2 = blai*blai;

// Search diagonal supernode

			j = _mtrau.ia[isup-isupbeg];
			if (blamxloc == 1) {
				ibs = j;
			} else {
				ibs = _mtrau.bsa[j];
			};

			ibsv = _gmtra.sprndr[isup]-_gmtra.sprndr[isupbeg];

			pa = _mtrau.a+ibs;

			for (kii=0;kii<blai;kii++) {
				caux = pa[kii*blai+kii];
				daux = caux;
				if (daux < 0.0e0) daux = -daux;
				dpivloc[nz+kii] = daux;
				daux = sqrt(daux);
				daux = 1.0e0 / daux;
				sclptloc[nz+kii] = daux;
			};

			nz += blai;

		};

	} else {
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = done;
		for (kii=0;kii<nlocblk;kii++) dpivloc[kii] = 1.0e0;
	};

// Update hystogram

	hyst2.UpdateHystogram (nlocblk,dpivloc);

// Init scaling data for the current block row

	int npt=0;

	nz = 0;

	for (isup=isupbeg;isup<=isupend;isup++) {

		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		blai_2 = blai*blai;

// Assign diagonal supernode

		j = _mtrau.ia[isup-isupbeg];
		if (blamxloc == 1) {
			ibs = j;
		} else {
			ibs = _mtrau.bsa[j];
		};

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

// Compute its singular value decomposition

		int blal = blai;
		int info = 0;

		if (blal != 1) {
			dgesvd_ ("A", "A", &blal, &blal,
						_fct.aloc, &blal, _fct.eig, 
						_fct.uloc, &blal, _fct.vloc, &blal, 
						_fct.work, &_fct.lwork, &info);
			if (info != 0) throw " Error in the Lapack routine DGESVD";
		} else {
			if (_fct.aloc[0] >= 0) {
				_fct.eig[0] = _fct.aloc[0];
				_fct.uloc[0] = 1.0e0;
				_fct.vloc[0] = 1.0e0;
			} else {
				_fct.eig[0] = -_fct.aloc[0];
				_fct.uloc[0] = -1.0e0;
				_fct.vloc[0] = 1.0e0;
			};
		};

		_fct.ops += 8 * blai_2*blai;

// Modify the pivots

		for (kii=0;kii<blai;kii++) {
			if (_fct.eig[kii] < _sclmin) _fct.eig[kii] = _sclmin;
		};

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			daux = _fct.eig[kii];
			daux = sqrt(daux);
			_fct.eig[kii] = daux;
			dpivloc[ibsv+kii] = _fct.eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			peig = _fct.eig;
			paloc = _fct.uloc+kii;
			pscl = sclinvlloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvl[ibsd+kjj*blai+kii] = _uloc[kjj*blai+kii] * _eig[kjj];
				*pscl = *paloc * *peig++;
				paloc += blai;
				pscl += blai;
			};
		};

		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			peig = _fct.eig;
			paloc = _fct.vloc+kiib;
			pscl = sclinvuloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvu[ibsd+kjj*blai+kii] = _vloc[kii*blai+kjj] * _eig[kjj];
				*pscl = *paloc++ * *peig++;
				pscl += blai;
			};
			kiib += blai;
		};

		for (kii=0;kii<blai;kii++) {
			daux = _fct.eig[kii];
			daux = 1.0e0 / daux;
			_fct.eig[kii] = daux;
		};

		peig = _fct.eig;
		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			paloc = _fct.uloc+kiib;
			pscl = scllloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scll[ibsd+kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//				*pscl = conj(*paloc++) * *peig;
				double temp = *paloc++;
				*pscl = temp * *peig;
				pscl += blai;
			};
			peig++;
			kiib += blai;
		};

		peig = _fct.eig;
		for (kii=0;kii<blai;kii++) {
			paloc = _fct.vloc+kii;
			pscl = scluloc+nz+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclu[ibsd+kjj*blai+kii] = _vloc[kjj*blai+kii] * _eig[kii];
//				*pscl = conj(*paloc) * *peig;
				double temp = *paloc;
				*pscl = temp * *peig;
				paloc += blai;
				pscl += blai;
			};
			peig++;
		};

		_fct.ops += 4*blai_2 + 2*blai;

		nz += blai*blai;
		npt += blai;

	};

// Init diagg

	for (kii=0;kii<nz;kii++) diaggloc[kii] = dzero;

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
// CFctDiagR: Store point scaling data in the vector array
//========================================================================================
void CFctDiagR::StorePointScaling (CGSMatrixR &_mtrl) const { // Store point scaling data in the vector array

//	const char *funcname = "StorePointScaling";

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

	CSVector temp (nloc, nblkloc, 1);

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
// CFctDiagR: Store point scaling data in the vector array
//========================================================================================
void CFctDiagR::StorePointScaling (CGSMatrixRS &_mtrl) const { // Store point scaling data in the vector array

//	const char *funcname = "StorePointScaling";

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

	CSVector temp (nloc, nblkloc, 1);

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
// CFctDiagR: Store block scaling data in the vector arrays
//========================================================================================
void CFctDiagR::StoreBlockScaling (char _mtrtype, CGSMatrixR &_mtrlu) const { // Store block scaling data in the vector arrays

//	const char *funcname = "StoreBlockScaling";

// Determine the size of the arrays to be stored

	int nblkloc = 0;
	int nloc = 0;
	int i, kii;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			nloc += diasize[i];
			nblkloc++;
		};
	};

// Allocate temporary structures

	CSVector temp1 (nloc, nblkloc, 1);
	CSVector temp2 (nloc, nblkloc, 1);

// Write scaling into SVectorC

	nblkloc = 0;
	nloc = 0;

	temp1.blksv[0] = 0;
	temp2.blksv[0] = 0;

	double *carr;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			if (_mtrtype == 'L') {

				if (nfiles>0) {

					int ifileloc  = blk2file[i*5];
					int ibsloc    = blk2bs  [i*5];

					FGet (files[ifileloc],diasize[i],temp1.vect+nloc,ibsloc);

					ifileloc  = blk2file[i*5+2];
					ibsloc    = blk2bs  [i*5+2];

					FGet (files[ifileloc],diasize[i],temp2.vect+nloc,ibsloc);

				} else {
					carr = scll[i];
					for (kii=0;kii<diasize[i];kii++) temp1.vect[nloc+kii] = carr[kii];
					carr = sclinvl[i];
					for (kii=0;kii<diasize[i];kii++) temp2.vect[nloc+kii] = carr[kii];
				};

			} else {

				if (nfiles>0) {

					int ifileloc  = blk2file[i*5+1];
					int ibsloc    = blk2bs  [i*5+1];

					FGet (files[ifileloc],diasize[i],temp1.vect+nloc,ibsloc);

					ifileloc  = blk2file[i*5+3];
					ibsloc    = blk2bs  [i*5+3];

					FGet (files[ifileloc],diasize[i],temp2.vect+nloc,ibsloc);

				} else {
					carr = sclu[i];
					for (kii=0;kii<diasize[i];kii++) temp1.vect[nloc+kii] = carr[kii];
					carr = sclinvu[i];
					for (kii=0;kii<diasize[i];kii++) temp2.vect[nloc+kii] = carr[kii];
				};

			};
			temp1.blksv[nblkloc+1] = temp1.blksv[nblkloc] + diasize[i];
			temp1.inda[nblkloc] = _mtrlu.bl2ndc[i];
			temp2.blksv[nblkloc+1] = temp2.blksv[nblkloc] + diasize[i];
			temp2.inda[nblkloc] = _mtrlu.bl2ndc[i];
			nloc += diasize[i];
			nblkloc++;
		};
	};

// Store

	_mtrlu.vector1 = temp1;
	_mtrlu.vector2 = temp2;

};

// Author: Kharchenko S.A.
// CFctDiagR: Store block scaling data in the vector arrays
//========================================================================================
void CFctDiagR::StoreBlockScaling (char _mtrtype, CGSMatrixRS &_mtrlu) const { // Store block scaling data in the vector arrays

//	const char *funcname = "StoreBlockScaling";

// Determine the size of the arrays to be stored

	int nblkloc = 0;
	int nloc = 0;
	int i, kii;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			nloc += diasize[i];
			nblkloc++;
		};
	};

// Allocate temporary structures

	CSVector temp1 (nloc, nblkloc, 1);
	CSVector temp2 (nloc, nblkloc, 1);

// Write scaling into SVectorC

	nblkloc = 0;
	nloc = 0;

	temp1.blksv[0] = 0;
	temp2.blksv[0] = 0;

	double *carr;

	for (i=0;i<nblks;i++) {
		if (diasize[i] != 0) {
			if (_mtrtype == 'L') {

				if (nfiles>0) {

					int ifileloc  = blk2file[i*5];
					int ibsloc    = blk2bs  [i*5];

					FGet (files[ifileloc],diasize[i],temp1.vect+nloc,ibsloc);

					ifileloc  = blk2file[i*5+2];
					ibsloc    = blk2bs  [i*5+2];

					FGet (files[ifileloc],diasize[i],temp2.vect+nloc,ibsloc);

				} else {
					carr = scll[i];
					for (kii=0;kii<diasize[i];kii++) temp1.vect[nloc+kii] = carr[kii];
					carr = sclinvl[i];
					for (kii=0;kii<diasize[i];kii++) temp2.vect[nloc+kii] = carr[kii];
				};

			} else {

				if (nfiles>0) {

					int ifileloc  = blk2file[i*5+1];
					int ibsloc    = blk2bs  [i*5+1];

					FGet (files[ifileloc],diasize[i],temp1.vect+nloc,ibsloc);

					ifileloc  = blk2file[i*5+3];
					ibsloc    = blk2bs  [i*5+3];

					FGet (files[ifileloc],diasize[i],temp2.vect+nloc,ibsloc);

				} else {
					carr = sclu[i];
					for (kii=0;kii<diasize[i];kii++) temp1.vect[nloc+kii] = carr[kii];
					carr = sclinvu[i];
					for (kii=0;kii<diasize[i];kii++) temp2.vect[nloc+kii] = carr[kii];
				};

			};
			temp1.blksv[nblkloc+1] = temp1.blksv[nblkloc] + diasize[i];
			temp1.inda[nblkloc] = _mtrlu.bl2ndc[i];
			temp2.blksv[nblkloc+1] = temp2.blksv[nblkloc] + diasize[i];
			temp2.inda[nblkloc] = _mtrlu.bl2ndc[i];
			nloc += diasize[i];
			nblkloc++;
		};
	};

// Store

	_mtrlu.vector1 = temp1;
	_mtrlu.vector2 = temp2;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Perform point scaling of the block row
//========================================================================================
void CSMatrixRS::PointScaling (int *_sprndr, int *_ibsscl, double *_scl) { // Perform point scaling of the block row

//	const char *funcname = "PointScaling";

// Scan all local supernodes

	int ilist, isup, jsup, j, blai, blaj, ibs, ibsv, jbsv, kii, kjj;
	double caux;

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
// CSMatrixRS: Perform supernode scaling of the block row
//========================================================================================
void CSMatrixRS::BlockScaling (int *_sprndr, int *_ibsscl, double *_scll, double *_sclu) { // Perform supernode scaling of the block row

	const char *funcname = "BlockScaling";

	double dzero = 0.0e0;

// Allocate work array

	double *aloc;
	double *vloc;

	aloc = new double [blamx*blamx];
	if (!aloc) MemoryFail (funcname);
	vloc = new double [blamx*blamx];
	if (!vloc) MemoryFail (funcname);

// Scan all local supernodes

	int ilist, isup, jsup, j, blai, blaj, ibs, ibsv, jbsv, kii, kjj;
	int kjjb, kkk, blaij;
	double aux;
	double *paloc, *pvloc, *pscloc;

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
					aux = dzero;
					for (kkk=0;kkk<blai;kkk++) {
//						aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//						aux += *paloc++ * *pscloc;
						aux += *paloc++ * *pscloc;
						pscloc += blai;
					};
					aloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					aux = dzero;
					paloc = aloc+kii;
					pscloc = _sclu+jbsv+kjj;
					for (kkk=0;kkk<blaj;kkk++) {
//						aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//						aux += *paloc * *pscloc;
						aux += *paloc * *pscloc;
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

// Author: Kharchenko S.A.
// CFctDiagR: Load diagonal data for the computations with the block row
//========================================================================================
void CFctDiagR::Ilu2LoadDiagonalData (int _iblkrow, // Load diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk,
										const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2LoadDiagonalData";

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

	int ilist, jblk, jsupbeg, jsupend;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		jsupbeg = _gmtra.blksc[jblk];
		jsupend = _gmtra.blksc[jblk+1]-1;
		nzscl += diasize[jblk];
		nzdpiv += _gmtra.sprndc[jsupend+1]-_gmtra.sprndc[jsupbeg];
	};

	nzfmask = nzdpiv * _fct.blamx;

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

	int kii, isup, blai;
	double *carr;

	nzscl = 0;
	nzdpiv = 0;

	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		jsupbeg = _gmtra.blksc[jblk];
		jsupend = _gmtra.blksc[jblk+1]-1;
		int nzloc = diasize[jblk];
		int nloc = _gmtra.sprndc[jsupend+1]-_gmtra.sprndc[jsupbeg];
		carr = sclpt[jblk];
		for (kii=0;kii<nloc;kii++) sclptloc[nzdpiv+kii] = carr[kii];
		if (nfiles == 0) {

			carr = scll[jblk];
			for (kii=0;kii<nzloc;kii++) scllloc[nzscl+kii] = carr[kii];
			carr = sclu[jblk];
			for (kii=0;kii<nzloc;kii++) scluloc[nzscl+kii] = carr[kii];
			carr = sclinvl[jblk];
			for (kii=0;kii<nzloc;kii++) sclinvlloc[nzscl+kii] = carr[kii];
			carr = sclinvu[jblk];
			for (kii=0;kii<nzloc;kii++) sclinvuloc[nzscl+kii] = carr[kii];
			carr = diagg[jblk];
			for (kii=0;kii<nzloc;kii++) diaggloc[nzscl+kii] = carr[kii];

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

	int nzdia=0;
	int ibsv=0;
	for (ilist=0;ilist<_nlstblk;ilist++) {
		jblk = _lstblk[ilist];
		for (isup=_gmtra.blksc[jblk];isup<_gmtra.blksc[jblk+1];isup++) {
			blai = _gmtra.sprndr[isup+1] - _gmtra.sprndr[isup];
			_fct.bsdia[isup] = nzdia;
			_fct.adrfmask[isup] = ibsv;
			nzdia += blai*blai;
			ibsv += blai;
		};
	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Unload diagonal data for the computations with the block row
//========================================================================================
void CFctDiagR::Ilu2UnloadDiagonalData (int _iblkrow, // Unload diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk,
										const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2UnloadDiagonalData";

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

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Filter Schur block row/column according to the indices
//========================================================================================
void CGSMatrixRS::FilterSchurBlock (int _iblk, // Filter Schur block row/column according to the indices
									int _isupbeg, int _isupend, int _isupbegrest, double _theta,
									CFctR &_fct) {

	const char *funcname = "FilterSchurBlock";

// Create new list of supernode rows

	int *listcount;

	int nlistloc = mtrarr[_iblk].nlist;

	listcount = new int [nlistloc];
	if (!listcount) MemoryFail (funcname);

	int ilist;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int nz = 0;
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jj = mtrarr[_iblk].ja[j];
			if (jj < 0) jj = -jj;
			if ((jj>=_isupbeg && jj <=_isupend) || jj>=_isupbegrest) nz++;
		};
		listcount[ilist] = nz;
	};

	int nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) nlistnew++;
	};

	int *listnew;

	listnew = new int [nlistnew];
	if (!listnew) MemoryFail (funcname);

	nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) {
			listnew[nlistnew] = mtrarr[_iblk].list[ilist];
			nlistnew++;
		};
	};

// Create zero sprndr and sprndc arrays

	int *sprndrloc;
	int *sprndcloc;

	sprndrloc = new int [nlistnew+1];
	if (!sprndrloc) MemoryFail (funcname);
	sprndcloc = new int [nlistnew+1];
	if (!sprndcloc) MemoryFail (funcname);

	int i;

	for (i=0;i<=nlistnew;i++) sprndrloc[i] = 0;
	for (i=0;i<=nlistnew;i++) sprndcloc[i] = 0;

// Create ia, ja and bsa arrays

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
	int *bsaloc;

	jaloc = new int [nzjanew];
	if (!jaloc) MemoryFail (funcname);
	bsaloc = new int [nzjanew];
	if (!bsaloc) MemoryFail (funcname);

	int nzjaloc = 0;
	int nzaloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = mtrarr[_iblk].list[ilist];
		int blai = sprndr[isup+1]-sprndr[isup];
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jsup = mtrarr[_iblk].ja[j];
			int jsupnew=jsup;
			if (jsupnew < 0) jsupnew = -jsupnew;
			int blaj = sprndr[jsupnew+1]-sprndr[jsupnew];
			if ((jsupnew>=_isupbeg && jsupnew <=_isupend) || jsupnew>=_isupbegrest) {
				jaloc[nzjaloc] = jsup;
				bsaloc[nzjaloc] = nzaloc;
				nzjaloc++;
				nzaloc += blai*blaj;
			};
		};
	};

// Create a and modify the corresponding diagg array

	double *aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	nzaloc = 0;

	int blamxloc = mtrarr[_iblk].blamx;

	int ibs;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = mtrarr[_iblk].list[ilist];
		int blai = sprndr[isup+1]-sprndr[isup];
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jsup = mtrarr[_iblk].ja[j];
			int jsupnew=jsup;
			if (jsupnew < 0) jsupnew = -jsupnew;
			int blaj = sprndr[jsupnew+1]-sprndr[jsupnew];
			if ((jsupnew>=_isupbeg && jsupnew <=_isupend) || jsupnew>=_isupbegrest) {
				if (blamxloc == 1) {
					ibs = j;
				} else {
					ibs = mtrarr[_iblk].bsa[j];
				};
				for (int kii=0;kii<blai*blaj;kii++) aloc[nzaloc+kii] = mtrarr[_iblk].a[ibs+kii];
				nzaloc += blai*blaj;
			} else {
				if (blamxloc == 1) {
					ibs = j;
				} else {
					ibs = mtrarr[_iblk].bsa[j];
				};
				int ibsd = _fct.bsdia[isup];
				int jbsd = _fct.bsdia[jsupnew];
				for (int kii=0;kii<blai;kii++) {
					for (int kjj=0;kjj<blaj;kjj++) {
						double auxloc = mtrarr[_iblk].a[ibs+kjj*blai+kii];
						if (auxloc < 0.0e0) auxloc = -auxloc;
						auxloc *= _theta;
						_fct.diagg[ibsd+kii*blai+kii] += auxloc;
						_fct.diagg[jbsd+kjj*blaj+kjj] += auxloc;
					};
				};
			};
		};
	};

// Free previous data of the block

	delete [] mtrarr[_iblk].list;
	delete [] mtrarr[_iblk].sprndr;
	delete [] mtrarr[_iblk].sprndc;
	delete [] mtrarr[_iblk].ia;
	delete [] mtrarr[_iblk].ja;
	delete [] mtrarr[_iblk].bsa;
	delete [] mtrarr[_iblk].a;

// Register new data

	mtrarr[_iblk].list   = listnew;
	mtrarr[_iblk].sprndr = sprndrloc;
	mtrarr[_iblk].sprndc = sprndcloc;
	mtrarr[_iblk].ia     = ialoc;
	mtrarr[_iblk].ja     = jaloc;
	mtrarr[_iblk].bsa    = bsaloc;
	mtrarr[_iblk].a      = aloc;

	mtrarr[_iblk].nsupc  = nlistnew;
	mtrarr[_iblk].nsupr  = nlistnew;
	mtrarr[_iblk].nlist  = nlistnew;
	mtrarr[_iblk].nzja   = nzjaloc;
	mtrarr[_iblk].nza    = nzaloc;
	mtrarr[_iblk].nzatot = nzaloc;

	if (mtrarr[_iblk].blamx == 1) mtrarr[_iblk].CleanBsa ();

// Free work memory

	delete [] listcount;

};

// Author: Kharchenko S.A.
// CFctR: Init current block row
//========================================================================================
void CFctR::Ilu2InitBlkRow (int _iblkrow, // Init current block row
							double _dshift,
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, 
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru) {

//	const char *funcname = "Ilu2InitBlkRow";

	double dzero = 0.0e0;

// Check the data on entry

	if (_iblkrow >= _gmtra.nblksc || _iblkrow < 0) {
		cout << " Wrong block row number " << endl;
		throw " Wrong block row number ";
	};

	int ilist, isup, isupbeg, isupend, nlistl;

	isupbeg = _gmtra.blksr[_iblkrow];
	isupend = _gmtra.blksr[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	if (nlistl != _mtral.nlist || nlistl != _mtrau.nlist) {
		cout << " Wrong number of supernode rows in the block row" << endl;
		throw " Wrong number of supernode rows in the block row";
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = _mtral.list[ilist];
		if (isup < isupbeg || isup > isupend) {
			cout << " Wrong supernode row in the block row" << endl;
			throw " Wrong supernode row in the block row";
		};
		isup = _mtrau.list[ilist];
		if (isup < isupbeg || isup > isupend) {
			cout << " Wrong supernode row in the block row" << endl;
			throw " Wrong supernode row in the block row";
		};
	};

// Copy the data

	_mtrl = _mtral;
	_mtru = _mtrau;

	_mtrl.nsupr = _mtrl.nlist;
	_mtrl.nsupc = _mtrl.nlist;
	_mtru.nsupr = _mtru.nlist;
	_mtru.nsupc = _mtru.nlist;

// Scale 

	int ibs, ibssc, jbssc, blai, blaj, blaij, blaj_2;
	int kii, kjj, kkk, kjjb, ibsv, jbsv;

	int blamxloc = _mtrl.blamx;

	double *paloc, *pscloc, *pvloc;
	double aux;

	for (ilist=0;ilist<nlistl;ilist++) {

		isup = _mtrl.list[ilist];
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];

		for (int j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {

			int jsup = _mtrl.ja[j];

			if (blamxloc == 1) {
				ibs = j;
			} else {
				ibs = _mtrl.bsa[j];
			};
			ibssc = bsdia[isup];
			jbssc = bsdia[jsup];
			blaj = _gmtra.sprndr[jsup+1]-_gmtra.sprndr[jsup];
			ibsv = adrfmask[isup];
			jbsv = adrfmask[jsup];

			blaij = blai * blaj;
			blaj_2 = blaj*blaj;

// Local point scaling (U part)

			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					aux = _mtru.a[ibs+kjj*blai+kii] * sclpt[ibsv+kii];
					vloc[kjj*blai+kii] = aux * sclpt[jbsv+kjj];
				};
			};

// Local supernode scaling (U part)

			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					paloc = vloc+kjjb;
					pscloc = scll+ibssc+kii;
					aux = dzero;
					for (kkk=0;kkk<blai;kkk++) {
//						aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//						aux += *paloc++ * *pscloc;
						aux += *paloc++ * *pscloc;
						pscloc += blai;
					};
					aloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					aux = dzero;
					paloc = aloc+kii;
					pscloc = sclu+jbssc+kjj;
					for (kkk=0;kkk<blaj;kkk++) {
//						aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//						aux += *paloc * *pscloc;
						aux += *paloc * *pscloc;
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

// Local point scaling (L part)

			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					aux = _mtrl.a[ibs+kjj*blai+kii] * sclpt[ibsv+kii];
					vloc[kjj*blai+kii] = aux * sclpt[jbsv+kjj];
				};
			};

// Local supernode scaling (L part)

			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					paloc = vloc+kjjb;
					pscloc = sclu+ibssc+kii;
					aux = dzero;
					for (kkk=0;kkk<blai;kkk++) {
//						aux += a[ibs+kjj*blai+kkk] * _sclu[ibssc+kkk*blai+kii];
//						aux += *paloc++ * *pscloc;
						aux += *paloc++ * *pscloc;
						pscloc += blai;
					};
					aloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			for (kii=0;kii<blai;kii++) {
				kjjb = 0;
				for (kjj=0;kjj<blaj;kjj++) {
					aux = dzero;
					paloc = aloc+kii;
					pscloc = scll+jbssc+kjj;
					for (kkk=0;kkk<blaj;kkk++) {
//						aux += _aloc[kkk*blai+kii] * _scll[jbssc+kkk*blaj+kjj];
//						aux += *paloc * *pscloc;
						aux += *paloc * *pscloc;
						pscloc += blaj;
						paloc += blai;
					};
					vloc[kjjb+kii] = aux;
					kjjb += blai;
				};
			};
			paloc = _mtrl.a+ibs;
			pvloc = vloc;
			for (kii=0;kii<blaij;kii++) *paloc++ = *pvloc++;

			ops += 2*blai*(blaij+blaj_2);

// Shift diagonal data if necessary

			paloc = _mtrl.a+ibs;

			if (isup == jsup) {
				for (kii=0;kii<blai;kii++) paloc[kii*blai+kii] += _dshift;
			};

		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Init current supernode row
//========================================================================================
void CFctR::Ilu2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau) {

//	const char *funcname = "Ilu2InitRow";

	double dzero = 0.0e0;

	int kii, ibs;
	int blamxloc, blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2;

	double *paloc;
	double *pfmask;
	double *pvloc;

	icycle++;

	blamxloc = _mtral.blamx;
	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	blai_2 = blai*blai;
	ibsd = bsdia[_irow];
	ibsv = adrfmask[_irow];

	nlist = 0;

	lstloc[nlist] = _irow;
	imask[_irow] = icycle;

	pfmask = fmaskl+ibsv*blai;
	pvloc = diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = czero;
	for (kii=0;kii<blai_2;kii++) *pfmask++ = dzero;

	int j = _mtrau.ia[_ilist];
	if (blamxloc == 1) {
		ibs = j;
	} else {
		ibs = _mtrau.bsa[j];
	};
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

	for (j=_mtral.ia[_ilist]+1;j<_mtral.ia[_ilist+1];j++) {
		int jj = _mtral.ja[j];
		lstloc[nlist] = jj;
		imask[jj] = icycle;

// Local copy

		if (blamxloc == 1) {
			ibs = j;
		} else {
			ibs = _mtral.bsa[j];
		};
		blaj = _gmtra.sprndr[jj+1]-_gmtra.sprndr[jj];
		jbsv = adrfmask[jj];
		blaij = blai * blaj;

		pfmask = fmasku+jbsv*blai;
		paloc = _mtrau.a+ibs;
		for (kii=0;kii<blaij;kii++) *pfmask++ = *paloc++;

		pfmask = fmaskl+jbsv*blai;
		paloc = _mtral.a+ibs;
		for (kii=0;kii<blaij;kii++) *pfmask++ = *paloc++;

		nlist++;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init current supernode row
//========================================================================================
void CFctR::Ilu2PointInitRow (char _type, int _ilist, int _irow, // Init current supernode row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau) {

//	const char *funcname = "Ilu2PointInitRow";

	double dzero = 0.0e0;

	int ibs;
	int ibsd, ibsv, jbsv;

	double *paloc;
	double *pfmask;
	double *pvloc;

	icycle++;

	ibsd = bsdia[_irow];
	ibsv = adrfmask[_irow];

	nlist = 0;

	lstloc[nlist] = _irow;
	imask[_irow] = icycle;

	pfmask = fmaskl+ibsv;
	pvloc = diagg+ibsd;
	*pfmask++ = dzero;

	int j = _mtrau.ia[_ilist];
	ibs = j;
	pfmask = fmasku+ibsv;
	pvloc = diagg+ibsd;
	paloc = _mtrau.a+ibs;
	if (_type == 'F') {
		*pfmask = *paloc + *pvloc;
	} else {
		*pfmask = *paloc;
	};

	nlist++;

	for (j=_mtral.ia[_ilist]+1;j<_mtral.ia[_ilist+1];j++) {
		int jj = _mtral.ja[j];
		lstloc[nlist] = jj;
		imask[jj] = icycle;

// Local copy

//		ibs = _mtral.bsa[j];
		jbsv = adrfmask[jj];

//		pfmask = fmasku+jbsv;
//		paloc = _mtrau.a+ibs;
//		*pfmask = *paloc;
		fmasku[jbsv] = _mtrau.a[j];

//		pfmask = fmaskl+jbsv;
//		paloc = _mtral.a+ibs;
//		*pfmask = *paloc;
		fmaskl[jbsv] = _mtral.a[j];

		nlist++;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2InitMask";

	double dzero = 0.0e0;

	int kii;
	int blai, blaj, jbsv, blaij;
	double *pfmask;

	int *addrgloc;
	double *glloc;
	double *guloc;
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
			pfmask = fmaskl+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = dzero;
			pfmask = fmasku+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = dzero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		addrgloc = addrg[_ilstprv];
		glloc = gl[_ilstprv];
		guloc = gu[_ilstprv];

		adrupdl[nupd] = glloc+addrgloc[kloc];
		adrupdu[nupd] = guloc+addrgloc[kloc];
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2PointInitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2PointInitMask";

	double dzero = 0.0e0;

	int jbsv;
//	double *pfmask;

	int *addrgloc;
	double *glloc;
	double *guloc;
	int k, kloc;
	int *jgloc;
	int jcolmn;

	nupd = 0;

	for (k=iv[_ilstprv];k<ig[_ilstprv+1];k++) {

		kloc = k-ig[_ilstprv];
		
		jgloc = jg[_ilstprv];
		jcolmn = jgloc[kloc];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (imask[jcolmn] != icycle) {
			imask[jcolmn] = icycle;
			lstloc[nlist] = jcolmn;
			jbsv = adrfmask[jcolmn];
//			pfmask = fmaskl+jbsv;
//			*pfmask++ = dzero;
			fmaskl[jbsv] = dzero;
//			pfmask = fmasku+jbsv;
//			*pfmask++ = dzero;
			fmasku[jbsv] = dzero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		addrgloc = addrg[_ilstprv];
		glloc = gl[_ilstprv];
		guloc = gu[_ilstprv];

		adrupdl[nupd] = glloc+addrgloc[kloc];
		adrupdu[nupd] = guloc+addrgloc[kloc];
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctR: Update current row
//========================================================================================
void CFctR::Ilu2UpdateRow (int _irow, int _irwprv, // Update current row
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2UpdateRow";

	int kii, kjj, ibs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv, blaik;

	double *pfmloc;
	double *plloc;
	double *pgloc;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	blak = _gmtra.sprndr[_irwprv+1]-_gmtra.sprndr[_irwprv];
	blaik = blai*blak;

	for (int k=0;k<nupd;k++) {
		int icol = lstupd[k];
		blaj = _gmtra.sprndr[icol+1]-_gmtra.sprndr[icol];
		jbsv = adrfmask[icol];
		ibs = jbsv*blai;

		double *elemg = adrupdu[k];
		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = fmasku+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = lelem+kikb;
				pgloc = elemg+kjkb;
				double *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
//					*pfmloc -= *plloc * *pgloc;
					(*pfmloc) -= *plloc * *pgloc;
				};
				pfmloc += blai;
				kjkb += blak;
			};
			kikb += blak;
		};

		elemg = adrupdl[k];
		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = fmaskl+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = uelem+kikb;
				pgloc = elemg+kjkb;
				double *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
//					*pfmloc -= *plloc * *pgloc;
					(*pfmloc) -= *plloc * *pgloc;
				};
				pfmloc += blai;
				kjkb += blak;
			};
			kikb += blak;
		};

		ops += 2 * blaik * blaj;

	};

};

// Author: Kharchenko S.A.
// CFctR: Update current row
//========================================================================================
void CFctR::Ilu2PointUpdateRow (int _irow, int _irwprv, // Update current row
							const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2UpdateRow";

//	int ibs;
	int jbsv;

//	double *pfmloc;
//	double *plloc;
//	double *pgloc;
	double *elemg; 
	int k, icol;

	for (k=0;k<nupd;k++) {
		icol = lstupd[k];
		jbsv = adrfmask[icol];

		elemg = adrupdu[k];
//		pfmloc = fmasku+ibs;
//		plloc = lelem;
//		pgloc = elemg;
//		(*pfmloc) -= *plloc * *pgloc;
		fmasku[jbsv] -= *lelem * *elemg;

		elemg = adrupdl[k];
//		pfmloc = fmaskl+ibs;
//		plloc = uelem;
//		pgloc = elemg;
//		(*pfmloc) -= *plloc * *pgloc;
		fmaskl[jbsv] -= *uelem * *elemg;

	};

	ops += 2*nupd;

};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend) { // Update arrays that support transposed structures search

//	const char *funcname = "Ilu2UpdateTransp";

	int ilstprv = ibegm[_ilist];

	while (ilstprv != -1) {

		int ilstprv1 = madj[ilstprv];

		if (iv[ilstprv] >= ig[ilstprv+1]-1) {

			madj[ilstprv] = -1;

		} else {

			int j = iv[ilstprv]+1;
			int *jgloc = jg[ilstprv];
			int jloc = j-ig[ilstprv];
			int icolmn = jgloc[jloc];
			if (icolmn < 0) icolmn = -icolmn;
			if (icolmn > _isupend) {
				madj[ilstprv] = -1;
			} else {
				int icolmnl=icolmn-_isupbeg;
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
void CFctR::Ilu2FiltrRow (int _ilist, int _irow, // Perform filtering of the current row
									int _fcttyp, double _theta, double _tau2, 
									int _ndiagfct2, int _isupbrd, double _tau2brd,
									const CGSMatrix &_gmtra, const CSMatrix &_mtraini) {

//	const char *funcname = "Ilu2FiltrRow";

	int nlist2 = 0;
	int kii, kjj, blai, blaj, ibsv, jbsv, jbsd, blaij;
	int ibeg1, iend1, ibeg2, iend2, ip1, ip2, jj1, jj2;
	double *pfmaskl, *pfmasku, *pfmaskuu, *pdiag;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	ibsv = adrfmask[_irow];

	if (_fcttyp == 0) {

		if (nlist > 0) qsort (lstloc+1, nlist-1, sizeof(int), compint);

		ibeg1 = 1;
		iend1 = nlist-1;
		ibeg2 = _mtraini.ia[_ilist]+1;
		iend2 = _mtraini.ia[_ilist+1]-1;
		ip1 = ibeg1;
		ip2 = ibeg2;
		while (ip1 <= iend1 && ip2 <= iend2) {
			jj1 = lstloc[ip1];
			jj2 = _mtraini.ja[ip2];
			if (jj1 == jj2) {
				lstupd[nlist2] = jj1;
				nlist2++;
				ip1++;
				ip2++;
			} else if (jj1 < jj2) {

				int j = lstloc[ip1];

				blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
				jbsv = adrfmask[j];
				jbsd = bsdia[j];
				blaij = blai*blaj;

				double aux = 0.0e0;

				pfmaskl = fmaskl+jbsv*blai;
				for (kii=0;kii<blaij;kii++) {
//					double auxloc = _fmaskl[jbsv*blai+kii];
					double auxloc = *pfmaskl;
					if (auxloc < 0.0e0) auxloc = -auxloc;
					pfmaskl++;
					if (auxloc > aux) aux = auxloc;
				};
				pfmasku = fmasku+jbsv*blai;
				for (kii=0;kii<blaij;kii++) {
//					double auxloc = _fmasku[jbsv*blai+kii];
					double auxloc = *pfmasku;
					if (auxloc < 0.0e0) auxloc = -auxloc;
					pfmasku++;
					if (auxloc > aux) aux = auxloc;
				};

				pfmaskuu = fmasku+ibsv*blai;
				for (kii=0;kii<blai;kii++) {
					pfmaskl = fmaskl+jbsv*blai+kii;
					pfmasku = fmasku+jbsv*blai+kii;
					pdiag = diagg+jbsd;
					for (kjj=0;kjj<blaj;kjj++) {
//						double auxloc = _fmaskl[jbsv*blai+kjj*blai+kii];
						double auxloc = *pfmaskl;
						if (auxloc < 0.0e0) auxloc = -auxloc;
//						double auxloc1 = _fmasku[jbsv*blai+kjj*blai+kii];
						double auxloc1 = *pfmasku;
						if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
						if (auxloc1 > auxloc) auxloc = auxloc1;
						auxloc *= _theta;
//						_fmasku[ibsv*blai+kii*blai+kii] += auxloc;
						*pfmaskuu += auxloc;
//						_diagg[jbsd+kjj*blaj+kjj] += auxloc;
						*pdiag += auxloc;
						pfmaskl += blai;
						pfmasku += blai;
						pdiag += blaj+1;
					};
					pfmaskuu += blai+1;
				};

				ip1++;
			} else if (jj1 > jj2) {
				ip2++;
			};
		};

	} else {

		for (int i=1;i<nlist;i++) {

			int j = lstloc[i];

			blaj = _gmtra.sprndr[j+1]-_gmtra.sprndr[j];
			jbsv = adrfmask[j];
			jbsd = bsdia[j];
			blaij = blai*blaj;

			double aux = 0.0e0;

			pfmaskl = fmaskl+jbsv*blai;
			for (kii=0;kii<blaij;kii++) {
//				double auxloc = _fmaskl[jbsv*blai+kii];
				double auxloc = *pfmaskl;
				if (auxloc < 0.0e0) auxloc = -auxloc;
				pfmaskl++;
				if (auxloc > aux) aux = auxloc;
			};
			pfmasku = fmasku+jbsv*blai;
			for (kii=0;kii<blaij;kii++) {
//				double auxloc = _fmasku[jbsv*blai+kii];
				double auxloc = *pfmasku;
				if (auxloc < 0.0e0) auxloc = -auxloc;
				pfmasku++;
				if (auxloc > aux) aux = auxloc;
			};

			if ((j < _isupbrd && aux >= _tau2) || (j >= _isupbrd && aux >= _tau2brd)) {
				lstupd[nlist2] = j;
				nlist2++;
			} else {
				pfmaskuu = fmasku+ibsv*blai;
				for (kii=0;kii<blai;kii++) {
					pfmaskl = fmaskl+jbsv*blai+kii;
					pfmasku = fmasku+jbsv*blai+kii;
					pdiag = diagg+jbsd;
					for (kjj=0;kjj<blaj;kjj++) {
//						double auxloc = _fmaskl[jbsv*blai+kjj*blai+kii];
						double auxloc = *pfmaskl;
						if (auxloc < 0.0e0) auxloc = -auxloc;
//						double auxloc1 = _fmasku[jbsv*blai+kjj*blai+kii];
						double auxloc1 = *pfmasku;
						if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
						if (auxloc1 > auxloc) auxloc = auxloc1;
						auxloc *= _theta;
//						_fmasku[ibsv*blai+kii*blai+kii] += auxloc;
						*pfmaskuu += auxloc;
//						_diagg[jbsd+kjj*blaj+kjj] += auxloc;
						*pdiag += auxloc;
						pfmaskl += blai;
						pfmasku += blai;
						pdiag += blaj+1;
					};
					pfmaskuu += blai+1;
				};
			};

		};

	};

	if (_fcttyp == 3 || _fcttyp == 4) throw " CFctR::Ilu2FiltrRow: diagonal elements limitation case is not implemented ";

	nlist = nlist2;

};

// Author: Kharchenko S.A.
// CFctR: Perform filtering of the current row
//========================================================================================
void CFctR::Ilu2PointFiltrRow (int _ilist, int _irow, // Perform filtering of the current row
											int _fcttyp, double _theta, double _tau2, 
											int _ndiagfct2, int _isupbrd, double _tau2brd,
											const CGSMatrix &_gmtra, const CSMatrix &_mtraini) {

//	const char *funcname = "Ilu2PointFiltrRow";

	int nlist2 = 0;
	int ibsv, jbsv, jbsd;
	int i, j;
	int ibeg1, iend1, ibeg2, iend2, ip1, ip2, jj1, jj2;
	double aux, auxloc, auxloc1;
//	double *pfmaskl, *pfmasku, *pfmaskuu, *pdiag;

	ibsv = adrfmask[_irow];

	if (_fcttyp == 0) {

		if (nlist > 0) qsort (lstloc+1, nlist-1, sizeof(int), compint);

		ibeg1 = 1;
		iend1 = nlist-1;
		ibeg2 = _mtraini.ia[_ilist]+1;
		iend2 = _mtraini.ia[_ilist+1]-1;
		ip1 = ibeg1;
		ip2 = ibeg2;
		while (ip1 <= iend1 && ip2 <= iend2) {
			jj1 = lstloc[ip1];
			jj2 = _mtraini.ja[ip2];
			if (jj1 == jj2) {
				lstupd[nlist2] = jj1;
				nlist2++;
				ip1++;
				ip2++;
			} else if (jj1 < jj2) {

				j = lstloc[ip1];

				jbsv = adrfmask[j];
				jbsd = bsdia[j];

				aux = 0.0e0;

				auxloc = fmaskl[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;
				auxloc = fmasku[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				auxloc = fmaskl[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc1 = fmasku[jbsv];
				if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
				if (auxloc1 > auxloc) auxloc = auxloc1;
				auxloc *= _theta;
				fmasku[ibsv] += auxloc;
				diagg[jbsd] += auxloc;

				ip1++;
			} else if (jj1 > jj2) {
				ip2++;
			};
		};

	} else {

		for (i=1;i<nlist;i++) {

			j = lstloc[i];

			jbsv = adrfmask[j];
			jbsd = bsdia[j];

			aux = 0.0e0;

			auxloc = fmaskl[jbsv];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsv];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;

			if ((j < _isupbrd && aux >= _tau2) || (j >= _isupbrd && aux >= _tau2brd)) {
				lstupd[nlist2] = j;
				nlist2++;
			} else {
				auxloc = fmaskl[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc1 = fmasku[jbsv];
				if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
				if (auxloc1 > auxloc) auxloc = auxloc1;
				auxloc *= _theta;
				fmasku[ibsv] += auxloc;
				diagg[jbsd] += auxloc;
			};

		};

	};

	if ((_fcttyp == 3 || _fcttyp == 4) && nlist2 > _ndiagfct2) {
		int nlist3 = 0;
		int nlistsrt = 0;
		for (i=0;i<nlist2;i++) {
			j = lstupd[i];
			jbsv = adrfmask[j];
			jbsd = bsdia[j];
			if (j == _irow || j >= _isupbrd) {
				lstloc[nlist3] = j;
				nlist3++;
			} else {
				aux = 0.0e0;
				auxloc = fmaskl[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;
				auxloc = fmasku[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				if (auxloc > aux) aux = auxloc;

				diarr[nlistsrt].dvalue = aux;
				diarr[nlistsrt].intvalue = j;
				nlistsrt++;
			};
		};

		if (nlistsrt > _ndiagfct2) {

			if (nlistsrt != 0) qsort (diarr, nlistsrt, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

			for (i=0;i<nlistsrt-_ndiagfct2;i++) {
				j = diarr[i].intvalue;
				jbsv = adrfmask[j];
				jbsd = bsdia[j];
				auxloc = fmaskl[jbsv];
				if (auxloc < 0.0e0) auxloc = -auxloc;
				auxloc1 = fmasku[jbsv];
				if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
				if (auxloc1 > auxloc) auxloc = auxloc1;
				auxloc *= _theta;
				fmasku[ibsv] += auxloc;
				diagg[jbsd] += auxloc;
			};
			for (i=nlistsrt-_ndiagfct2;i<nlistsrt;i++) {
				lstloc[nlist3] = diarr[i].intvalue;
				nlist3++;
			};

			nlist2 = nlist3;
			for (i=0;i<nlist2;i++) lstupd[i] = lstloc[i];

		};
	};

	nlist = nlist2;

};

// Author: Kharchenko S.A.
// CFctR: Compute diagonal pivot and scale current row
//========================================================================================
void CFctR::Ilu2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
								const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2PivScaleRow";

	double czero = 0.0e0;

	int kii, kjj, kkk, ibs;
	int blai, blaj, ibsv, jbsv;
	double caux;

// Copy diagonal supernode

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	ibsv = adrfmask[_irow];
	ibs = ibsv * blai;

	for (kii=0;kii<blai*blai;kii++) aloc[kii] = fmasku[ibs+kii];

// Compute singular value decomposition of the diagonal supernode

	int blal = blai;
	int info = 0;

	dgesvd_ ("A", "A", &blal, &blal,
				aloc, &blal, eig, 
				uloc, &blal, vloc, &blal, 
				work, &lwork, &info);

	if (info != 0) throw " Error in the Lapack routine DGESVD";

	ops += 8 * blai * blai * blai;

//	ofstream fout ("CheckPivIlu.dat");
//	fout << " irow = " << _irow << endl;
//	OutArr(fout," Eig ",blal,eig);

// Modify singular values if necessary

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
	double *pdeloc;
	double *pfmloc;
	double *palloc;

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = deleml+kii;
		palloc = uloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_deleml[kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//			*pdeloc = eig[kii] * conj(*palloc++);
			double temp = *palloc++;
			*pdeloc = eig[kii] * temp;
			pdeloc += blai;
		};
		kiib += blai;
	};

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = delemu+kii;
		palloc = vloc+kii;
//		palloc = _vloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_delemu[kjj*blai+kii] = _vloc[kjj*blai+kii] * _eig[kii];
//			*pdeloc = eig[kii] * conj(*palloc);
			double temp = *palloc;
			*pdeloc = eig[kii] * temp;
			pdeloc += blai;
			palloc += blai;
//			palloc++;
		};
		kiib += blai;
	};

	ops += 2 * blai * blai;

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
				pfmloc = fmaskl+ibs+kjjb;
				pdeloc = delemu+kii;
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _delemu[kkk*blai+kii] * _fmaskl[ibs+kjj*blai+kkk];
//					caux += *pdeloc * *pfmloc++;
					caux += *pdeloc * *pfmloc++;
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = caux;
				palloc += blai;
			};
		};

		pfmloc = fmaskl+ibs;
		palloc = aloc;
//		for (kii=0;kii<blaij;kii++) _fmaskl[ibs+kii] = _aloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmloc++ = *palloc++;

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
					caux += *pdeloc * *pfmloc++;
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

		ops += 2 * blai * blaij;

	};

};

// Author: Kharchenko S.A.
// CFctR: Compute diagonal pivot and scale current row
//========================================================================================
void CFctR::Ilu2PointPivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
								const CGSMatrix &_gmtra) {

//	const char *funcname = "Ilu2PointPivScaleRow";

	int ibs;
	int ibsv, jbsv;

// Copy diagonal supernode

	ibsv = adrfmask[_irow];
	ibs = ibsv;

	aloc[0] = fmasku[ibs];

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
	dpiv[ibsv] = eig[0];
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

	int i, j;

	for (i=0;i<nlist;i++) {

		j = lstupd[i];
		jbsv = adrfmask[j];

		ibs = jbsv;

		fmaskl[ibs] = *delemu * fmaskl[ibs];

		fmasku[ibs] = *deleml * fmasku[ibs];

	};

	ops += 2*nlist;

};

// Author: Kharchenko S.A.
// CFctR: Store current row
//========================================================================================
void CFctR::Ilu2StoreRow (int _ilist, int _irow, // Store current row
									int _fcttyp, double _tau1, int _ndiagfct, 
									int _isupbrd, double _tau1brd,
									int &_nzrtot, 
									const CGSMatrix &_gmtra) {

	const char *funcname = "Ilu2StoreRow";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	double *pg, *pfmask, *pd;

	if (_fcttyp == 3 || _fcttyp == 4) throw " CFctR::Ilu2StoreRow: diagonal elements limitation case is not implemented ";

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
	double *glloc;
	double *guloc;

	jgloc = new int [nlist+1];
	if (!jgloc) MemoryFail (funcname);
	adrgloc = new int [nlist+1];
	if (!adrgloc) MemoryFail (funcname);
	glloc = new double [nz];
	if (!glloc) MemoryFail (funcname);
	guloc = new double [nz];
	if (!guloc) MemoryFail (funcname);

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gl[_ilist] = glloc;
	gu[_ilist] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

// Store current diagonal supernode

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];
	jgloc[nzjgloc] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = deleml+kiib;
		pg = glloc+nzloc+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gl [_nznew+kjj*blai+kii] = _deleml[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
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
		pfmask = fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = *pfmask;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			pfmask++;
			if (auxloc > aux) aux = auxloc;
		};
		pfmask = fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = *pfmask;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			pfmask++;
			if (auxloc > aux) aux = auxloc;
		};
		adrgloc[nzjgloc] = nzloc;
		if (_fcttyp == 0 || (j < _isupbrd && aux >= _tau1) || (j >= _isupbrd && aux >= _tau1brd)) {
			jgloc[nzjgloc] = j;
		} else {
			jgloc[nzjgloc] = -j;
			_nzrtot += blaij;
		};
		pg = glloc+nzloc;
		pfmask = fmaskl+jbsv*blai;
//		for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
		for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
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
// CFctR: Store current row
//========================================================================================
void CFctR::Ilu2PointStoreRow (int _ilist, int _irow, // Store current row
											int _fcttyp, double _tau1, int _ndiagfct, 
											int _isupbrd, double _tau1brd, 
											int &_nzrtot, 
											const CGSMatrix &_gmtra) {

	const char *funcname = "Ilu2PointStoreRow";

	int jbsv;
	double *pg, *pd;

// Count the number of elements to be stored

	int nz = 1;

	int i, j;
	for (i=0;i<nlist;i++) {
		j = lstupd[i];
		nz += 1;
	};

// Allocate the memory

	int *jgloc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		adrgloc = new int [nzalloc];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [nzalloc];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc]+nzused;
		glloc   = glblk[iblkalloc]+nzused;
		guloc   = gublk[iblkalloc]+nzused;
	};

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gl[_ilist] = glloc;
	gu[_ilist] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

// Store current diagonal supernode

	jgloc[nzjgloc] = _irow;
	pd = deleml;
	pg = glloc+nzloc;
	*pg = *pd;
	pd = delemu;
	pg = guloc+nzloc;
	*pg = *pd;
	adrgloc[nzjgloc] = nzloc;
	nzjgloc++;
	nzloc += 1;

	double aux, auxloc;

	if ((_fcttyp != 3 && _fcttyp != 4) || nlist < _ndiagfct) {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jbsv = adrfmask[j];
			aux = 0.0e0;
			auxloc = fmaskl[jbsv];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsv];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			adrgloc[nzjgloc] = nzloc;
			if (_fcttyp == 0 || (j < _isupbrd && aux >= _tau1) || (j >= _isupbrd && aux >= _tau1brd)) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j;
				_nzrtot += 1;
			};
			glloc[nzloc] = fmaskl[jbsv];
			guloc[nzloc] = fmasku[jbsv];

			nzjgloc++;
			nzloc += 1;
		};

	} else {

		for (i=0;i<nlist;i++) {
			j = lstupd[i];
			jbsv = adrfmask[j];
			aux = 0.0e0;
			auxloc = fmaskl[jbsv];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
			auxloc = fmasku[jbsv];
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
			iord = lstloc[i];
			i1 = nlist-iord-1;
			jbsv = adrfmask[j];
			adrgloc[nzjgloc] = nzloc;
			if (i1 < _ndiagfct) {
				jgloc[nzjgloc] = j;
			} else {
				jgloc[nzjgloc] = -j;
				_nzrtot += 1;
			};
			glloc[nzloc] = fmaskl[jbsv];
			guloc[nzloc] = fmasku[jbsv];

			nzjgloc++;
			nzloc += 1;

		};

	};

	ig[_ilist+1] = ig[_ilist]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
								const CGSMatrix &_gmtra, 
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
								CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2StoreScaleU";

	double dzero = 0.0e0;

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	double caux;

	double *pscloc;
	double *pgloc;
	double *paloc;

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

	CSMatrixRS mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixRS temp (nlistl,nzjlu,nzlu);

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

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	double *glloc;
	double *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jj >= 0) {
				ibs = adrgloc[jloc];
				for (kii=0;kii<blai*blaj;kii++) _mtrl.a[nzlu+kii] = glloc[ibs+kii];
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
			for (int j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {
				jbs = _mtrl.bsa[j];
				int jj = _mtrl.ja[j];
				blaj = _gmtra.sprndr[jj+1]-_gmtra.sprndr[jj];
				blaij = blai*blaj;
				jjbssc = bsdia[jj];
				if (jj == isup) {

					kiib = 0;
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = aloc+kii;
						for (kjj=0;kjj<blai;kjj++) {
							pgloc = _mtrl.a+jbs+kjjb;
							pscloc = scll+ibssc+kiib;
							caux = dzero;
							for (kkk=0;kkk<blai;kkk++) {
//								caux += *pscloc * *pgloc;
								caux += *pscloc * *pgloc;
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
					pgloc = _mtrl.a+jbs;
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;

					kiib = 0;
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = aloc+kii;
						for (kjj=0;kjj<blai;kjj++) {
							pgloc = _mtru.a+jbs+kjjb;
							pscloc = sclu+ibssc+kiib;
							caux = dzero;
							for (kkk=0;kkk<blai;kkk++) {
//								caux += *pscloc * *pgloc;
								caux += *pscloc * *pgloc;
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
							caux = dzero;
							pgloc = _mtrl.a+jbs+kii;
							pscloc = sclinvl+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								caux += *pgloc * *pscloc;
								caux += *pgloc * *pscloc;
								pgloc += blai;
								pscloc += blaj;
							};
							*paloc = caux;
							paloc += blai;
							kjjb += blai;
						};
					};
					paloc = aloc;
					pgloc = _mtrl.a+jbs;
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;

					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = aloc+kii;
						for (kjj=0;kjj<blaj;kjj++) {
							caux = dzero;
							pgloc = _mtru.a+jbs+kii;
							pscloc = sclinvu+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								caux += *pgloc * *pscloc;
								caux += *pgloc * *pscloc;
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

					ops += 2*blaij*blaj;
				};
			};
		};

	};

// Allocate and init second order matrices

	_mtrl2 = mtrdummy;
	_mtru2 = mtrdummy;

	CSMatrixRS temp2 (nlistl,nzjlu2,nzlu2);

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

	_mtrl2 = temp2;
	_mtru2 = temp2;

	temp2 = mtrdummy;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			if (jjloc > isupend) {
				ibs = adrgloc[jloc];
				for (kii=0;kii<blai*blaj;kii++) _mtrl2.a[nzlu2+kii] = glloc[ibs+kii];
				for (kii=0;kii<blai*blaj;kii++) _mtru2.a[nzlu2+kii] = guloc[ibs+kii];
				nzjlu2++;
				nzlu2 += blai*blaj;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2PointStoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
								const CGSMatrix &_gmtra, 
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
								CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2PointStoreScaleU";

	int ibssc, jbs, jjbssc;

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
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jj >= 0) {
				nzjlu++;
				nzlu++;
			};
		};
	};

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jjloc > isupend) {
				nzjlu2++;
				nzlu2++;
			};
		};
	};

// Allocate and init first order matrices

	CSMatrixRS mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixRS temp (nlistl,nzjlu,nzlu);

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
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jj >= 0) {
				temp.ja[nzjlu] = jjloc;
				temp.bsa[nzjlu] = nzlu;
				nzjlu++;
				nzlu++;
			};
		};
		temp.list[ilist] = isup;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nzja = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	temp.CleanBsa ();

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	double *glloc;
	double *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jj >= 0) {
				ibs = adrgloc[jloc];
				_mtrl.a[nzlu] = glloc[ibs];
				_mtru.a[nzlu] = guloc[ibs];
				nzjlu++;
				nzlu++;
			};
		};
	};

// Scale

	if (_param.sclpttype < 2) {

		for (ilist=0;ilist<nlistl;ilist++) {
			isup = isupbeg+ilist;
			ibssc = bsdia[isup];
			for (int j=_mtrl.ia[ilist];j<_mtrl.ia[ilist+1];j++) {
				jbs = j;
				jj = _mtrl.ja[j];
				jjbssc = bsdia[jj];
				if (jj == isup) {

					_mtrl.a[jbs] = _mtrl.a[jbs] * scll[ibssc];
					_mtru.a[jbs] = _mtru.a[jbs] * sclu[ibssc];

	 				ops += 2;

				} else {

					_mtrl.a[jbs] = _mtrl.a[jbs] * sclinvl[jjbssc];
					_mtru.a[jbs] = _mtru.a[jbs] * sclinvu[jjbssc];

					ops += 2;
				};
			};
		};

	};

// Allocate and init second order matrices

	_mtrl2 = mtrdummy;
	_mtru2 = mtrdummy;

	CSMatrixRS temp2 (nlistl,nzjlu2,nzlu2);

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
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jjloc > isupend) {
				temp2.ja[nzjlu2] = jj;
				temp2.bsa[nzjlu2] = nzlu2;
				nzjlu2++;
				nzlu2++;
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

	temp2.CleanBsa ();

	_mtrl2 = temp2;
	_mtru2 = temp2;

	temp2 = mtrdummy;

	nzjlu2 = 0;
	nzlu2 = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			if (jjloc > isupend) {
				ibs = adrgloc[jloc];
				_mtrl2.a[nzlu2] = glloc[ibs];
				_mtru2.a[nzlu2] = guloc[ibs];
				nzjlu2++;
				nzlu2++;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform factorization of the current block row
//========================================================================================
void CFctR::Ilu2FctBlkRow (int _iblkrow, const CSlvParam &_param, int _isupbrd, // Perform factorization of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, const CSMatrixRS &_mtraini,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
							CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2FctBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int nzrtot;
	int irwprv, ilstprv, icolmn, j, jloc, ibs;
	int isupbeg, isupend, nlistl, blai, blaj, kii;
	bool elmisr;
	int *jgloc;
	int *adrgloc;
	double *glloc, *guloc;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	nzrtot = 0;

	for (ilist=0;ilist<=_mtral.nlist;ilist++) ibegm[ilist] = -1;

	ig[0] = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];

// Init current scaled row

		Ilu2InitRow ('F', ilist, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			jgloc = jg[ilstprv];
			adrgloc = addrg[ilstprv];
			glloc = gl[ilstprv];
			guloc = gu[ilstprv];

			j = iv[ilstprv];

			jloc = j-ig[ilstprv];

			icolmn = jgloc[jloc];
			ibs = adrgloc[jloc];

			irwprv = _mtral.list[ilstprv];
			int blailoc = _gmtra.sprndr[irwprv+1]-_gmtra.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _gmtra.sprndr[icolmnl+1]-_gmtra.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) lelem[kii] = glloc[ibs+kii];
			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = guloc[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2InitMask (elmisr, irow, ilstprv, _gmtra);

// Perform updates

			Ilu2UpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend);

// Perform filtering of the row with diagonal modification if necessary 

		Ilu2FiltrRow (ilist, irow, 
								_param.fcttyp, _param.theta, _param.tau2, _param.ndiagfct2, 
								_isupbrd, _param.tau2bord, 
								_gmtra, _mtraini);

// Compute and store the pivot value and scale current supernode row

		Ilu2PivScaleRow (irow, _param.pivmin,
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

		Ilu2StoreRow (ilist, irow, 
							_param.fcttyp, _param.tau1, _param.ndiagfct2, 
							_isupbrd, _param.tau1bord, 
							nzrtot,
							_gmtra);

	};

// Store the resulting U objects

	Ilu2StoreScaleU (_iblkrow, _param, 
						_gmtra,
						_mtrl, _mtru,
						_mtrl2, _mtru2);

	nzr += nzrtot;

// Store the hystogram of values

	int nloc = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {
		irow = _mtral.list[ilist];
		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];
		nloc += blai;
	};

	double *dpivarr;

	dpivarr = new double [nloc];

	nloc = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {
		irow = _mtral.list[ilist];
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
		delete [] gl[ilist];
		delete [] gu[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform factorization of the current block row
//========================================================================================
void CFctR::Ilu2PointFctBlkRow (int _iblkrow, const CSlvParam &_param, int _isupbrd, // Perform factorization of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, const CSMatrixRS &_mtraini,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
							CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2PointFctBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int nzrtot;
	int irwprv, ilstprv, icolmn, j, jloc, ibs;
	int isupbeg, isupend, nlistl, blai, blaj;
	bool elmisr;
	int *jgloc;
	int *adrgloc;
	double *glloc, *guloc;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

	nzrtot = 0;

	for (ilist=0;ilist<=_mtral.nlist;ilist++) ibegm[ilist] = -1;

	ig[0] = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

// Init current scaled row

		Ilu2PointInitRow ('F', ilist, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			jgloc = jg[ilstprv];
			adrgloc = addrg[ilstprv];
			glloc = gl[ilstprv];
			guloc = gu[ilstprv];

			j = iv[ilstprv];

			jloc = j-ig[ilstprv];

			icolmn = jgloc[jloc];
			ibs = adrgloc[jloc];

			irwprv = _mtral.list[ilstprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = 1;

			lelem[0] = glloc[ibs];
			uelem[0] = guloc[ibs];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2PointInitMask (elmisr, irow, ilstprv, _gmtra);

// Perform updates

			Ilu2PointUpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend);

// Perform filtering of the row with diagonal modification if necessary 

		Ilu2PointFiltrRow (ilist, irow, 
									_param.fcttyp, _param.theta, _param.tau2, _param.ndiagfct2,
									_isupbrd, _param.tau2bord, 
									_gmtra, _mtraini);

// Compute and store the pivot value and scale current supernode row

		Ilu2PointPivScaleRow (irow, _param.pivmin,
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

		Ilu2PointStoreRow (ilist, irow, 
									_param.fcttyp, _param.tau1, _param.ndiagfct, 
									_isupbrd, _param.tau1bord, 
									nzrtot,
									_gmtra);

	};

// Store the resulting U objects

	Ilu2PointStoreScaleU (_iblkrow, _param, 
						_gmtra,
						_mtrl, _mtru,
						_mtrl2, _mtru2);

	nzr += nzrtot;

// Store the hystogram of values

	int nloc = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {
		irow = _mtral.list[ilist];
		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];
		nloc += blai;
	};

	double *dpivarr;

	dpivarr = new double [nloc];

	nloc = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {
		irow = _mtral.list[ilist];
		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];
		int ibsv = adrfmask[irow];;
		dpivarr[nloc] = dpiv[ibsv];
		nloc += blai;
	};

	hyst.UpdateHystogram (nloc, dpivarr);

	delete [] dpivarr;

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

/*
	for (ilist=0;ilist<nlistl;ilist++) {
		delete [] jg[ilist];
		delete [] addrg[ilist];
		delete [] gl[ilist];
		delete [] gu[ilist];
	};
*/
};

// Author: Kharchenko S.A.
// CFctR: Init arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2InitTransp (int _isupbeg, int _isupend, // Init arrays that support transposed structures search
							const CSMatrixRS &_mtrl2) { 

//	const char *funcname = "Ilu2InitTransp";

	int ilist, isup, j, jj, nlistl, ilistn;

	for (ilist=0;ilist<_mtrl2.nlist;ilist++) {
		isup = _mtrl2.list[ilist];
		iv[ilist] = _mtrl2.ia[ilist+1];
		for (j=_mtrl2.ia[ilist];j<_mtrl2.ia[ilist+1];j++) {
			jj = _mtrl2.ja[j];
			if (jj < 0) jj = -jj;
			if (jj >= _isupbeg) {
				iv[ilist] = j;
				break;
			};
		};
	};

	nlistl = _isupend-_isupbeg+1;

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
			if (jj < 0) jj = -jj;
			if (jj >= _isupbeg && jj <= _isupend) {
				ilistn = jj-_isupbeg;
				madj[ilist] = ibegm[ilistn];
				ibegm[ilistn] = ilist;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra, 
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2InitMask";

	double dzero = 0.0e0;

	int kii;
	int blai, blaj, jbsv, blaij, ibs;
	double *pfmask;

	int k;
	int jcolmn;

	nupd = 0;

	blai = _gmtra.sprndr[_irow+1]-_gmtra.sprndr[_irow];

	for (k=iv[_ilstprv];k<_mtrl2.ia[_ilstprv+1];k++) {

		jcolmn = _mtrl2.ja[k];
		ibs = _mtrl2.bsa[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (imask[jcolmn] != icycle) {
			imask[jcolmn] = icycle;
			lstloc[nlist] = jcolmn;
			blaj = _gmtra.sprndr[jcolmn+1]-_gmtra.sprndr[jcolmn];
			jbsv = adrfmask[jcolmn];
			blaij = blai*blaj;
			pfmask = fmaskl+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = dzero;
			pfmask = fmasku+jbsv*blai;
			for (kii=0;kii<blaij;kii++) *pfmask++ = dzero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		ibs = _mtrl2.bsa[k];

		adrupdl[nupd] = _mtrl2.a+ibs;
		adrupdu[nupd] = _mtru2.a+ibs;
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctR: Init mask and update arrays
//========================================================================================
void CFctR::Ilu2PointInitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra, 
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2) {

//	const char *funcname = "Ilu2PointInitMask";

	double dzero = 0.0e0;

	int jbsv, ibs;

	int k;
	int jcolmn;

	nupd = 0;

	for (k=iv[_ilstprv];k<_mtrl2.ia[_ilstprv+1];k++) {

		jcolmn = _mtrl2.ja[k];
		ibs = k;

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (imask[jcolmn] != icycle) {
			imask[jcolmn] = icycle;
			lstloc[nlist] = jcolmn;
			jbsv = adrfmask[jcolmn];
			fmaskl[jbsv] = dzero;
			fmasku[jbsv] = dzero;
			nlist++;
		};
		lstupd[nupd] = jcolmn;
		ibs = k;

		adrupdl[nupd] = _mtrl2.a+ibs;
		adrupdu[nupd] = _mtru2.a+ibs;
		nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CFctR: Update arrays that support transposed structures search
//========================================================================================
void CFctR::Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend, // Update arrays that support transposed structures search
								const CSMatrixRS &_mtrl2) {

//	const char *funcname = "Ilu2UpdateTransp";

	int ilstprv = ibegm[_ilist];

	while (ilstprv != -1) {

		int ilstprv1 = madj[ilstprv];

		if (iv[ilstprv] >= _mtrl2.ia[ilstprv+1]-1) {

			madj[ilstprv] = -1;

		} else {

			int j = iv[ilstprv]+1;
			int icolmn = _mtrl2.ja[j];
			if (icolmn < 0) icolmn = -icolmn;
			if (icolmn > _isupend) {
				madj[ilstprv] = -1;
			} else {
				int icolmnl=icolmn-_isupbeg;
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
void CFctR::Ilu2StoreRow (int _ilist, int _irow, // Store current row
							const CGSMatrix &_gmtra) {

	const char *funcname = "Ilu2StoreRow";

	int kii;
	int blai, blaj, jbsv, blaij;
	double *pg, *pfmask;

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
	double *glloc;
	double *guloc;

	jgloc = new int [nlist+1];
	if (!jgloc) MemoryFail (funcname);
	adrgloc = new int [nlist+1];
	if (!adrgloc) MemoryFail (funcname);
	glloc = new double [nz];
	if (!glloc) MemoryFail (funcname);
	guloc = new double [nz];
	if (!guloc) MemoryFail (funcname);

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gl[_ilist] = glloc;
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
		pg = glloc+nzloc;
		pfmask = fmaskl+jbsv*blai;
//		for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
		for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
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
// CFctR: Store current row
//========================================================================================
void CFctR::Ilu2PointStoreRow (int _ilist, int _irow, // Store current row
							const CGSMatrix &_gmtra) {

	const char *funcname = "Ilu2PointStoreRow";

	int jbsv;
	double *pg, *pfmask;

// Allocate the memory

	int *jgloc;
	int *adrgloc;
	double *glloc;
	double *guloc;

	if (nlist+1 > nzalloc-nzused) {
		nzalloc = nzallocmax;
		if (nzalloc < nlist+1) nzalloc = nlist+1;
		jgloc = new int [nzalloc];
		if (!jgloc) MemoryFail (funcname);
		adrgloc = new int [nzalloc];
		if (!adrgloc) MemoryFail (funcname);
		glloc = new double [nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [nzalloc];
		if (!guloc) MemoryFail (funcname);
		iblkalloc++;
		jgblk[iblkalloc]   = jgloc;
		adrgblk[iblkalloc] = adrgloc;
		glblk[iblkalloc]   = glloc;
		gublk[iblkalloc]   = guloc;
		nzused = 0;
	} else {
		jgloc   = jgblk[iblkalloc]+nzused;
		adrgloc = adrgblk[iblkalloc]+nzused;
		glloc   = glblk[iblkalloc]+nzused;
		guloc   = gublk[iblkalloc]+nzused;
	};

	jg[_ilist] = jgloc;
	addrg[_ilist] = adrgloc;
	gl[_ilist] = glloc;
	gu[_ilist] = guloc;

	int nzjgloc = 0;
	int nzloc = 0;

	int i, j;

	for (i=0;i<nlist;i++) {
		j = lstloc[i];
		jbsv = adrfmask[j];
		adrgloc[nzjgloc] = nzloc;
		jgloc[nzjgloc] = j;
		pg = glloc+nzloc;
		pfmask = fmaskl+jbsv;
		*pg = *pfmask;
		pg = guloc+nzloc;
		pfmask = fmasku+jbsv;
		*pg = *pfmask;

		nzjgloc++;
		nzloc++;
	};

	ig[_ilist+1] = ig[_ilist]+nzjgloc;

	nzused += nzjgloc;

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2StoreU (int _iblkrow, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru) {

//	const char *funcname = "Ilu2StoreU";

//	double dzero = 0.0e0;

	int kii;
	int blai, blaj;

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

	CSMatrixRS mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixRS temp (nlistl,nzjlu,nzlu);

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

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	double *glloc;
	double *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = _gmtra.sprndr[jjloc+1]-_gmtra.sprndr[jjloc];
			ibs = adrgloc[jloc];
			for (kii=0;kii<blai*blaj;kii++) _mtrl.a[nzlu+kii] = glloc[ibs+kii];
			for (kii=0;kii<blai*blaj;kii++) _mtru.a[nzlu+kii] = guloc[ibs+kii];
			nzjlu++;
			nzlu += blai*blaj;
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Store the resulting U objects
//========================================================================================
void CFctR::Ilu2PointStoreU (int _iblkrow, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru) {

//	const char *funcname = "Ilu2PointStoreU";

//	double dzero = 0.0e0;

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
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			nzjlu++;
			nzlu++;
		};
	};

// Allocate and init matrices

	CSMatrixRS mtrdummy;

	_mtrl = mtrdummy;
	_mtru = mtrdummy;

	CSMatrixRS temp (nlistl,nzjlu,nzlu);

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
		jgloc = jg[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			temp.ja[nzjlu] = jjloc;
			temp.bsa[nzjlu] = nzlu;
			nzjlu++;
			nzlu++;
		};
		temp.list[ilist] = isup;
		temp.ia[ilist+1] = nzjlu;
	};

	temp.nzja = nzjlu;
	temp.nza = nzlu;
	temp.nzatot = nzlu;

	if (blamx == 1) temp.CleanBsa ();

	_mtrl = temp;
	_mtru = temp;

	temp = mtrdummy;

	int *adrgloc;
	double *glloc;
	double *guloc;

	nzjlu = 0;
	nzlu = 0;

	for (ilist=0;ilist<nlistl;ilist++) {
		isup = isupbeg+ilist;
		jgloc = jg[ilist];
		adrgloc = addrg[ilist];
		glloc = gl[ilist];
		guloc = gu[ilist];
		for (j=ig[ilist];j<ig[ilist+1];j++) {
			jloc = j-ig[ilist];
			jj = jgloc[jloc];
			jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			ibs = adrgloc[jloc];
			_mtrl.a[nzlu] = glloc[ibs];
			_mtru.a[nzlu] = guloc[ibs];
			nzjlu++;
			nzlu++;
		};
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row
//========================================================================================
void CFctR::Ilu2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru) {

//	const char *funcname = "Ilu2UpdateBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int irwprv, ilstprv, icolmn, j, ibs;
	int isupbeg, isupend, nlistl, blai, blaj, kii;
	bool elmisr;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

// Init arrays that support transposed structures search

	Ilu2InitTransp (isupbeg, isupend, _mtrl2);

	ig[0] = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];

// Init current row

		Ilu2InitRow ('U', ilist, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			j = iv[ilstprv];

			icolmn = _mtrl2.ja[j];
			ibs = _mtrl2.bsa[j];

			irwprv = _mtrl2.list[ilstprv];
			int blailoc = _gmtra.sprndr[irwprv+1]-_gmtra.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _gmtra.sprndr[icolmnl+1]-_gmtra.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) lelem[kii] = _mtrl2.a[ibs+kii];
			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = _mtru2.a[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2InitMask (elmisr, irow, ilstprv, _gmtra,
								_mtrl2, _mtru2);

// Perform updates

			Ilu2UpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend, _mtrl2);

// Sort the list of indices

		if (nlist != 0) qsort (lstloc, nlist, sizeof(int), compint);

// Store the data in the global arrays

		Ilu2StoreRow (ilist, irow,
						_gmtra);

	};

// Store the resulting U objects

	Ilu2StoreU (_iblkrow, 
					_gmtra,
					_mtrl, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<nlistl;ilist++) {
		delete [] jg[ilist];
		delete [] addrg[ilist];
		delete [] gl[ilist];
		delete [] gu[ilist];
	};

};

// Author: Kharchenko S.A.
// CFctR: Perform updates of the current block row
//========================================================================================
void CFctR::Ilu2PointUpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru) {

//	const char *funcname = "Ilu2PointUpdateBlkRow";

// Main cycle over the rows

	int ilist, irow;
	int irwprv, ilstprv, icolmn, j, ibs;
	int isupbeg, isupend, nlistl, blai, blaj, kii;
	bool elmisr;

	int nzaloc = _mtral.nza;
	nzallocmax = nzaloc * 3 / 10;
	if (nzallocmax < 1) nzallocmax = 1;

	nzalloc = 0;
	iblkalloc = -1;
	nzused = 0;

	isupbeg = _gmtra.blksc[_iblkrow];
	isupend = _gmtra.blksc[_iblkrow+1]-1;

	nlistl = isupend-isupbeg+1;

// Init arrays that support transposed structures search

	Ilu2InitTransp (isupbeg, isupend, _mtrl2);

	ig[0] = 0;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {

		irow = _mtral.list[ilist];

		blai = _gmtra.sprndr[irow+1]-_gmtra.sprndr[irow];

// Init current row

		Ilu2InitRow ('U', ilist, irow, _gmtra, _mtral, _mtrau);

// Update current row by the previous ones

		ilstprv = ibegm[ilist];

		while (ilstprv != -1) {

			j = iv[ilstprv];

			icolmn = _mtrl2.ja[j];
			ibs = j;
//			ibs = _mtrl2.bsa[j];

			irwprv = _mtrl2.list[ilstprv];
			int blailoc = _gmtra.sprndr[irwprv+1]-_gmtra.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _gmtra.sprndr[icolmnl+1]-_gmtra.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) lelem[kii] = _mtrl2.a[ibs+kii];
			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = _mtru2.a[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ilu2PointInitMask (elmisr, irow, ilstprv, _gmtra,
									_mtrl2, _mtru2);

// Perform updates

			Ilu2PointUpdateRow (irow, irwprv, _gmtra);

			ilstprv = madj[ilstprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ilu2UpdateTransp (ilist, isupbeg, isupend, _mtrl2);

// Sort the list of indices

		if (nlist != 0) qsort (lstloc, nlist, sizeof(int), compint);

// Store the data in the global arrays

		Ilu2PointStoreRow (ilist, irow,
						_gmtra);

	};

// Store the resulting U objects

	Ilu2PointStoreU (_iblkrow, 
					_gmtra,
					_mtrl, _mtru);

// Free local factorization arrays 

	for (ilist=0;ilist<iblkalloc+1;ilist++) {
		delete [] jgblk[ilist];
		delete [] adrgblk[ilist];
		delete [] glblk[ilist];
		delete [] gublk[ilist];
	};

//	for (ilist=0;ilist<nlistl;ilist++) {
//		delete [] jg[ilist];
//		delete [] addrg[ilist];
//		delete [] gl[ilist];
//		delete [] gu[ilist];
//	};

};

// Author: Kharchenko S.A.
// CFctR: Add two block rows/columns assuming that the list of indices coincide
//========================================================================================
CSMatrixRS CFctR::AddBlocks (const CGSMatrix &_gmtra, // Add two block rows/columns assuming that the list of indices coincide
								CSMatrixRS &_mtru, CSMatrixRS &_mtru2) {

	const char *funcname = "AddBlocks";

//	double dzero = 0.0e0;

// Check that the lists of indices coincide

	if (_mtru.nlist != _mtru2.nlist) {
		throw " The lists of indices do not coincide";
	};
	int i;
	for (i=0;i<_mtru.nlist;i++) {
		if (_mtru.list[i] != _mtru2.list[i]) {
			throw " The lists of indices do not coincide";
		};
	};

// Count the total number of nonzero supernodes

	int *ialoc;

	ialoc = new int [_mtru.nlist+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	int nzjaloc = 0;

	int ilist;
	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		icycle++;
		int j;
		for (j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {
			int jsup = _mtru.ja[j];
			if (imask[jsup] != icycle) {
				imask[jsup] = icycle;
				nzjaloc++;
			};
		};
		for (j=_mtru2.ia[ilist];j<_mtru2.ia[ilist+1];j++) {
			int jsup = _mtru2.ja[j];
			if (imask[jsup] != icycle) {
				imask[jsup] = icycle;
				nzjaloc++;
			};
		};
		ialoc[ilist+1] = nzjaloc;
	};

// Allocate and init array jaloc

	int *jaloc;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);

	nzjaloc = 0;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		icycle++;
		int j;
		for (j=_mtru.ia[ilist];j<_mtru.ia[ilist+1];j++) {
			int jsup = _mtru.ja[j];
			if (imask[jsup] != icycle) {
				imask[jsup] = icycle;
				jaloc[nzjaloc] = jsup;
				nzjaloc++;
			};
		};
		for (j=_mtru2.ia[ilist];j<_mtru2.ia[ilist+1];j++) {
			int jsup = _mtru2.ja[j];
			if (imask[jsup] != icycle) {
				imask[jsup] = icycle;
				jaloc[nzjaloc] = jsup;
				nzjaloc++;
			};
		};
	};

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		int nz = ialoc[ilist+1]-ialoc[ilist];
		int ibs = ialoc[ilist];
		qsort (jaloc+ibs, nz, sizeof(int), compint);
	};

// Compute array bsa

	int *bsaloc;

	bsaloc = new int [nzjaloc];
	if (!bsaloc) MemoryFail (funcname);

	int nzaloc = 0;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		int isup = _mtru.list[ilist];
		int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		for (int j=ialoc[ilist];j<ialoc[ilist+1];j++) {
			int jsup = jaloc[j];
			int blaj = _gmtra.sprndr[jsup+1]-_gmtra.sprndr[jsup];
			bsaloc[j] = nzaloc;
			nzaloc += blai*blaj;
		};
	};

// Store computed data

	CSMatrixRS temp (_mtru.nlist,nzjaloc,nzaloc);

	for (i=0;i< _mtru.nlist;i++) temp.list[i] = _mtru.list[i];
	for (i=0;i<=_mtru.nlist;i++) temp.sprndr[i] = 0;
	for (i=0;i<=_mtru.nlist;i++) temp.sprndc[i] = 0;
	for (i=0;i<=_mtru.nlist;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];
	for (i=0;i<nzjaloc;i++) temp.bsa[i] = bsaloc[i];

// Free work memory

	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;

// Compute array a

	nzaloc = 0;
	nzjaloc = 0;

	int blamx1 = _mtru.blamx;
	int blamx2 = _mtru2.blamx;
	int ibs1, ibs2;

	for (ilist=0;ilist<_mtru.nlist;ilist++) {
		int isup = _mtru.list[ilist];
		int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		int ibeg1 = _mtru.ia[ilist];
		int iend1 = _mtru.ia[ilist+1]-1;
		int ibeg2 = _mtru2.ia[ilist];
		int iend2 = _mtru2.ia[ilist+1]-1;
		int ind1 = ibeg1;
		int ind2 = ibeg2;
		while (ind1 <= iend1 || ind2 <= iend2) {
			if (ind1 <= iend1 && ind2 <= iend2) {
				int jsup1 = _mtru.ja[ind1];
				int jsup2 = _mtru2.ja[ind2];
				if (jsup1 == jsup2) {
					if (blamx1 == 1) {
						ibs1 = ind1;
					} else {
						ibs1 = _mtru.bsa[ind1];
					};
					if (blamx2 == 1) {
						ibs2 = ind2;
					} else {
						ibs2 = _mtru2.bsa[ind2];
					};
					int blaj = _gmtra.sprndr[jsup1+1]-_gmtra.sprndr[jsup1];
					for (int kii=0;kii<blai*blaj;kii++) temp.a[nzaloc+kii] = _mtru.a[ibs1+kii]+_mtru2.a[ibs2+kii];
					nzjaloc++;
					nzaloc += blai*blaj;
					ind1++;
					ind2++;
				} else if (jsup1 < jsup2) {
					if (blamx1 == 1) {
						ibs1 = ind1;
					} else {
						ibs1 = _mtru.bsa[ind1];
					};
					int blaj = _gmtra.sprndr[jsup1+1]-_gmtra.sprndr[jsup1];
					for (int kii=0;kii<blai*blaj;kii++) temp.a[nzaloc+kii] = _mtru.a[ibs1+kii];
					nzjaloc++;
					nzaloc += blai*blaj;
					ind1++;
				} else if (jsup2 < jsup1) {
					if (blamx2 == 1) {
						ibs2 = ind2;
					} else {
						ibs2 = _mtru2.bsa[ind2];
					};
					int blaj = _gmtra.sprndr[jsup2+1]-_gmtra.sprndr[jsup2];
					for (int kii=0;kii<blai*blaj;kii++) temp.a[nzaloc+kii] = _mtru2.a[ibs2+kii];
					nzjaloc++;
					nzaloc += blai*blaj;
					ind2++;
				};
			} else if (ind1 <= iend1) {
				int jsup1 = _mtru.ja[ind1];
				if (blamx1 == 1) {
					ibs1 = ind1;
				} else {
					ibs1 = _mtru.bsa[ind1];
				};
				int blaj = _gmtra.sprndr[jsup1+1]-_gmtra.sprndr[jsup1];
				for (int kii=0;kii<blai*blaj;kii++) temp.a[nzaloc+kii] = _mtru.a[ibs1+kii];
				nzjaloc++;
				nzaloc += blai*blaj;
				ind1++;
			} else if (ind2 <= iend2) {
				int jsup2 = _mtru2.ja[ind2];
				if (blamx2 == 1) {
					ibs2 = ind2;
				} else {
					ibs2 = _mtru2.bsa[ind2];
				};
				int blaj = _gmtra.sprndr[jsup2+1]-_gmtra.sprndr[jsup2];
				for (int kii=0;kii<blai*blaj;kii++) temp.a[nzaloc+kii] = _mtru2.a[ibs2+kii];
				nzjaloc++;
				nzaloc += blai*blaj;
				ind2++;
			};
		};
	};

// Fill control data

	temp.m      = _mtru.m;
	temp.n      = _mtru.n;
	temp.nsupr  = _mtru.nlist;
	temp.nsupc  = _mtru.nlist;
	temp.nlist  = _mtru.nlist;
	temp.nzja   = nzjaloc;
	temp.nza    = nzaloc;
	temp.nzatot = nzaloc;
	temp.blamx  = _mtru.blamx;

	if (temp.blamx == 1) temp.CleanBsa ();

	return temp;

};

// Author: Kharchenko S.A.
// CFctR: Init zero data diagonal block
//========================================================================================
CSMatrixRS CFctR::CreateBlock (int _isupbeg, int _isupend, const CGSMatrix &_gmtra) { // Init zero data diagonal block

//	const char *funcname = "CreateBlock";

	double dzero = 0.0e0;

// Count the total number of nonzero elements

	int nzaloc = 0;

	int isup;
	for (isup=_isupbeg;isup<=_isupend;isup++) {
		int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		nzaloc += blai*blai;
	};

	int nlistloc = _isupend-_isupbeg+1;

// Store the data

	CSMatrixRS temp (nlistloc,nlistloc,nzaloc);

	int i;
	for (i=0;i<nlistloc;i++) temp.list[i] = _isupbeg+i;
	for (i=0;i<=nlistloc;i++) temp.sprndr[i] = 0;
	for (i=0;i<=nlistloc;i++) temp.sprndc[i] = 0;
	for (i=0;i<=nlistloc;i++) temp.ia[i] = i;
	for (i=0;i<nlistloc;i++) temp.ja[i] = _isupbeg+i;

	nzaloc = 0;

	for (isup=_isupbeg;isup<=_isupend;isup++) {
		int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
		temp.bsa[isup-_isupbeg] = nzaloc;
		nzaloc += blai*blai;
	};

	for (i=0;i<nzaloc;i++) temp.a[i] = dzero;

// Fill control data

	temp.m      = _gmtra.m;
	temp.n      = _gmtra.n;
	temp.nsupr  = nlistloc;
	temp.nsupc  = nlistloc;
	temp.nlist  = nlistloc;
	temp.nzja   = nlistloc;
	temp.nza    = nzaloc;
	temp.nzatot = nzaloc;
	temp.blamx  = blamx;

	return temp;

};
