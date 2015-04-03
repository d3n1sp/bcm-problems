//------------------------------------------------------------------------------------------------
// File: gsmatrixc.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "globals.h"
#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CGSMatrixCS: Memory allocation zero data constructor
//========================================================================================
CGSMatrixCS::CGSMatrixCS (): CGSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixCS_00";

	mtrarr = new CSMatrixCS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

	ivect = -1;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Memory allocation zero data constructor
//========================================================================================
CGSMatrixCS::CGSMatrixCS (int _nblks, int _nsups, int _nlist): CGSMatrix (_nblks, _nsups, _nlist) { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixCS_01";

	mtrarr = new CSMatrixCS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

	ivect = -1;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Destructor
//========================================================================================
CGSMatrixCS::~CGSMatrixCS () { // Destructor
//	std::cout << " On entry to CGSMatrixCS destructor " << std::endl;
	delete [] mtrarr;
	ivect = -1;
	CSVectorC vectdummy;
	vector1 = vectdummy;
	vector2 = vectdummy;
//	std::cout << " On return from CGSMatrixCS destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CGSMatrixCS: Equality operator
//========================================================================================
CGSMatrixCS &CGSMatrixCS::operator= (const CGSMatrixCS &_aa) { // Equality operator

	const char *funcname = "CGSMatrixCS_=";

	delete [] mtrarr;

	this->CGSMatrix::operator= ((const CGSMatrix &) _aa);

	mtrarr = new CSMatrixCS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);
//	for (int iblk=0;iblk<nblksr;iblk++) {
//		mtrarr[iblk] = _aa.mtrarr[iblk];
//	};

	ivect = _aa.ivect;
	vector1 = _aa.vector1;
	vector2 = _aa.vector2;

	return *this;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Filter block row/column according to the index
//========================================================================================
void CGSMatrixCS::FilterBlock (int _iblk, int _isupend) { // Filter block row/column according to the index

	const char *funcname = "FilterBlock";

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
			if (jj>_isupend) nz++;
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
			if (jsupnew>_isupend) {
				jaloc[nzjaloc] = jsup;
				bsaloc[nzjaloc] = nzaloc;
				nzjaloc++;
				nzaloc += blai*blaj;
			};
		};
	};

// Create a

	dcmplx *aloc;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	nzaloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = mtrarr[_iblk].list[ilist];
		int blai = sprndr[isup+1]-sprndr[isup];
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jsup = mtrarr[_iblk].ja[j];
			int jsupnew=jsup;
			if (jsupnew < 0) jsupnew = -jsupnew;
			int blaj = sprndr[jsupnew+1]-sprndr[jsupnew];
			if (jsupnew>_isupend) {
				int ibs = mtrarr[_iblk].bsa[j];
				for (int kii=0;kii<blai*blaj;kii++) aloc[nzaloc+kii] = mtrarr[_iblk].a[ibs+kii];
				nzaloc += blai*blaj;
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

// Free work memory

	delete [] listcount;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Write block row/column into the disk
//========================================================================================
void CGSMatrixCS::WriteBlock (int _iblk) { // Write block row/column into the disk

//	const char *funcname = "WriteBlock";

	int nzaloc;
	int ibsloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	ibsloc = ibsfile[ifile];

	FPut (files[ifile], nzaloc, aloc, ibsloc);

	bl2ibsfile[_iblk] = ibsloc;
	bl2file[_iblk] = ifile;

	ibsloc += nzaloc;

	ibsfile[ifile] = ibsloc;

	ifile++;

	if (ifile >= nfiles) ifile = 0;

Exit:;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Rewrite block row/column into the disk
//========================================================================================
void CGSMatrixCS::ReWriteBlock (int _iblk) { // Rewrite block row/column into the disk

//	const char *funcname = "ReWriteBlock";

	int nzaloc;
	int ibsloc, ifileloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FPut (files[ifileloc], nzaloc, aloc, ibsloc);

Exit:;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Read block row/column from the disk into the main memory if necessary
//========================================================================================
void CGSMatrixCS::ReadBlock (int _iblk) { // Read block row/column from the disk into the main memory if necessary

	const char *funcname = "ReadBlock";

	int nzaloc;
	int ibsloc;
	int ifileloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FGet (files[ifileloc], nzaloc, aloc, ibsloc);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Allocate block row/column
//========================================================================================
void CGSMatrixCS::AllocateBlock (int _iblk) { // Allocate block row/column

	const char *funcname = "AllocateBlock";

	int nzaloc;
	dcmplx *aloc;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Free block row/column from main mamory
//========================================================================================
void CGSMatrixCS::FreeBlock (int _iblk) { // Free block row/column from main mamory

	const char *funcname = "FreeBlock";

	int nzaloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	nzaloc = 0;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Create GSMatrixCS structure via SMatrixCS
//========================================================================================
void CGSMatrixCS::CS2GCS (const CSMatrixCS &_mtra, int _nblks, int *_blks) { // Create GSMatrixCS structure via SMatrixCS

	const char *funcname = "CS2GCS";

// Create temporary bl2cpu data

	int *bl2cpu;

	bl2cpu = new int [_nblks];
	if (!bl2cpu) MemoryFail (funcname);

	for (int iblk=0;iblk<_nblks;iblk++) bl2cpu[iblk] = 0;

	int myid = 0;

// Create the structure

	CS2GCS (_mtra, _nblks, _blks, myid, bl2cpu);

// Free work arrays

	delete [] bl2cpu;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Create GSMatrixCS structure via SMatrixCS for specified indices only
//========================================================================================
void CGSMatrixCS::CS2GCS (const CSMatrixCS &_mtra, int _nblks, int *_blks, // Create GSMatrixCS structure via SMatrixCS for specified indices only
							int _myid, int *_bl2cpu) {

//	const char *funcname = "CS2GCS";

// Create temporary CGSMatrixCS data

	CGSMatrixCS tempgcs (_nblks, _mtra.nsupr, _nblks);

	*this = tempgcs;

	m = _mtra.m;
	n = _mtra.n;
	nblksr = _nblks;
	nblksc = _nblks;
	nsupr = _mtra.nsupr;
	nsupc = _mtra.nsupc;

	int iblk, isup;

	for (iblk=0;iblk<=_nblks;iblk++) blksr[iblk] = _blks[iblk];
	for (iblk=0;iblk<=_nblks;iblk++) blksc[iblk] = _blks[iblk];
	for (isup=0;isup<=nsupr;isup++) sprndr[isup] = _mtra.sprndr[isup];
	for (isup=0;isup<=nsupc;isup++) sprndc[isup] = _mtra.sprndc[isup];

	bl2ndr[0] = 0;
	for (iblk=0;iblk<nblksr;iblk++) {
		int isupbeg=blksr[iblk];
		int isupend=blksr[iblk+1]-1;
		int ni=sprndr[isupend+1]-sprndr[isupbeg];
		bl2ndr[iblk+1] = bl2ndr[iblk] + ni;
	};

	bl2ndc[0] = 0;
	for (iblk=0;iblk<nblksc;iblk++) {
		int isupbeg=blksc[iblk];
		int isupend=blksc[iblk+1]-1;
		int ni=sprndc[isupend+1]-sprndc[isupbeg];
		bl2ndc[iblk+1] = bl2ndc[iblk] + ni;
	};

	nlist=0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_bl2cpu[iblk]==_myid) listb[nlist++] = iblk;
	};

// Create the set of CSMatrixCS block rows

	CSMatrixCS mtrdummy;

	int isupbeg, isupend, nsupl;
	int nzjal, nzal, blai, blaj, j, jsup, kii;

	for (iblk=0;iblk<_nblks;iblk++) {

		if (_bl2cpu[iblk] == _myid) {

			isupbeg = _blks[iblk];
			isupend = _blks[iblk+1]-1;

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup];j<_mtra.ia[isup+1];j++) {
					jsup = _mtra.ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					nzjal++;
					nzal += blai*blaj;
				};
			};

// Allocate the data

			CSMatrixCS temp (nsupl, nzjal, nzal);

// Init the matrix

			temp.m = _mtra.m;
			temp.n = _mtra.n;
			temp.nsupr = nsupl;
			temp.nsupc = nsupl;
			temp.nlist = nsupl;
			temp.nzja = nzjal;
			temp.nza = nzal;
			temp.nzatot = nzal;
			temp.blamx = _mtra.blamx;

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
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup];j<_mtra.ia[isup+1];j++) {
					jsup = _mtra.ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					temp.ja[nzjal] = jsup;
					temp.bsa[nzjal] = nzal;
					int ibs = _mtra.bsa[j];
					for (kii=0;kii<blai*blaj;kii++) temp.a[nzal+kii] = _mtra.a[ibs+kii];
					nzjal++;
					nzal += blai*blaj;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

// Store the matrix

			mtrarr[iblk] = temp;

			temp = mtrdummy;

		};

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Combine matrices into one
//========================================================================================
void CGSMatrixCS::CombineMatrices (int _nlist, int *_list, // Combine matrices into one
									const CSMatrixCS *_mtrarr,
									CSMatrixCS &_mtrnew) const {

//	const char *funcname = "CombineMatrices";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	int i, j, isup, jsup, blai, blaj, ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
		nzatotal += _mtrarr[iblk].nza;
	};

// Allocate and fill

	CSMatrixCS temp(nlisttotal,nzjatotal,nzatotal);

	int nlistloc=0, nzjaloc=0, nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nza;kii++) temp.a[nzaloc+kii] = _mtrarr[iblk].a[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
		nzaloc += _mtrarr[iblk].nza;
	};

// Assign bsa array

	nzaloc = 0;

	for (ilist=0;ilist<nlisttotal;ilist++) {
		isup = temp.list[ilist];
		blai = sprndc[isup+1]-sprndc[isup];
		for (j=temp.ia[ilist];j<temp.ia[ilist+1];j++) {
			jsup = temp.ja[j];
			if (jsup < 0) jsup = -jsup;
			blaj = sprndr[jsup+1]-sprndr[jsup];
			temp.bsa[j] = nzaloc;
			nzaloc += blai*blaj;
		};
	};

// Fill control data

	for (i=0;i<=nlisttotal;i++) temp.sprndr[i] = 0;
	for (i=0;i<=nlisttotal;i++) temp.sprndc[i] = 0;

	temp.m = m;
	temp.n = n;
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nzatot = nzatotal;
	temp.nza = nzatotal;
	iblk = _list[0];
	temp.blamx = _mtrarr[iblk].blamx;

// Return the result

	_mtrnew = temp;

	CSMatrixCS mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Split matrix into the set of matrices according to the list
//========================================================================================
void CGSMatrixCS::SplitMatrix (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
								const CSMatrixCS &_mtr, 
								CSMatrixCS *_mtrarr) {

//	const char *funcname = "SplitMatrix";

// Check the data on entry

	int i, j, isup, jsup, blai, blaj, ilist, ilistblk, iblk, kii;
	int nlistloc=0, nzjaloc=0, nzaloc=0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlistloc += blksc[iblk+1]-blksc[iblk];
	};

	if (nlistloc != _mtr.nlist) {
		cout << " CGSMatrixCS::SplitMatrix: Error in the list of indices when splitting the matrix" << endl;
		throw " CGSMatrixCS::SplitMatrix: Error in the list of indices when splitting the matrix ";
	};

// Compute the set of matrices

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	for (ilistblk=0;ilistblk<_nlist;ilistblk++) {

		iblk = _list[ilistblk];

// Determine the sizes of the data to be stored

		nlistloc = blksc[iblk+1]-blksc[iblk];

		nzjaloc = _mtr.ia[nlisttotal+nlistloc]-_mtr.ia[nlisttotal];

		nzaloc = 0;

		for (ilist=nlisttotal;ilist<nlisttotal+nlistloc;ilist++) {
			isup = _mtr.list[ilist];
			blai = sprndc[isup+1]-sprndc[isup];
			for (j=_mtr.ia[ilist];j<_mtr.ia[ilist+1];j++) {
				jsup = _mtr.ja[j];
				if (jsup < 0) jsup = -jsup;
				blaj = sprndr[jsup+1]-sprndr[jsup];
				nzaloc += blai*blaj;
			};
		};

// Allocate and fill

		CSMatrixCS temp(nlistloc,nzjaloc,nzaloc);

		temp.ia[0] = 0;

		for (kii=0;kii<nlistloc;kii++) temp.list[kii] = _mtr.list[nlisttotal+kii];
		for (kii=0;kii<nlistloc;kii++) {
			temp.ia[kii+1] = temp.ia[kii] + (_mtr.ia[nlisttotal+kii+1] - _mtr.ia[nlisttotal+kii]);
		};
		for (kii=0;kii<nzjaloc;kii++) temp.ja[kii] = _mtr.ja[nzjatotal+kii];
		for (kii=0;kii<nzaloc;kii++) temp.a[kii] = _mtr.a[nzatotal+kii];

		nlisttotal += nlistloc;
		nzjatotal += nzjaloc;
		nzatotal += nzaloc;

// Assign bsa array

		nzaloc = 0;

		for (ilist=0;ilist<nlistloc;ilist++) {
			isup = temp.list[ilist];
			blai = sprndc[isup+1]-sprndc[isup];
			for (j=temp.ia[ilist];j<temp.ia[ilist+1];j++) {
				jsup = temp.ja[j];
				if (jsup < 0) jsup = -jsup;
				blaj = sprndr[jsup+1]-sprndr[jsup];
				temp.bsa[j] = nzaloc;
				nzaloc += blai*blaj;
			};
		};

// Fill control data

		for (i=0;i<=nlistloc;i++) temp.sprndr[i] = 0;
		for (i=0;i<=nlistloc;i++) temp.sprndc[i] = 0;

		temp.m = m;
		temp.n = n;
		temp.nsupr = nlistloc;
		temp.nsupc = nlistloc;
		temp.nlist = nlistloc;
		temp.nzja = nzjaloc;
		temp.nzatot = nzaloc;
		temp.nza = nzaloc;
		temp.blamx = _mtr.blamx;

// Store the result

		_mtrarr[iblk] = temp;

		CSMatrixCS mtrdummy;

		temp = mtrdummy;

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Split matrix into the set of matrices according to the list
//========================================================================================
void CGSMatrixCS::SplitMatrixList (int &_nlistblk, int *_listblk, int *_sp2blk, // Split matrix into the set of matrices according to the list
												const CSMatrixCS &_mtr, 
												CSMatrixCS *_mtrarr) {

	const char *funcname = "SplitMatrixList";

// Compute the list of blocks

	int i;

	int *imaskblk;

	imaskblk = new int [nblksc];
	if (!imaskblk) MemoryFail (funcname);

	for (i=0;i<nblksc;i++) imaskblk[i] = -1;

	int icycle = -1;

	icycle++;

	_nlistblk = 0;

	int isup, iblk;

	for (i=0;i<_mtr.nlist;i++) {
		isup = _mtr.list[i];
		iblk = _sp2blk[isup];
		if (imaskblk[iblk] != icycle) {
			_listblk[_nlistblk] = iblk;
			_nlistblk++;
			imaskblk[iblk] = icycle;
		};
	};

	std::sort (_listblk, _listblk+_nlistblk);

// Split the matrix

	SplitMatrix (_nlistblk, _listblk,
						_mtr, _mtrarr);

// Free work arrays

	delete [] imaskblk;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Prepare the matrix for exchange
//========================================================================================
void CGSMatrixCS::PrepareSendMatrix (int _inode, const CTree &_tree, // Prepare the matrix for exchange
										CGSMatrixCS &_gmtru2, const CSMatrixCS *_mtru2charr,
										CSMatrixCS &_mtru2send) const {

//	const char *funcname = "PrepareSendMatrix";

// Determine control data

	int iblkbeg, iblkend, isupbeg, nchilds;

	iblkbeg = _tree.nodes[_inode].indbeg;
	iblkend = _tree.nodes[_inode].indend;
	nchilds = _tree.nodes[_inode].nchilds;

	if (nchilds == 1) nchilds = 0;

	isupbeg = blksr[iblkend+1];

// Count the numbers

	int nlistloc = 0;
	int nzjaloc = 0;
	int nzaloc = 0;

	int ichild;

	for (ichild=0;ichild<nchilds;ichild++) {
		for (int ilist=0;ilist<_mtru2charr[ichild].nlist;ilist++) {
			int isup = _mtru2charr[ichild].list[ilist];
			int blai = sprndr[isup+1]-sprndr[isup];
			int icount = 0;
			for (int j=_mtru2charr[ichild].ia[ilist];j<_mtru2charr[ichild].ia[ilist+1];j++) {
				int jsup = _mtru2charr[ichild].ja[j];
				if (jsup < 0) jsup = -jsup;
				int blaj = sprndr[jsup+1]-sprndr[jsup];
				if (jsup >= isupbeg) {
					if (icount==0) {
						nlistloc++;
						icount = 1;
					};
					nzjaloc++;
					nzaloc += blai*blaj;
				};
			};
		};
	};

	int iblk;

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {
		for (int ilist=0;ilist<_gmtru2.mtrarr[iblk].nlist;ilist++) {
			int isup = _gmtru2.mtrarr[iblk].list[ilist];
			int blai = sprndr[isup+1]-sprndr[isup];
			int icount = 0;
			for (int j=_gmtru2.mtrarr[iblk].ia[ilist];j<_gmtru2.mtrarr[iblk].ia[ilist+1];j++) {
				int jsup = _gmtru2.mtrarr[iblk].ja[j];
				if (jsup < 0) jsup = -jsup;
				int blaj = sprndr[jsup+1]-sprndr[jsup];
				if (jsup >= isupbeg) {
					if (icount==0) {
						nlistloc++;
						icount = 1;
					};
					nzjaloc++;
					nzaloc += blai*blaj;
				};
			};
		};
	};

// Allocate the data to be stored

	CSMatrixCS temp(nlistloc,nzjaloc,nzaloc);

// Store the data

	nlistloc = 0;
	nzjaloc = 0;
	nzaloc = 0;

	temp.ia[0] = 0;

	for (ichild=0;ichild<nchilds;ichild++) {
		for (int ilist=0;ilist<_mtru2charr[ichild].nlist;ilist++) {
			int isup = _mtru2charr[ichild].list[ilist];
			int blai = sprndr[isup+1]-sprndr[isup];
			int icount = 0;
			for (int j=_mtru2charr[ichild].ia[ilist];j<_mtru2charr[ichild].ia[ilist+1];j++) {
				int jsup = _mtru2charr[ichild].ja[j];
				int jsupini=jsup;
				if (jsup < 0) jsup = -jsup;
				int blaj = sprndr[jsup+1]-sprndr[jsup];
				if (jsup >= isupbeg) {
					if (icount==0) {
						temp.list[nlistloc] = isup;
						nlistloc++;
						icount = 1;
					};
					temp.bsa[nzjaloc] = nzaloc;
					temp.ja[nzjaloc] = jsupini;
					int ibs = _mtru2charr[ichild].bsa[j];
					int blaij = blai*blaj;
					for (int kii=0;kii<blaij;kii++) temp.a[nzaloc+kii] = _mtru2charr[ichild].a[ibs+kii];
					nzjaloc++;
					nzaloc += blaij;
				};
			};
			if (icount == 1) temp.ia[nlistloc] = nzjaloc;
		};
	};

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		_gmtru2.ReadBlock(iblk);

		for (int ilist=0;ilist<_gmtru2.mtrarr[iblk].nlist;ilist++) {
			int isup = _gmtru2.mtrarr[iblk].list[ilist];
			int blai = sprndr[isup+1]-sprndr[isup];
			int icount = 0;
			for (int j=_gmtru2.mtrarr[iblk].ia[ilist];j<_gmtru2.mtrarr[iblk].ia[ilist+1];j++) {
				int jsup = _gmtru2.mtrarr[iblk].ja[j];
				int jsupini=jsup;
				if (jsup < 0) jsup = -jsup;
				int blaj = sprndr[jsup+1]-sprndr[jsup];
				if (jsup >= isupbeg) {
					if (icount==0) {
						temp.list[nlistloc] = isup;
						nlistloc++;
						icount = 1;
					};
					temp.bsa[nzjaloc] = nzaloc;
					temp.ja[nzjaloc] = jsupini;
					int ibs = _gmtru2.mtrarr[iblk].bsa[j];
					int blaij = blai*blaj;
					for (int kii=0;kii<blaij;kii++) temp.a[nzaloc+kii] = _gmtru2.mtrarr[iblk].a[ibs+kii];
					nzjaloc++;
					nzaloc += blaij;
				};
			};
			if (icount == 1) temp.ia[nlistloc] = nzjaloc;
		};

		_gmtru2.FreeBlock(iblk);

	};

// Fill control data

	int i;

	for (i=0;i<=nlistloc;i++) temp.sprndr[i] = 0;
	for (i=0;i<=nlistloc;i++) temp.sprndc[i] = 0;

	temp.m = m;
	temp.n = n;
	temp.nsupr = nlistloc;
	temp.nsupc = nlistloc;
	temp.nlist = nlistloc;
	temp.nzja = nzjaloc;
	temp.nzatot = nzaloc;
	temp.nza = nzaloc;
	temp.blamx = _gmtru2.mtrarr[iblkbeg].blamx;

// Return the result

	_mtru2send = temp;

	CSMatrixCS mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Copy a part of data of the block into the new one
//========================================================================================
void CGSMatrixCS::CopyBlk2Blk (char _diatype, const CSMatrixCS &_a, CSMatrixCS &_anew) const { // Copy a part of data of the block into the new one

//	const char *funcname = "CopyBlk2Blk";

// Scan both block rows

	int ilist1, ilist2, isup1, isup2, blai, ibeg1, ibeg2, iend1, iend2;
	int jind1, jind2, jj1, jj2, blaj, ibs1, ibs2, kii;

	ilist1 = 0;
	ilist2 = 0;

	while (ilist1 < _a.nlist && ilist2 < _anew.nlist) {
		isup1 = _a.list[ilist1];
		isup2 = _anew.list[ilist2];
		if (isup1 == isup2) {
			blai = sprndr[isup1+1]-sprndr[isup1];
			ibeg1 = _a.ia[ilist1]; 
			iend1 = _a.ia[ilist1+1]-1; 
			ibeg2 = _anew.ia[ilist2]; 
			iend2 = _anew.ia[ilist2+1]-1; 
			jind1 = ibeg1;
			jind2 = ibeg2;
			while (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = _a.ja[jind1];
				jj2 = _anew.ja[jind2];
				if (jj1 == jj2) {
					blaj = sprndr[jj1+1]-sprndr[jj1];
					ibs1 = _a.bsa[jind1];
					ibs2 = _anew.bsa[jind2];
					if ((isup1 == jj1 && _diatype == 'D') || isup1 != jj1) {
						for (kii=0;kii<blai*blaj;kii++) _anew.a[ibs2+kii] = _a.a[ibs1+kii];
					};
				} else if (jj1 < jj2) {
					jind1++;
				} else if (jj1 > jj2) {
					jind2++;
				};
			};
		} else if (isup1 < isup2) {
			ilist1++;
		} else if (isup1 > isup2) {
			ilist2++;
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Copy a part of data of the block into the new one
//========================================================================================
void CGSMatrixCS::CopyBlkT2Blk (char _diatype, int &_icycle, int *_imask, 
								const CSMatrixCS &_a, CSMatrixCS &_anew) const { // Copy a part of data of the block into the new one

	const char *funcname = "CopyBlkT2Blk";

// Compute transposed structure of the given block

	int ilist1, ilist2, isup1, isup2, blai, ibeg1, ibeg2, iend1, iend2;
	int jind1, jind2, jj1, jj2, blaj, ibs1, ibs2, kii, kjj;
	int j, k, isup, jsup;

	int *listloc, *iptr, *iptr1, *ialoc, *jaloc, *bsaloc;

	listloc = new int [nsupr];
	if (!listloc) MemoryFail (funcname);
	iptr = new int [nsupr];
	if (!iptr) MemoryFail (funcname);
	iptr1 = new int [nsupr];
	if (!iptr1) MemoryFail (funcname);
	ialoc = new int [nsupr+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [_a.nzja];
	if (!jaloc) MemoryFail (funcname);
	bsaloc = new int [_a.nzja];
	if (!bsaloc) MemoryFail (funcname);

	_icycle++;

	int nlistloc = 0;

	for (ilist1=0;ilist1<_a.nlist;ilist1++) {
		for (j=_a.ia[ilist1];j<_a.ia[ilist1+1];j++) {
			jsup=_a.ja[j];
			if (_imask[jsup] != _icycle) {
				listloc[nlistloc] = jsup;
				nlistloc++;
				_imask[jsup] = _icycle;
				iptr[jsup] = 1;
			} else {
				iptr[jsup]++;
			};
		};
	};

	if (nlistloc != 0) qsort (listloc, nlistloc, sizeof(int), compint);

	ialoc[0] = 0;

	for (ilist1=0;ilist1<nlistloc;ilist1++) {
		jsup = listloc[ilist1];
		ialoc[ilist1+1] = ialoc[ilist1] + iptr[jsup];
		iptr1[jsup] = ilist1;
	};

	for (ilist1=0;ilist1<nlistloc;ilist1++) iptr[ilist1] = ialoc[ilist1];

	for (ilist1=0;ilist1<_a.nlist;ilist1++) {
		isup=_a.list[ilist1];
		for (j=_a.ia[ilist1];j<_a.ia[ilist1+1];j++) {
			jsup=_a.ja[j];
			ilist2 = iptr1[jsup];
			k = iptr[ilist2];
			jaloc[k] = isup;
			bsaloc[k] = _a.bsa[j];
			iptr[ilist2]++;
		};
	};

// Scan both block rows

	ilist1 = 0;
	ilist2 = 0;

	while (ilist1 < nlistloc && ilist2 < _anew.nlist) {
		isup1 = listloc[ilist1];
		isup2 = _anew.list[ilist2];
		if (isup1 == isup2) {
			blai = sprndr[isup1+1]-sprndr[isup1];
			ibeg1 = ialoc[ilist1]; 
			iend1 = ialoc[ilist1+1]-1; 
			ibeg2 = _anew.ia[ilist2]; 
			iend2 = _anew.ia[ilist2+1]-1; 
			jind1 = ibeg1;
			jind2 = ibeg2;
			while (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = jaloc[jind1];
				jj2 = _anew.ja[jind2];
				if (jj1 == jj2) {
					blaj = sprndr[jj1+1]-sprndr[jj1];
					ibs1 = bsaloc[jind1];
					ibs2 = _anew.bsa[jind2];
					if ((isup1 == jj1 && _diatype == 'D') || isup1 != jj1) {
						for (kii=0;kii<blai;kii++) {
							for (kjj=0;kjj<blaj;kjj++) {
								_anew.a[ibs2+kjj*blai+kii] = _a.a[ibs1+kii*blaj+kjj];
							};
						};
					};
				} else if (jj1 < jj2) {
					jind1++;
				} else if (jj1 > jj2) {
					jind2++;
				};
			};
		} else if (isup1 < isup2) {
			ilist1++;
		} else if (isup1 > isup2) {
			ilist2++;
		};
	};

// Free work arrays

	delete [] listloc;
	delete [] iptr;
	delete [] iptr1;
	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Write current block into full format
//========================================================================================
void CGSMatrixCS::FullFormat (int _iblk, int _ntotal, dcmplx *_a) { // Write current block into full format

//	const char *funcname = "FullFormat";

// Check sizes

	if (m != _ntotal || n != _ntotal) {
		throw " FullFormat: incorrect matrix sizes ";
	};

	if (_iblk < 0 || _iblk >= nblksr) {
		throw " FullFormat: wrong block number ";
	};

// Read block data if necessary

	ReadBlock (_iblk);

// Cycle over the supernodes

	int ilist, isup, blai, ibegi, j, jsup, ibegj, blaj, ibs, kii, kjj;

	for (ilist=0;ilist<mtrarr[_iblk].nlist;ilist++) {
		isup = mtrarr[_iblk].list[ilist];
		blai = sprndr[isup+1]-sprndr[isup];
		ibegi = sprndr[isup];
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jsup = mtrarr[_iblk].ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			ibegj = sprndr[jsup];
			ibs = mtrarr[_iblk].bsa[j];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					_a[(ibegj+kjj)*_ntotal+ibegi+kii] = mtrarr[_iblk].a[ibs+kjj*blai+kii];
				};
			};
		};
	};

// Free block data if necessary

	FreeBlock (_iblk);

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Compute the sparsity structure of all local block rows
//========================================================================================
CSMatrix CGSMatrixCS::GlobalSparsity () const { // Compute the sparsity structure of all local block rows

//	const char *funcname = "GlobalSparsity";

// Count the number of elements in the matrix

	int ilist, iblk, i, j;
	int nlistloc=0, nzloc=0;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		nlistloc += mtrarr[iblk].nlist;
		nzloc += mtrarr[iblk].nzja;
	};

// Allocate and init the matrix

	CSMatrix aloc(nlistloc,nzloc);

	aloc.m = nsupr;
	aloc.n = nsupr;
	aloc.nsupr = nsupr;
	aloc.nsupc = nsupr;
	aloc.nlist = nlistloc;
	aloc.nzja = nzloc;

	nlistloc=0;
	nzloc=0;

	aloc.ia[0] = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		for (i=0;i<mtrarr[iblk].nlist;i++) {
			aloc.list[nlistloc] = mtrarr[iblk].list[i];
			for (j=mtrarr[iblk].ia[i];j<mtrarr[iblk].ia[i+1];j++) {
				aloc.ja[nzloc] = mtrarr[iblk].ja[j];
				nzloc++;
			};
			aloc.ia[nlistloc+1] = nzloc;
			nlistloc++;
		};
	};

// Return the result

	return aloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Compute the point sparse matrix via very dense rectangular matrix
//========================================================================================
CSMatrixCS CGSMatrixCS::Sparsify (int _nitermax, int _nza) { // Compute the point sparse matrix via very dense rectangular matrix

	const char *funcname = "Sparsify";

// Initial scan for scaling and maximal values

	int nloc = n;
	int mloc = m;

	int *jamax;
	double *valmax, *valmin;
	dcmplx *amtrmax;
	double *scl;

	jamax = new int [nloc];
	if (!jamax) MemoryFail (funcname);
	valmax = new double [nloc];
	if (!valmax) MemoryFail (funcname);
	valmin = new double [nloc];
	if (!valmin) MemoryFail (funcname);
	amtrmax = new dcmplx [nloc];
	if (!amtrmax) MemoryFail (funcname);
	scl = new double [nloc];
	if (!scl) MemoryFail (funcname);

	int iblk, i, j, nz, isup, blai, ibeg, jsup, blaj, jbeg;
	int kii, kjj, ind, nztot;
	double aux;

	int nlistloc;
	int *listloc, *ialoc, *jaloc, *bsaloc;
	dcmplx *aaloc;

	for (i=0;i<nloc;i++) jamax[i] = -1;
	for (i=0;i<nloc;i++) valmax[i] = -1.0e0;
	for (i=0;i<nloc;i++) valmin[i] = 1.0e10;
	for (i=0;i<nloc;i++) scl[i] = 0.0e0;

	for (iblk=0;iblk<nblksc;iblk++) {

		ReadBlock (iblk);

		nlistloc = mtrarr[iblk].GetNlist ();

		mtrarr[iblk].GetSparsity (listloc, ialoc, jaloc);
		bsaloc = mtrarr[iblk].GetBsa ();
		aaloc = mtrarr[iblk].GetA ();

		nz = 0;
		for (i=0;i<nlistloc;i++) {
			isup = listloc[i];
			blai = sprndc[isup+1]-sprndc[isup];
			ibeg = sprndc[isup];
			for (j=ialoc[i];j<ialoc[i+1];j++) {
				jsup = jaloc[j];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				jbeg = sprndr[jsup];
				for (kjj=0;kjj<blaj;kjj++) {
					for (kii=0;kii<blai;kii++) {
						ind = kii+ibeg;
						aux = aaloc[nz].x * aaloc[nz].x + aaloc[nz].y * aaloc[nz].y;
						scl[ind] += aux;
						if (aux < valmin[ind]) {
							valmin[ind] = aux;
						};
						if (aux > valmax[ind]) {
							jamax[ind] = kjj+jbeg;
							valmax[ind] = aux;
							amtrmax[ind] = aaloc[nz];
						};
						nz++;
					};
				};
			};
		};

		FreeBlock (iblk);

	};

	double threshmax = 0.0e0;
	double threshmin = 10.0e0;

	for (i=0;i<nloc;i++) {
		aux = scl[i];
		aux = 1.0e0 / sqrt(aux);
		scl[i] = aux;
		valmin[i] *= scl[i]*scl[i];
		valmax[i] *= scl[i]*scl[i];
		if (valmin[i] < threshmin) threshmin = valmin[i];
		if (valmax[i] > threshmax) threshmax = valmax[i];
	};

	threshmin = sqrt(threshmin) * 1.0e-1;
	threshmax = sqrt(threshmax);

// Iterative cycle on threshold value

	double thresh = sqrt(threshmax*threshmin);
	double thresh_2;

	int iter;

	for (iter=0;iter<_nitermax;iter++) {

		thresh_2 = thresh * thresh;

		nztot = 0;

		for (iblk=0;iblk<nblksc;iblk++) {

			ReadBlock (iblk);

			nlistloc = mtrarr[iblk].GetNlist ();

			mtrarr[iblk].GetSparsity (listloc, ialoc, jaloc);
			bsaloc = mtrarr[iblk].GetBsa ();
			aaloc = mtrarr[iblk].GetA ();

			nz = 0;
			for (i=0;i<nlistloc;i++) {
				isup = listloc[i];
				blai = sprndc[isup+1]-sprndc[isup];
				ibeg = sprndc[isup];
				for (j=ialoc[i];j<ialoc[i+1];j++) {
					jsup = jaloc[j];
					blaj = sprndr[jsup+1]-sprndr[jsup];
					for (kjj=0;kjj<blaj;kjj++) {
						for (kii=0;kii<blai;kii++) {
							ind = kii+ibeg;
							aux = aaloc[nz].x * aaloc[nz].x + aaloc[nz].y * aaloc[nz].y;
							aux *= scl[ind] * scl[ind];
							if (aux > thresh_2) nztot++;
							nz++;
						};
					};
				};
			};

			FreeBlock (iblk);

		};

		cout << " Iter = " << iter << " ThrMin = " << threshmin << " ThrMax = " << threshmax << " Thr = " << thresh << endl;
		cout << "    NzaOpt = " << _nza << " Nza = " << nztot << endl;

		if (iter < _nitermax-1) {
			if (nztot > _nza) {
				threshmin = thresh;
			} else {
				threshmax = thresh;
			};
			thresh = sqrt(threshmax*threshmin);
		};

	};

// Compute sparse matrix according to the threshold

	int *iasave;

	iasave = new int [nloc+1];
	if (!iasave) MemoryFail (funcname);

	for (i=0;i<=nloc;i++) iasave[i] = 0;

	thresh_2 = thresh * thresh;

	for (iblk=0;iblk<nblksc;iblk++) {

		ReadBlock (iblk);

		nlistloc = mtrarr[iblk].GetNlist ();

		mtrarr[iblk].GetSparsity (listloc, ialoc, jaloc);
		bsaloc = mtrarr[iblk].GetBsa ();
		aaloc = mtrarr[iblk].GetA ();

		nz = 0;
		for (i=0;i<nlistloc;i++) {
			isup = listloc[i];
			blai = sprndc[isup+1]-sprndc[isup];
			ibeg = sprndc[isup];
			for (j=ialoc[i];j<ialoc[i+1];j++) {
				jsup = jaloc[j];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				for (kjj=0;kjj<blaj;kjj++) {
					for (kii=0;kii<blai;kii++) {
						ind = kii+ibeg;
						aux = aaloc[nz].x * aaloc[nz].x + aaloc[nz].y * aaloc[nz].y;
						aux *= scl[ind] * scl[ind];
						if (aux > thresh_2) iasave[ind+1]++;
						nz++;
					};
				};
			};
		};

		FreeBlock (iblk);

	};

	for (i=0;i<nloc;i++) {
		if (iasave[i+1] == 0) iasave[i+1] = 1;
	};

	for (i=0;i<nloc;i++) iasave[i+1] = iasave[i] + iasave[i+1];

	nztot = iasave[nloc];

	int *jasave;
	int *iptr;
	dcmplx *asave;

	jasave = new int [nztot];
	if (!jasave) MemoryFail (funcname);
	asave = new dcmplx [nztot];
	if (!asave) MemoryFail (funcname);
	iptr = new int [nloc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nloc;i++) iptr[i] = iasave[i];

	int k;

	for (iblk=0;iblk<nblksc;iblk++) {

		ReadBlock (iblk);

		nlistloc = mtrarr[iblk].GetNlist ();

		mtrarr[iblk].GetSparsity (listloc, ialoc, jaloc);
		bsaloc = mtrarr[iblk].GetBsa ();
		aaloc = mtrarr[iblk].GetA ();

		nz = 0;
		for (i=0;i<nlistloc;i++) {
			isup = listloc[i];
			blai = sprndc[isup+1]-sprndc[isup];
			ibeg = sprndc[isup];
			for (j=ialoc[i];j<ialoc[i+1];j++) {
				jsup = jaloc[j];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				jbeg = sprndr[jsup];
				for (kjj=0;kjj<blaj;kjj++) {
					for (kii=0;kii<blai;kii++) {
						ind = kii+ibeg;
						aux = aaloc[nz].x * aaloc[nz].x + aaloc[nz].y * aaloc[nz].y;
						aux *= scl[ind] * scl[ind];
						if (aux > thresh_2 && iptr[ind] < iasave[ind+1]) {
							k = iptr[ind];
							jasave[k] = jbeg+kjj;
							asave[k] = aaloc[nz];
							iptr[ind]++;
						};
						nz++;
					};
				};
			};
		};

		FreeBlock (iblk);

	};

	for (i=0;i<nloc;i++) {
		if (iptr[i] < iasave[i+1]-1) {
			throw " CGSMatrixCS::Sparsify: error when filtering matrix entries";
		} else if (iptr[i] == iasave[i+1]-1) {
			k = iptr[i];
			jasave[k] = jamax[i];
			asave[k] = amtrmax[i];
		};
	};

// Allocate and init the matrix

	CSMatrixCS aloc(mloc,nztot,nztot);

	aloc.m = mloc;
	aloc.n = nloc;
	aloc.nsupr = mloc;
	aloc.nsupc = nloc;
	aloc.nlist = nloc;
	aloc.nzja = nztot;
	aloc.nza = nztot;
	aloc.nzatot = nztot;

	for (i=0;i<nloc;i++) {
		aloc.list[i] = i;
	};

	for (i=0;i<=nloc;i++) aloc.ia[i] = iasave[i];
	for (i=0;i<nztot;i++) aloc.ja[i] = jasave[i];
	for (i=0;i<nztot;i++) aloc.a[i] = asave[i];

	aloc.blamx = 1;

	for (i=0;i<=nloc;i++) aloc.sprndc[i] = i;
	for (i=0;i<=nloc;i++) aloc.sprndr[i] = i;
	for (i=0;i<nztot;i++) aloc.bsa[i] = i;

// Free work arrays

	delete [] jamax;
	delete [] valmax;
	delete [] valmin;
	delete [] amtrmax;
	delete [] scl;
	delete [] iasave;
	delete [] jasave;
	delete [] asave;
	delete [] iptr;

// Return the result

	return aloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CGSMatrixCS &_a) { // Output matrix

	_stream << " CGSMatrixCS:" << endl;

	_stream << (const CGSMatrix &) _a << endl;

	for (int iblk=0;iblk<_a.nblksr;iblk++) {
		_stream << " IblkRow = " << iblk << endl;
		_stream << " BlockRow = " << _a.mtrarr[iblk] << endl;
	};

	_stream << " Ivect = " << _a.ivect << endl;
	_stream << " Vector1 = " << _a.vector1 << endl;
	_stream << " Vector2 = " << _a.vector2 << endl;

	return _stream;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
