//------------------------------------------------------------------------------------------------
// File: gsmatrix.cpp
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
// Description: Memory allocation zero data constructor
// CGSMatrix::CGSMatrix()
//========================================================================================
CGSMatrix::CGSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrix_00";

	m = 0;
	n = 0;
	nblksr = 0;
	nblksc = 0;
	nsupr = 0;
	nsupc = 0;
	nlist = 0;
	nlistcnd = 0;
	nfiles = 0;
	ifile = 0;

	int i;

	blksr = new int [nblksr+1];
	if (!blksr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blksr[i] = 0;
	};
	blksc = new int [nblksc+1];
	if (!blksc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		blksc[i] = 0;
	};
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	for (i=0;i<nsupr+1;i++) {
		sprndr[i] = 0;
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<nsupc+1;i++) {
		sprndc[i] = 0;
	};
	bl2ndr = new int [nblksr+1];
	if (!bl2ndr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		bl2ndr[i] = 0;
	};
	bl2ndc = new int [nblksc+1];
	if (!bl2ndc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		bl2ndc[i] = 0;
	};
	listb = new int [nlist];
	if (!listb) MemoryFail (funcname);
	for (i=0;i<nlist;i++) {
		listb[i] = i;
	};
	blkscnd = new int [nblksr+1];
	if (!blkscnd) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blkscnd[i] = 0;
	};
	listcnd = new int [nlistcnd];
	if (!listcnd) MemoryFail (funcname);
	for (i=0;i<nlistcnd;i++) {
		listcnd[i] = 0;
	};
	iblkcopymask = new int [nblksr+1];
	if (!iblkcopymask) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		iblkcopymask[i] = 0;
	};
	ptrsparr = new int * [6*nblksr];
	if (!ptrsparr) MemoryFail (funcname);
	bl2file = new int [nblksc];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [nblksc];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};
};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CGSMatrix::CGSMatrix()
//========================================================================================
CGSMatrix::CGSMatrix (int _nblks, int _nsups, int _nlist) { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrix_01";

	m = 0;
	n = 0;
	nblksr = _nblks;
	nblksc = _nblks;
	nsupr = _nsups;
	nsupc = _nsups;
	nlist = _nlist;
	nlistcnd = 0;
	nfiles = 0;
	ifile = 0;

	blksr = new int [nblksr+1];
	if (!blksr) MemoryFail (funcname);

	int i;

	for (i=0;i<nblksr+1;i++) {
		blksr[i] = 0;
	};
	blksc = new int [nblksc+1];
	if (!blksc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		blksc[i] = 0;
	};
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	for (i=0;i<nsupr+1;i++) {
		sprndr[i] = 0;
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<nsupc+1;i++) {
		sprndc[i] = 0;
	};
	bl2ndr = new int [nblksr+1];
	if (!bl2ndr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		bl2ndr[i] = 0;
	};
	bl2ndc = new int [nblksc+1];
	if (!bl2ndc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		bl2ndc[i] = 0;
	};
	listb = new int [nlist];
	if (!listb) MemoryFail (funcname);
	for (i=0;i<nlist;i++) {
		listb[i] = i;
	};
	blkscnd = new int [nblksr+1];
	if (!blkscnd) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blkscnd[i] = 0;
	};
	listcnd = new int [nlistcnd];
	if (!listcnd) MemoryFail (funcname);
	for (i=0;i<nlistcnd;i++) {
		listcnd[i] = 0;
	};
	iblkcopymask = new int [nblksr+1];
	if (!iblkcopymask) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		iblkcopymask[i] = 0;
	};
	ptrsparr = new int * [6*nblksr];
	if (!ptrsparr) MemoryFail (funcname);
	bl2file = new int [nblksc];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [nblksc];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

};

// Author: Kharchenko S.A.
// Description: Destructor
// CGSMatrix::~CGSMatrix()
//========================================================================================
CGSMatrix::~CGSMatrix () { // Destructor
//	std::cout << " On entry to CGSMatrix destructor " << std::endl;
	m      = 0;
	n      = 0;
	nblksr = 0;
	nblksc = 0;
	nsupr  = 0;
	nsupc  = 0;
	nlist  = 0;
	nlistcnd = 0;
	delete [] blksr;
	delete [] blksc;
	delete [] sprndr;
	delete [] sprndc;
	delete [] bl2ndr;
	delete [] bl2ndc;
	delete [] listb;
	delete [] blkscnd;
	delete [] listcnd;
	delete [] iblkcopymask;
	delete [] ptrsparr;
	int i;
	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	nfiles = 0;
	ifile  = 0;
	delete [] bl2file;
	delete [] bl2ibsfile;
	delete [] ibsfile;
	delete [] files;
	delete [] fnames;
//	std::cout << " On return from CGSMatrix destructor " << std::endl;
};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CGSMatrix::CGSMatrix()
//========================================================================================
CGSMatrix::CGSMatrix (const CGSMatrix &_aa) { // Copy constructor

	const char *funcname = "CGSMatrix_copy";

	m = _aa.m;
	n = _aa.n;
	nblksr = _aa.nblksr;
	nblksc = _aa.nblksc;
	nsupr = _aa.nsupr;
	nsupc = _aa.nsupc;
	nlist = _aa.nlist;
	nlistcnd = _aa.nlistcnd;
	nfiles = 0;
	ifile = 0;

	int i;

	blksr = new int [nblksr+1];
	if (!blksr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blksr[i] = _aa.blksr[i];
	};
	blksc = new int [nblksc+1];
	if (!blksc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		blksc[i] = _aa.blksc[i];
	};
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	for (i=0;i<nsupr+1;i++) {
		sprndr[i] = _aa.sprndr[i];
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<nsupc+1;i++) {
		sprndc[i] = _aa.sprndc[i];
	};
	bl2ndr = new int [nblksr+1];
	if (!bl2ndr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		bl2ndr[i] = _aa.bl2ndr[i];
	};
	bl2ndc = new int [nblksc+1];
	if (!bl2ndc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		bl2ndc[i] = _aa.bl2ndc[i];
	};
	listb = new int [nlist];
	if (!listb) MemoryFail (funcname);
	for (i=0;i<nlist;i++) {
		listb[i] = _aa.listb[i];
	};
	blkscnd = new int [nblksr+1];
	if (!blkscnd) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blkscnd[i] = _aa.blkscnd[i];
	};
	listcnd = new int [nlistcnd];
	if (!listcnd) MemoryFail (funcname);
	for (i=0;i<nlistcnd;i++) {
		listcnd[i] = _aa.listcnd[i];
	};
	iblkcopymask = new int [nblksr+1];
	if (!iblkcopymask) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		iblkcopymask[i] = 0;
	};
	ptrsparr = new int * [6*nblksr];
	if (!ptrsparr) MemoryFail (funcname);
	bl2file = new int [nblksc];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [nblksc];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CGSMatrix::operator=()
//========================================================================================
CGSMatrix &CGSMatrix::operator= (const CGSMatrix &_aa) { // Equality operator

	const char *funcname = "CGSMatrix_=";

	delete [] blksr;
	delete [] blksc;
	delete [] sprndr;
	delete [] sprndc;
	delete [] bl2ndr;
	delete [] bl2ndc;
	delete [] listb;
	delete [] blkscnd;
	delete [] listcnd;
	delete [] iblkcopymask;
	delete [] ptrsparr;
	int i;
	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	nfiles = 0;
	ifile  = 0;
	delete [] bl2file;
	delete [] bl2ibsfile;
	delete [] ibsfile;
	delete [] files;
	delete [] fnames;

	m = _aa.m;
	n = _aa.n;
	nblksr = _aa.nblksr;
	nblksc = _aa.nblksc;
	nsupr = _aa.nsupr;
	nsupc = _aa.nsupc;
	nlist = _aa.nlist;
	nfiles = 0;
	ifile = 0;

	blksr = new int [nblksr+1];
	if (!blksr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blksr[i] = _aa.blksr[i];
	};
	blksc = new int [nblksc+1];
	if (!blksc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		blksc[i] = _aa.blksc[i];
	};
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	for (i=0;i<nsupr+1;i++) {
		sprndr[i] = _aa.sprndr[i];
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<nsupc+1;i++) {
		sprndc[i] = _aa.sprndc[i];
	};
	bl2ndr = new int [nblksr+1];
	if (!bl2ndr) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		bl2ndr[i] = _aa.bl2ndr[i];
	};
	bl2ndc = new int [nblksc+1];
	if (!bl2ndc) MemoryFail (funcname);
	for (i=0;i<nblksc+1;i++) {
		bl2ndc[i] = _aa.bl2ndc[i];
	};
	listb = new int [nlist];
	if (!listb) MemoryFail (funcname);
	for (i=0;i<nlist;i++) {
		listb[i] = _aa.listb[i];
	};
	blkscnd = new int [nblksr+1];
	if (!blkscnd) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		blkscnd[i] = _aa.blkscnd[i];
	};
	listcnd = new int [nlistcnd];
	if (!listcnd) MemoryFail (funcname);
	for (i=0;i<nlistcnd;i++) {
		listcnd[i] = _aa.listcnd[i];
	};
	iblkcopymask = new int [nblksr+1];
	if (!iblkcopymask) MemoryFail (funcname);
	for (i=0;i<nblksr+1;i++) {
		iblkcopymask[i] = 0;
	};
	ptrsparr = new int * [6*nblksr];
	if (!ptrsparr) MemoryFail (funcname);
	bl2file = new int [nblksc];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [nblksc];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

	return *this;
};

// Author: Kharchenko S.A.
// Description: Setup files
// CGSMatrix::SetupFiles()
//========================================================================================
void CGSMatrix::SetupFiles (int _nfiles, char *_name) { // Setup files

	const char *funcname = "SetupFiles";

	int i;

	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	delete [] ibsfile;
	delete [] files;
	delete [] fnames;

	nfiles = _nfiles;
	ifile = 0;

	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

	OpenDIOFiles (_nfiles, files, fnames, "", _name, ".bin");

};

// Author: Kharchenko S.A.
// Description: Close files
// CGSMatrix::CloseFiles()
//========================================================================================
void CGSMatrix::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

//	CloseDIOFiles (nfiles, files);
	CloseDIOFiles (nfiles, files, fnames);

};

// Author: Kharchenko S.A.
// Description: Free files
// CGSMatrix::FreeFiles()
//========================================================================================
void CGSMatrix::FreeFiles () { // Free files

//	const char *funcname = "FreeFiles";

	ifile = 0;

	int i;

	for (i=0;i<nblksc;i++) {
		bl2file[i] = -1;
	};
	for (i=0;i<nblksc;i++) {
		bl2ibsfile[i] = -1;
	};
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the block row
// CGSMatrix::BlockSparsity2Index()
//========================================================================================
void CGSMatrix::BlockSparsity2Index (const CSMatrix &_mtra, // Compute block sparsity of the block row
								int &_nlstblk, int *_lstblk,
								int &_icycle, int *_imask) const {

//	const char *funcname = "BlockSparsity";

	_icycle++;

	_nlstblk = 0;

	for (int ilist=0;ilist<_mtra.nlist;ilist++) {
		for (int j=_mtra.ia[ilist];j<_mtra.ia[ilist+1];j++) {
			int jblk = _mtra.ja2[j];
			if (_imask[jblk] != _icycle) {
				_lstblk[_nlstblk] = jblk;
				_imask[jblk] = _icycle;
				_nlstblk++;
			};
		};
	};

	if (_nlstblk != 0) qsort (_lstblk, _nlstblk, sizeof(int), compint);

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the block row
// CGSMatrix::BlockSparsity()
//========================================================================================
void CGSMatrix::BlockSparsity (const CSMatrix &_mtra, // Compute block sparsity of the block row
								int &_nlstblk, int *_lstblk,
								int &_icycle, int *_imask, 
								const int *_sp2blk) const {

//	const char *funcname = "BlockSparsity";

	_icycle++;

	_nlstblk = 0;

	for (int ilist=0;ilist<_mtra.nlist;ilist++) {
		for (int j=_mtra.ia[ilist];j<_mtra.ia[ilist+1];j++) {
			int jsup = _mtra.ja[j];
			if (jsup < 0) jsup = -jsup;
			int jblk = _sp2blk[jsup];
			if (_imask[jblk] != _icycle) {
				_lstblk[_nlstblk] = jblk;
				_imask[jblk] = _icycle;
				_nlstblk++;
			};
		};
	};

	if (_nlstblk != 0) qsort (_lstblk, _nlstblk, sizeof(int), compint);

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the block row
// CGSMatrix::MtrBlockSparsity2Index()
//========================================================================================
void CGSMatrix::MtrBlockSparsity2Index (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
										int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity2Index";

// Allocate and init working arrays

	int *lstblk, *imask;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		BlockSparsity2Index (*_mtrarr[iblk], 
						nlstblk, lstblk,
						icycle, imask);

		nzjab += nlstblk;

	};

// Allocate and compute arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjab];
	if (!_jab) MemoryFail (funcname);

	_iab[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		BlockSparsity2Index (*_mtrarr[iblk], 
							nlstblk, _jab+nzjab,
							icycle, imask);

		nzjab += nlstblk;

		_iab[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the matrix
// CGSMatrix::MtrBlockSparsity()
//========================================================================================
void CGSMatrix::MtrBlockSparsity (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
									int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity";

// Allocate and init working arrays

	int *lstblk, *imask, *sp2blk;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);
	sp2blk = new int [nsupr];
	if (!sp2blk) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};
	for (iblk=0;iblk<nblksr;iblk++) {
		for (int isup=blksr[iblk];isup<blksr[iblk+1];isup++) {
			sp2blk[isup] = iblk;
		};
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		BlockSparsity (*_mtrarr[iblk], 
						nlstblk, lstblk,
						icycle, imask, sp2blk);

		nzjab += nlstblk;

	};

// Allocate and compute arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjab];
	if (!_jab) MemoryFail (funcname);

	_iab[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		BlockSparsity (*_mtrarr[iblk], 
						nlstblk, _jab+nzjab,
						icycle, imask, sp2blk);

		nzjab += nlstblk;

		_iab[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;
	delete [] sp2blk;

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the matrix
// CGSMatrix::MtrBlockSparsity2Index()
//========================================================================================
void CGSMatrix::MtrBlockSparsity2Index (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
										int _iblkbeg, int _iblkend, 
										int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity2Index";

// Allocate and init working arrays

	int *lstblk, *imask;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {

		BlockSparsity2Index (*_mtrarr[iblk], 
						nlstblk, lstblk,
						icycle, imask);

		nzjab += nlstblk;

	};

// Allocate and compute arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjab];
	if (!_jab) MemoryFail (funcname);

	_iab[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (iblk >= _iblkbeg && iblk <= _iblkend) {

			BlockSparsity2Index (*_mtrarr[iblk], 
								nlstblk, _jab+nzjab,
								icycle, imask);

			nzjab += nlstblk;

		};

		_iab[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the matrix
// CGSMatrix::MtrBlockSparsity()
//========================================================================================
void CGSMatrix::MtrBlockSparsity (CSMatrix **_mtrarr, int _iblkbeg, int _iblkend, // Compute block sparsity of the matrix
									int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity";

// Allocate and init working arrays

	int *lstblk, *imask, *sp2blk;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);
	sp2blk = new int [nsupr];
	if (!sp2blk) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};
	for (iblk=0;iblk<nblksr;iblk++) {
		for (int isup=blksr[iblk];isup<blksr[iblk+1];isup++) {
			sp2blk[isup] = iblk;
		};
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {

		BlockSparsity (*_mtrarr[iblk], 
						nlstblk, lstblk,
						icycle, imask, sp2blk);

		nzjab += nlstblk;

	};

// Allocate and compute arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjab];
	if (!_jab) MemoryFail (funcname);

	_iab[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (iblk >= _iblkbeg && iblk <= _iblkend) {

			BlockSparsity (*_mtrarr[iblk], 
							nlstblk, _jab+nzjab,
							icycle, imask, sp2blk);

			nzjab += nlstblk;

		};

		_iab[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;
	delete [] sp2blk;

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the matrix in the parallel mode
// CGSMatrix::MtrBlockSparsity2Index()
//========================================================================================
void CGSMatrix::MtrBlockSparsity2Index (bool _check, const CTree &_tree, int *_bl2cpu, // Compute block sparsity of the matrix in the parallel mode
									CSMatrix **_mtrarr, int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity2Index";

// Allocate and init working arrays

	int *lstblk, *imask;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (_bl2cpu[iblk] == _tree.myid) {

			BlockSparsity2Index (*_mtrarr[iblk], 
							nlstblk, lstblk,
							icycle, imask);

			if (nlstblk == 0) {
				lstblk[0] = iblk;
				nlstblk = 1;
			};

			nzjab += nlstblk;

		};

	};

// Allocate and compute arrays that describe local block sparsity

	int *iabloc, *jabloc;

	iabloc = new int [nblksr+1];
	if (!iabloc) MemoryFail (funcname);
	jabloc = new int [nzjab];
	if (!jabloc) MemoryFail (funcname);

	iabloc[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (_bl2cpu[iblk] == _tree.myid) {

			BlockSparsity2Index (*_mtrarr[iblk], 
							nlstblk, jabloc+nzjab,
							icycle, imask);

			if (nlstblk == 0) {
				jabloc[nzjab] = iblk;
				nlstblk = 1;
			};

			nzjab += nlstblk;

		};

		iabloc[iblk+1] = nzjab;

	};

// Compute the total number of nonzero blocks in the matrix

	int nzjabtot = 0;

	int *nzarr, *nzsum;

	nzarr = new int [_tree.nproc];
	if (!nzarr) MemoryFail (funcname);
	nzsum = new int [_tree.nproc+1];
	if (!nzsum) MemoryFail (funcname);

	int i;

	for (i=0;i<_tree.nproc;i++) nzarr[i] = 0;
	nzarr[_tree.myid] = nzjab;
	nzsum[0] = 0;

	if (_tree.nproc == 1) {
		nzsum[1] = nzjab;
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													_tree.nproc, nzarr, nzsum+1);
		for (i=0;i<_tree.nproc;i++) nzsum[i+1] += nzsum[i];
	};

	nzjabtot = nzsum[_tree.nproc];

// Allocate arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjabtot];
	if (!_jab) MemoryFail (funcname);

	for (i=0;i<=nblksr;i++) _iab[i] = 0;

// Exchange IA arrays

	int *iasnd;

	iasnd = new int [nblksr];
	if (!iasnd) MemoryFail (funcname);

	for (iblk=0;iblk<nblksr;iblk++) {
		if (_bl2cpu[iblk] == _tree.myid) {
			iasnd[iblk] = iabloc[iblk+1]-iabloc[iblk];
		} else {
			iasnd[iblk] = 0;
		};
	};

	if (_tree.nproc == 1) {
		for (iblk=0;iblk<=nblksr;iblk++) _iab[iblk] = iabloc[iblk];
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nblksr, iasnd, _iab+1);
		for (i=0;i<nblksr;i++) _iab[i+1] += _iab[i];
	};

// Exchange JA arrays

	int *jasnd;

	jasnd = new int [nzjabtot];
	if (!jasnd) MemoryFail (funcname);

	for (i=0;i<nzjabtot;i++) jasnd[i] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {
		if (_bl2cpu[iblk] == _tree.myid) {
			for (int j=_iab[iblk];j<_iab[iblk+1];j++) {
				jasnd[j] = jabloc[nzjab];
				nzjab++;
			};
		};
	};

	if (_tree.nproc == 1) {
		for (i=0;i<nzjabtot;i++) _jab[i] = jabloc[i];
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nzjabtot, jasnd, _jab);
	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;
	delete [] nzarr;
	delete [] nzsum;
	delete [] iabloc;
	delete [] jabloc;
	delete [] iasnd;
	delete [] jasnd;

};

// Author: Kharchenko S.A.
// Description: Compute block sparsity of the matrix in the parallel mode
// CGSMatrix::MtrBlockSparsity()
//========================================================================================
void CGSMatrix::MtrBlockSparsity (const CTree &_tree, int *_bl2cpu, // Compute block sparsity of the matrix in the parallel mode
									CSMatrix **_mtrarr, int *&_iab, int *&_jab) const {

	const char *funcname = "MtrBlockSparsity";

// Allocate and init working arrays

	int *lstblk, *imask, *sp2blk;

	lstblk = new int [nblksr];
	if (!lstblk) MemoryFail (funcname);
	imask = new int [nblksr];
	if (!imask) MemoryFail (funcname);
	sp2blk = new int [nsupr];
	if (!sp2blk) MemoryFail (funcname);

	int icycle = -1;

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {
		imask[iblk] = icycle;
	};
	for (iblk=0;iblk<nblksr;iblk++) {
		for (int isup=blksr[iblk];isup<blksr[iblk+1];isup++) {
			sp2blk[isup] = iblk;
		};
	};

// Count the number of elements

	int nzjab = 0;
	int nlstblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (_bl2cpu[iblk] == _tree.myid) {

			BlockSparsity (*_mtrarr[iblk], 
							nlstblk, lstblk,
							icycle, imask, sp2blk);

			if (nlstblk == 0) {
				lstblk[0] = iblk;
				nlstblk = 1;
			};

			nzjab += nlstblk;

		};

	};

// Allocate and compute arrays that describe local block sparsity

	int *iabloc, *jabloc;

	iabloc = new int [nblksr+1];
	if (!iabloc) MemoryFail (funcname);
	jabloc = new int [nzjab];
	if (!jabloc) MemoryFail (funcname);

	iabloc[0] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {

		if (_bl2cpu[iblk] == _tree.myid) {

			BlockSparsity (*_mtrarr[iblk], 
							nlstblk, jabloc+nzjab,
							icycle, imask, sp2blk);

			if (nlstblk == 0) {
				jabloc[nzjab] = iblk;
				nlstblk = 1;
			};

			nzjab += nlstblk;

		};

		iabloc[iblk+1] = nzjab;

	};

// Compute the total number of nonzero blocks in the matrix

	int nzjabtot = 0;

	int *nzarr, *nzsum;

	nzarr = new int [_tree.nproc];
	if (!nzarr) MemoryFail (funcname);
	nzsum = new int [_tree.nproc+1];
	if (!nzsum) MemoryFail (funcname);

	int i;

	for (i=0;i<_tree.nproc;i++) nzarr[i] = 0;
	nzarr[_tree.myid] = nzjab;
	nzsum[0] = 0;

	if (_tree.nproc == 1) {
		nzsum[1] = nzjab;
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													_tree.nproc, nzarr, nzsum+1);
		for (i=0;i<_tree.nproc;i++) nzsum[i+1] += nzsum[i];
	};

	nzjabtot = nzsum[_tree.nproc];

// Allocate arrays to be returned

	_iab = new int [nblksr+1];
	if (!_iab) MemoryFail (funcname);
	_jab = new int [nzjabtot];
	if (!_jab) MemoryFail (funcname);

	for (i=0;i<=nblksr;i++) _iab[i] = 0;

// Exchange IA arrays

	int *iasnd;

	iasnd = new int [nblksr];
	if (!iasnd) MemoryFail (funcname);

	for (iblk=0;iblk<nblksr;iblk++) {
		if (_bl2cpu[iblk] == _tree.myid) {
			iasnd[iblk] = iabloc[iblk+1]-iabloc[iblk];
		} else {
			iasnd[iblk] = 0;
		};
	};

	if (_tree.nproc == 1) {
		for (iblk=0;iblk<=nblksr;iblk++) _iab[iblk] = iabloc[iblk];
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nblksr, iasnd, _iab+1);
		for (i=0;i<nblksr;i++) _iab[i+1] += _iab[i];
	};

// Exchange JA arrays

	int *jasnd;

	jasnd = new int [nzjabtot];
	if (!jasnd) MemoryFail (funcname);

	for (i=0;i<nzjabtot;i++) jasnd[i] = 0;

	nzjab = 0;

	for (iblk=0;iblk<nblksr;iblk++) {
		if (_bl2cpu[iblk] == _tree.myid) {
			for (int j=_iab[iblk];j<_iab[iblk+1];j++) {
				jasnd[j] = jabloc[nzjab];
				nzjab++;
			};
		};
	};

	if (_tree.nproc == 1) {
		for (i=0;i<nzjabtot;i++) _jab[i] = jabloc[i];
	} else {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													INTEGERVALUE, ADD,
													nzjabtot, jasnd, _jab);
	};

// Free working arrays

	delete [] lstblk;
	delete [] imask;
	delete [] sp2blk;
	delete [] nzarr;
	delete [] nzsum;
	delete [] iabloc;
	delete [] jabloc;
	delete [] iasnd;
	delete [] jasnd;

};

// Author: Kharchenko S.A.
// Description: Compute extended block sparsity with respect to triangular factorization
// CGSMatrix::ExtendBlockSparsity()
//========================================================================================
void CGSMatrix::ExtendBlockSparsity (int _nblks, // Compute extended block sparsity with respect to triangular factorization
										int *_iab, int *_jab, int *&_iabext, int *&_jabext) const {

	const char *funcname = "ExtendBlockSparsity";

// Allocate and init working arrays

	int *lstblk, *lstblk1;

	lstblk = new int [_nblks];
	if (!lstblk) MemoryFail (funcname);
	lstblk1 = new int [_nblks];
	if (!lstblk1) MemoryFail (funcname);

// Count the number of elements

	_iabext = new int [_nblks+1];
	if (!_iabext) MemoryFail (funcname);

	_iabext[0] = 0;

	int nzjab = 0;
	int nlstblk = 0;
	int nlstblk1 = 0;

	int iblk;

	for (iblk=0;iblk<_nblks;iblk++) {

		int nzjabloc = _iab[iblk+1]-_iab[iblk];
		int ibs = _iab[iblk];

		AddLists (nlstblk, lstblk,
					nzjabloc, _jab+ibs,
					nlstblk1, lstblk1);

		nlstblk = nlstblk1-1;

		for (int j=0;j<nlstblk;j++) lstblk[j] = lstblk1[j+1];

		nzjab += nlstblk1;

		_iabext[iblk+1] = nzjab;

	};

// Store jabext

	_jabext = new int [nzjab];
	if (!_jabext) MemoryFail (funcname);

	nzjab = 0;
	nlstblk = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		int nzjabloc = _iab[iblk+1]-_iab[iblk];
		int ibs = _iab[iblk];

		AddLists (nlstblk, lstblk,
					nzjabloc, _jab+ibs,
					nlstblk1, lstblk1);

		int j;

		for (j=0;j<nlstblk1;j++) _jabext[nzjab+j] = lstblk1[j];

		nlstblk = nlstblk1-1;

		for (j=0;j<nlstblk;j++) lstblk[j] = lstblk1[j+1];

		nzjab += nlstblk1;

		_iabext[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] lstblk1;

};

// Author: Kharchenko S.A.
// Description: Compute extended block sparsity with respect to triangular factorization
// CGSMatrix::ExtendBlockSparsityTri()
//========================================================================================
void CGSMatrix::ExtendBlockSparsityTri (int _nblks, // Compute extended block sparsity with respect to triangular factorization
										int *_iab, int *_jab, int *&_iabext, int *&_jabext) const {

	const char *funcname = "ExtendBlockSparsity";

// Allocate and init working arrays

	int *lstblk, *lstblk1;

	lstblk = new int [_nblks];
	if (!lstblk) MemoryFail (funcname);
	lstblk1 = new int [_nblks];
	if (!lstblk1) MemoryFail (funcname);

// Count the number of elements

	_iabext = new int [_nblks+1];
	if (!_iabext) MemoryFail (funcname);

	_iabext[0] = 0;

	int nzjab = 0;
	int nlstblk = 0;
	int nlstblk1 = 0;

	int iblk, j;

	for (iblk=0;iblk<_nblks;iblk++) {

		int nzjabloc = _iab[iblk+1]-_iab[iblk];
		int ibs = _iab[iblk];

		AddLists (nlstblk, lstblk,
					nzjabloc, _jab+ibs,
					nlstblk1, lstblk1);

		nlstblk = 0;

		for (j=0;j<nlstblk1;j++) {
			if (lstblk1[j] >= iblk) {
				lstblk[nlstblk] = lstblk1[j];
				nlstblk++;
			};
		};

		nzjab += nlstblk;

		_iabext[iblk+1] = nzjab;

	};

// Store jabext

	_jabext = new int [nzjab];
	if (!_jabext) MemoryFail (funcname);

	nzjab = 0;
	nlstblk = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		int nzjabloc = _iab[iblk+1]-_iab[iblk];
		int ibs = _iab[iblk];

		AddLists (nlstblk, lstblk,
					nzjabloc, _jab+ibs,
					nlstblk1, lstblk1);

		nlstblk = 0;

		for (j=0;j<nlstblk1;j++) {
			if (lstblk1[j] >= iblk) {
				lstblk[nlstblk] = lstblk1[j];
				_jabext[nzjab+nlstblk] = lstblk1[j];
				nlstblk++;
			};
		};

		nzjab += nlstblk;

		_iabext[iblk+1] = nzjab;

	};

// Free working arrays

	delete [] lstblk;
	delete [] lstblk1;

};

// Author: Kharchenko S.A.
// Description: Compute extended block sparsity with respect to triangular factorization
// CGSMatrix::ExtendBlockSparsityOpt()
//========================================================================================
void CGSMatrix::ExtendBlockSparsityOpt (int _nblks, // Compute extended block sparsity with respect to triangular factorization
										int *_iab, int *_jab, int *&_iabext, int *&_jabext) const {

	const char *funcname = "ExtendBlockSparsityOpt";

// Compute initial extended structure

	int *iabextini, *jabextini;

	ExtendBlockSparsityTri (_nblks,
							_iab, _jab, iabextini, jabextini);

// Create transposed structures

	int nzjaextini = iabextini[_nblks];

	int *iat, *jat, *ja2at, *jamask, *jatmask, *iptr, *listloc, *imask;

	iat = new int [_nblks+1];
	if (!iat) MemoryFail (funcname);
	jat = new int [nzjaextini];
	if (!jat) MemoryFail (funcname);
	ja2at = new int [nzjaextini];
	if (!ja2at) MemoryFail (funcname);
	jamask = new int [nzjaextini];
	if (!jamask) MemoryFail (funcname);
	jatmask = new int [nzjaextini];
	if (!jatmask) MemoryFail (funcname);
	iptr = new int [_nblks+1];
	if (!iptr) MemoryFail (funcname);
	listloc = new int [_nblks+1];
	if (!listloc) MemoryFail (funcname);
	imask = new int [_nblks+1];
	if (!imask) MemoryFail (funcname);

	int i, j, jj;

	for (i=0;i<=_nblks;i++) iat[i] = 0;

	for (i=0;i<_nblks;i++) {
		for (j=iabextini[i];j<iabextini[i+1];j++) {
			jj = jabextini[j];
			iat[jj+1]++;
		};
	};

	for (i=0;i<_nblks;i++) iat[i+1] = iat[i]+iat[i+1];

	for (i=0;i<_nblks;i++) iptr[i] = iat[i];

	int k;

	for (i=0;i<_nblks;i++) {
		for (j=iabextini[i];j<iabextini[i+1];j++) {
			jj = jabextini[j];
			k = iptr[jj];
			ja2at[j] = k;
			jat[k] = i;
			iptr[jj]++;
		};
	};

	for (i=0;i<nzjaextini;i++) jamask[i] = 0;
	for (i=0;i<nzjaextini;i++) jatmask[i] = 0;

	int icycle = -1;

	for (i=0;i<_nblks;i++) imask[i] = -1;

// Compute extended arrays

	_iabext = new int [_nblks+1];
	if (!_iabext) MemoryFail (funcname);

	int *jabloc;

	jabloc = new int [nzjaextini];
	if (!jabloc) MemoryFail (funcname);

	_iabext[0] = 0;

	int nzjab = 0;

	int nlistloc, kk, ind, indat;

	for (i=0;i<_nblks;i++) {

		icycle++;

// Init list by the current row

		nlistloc = 0;

		for (j=_iab[i];j<_iab[i+1];j++) {
			jj = _jab[j];
			if (jj >= i) {
				listloc[nlistloc] = jj;
				imask[jj] = icycle;
				nlistloc++;
			};
		};

// Update row by the previous nonzero ones

		for (j=iat[i];j<iat[i+1];j++) {
			jj = jat[j];
			if (jatmask[j] != 0) {
				for (k=iabextini[jj];k<iabextini[jj+1];k++) {
					kk = jabextini[k];
					if (jamask[k] != 0 && imask[kk] != icycle && kk > i) {
						listloc[nlistloc] = kk;
						imask[kk] = icycle;
						nlistloc++;
					};
				};
			};
		};

// Sort the row

		qsort (listloc, nlistloc, sizeof(int), compint);

// Mark nonzero data in the structures

		ind = 0;

		for (j=iabextini[i];j<iabextini[i+1];j++) {
			jj = jabextini[j];
			if (ind < nlistloc) {
				kk = listloc[ind];
				if (jj == kk) {
					jamask[j] = 1;
					indat = ja2at[j];
					jatmask[indat] = 1;
					ind++;
				};
			};
		};

// Store the result and mark nonzero data in the structures

		for (j=0;j<nlistloc;j++) jabloc[nzjab+j] = listloc[j];

		nzjab += nlistloc;

		_iabext[i+1] = nzjab;

	};

// Compress the result

	_jabext = new int [nzjab];
	if (!_jabext) MemoryFail (funcname);

	for (j=0;j<nzjab;j++) _jabext[j] = jabloc[j];

// Free working arrays

	delete [] iabextini;
	delete [] jabextini;
	delete [] iat;
	delete [] jat;
	delete [] ja2at;
	delete [] jamask;
	delete [] jatmask;
	delete [] iptr;
	delete [] listloc;
	delete [] imask;
	delete [] jabloc;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result
// CGSMatrix::CombineMatricesSort()
//========================================================================================
void CGSMatrix::CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												const CSMatrix *_mtrarr,
												CSMatrix &_mtrnew) {

	const char *funcname = "CombineMatricesSort";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrix temp(nlisttotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;
//	int nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;

// Reorder the rows

	CIntInt *iiarr;

	iiarr = new CIntInt [nlisttotal];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlisttotal;i++) {
		iiarr[i].intvalue = temp.list[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr,nlisttotal,sizeof(CIntInt),CIntInt::CompareIntInt);

	int *order, *iorder;

	order = new int [nlisttotal];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlisttotal];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		iorder[i] = iiarr[i].int2value;
	};

	for (i=0;i<nlisttotal;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *ialoc;
	int *jaloc;

	listloc = new int [nlisttotal];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlisttotal+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzjatotal];
	if (!jaloc) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		listloc[j] = temp.list[i];
		ialoc[j+1] = temp.ia[i+1]-temp.ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlisttotal;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		for (k=temp.ia[i];k<temp.ia[i+1];k++) {
			jaloc[k-temp.ia[i]+ialoc[j]] = temp.ja[k];
		};
	};

	for (i=0;i<nlisttotal;i++) temp.list[i] = listloc[i];
	for (i=0;i<=nlisttotal;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjatotal;i++) temp.ja[i] = jaloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;

	delete [] order;
	delete [] iorder;

// Return the result

	_mtrnew = temp;

	CSMatrix mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result (2index version)
// CGSMatrix::CombineMatricesSort2Index()
//========================================================================================
void CGSMatrix::CombineMatricesSort2Index (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result (2index version)
															const CSMatrix *_mtrarr,
															CSMatrix &_mtrnew) {

	const char *funcname = "CombineMatricesSort2Index";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrix temp(nlisttotal,nlisttotal,nzjatotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;
//	int nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list2[nlistloc+kii] = _mtrarr[iblk].list2[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja2[nzjaloc+kii] = _mtrarr[iblk].ja2[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nlist2 = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nzja2 = nzjatotal;

// Reorder the rows

	CInd2Int *iiarr;

	iiarr = new CInd2Int [nlisttotal];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlisttotal;i++) {
		iiarr[i].indx = temp.list2[i];
		iiarr[i].indy = temp.list[i];
		iiarr[i].intvalue = i;
	};

	qsort (iiarr,nlisttotal,sizeof(CInd2Int),CInd2Int::CompareInd2Int);

	int *order, *iorder;

	order = new int [nlisttotal];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlisttotal];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		iorder[i] = iiarr[i].intvalue;
	};

	for (i=0;i<nlisttotal;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *list2loc;
	int *ialoc;
	int *jaloc;
	int *ja2loc;

	listloc = new int [nlisttotal];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlisttotal];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlisttotal+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzjatotal];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjatotal];
	if (!ja2loc) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		listloc[j] = temp.list[i];
		list2loc[j] = temp.list2[i];
		ialoc[j+1] = temp.ia[i+1]-temp.ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlisttotal;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		for (k=temp.ia[i];k<temp.ia[i+1];k++) {
			jaloc[k-temp.ia[i]+ialoc[j]] = temp.ja[k];
			ja2loc[k-temp.ia[i]+ialoc[j]] = temp.ja2[k];
		};
	};

	for (i=0;i<nlisttotal;i++) temp.list[i] = listloc[i];
	for (i=0;i<nlisttotal;i++) temp.list2[i] = list2loc[i];
	for (i=0;i<=nlisttotal;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjatotal;i++) temp.ja[i] = jaloc[i];
	for (i=0;i<nzjatotal;i++) temp.ja2[i] = ja2loc[i];

	delete [] listloc;
	delete [] list2loc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] ja2loc;

	delete [] order;
	delete [] iorder;

// Return the result

	_mtrnew = temp;

	CSMatrix mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Symmetrize the matrix
// CGSMatrix::SymmetrizeMatrix2Index()
//========================================================================================
void CGSMatrix::SymmetrizeMatrix2Index (int _nblks, int *_blks, // Symmetrize the matrix
														int _nlist, int *_list, const CSMatrix *_mtrarr, 
														int &_nlistnew, int *_listnew, CSMatrix *_mtrarrnew) {

	const char *funcname = "SymmetrizeMatrix2Index";

// Check that matrix on entry satisfies assumptions

	int i, j, jjblk, iblk;

	for (i=0;i<_nlist-1;i++) {
		if (_list[i] > _list[i+1]) {
			throw " CGSMatrix::SymmetrizeMatrix2Index: nonmonotone list ";
		};
	};

	int ilist;
	int nlist1;
	int *plist, *plist2;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlist1 = _mtrarr[iblk].nlist-1;
		plist = _mtrarr[iblk].list;
		plist2 = _mtrarr[iblk].list2;
		for (i=0;i<nlist1;i++) {
			if (plist2[i] > plist2[i+1]) {
				throw " CGSMatrix::SymmetrizeMatrix2Index: nonmonotone list ";
			};
			if (plist2[i] == plist2[i+1]) {
				if (plist[i] > plist[i+1]) {
					throw " CGSMatrix::SymmetrizeMatrix2Index: nonmonotone list ";
				};
			};
		};
	};

// Create the total symmetrized list of matrix blocks

	int *imaskblk;

	imaskblk = new int [_nblks];
	if (!imaskblk) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<_nblks;i++) imaskblk[i] = icycle;

	icycle++;

	_nlistnew = 0;

	int nzjaloc;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		if (imaskblk[iblk] != icycle) {
			_listnew[_nlistnew] = iblk;
			_nlistnew++;
			imaskblk[iblk] = icycle;
		};
		nlist1 = _mtrarr[iblk].nlist;
		nzjaloc = _mtrarr[iblk].nzja;
		for (j=0;j<nzjaloc;j++) {
			jjblk = _mtrarr[iblk].ja2[j];
			if (imaskblk[jjblk] != icycle) {
				_listnew[_nlistnew] = jjblk;
				_nlistnew++;
				imaskblk[jjblk] = icycle;
			};
		};
	};

	qsort (_listnew, _nlistnew, sizeof(int), compint);

	for (i=0;i<_nlistnew;i++) {
		iblk = _listnew[i];
		imaskblk[iblk] = i;
	};

// Compute base and masks arrays

	int *ibsblks;

	ibsblks = new int [_nblks+1];
	if (!ibsblks) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) ibsblks[i] = -1;

	int nz = 0;

	for (i=0;i<_nlistnew;i++) {
		jjblk = _listnew[i];
		ibsblks[jjblk] = nz;
		nz += _blks[jjblk+1]-_blks[jjblk];
	};

// Compute estimate from above of the number of elements in each row

	int *ialoc;
	int *iptr;
	int *ialocnew;

	ialoc = new int [nz+1];
	if (!ialoc) MemoryFail (funcname);
	iptr = new int [nz+1];
	if (!iptr) MemoryFail (funcname);
	ialocnew = new int [nz+1];
	if (!ialocnew) MemoryFail (funcname);

	for (i=0;i<=nz;i++) ialoc[i] = 0;

	int ibs, jbs, iiblk, iirow, jjcol;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			iirow = _mtrarr[iblk].list[i];
			iiblk = _mtrarr[iblk].list2[i];
			ibs = ibsblks[iiblk];
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jjcol = _mtrarr[iblk].ja[j];
				jjblk = _mtrarr[iblk].ja2[j];
				jbs = ibsblks[jjblk];
				ialoc[ibs+iirow+1]++;
				ialoc[jbs+jjcol+1]++;
			};
		};
	};

	for (i=0;i<nz;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	for (i=0;i<nz;i++) iptr[i] = ialoc[i];

// Allocate and init data in each block row

	int **pjablk, **pja2blk;

	pjablk = new int * [_nlistnew];
	if (!pjablk) MemoryFail (funcname);
	pja2blk = new int * [_nlistnew];
	if (!pja2blk) MemoryFail (funcname);

	int niloc, nzj;

	for (ilist=0;ilist<_nlistnew;ilist++) {
		iblk = _listnew[ilist];
		ibs = ibsblks[iblk];
		niloc = _blks[iblk+1]-_blks[iblk];
		nzj = ialoc[ibs+niloc]-ialoc[ibs];
		pjablk[ilist] = new int [nzj];
		if (!pjablk[ilist]) MemoryFail (funcname);
		pja2blk[ilist] = new int [nzj];
		if (!pja2blk[ilist]) MemoryFail (funcname);
	};

	int ibs0, ind, ind0, jbs0, jind, ilisti, ilistj;
	int *jablkloc, *ja2blkloc;
	int *jjablkloc, *jja2blkloc;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			iirow = _mtrarr[iblk].list[i];
			iiblk = _mtrarr[iblk].list2[i];
			ibs0 = ibsblks[iiblk];
			ibs = ibsblks[iiblk]+iirow;
			ilisti = imaskblk[iiblk];
			jablkloc = pjablk[ilisti];
			ja2blkloc = pja2blk[ilisti];
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jjcol = _mtrarr[iblk].ja[j];
				jjblk = _mtrarr[iblk].ja2[j];
				ilistj = imaskblk[jjblk];
				jjablkloc = pjablk[ilistj];
				jja2blkloc = pja2blk[ilistj];
				jbs0 = ibsblks[jjblk];
				jbs = ibsblks[jjblk]+jjcol;
				ind = iptr[ibs]-ialoc[ibs0];
				jablkloc[ind] = jjcol;
				ja2blkloc[ind] = jjblk;
				iptr[ibs]++;
				jind = iptr[jbs]-ialoc[jbs0];
				jjablkloc[jind] = iirow;
				jja2blkloc[jind] = iiblk;
				iptr[jbs]++;
			};
		};
	};

// Perform sorting and filtering of the elements in each row

	int nzmax = 0;

	for (i=0;i<nz;i++) {
		niloc = ialoc[i+1]-ialoc[i];
		if (niloc > nzmax) nzmax = niloc;
	};

	CInd2Int *i2arr;

	i2arr = new CInd2Int [nzmax];
	if (!i2arr) MemoryFail (funcname);

	ialocnew[0] = 0;

	bool insert;

	int indnew, niblk;

	for (ilist=0;ilist<_nlistnew;ilist++) {
		iblk = _listnew[ilist];
		jablkloc = pjablk[ilist];
		ja2blkloc = pja2blk[ilist];
		niblk = _blks[iblk+1]-_blks[iblk];
		ibs0 = ibsblks[iblk];
		indnew = 0;
		for (i=0;i<niblk;i++) {
			ibs = ibsblks[iblk]+i;
			niloc = ialoc[ibs+1]-ialoc[ibs];
			ind0 = ialoc[ibs]-ialoc[ibs0];
			for (j=0;j<niloc;j++) {
				ind = ind0+j;
				iirow = jablkloc[ind];
				iiblk = ja2blkloc[ind];
				i2arr[j].indx = iiblk;
				i2arr[j].indy = iirow;
				i2arr[j].intvalue = j;
			};
			qsort (i2arr, niloc, sizeof (CInd2Int), CInd2Int::CompareInd2Int);
			for (j=0;j<niloc;j++) {
				insert = false;
				if (j==0) {
					insert = true;
				} else {
					if (i2arr[j].indx != i2arr[j-1].indx || i2arr[j].indy != i2arr[j-1].indy) insert = true;
				};
				if (insert) {
					iiblk = i2arr[j].indx;
					iirow = i2arr[j].indy;
					jablkloc[indnew] = iirow;
					ja2blkloc[indnew] = iiblk;
					indnew++;
				};
			};
			ialocnew[ibs+1] = ialocnew[ibs0]+indnew;
		};
	};

// Store the result

	int ntot = _blks[_nblks];

	int *pia, *pja, *pja2;

	for (ilist=0;ilist<_nlistnew;ilist++) {

		iblk = _listnew[ilist];
		jablkloc = pjablk[ilist];
		ja2blkloc = pja2blk[ilist];
		niblk = _blks[iblk+1]-_blks[iblk];
		ibs0 = ibsblks[iblk];
		nzjaloc = ialocnew[ibs0+niblk]-ialocnew[ibs0];

		CSMatrix aloc (niblk,niblk,nzjaloc,nzjaloc);

		plist = aloc.GetList ();
		plist2 = aloc.GetList2 ();
		pia = aloc.GetIa ();
		pja = aloc.GetJa ();
		pja2 = aloc.GetJa2 ();

		for (i=0;i<niblk;i++) {
			plist[i] = i;
			plist2[i] = iblk;
		};
		for (i=0;i<=niblk;i++) {
			pia[i] = ialocnew[ibs0+i]-ialocnew[ibs0];
		};
		for (i=0;i<nzjaloc;i++) {
			pja[i] = jablkloc[i];
			pja2[i] = ja2blkloc[i];
		};

		aloc.n = ntot;
		aloc.m = ntot;
		aloc.nsupr = ntot;
		aloc.nsupc = ntot;
		aloc.nlist = niblk;
		aloc.nlist2 = niblk;
		aloc.nzja = nzjaloc;
		aloc.nzja2 = nzjaloc;

		_mtrarrnew[iblk] = aloc;

	};

// Free work memory

	delete [] imaskblk;
	delete [] ibsblks;
	delete [] ialoc;
	delete [] iptr;
	delete [] ialocnew;

	for (i=0;i<_nlistnew;i++) {
		delete [] pjablk[i];
		delete [] pja2blk[i];
	};
	delete [] pjablk;
	delete [] pja2blk;
	delete [] i2arr;

};

// Author: Kharchenko S.A.
// Description: Output matrix
// CGSMatrix::operator<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CGSMatrix &_a) { // Output matrix

	_stream << " CGSMatrix:" << endl;
	_stream << " M = " << _a.m << " N = " << _a.n;
	_stream << " NblksR = " << _a.nblksr << " NblksC = " << _a.nblksc;
	_stream << " NsupR = " << _a.nsupr << " NsupC = " << _a.nsupc;
	_stream << " Nlist = " << _a.nlist << endl;

	OutArr (_stream, " BlksR =", _a.nblksr+1, _a.blksr);
	OutArr (_stream, " BlksC =", _a.nblksc+1, _a.blksc);
	OutArr (_stream, " SprndR =", _a.nsupr+1, _a.sprndr);
	OutArr (_stream, " SprndC =", _a.nsupc+1, _a.sprndc);
	OutArr (_stream, " Bl2ndR =", _a.nblksr+1, _a.bl2ndr);
	OutArr (_stream, " Bl2ndC =", _a.nblksc+1, _a.bl2ndc);
	OutArr (_stream, " Listb =", _a.nlist, _a.listb);

	return _stream;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CGSMatrixR::CGSMatrixR()
//========================================================================================
CGSMatrixR::CGSMatrixR (): CGSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixR_00";

	mtrarr = new CSMatrixRS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CGSMatrixR::CGSMatrixR()
//========================================================================================
CGSMatrixR::CGSMatrixR (int _nblks, int _nsups, int _nlist): CGSMatrix (_nblks, _nsups, _nlist) { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixR_01";

	mtrarr = new CSMatrixR [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destructor
// CGSMatrixR::~CGSMatrixR()
//========================================================================================
CGSMatrixR::~CGSMatrixR () { // Destructor
//	std::cout << " On entry to CGSMatrixR destructor " << std::endl;
	int i;
	for (i=0;i<nblksr;i++) {
		if (iblkcopymask[i] != 0) {
			RestoreSparsityPointers (i);
		};
	};
	delete [] mtrarr;
//	std::cout << " On return from CGSMatrixR destructor " << std::endl;
};

// Author: Kharchenko S.A.
// Description: Equality operator
// CGSMatrixR::operator=()
//========================================================================================
CGSMatrixR &CGSMatrixR::operator= (const CGSMatrixR &_aa) { // Equality operator

	const char *funcname = "CGSMatrixR_=";

	int i;
	for (i=0;i<nblksr;i++) {
		if (iblkcopymask[i] != 0) {
			RestoreSparsityPointers (i);
		};
	};

	delete [] mtrarr;

	this->CGSMatrix::operator= ((const CGSMatrix &) _aa);

	mtrarr = new CSMatrixR [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);
//	for (int iblk=0;iblk<nblksr;iblk++) {
//		mtrarr[iblk] = _aa.mtrarr[iblk];
//	};

	return *this;

};

// Author: Kharchenko S.A.
// Description: Write block row/column into the disk
// CGSMatrixR::WriteBlock()
//========================================================================================
void CGSMatrixR::WriteBlock (int _iblk) { // Write block row/column into the disk

//	const char *funcname = "WriteBlock";

	int nzaloc;
	int ibsloc;
	double *aloc;

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
// Description: Rewrite block row/column into the disk
// CGSMatrixR::ReWriteBlock()
//========================================================================================
void CGSMatrixR::ReWriteBlock (int _iblk) { // Rewrite block row/column into the disk

//	const char *funcname = "ReWriteBlock";

	int nzaloc;
	int ibsloc, ifileloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FPut (files[ifileloc], nzaloc, aloc, ibsloc);

Exit:;

};

// Author: Kharchenko S.A.
// Description: Read block row/column from the disk into the main memory if necessary
// CGSMatrixR::ReadBlock()
//========================================================================================
void CGSMatrixR::ReadBlock (int _iblk) { // Read block row/column from the disk into the main memory if necessary

	const char *funcname = "ReadBlock";

	int nzaloc;
	int ibsloc;
	int ifileloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FGet (files[ifileloc], nzaloc, aloc, ibsloc);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// Description: Allocate block row/column
// CGSMatrixR::AllocateBlock()
//========================================================================================
void CGSMatrixR::AllocateBlock (int _iblk) { // Allocate block row/column

	const char *funcname = "AllocateBlock";

	int nzaloc;
	double *aloc;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

};

// Author: Kharchenko S.A.
// Description: Free block row/column from main mamory
// CGSMatrixR::FreeBlock()
//========================================================================================
void CGSMatrixR::FreeBlock (int _iblk) { // Free block row/column from main mamory

	const char *funcname = "FreeBlock";

	int nzaloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	nzaloc = 0;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// Description: Filter block row/column according to the index
// CGSMatrixR::FilterBlock2Index()
//========================================================================================
void CGSMatrixR::FilterBlock2Index (int _iblk, int _iblkend) { // Filter block row/column according to the index

	const char *funcname = "FilterBlock2Index";

// Create new list of supernode rows

	int *listcount;

	int nlistloc = mtrarr[_iblk].nlist;

	listcount = new int [nlistloc];
	if (!listcount) MemoryFail (funcname);

	int ilist, j, jj, jj2, nz;

	for (ilist=0;ilist<nlistloc;ilist++) {
		nz = 0;
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj2 = mtrarr[_iblk].ja2[j];
			if (jj2>_iblkend) nz++;
		};
		listcount[ilist] = nz;
	};

	int nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) nlistnew++;
	};

	int *listnew, *list2new;

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

	jaloc = new int [nzjanew];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjanew];
	if (!ja2loc) MemoryFail (funcname);

	int nzjaloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj = mtrarr[_iblk].ja[j];
			jj2 = mtrarr[_iblk].ja2[j];
			int jjnew=jj;
			if (jjnew < 0) jjnew = -jjnew-1;
			if (jj2 > _iblkend) {
				jaloc[nzjaloc] = jj;
				ja2loc[nzjaloc] = jj2;
				nzjaloc++;
			};
		};
	};

// Create ja3 and a

	int *ja3loc;

	ja3loc = new int [nzjaloc];
	if (!ja3loc) MemoryFail (funcname);

	double *aloc;

	aloc = new double [nzjaloc];
	if (!aloc) MemoryFail (funcname);

	int nzaloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		for (j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			jj2 = mtrarr[_iblk].ja2[j];
			if (jj2 > _iblkend) {
				ja3loc[nzaloc] = mtrarr[_iblk].ja3[j];
				aloc[nzaloc] = mtrarr[_iblk].a[j];
				nzaloc++;
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
// Description: Transform secondary matrix into the hypersparse format
// CGSMatrixR::TransformIntoHyperSparseFormat()
//========================================================================================
void CGSMatrixR::TransformIntoHyperSparseFormat (CTree &_tree) { // Transform secondary matrix into the hypersparse format

	const char *funcname = "TransformIntoHyperSparseFormat";

// Get cpu ID

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Compute blk2cpu array

	int nblksloc = nblksr;

	int *blk2cpu;

	blk2cpu = new int [nblksloc];
	if (!blk2cpu) MemoryFail (funcname);

	_tree.Block2Cpu (blk2cpu);

// Two times scan local block rows and create the local lists

	int *nd2blk;

	nd2blk = new int [nblksloc+1];
	if (!nd2blk) MemoryFail (funcname);

	int i, iblk, nlistloc, iiblk, ilist, j, jj, jjblk;
	int *plist, *plist2, *pia, *pja, *pja2;

	for (i=0;i<nblksloc+1;i++) nd2blk[i] = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			nlistloc = mtrarr[iblk].GetNlist ();
			plist = mtrarr[iblk].GetList ();
			plist2 = mtrarr[iblk].GetList2 ();
			pia = mtrarr[iblk].GetIa ();
			pja = mtrarr[iblk].GetJa ();
			pja2 = mtrarr[iblk].GetJa2 ();
			for (ilist=0;ilist<nlistloc;ilist++) {
				i = plist[ilist];
				iiblk = plist2[ilist];
				nd2blk[iiblk+1]++;
				for (j=pia[ilist];j<pia[ilist+1];j++) {
					jj = pja[j];
					jjblk = pja2[j];
					nd2blk[jjblk+1]++;
				};
			};
		};
	};

	for (i=0;i<nblksloc;i++) nd2blk[i+1] = nd2blk[i]+nd2blk[i+1];

	int nzlist = nd2blk[nblksloc];
	int *iptr;
	int *listtot;

	iptr = new int [nblksloc+1];
	if (!iptr) MemoryFail (funcname);
	listtot = new int [nzlist];
	if (!listtot) MemoryFail (funcname);

	for (i=0;i<nblksloc;i++) iptr[i] = nd2blk[i];

	int k;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			nlistloc = mtrarr[iblk].GetNlist ();
			plist = mtrarr[iblk].GetList ();
			plist2 = mtrarr[iblk].GetList2 ();
			pia = mtrarr[iblk].GetIa ();
			pja = mtrarr[iblk].GetJa ();
			pja2 = mtrarr[iblk].GetJa2 ();
			for (ilist=0;ilist<nlistloc;ilist++) {
				i = plist[ilist];
				iiblk = plist2[ilist];
				k = iptr[iiblk];
				listtot[k] = i;
				iptr[iiblk]++;
				for (j=pia[ilist];j<pia[ilist+1];j++) {
					jj = pja[j];
					jjblk = pja2[j];
					k = iptr[jjblk];
					listtot[k] = jj;
					iptr[jjblk]++;
				};
			};
		};
	};

	int ni, ibs;

	for (i=0;i<nblksloc;i++) {
		ni = nd2blk[i+1]-nd2blk[i];
		ibs = nd2blk[i];
		if (ni > 0) qsort (listtot+ibs, ni, sizeof(int), compint);
	};

	iptr[0] = 0;

	int nz = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=nd2blk[iblk];j<nd2blk[iblk+1];j++) {
			if (j==nd2blk[iblk]) {
				listtot[nz] = listtot[j];
				nz++;
			} else {
				if (listtot[j] != listtot[nz]) {
					listtot[nz] = listtot[j];
					nz++;
				};
			};
		};
		iptr[iblk+1] = nz;
	};

	for (i=0;i<nblksloc+1;i++) nd2blk[i] = iptr[i];

// Exchange nd2blk arrays

	int *nd2blkcpu;

	nd2blkcpu = new int [nproc*(nblksloc+1)];
	if (!nd2blkcpu) MemoryFail (funcname);

	for (i=0;i<nproc*(nblksloc+1);i++) nd2blkcpu[i] = 0;
	for (i=0;i<nblksloc+1;i++) nd2blkcpu[i+myid*(nblksloc+1)] = nd2blk[i];

	CMPIComm pcomm = _tree.GetComm ();

	CMPIExchange::ExchangeArrayMPI (pcomm,
									INTEGERVALUE, ADD,
									nproc*(nblksloc+1), nd2blkcpu, nd2blkcpu);

// Exchange the lists among processors

	int NObjSend = nproc;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	int iproc;

	for (iproc=0;iproc<nproc;iproc++) {
		ObjTypeSend[iobj] = 1;
		ObjIDSend[iobj] = myid;
		CpuIDSend[iobj] = iproc;
		ObjSizeSend[iobj] = nz * sizeof(int);
		ObjSend[iobj] = (char *) (listtot);
		iobj++;
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info = CMPIExchange::DataExchangeMPI (pcomm,
												NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
												ObjSizeSend, ObjSend,
												NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
												ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixR::TransformIntoHyperSparseFormat";

// Free send arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Free local data array

	delete [] listtot;

// Scan received data and create new list array

	for (i=0;i<nblksloc+1;i++) nd2blk[i] = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (iproc=0;iproc<nproc;iproc++) {
			nd2blk[iblk+1] += (nd2blkcpu[iproc*(nblksloc+1)+iblk+1]-nd2blkcpu[iproc*(nblksloc+1)+iblk]);
		};
	};

	for (i=0;i<nblksloc;i++) nd2blk[i+1] = nd2blk[i]+nd2blk[i+1];

	nzlist = nd2blk[nblksloc];

	listtot = new int [nzlist];
	if (!listtot) MemoryFail (funcname);

	for (i=0;i<nblksloc;i++) iptr[i] = nd2blk[i];

	int iiproc;
	int *iarr;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (iiproc=0;iiproc<NObjRecv;iiproc++) {
			iproc = CpuIDRecv[iiproc];
			iarr = (int *) ObjRecv[iiproc];
			for (i=nd2blkcpu[iproc*(nblksloc+1)+iblk];i<nd2blkcpu[iproc*(nblksloc+1)+iblk+1];i++) {
				k = iptr[iblk];
				listtot[k] = iarr[i];
				iptr[iblk]++;
			};
		};
	};

	for (i=0;i<nblksloc;i++) {
		ni = nd2blk[i+1]-nd2blk[i];
		ibs = nd2blk[i];
		if (ni > 0) qsort (listtot+ibs, ni, sizeof(int), compint);
	};

	iptr[0] = 0;

	nz = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=nd2blk[iblk];j<nd2blk[iblk+1];j++) {
			if (j==nd2blk[iblk]) {
				listtot[nz] = listtot[j];
				nz++;
			} else {
				if (listtot[j] != listtot[nz]) {
					listtot[nz] = listtot[j];
					nz++;
				};
			};
		};
		iptr[iblk+1] = nz;
	};

	for (i=0;i<nblksloc+1;i++) nd2blk[i] = iptr[i];

// Free receive arrays

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Recompute the data

	delete [] listcnd;

	listcnd = new int [nz];
	if (!listcnd) MemoryFail (funcname);

	for (i=0;i<nz;i++) listcnd[i] = listtot[i];

	nlistcnd = nz;

	for (i=0;i<nblksr+1;i++) blkscnd[i] = nd2blk[i];

	int nlistl, ilstblk, ii, ii3;
	int *plist3, *pja3;

	for (ilstblk=0;ilstblk<nlist;ilstblk++) {

		iblk = listb[ilstblk];

		nlistl = mtrarr[iblk].GetNlist ();
		plist = mtrarr[iblk].GetList ();
		plist2 = mtrarr[iblk].GetList2 ();
		plist3 = mtrarr[iblk].GetList3 ();

		delete [] plist3;

		plist3 = new int [nlistl];
		if (!plist3) MemoryFail (funcname);

		for (i=0;i<nlistl;i++) {
			ii = plist[i];
			iiblk = plist2[i];
			ii3 = -1;
			for (j=blkscnd[iiblk];j<blkscnd[iiblk+1];j++) {
				if (listcnd[j] == ii) ii3 = j-blkscnd[iiblk];
			};
			if (ii3 < 0) throw " CGSMatrixR::TransformIntoHyperSparseFormat: index not found";
			plist3[i] = ii3;
		};

		mtrarr[iblk].SetNlist3 (nlistl);
		mtrarr[iblk].SetList3 (plist3);

		pia = mtrarr[iblk].GetIa ();
		pja = mtrarr[iblk].GetJa ();
		pja2 = mtrarr[iblk].GetJa2 ();
		pja3 = mtrarr[iblk].GetJa3 ();

		nz = pia[nlistl];

		delete [] pja3;

		pja3 = new int [nz];
		if (!pja3) MemoryFail (funcname);

		for (i=0;i<nz;i++) {
			ii = pja[i];
			iiblk = pja2[i];
			ii3 = -1;
			for (j=blkscnd[iiblk];j<blkscnd[iiblk+1];j++) {
				if (listcnd[j] == ii) ii3 = j-blkscnd[iiblk];
			};
			if (ii3 < 0) throw " CGSMatrixR::TransformIntoHyperSparseFormat: index not found";
			pja3[i] = ii3;
		};

		mtrarr[iblk].SetNzja3 (nz);
		mtrarr[iblk].SetJa3 (pja3);

	};

// Free work arrays

	delete [] blk2cpu;
	delete [] nd2blk;
	delete [] iptr;
	delete [] nd2blkcpu;
	delete [] listtot;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result
// CGSMatrixR::CombineMatricesSort()
//========================================================================================
void CGSMatrixR::CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
													const CSMatrixR *_mtrarr,
													CSMatrixR &_mtrnew) {

	const char *funcname = "CombineMatricesSort";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrixR temp(nlisttotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;
//	int nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.a[nzjaloc+kii] = _mtrarr[iblk].a[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;

// Reorder the rows

	CIntInt *iiarr;

	iiarr = new CIntInt [nlisttotal];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlisttotal;i++) {
		iiarr[i].intvalue = temp.list[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr,nlisttotal,sizeof(CIntInt),CIntInt::CompareIntInt);

	int *order, *iorder;

	order = new int [nlisttotal];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlisttotal];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		iorder[i] = iiarr[i].int2value;
	};

	for (i=0;i<nlisttotal;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *ialoc;
	int *jaloc;
	double *aloc;

	listloc = new int [nlisttotal];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlisttotal+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzjatotal];
	if (!jaloc) MemoryFail (funcname);
	aloc = new double [nzjatotal];
	if (!aloc) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		listloc[j] = temp.list[i];
		ialoc[j+1] = temp.ia[i+1]-temp.ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlisttotal;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		for (k=temp.ia[i];k<temp.ia[i+1];k++) {
			jaloc[k-temp.ia[i]+ialoc[j]] = temp.ja[k];
			aloc[k-temp.ia[i]+ialoc[j]] = temp.a[k];
		};
	};

	for (i=0;i<nlisttotal;i++) temp.list[i] = listloc[i];
	for (i=0;i<=nlisttotal;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjatotal;i++) temp.ja[i] = jaloc[i];
	for (i=0;i<nzjatotal;i++) temp.a[i] = aloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] aloc;

	delete [] order;
	delete [] iorder;

// Return the result

	_mtrnew = temp;

	CSMatrixR mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one
// CGSMatrixR::CombineMatrices2Index()
//========================================================================================
CSMatrixR CGSMatrixR::CombineMatrices2Index () { // Combine matrices into one

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		nlisttotal += mtrarr[iblk].nlist;
		nzjatotal += mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrixR temp(nlisttotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;
//	int nzaloc=0;

	temp.ia[0] = 0;

	int ishft, jblk, jshft;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		ishft = bl2ndr[iblk];
		for (kii=0;kii<mtrarr[iblk].nlist;kii++) {
			temp.list[nlistloc+kii] = ishft+mtrarr[iblk].list[kii];
		};
		for (kii=0;kii<mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (mtrarr[iblk].ia[kii+1] - mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<mtrarr[iblk].nzja;kii++) {
			jblk = mtrarr[iblk].ja2[kii];
			jshft = bl2ndr[jblk];
			temp.ja[nzjaloc+kii] = jshft+mtrarr[iblk].ja[kii];
			temp.a[nzjaloc+kii] = mtrarr[iblk].a[kii];
		};
		nlistloc += mtrarr[iblk].nlist;
		nzjaloc += mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = bl2ndr[nblksr];
	temp.n = bl2ndc[nblksc];
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;

	return temp;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one
// CGSMatrixR::CombineMatrices2Index()
//========================================================================================
CSMatrixR CGSMatrixR::CombineMatrices2IndexInto2Index () { // Combine matrices into one

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		nlisttotal += mtrarr[iblk].nlist;
		nzjatotal += mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrixR temp(nlisttotal,nlisttotal,nzjatotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;
//	int nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		for (kii=0;kii<mtrarr[iblk].nlist;kii++) {
			temp.list[nlistloc+kii] = mtrarr[iblk].list[kii];
			temp.list2[nlistloc+kii] = mtrarr[iblk].list2[kii];
		};
		for (kii=0;kii<mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (mtrarr[iblk].ia[kii+1] - mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<mtrarr[iblk].nzja;kii++) {
			temp.ja[nzjaloc+kii] = mtrarr[iblk].ja[kii];
			temp.ja2[nzjaloc+kii] = mtrarr[iblk].ja2[kii];
			temp.a[nzjaloc+kii] = mtrarr[iblk].a[kii];
		};
		nlistloc += mtrarr[iblk].nlist;
		nzjaloc += mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = bl2ndr[nblksr];
	temp.n = bl2ndc[nblksc];
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nlist2 = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nzja2 = nzjatotal;
	temp.nza = nzjatotal;
	temp.nzatot = nzjatotal;

	return temp;

};

// Author: Kharchenko S.A.
// Description: Reorder diagonal part of the matrix
// CGSMatrixR::ReorderDiagonalPart2Index()
//========================================================================================
void CGSMatrixR::ReorderDiagonalPart2Index (int _nblks, int *_blks, // Reorder diagonal part of the matrix
													int _nlist, int *_list, const CSMatrixR *_mtrarr, 
													int _iblkbeg, int _iblkend, int *_order,
													CSMatrixR *_mtrarrnew) {

	const char *funcname = "ReorderDiagonalPart2Index";

// Find block interval in the list

	int ilistbeg = -1;
	int ilistend = -1;

	int ilist;

	for (ilist=0;ilist<_nlist;ilist++) {
		if (_list[ilist] == _iblkbeg) ilistbeg = ilist;
		if (_list[ilist] == _iblkend) ilistend = ilist;
	};

	if (ilistbeg == -1 || ilistend == -1) throw " CGSMatrixR::ReorderDiagonalPart: diagonal block indices not found in the list ";

	if (ilistend-ilistbeg != _iblkend-_iblkbeg) throw " CGSMatrixR::ReorderDiagonalPart: not all blocks of the diagonal block are included into the list ";

// Check that all diagonal blocks have all rows

	int iblk;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {
		if (_mtrarr[iblk].nlist != _blks[iblk+1]-_blks[iblk]) throw " CGSMatrixR::ReorderDiagonalPart: incomplete diagonal block ";
	};

// Compute arrays that describe ordering

	int nz = 0;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {
		nz += _blks[iblk+1]-_blks[iblk];
	};

	int *iorder, *iind2blk, *iindarr, *iblkarr, *iindoldarr, *iblkoldarr;

	iorder = new int [nz];
	if (!iorder) MemoryFail (funcname);
	iind2blk = new int [nz];
	if (!iind2blk) MemoryFail (funcname);
	iindarr = new int [nz];
	if (!iindarr) MemoryFail (funcname);
	iblkarr = new int [nz];
	if (!iblkarr) MemoryFail (funcname);
	iindoldarr = new int [nz];
	if (!iindoldarr) MemoryFail (funcname);
	iblkoldarr = new int [nz];
	if (!iblkoldarr) MemoryFail (funcname);

	int i, j;

	for (i=0;i<nz;i++) iorder[_order[i]] = i;

	nz = 0;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			iind2blk[nz] = iblk;
			nz++;
		};
	};

	int inew, iblknew, indnew, iold, iblkold, indold;

	for (i=0;i<nz;i++) {
		inew = _order[i];
		iblknew = iind2blk[inew];
		indnew = inew - (_blks[iblknew]-_blks[_iblkbeg]);
		iindarr[i] = indnew;
		iblkarr[i] = iblknew;
		iold = iorder[i];
		iblkold = iind2blk[iold];
		indold = iold - (_blks[iblkold]-_blks[_iblkbeg]);
		iindoldarr[i] = indold;
		iblkoldarr[i] = iblkold;
	};

// Create the new ia array

	int *ialocnew;

	ialocnew = new int [nz+1];
	if (!ialocnew) MemoryFail (funcname);

	for (i=0;i<=nz;i++) ialocnew[i] = 0;

	int ind;

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {

		for (i=0;i<_mtrarr[iblk].nlist;i++) {

			ind = _blks[iblk]+i-_blks[_iblkbeg];
			indnew = _order[ind];
			ialocnew[indnew+1] = _mtrarr[iblk].ia[i+1]-_mtrarr[iblk].ia[i];

		};

	};

	for (i=0;i<nz;i++) ialocnew[i+1] = ialocnew[i]+ialocnew[i+1];

// Copy the data into the new matrix and reorder rows of the diagonal part

	for (ilist=0;ilist<ilistbeg;ilist++) {
		iblk = _list[ilist];
		_mtrarrnew[iblk] = _mtrarr[iblk];
	};
	for (ilist=ilistend+1;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		_mtrarrnew[iblk] = _mtrarr[iblk];
	};

	int indbeg, indend, iindold, niloc, nzjaloc;
	int *plistloc, *plist2loc, *pialoc, *pjaloc, *pja2loc;
	double *paloc;

	int ntot = _blks[_nblks];

	for (iblk=_iblkbeg;iblk<=_iblkend;iblk++) {

// Allocate the matrix

		indbeg = _blks[iblk]-_blks[_iblkbeg];
		indend = _blks[iblk+1]-1-_blks[_iblkbeg];

		niloc = _blks[iblk+1]-_blks[iblk];
		nzjaloc = ialocnew[indend+1]-ialocnew[indbeg];

		CSMatrixR aloc(niloc,niloc,nzjaloc,nzjaloc,nzjaloc);

		plistloc = aloc.GetList ();
		plist2loc = aloc.GetList2 ();
		pialoc = aloc.GetIa ();
		pjaloc = aloc.GetJa ();
		pja2loc = aloc.GetJa2 ();
		paloc = aloc.GetA ();

		for (i=0;i<niloc;i++) plistloc[i] = i;
		for (i=0;i<niloc;i++) plist2loc[i] = iblk;
		for (i=0;i<=niloc;i++) pialoc[i] = ialocnew[indbeg+i]-ialocnew[indbeg];

		nzjaloc = 0;

		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			ind = _blks[iblk]+i-_blks[_iblkbeg];
			indold = iorder[ind];
			iindold = iindoldarr[ind];
			iblkold = iblkoldarr[ind];
			for (j=_mtrarr[iblkold].ia[iindold];j<_mtrarr[iblkold].ia[iindold+1];j++) {
				pjaloc[nzjaloc] = _mtrarr[iblkold].ja[j];
				pja2loc[nzjaloc] = _mtrarr[iblkold].ja2[j];
				paloc[nzjaloc] = _mtrarr[iblkold].a[j];
				nzjaloc++;
			};
		};

		aloc.m = ntot;
		aloc.n = ntot;
		aloc.nsupr = ntot;
		aloc.nsupc = ntot;
		aloc.nlist = niloc;
		aloc.nlist2 = niloc;
		aloc.nzja = nzjaloc;
		aloc.nzja2 = nzjaloc;
		aloc.nza = nzjaloc;
		aloc.nzatot = nzjaloc;

		_mtrarrnew[iblk] = aloc;

	};

// Compute the maximal number of elements in the rows

	int nzjmax = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			niloc = _mtrarr[iblk].ia[i+1]-_mtrarr[iblk].ia[i];
			if (niloc > nzjmax) nzjmax = niloc;
		};
	};

// Reorder in-place columns in the whole list

	double *darr;
	CInd2Int *i2arr;

	darr = new double [nzjmax];
	if (!darr) MemoryFail (funcname);
	i2arr = new CInd2Int [nzjmax];
	if (!i2arr) MemoryFail (funcname);

	int ibeg, jj, jjblk, jjblknew, jjnew;
	int jloc, jlocold;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (i=0;i<_mtrarrnew[iblk].nlist;i++) {

			niloc = _mtrarrnew[iblk].ia[i+1]-_mtrarrnew[iblk].ia[i];

			ibeg = _mtrarrnew[iblk].ia[i];

// Modify column indices according to the ordering

			for (j=_mtrarrnew[iblk].ia[i];j<_mtrarrnew[iblk].ia[i+1];j++) {
				jloc = j-ibeg;
				jj = _mtrarrnew[iblk].ja[j];
				jjblk = _mtrarrnew[iblk].ja2[j];
				if (jjblk >= _iblkbeg && jjblk <= _iblkend) {
					ind = _blks[jjblk]+jj-_blks[_iblkbeg];
					indnew = _order[ind];
					jjblknew = iblkarr[ind];
					jjnew = iindarr[ind];
					_mtrarrnew[iblk].ja[j] = jjnew;
					_mtrarrnew[iblk].ja2[j] = jjblknew;
				};
			};

// Sort indices and values

			for (j=_mtrarrnew[iblk].ia[i];j<_mtrarrnew[iblk].ia[i+1];j++) {
				jloc = j-ibeg;
				jj = _mtrarrnew[iblk].ja[j];
				jjblk = _mtrarrnew[iblk].ja2[j];
				darr[jloc] = _mtrarrnew[iblk].a[j];
				i2arr[jloc].indx = jjblk;
				i2arr[jloc].indy = jj;
				i2arr[jloc].intvalue = jloc;
			};

			qsort (i2arr, niloc, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

			for (j=_mtrarrnew[iblk].ia[i];j<_mtrarrnew[iblk].ia[i+1];j++) {
				jloc = j-ibeg;
				jjblk = i2arr[jloc].indx;
				jj = i2arr[jloc].indy;
				jlocold = i2arr[jloc].intvalue;
				_mtrarrnew[iblk].ja[j]= jj;
				_mtrarrnew[iblk].ja2[j] = jjblk;
				_mtrarrnew[iblk].a[j] = darr[jlocold];
			};

		};

	};

// Free work arrays

	delete [] iorder;
	delete [] iind2blk;
	delete [] iindarr;
	delete [] iblkarr;
	delete [] iindoldarr;
	delete [] iblkoldarr;
	delete [] ialocnew;
	delete [] darr;
	delete [] i2arr;

};

// Author: Kharchenko S.A.
// Description: Split matrix into L and U parts
// CGSMatrixR::SplitLU2Index()
//========================================================================================
void CGSMatrixR::SplitLU2Index (int _nblks, int *_blks, // Split matrix into L and U parts
											int _nlist, int *_list, const CSMatrixR *_mtrarr, 
											CSMatrixR *_mtrlarr, CSMatrixR *_mtruarr) {

	const char *funcname = "SplitLU2Index";

// Count the number of rows/columns in the submatrix

	int ilist, iblk;
	int *ibsblk, *imaskblk;

	imaskblk = new int [_nblks];
	if (!imaskblk) MemoryFail (funcname);
	ibsblk = new int [_nblks];
	if (!ibsblk) MemoryFail (funcname);

	for (iblk=0;iblk<_nblks;iblk++) {
		imaskblk[iblk] = -1;
		ibsblk[iblk] = -1;
	};

	int nz = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		imaskblk[iblk] = ilist;
		ibsblk[iblk] = nz;
		nz += (_blks[iblk+1]-_blks[iblk]);
	};

// Count the numbers of elements in L and U parts

	int *ialoc;
	int *iptr;

	ialoc = new int [nz+1];
	if (!ialoc) MemoryFail (funcname);
	iptr = new int [nz];
	if (!iptr) MemoryFail (funcname);

	int i, j, jj, jjblk;

	for (i=0;i<=nz;i++) ialoc[i] = 0;

	int irow = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jj = _mtrarr[iblk].ja[j];
				jjblk = _mtrarr[iblk].ja2[j];
				if (jjblk > iblk || (jjblk == iblk && jj >= i)) {
					ialoc[irow+1]++;
				};
			};
			irow++;
		};
	};

	for (i=0;i<nz;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];
	for (i=0;i<nz;i++) iptr[i] = ialoc[i];

// Create the set of U block rows

	int niloc, nzjaloc;
	int *plistloc, *plist2loc, *pialoc, *pjaloc, *pja2loc;
	double *paloc;

	int ntot = _blks[_nblks];

	irow = 0;

	for (ilist=0;ilist<_nlist;ilist++) {

		iblk = _list[ilist];
		niloc = _blks[iblk+1]-_blks[iblk];
		nzjaloc = ialoc[irow+niloc]-ialoc[irow];

		CSMatrixR aloc(niloc,niloc,nzjaloc,nzjaloc,nzjaloc);

		plistloc = aloc.GetList ();
		plist2loc = aloc.GetList2 ();
		pialoc = aloc.GetIa ();
		pjaloc = aloc.GetJa ();
		pja2loc = aloc.GetJa2 ();
		paloc = aloc.GetA ();

		for (i=0;i<niloc;i++) plistloc[i] = i;
		for (i=0;i<niloc;i++) plist2loc[i] = iblk;
		for (i=0;i<=niloc;i++) pialoc[i] = ialoc[irow+i]-ialoc[irow];

		nzjaloc = 0;

		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jj = _mtrarr[iblk].ja[j];
				jjblk = _mtrarr[iblk].ja2[j];
				if (jjblk > iblk || (jjblk == iblk && jj >= i)) {
					pjaloc[nzjaloc] = jj;
					pja2loc[nzjaloc] = jjblk;
					paloc[nzjaloc] = _mtrarr[iblk].a[j];
					nzjaloc++;
				};
			};
			irow++;
		};

		aloc.m = ntot;
		aloc.n = ntot;
		aloc.nsupr = ntot;
		aloc.nsupc = ntot;
		aloc.nlist = niloc;
		aloc.nlist2 = niloc;
		aloc.nzja = nzjaloc;
		aloc.nzja2 = nzjaloc;
		aloc.nza = nzjaloc;
		aloc.nzatot = nzjaloc;

		_mtruarr[iblk] = aloc;

	};

// Init L data (sparsity structure and zero elements)

	int nzaloc;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];

		_mtrlarr[iblk] = _mtruarr[iblk];

		nzaloc = _mtrlarr[iblk].GetNza ();
		paloc = _mtrlarr[iblk].GetA ();

		for (i=0;i<nzaloc;i++) paloc[i] = 0.0e0;

	};

// Scan matrix data and fill in the elements of L

	int ibs, k;

	irow = 0;

	for (ilist=0;ilist<_nlist;ilist++) {

		iblk = _list[ilist];
		niloc = _blks[iblk+1]-_blks[iblk];
		nzjaloc = ialoc[irow+niloc]-ialoc[irow];

		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jj = _mtrarr[iblk].ja[j];
				jjblk = _mtrarr[iblk].ja2[j];
				if (jjblk < iblk || (jjblk == iblk && jj <= i)) {
					if (imaskblk[jjblk] < 0) {
						throw " CGSMatrixR::SplitLU2Index: block is not found in the list ";
					};
					ibs = ibsblk[jjblk];
					pjaloc = _mtrlarr[jjblk].GetJa ();
					pja2loc = _mtrlarr[jjblk].GetJa2 ();
					paloc = _mtrlarr[jjblk].GetA ();
					irow = ibs+jj;
					k = iptr[irow]-ialoc[ibs];
					if (pjaloc[k] != i || pja2loc[k] != iblk) {
						throw " CGSMatrixR::SplitLU2Index: indices not found in the sparsity structure ";
					};
					paloc[k] = _mtrarr[iblk].a[j];
					iptr[irow]++;
				};
			};
		};

	};

// Free work arrays

	delete [] imaskblk;
	delete [] ibsblk;
	delete [] ialoc;
	delete [] iptr;

};

// Author: Kharchenko S.A.
// Description: Output matrix
// CGSMatrixR::operator<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CGSMatrixR &_a) { // Output matrix

	_stream << " CGSMatrixR:" << endl;
	_stream << " Control = " << (CGSMatrix &)_a << endl;

	int nlist = _a.GetNlist();
	int *plist = _a.GetListb ();
	CSMatrixR *pmtrarr = _a.GetMtrarr ();

	for (int i=0;i<nlist;i++) {
		int iblk = plist[i];
		_stream << " Iblk = " << iblk << endl;
		_stream << " Block = " << pmtrarr[iblk] << endl;
	};

	return _stream;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CGSMatrixRS::CGSMatrixRS()
//========================================================================================
CGSMatrixRS::CGSMatrixRS (): CGSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixRS_00";

	mtrarr = new CSMatrixRS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

	ivect = -1;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CGSMatrixRS::CGSMatrixRS()
//========================================================================================
CGSMatrixRS::CGSMatrixRS (int _nblks, int _nsups, int _nlist): CGSMatrix (_nblks, _nsups, _nlist) { // Memory allocation zero data constructor

	const char *funcname = "CGSMatrixRS_01";

	mtrarr = new CSMatrixRS [nblksr+1];
	if (!mtrarr) MemoryFail (funcname);

	ivect = -1;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CGSMatrixRS::~CGSMatrixRS()
//========================================================================================
CGSMatrixRS::~CGSMatrixRS () { // Destructor
//	std::cout << " On entry to CGSMatrixRS destructor " << std::endl;
	delete [] mtrarr;
	ivect = -1;
	CSVector vectdummy;
	vector1 = vectdummy;
	vector2 = vectdummy;
//	std::cout << " On return from CGSMatrixRS destructor " << std::endl;
};

// Author: Kharchenko S.A.
// Description: Equality operator
// CGSMatrixRS::operator=()
//========================================================================================
CGSMatrixRS &CGSMatrixRS::operator= (const CGSMatrixRS &_aa) { // Equality operator

	const char *funcname = "CGSMatrixRS_=";

	delete [] mtrarr;

	this->CGSMatrix::operator= ((const CGSMatrix &) _aa);

	mtrarr = new CSMatrixRS [nblksr+1];
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
// Description: Filter block row/column according to the index
// CGSMatrixRS::FilterBlock()
//========================================================================================
void CGSMatrixRS::FilterBlock (int _iblk, int _isupend) { // Filter block row/column according to the index

	const char *funcname = "FilterBlock";

	int blamxloc = mtrarr[_iblk].blamx;

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
		int blai = blamxloc;
		if (blamxloc > 1) blai = sprndr[isup+1]-sprndr[isup];
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jsup = mtrarr[_iblk].ja[j];
			int jsupnew=jsup;
			if (jsupnew < 0) jsupnew = -jsupnew;
			int blaj = blamxloc;
			if (blamxloc > 1) blaj = sprndr[jsupnew+1]-sprndr[jsupnew];
			if (jsupnew>_isupend) {
				jaloc[nzjaloc] = jsup;
				bsaloc[nzjaloc] = nzaloc;
				nzjaloc++;
				nzaloc += blai*blaj;
			};
		};
	};

// Create a

	double *aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	nzaloc = 0;

	int ibs;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = mtrarr[_iblk].list[ilist];
		int blai = blamxloc;
		if (blamxloc > 1) blai = sprndr[isup+1]-sprndr[isup];
		for (int j=mtrarr[_iblk].ia[ilist];j<mtrarr[_iblk].ia[ilist+1];j++) {
			int jsup = mtrarr[_iblk].ja[j];
			int jsupnew=jsup;
			if (jsupnew < 0) jsupnew = -jsupnew;
			int blaj = blamxloc;
			if (blamxloc > 1) blaj = sprndr[jsupnew+1]-sprndr[jsupnew];
			if (jsupnew>_isupend) {
				if (blamxloc == 1) {
					ibs = j;
				} else {
					ibs = mtrarr[_iblk].bsa[j];
				};
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

	if (blamxloc == 1) mtrarr[_iblk].CleanBsa ();

// Free work memory

	delete [] listcount;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result
// CGSMatrixRS::CombineMatricesSort()
//========================================================================================
void CGSMatrixRS::CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												CSMatrixR *_mtrarr,
												CSMatrixR &_mtrnew) {

	const char *funcname = "CombineMatricesSort";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
	};

// Allocate and fill

	CSMatrixR temp(nlisttotal,nzjatotal);

	int nlistloc=0, nzjaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.a[nzjaloc+kii] = _mtrarr[iblk].a[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = nlisttotal;
	temp.nsupc = nlisttotal;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nza = nzjatotal;
	temp.nzatot = nzjatotal;

// Reorder the rows

	CIntInt *iiarr;

	iiarr = new CIntInt [nlisttotal];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlisttotal;i++) {
		iiarr[i].intvalue = temp.list[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr,nlisttotal,sizeof(CIntInt),CIntInt::CompareIntInt);

	int *order, *iorder;

	order = new int [nlisttotal];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlisttotal];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		iorder[i] = iiarr[i].int2value;
	};

	for (i=0;i<nlisttotal;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *ialoc;
	int *jaloc;
	double *aloc;

	listloc = new int [nlisttotal];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlisttotal+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzjatotal];
	if (!jaloc) MemoryFail (funcname);
	aloc = new double [nzjatotal];
	if (!aloc) MemoryFail (funcname);

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		listloc[j] = temp.list[i];
		ialoc[j+1] = temp.ia[i+1]-temp.ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlisttotal;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	int kloc;

	for (i=0;i<nlisttotal;i++) {
		j = order[i];
		for (k=temp.ia[i];k<temp.ia[i+1];k++) {
			kloc = k-temp.ia[i]+ialoc[j];
			jaloc[kloc] = temp.ja[k];
			aloc[kloc] = temp.a[k];
		};
	};

	for (i=0;i<nlisttotal;i++) temp.list[i] = listloc[i];
	for (i=0;i<=nlisttotal;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjatotal;i++) temp.ja[i] = jaloc[i];
	for (i=0;i<nzjatotal;i++) temp.a[i] = aloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] aloc;

	delete [] order;
	delete [] iorder;

// Return the result

	_mtrnew = temp;

	CSMatrixR mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result
// CGSMatrixRS::CombineMatricesSort()
//========================================================================================
void CGSMatrixRS::CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												int _nsup, int *_sprnds,
												CSMatrixRS *_mtrarr,
												CSMatrixRS &_mtrnew) {

	const char *funcname = "CombineMatricesSort";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
		nzatotal += _mtrarr[iblk].nza;
	};

// Allocate and fill

	CSMatrixRS temp(nlisttotal,nzjatotal,nzatotal);

	int nlistloc=0, nzjaloc=0, nzaloc=0;

	temp.ia[0] = 0;

	int i, irow, blai, j, jj, blaj;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nza;kii++) temp.a[nzaloc+kii] = _mtrarr[iblk].a[kii];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			irow = _mtrarr[iblk].list[i];
			blai = _sprnds[irow+1]-_sprnds[irow];
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jj = _mtrarr[iblk].ja[j];
				blaj = _sprnds[jj+1]-_sprnds[jj];
				temp.bsa[nzjaloc+j] = nzaloc;
				nzaloc += blai*blaj;
			};
		};
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	delete [] temp.sprndr;
	delete [] temp.sprndc;

	temp.sprndr = new int [_nsup+1];
	if (!temp.sprndr) MemoryFail (funcname);
	temp.sprndc = new int [_nsup+1];
	if (!temp.sprndc) MemoryFail (funcname);

	for (i=0;i<=_nsup;i++) temp.sprndr[i] = _sprnds[i];
	for (i=0;i<=_nsup;i++) temp.sprndc[i] = _sprnds[i];

	int blamxloc = 0;

	for (i=0;i<_nsup;i++) {
		blai = _sprnds[i+1]-_sprnds[i];
		if (blai > blamxloc) blamxloc = blai;
	};

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = _nsup;
	temp.nsupc = _nsup;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nza = nzatotal;
	temp.nzatot = nzatotal;
	temp.blamx = blamxloc;

// Reorder the rows

	temp.SortListAndColumns ();

// Return the result

	_mtrnew = temp;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one and sort the result
// CGSMatrixCS::CombineMatricesSort()
//========================================================================================
void CGSMatrixCS::CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												int _nsup, int *_sprnds,
												CSMatrixCS *_mtrarr,
												CSMatrixCS &_mtrnew) {

	const char *funcname = "CombineMatricesSort";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	int ilist, iblk, kii;

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

	int i, irow, blai, j, jj, blaj;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nza;kii++) temp.a[nzaloc+kii] = _mtrarr[iblk].a[kii];
		for (i=0;i<_mtrarr[iblk].nlist;i++) {
			irow = _mtrarr[iblk].list[i];
			blai = _sprnds[irow+1]-_sprnds[irow];
			for (j=_mtrarr[iblk].ia[i];j<_mtrarr[iblk].ia[i+1];j++) {
				jj = _mtrarr[iblk].ja[j];
				blaj = _sprnds[jj+1]-_sprnds[jj];
				temp.bsa[nzjaloc+j] = nzaloc;
				nzaloc += blai*blaj;
			};
		};
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
	};

// Fill control data

	delete [] temp.sprndr;
	delete [] temp.sprndc;

	temp.sprndr = new int [_nsup+1];
	if (!temp.sprndr) MemoryFail (funcname);
	temp.sprndc = new int [_nsup+1];
	if (!temp.sprndc) MemoryFail (funcname);

	for (i=0;i<=_nsup;i++) temp.sprndr[i] = _sprnds[i];
	for (i=0;i<=_nsup;i++) temp.sprndc[i] = _sprnds[i];

	int blamxloc = 0;

	for (i=0;i<_nsup;i++) {
		blai = _sprnds[i+1]-_sprnds[i];
		if (blai > blamxloc) blamxloc = blai;
	};

	temp.m = _n;
	temp.n = _n;
	temp.nsupr = _nsup;
	temp.nsupc = _nsup;
	temp.nlist = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nza = nzatotal;
	temp.nzatot = nzatotal;
	temp.blamx = blamxloc;

// Reorder the rows

	temp.SortListAndColumns ();

// Return the result

	_mtrnew = temp;

};

// Author: Kharchenko S.A.
// Description: Write block row/column into the disk
// CGSMatrixRS::WriteBlock()
//========================================================================================
void CGSMatrixRS::WriteBlock (int _iblk) { // Write block row/column into the disk

//	const char *funcname = "WriteBlock";

	int nzaloc;
	int ibsloc;
	double *aloc;

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
// Description: Rewrite block row/column into the disk
// CGSMatrixRS::ReWriteBlock()
//========================================================================================
void CGSMatrixRS::ReWriteBlock (int _iblk) { // Rewrite block row/column into the disk

//	const char *funcname = "ReWriteBlock";

	int nzaloc;
	int ibsloc, ifileloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FPut (files[ifileloc], nzaloc, aloc, ibsloc);

Exit:;

};

// Author: Kharchenko S.A.
// Description: Read block row/column from the disk into the main memory if necessary
// CGSMatrixRS::ReadBlock()
//========================================================================================
void CGSMatrixRS::ReadBlock (int _iblk) { // Read block row/column from the disk into the main memory if necessary

	const char *funcname = "ReadBlock";

	int nzaloc;
	int ibsloc;
	int ifileloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	ibsloc = bl2ibsfile[_iblk];
	ifileloc = bl2file[_iblk];

	FGet (files[ifileloc], nzaloc, aloc, ibsloc);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// Description: Allocate block row/column
// CGSMatrixRS::AllocateBlock()
//========================================================================================
void CGSMatrixRS::AllocateBlock (int _iblk) { // Allocate block row/column

	const char *funcname = "AllocateBlock";

	int nzaloc;
	double *aloc;

	nzaloc = mtrarr[_iblk].nzatot;
	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

};

// Author: Kharchenko S.A.
// Description: Free block row/column from main mamory
// CGSMatrixRS::FreeBlock()
//========================================================================================
void CGSMatrixRS::FreeBlock (int _iblk) { // Free block row/column from main mamory

	const char *funcname = "FreeBlock";

	int nzaloc;
	double *aloc;

	if (nfiles == 0) goto Exit;

	aloc = mtrarr[_iblk].a;

	delete [] aloc;

	nzaloc = 0;

	aloc = new double [nzaloc];
	if (!aloc) MemoryFail (funcname);

	mtrarr[_iblk].a = aloc;
	mtrarr[_iblk].nza = nzaloc;

Exit:;

};

// Author: Kharchenko S.A.
// Description: Supersparsify block rows
// CGSMatrixRS::SupersparsifyBlockRows()
//========================================================================================
void CGSMatrixRS::SupersparsifyBlockRows () { // Supersparsify block rows

//	const char *funcname = "SupersparsifyBlockRows";

	int i, iblk;

	for (i=0;i<nlist;i++) {
		iblk = listb[i];
		mtrarr[iblk].SuperSparsify ();
	};

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixRS structure via SMatrixRS
// CGSMatrixRS::RS2GRS()
//========================================================================================
void CGSMatrixRS::RS2GRS (const CSMatrixRS &_mtra, int _nblks, int *_blks) { // Create GSMatrixCS structure via SMatrixCS

	const char *funcname = "RS2GRS";

// Create temporary bl2cpu data

	int *bl2cpu;

	bl2cpu = new int [_nblks];
	if (!bl2cpu) MemoryFail (funcname);

	for (int iblk=0;iblk<_nblks;iblk++) bl2cpu[iblk] = 0;

	int myid = 0;

// Create the structure

	RS2GRS (true, _mtra, _nblks, _blks, myid, bl2cpu);

// Free work arrays

	delete [] bl2cpu;

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixR2Ind structure via SMatrixR for specified indices only
// CGSMatrixR::R2GR2Ind()
//========================================================================================
void CGSMatrixR::R2GR2Ind (const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only
							int _myid, int *_bl2cpu) {

	const char *funcname = "R2GR2Ind";

// Create node 2 blk array

	int *nd2blk;

	int nloc = _blks[_nblks];

	nd2blk = new int [nloc];
	if (!nd2blk) MemoryFail (funcname);

	int i, iblk;

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) nd2blk[i] = iblk;
	};

// Create temporary CGSMatrixR data

	CGSMatrixR tempgcs (_nblks, 0, _nblks);

	*this = tempgcs;

	m = _mtra.m;
	n = _mtra.n;
	nblksr = _nblks;
	nblksc = _nblks;
//	nsupr = _mtra.nsupr;
//	nsupc = _mtra.nsupc;
	nsupr = 0;
	nsupc = 0;

	int isup;

	for (iblk=0;iblk<=_nblks;iblk++) blksr[iblk] = _blks[iblk];
	for (iblk=0;iblk<=_nblks;iblk++) blksc[iblk] = _blks[iblk];

	bl2ndr[0] = 0;
	for (iblk=0;iblk<nblksr;iblk++) {
		int isupbeg=blksr[iblk];
		int isupend=blksr[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndr[iblk+1] = bl2ndr[iblk] + ni;
	};

	bl2ndc[0] = 0;
	for (iblk=0;iblk<nblksc;iblk++) {
		int isupbeg=blksc[iblk];
		int isupend=blksc[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndc[iblk+1] = bl2ndc[iblk] + ni;
	};

	nlist=0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_bl2cpu[iblk]==_myid) listb[nlist++] = iblk;
	};

// Create the set of CSMatrixCS block rows

	CSMatrixR mtrdummy;

	int isupbeg, isupend, nsupl;
	int nzjal, nzal, j, jsup, jblk;

	for (iblk=0;iblk<_nblks;iblk++) {

		if (_bl2cpu[iblk] == _myid) {

			isupbeg = _blks[iblk];
			isupend = _blks[iblk+1]-1;

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				for (j=_mtra.ia[isup];j<_mtra.ia[isup+1];j++) {
					jsup = _mtra.ja[j];
					nzjal++;
					nzal++;
				};
			};

// Allocate the data

			CSMatrixR temp (nsupl, nsupl, nzjal, nzjal, nzal);

// Init the matrix

			temp.m = _mtra.m;
			temp.n = _mtra.n;
			temp.nsupr = nsupl;
			temp.nsupc = nsupl;
			temp.nlist = nsupl;
			temp.nlist2 = nsupl;
			temp.nzja = nzjal;
			temp.nzja2 = nzjal;
			temp.nza = nzal;
			temp.nzatot = nzal;

			for (isup=isupbeg;isup<=isupend;isup++) {
				temp.list[isup-isupbeg] = isup-isupbeg;
				temp.list2[isup-isupbeg] = iblk;
			};

			temp.ia[0] = 0;

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				for (j=_mtra.ia[isup];j<_mtra.ia[isup+1];j++) {
					jsup = _mtra.ja[j];
					jblk = nd2blk[jsup];
					temp.ja[nzjal] = jsup-_blks[jblk];
					temp.ja2[nzjal] = jblk;
					temp.a[nzal] = _mtra.a[j];
					nzjal++;
					nzal++;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

// Store the matrix

			mtrarr[iblk] = temp;

			temp = mtrdummy;

		};

	};

// Free work arrays

	delete [] nd2blk;

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixR2Ind structure via SMatrixR for specified indices only, matrix contains only necessary rows
// CGSMatrixR::R2GRSubmatrix2Ind()
//========================================================================================
void CGSMatrixR::R2GRSubmatrix2Ind (const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only, matrix contains only necessary rows
							int _myid, int *_bl2cpu) {

	const char *funcname = "R2GRSubmatrix2Ind";

// Create node 2 blk array

	int *nd2blk;

	int nloc = _blks[_nblks];

	nd2blk = new int [nloc];
	if (!nd2blk) MemoryFail (funcname);

	int i, iblk;

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) nd2blk[i] = iblk;
	};

// Create temporary CGSMatrixR data

	CGSMatrixR tempgcs (_nblks, 0, _nblks);

	*this = tempgcs;

	m = _mtra.m;
	n = _mtra.n;
	nblksr = _nblks;
	nblksc = _nblks;
//	nsupr = _mtra.nsupr;
//	nsupc = _mtra.nsupc;
	nsupr = 0;
	nsupc = 0;

	int isup;

	for (iblk=0;iblk<=_nblks;iblk++) blksr[iblk] = _blks[iblk];
	for (iblk=0;iblk<=_nblks;iblk++) blksc[iblk] = _blks[iblk];

	bl2ndr[0] = 0;
	for (iblk=0;iblk<nblksr;iblk++) {
		int isupbeg=blksr[iblk];
		int isupend=blksr[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndr[iblk+1] = bl2ndr[iblk] + ni;
	};

	bl2ndc[0] = 0;
	for (iblk=0;iblk<nblksc;iblk++) {
		int isupbeg=blksc[iblk];
		int isupend=blksc[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndc[iblk+1] = bl2ndc[iblk] + ni;
	};

	nlist=0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_bl2cpu[iblk]==_myid) listb[nlist++] = iblk;
	};

// Create the set of CSMatrixCS block rows

	CSMatrixR mtrdummy;

	int isupbeg, isupend, nsupl;
	int nzjal, nzal, j, jsup, jblk, isup0;

	int isupbeg0 = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		if (_bl2cpu[iblk] == _myid) {

			isupbeg = _blks[iblk];
			isupend = _blks[iblk+1]-1;

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isupbeg0 + (isup-isupbeg);
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
					jsup = _mtra.ja[j];
					nzjal++;
					nzal++;
				};
			};

// Allocate the data

			CSMatrixR temp (nsupl, nsupl, nzjal, nzjal, nzal);

// Init the matrix

			temp.m = _mtra.m;
			temp.n = _mtra.n;
			temp.nsupr = _mtra.m;
			temp.nsupc = _mtra.n;
			temp.nlist = nsupl;
			temp.nlist2 = nsupl;
			temp.nzja = nzjal;
			temp.nzja2 = nzjal;
			temp.nza = nzal;
			temp.nzatot = nzal;

			for (isup=isupbeg;isup<=isupend;isup++) {
				temp.list[isup-isupbeg] = isup-isupbeg;
				temp.list2[isup-isupbeg] = iblk;
			};

			temp.ia[0] = 0;

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isupbeg0 + (isup-isupbeg);
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
					jsup = _mtra.ja[j];
					jblk = nd2blk[jsup];
					temp.ja[nzjal] = jsup-_blks[jblk];
					temp.ja2[nzjal] = jblk;
					temp.a[nzal] = _mtra.a[j];
					nzjal++;
					nzal++;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

// Store the matrix

			mtrarr[iblk] = temp;

			temp = mtrdummy;

			isupbeg0 += nsupl;

		};

	};

// Free work arrays

	delete [] nd2blk;

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixR2Ind structure via SMatrixR for specified indices only, matrix contains only necessary rows
// CGSMatrixR::R2GRSubmatrix2IndFastSearch()
//========================================================================================
void CGSMatrixR::R2GRSubmatrix2IndFastSearch (bool _is_all_rows, const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only, matrix contains only necessary rows
							int _myid, int *_bl2cpu) {

//	const char *funcname = "R2GRSubmatrix2IndFastSearch";

// Create temporary CGSMatrixR data

	CGSMatrixR tempgcs (_nblks, 0, _nblks);

	*this = tempgcs;

	m = _mtra.m;
	n = _mtra.n;
	nblksr = _nblks;
	nblksc = _nblks;
//	nsupr = _mtra.nsupr;
//	nsupc = _mtra.nsupc;
	nsupr = 0;
	nsupc = 0;

	int iblk;
	int isup;

	for (iblk=0;iblk<=_nblks;iblk++) blksr[iblk] = _blks[iblk];
	for (iblk=0;iblk<=_nblks;iblk++) blksc[iblk] = _blks[iblk];

	bl2ndr[0] = 0;
	for (iblk=0;iblk<nblksr;iblk++) {
		int isupbeg=blksr[iblk];
		int isupend=blksr[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndr[iblk+1] = bl2ndr[iblk] + ni;
	};

	bl2ndc[0] = 0;
	for (iblk=0;iblk<nblksc;iblk++) {
		int isupbeg=blksc[iblk];
		int isupend=blksc[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		bl2ndc[iblk+1] = bl2ndc[iblk] + ni;
	};

	nlist=0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_bl2cpu[iblk]==_myid) listb[nlist++] = iblk;
	};

// Create the set of CSMatrixR block rows

	CSMatrixR mtrdummy;

	int isupbeg, isupend, nsupl, jblkprev, jbegblk, jendblk;
	int nzjal, nzal, j, jsup, jblk, isup0;

	int isupbeg0 = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		if (_is_all_rows) isupbeg0 = _blks[iblk];

		if (_bl2cpu[iblk] == _myid) {

			isupbeg = _blks[iblk];
			isupend = _blks[iblk+1]-1;

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isupbeg0 + (isup-isupbeg);
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
					jsup = _mtra.ja[j];
					nzjal++;
					nzal++;
				};
			};

// Allocate the data

			CSMatrixR temp (nsupl, nsupl, nzjal, nzjal, nzal);

// Init the matrix

			temp.m = _mtra.m;
			temp.n = _mtra.n;
			temp.nsupr = _mtra.m;
			temp.nsupc = _mtra.n;
			temp.nlist = nsupl;
			temp.nlist2 = nsupl;
			temp.nzja = nzjal;
			temp.nzja2 = nzjal;
			temp.nza = nzal;
			temp.nzatot = nzal;

			for (isup=isupbeg;isup<=isupend;isup++) {
				temp.list[isup-isupbeg] = isup-isupbeg;
				temp.list2[isup-isupbeg] = iblk;
			};

			temp.ia[0] = 0;

			nzjal = 0;
			nzal = 0;

			jblkprev = iblk;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isupbeg0 + (isup-isupbeg);
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
					jsup = _mtra.ja[j];
					if (jsup >= _blks[jblkprev] && jsup < _blks[jblkprev+1]) {
						jblk = jblkprev;
					} else {
						jbegblk = 0;
						jendblk = _nblks-1;
						while (jbegblk != jendblk) {
							jblk = (jbegblk+jendblk)/2;
							if (jsup >= _blks[jblk] && jsup < _blks[jblk+1]) {
								jbegblk = jblk;
								jendblk = jblk;
							} else if (jsup < _blks[jblk]) {
								jendblk = jblk-1;
							} else if (jsup >= _blks[jblk+1]) {
								jbegblk = jblk+1;
							};
						};
						jblk = jbegblk;
					};
					jblkprev = jblk;
					temp.ja[nzjal] = jsup-_blks[jblk];
					temp.ja2[nzjal] = jblk;
					temp.a[nzal] = _mtra.a[j];
					nzjal++;
					nzal++;
				};
				temp.ia[isup-isupbeg+1] = nzjal;
			};

// Store the matrix

			mtrarr[iblk] = temp;

			temp = mtrdummy;

			isupbeg0 += nsupl;

		};

	};

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixRS structure via SMatrixRS for specified indices only
// CGSMatrixRS::RS2GRS()
//========================================================================================
void CGSMatrixRS::RS2GRS (bool _is_all_rows, const CSMatrixRS &_mtra, int _nblks, int *_blks, // Create GSMatrixCS structure via SMatrixCS for specified indices only
							int _myid, int *_bl2cpu) {

//	const char *funcname = "RS2GRS";

// Create temporary CGSMatrixCS data

	CGSMatrixRS tempgcs (_nblks, _mtra.nsupr, _nblks);

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

	CSMatrixRS mtrdummy;

	int isupbeg, isupend, nsupl;
	int nzjal, nzal, blai, blaj, j, jsup, kii, isupbeg0, isup0;

	isupbeg0 = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		isupbeg = _blks[iblk];
		isupend = _blks[iblk+1]-1;

		if (_is_all_rows) isupbeg0 = isupbeg;

		if (_bl2cpu[iblk] == _myid) {

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isup-isupbeg+isupbeg0;
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
					jsup = _mtra.ja[j];
					blaj = _mtra.sprndr[jsup+1]-_mtra.sprndr[jsup];
					nzjal++;
					nzal += blai*blaj;
				};
			};

// Allocate the data

			CSMatrixRS temp (nsupl, nzjal, nzal);

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
				isup0 = isup-isupbeg+isupbeg0;
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
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

			isupbeg0 += nsupl;

		};

	};

};

// Author: Kharchenko S.A.
// Description: Create GSMatrixRS structure via SMatrixRS for specified indices only
// CGSMatrixCS::CS2GCS()
//========================================================================================
void CGSMatrixCS::CS2GCS (bool _is_all_rows, const CSMatrixCS &_mtra, int _nblks, int *_blks, // Create GSMatrixCS structure via SMatrixCS for specified indices only
							int _myid, int *_bl2cpu) {

//	const char *funcname = "RS2GRS";

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
	int nzjal, nzal, blai, blaj, j, jsup, kii, isupbeg0, isup0;

	isupbeg0 = 0;

	for (iblk=0;iblk<_nblks;iblk++) {

		isupbeg = _blks[iblk];
		isupend = _blks[iblk+1]-1;

		if (_is_all_rows) isupbeg0 = isupbeg;

		if (_bl2cpu[iblk] == _myid) {

			nsupl = isupend-isupbeg+1;

// Count the number of elements in the block row

			nzjal = 0;
			nzal = 0;

			for (isup=isupbeg;isup<=isupend;isup++) {
				isup0 = isup-isupbeg+isupbeg0;
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
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
				isup0 = isup-isupbeg+isupbeg0;
				blai = _mtra.sprndr[isup+1]-_mtra.sprndr[isup];
				for (j=_mtra.ia[isup0];j<_mtra.ia[isup0+1];j++) {
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

			isupbeg0 += nsupl;

		};

	};

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one
// CGSMatrixR::CombineMatrices2Index()
//========================================================================================
void CGSMatrixR::CombineMatrices2Index (bool _cpyindex3, int _nlist, int *_list, // Combine matrices into one
										const CSMatrixR *_mtrarr,
										CSMatrixR &_mtrnew) const {

//	const char *funcname = "CombineMatrices";

// Compute the sizes

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	int ilist, iblk, kii;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlisttotal += _mtrarr[iblk].nlist;
		nzjatotal += _mtrarr[iblk].nzja;
		nzatotal += _mtrarr[iblk].nza;
	};

// Allocate and fill

	CSMatrixR temp(nlisttotal,nlisttotal,nzjatotal,nzjatotal,nzatotal);

	temp.AllocateJa3 (nzjatotal);

	int nlistloc=0, nzjaloc=0, nzaloc=0;

	temp.ia[0] = 0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list[nlistloc+kii] = _mtrarr[iblk].list[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) temp.list2[nlistloc+kii] = _mtrarr[iblk].list2[kii];
		for (kii=0;kii<_mtrarr[iblk].nlist;kii++) {
			temp.ia[nlistloc+1+kii] = temp.ia[nlistloc+kii] + (_mtrarr[iblk].ia[kii+1] - _mtrarr[iblk].ia[kii]);
		};
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja[nzjaloc+kii] = _mtrarr[iblk].ja[kii];
		for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja2[nzjaloc+kii] = _mtrarr[iblk].ja2[kii];
		if (_cpyindex3) {
			for (kii=0;kii<_mtrarr[iblk].nzja;kii++) temp.ja3[nzjaloc+kii] = _mtrarr[iblk].ja3[kii];
		};
		for (kii=0;kii<_mtrarr[iblk].nza;kii++) temp.a[nzaloc+kii] = _mtrarr[iblk].a[kii];
		nlistloc += _mtrarr[iblk].nlist;
		nzjaloc += _mtrarr[iblk].nzja;
		nzaloc += _mtrarr[iblk].nza;
	};

// Fill control data

	temp.m = m;
	temp.n = n;
	temp.nsupr = m;
	temp.nsupc = n;
	temp.nlist = nlisttotal;
	temp.nlist2 = nlisttotal;
	temp.nzja = nzjatotal;
	temp.nzja2 = nzjatotal;
	temp.nzja3 = nzjatotal;
	temp.nzatot = nzatotal;
	temp.nza = nzatotal;

// Return the result

	_mtrnew = temp;

	CSMatrixR mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Combine matrices into one
// CGSMatrixRS::CombineMatrices()
//========================================================================================
void CGSMatrixRS::CombineMatrices (int _nlist, int *_list, // Combine matrices into one
									const CSMatrixRS *_mtrarr,
									CSMatrixRS &_mtrnew) const {

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

	CSMatrixRS temp(nlisttotal,nzjatotal,nzatotal);

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
		blai = sprndr[isup+1]-sprndr[isup];
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

	if (temp.blamx == 1) temp.CleanBsa ();

// Return the result

	_mtrnew = temp;

	CSMatrixRS mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Split matrix into the set of matrices
// CGSMatrixR::SplitMatrixList2Index()
//========================================================================================
void CGSMatrixR::SplitMatrixList2Index (int &_nlistblk, int *_listblk, // Split matrix into the set of matrices
														const CSMatrixR &_mtr, 
														CSMatrixR *_mtrarr) {

	const char *funcname = "SplitMatrixList2Index";

// Create the list

	int *imaskblk;

	imaskblk = new int [nblksr];
	if (!imaskblk) MemoryFail (funcname);

	int i, iblk;

	for (i=0;i<nblksr;i++) imaskblk[i] = -1;

	_nlistblk = 0;

	int icycle = -1;

	icycle++;

	for (i=0;i<_mtr.nlist;i++) {
		iblk = _mtr.list2[i];
		if (imaskblk[iblk] != icycle) {
			_listblk[_nlistblk] = iblk;
			_nlistblk++;
			imaskblk[iblk] = icycle;
		};
	};

	qsort (_listblk, _nlistblk, sizeof(int), compint);

// Compute numbers of list values and nonzero entries

	int *nlistblk, *nzjablk;

	nlistblk = new int [nblksr];
	if (!nlistblk) MemoryFail (funcname);
	nzjablk = new int [nblksr];
	if (!nzjablk) MemoryFail (funcname);

	for (i=0;i<_nlistblk;i++) {
		iblk = _listblk[i];
		nlistblk[iblk] = 0;
		nzjablk[iblk] = 0;
	};

	for (i=0;i<_mtr.nlist;i++) {
		iblk = _mtr.list2[i];
		nlistblk[iblk]++;
		nzjablk[iblk] += (_mtr.ia[i+1]-_mtr.ia[i]);
	};

// Check the data on entry

	int ilistblk, kii;

// Compute the set of matrices

	CSMatrixR mtrdummy;

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	for (ilistblk=0;ilistblk<_nlistblk;ilistblk++) {

		iblk = _listblk[ilistblk];

// Determine the sizes of the data to be stored

		int nlistloc = nlistblk[iblk];
		int nzjaloc = nzjablk[iblk];
		int nzaloc = nzjaloc;

// Allocate and fill

		CSMatrixR temp(nlistloc,nlistloc,nzjaloc,nzjaloc,nzaloc);

		temp.AllocateJa3 (nzjaloc);

		temp.ia[0] = 0;

		for (kii=0;kii<nlistloc;kii++) temp.list[kii] = _mtr.list[nlisttotal+kii];
		for (kii=0;kii<nlistloc;kii++) temp.list2[kii] = _mtr.list2[nlisttotal+kii];
		for (kii=0;kii<nlistloc;kii++) {
			temp.ia[kii+1] = temp.ia[kii] + (_mtr.ia[nlisttotal+kii+1] - _mtr.ia[nlisttotal+kii]);
		};
		for (kii=0;kii<nzjaloc;kii++) temp.ja[kii] = _mtr.ja[nzjatotal+kii];
		for (kii=0;kii<nzjaloc;kii++) temp.ja2[kii] = _mtr.ja2[nzjatotal+kii];
		for (kii=0;kii<nzjaloc;kii++) temp.ja3[kii] = _mtr.ja3[nzjatotal+kii];
		for (kii=0;kii<nzaloc;kii++) temp.a[kii] = _mtr.a[nzatotal+kii];

		nlisttotal += nlistloc;
		nzjatotal += nzjaloc;
		nzatotal += nzaloc;

// Fill control data

		temp.m = m;
		temp.n = n;
		temp.nsupr = nlistloc;
		temp.nsupc = nlistloc;
		temp.nlist = nlistloc;
		temp.nlist2 = nlistloc;
		temp.nzja = nzjaloc;
		temp.nzja2 = nzjaloc;
		temp.nzja3 = nzjaloc;
		temp.nzatot = nzaloc;
		temp.nza = nzaloc;

// Store the result

		_mtrarr[iblk] = temp;

		temp = mtrdummy;

	};

// Free work arrays

	delete [] imaskblk;
	delete [] nlistblk;
	delete [] nzjablk;

};

// Author: Kharchenko S.A.
// Description: Split matrix into the set of matrices according to the list
// CGSMatrixR::SplitMatrix2Index()
//========================================================================================
void CGSMatrixR::SplitMatrix2Index (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
								const CSMatrixR &_mtr, 
								CSMatrixR *_mtrarr) {

//	const char *funcname = "SplitMatrix2Index";

// Check the data on entry

	int j, isup, jsup, ilist, ilistblk, iblk, kii;
	int nlistloc=0, nzjaloc=0, nzaloc=0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlistloc += blksr[iblk+1]-blksr[iblk];
	};

	if (nlistloc != _mtr.nlist) {
		cout << " Error in the list of indices when splitting the matrix" << endl;
		throw " Error in the list of indices when splitting the matrix ";
	};

// Compute the set of matrices

	CSMatrixR mtrdummy;

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	for (ilistblk=0;ilistblk<_nlist;ilistblk++) {

		iblk = _list[ilistblk];

// Determine the sizes of the data to be stored

		nlistloc = blksr[iblk+1]-blksr[iblk];

		nzjaloc = _mtr.ia[nlisttotal+nlistloc]-_mtr.ia[nlisttotal];

		nzaloc = 0;

		for (ilist=nlisttotal;ilist<nlisttotal+nlistloc;ilist++) {
			isup = _mtr.list[ilist];
			for (j=_mtr.ia[ilist];j<_mtr.ia[ilist+1];j++) {
				jsup = _mtr.ja[j];
				if (jsup < 0) jsup = -jsup-1;
				nzaloc++;
			};
		};

// Allocate and fill

		CSMatrixR temp(nlistloc,nlistloc,nzjaloc,nzjaloc,nzaloc);

		temp.AllocateJa3 (nzjaloc);

		temp.ia[0] = 0;

		for (kii=0;kii<nlistloc;kii++) temp.list[kii] = _mtr.list[nlisttotal+kii];
		for (kii=0;kii<nlistloc;kii++) temp.list2[kii] = _mtr.list2[nlisttotal+kii];
		for (kii=0;kii<nlistloc;kii++) {
			temp.ia[kii+1] = temp.ia[kii] + (_mtr.ia[nlisttotal+kii+1] - _mtr.ia[nlisttotal+kii]);
		};
		for (kii=0;kii<nzjaloc;kii++) temp.ja[kii] = _mtr.ja[nzjatotal+kii];
		for (kii=0;kii<nzjaloc;kii++) temp.ja2[kii] = _mtr.ja2[nzjatotal+kii];
		for (kii=0;kii<nzjaloc;kii++) temp.ja3[kii] = _mtr.ja3[nzjatotal+kii];
		for (kii=0;kii<nzaloc;kii++) temp.a[kii] = _mtr.a[nzatotal+kii];

		nlisttotal += nlistloc;
		nzjatotal += nzjaloc;
		nzatotal += nzaloc;

// Fill control data

		temp.m = m;
		temp.n = n;
		temp.nsupr = nlistloc;
		temp.nsupc = nlistloc;
		temp.nlist = nlistloc;
		temp.nlist2 = nlistloc;
		temp.nzja = nzjaloc;
		temp.nzja2 = nzjaloc;
		temp.nzja3 = nzjaloc;
		temp.nzatot = nzaloc;
		temp.nza = nzaloc;

// Store the result

		_mtrarr[iblk] = temp;

		temp = mtrdummy;

	};

};

// Author: Kharchenko S.A.
// Description: Free pointers of the current block and copy them from given matrix
// CGSMatrixR::FreeAndCopySparsityPointers()
//========================================================================================
void CGSMatrixR::FreeAndCopySparsityPointers (int _iblk, CSMatrix &_mtr) { // Free pointers of the current block and copy them from given matrix

	const char *funcname = "FreeAndCopySparsityPointers";

	delete [] mtrarr[_iblk].list;
	delete [] mtrarr[_iblk].list2;
	delete [] mtrarr[_iblk].ia;
	delete [] mtrarr[_iblk].ja;
	delete [] mtrarr[_iblk].ja2;
	delete [] mtrarr[_iblk].ja3;

	ptrsparr[_iblk*6+0] = new int [0];
	if (!ptrsparr[_iblk*6+0]) MemoryFail (funcname);
	ptrsparr[_iblk*6+1] = new int [0];
	if (!ptrsparr[_iblk*6+1]) MemoryFail (funcname);
	ptrsparr[_iblk*6+2] = new int [0];
	if (!ptrsparr[_iblk*6+2]) MemoryFail (funcname);
	ptrsparr[_iblk*6+3] = new int [0];
	if (!ptrsparr[_iblk*6+3]) MemoryFail (funcname);
	ptrsparr[_iblk*6+4] = new int [0];
	if (!ptrsparr[_iblk*6+4]) MemoryFail (funcname);
	ptrsparr[_iblk*6+5] = new int [0];
	if (!ptrsparr[_iblk*6+5]) MemoryFail (funcname);

	mtrarr[_iblk].list = _mtr.list;
	mtrarr[_iblk].list2 = _mtr.list2;
	mtrarr[_iblk].ia = _mtr.ia;
	mtrarr[_iblk].ja = _mtr.ja;
	mtrarr[_iblk].ja2 = _mtr.ja2;
	mtrarr[_iblk].ja3 = _mtr.ja3;

	iblkcopymask[_iblk] = 1;

};

// Author: Kharchenko S.A.
// Description: Restore sparsity pointers of the current block
// CGSMatrixR::RestoreSparsityPointers()
//========================================================================================
void CGSMatrixR::RestoreSparsityPointers (int _iblk) { // Restore sparsity pointers of the current block

	mtrarr[_iblk].list  = ptrsparr[_iblk*6+0];
	mtrarr[_iblk].list2 = ptrsparr[_iblk*6+1];
	mtrarr[_iblk].ia    = ptrsparr[_iblk*6+2];
	mtrarr[_iblk].ja    = ptrsparr[_iblk*6+3];
	mtrarr[_iblk].ja2   = ptrsparr[_iblk*6+4];
	mtrarr[_iblk].ja3   = ptrsparr[_iblk*6+5];

	iblkcopymask[_iblk] = 0;

};

// Author: Kharchenko S.A.
// Description: Split matrix into the set of matrices according to the list
// CGSMatrixRS::SplitMatrix()
//========================================================================================
void CGSMatrixRS::SplitMatrix (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
								const CSMatrixRS &_mtr, 
								CSMatrixRS *_mtrarr) {

//	const char *funcname = "SplitMatrix";

// Check the data on entry

	int i, j, isup, jsup, blai, blaj, ilist, ilistblk, iblk, kii;
	int nlistloc=0, nzjaloc=0, nzaloc=0;

	for (ilist=0;ilist<_nlist;ilist++) {
		iblk = _list[ilist];
		nlistloc += blksr[iblk+1]-blksr[iblk];
	};

	if (nlistloc != _mtr.nlist) {
		cout << " Error in the list of indices when splitting the matrix" << endl;
		throw " Error in the list of indices when splitting the matrix ";
	};

// Compute the set of matrices

	int nlisttotal=0, nzjatotal=0, nzatotal=0;

	for (ilistblk=0;ilistblk<_nlist;ilistblk++) {

		iblk = _list[ilistblk];

// Determine the sizes of the data to be stored

		nlistloc = blksr[iblk+1]-blksr[iblk];

		nzjaloc = _mtr.ia[nlisttotal+nlistloc]-_mtr.ia[nlisttotal];

		nzaloc = 0;

		for (ilist=nlisttotal;ilist<nlisttotal+nlistloc;ilist++) {
			isup = _mtr.list[ilist];
			blai = sprndr[isup+1]-sprndr[isup];
			for (j=_mtr.ia[ilist];j<_mtr.ia[ilist+1];j++) {
				jsup = _mtr.ja[j];
				if (jsup < 0) jsup = -jsup;
				blaj = sprndr[jsup+1]-sprndr[jsup];
				nzaloc += blai*blaj;
			};
		};

// Allocate and fill

		CSMatrixRS temp(nlistloc,nzjaloc,nzaloc);

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
			blai = sprndr[isup+1]-sprndr[isup];
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

		if (temp.blamx == 1) temp.CleanBsa ();

// Store the result

		_mtrarr[iblk] = temp;

		CSMatrixRS mtrdummy;

		temp = mtrdummy;

	};

};

// Author: Kharchenko S.A.
// Description: Prepare the matrix for exchange
// CGSMatrixRS::PrepareSendMatrix()
//========================================================================================
void CGSMatrixRS::PrepareSendMatrix (int _inode, const CTree &_tree, // Prepare the matrix for exchange
										CGSMatrixRS &_gmtru2, const CSMatrixRS *_mtru2charr,
										CSMatrixRS &_mtru2send) const {

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

	CSMatrixRS temp(nlistloc,nzjaloc,nzaloc);

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

	CSMatrixRS mtrdummy;

	temp = mtrdummy;

};

// Author: Kharchenko S.A.
// Description: Copy a part of data of the block into the new one
// CGSMatrixRS::CopyBlk2Blk()
//========================================================================================
void CGSMatrixRS::CopyBlk2Blk (char _diatype, const CSMatrixRS &_a, CSMatrixRS &_anew) const { // Copy a part of data of the block into the new one

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
// CGSMatrixRS: Copy a part of data of the block into the new one
//========================================================================================
void CGSMatrixRS::CopyBlkT2Blk (char _diatype, int &_icycle, int *_imask, 
								const CSMatrixRS &_a, CSMatrixRS &_anew) const { // Copy a part of data of the block into the new one

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
// CGSMatrixRS: Write current block into full format
//========================================================================================
void CGSMatrixRS::FullFormat (int _iblk, int _ntotal, double *_a) { // Write current block into full format

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
// CGSMatrixRS: Compute the sparsity structure of all local block rows
//========================================================================================
CSMatrix CGSMatrixRS::GlobalSparsity () const { // Compute the sparsity structure of all local block rows

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
// CGSMatrixRS: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CGSMatrixRS &_a) { // Output matrix

	_stream << " CGSMatrixRS:" << endl;

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
