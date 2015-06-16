//------------------------------------------------------------------------------------------------
// File: fct.cpp
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

#include "fct.h"
#include "fctdiag.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFct::CFct()
//========================================================================================
CFct::CFct () { // Memory allocation zero data constructor

	const char *funcname = "CFct_00";

	nblks   = 0;
	nsups   = 0;
	icycle  = -1;
	icyclegl = -1;
	nlist   = 0;
	nupd    = 0;
	nzr     = 0;

	int i;

	ibsdia = new int [nblks+1];
	if (!ibsdia) MemoryFail (funcname);
	ibsimask = new int [nblks+1];
	if (!ibsimask) MemoryFail (funcname);
	ibsfmask = new int [nblks+1];
	if (!ibsfmask) MemoryFail (funcname);
	imask = new int [nsups+1];
	if (!imask) MemoryFail (funcname);
	imasklev = new int [nsups+1];
	if (!imasklev) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		imask[i] = -1;
	};
	lstloc = new int [nsups+1];
	if (!lstloc) MemoryFail (funcname);
	lstloc2 = new int [nsups+1];
	if (!lstloc2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		lstloc[i] = -1;
	};
	iv = new int [nsups+1];
	if (!iv) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		iv[i] = -1;
	};
	madj = new int [nsups+1];
	if (!madj) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj[i] = -1;
	};
	madj2 = new int [nsups+1];
	if (!madj2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj2[i] = -1;
	};
	ibegm = new int [nsups+1];
	if (!ibegm) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm[i] = -1;
	};
	ibegm2 = new int [nsups+1];
	if (!ibegm2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm2[i] = -1;
	};
	ibsblk = new int [nblks+1];
	if (!ibsblk) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		ibsblk[i] = -1;
	};
	irw2lst = new int [nsups+1];
	if (!irw2lst) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		irw2lst[i] = -1;
	};
	ivgl = new int [nsups+1];
	if (!ivgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ivgl[i] = -1;
	};
	madjgl = new int [nsups+1];
	if (!madjgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madjgl[i] = -1;
	};
	madj2gl = new int [nsups+1];
	if (!madj2gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj2gl[i] = -1;
	};
	ibegmgl = new int [nsups+1];
	if (!ibegmgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegmgl[i] = -1;
	};
	ibegm2gl = new int [nsups+1];
	if (!ibegm2gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm2gl[i] = -1;
	};
	ibsblkgl = new int [nblks+1];
	if (!ibsblkgl) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		ibsblkgl[i] = -1;
	};
	imaskgl = new int [nsups+1];
	if (!imaskgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		imaskgl[i] = -1;
	};
	lstupd = new int [nsups+1];
	if (!lstupd) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		lstupd[i] = -1;
	};
	lstupd2 = new int [nsups+1];
	if (!lstupd2) MemoryFail (funcname);
	bsdia = new int [nsups+1];
	if (!bsdia) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		bsdia[i] = -1;
	};
	adrfmask = new int [nsups+1];
	if (!adrfmask) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		adrfmask[i] = -1;
	};
	ig = new int [nsups+1];
	if (!ig) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ig[i] = -1;
	};
	jg = new int * [nsups+1];
	if (!jg) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		jg[i] = 0;
	};
	jg2 = new int * [nsups+1];
	if (!jg2) MemoryFail (funcname);
	jg3 = new int * [nsups+1];
	if (!jg3) MemoryFail (funcname);
	addrg = new int * [nsups+1];
	if (!addrg) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		addrg[i] = 0;
	};
	diarr = new CDoubleInt [nsups+1];
	if (!diarr) MemoryFail (funcname);
	ind2arr = new CInd2Int [nsups+1];
	if (!ind2arr) MemoryFail (funcname);
	ind2arrgl = new CInd2Int [nsups+1];
	if (!ind2arrgl) MemoryFail (funcname);

	nmodsv = 0;
	ops = 0.0e0;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFct::CFct()
//========================================================================================
CFct::CFct (int _nblks, int _nsups) { // Memory allocation zero data constructor

	const char *funcname = "CFct_01";

	nblks   = _nblks;
	nsups   = _nsups;
	icycle  = -1;
	icyclegl = -1;
	nlist   = 0;
	nupd    = 0;
	nzr     = 0;

	int i;

	ibsdia = new int [nblks+1];
	if (!ibsdia) MemoryFail (funcname);
	ibsimask = new int [nblks+1];
	if (!ibsimask) MemoryFail (funcname);
	ibsfmask = new int [nblks+1];
	if (!ibsfmask) MemoryFail (funcname);

	imask = new int [nsups+1];
	if (!imask) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		imask[i] = -1;
	};
	imasklev = new int [nsups+1];
	if (!imasklev) MemoryFail (funcname);
	lstloc = new int [nsups+1];
	if (!lstloc) MemoryFail (funcname);
	lstloc2 = new int [nsups+1];
	if (!lstloc2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		lstloc[i] = -1;
	};
	iv = new int [nsups+1];
	if (!iv) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		iv[i] = -1;
	};
	madj = new int [nsups+1];
	if (!madj) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj[i] = -1;
	};
	madj2 = new int [nsups+1];
	if (!madj2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj2[i] = -1;
	};
	ibegm = new int [nsups+1];
	if (!ibegm) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm[i] = -1;
	};
	ibegm2 = new int [nsups+1];
	if (!ibegm2) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm2[i] = -1;
	};
	ibsblk = new int [nblks+1];
	if (!ibsblk) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		ibsblk[i] = -1;
	};
	irw2lst = new int [nsups+1];
	if (!irw2lst) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		irw2lst[i] = -1;
	};
	ivgl = new int [nsups+1];
	if (!ivgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ivgl[i] = -1;
	};
	madjgl = new int [nsups+1];
	if (!madjgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madjgl[i] = -1;
	};
	madj2gl = new int [nsups+1];
	if (!madj2gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		madj2gl[i] = -1;
	};
	ibegmgl = new int [nsups+1];
	if (!ibegmgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegmgl[i] = -1;
	};
	ibegm2gl = new int [nsups+1];
	if (!ibegm2gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ibegm2gl[i] = -1;
	};
	ibsblkgl = new int [nblks+1];
	if (!ibsblkgl) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		ibsblkgl[i] = -1;
	};
	imaskgl = new int [nsups+1];
	if (!imaskgl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		imaskgl[i] = -1;
	};
	lstupd = new int [nsups+1];
	if (!lstupd) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		lstupd[i] = -1;
	};
	lstupd2 = new int [nsups+1];
	if (!lstupd2) MemoryFail (funcname);
	bsdia = new int [nsups+1];
	if (!bsdia) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		bsdia[i] = -1;
	};
	adrfmask = new int [nsups+1];
	if (!adrfmask) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		adrfmask[i] = -1;
	};
	ig = new int [nsups+1];
	if (!ig) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		ig[i] = -1;
	};
	jg = new int * [nsups+1];
	if (!jg) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		jg[i] = 0;
	};
	jg2 = new int * [nsups+1];
	if (!jg2) MemoryFail (funcname);
	jg3 = new int * [nsups+1];
	if (!jg3) MemoryFail (funcname);
	addrg = new int * [nsups+1];
	if (!addrg) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) {
		addrg[i] = 0;
	};
	diarr = new CDoubleInt [nsups+1];
	if (!diarr) MemoryFail (funcname);
	ind2arr = new CInd2Int [nsups+1];
	if (!ind2arr) MemoryFail (funcname);
	ind2arrgl = new CInd2Int [nsups+1];
	if (!ind2arrgl) MemoryFail (funcname);

	nmodsv = 0;
	ops = 0.0e0;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFct::~CFct()
//========================================================================================
CFct::~CFct () { // Destructor
//	cout << " On entry to CFct destructor " << endl;
	nblks   = 0;
	nsups   = 0;
	icycle  = -1;
	icyclegl = -1;
	nlist   = 0;
	nupd    = 0;
	nzr     = 0;
	delete [] ibsdia;
	delete [] ibsimask;
	delete [] ibsfmask;
	delete [] imask;
	delete [] imasklev;
	delete [] lstloc;
	delete [] lstloc2;
	delete [] iv;
	delete [] madj;
	delete [] madj2;
	delete [] ibegm;
	delete [] ibegm2;
	delete [] ibsblk;
	delete [] irw2lst;
	delete [] ivgl;
	delete [] madjgl;
	delete [] madj2gl;
	delete [] ibegmgl;
	delete [] ibegm2gl;
	delete [] ibsblkgl;
	delete [] imaskgl;
	delete [] lstupd;
	delete [] lstupd2;
	delete [] bsdia;
	delete [] adrfmask;
	delete [] ig;
	delete [] jg;
	delete [] jg2;
	delete [] jg3;
	delete [] addrg;
	delete [] diarr;
	delete [] ind2arr;
	delete [] ind2arrgl;
	nmodsv = 0;
	ops = 0.0e0;
//	cout << " On return from CFct destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctR::CFctR()
//========================================================================================
CFctR::CFctR (int _nblks, int _nsups, int _blamx): CFct (_nblks,_nsups) { // Memory allocation zero data constructor

	const char *funcname = "CFctR_01";

	blamx = _blamx;
	lwork = blamx * 10;

	int i;

	iblkalloc = -1;
	nzused = 0;
	nzalloc = 0;
	nzallocmax = 0;

	jgblk = new int * [nsups+1];
	if (!jgblk) MemoryFail (funcname);
	jg2blk = new int * [nsups+1];
	if (!jg2blk) MemoryFail (funcname);
	jg3blk = new int * [nsups+1];
	if (!jg3blk) MemoryFail (funcname);
	adrgblk = new int * [nsups+1];
	if (!adrgblk) MemoryFail (funcname);
	glblk = new double * [nsups+1];
	if (!glblk) MemoryFail (funcname);
	gublk = new double * [nsups+1];
	if (!gublk) MemoryFail (funcname);
	dpiv = new double [1];
	if (!dpiv) MemoryFail (funcname);
	sclpt = new double [1];
	if (!sclpt) MemoryFail (funcname);
	scll = new double [1];
	if (!scll) MemoryFail (funcname);
	sclu = new double [1];
	if (!sclu) MemoryFail (funcname);
	sclinvl = new double [1];
	if (!sclinvl) MemoryFail (funcname);
	sclinvu = new double [1];
	if (!sclinvu) MemoryFail (funcname);
	diagg = new double [1];
	if (!diagg) MemoryFail (funcname);
	fmaskl = new double [1];
	if (!fmaskl) MemoryFail (funcname);
	fmasku = new double [1];
	if (!fmasku) MemoryFail (funcname);

	adrupdl = new double * [nsups+1];
	if (!adrupdl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) adrupdl[i] = 0;
	adrupdu = new double * [nsups+1];
	if (!adrupdu) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) adrupdu[i] = 0;
	gl = new double * [nsups+1];
	if (!gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) gl[i] = 0;
	gu = new double * [nsups+1];
	if (!gu) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) gu[i] = 0;

	int blamx_2 = blamx*blamx;

	lelem = new double [blamx_2];
	if (!lelem) MemoryFail (funcname);
	uelem = new double [blamx_2];
	if (!uelem) MemoryFail (funcname);
	deleml = new double [blamx_2];
	if (!deleml) MemoryFail (funcname);
	delemu = new double [blamx_2];
	if (!delemu) MemoryFail (funcname);
	aloc = new double [blamx_2];
	if (!aloc) MemoryFail (funcname);
	uloc = new double [blamx_2];
	if (!uloc) MemoryFail (funcname);
	vloc = new double [blamx_2];
	if (!vloc) MemoryFail (funcname);
	eig = new double [blamx];
	if (!eig) MemoryFail (funcname);
	work = new double [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctR::~CFctR()
//========================================================================
CFctR::~CFctR () { // Destructor
//	cout << " On entry to CFctR destructor " << endl;
	blamx  = 0;
	lwork  = 0;
	iblkalloc = -1;
	nzused = 0;
	nzalloc = 0;
	nzallocmax = 0;

	delete [] jgblk;
	delete [] jg2blk;
	delete [] jg3blk;
	delete [] adrgblk;
	delete [] glblk;
	delete [] gublk;
	delete [] dpiv;
	delete [] sclpt;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] diagg;
	delete [] fmaskl;
	delete [] fmasku;
	delete [] adrupdl;
	delete [] adrupdu;
	delete [] gl;
	delete [] gu;
	delete [] lelem;
	delete [] uelem;
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] eig;
	delete [] work;
	delete [] rwork;
//	cout << " On return from CFctR destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctDiagR::CFctDiagR()
//========================================================================================
CFctDiagR::CFctDiagR (int _nblks) { // Memory allocation zero data constructor

	const char *funcname = "CFctDiagR_00";

	nblks  = _nblks;
	nfiles = 0;
	ifile = 0;

	int i;

	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	blk2file = new int [nblks*5];
	if (!blk2file) MemoryFail (funcname);
	for (i=0;i<nblks*5;i++) {
		blk2file[i] = -1;
	};
	blk2bs = new int [nblks*5];
	if (!blk2bs) MemoryFail (funcname);
	for (i=0;i<nblks*5;i++) {
		blk2bs[i] = -1;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);

	diasize = new int [nblks+1];
	if (!diasize) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		diasize[i] = 0;
	};
	sclpt = new double * [nblks];
	if (!sclpt) MemoryFail (funcname);
	scll = new double * [nblks];
	if (!scll) MemoryFail (funcname);
	sclu = new double * [nblks];
	if (!sclu) MemoryFail (funcname);
	sclinvl = new double * [nblks];
	if (!sclinvl) MemoryFail (funcname);
	sclinvu = new double * [nblks];
	if (!sclinvu) MemoryFail (funcname);
	diagg = new double * [nblks];
	if (!diagg) MemoryFail (funcname);
	dpiv = new double * [nblks];
	if (!dpiv) MemoryFail (funcname);

	for (i=0;i<nblks;i++) {
		sclpt[i] = new double [0];
		if (!sclpt[i]) MemoryFail (funcname);
		scll[i] = new double [0];
		if (!scll[i]) MemoryFail (funcname);
		sclu[i] = new double [0];
		if (!sclu[i]) MemoryFail (funcname);
		sclinvl[i] = new double [0];
		if (!sclinvl[i]) MemoryFail (funcname);
		sclinvu[i] = new double [0];
		if (!sclinvu[i]) MemoryFail (funcname);
		diagg[i] = new double [0];
		if (!diagg[i]) MemoryFail (funcname);
		dpiv[i] = new double [0];
		if (!dpiv[i]) MemoryFail (funcname);
	};

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctDiagR::~CFctDiagR()
//========================================================================================
CFctDiagR::~CFctDiagR () { // Destructor
//	cout << " On entry to CFctDiagR destructor " << endl;

	for (int iblk=0;iblk<nblks;iblk++) {
//		if (diasize[iblk] != 0) {
		delete [] sclpt[iblk];
		delete [] scll[iblk];
		delete [] sclu[iblk];
		delete [] sclinvl[iblk];
		delete [] sclinvu[iblk];
		delete [] diagg[iblk];
		delete [] dpiv[iblk];
//		};
	};

	delete [] diasize;
	delete [] ibsfile;
	delete [] blk2file;
	delete [] blk2bs;
	delete [] files;
	delete [] sclpt;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] diagg;
	delete [] dpiv;

	nblks = 0;
	nfiles = 0;
	ifile = 0;
//	cout << " On return from CFctDiagR destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Setup files
// CFctDiagR::SetupFiles()
//========================================================================================
void CFctDiagR::SetupFiles (int _nfiles, char *_name) { // Setup files

	const char *funcname = "SetupFiles";

	delete [] ibsfile;
	delete [] files;
	
	int i;

	nfiles = _nfiles;
	ifile = 0;
	for (i=0;i<nblks*5;i++) {
		blk2file[i] = -1;
	};
	for (i=0;i<nblks*5;i++) {
		blk2bs[i] = -1;
	};
	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);

	OpenDIOFiles (_nfiles, files, "", _name, ".bin");

};

// Author: Kharchenko S.A.
// Description: Close files
// CFctDiagR::CloseFiles()
//========================================================================================
void CFctDiagR::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

	CloseDIOFiles (nfiles, files);

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctC::CFctC()
//========================================================================================
CFctC::CFctC (int _nblks, int _nsups, int _blamx): CFct (_nblks,_nsups) { // Memory allocation zero data constructor

	const char *funcname = "CFctC_01";

	blamx = _blamx;
	lwork = blamx * 10;

	int i;

	sp2blk = new int [_nsups];
	if (!sp2blk) MemoryFail (funcname);

	dpiv = new double [1];
	if (!dpiv) MemoryFail (funcname);
	sclpt = new dcmplx [1];
	if (!sclpt) MemoryFail (funcname);
	scll = new dcmplx [1];
	if (!scll) MemoryFail (funcname);
	sclu = new dcmplx [1];
	if (!sclu) MemoryFail (funcname);
	sclinvl = new dcmplx [1];
	if (!sclinvl) MemoryFail (funcname);
	sclinvu = new dcmplx [1];
	if (!sclinvu) MemoryFail (funcname);
	diagg = new dcmplx [1];
	if (!diagg) MemoryFail (funcname);
	fmaskl = new dcmplx [1];
	if (!fmaskl) MemoryFail (funcname);
	fmasku = new dcmplx [1];
	if (!fmasku) MemoryFail (funcname);

	adrupdl = new dcmplx * [nsups+1];
	if (!adrupdl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) adrupdl[i] = 0;
	adrupdu = new dcmplx * [nsups+1];
	if (!adrupdu) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) adrupdu[i] = 0;
	gl = new dcmplx * [nsups+1];
	if (!gl) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) gl[i] = 0;
	gu = new dcmplx * [nsups+1];
	if (!gu) MemoryFail (funcname);
	for (i=0;i<nsups+1;i++) gu[i] = 0;

	int blamx_2 = blamx*blamx;

	lelem = new dcmplx [blamx_2];
	if (!lelem) MemoryFail (funcname);
	uelem = new dcmplx [blamx_2];
	if (!uelem) MemoryFail (funcname);
	deleml = new dcmplx [blamx_2];
	if (!deleml) MemoryFail (funcname);
	delemu = new dcmplx [blamx_2];
	if (!delemu) MemoryFail (funcname);
	aloc = new dcmplx [blamx_2];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [blamx_2];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [blamx_2];
	if (!vloc) MemoryFail (funcname);
	eig = new double [blamx];
	if (!eig) MemoryFail (funcname);
	work = new dcmplx [lwork];
	if (!work) MemoryFail (funcname);
	rwork = new double [lwork];
	if (!rwork) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctC::~CFctC()
//========================================================================================
CFctC::~CFctC () { // Destructor
//	cout << " On entry to CFctC destructor " << endl;
	blamx  = 0;
	lwork  = 0;
	delete [] sp2blk;
	delete [] dpiv;
	delete [] sclpt;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] diagg;
	delete [] fmaskl;
	delete [] fmasku;
	delete [] adrupdl;
	delete [] adrupdu;
	delete [] gl;
	delete [] gu;
	delete [] lelem;
	delete [] uelem;
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] eig;
	delete [] work;
	delete [] rwork;
//	cout << " On return from CFctC destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctDiagC::CFctDiagC()
//========================================================================================
CFctDiagC::CFctDiagC (int _nblks) { // Memory allocation zero data constructor

	const char *funcname = "CFctDiagC_00";

	nblks  = _nblks;
	nfiles = 0;
	ifile = 0;

	int i;

	ibsfile = new int [nfiles];
	if (!ibsfile) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};
	blk2file = new int [nblks*5];
	if (!blk2file) MemoryFail (funcname);
	for (i=0;i<nblks*5;i++) {
		blk2file[i] = -1;
	};
	blk2bs = new int [nblks*5];
	if (!blk2bs) MemoryFail (funcname);
	for (i=0;i<nblks*5;i++) {
		blk2bs[i] = -1;
	};
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

	diasize = new int [nblks+1];
	if (!diasize) MemoryFail (funcname);
	for (i=0;i<nblks+1;i++) {
		diasize[i] = 0;
	};
	sclpt = new dcmplx * [nblks+1];
	if (!sclpt) MemoryFail (funcname);
	scll = new dcmplx * [nblks+1];
	if (!scll) MemoryFail (funcname);
	sclu = new dcmplx * [nblks+1];
	if (!sclu) MemoryFail (funcname);
	sclinvl = new dcmplx * [nblks+1];
	if (!sclinvl) MemoryFail (funcname);
	sclinvu = new dcmplx * [nblks+1];
	if (!sclinvu) MemoryFail (funcname);
	diagg = new dcmplx * [nblks+1];
	if (!diagg) MemoryFail (funcname);
	dpiv = new double * [nblks+1];
	if (!dpiv) MemoryFail (funcname);

	for (i=0;i<nblks;i++) {
/*
		sclpt[i] = new dcmplx [0];
		if (!sclpt[i]) MemoryFail (funcname);
		scll[i] = new dcmplx [0];
		if (!scll[i]) MemoryFail (funcname);
		sclu[i] = new dcmplx [0];
		if (!sclu[i]) MemoryFail (funcname);
		sclinvl[i] = new dcmplx [0];
		if (!sclinvl[i]) MemoryFail (funcname);
		sclinvu[i] = new dcmplx [0];
		if (!sclinvu[i]) MemoryFail (funcname);
		diagg[i] = new dcmplx [0];
		if (!diagg[i]) MemoryFail (funcname);
		dpiv[i] = new double [0];
		if (!dpiv[i]) MemoryFail (funcname);
*/
		sclpt[i] = 0;
		scll[i] = 0;
		sclu[i] = 0;
		sclinvl[i] = 0;
		sclinvu[i] = 0;
		diagg[i] = 0;
		dpiv[i] = 0;
	};

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctDiagC::~CFctDiagC()
//========================================================================================
CFctDiagC::~CFctDiagC () { // Destructor
//	cout << " On entry to CFctDiagC destructor " << endl;

	for (int iblk=0;iblk<nblks;iblk++) {
		if (diasize[iblk] != 0) {
			delete [] sclpt[iblk];
			delete [] scll[iblk];
			delete [] sclu[iblk];
			delete [] sclinvl[iblk];
			delete [] sclinvu[iblk];
			delete [] diagg[iblk];
			delete [] dpiv[iblk];
		};
	};

	int i;

	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};

	delete [] diasize;
	delete [] ibsfile;
	delete [] blk2file;
	delete [] blk2bs;
	delete [] files;
	delete [] fnames;
	delete [] sclpt;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] diagg;
	delete [] dpiv;

	nblks = 0;
	nfiles = 0;
	ifile = 0;
//	cout << " On return from CFctDiagC destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Setup files
// CFctDiagC::SetupFiles()
//========================================================================================
void CFctDiagC::SetupFiles (int _nfiles, char *_name) { // Setup files

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
	for (i=0;i<nblks*5;i++) {
		blk2file[i] = -1;
	};
	for (i=0;i<nblks*5;i++) {
		blk2bs[i] = -1;
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
// CFctDiagC::CloseFiles()
//========================================================================================
void CFctDiagC::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

	CloseDIOFiles (nfiles, files, fnames);

};

// Author: Kharchenko S.A.
// CFctSchur: Implement abstract block Fct for ILU2 with Schur complement parallelization
//========================================================================================
void CFctSchur::AbstractBlockIlu2 (CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu) { // Implement abstract block Fct for ILU2 with Schur complement parallelization

	const char *funcname = "AbstractBlockIlu2";

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",myid,".dat");
//	ofstream fout (strbuff);

// Prepare all work data

	PrepareWorkData ();

// Open the tree structure

	int nlevtot = _tree.GetNlev ();
//	int nnodestot = _tree.GetNnodes ();
	int *pcpuidend = _tree.GetCpuidend ();
	CNode *pnodes = _tree.GetNode ();

// Open the maximal factorization block sparsity structure

	int *piablksp = _ablksp.GetIa ();
	int *pjablksp = _ablksp.GetJa ();

// Get description of the tree for the current cpu

	int *i_lev_arr;

	i_lev_arr = new int [nlevtot];
	if (!i_lev_arr) MemoryFail (funcname);

	int ilev, ifather, inodeloc;

	inodeloc = pcpuidend[myid];

	while (true) {
		ilev = pnodes[inodeloc].GetNodelv ();
		i_lev_arr[ilev] = inodeloc;
		ifather = pnodes[inodeloc].FatherId ();
		if (ifather == inodeloc) break;
		inodeloc = ifather;
	};

// Allocate the work arrays

	int *iagroup;
	int *grp2cpu;
	int *grptype;
	int *iagrp2cpu;
	int *jagrp2cpu;
	int *imaskcpu;
	int *listcpu;
	int *iaset2grp;
	int *listfct;
	int *listschur;
	int *imaskschur;

	iagroup = new int [_nblks+1];
	if (!iagroup) MemoryFail (funcname);
	grp2cpu = new int [_nblks];
	if (!grp2cpu) MemoryFail (funcname);
	grptype = new int [_nblks];
	if (!grptype) MemoryFail (funcname);
	iagrp2cpu = new int [_nblks+1];
	if (!iagrp2cpu) MemoryFail (funcname);
	jagrp2cpu = new int [_nblks*nproc];
	if (!jagrp2cpu) MemoryFail (funcname);
	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	iaset2grp = new int [_nblks+1];
	if (!iaset2grp) MemoryFail (funcname);
	listfct = new int [_nblks];
	if (!listfct) MemoryFail (funcname);
	listschur = new int [_nblks];
	if (!listschur) MemoryFail (funcname);
	imaskschur = new int [_nblks];
	if (!imaskschur) MemoryFail (funcname);

	int icyclecpu = -1;

	int i;

	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	int icycleschur = -1;

	for (i=0;i<_nblks;i++) imaskschur[i] = icycleschur;

// Init Schur blocks

	int j, iblk, jjblk, k, kkblk;

	icycleschur++;

	int nlistschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			if (imaskschur[iblk] != icycleschur) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskschur[iblk] = icycleschur;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskschur[jjblk] != icycleschur) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskschur[jjblk] = icycleschur;
					};
				} else if (jjblk < iblk) {
					for (k=piablksp[jjblk];k<piablksp[jjblk+1];k++) {
						kkblk = pjablksp[k];
						if (kkblk >= iblk) {
							if (imaskschur[kkblk] != icycleschur) {
								listschur[nlistschur] = kkblk;
								nlistschur++;
								imaskschur[kkblk] = icycleschur;
							};
						};
					};
				};
			};
		};
	};

// Add tree data

	int inode = pcpuidend[myid];

	int indbeg, indend;

	while (true) {
		indbeg = pnodes[inode].GetIndbeg ();
		indend = pnodes[inode].GetIndend ();
		for (i=indbeg;i<=indend;i++) {
			if (imaskschur[i] != icycleschur) {
				listschur[nlistschur] = i;
				nlistschur++;
				imaskschur[i] = icycleschur;
			};
		};
		ifather = pnodes[inode].FatherId ();
		if (ifather == inode) break;
		inode = ifather;
	};

// Sort the list

	qsort (listschur, nlistschur, sizeof(int), compint);

// Prepare the scaling

	PrepareScaling (nlistschur, listschur);

// Init Schur blocks

	InitSchurBlocks (nlistschur, listschur);

// Main cycle over the nodes

	int ibegblk, iendblk, ngroups;
	int iproc, jproc, jbegblk, jendblk, icheck, nz, nlistcpu, nsets, igroup;
	int jgroup, jjbegblk, jjendblk, iset, ibeggrp, iendgrp, myidloc;
	int nchildsloc, ichild0, ichild1, iproc0, iproc1, nlistfct;
	int *pchildsloc;

	for (ilev=nlevtot-1;ilev>=0;ilev--) {

		inode = i_lev_arr[ilev];

		ibegblk = pnodes[inode].GetIndbeg ();
		iendblk = pnodes[inode].GetIndend ();
//		cout << " IbegBlk = " << ibegblk << " IendBlk = " << iendblk << endl;

// Create the whole set of continuous groups of blocks (set of consecutive blocks on the same cpu)

		iblk = ibegblk;
		iproc = _blk2cpu[iblk];
		ngroups = 0;
		iagroup[0] = 0;
		grp2cpu[ngroups] = iproc;
		iagroup[ngroups+1] = iblk+1-ibegblk;
		while (iblk < iendblk) {
			iblk++;
			jproc = _blk2cpu[iblk];
			if (jproc == iproc) {
				iagroup[ngroups+1] = iblk+1-ibegblk;
			} else {
				ngroups++;
				iproc = jproc;
				grp2cpu[ngroups] = iproc;
				iagroup[ngroups+1] = iblk+1-ibegblk;
			};
		};
		ngroups++;

// Mark all groups (not used, Schur or compute)

		for (i=0;i<ngroups;i++) {
			iproc = grp2cpu[i];
			if (iproc == myid) {
				grptype[i] = 1;
			} else {
				jbegblk = iagroup[i]+ibegblk;
				jendblk = iagroup[i+1]-1+ibegblk;
				icheck = -1;
				for (iblk=jbegblk;iblk<=jendblk;iblk++) {
					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk <= iblk) {
							if (_blk2cpu[jjblk] == myid) icheck = 1;
						};
					};
				};
				if (icheck > 0) {
					grptype[i] = 0;
				} else {
					grptype[i] = -1;
				};
			};
		};

// For each working group determine the list of cpu's

		iagrp2cpu[0] = 0;
		nz = 0;

		for (i=0;i<ngroups;i++) {
			icyclecpu++;
			jbegblk = iagroup[i]+ibegblk;
			jendblk = iagroup[i+1]-1+ibegblk;
			iproc = grp2cpu[i];
			listcpu[0] = iproc;
			imaskcpu[iproc] = icyclecpu;
			nlistcpu = 1;
			for (iblk=jbegblk;iblk<=jendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk <= iblk) {
						jproc = _blk2cpu[jjblk];
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[nlistcpu] = jproc;
							nlistcpu++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
			};
			if (nlistcpu > 1) qsort (listcpu+1,nlistcpu-1,sizeof(int),compint);
			for (j=0;j<nlistcpu;j++) jagrp2cpu[nz+j] = listcpu[j];
			nz += nlistcpu;
			iagrp2cpu[i+1] = nz;
		};

// Find the sets of sequential groups that should be processed in parallel

		nsets = 0;

		igroup = 0;
		iaset2grp[0] = 0;
		iaset2grp[nsets+1] = igroup+1;
		jbegblk = iagroup[igroup]+ibegblk;
		jendblk = iagroup[igroup+1]-1+ibegblk;
		while (igroup < ngroups-1) {
			jgroup = igroup+1;
			jjbegblk = iagroup[jgroup]+ibegblk;
			jjendblk = iagroup[jgroup+1]-1+ibegblk;
			icheck = 1;
			for (iblk=jjbegblk;iblk<=jjendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk >= jbegblk && jjblk <= jendblk) icheck = -1;
				};
			};
			if (icheck == 1) {
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			} else {
				nsets++;
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jbegblk = iagroup[igroup]+ibegblk;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			};
		};
		nsets++;

// Main cycle over the sets

		for (iset=0;iset<nsets;iset++) {

			ibeggrp = iaset2grp[iset];
			iendgrp = iaset2grp[iset+1]-1;

// Make binary tree type Schur exchanges and summation for all groups of blocks in the set

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {
				if (grptype[igroup] >= 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

// Plan sends and receives

					myidloc = -1;

					nlistcpu = 0;
					for (j=iagrp2cpu[igroup];j<iagrp2cpu[igroup+1];j++) {
						jproc = jagrp2cpu[j];
						if (jproc == myid) myidloc = nlistcpu;
						listcpu[nlistcpu] = jproc;
						nlistcpu++;
					};

// Do exchanges and summation according to the plan

					inodeloc = pcpuidend[myidloc];

					while (true) {
						nchildsloc = pnodes[inodeloc].GetNchilds ();
						pchildsloc = pnodes[inodeloc].GetChilds ();
						iproc = pnodes[inodeloc].GetNodecpu ();
						if (iproc != myidloc) break;
						if (nchildsloc > 1) {
							ichild0 = pchildsloc[0];
							ichild1 = pchildsloc[1];
							iproc0 = pnodes[ichild0].GetNodecpu ();
							iproc1 = pnodes[ichild1].GetNodecpu ();
							if (iproc0 == myidloc) {
								if (iproc1 < nlistcpu) {
									jproc = listcpu[iproc1];
									ReceiveSchurBlocks (jproc, nlistschur, listschur);
									AddSchurBlocks (nlistschur, listschur);
								};
							} else if (iproc1 == myidloc) {
								if (iproc0 < nlistcpu) {
									jproc = listcpu[iproc0];
									ReceiveSchurBlocks (jproc, nlistschur, listschur);
									AddSchurBlocks (nlistschur, listschur);
								};
							};
						};
						ifather = pnodes[inodeloc].FatherId ();
						if (ifather == inodeloc) break;
						iproc0 = pnodes[ifather].GetNodecpu ();
						if (iproc0 != myidloc) {
							jproc = listcpu[iproc0];
							FilterSchurBlocks (nlistschur, listschur);
							SendSchurBlocks (jproc, nlistschur, listschur);
						};
						inodeloc = ifather;
					};

				};
			};

// Wait for completion of sends

			WaitSends ();

// Perform factorization for own blocks

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {

				if (grptype[igroup] == 1) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistfct = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listfct[nlistfct] = i;
						nlistfct++;
					};

					FctAndFilterBlocks (nlistfct, listfct);

// Update Schur complement for all necessary blocks

					icycleschur++;
					nlistschur = 0;
					for (iblk=jbegblk;iblk<=jendblk;iblk++) {
						for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
							jjblk = pjablksp[j];
							if (jjblk > jendblk) {
								if (imaskschur[jjblk] != icycleschur) {
									listschur[nlistschur] = jjblk;
									nlistschur++;
									imaskschur[jjblk] = icycleschur;
								};
							};
						};
					};

					qsort (listschur, nlistschur, sizeof(int), compint);

					UpdateSchurBlocks (nlistfct, listfct, nlistschur, listschur);

// Free second order matrices for fct blocks

					FreeSecondOrderBlocks (nlistfct, listfct);

				} else if (grptype[igroup] == 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

					FreeSchurBlocks (nlistschur, listschur);

				};

			};

		};

	};

// Free work arrays

	delete [] i_lev_arr;
	delete [] iagroup;
	delete [] grp2cpu;
	delete [] grptype;
	delete [] iagrp2cpu;
	delete [] jagrp2cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] iaset2grp;
	delete [] listfct;
	delete [] listschur;
	delete [] imaskschur;

// Destroy all work data

	DestroyWorkData ();

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctSchurR::CFctSchurR()
//========================================================================================
CFctSchurR::CFctSchurR (): CFctSchur () { // Memory allocation zero data constructor

	const char *funcname = "CFctSchurR_00";

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);

	pgmtrlschr = new CGSMatrixR;
	if (!pgmtrlschr) MemoryFail (funcname);
	pgmtruschr = new CGSMatrixR;
	if (!pgmtruschr) MemoryFail (funcname);
	pgmtrladd = new CGSMatrixR;
	if (!pgmtrladd) MemoryFail (funcname);
	pgmtruadd = new CGSMatrixR;
	if (!pgmtruadd) MemoryFail (funcname);

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctSchurR::CFctSchurR()
//========================================================================================
CFctSchurR::CFctSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
								int _nblks, int *_blks, int *_blk2cpu,
								CTree *_ptree, CSMatrix *_pablkstr, 
								void *_pgmtral, void *_pgmtrau, 
								void *_pgmtrl, void *_pgmtru,
								CSlvParam *_pparam
								): CFctSchur () {

	const char *funcname = "CFctSchurR_01";

	pfout = _pfout;
	nblks = _nblks;
	pblks = _blks;
	pblk2cpu = _blk2cpu;
	ptree = _ptree;
	pablkstr = _pablkstr;
	pgmtral = _pgmtral;
	pgmtrau = _pgmtrau;
	pgmtrl = _pgmtrl;
	pgmtru = _pgmtru;
	pparam = _pparam;

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);

	pgmtrlschr = new CGSMatrixR;
	if (!pgmtrlschr) MemoryFail (funcname);
	pgmtruschr = new CGSMatrixR;
	if (!pgmtruschr) MemoryFail (funcname);
	pgmtrladd = new CGSMatrixR;
	if (!pgmtrladd) MemoryFail (funcname);
	pgmtruadd = new CGSMatrixR;
	if (!pgmtruadd) MemoryFail (funcname);

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctSchurR::~CFctSchurR()
//========================================================================================
CFctSchurR::~CFctSchurR () { // Destructor
//	cout << " On entry to CFctSchurR destructor " << endl;

	delete [] imaskflt;

	delete ((CGSMatrix *)pgmtrlschr);
	delete ((CGSMatrix *)pgmtruschr);
	delete ((CGSMatrix *)pgmtrladd);
	delete ((CGSMatrix *)pgmtruadd);

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

//	cout << " On return from CFctSchurR destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctSymmSchurR::CFctSymmSchurR()
//========================================================================================
CFctSymmSchurR::CFctSymmSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
								int _nblks, int *_blks, int *_blk2cpu,
								CTree *_ptree, CSMatrix *_pablkstr, 
								void *_pgmtral, 
								void *_pgmtrl, 
								CSlvParam *_pparam
								): CFctSchurR (_pfout, _nblks, _blks, _blk2cpu,
								_ptree, _pablkstr, 
								_pgmtral, _pgmtral, 
								_pgmtrl, _pgmtrl,
								_pparam) {

//	const char *funcname = "CFctSymmSchurR_00";

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctSchurC::CFctSchurC()
//========================================================================================
CFctSchurC::CFctSchurC (): CFctSchur () { // Memory allocation zero data constructor

	const char *funcname = "CFctSchurC_00";

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);

	pgmtrlschr = new CGSMatrixCS;
	if (!pgmtrlschr) MemoryFail (funcname);
	pgmtruschr = new CGSMatrixCS;
	if (!pgmtruschr) MemoryFail (funcname);
	pgmtrladd = new CGSMatrixCS;
	if (!pgmtrladd) MemoryFail (funcname);
	pgmtruadd = new CGSMatrixCS;
	if (!pgmtruadd) MemoryFail (funcname);

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFctSchurC::CFctSchurC()
//========================================================================================
CFctSchurC::CFctSchurC (std::ofstream *_pfout, // Memory allocation zero data constructor
								int _nblks, int *_blks, int *_blk2cpu,
								CTree *_ptree, CSMatrix *_pablkstr, 
								void *_pgmtral, void *_pgmtrau, 
								void *_pgmtrl, void *_pgmtru,
								CSlvParam *_pparam
								): CFctSchur () {

	const char *funcname = "CFctSchurC_01";

	pfout = _pfout;
	nblks = _nblks;
	pblks = _blks;
	pblk2cpu = _blk2cpu;
	ptree = _ptree;
	pablkstr = _pablkstr;
	pgmtral = _pgmtral;
	pgmtrau = _pgmtrau;
	pgmtrl = _pgmtrl;
	pgmtru = _pgmtru;
	pparam = _pparam;

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);

	pgmtrlschr = new CGSMatrixCS;
	if (!pgmtrlschr) MemoryFail (funcname);
	pgmtruschr = new CGSMatrixCS;
	if (!pgmtruschr) MemoryFail (funcname);
	pgmtrladd = new CGSMatrixCS;
	if (!pgmtrladd) MemoryFail (funcname);
	pgmtruadd = new CGSMatrixCS;
	if (!pgmtruadd) MemoryFail (funcname);

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destructor
// CFctSchurC::~CFctSchurC()
//========================================================================================
CFctSchurC::~CFctSchurC () { // Destructor
//	cout << " On entry to CFctSchurC destructor " << endl;

	delete [] imaskflt;

	delete ((CGSMatrix *)pgmtrlschr);
	delete ((CGSMatrix *)pgmtruschr);
	delete ((CGSMatrix *)pgmtrladd);
	delete ((CGSMatrix *)pgmtruadd);

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

//	cout << " On return from CFctSchurC destructor " << endl;
};
