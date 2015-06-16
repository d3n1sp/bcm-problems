//------------------------------------------------------------------------------------------------
// File: corr.cpp
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

#include "corr.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A. 
// Description: Memory allocation zero data constructor
// CCorrectorC::CCorrectorC()
//========================================================================================
CCorrectorC::CCorrectorC () { // Memory allocation zero data constructor

	const char *funcname = "CCorrectorC_00";

	ncorr    = 0;
	ncorrmax = 1;
	nfiles = 0;
	ifile = 0;

	int i;

	ldcorrarr = new int [ncorrmax];
	if (!ldcorrarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		ldcorrarr[i] = 0;
	};
	ldcorrarr2 = new int [ncorrmax];
	if (!ldcorrarr2) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		ldcorrarr2[i] = 0;
	};

	corrsizes = new int [ncorrmax];
	if (!corrsizes) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		corrsizes[i] = 0;
	};

	yarr = new dcmplx * [ncorrmax];
	if (!yarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		yarr[i] = new dcmplx [0];
		if (!yarr[i]) MemoryFail (funcname);
	};
	warr = new dcmplx * [ncorrmax];
	if (!warr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		warr[i] = new dcmplx [0];
		if (!warr[i]) MemoryFail (funcname);
	};
	darr = new dcmplx * [ncorrmax];
	if (!darr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		darr[i] = new dcmplx [0];
		if (!darr[i]) MemoryFail (funcname);
	};
	tauarr = new dcmplx * [ncorrmax];
	if (!tauarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		tauarr[i] = new dcmplx [0];
		if (!tauarr[i]) MemoryFail (funcname);
	};

	bl2file = new int [ncorrmax];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [ncorrmax];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
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
// CCorrectorC::CCorrectorC()
//========================================================================================
CCorrectorC::CCorrectorC (int _ncorrmax) { // Memory allocation zero data constructor

	const char *funcname = "CCorrectorC_01";

	ncorr    = 0;
	ncorrmax = _ncorrmax;
	nfiles = 0;
	ifile = 0;

	int i;

	ldcorrarr = new int [ncorrmax];
	if (!ldcorrarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		ldcorrarr[i] = 0;
	};
	ldcorrarr2 = new int [ncorrmax];
	if (!ldcorrarr2) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		ldcorrarr2[i] = 0;
	};

	corrsizes = new int [ncorrmax];
	if (!corrsizes) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		corrsizes[i] = 0;
	};

	yarr = new dcmplx * [ncorrmax];
	if (!yarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		yarr[i] = new dcmplx [0];
		if (!yarr[i]) MemoryFail (funcname);
	};
	warr = new dcmplx * [ncorrmax];
	if (!warr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		warr[i] = new dcmplx [0];
		if (!warr[i]) MemoryFail (funcname);
	};
	darr = new dcmplx * [ncorrmax];
	if (!darr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		darr[i] = new dcmplx [0];
		if (!darr[i]) MemoryFail (funcname);
	};
	tauarr = new dcmplx * [ncorrmax];
	if (!tauarr) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		tauarr[i] = new dcmplx [0];
		if (!tauarr[i]) MemoryFail (funcname);
	};

	bl2file = new int [ncorrmax];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [ncorrmax];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
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
// CCorrectorC::~CCorrectorC()
//========================================================================================
CCorrectorC::~CCorrectorC () { // Destructor
//	std::cout << " CCorrectorC destructor: " << endl;
	ncorr = 0;
	ifile  = 0;
	delete [] ldcorrarr;
	delete [] ldcorrarr2;
	delete [] corrsizes;
	int i;
	for (i=0;i<ncorrmax;i++) {
		delete [] yarr[i];
		delete [] warr[i];
		delete [] darr[i];
		delete [] tauarr[i];
	};
	ncorrmax = 0;
	delete [] yarr;
	delete [] warr;
	delete [] darr;
	delete [] tauarr;
	delete [] bl2file;
	delete [] bl2ibsfile;
	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	nfiles = 0;
	delete [] ibsfile;
	delete [] files;
	delete [] fnames;
//	std::cout << " On return from CCorrectorC destructor " << endl;
};

// Author: Kharchenko S.A. 
// Description: Equality operator
// CCorrectorC &CCorrectorC::operator=()
//========================================================================================
CCorrectorC &CCorrectorC::operator= (const CCorrectorC &_corr) { // Equality operator

	const char *funcname = "CCorrectorC_=";

	if (_corr.nfiles > 0) {
		throw " CCorrectorC::operator=: the case nfiles > 0 is not implemented yet ";
	};

	ncorr = 0;
	ifile  = 0;

	delete [] ldcorrarr;
	delete [] ldcorrarr2;
	delete [] corrsizes;
	int i;
	for (i=0;i<ncorrmax;i++) {
		delete [] yarr[i];
		delete [] warr[i];
		delete [] darr[i];
		delete [] tauarr[i];
	};

	ncorrmax = 0;

	delete [] yarr;
	delete [] warr;
	delete [] darr;
	delete [] tauarr;
	delete [] bl2file;
	delete [] bl2ibsfile;

	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	nfiles = 0;
	delete [] ibsfile;
	delete [] files;
	delete [] fnames;

	ncorr    = _corr.ncorr;
	ncorrmax = _corr.ncorrmax;
	nfiles = _corr.nfiles;
	ifile = _corr.ifile;

	int j;

	ldcorrarr = new int [ncorrmax];
	if (!ldcorrarr) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		ldcorrarr[i] = _corr.ldcorrarr[i];
	};
	for (i=ncorr;i<ncorrmax;i++) {
		ldcorrarr[i] = 0;
	};

	ldcorrarr2 = new int [ncorrmax];
	if (!ldcorrarr2) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		ldcorrarr2[i] = _corr.ldcorrarr2[i];
	};
	for (i=ncorr;i<ncorrmax;i++) {
		ldcorrarr2[i] = 0;
	};

	corrsizes = new int [ncorrmax];
	if (!corrsizes) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		corrsizes[i] = _corr.corrsizes[i];
	};
	for (i=ncorr;i<ncorrmax;i++) {
		corrsizes[i] = 0;
	};

	yarr = new dcmplx * [ncorrmax];
	if (!yarr) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		yarr[i] = new dcmplx [corrsizes[i]*ldcorrarr[i]];
		if (!yarr[i]) MemoryFail (funcname);
		for (j=0;j<corrsizes[i]*ldcorrarr[i];j++) {
			yarr[i][j] = _corr.yarr[i][j];
		};
	};
	for (i=ncorr;i<ncorrmax;i++) {
		yarr[i] = new dcmplx [0];
	};
	warr = new dcmplx * [ncorrmax];
	if (!warr) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		warr[i] = new dcmplx [corrsizes[i]*ldcorrarr2[i]];
		if (!warr[i]) MemoryFail (funcname);
		for (j=0;j<corrsizes[i]*ldcorrarr2[i];j++) {
			warr[i][j] = _corr.warr[i][j];
		};
	};
	for (i=ncorr;i<ncorrmax;i++) {
		warr[i] = new dcmplx [0];
	};
	darr = new dcmplx * [ncorrmax];
	if (!darr) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		darr[i] = new dcmplx [4*corrsizes[i]*corrsizes[i]];
		if (!darr[i]) MemoryFail (funcname);
		for (j=0;j<4*corrsizes[i]*corrsizes[i];j++) {
			darr[i][j] = _corr.darr[i][j];
		};
	};
	for (i=ncorr;i<ncorrmax;i++) {
		darr[i] = new dcmplx [0];
	};
	tauarr = new dcmplx * [ncorrmax];
	if (!tauarr) MemoryFail (funcname);
	for (i=0;i<ncorr;i++) {
		tauarr[i] = new dcmplx [2*corrsizes[i]];
		if (!tauarr[i]) MemoryFail (funcname);
		for (j=0;j<2*corrsizes[i];j++) {
			tauarr[i][j] = _corr.tauarr[i][j];
		};
	};
	for (i=ncorr;i<ncorrmax;i++) {
		tauarr[i] = new dcmplx [0];
	};

	bl2file = new int [ncorrmax];
	if (!bl2file) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
		bl2file[i] = -1;
	};
	bl2ibsfile = new int [ncorrmax];
	if (!bl2ibsfile) MemoryFail (funcname);
	for (i=0;i<ncorrmax;i++) {
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
// Description: Store current corrector
//========================================================================================
void CCorrectorC::StoreCorr (int _icorr, int _ld, int _ld2, int _ncol, // Store current corrector
										dcmplx *_y, dcmplx *_w, dcmplx *_d, dcmplx *_tau) {

	const char *funcname = "StoreCorr";

	if (nfiles > 0) {
		throw " CCorrectorC::StoreCorr: the case nfiles > 0 is not implemented yet ";
	};

	if (_icorr >= ncorrmax) {
		throw " CCorrectorC::StoreCorr: insufficient number of correctors allocated";
	};

	delete [] yarr[_icorr];
	delete [] warr[_icorr];
	delete [] darr[_icorr];
	delete [] tauarr[_icorr];

	ldcorrarr[_icorr] = _ld;
	ldcorrarr2[_icorr] = _ld2;
	corrsizes[_icorr] = _ncol;

	int j;

	yarr[_icorr] = new dcmplx [corrsizes[_icorr]*ldcorrarr[_icorr]];
	if (!yarr[_icorr]) MemoryFail (funcname);
	for (j=0;j<corrsizes[_icorr]*ldcorrarr[_icorr];j++) {
		yarr[_icorr][j] = _y[j];
	};
	warr[_icorr] = new dcmplx [corrsizes[_icorr]*ldcorrarr2[_icorr]];
	if (!warr[_icorr]) MemoryFail (funcname);
	for (j=0;j<corrsizes[_icorr]*ldcorrarr2[_icorr];j++) {
		warr[_icorr][j] = _w[j];
	};
	darr[_icorr] = new dcmplx [4*corrsizes[_icorr]*corrsizes[_icorr]];
	if (!darr[_icorr]) MemoryFail (funcname);
	for (j=0;j<4*corrsizes[_icorr]*corrsizes[_icorr];j++) {
		darr[_icorr][j] = _d[j];
	};
	tauarr[_icorr] = new dcmplx [2*corrsizes[_icorr]];
	if (!tauarr[_icorr]) MemoryFail (funcname);
	for (j=0;j<2*corrsizes[_icorr];j++) {
		tauarr[_icorr][j] = _tau[j];
	};

	if (ncorr < _icorr+1) ncorr = _icorr+1;

};

// Author: Kharchenko S.A. 
// Description: Setup files
// CCorrectorC::SetupFiles()
//========================================================================================
void CCorrectorC::SetupFiles (int _nfiles, char *_name) { // Setup files

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

	for (i=0;i<ncorrmax;i++) {
		bl2file[i] = -1;
	};
	for (i=0;i<ncorrmax;i++) {
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
// CCorrectorC::CloseFiles()
//========================================================================================
void CCorrectorC::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

	CloseDIOFiles (nfiles, files, fnames);

};

// Author: Kharchenko S.A. 
// Description: Free files
// CCorrectorC::FreeFiles()
//========================================================================================
void CCorrectorC::FreeFiles () { // Free files

//	const char *funcname = "FreeFiles";

	ifile = 0;

	int i;

	for (i=0;i<ncorrmax;i++) {
		bl2file[i] = -1;
	};
	for (i=0;i<ncorrmax;i++) {
		bl2ibsfile[i] = -1;
	};
	for (i=0;i<nfiles;i++) {
		ibsfile[i] = 0;
	};

};

// Author: Kharchenko S.A. 
// Description: Write corrector into the disk
// CCorrectorC::WriteCorr()
//========================================================================================
void CCorrectorC::WriteCorr (int _icorr) { // Write corrector into the disk

//	const char *funcname = "WriteCorr";

	int nzaloc, nzaloc2;
	int ibsloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = ldcorrarr[_icorr]*corrsizes[_icorr];
	nzaloc2 = ldcorrarr2[_icorr]*corrsizes[_icorr];

	ibsloc = ibsfile[ifile];

	aloc = yarr[_icorr];

	FPut (files[ifile], nzaloc, aloc, ibsloc);

	bl2ibsfile[_icorr] = ibsloc;
	bl2file[_icorr] = ifile;

	ibsloc += nzaloc;

	aloc = warr[_icorr];

	FPut (files[ifile], nzaloc2, aloc, ibsloc);

	ibsloc += nzaloc;

	ibsfile[ifile] = ibsloc;

	ifile++;
	if (ifile >= nfiles) ifile = 0;

Exit:;

};

// Author: Kharchenko S.A. 
// Description: Rewrite corrector into the disk
// CCorrectorC::ReWriteCorr()
//========================================================================================
void CCorrectorC::ReWriteCorr (int _icorr) { // Rewrite corrector into the disk

//	const char *funcname = "ReWriteCorr";

	int nzaloc, nzaloc2;
	int ibsloc, ifileloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = ldcorrarr[_icorr]*corrsizes[_icorr];
	nzaloc2 = ldcorrarr2[_icorr]*corrsizes[_icorr];
	aloc = yarr[_icorr];

	ibsloc = bl2ibsfile[_icorr];
	ifileloc = bl2file[_icorr];

	FPut (files[ifileloc], nzaloc, aloc, ibsloc);

	ibsloc += nzaloc;

	aloc = warr[_icorr];

	FPut (files[ifileloc], nzaloc2, aloc, ibsloc);

Exit:;

};

// Author: Kharchenko S.A. 
// Description: Read corrector from the disk into the main memory if necessary
// CCorrectorC::ReadCorr()
//========================================================================================
void CCorrectorC::ReadCorr (int _icorr) { // Read corrector from the disk into the main memory if necessary

	const char *funcname = "ReadCorr";

	int nzaloc, nzaloc2;
	int ibsloc;
	int ifileloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	nzaloc = ldcorrarr[_icorr]*corrsizes[_icorr];
	nzaloc2 = ldcorrarr2[_icorr]*corrsizes[_icorr];

	ibsloc = bl2ibsfile[_icorr];
	ifileloc = bl2file[_icorr];

	aloc = yarr[_icorr];

	delete [] aloc;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	FGet (files[ifileloc], nzaloc, aloc, ibsloc);

	ibsloc += nzaloc;

	yarr[_icorr] = aloc;

	aloc = warr[_icorr];

	delete [] aloc;

	aloc = new dcmplx [nzaloc2];
	if (!aloc) MemoryFail (funcname);

	FGet (files[ifileloc], nzaloc2, aloc, ibsloc);

	warr[_icorr] = aloc;

Exit:;

};

// Author: Kharchenko S.A. 
// Description: Allocate corrector
// CCorrectorC::AllocateCorr()
//========================================================================================
void CCorrectorC::AllocateCorr (int _icorr) { // Allocate corrector

	const char *funcname = "AllocateCorr";

	dcmplx czero (0.0e0,0.0e0);

	int nzaloc, nzaloc2;
	dcmplx *aloc;

	nzaloc = ldcorrarr[_icorr]*corrsizes[_icorr];
	nzaloc2 = ldcorrarr2[_icorr]*corrsizes[_icorr];

	aloc = yarr[_icorr];

	delete [] aloc;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	int i;

	for (i=0;i<nzaloc;i++) aloc[i] = czero;

	yarr[_icorr] = aloc;

	aloc = warr[_icorr];

	delete [] aloc;

	aloc = new dcmplx [nzaloc2];
	if (!aloc) MemoryFail (funcname);

	for (i=0;i<nzaloc2;i++) aloc[i] = czero;

	warr[_icorr] = aloc;

};

// Author: Kharchenko S.A. 
// Description: Free corrector from main mamory
// CCorrectorC::FreeCorr()
//========================================================================================
void CCorrectorC::FreeCorr (int _icorr) { // Free corrector from main mamory

	const char *funcname = "FreeCorr";

	int nzaloc;
	dcmplx *aloc;

	if (nfiles == 0) goto Exit;

	aloc = yarr[_icorr];

	delete [] aloc;

	nzaloc = 0;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	yarr[_icorr] = aloc;

	aloc = warr[_icorr];

	delete [] aloc;

	nzaloc = 0;

	aloc = new dcmplx [nzaloc];
	if (!aloc) MemoryFail (funcname);

	warr[_icorr] = aloc;

Exit:;

};

// Author: Kharchenko S.A. 
// Description: Output correctors
// CCorrectorC::<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CCorrectorC &_corr) { // Output correctors

	_stream << " CCorrectorC:" << endl;

	_stream << " NCorr = " << _corr.ncorr << endl;
	_stream << " NCorrMax = " << _corr.ncorrmax << endl;

	_stream << " List of correctors:" << endl;

	for (int i=0;i<_corr.ncorr;i++) {
		_stream << " Icorr = " << i << endl;
		_stream << " Ldsize = " << _corr.ldcorrarr[i] << endl;
		_stream << " Ldsize2 = " << _corr.ldcorrarr2[i] << endl;
		_stream << " Isize = " << _corr.corrsizes[i] << endl;
		OutArr (_stream, " Y =", _corr.corrsizes[i]*_corr.ldcorrarr[i], _corr.yarr[i]);
		OutArr (_stream, " W =", _corr.corrsizes[i]*_corr.ldcorrarr2[i], _corr.warr[i]);
		OutArr (_stream, " D =", 4*_corr.corrsizes[i]*_corr.corrsizes[i], _corr.darr[i]);
		OutArr (_stream, " Tau =", 2*_corr.corrsizes[i], _corr.tauarr[i]);
	};

	return _stream;

};

