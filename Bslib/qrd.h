//------------------------------------------------------------------------------------------------
// File: qrd.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Qrd.h: Description of the QR decomposition support data
//
//////////////////////////////////////////////////////////////////////////////

#ifndef __Qrd
#define __Qrd

// Preliminary declarations

class CTree;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CQrd
{
protected:
// Data
	int nlev;          // nlev is the total number of levels
	int nupdmax;       // nupdmax is the maximal number of updates allowed
	int blksiz;        // blksiz is the block size
	int nchildsmax;    // nchildsmax is the maximal number of child nodes
	int nfiles;        // nfiles is the number of files used to store data
	int **mlocarr;     // mlocarr [nlev] is qrd working array
	int **bsqarr;      // bsqarr  [nlev] is qrd working array
	double **qarr;     // qarr    [nlev] is qrd working array
	double **tauarr;   // tauarr  [nlev] is qrd working array
	double *qsendarr;  // qsendarr[    ] is qrd working array
	double *qrecvarr;  // qrecvarr[    ] is qrd working array
	int *bl2file;      // bl2file[nupdmax+1+nfiles] array contains block to file correspondence
	int *bsblk;        // bsblk[nupdmax+1+nfiles] array contains pointers to blocks in files
	FILE **files;      // files[nfiles] array contains description of the files for I/O
	char **fnames;//fnames[nfiles]   array describes the file names used to store the matrix
public:
// Functions
// Constructors and destructor
	CQrd () {throw " Default constructor for Qrd class called";}; // Memory allocation zero data constructor
	CQrd (const CTree &_tree, int _nupdmax, int _blksiz); // Memory allocation zero data constructor
	CQrd (const CQrd &_qrd) {throw " Copy constructor for Qrd class called";}; // Copy constructor
	~CQrd () { // Destructor
//		cout << " On entry to CQrd destructor " << endl;
		nchildsmax = 0;
		nupdmax = 0;
		blksiz = 0;
		for (int ilev=0;ilev<nlev;ilev++) {
			delete [] mlocarr[ilev];
			delete [] bsqarr[ilev];
			delete [] qarr[ilev];
			delete [] tauarr[ilev];
		};
		delete [] mlocarr;
		delete [] bsqarr;
		delete [] tauarr;
		delete [] qarr;
		nlev = 0;
		delete [] qrecvarr;
		delete [] qsendarr;
		delete [] bl2file;
		delete [] bsblk;
		delete [] files;

		int i;
		for (i=0;i<nfiles;i++) {
			delete [] fnames[i];
		};
		delete [] fnames;

		nfiles = 0;
//		cout << " On return from CQrd destructor " << endl;
	};
// Operator functions
	CQrd &operator= (const CQrd &_qrd) {throw " Equality operator for Qrd class called";}; // Equality operator
// Functions that support QR decomposition
	void UpdateQrdBlk (int _ilev, int _iblk, // Update QR decomposition of the block
						int _ldq, double *_qblk);
	void UpdateQrd (const CTree &_tree, // Update QR decomposition of the block in parallel mode
					int _m, int _ibeg, int _iend, double *_a, int _lda, double *_tau);
	void MvmQBlk (int _ilev, int _iblk, int _nrhs, // Multiply Q factor by the current block
					double *_x, int _ldx , double *_qx, int _ldqx);
	void MvmQBlk (const CTree &_tree, // Multiply Q factor by the current block in the parallel mode
					int _m, int _n, int _nrhs,
					double *_q, int _ldq, double *_tau, 
					double *_x, int _ldx, double *_qx, int _ldqx);
	void GetRPartQrd (const CTree &_tree, int _ilev, int _ibeg, int _iend, // Get R part of the QR decomposition
						double *_r, int _ldr);
// Input/Output functions
	void SetupFiles (int _nfiles, char *_name, int _ldq); // Setup files
	void CloseFiles (); // Close files
// Friend classes
};

#endif

#ifndef __QrdC
#define __QrdC

// Preliminary declarations

class CTree;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CQrdC
{
protected:
// Data
	int nlev;          // nlev is the total number of levels
	int nupdmax;       // nupdmax is the maximal number of updates allowed
	int blksiz;        // blksiz is the block size
	int nchildsmax;    // nchildsmax is the maximal number of child nodes
	int nfiles;        // nfiles is the number of files used to store data
	int **mlocarr;     // mlocarr [nlev] is qrd working array
	int **bsqarr;      // bsqarr  [nlev] is qrd working array
	dcmplx **qarr;     // qarr    [nlev] is qrd working array
	dcmplx **tauarr;   // tauarr  [nlev] is qrd working array
	dcmplx *qsendarr;  // qsendarr[    ] is qrd working array
	dcmplx *qrecvarr;  // qrecvarr[    ] is qrd working array
	int *bl2file;      // bl2file[nupdmax+1+nfiles] array contains block to file correspondence
	int *bsblk;        // bsblk[nupdmax+1+nfiles] array contains pointers to blocks in files
	FILE **files;      // files[nfiles] array contains description of the files for I/O
	char **fnames;//fnames[nfiles]   array describes the file names used to store the matrix
public:
// Functions
// Constructors and destructor
	CQrdC () {throw " Default constructor for QrdC class called";}; // Memory allocation zero data constructor
	CQrdC (const CTree &_tree, int _nupdmax, int _blksiz); // Memory allocation zero data constructor
	CQrdC (const CQrdC &_qrd) {throw " Copy constructor for QrdC class called";}; // Copy constructor
	~CQrdC () { // Destructor
//		cout << " On entry to CQrdC destructor " << endl;
		nchildsmax = 0;
		nupdmax = 0;
		blksiz = 0;
		for (int ilev=0;ilev<nlev;ilev++) {
			delete [] mlocarr[ilev];
			delete [] bsqarr[ilev];
			delete [] qarr[ilev];
			delete [] tauarr[ilev];
		};
		delete [] mlocarr;
		delete [] bsqarr;
		delete [] tauarr;
		delete [] qarr;
		nlev = 0;
		delete [] qrecvarr;
		delete [] qsendarr;
		delete [] bl2file;
		delete [] bsblk;
		delete [] files;

		int i;
		for (i=0;i<nfiles;i++) {
			delete [] fnames[i];
		};
		delete [] fnames;

		nfiles = 0;
//		cout << " On return from CQrdC destructor " << endl;
	};
// Operator functions
	CQrdC &operator= (const CQrdC &_qrd) {throw " Equality operator for QrdC class called";}; // Equality operator
// Functions that support QR decomposition
	void UpdateQrdBlk (int _ilev, int _iblk, // Update QR decomposition of the block
						int _ldq, dcmplx *_qblk);
	void UpdateQrd (const CTree &_tree, // Update QR decomposition of the block in parallel mode
					int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau);
	void UpdateQrdRPart (const CTree &_tree, // Update QR decomposition of the block in parallel mode
					int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau);
	void MvmQBlk (int _ilev, int _iblk, int _nrhs, // Multiply Q factor by the current block
					dcmplx *_x, int _ldx , dcmplx *_qx, int _ldqx);
	void MvmQBlk (const CTree &_tree, // Multiply Q factor by the current block in the parallel mode
					int _m, int _n, int _nrhs,
					dcmplx *_q, int _ldq, dcmplx *_tau, 
					dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx);
	void GetRPartQrd (const CTree &_tree, int _ilev, int _ibeg, int _iend, // Get R part of the QR decomposition
						dcmplx *_r, int _ldr);
// Input/Output functions
	void SetupFiles (int _nfiles, char *_name, int _ldq); // Setup files
	void CloseFiles (); // Close files
// Friend classes
};

#endif

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
