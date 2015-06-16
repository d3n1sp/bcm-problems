//------------------------------------------------------------------------------------------------
// File: corr.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include <iostream>

#include "globals.h"

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// corr.h: Description of complex corrector
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __CorrectorC
#define __CorrectorC

// Preliminary declarations

class CCorrectorC
{
// Data
	int ncorr;       // ncorr is the total number of correctors
	int ncorrmax;    // ncorrmax is the maximal allowable number of correctors
	int *ldcorrarr;  // ldcorrarr [ncorrmax] array is the leading dimension of the corrrector blocks (Y and W)
	int *ldcorrarr2; // ldcorrarr2 [ncorrmax] array is the leading dimension of the corrrector blocks (Y and W)
	int *corrsizes;  // corrsizes [ncorrmax] array describes the sizes of correctors
	dcmplx **yarr;   // yarr [ncorrmax] array describes Y part of each corrector
	dcmplx **warr;   // warr [ncorrmax] array describes W part of each corrector
	dcmplx **darr;   // darr [ncorrmax] array describes D part of each corrector
	dcmplx **tauarr; // tauarr [ncorrmax] array describes \tau part of each corrector
	int nfiles;      // nfiles is the total number of files used to store the matrix
	int ifile;       // ifile is the current file number
	int *bl2file;    // bl2file[ncorrmax]  array describes the block to file correspondence
	int *bl2ibsfile; // bl2ibsfile[ncorrmax] array describes the base address of the block in file
	int *ibsfile;    // ibsfile[nfiles] array describes the current free base addresses in file
	FILE **files;    // files[nfiles]   array describes the files used to store the matrix
	char **fnames;   // fnames[nfiles]  array describes the file names used to store the data
public:
// Functions
// Constructors and destructor
	CCorrectorC (); // Memory allocation constructor with zero data
	CCorrectorC (int _ncorrmax); // Memory allocation constructor with zero data
	CCorrectorC (const CCorrectorC &_corr) {throw " Copy constructor for CCorrectorC class called";}; // Copy constructor
	~CCorrectorC (); // Destructor
// Get/set functions
	int GetNcorr () const {return ncorr;}; // Get ncorr
	int * GetLdcorrarr () const {return ldcorrarr;}; // Get ldcorrarr
	int * GetCorrsizes () const {return corrsizes;}; // Get corrsizes
	dcmplx ** GetYarr  () const {return yarr;};      // Get yarr
	dcmplx ** GetWarr  () const {return warr;};      // Get warr
	void GetCorr (int _icorr, int &_ld, int &_ld2, int &_ncol, // Get current corrector
						dcmplx *&_y, dcmplx *&_w, dcmplx *&_d, dcmplx *&_tau) {
		_ld = ldcorrarr[_icorr];
		_ld2 = ldcorrarr2[_icorr];
		_ncol = corrsizes[_icorr];
		_y = yarr[_icorr];
		_w = warr[_icorr];
		_d = darr[_icorr];
		_tau = tauarr[_icorr];
	};
// Operator functions
	CCorrectorC &operator= (const CCorrectorC &_corr); // Equality operator
// Input and output
	void StoreCorr (int _icorr, int _ld, int _ld2, int _ncol, // Store current corrector
							dcmplx *_y, dcmplx *_w, dcmplx *_d, dcmplx *_tau); 
	void SetupFiles (int _nfiles, char *_name); // Setup files
	void CloseFiles (); // Close files
	void FreeFiles (); // Free files
	void WriteCorr    (int _icorr); // Write corrector into the disk
	void ReWriteCorr  (int _icorr); // Rewrite corrector into the disk
	void ReadCorr     (int _icorr); // Read corrector from the disk into the main memory if necessary
	void AllocateCorr (int _icorr); // Allocate corrector
	void FreeCorr     (int _icorr); // Free corrector from main mamory
	friend std::ostream &operator<< (std::ostream &_stream, const CCorrectorC &_corr); // Output correctors
// Friend classes
	friend class CSVectorC;
};

#endif
