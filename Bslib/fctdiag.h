//------------------------------------------------------------------------------------------------
// File: fctdiag.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"
#include "hyst.h"

// FctDiag.h: Description of the diagonal part of the factorization support data
//
//////////////////////////////////////////////////////////////////////////////

#ifndef __FctDiag
#define __FctDiag

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctDiag
{
// Data
public:
// Functions
// Constructors and destructor
	CFctDiag () {}; // Memory allocation zero data constructor
	virtual ~CFctDiag () {}; // Destructor
};

#endif

#ifndef __FctDiagR
#define __FctDiagR

// Preliminary declarations

class CTree;
class CMvm;
class CGSMatrix;
class CGSMatrixR;
class CGSMatrixRS;
class CSMatrixRS;
class CSlvParam;
class CHyst;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctDiagR: CFctDiag
{
// Data
	int nblks;        // nblks  is the total number of blocks
	int nfiles;       // nfiles is the total number of files used to store diagonal data
	int ifile;        // ifile is the current number of the file used to store data
	int *diasize;     // diasize  [nblks]   is fct working array
	int *ibsfile;     // ibsfile  [nfiles]  is fct working array
	int *blk2file;    // blk2file [nblks*5] is fct working array
	int *blk2bs;      // blk2bs   [nblks*5] is fct working array
	FILE **files;     // files    [nfiles]  is fct working array
	double **sclpt;   // sclpt    [nblks]   is fct working array
	double **scll;    // scll     [nblks]   is fct working array
	double **sclu;    // sclu     [nblks]   is fct working array
	double **sclinvl; // sclinvl  [nblks]   is fct working array
	double **sclinvu; // sclinvu  [nblks]   is fct working array
	double **diagg;   // diagg    [nblks]   is fct working array
	double **dpiv;    // dpiv     [nblks]   is fct working array
	CHyst hyst;       // hyst is the hystogram of values
	CHyst hyst2;      // hyst2 is the hystogram of values
public:
// Functions
// Constructors and destructor
	CFctDiagR (int _nblks); // Memory allocation zero data constructor
	CFctDiagR (const CFctDiagR &_aa) {throw " Copy constructor for FctDiagR class called";}; // Copy constructor
	~CFctDiagR (); // Destructor
// Operator functions
	CFctDiagR &operator= (const CFctDiagR &_aa) {throw " Equality operator for FctDiagR class called";}; // Equality operator
// Get/set functions
	CHyst * GetHyst () {return &hyst;}; // Get hyst
	CHyst * GetHyst2 () {return &hyst2;}; // Get hyst2
// Up level incomplete factorization functions
	void StorePointScaling (CGSMatrixR &_mtrl) const; // Store point scaling data in the vector array
	void StorePointScaling (CGSMatrixRS &_mtrl) const; // Store point scaling data in the vector array
	void StoreBlockScaling (char _mtrtype, CGSMatrixR &_mtrlu) const; // Store block scaling data in the vector arrays
	void StoreBlockScaling (char _mtrtype, CGSMatrixRS &_mtrlu) const; // Store block scaling data in the vector arrays
	void FindSups (int _iblkrow, int _nthresh, double *_thresh, // Find supernodes with small singular values
					const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau,
					int *_nlistarr, int **_listarr,
					CFctR &_fct);
	void Ich2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init current scaling for the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau,
							CFctR &_fct);
	void Ilu2ScaleBlkRow2Index (int _iblkrow, int _sclpttype, double _sclmin, // Init scaling for the current block row
								const CGSMatrix &_gmtra, const CSMatrixR &_mtrau,
								CFctR &_fct);
	void Ilu2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init current scaling for the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau,
							CFctR &_fct);
	void ExchangeScaling (bool _is2index, const CTree &_tree, const CMvm &_mvm, const CGSMatrix &_gmtra); // Exchange the scaling
	void ExchangeScalingGeneral (bool _is2index, const CTree &_tree, int *_blk2cpu, int _nlistschur, int *_listschur,
											int *_iab, int *_jab, const CGSMatrix &_gmtra); // Exchange the scaling in the general case
	void Ilu2LoadDiagonalData (int _iblkrow, // Load diagonal data for the computations with the block row
								int _nlstblk, int *_lstblk,
								const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2LoadDiagonalData2Index (int _iblkrow, // Load diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk, int _nlstblkcol, int *_lstblkcol,
										const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2LoadDiagonalData (int _inode, // Load diagonal data for the computations with the set of block rows in parallel mode
								const CTree &_tree, const CMvm &_mvm,
								const int *_iabext, const int *_jabext,
								const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2UnloadDiagonalData (int _iblkrow, // Unload back diagonal data after computations with the block row
									int _nlstblk, int *_lstblk,
									const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2UnloadDiagonalData2Index (int _iblkrow, // Unload diagonal data for the computations with the block row
										int _nlstblk, int *_lstblk,
										const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2UnloadDiagonalData (int _inode, // Unload diagonal data for the computations with the set of block rows in parallel mode
									const CTree &_tree, const CMvm &_mvm,
									const int *_iabext, const int *_jabext,
									const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2LoadIMask (// Load imask data
						int _nlstblk, int *_lstblk, int *_iab, int *_jab,
						const CGSMatrix &_gmtra, CFctR &_fct);
	void Ilu2UnLoadIMask (CFctR &_fct); // Unload imask data
// Low-level incomplete factorization functions
	void SendDiag (const CMPIComm &_comm, int _nlistloc, int* _listloc, int _iproc, int _msgtag, // Send diagg array
					const CGSMatrix &_gmtra);
	void ReceiveDiag (const CMPIComm &_comm, char _optype, int _nlistloc, int* _listloc, int _iproc, int _msgtag); // Receive diagg array
	void AllocateScaling (bool _is2index, int _nlistloc, int *_listloc, // Allocate the memory for scaling according to the prescribed list
									const CGSMatrix &_gmtra);
	void FreeScaling     (int _nlistloc, int *_listloc); // Free the scaling
	void SendScaling     (bool _is2index, const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Send scaling data
									const CGSMatrix &_gmtra);
	void CreateSetOfScalings (bool _is2index, int _nlistloc, int *_listloc, // Create scaling data for set of rows
										char *&_pdata, int &_length, const CGSMatrix &_gmtra);
	void ReceiveScaling  (bool _is2index, const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Receive scaling data
									const CGSMatrix &_gmtra);
	void RestoreSetOfScalings (bool _is2index, int _nlistloc, int *_listloc, // Restore received scaling data
										char *_pdata, int _length,
										const CGSMatrix &_gmtra);
	void WriteScaling    (int _nlistloc, int *_listloc); // Write scaling data to the disk
	void ReadScaling     (int _nlistloc, int *_listloc); // Read  scaling data from the disk
// Input/Output functions
	void SetupFiles (int _nfiles, char *_name); // Setup files
	void CloseFiles (); // Close files
// Friend classes
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
};

#endif

#ifndef __FctDiagC
#define __FctDiagC

// Preliminary declarations

class CTree;
class CMvm;
class CGSMatrix;
class CGSMatrixCS;
class CSMatrixCS;
class CSlvParam;
class CHyst;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctDiagC: CFctDiag
{
// Data
	int nblks;        // nblks  is the total number of blocks
	int nfiles;       // nfiles is the total number of files used to store diagonal data
	int ifile;        // ifile is the current number of the file used to store data
	int *diasize;     // diasize  [nblks]   is fct working array
	int *ibsfile;     // ibsfile  [nfiles]  is fct working array
	int *blk2file;    // blk2file [nblks*5] is fct working array
	int *blk2bs;      // blk2bs   [nblks*5] is fct working array
	FILE **files;     // files    [nfiles]  is fct working array
	char **fnames;    // fnames   [nfiles] array describes the file names used to store the data
	dcmplx **sclpt;   // sclpt    [nblks]   is fct working array
	dcmplx **scll;    // scll     [nblks]   is fct working array
	dcmplx **sclu;    // sclu     [nblks]   is fct working array
	dcmplx **sclinvl; // sclinvl  [nblks]   is fct working array
	dcmplx **sclinvu; // sclinvu  [nblks]   is fct working array
	dcmplx **diagg;   // diagg    [nblks]   is fct working array
	double **dpiv;    // dpiv     [nblks]   is fct working array
	CHyst hyst;       // hyst is the hystogram of values
	CHyst hyst2;      // hyst2 is the hystogram of values
public:
// Functions
// Constructors and destructor
	CFctDiagC (int _nblks); // Memory allocation zero data constructor
	CFctDiagC (const CFctDiagC &_aa) {throw " Copy constructor for FctDiagC class called";}; // Copy constructor
	~CFctDiagC (); // Destructor
// Operator functions
	CFctDiagC &operator= (const CFctDiagC &_aa) {throw " Equality operator for FctDiagC class called";}; // Equality operator
// Get/set functions
	CHyst * GetHyst () {return &hyst;}; // Get hyst
	CHyst * GetHyst2 () {return &hyst2;}; // Get hyst2
// Up level incomplete factorization functions
	void StorePointScaling (CGSMatrixCS &_mtrl) const; // Store point scaling data in the vector array
	void StoreBlockScaling (char _mtrtype, CGSMatrixCS &_mtrlu) const; // Store block scaling data in the vector arrays
	void FindSups (int _iblkrow, int _nthresh, double *_thresh, // Find supernodes with small singular values
					const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau,
					int *_nlistarr, int **_listarr,
					CFctC &_fct);
	void Ich2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init current scaling for the current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau,
							CFctC &_fct);
	void Ilu2ScaleBlkRow (int _iblkrow, int _sclpttype, double _sclmin, // Init current scaling for the current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau,
							CFctC &_fct);
	void ExchangeScaling (const CTree &_tree, const CMvm &_mvm, const CGSMatrix &_gmtra); // Exchange the scaling
	void ExchangeScalingGeneral (const CTree &_tree, int *_blk2cpu, int _nlistschur, int *_listschur, // Exchange the scaling in the general case
											int *_iab, int *_jab, const CGSMatrix &_gmtra);
	void Ilu2LoadDiagonalData (int _iblkrow, // Load diagonal data for the computations with the block row
								int _nlstblk, int *_lstblk,
								const CGSMatrix &_gmtra, CFctC &_fct);
	void Ilu2LoadDiagonalData (int _inode, // Load diagonal data for the computations with the set of block rows in parallel mode
								const CTree &_tree, const CMvm &_mvm,
								const int *_iabext, const int *_jabext,
								const CGSMatrix &_gmtra, CFctC &_fct);
	void Ilu2UnloadDiagonalData (int _iblkrow, // Unload back diagonal data after computations with the block row
									int _nlstblk, int *_lstblk,
									const CGSMatrix &_gmtra, CFctC &_fct);
	void Ilu2UnloadDiagonalData (int _inode, // Unload diagonal data for the computations with the set of block rows in parallel mode
									const CTree &_tree, const CMvm &_mvm,
									const int *_iabext, const int *_jabext,
									const CGSMatrix &_gmtra, CFctC &_fct);
// Low-level incomplete factorization functions
	void SendDiag    (const CMPIComm &_comm, int _nlistloc, int* _listloc, int _iproc, int _msgtag, // Send diagg array
					const CGSMatrix &_gmtra);
	void ReceiveDiag (const CMPIComm &_comm, char _optype, int _nlistloc, int* _listloc, int _iproc, int _msgtag); // Receive diagg array
	void AllocateScaling (int _nlistloc, int *_listloc, // Allocate the memory for scaling according to the prescribed list
									const CGSMatrix &_gmtra);
	void FreeScaling     (int _nlistloc, int *_listloc); // Free the scaling
	void CreateSetOfScalings (int _nlistloc, int *_listloc, // Create scaling data for set of rows
										char *&_pdata, int &_length, const CGSMatrix &_gmtra);
	void SendScaling     (const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Send scaling data
									const CGSMatrix &_gmtra);
	void ReceiveScaling  (const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Receive scaling data
									const CGSMatrix &_gmtra);
	void RestoreSetOfScalings (int _nlistloc, int *_listloc, // Restore received scaling data
										char *_pdata, int _length,
										const CGSMatrix &_gmtra);
	void WriteScaling    (int _nlistloc, int *_listloc); // Write scaling data to the disk
	void ReadScaling     (int _nlistloc, int *_listloc); // Read  scaling data from the disk
// Input/Output functions
	void SetupFiles (int _nfiles, char *_name); // Setup files
	void CloseFiles (); // Close files
// Friend classes
	friend class CGSMatrixCS;
};

#endif
