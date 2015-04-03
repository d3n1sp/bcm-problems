//------------------------------------------------------------------------------------------------
// File: globals.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------

#ifndef __Globals
#define __Globals

#include <cstdio>
#include <fstream>
#include <assert.h>
#include "complex.h"

//#ifndef __FV_2004__
//#define __FV_2004__
//#endif

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

//#define __FlowVisionOut

// Preliminary declarations

// Preliminary declarations

const int CHAR_BUF_SIZE = 2000;
const double DPI = 3.141592653589793238e0;


/////////////////////////////////////////////////////////////////////////////
// Type names declared in this file

typedef dcmplx_my dcmplx;

/////////////////////////////////////////////////////////////////////////////
// Lapack subroutines

extern "C" {

int dlarfg_ (int *n, double *alpha, double *x, 
				int *incx, double *tau);

int dsyev_ (const char *jobz, const char *uplo, int *n, double *a,
			int *lda, double *w, double *work, int *lwork, int *info);

int dgesvd_ (const char *jobu, const char *jobvt, int *m, int *n, 
				double *a, int *lda, double *s, double *u, int *ldu, 
				double *vt, int *ldvt, double *work, int *lwork, int *info);

int zlarfg_ (int *n, dcmplx *alpha, dcmplx *x, 
				int *incx, dcmplx *tau);

int zheev_ (const char *jobz, const char *uplo, int *n, dcmplx *a, int *lda, double *w, dcmplx *work, int *lwork, 
			double *rwork, int *info);

int zgesvd_ (const char *jobu, const char *jobvt, int *m, int *n, 
				dcmplx *a, int *lda, double *s, dcmplx *u, 
				int *ldu, dcmplx *vt, int *ldvt, dcmplx *work, 
				int *lwork, double *rwork, int *info);

};

// Preliminary declaration of the METIS and ordering routines

extern "C" {

void METIS_PartGraphRecursive(int *n, int *xadj, int *adjncy, int *vwgt, 
								int *adjwgt, int *wgtflag, int *numflag, int *nparts, 
								int *options, int *edgecut, int *part);

void METIS_WPartGraphRecursive(int *n, int *xadj, int *adjncy, int *vwgt, 
								int *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts,
								int *options, int *edgecut, int *part);

void METIS_PartGraphVKway(int *nvtxs, int *xadj, int *adjncy, int *vwgt, 
							int *vsize, int *wgtflag, int *numflag, int *nparts, 
							int *options, int *volume, int *part);

void METIS_PartGraphKway(int *nvtxs, int *xadj, int *adjncy, int *vwgt, 
							int *vsize, int *wgtflag, int *numflag, int *nparts, 
							int *options, int *edgecut, int *part);

void METIS_NodeND (int *n, int *xadj, int *adjncy, int *numflag, int *options, int *perm, int *iperm);

void METIS_EdgeND (int *n, int *xadj, int *adjncy, int *numflag, int *options, int *perm, int *iperm);

};

void GPSKCA(int N, int * DEGREE, int * RSTART, int * CONNEC, 
			int OPTPRO, int WRKLEN, int * PERMUT, int * WORK, 
			int & BANDWD, int & PROFIL, int & ERGPS, int & SPACE);

void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask);
void subrcm(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n);

/////////////////////////////////////////////////////////////////////////////
// Functions declared in this file

class CSlvParam;

void CreateLsqObjC (void *&_obj); // Create object of the LsqObjC class
void FreeLsqObjC (void *_obj); // Free data in the LsqObjC object
void DeleteLsqObjC (void *_obj); // Delete object of the LsqObjC class

void LSQDriver (int _slvtype, double _alpha, // Solve complex LSQ problem
						int _m, int _n, dcmplx *_amatr,
						int _nrhs, dcmplx *_rhs, dcmplx *_sol);
void SolveLsqQrd (int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via QRD
						int _nrhs, dcmplx *_rhs, dcmplx *_sol);
void SolveLsqSvd (double _alpha, int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via QRD
						int _nrhs, dcmplx *_rhs, dcmplx *_sol);
void SolveLsqSparse_old (std::ofstream &_fout, CSlvParam &_param,
							int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
							int _nrhs, dcmplx *_rhs, dcmplx *_sol);
void SolveLsqSparse (std::ofstream &_fout, CSlvParam &_param, void *_obj,
							int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
							int _nrhs, dcmplx *_rhs, dcmplx *_sol);
void SolveLsqSparse (int _m, int _n, dcmplx *_amatr, // Solve complex LSQ problem via sparsification
							int _nrhs, dcmplx *_rhs, dcmplx *_sol, int *_iparam, double *_dparam);

void MyZgesvd (int _m, int _n, dcmplx *_amatr, // Compute rectangular SVD
					double *_sv, dcmplx *_u, dcmplx *_v);
void MyDgesvd (int _m, int _n, double *_amatr, // Compute rectangular SVD
					double *_sv, double *_u, double *_v);

void SolveDenseBySVD (int _n, int _nrhs, // Solve system with dense matrix by SVD
								dcmplx *_rmatr, int _ldrmatr,
								dcmplx *_rhs, int _ldrhs,
								dcmplx *_sol, int _ldsol);

void SolveDenseBySVDTikhonov (double _alpha, int _n, int _nrhs, // Solve system with dense matrix by SVD
								dcmplx *_rmatr, int _ldrmatr,
								dcmplx *_rhs, int _ldrhs,
								dcmplx *_sol, int _ldsol);

void SvdBlk (int _m, int _n, // Compute SVD of the current block
				double *_a, double *_sv, double *_u, double *_v, 
				double *_vwork, double *_tau, 
				int _lwork, double *_work);

void SvdBlk (int _m, int _n, // Compute SVD of the current block
				dcmplx *_a, double *_sv, dcmplx *_u, dcmplx *_v, 
				dcmplx *_vwork, dcmplx *_tau, 
				int _lwork, dcmplx *_cwork, double *_work);

void QrdBlk (int _m, int _n, double *_a, double *_tau); // Compute QR decomposition of the block
void QrdBlk (int _m, int _n, double *_a, int _lda, double *_tau); // Compute QR decomposition of the block

void UpdateQrd (int _m, int _ibeg, int _iend, double *_a, int _lda, double *_tau); // Update QR decomposition of the block
void UpdateQrd (int _m, int _iblk, int _blksiz, // Update QR decomposition of the block in the out-of-core mode
				FILE **_afiles, int *_bl2file, int *_bsblk, int _lda, double *_tau,
				double *_abuff);

void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				double *_q, double *_tau, 
				double *_x, double *_qx);
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				double *_q, int _ldq, double *_tau, 
				double *_x, int _ldx, double *_qx, int _ldqx);
void MvmQBlk (int _m, int _nblk, int _blksiz, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, double *_tau, 
				double *_x, int _ldx, double *_qx, int _ldqx,
				double *_qbuff);

void GetRPartQrd (int _m, int _n, // Get R part of the QR decomposition
					double *_q, double *_r);
void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition
					double *_q, int _ldq, double *_r, int _ldr);
void GetRPartQrd (int _m, int _iblkbeg, int _iblkend, int _blksiz, // Get R part of the QR decomposition in the out-of-core mode
					FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, double *_r, int _ldr);

void QrdBlk (int _m, int _n, dcmplx *_a, int _lda, dcmplx *_tau); // Compute QR decomposition of the block

void UpdateQrd (int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau); // Update QR decomposition of the block
void UpdateQrd (int _m, int _ibeg, int _iend, // Update QR decomposition of the block in the out-of-core mode
				FILE *_afile, int _lda, dcmplx *_tau,
				dcmplx *_abuff, int _ncolbuff);
void UpdateQrd (int _m, int _iblk, int _blksiz, // Update QR decomposition of the block in the out-of-core mode
				FILE **_afiles, int *_bl2file, int *_bsblk, int _lda, dcmplx *_tau,
				dcmplx *_abuff);

void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				dcmplx *_q, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx);
void MvmQHBlk (int _m, int _n, int _nrhs, // Multiply QH factor by the current block
					dcmplx *_q, int _ldq, dcmplx *_tau, 
					dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx);
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE *_qfile, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx,
				dcmplx *_qbuff, int _ncolbuff);
void MvmQBlk (int _m, int _nblk, int _blksiz, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx,
				dcmplx *_qbuff);

void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition
					dcmplx *_q, int _ldq, dcmplx *_r, int _ldr);
void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition in the out-of-core mode
					FILE *_qfile, int _ldq, dcmplx *_r, int _ldr);
void GetRPartQrd (int _m, int _iblkbeg, int _iblkend, int _blksiz, // Get R part of the QR decomposition in the out-of-core mode
					FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, dcmplx *_r, int _ldr);

void TransposeMatrix (int _n, int *_ia, int *_ja, int *_iat, int *_jat); // Compute transposed sparsity

void csvhyst (std::ostream &_fout, const char *_name, int _n, double *_sv); // Create the histogram of values
void csvhyst (const char *_name, int _n, double *_sv); // Create the histogram of values

void CollectBool (int _myid, int _nproc, int *_blks, bool *_barr, bool *_barrglob); // Collect bool data

void ReadBool (std::istream &_fin, int _n, bool *_barr); // Read boolean array from the disk

void ReadBoolBin (std::istream &_fin, int _n, bool *_barr); // Read boolean array from disk in binary form

void OrdBool (char _ordtyp, int _n, int *_order, bool *_barr, bool *_barro); // Reorder boolean array

void AddSparseIArray (int _myid, int _nproc, int _nelem, int *_iarr); // Add sparse integer arrays with no common entries

//void AllocateRequests (int _nsends, int _nrecvs, MPI_Request *&_sndreqs, MPI_Request *&_rcvreqs, // Perform requests allocation
//						MPI_Status *&_sndstat, MPI_Status *&_rcvstat, char *_funcname);

//void WaitReceive (int _nrecvs, MPI_Request *_rcvreqs, int &_index, MPI_Status *_rcvstat, // Wait for any receive from the list
//					char *_funcname);

//void WaitSendFreeRequests (int _nsends, int _nrecvs, MPI_Request *_sndreqs, MPI_Request *_rcvreqs, // Wait all sends for completion and free all requests
//							MPI_Status *_sndstat, MPI_Status *_rcvstat, char *_funcname);

void MemoryFail (const char *_funcname); // Memory allocation failed

void IntCompDummy (int &_n);

void StartTimer();
void StopTimer(double& _time0, double& _time1);
void StopTimer(const char* _str);

int compint (const void *_arg1, const void *_arg2); // Compare two integer values

void OutArr (std::ostream &_stream, const char *_name, int _isize, int *_iarr); // Print integer array
void OutArr (std::ostream &_stream, const char *_name, int _isize, float *_farr); // Print float array
void OutArr (std::ostream &_stream, const char *_name, int _isize, double *_darr); // Print double array
void OutArrThresh (std::ostream &_stream, const char *_name, int _isize, double *_darr, double _thresh); // Print double complex array with small value threshold
void OutArr (std::ostream &_stream, const char *_name, int _isize, dcmplx *_carr); // Print double complex array
void OutArrThresh (std::ostream &_stream, const char *_name, int _isize, dcmplx *_carr, double _thresh); // Print double complex array with small value threshold
void OutArr (std::ostream &_stream, const char *_name, int _isize, bool *_barr); // Print boolean array
void OutArr (std::ostream &_stream, const char *_name, int _isize, char *_charr); // Print char array

void OutArrAcc (std::ostream &_stream, const char *_name, int _isize, int *_iarr); // Store integer array
void OutArrAcc (std::ostream &_stream, const char *_name, int _isize, double *_darr); // Store double array

void OutArrAccV (std::ostream &_stream, const char *_name, int _isize, int *_iarr); // Store integer array
void OutArrAccV (std::ostream &_stream, const char *_name, int _isize, double *_darr); // Store double array
void OutArrAccV (std::ostream &_stream, const char *_name, int _isize, bool *_barr); // Store boolean array

void WriteBoolBin (std::ostream &_stream, const char *_name, int _n, bool *_barr); // Read boolean array from disk in binary form

void ReadArrAcc (std::istream &_stream, int _isize, int *_iarr); // Read integer array
void ReadArrAcc (std::istream &_stream, int _isize, double *_darr); // Read double array

void ReadArrAccV (std::istream &_stream, int _isize, int *_iarr); // Read integer array
void ReadArrAccV (std::istream &_stream, int _isize, double *_darr); // Read double array
void ReadArrAccV (std::istream &_stream, int _isize, bool *_barr); // Read boolean array

void FPut (std::fstream &_stream, int _size, double *_darr, int _offset); // Store double array (direct access)
void FGet (std::fstream &_stream, int _size, double *_darr, int _offset); // Read double array (direct access)

void OpenDIOFiles (int _nfiles, std::fstream *_farr, // Open the set of direct IO files according to the path, name and extension
					const char *_path, const char *_name, const char *_extention);

void CloseDIOFiles (int _nfiles, std::fstream *_farr); // Close the set of direct IO files

void FPutG(FILE *_file, int _size, void *_varr, int _offset); // Store void array (direct access)
void FPut (FILE *_file, int _size, int *_iarr, int _offset); // Store integer array (direct access)
void FPut (FILE *_file, int _size, double *_darr, int _offset); // Store double array (direct access)
void FPut (FILE *_file, int _size, dcmplx *_darr, int _offset); // Store double complex array (direct access)

void FGetG(FILE *_file, int _size, void *_varr, int _offset); // Read void array (direct access)
void FGet (FILE *_file, int _size, double *_darr, int _offset); // Read double array (direct access)
void FGet (FILE *_file, int _size, dcmplx *_darr, int _offset); // Read double complex array (direct access)
void FGet (FILE *_file, int _size, char *_charr, int _offset); // Read char array (direct access)

void OpenDIOFiles (int _nfiles, FILE **_farr, // Open the set of direct IO files according to the path, name and extension
					const char *_path, const char *_name, const char *_extention);
void OpenDIOFiles (int _nfiles, FILE **_farr, char **_fnames, // Open the set of direct IO files according to the path, name and extension and store file names
					const char *_path, const char *_name, const char *_extention);

void FlushDIOFiles (int _nfiles, FILE **_farr); // Flush the set of direct IO files

void CloseDIOFiles (int _nfiles, FILE **_farr); // Close the set of direct IO files
void CloseDIOFiles (int _nfiles, FILE **_farr, char **_fnames); // Close the set of direct IO files and remove them

int ReadSymbol (char *_file, int &_k, int _fend, char _symbol); // Read prescribed symbol from file

void WriteSymbol (char _symbol, char *_file, int &_k, int _fend); // Write prescribed symbol into file

char ReadChar (char *_file, int &_k, int _fend); // Read current char from file
char ReadCharAsIs (char *_file, int &_k, int _fend); // Read current char from file
int ReadWordLen (char *_buf, char *_file, int &_k, int _fend, int _m); // Read current word from file according to length

void WriteWordLen (const char *_buf, char *_file, int &_k, int _fend, int _m); // Write current word from file according to length
void WriteIntLen (int _j, char *_file, int &_k, int _fend, int _m); // Write int word into file according to length
void WriteDblLen (double _f, char *_file, int &_k, int _fend, int _m); // Write double word into file according to length

int ReadWord (char *_buf, char *_file, int &_k, int _fend); // Read current word from file
int ReadWordAsIs (char *_buf, char *_file, int &_k, int _fend); // Read current word from file
int ReadWordInBrakets (char *_buf, char *_file, int &_k, int _fend); // Read current word in brakets from file

int LengthOfAsciiFile (char *_filename); // Determine the length of ascii file

char * ReadAsciiFile (char *_filename, int &_length); // Read current file into char variable

char * RemoveComments (char *_filename, int &_length); // Remove comments

void AddLists (int _nlist1, int *_iarr1, // Add two sorted lists
				int _nlist2, int *_iarr2,
				int &_nlist, int *_iarr);

void AddLists (int _nlist1, int *_iarr1, int *_iarrblk1, // Add two sorted 2 index lists
				int _nlist2, int *_iarr2, int *_iarrblk2,
				int &_nlist, int *_iarr, int *_iarrblk);

void SortElm (int _n, int *_iarr, int *_iarr2); // Sort elements of the array
void SortElm (int _n, int *_iarr, double *_darr); // Sort elements of the array
void SortPtr (int _n, int *_iarr, double **_pdarr); // Sort pointers to the double array
void SortPtr (int _n, int *_iarr, dcmplx **_pdarr); // Sort pointers to the double complex array

void Round (double &dx, double &dy); // Round two double values to three digits

std::ostream &SetPw (std::ostream &stream); // Setup output manipulator

void ListBlocks3D (int _nelem, double *_xyzrelem, int _nballs, double *_xyzrballs, // Compute the list of blocks to be checked
					int *&_ilinks, int *&_jlinks);

void InitTracing (double *_times); // Init timing tracing

void FinTracing (std::ofstream &_fout, char *_fname, double *_times); // Finalize timing tracing

class CDoubleInt
{
public:
// Data
	double dvalue; // dvalue contains double value to be sorted
	int intvalue; // intvalue contains integer value to be sorted
// Functions
// Constructors and destructor
	CDoubleInt () {dvalue = 0.0e0; intvalue = 0;}; // Default constructor
	CDoubleInt (double _dvalue, int _intvalue) { // Memory allocation constructor with zero data
		dvalue = _dvalue;
		intvalue = _intvalue;
	};
	~CDoubleInt () {dvalue = 0.0e0; intvalue = 0;}; // Destructor
// Operator functions
	CDoubleInt& operator= (const CDoubleInt &_di2) { // Equality operator
		dvalue = _di2.dvalue;
		intvalue = _di2.intvalue;
		return *this;
	};
// Comparison function
	static int CompareDoubleInt (const void *_arg1, const void *_arg2); // Compare double-integer values
// Input/output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CDoubleInt &_di) { // Output double-int value
		_stream << " CDoubleInt:" << std::endl;
		_stream << " DValue = " << _di.dvalue << " IValue = " << _di.intvalue << std::endl;
		return _stream;
	};
};

//========================================================================
// Compare two double integer values
inline int CDoubleInt::CompareDoubleInt (const void *arg1, const void *arg2) { // Compare two double integer values
	CDoubleInt *diarg1= (CDoubleInt *) arg1;
	CDoubleInt *diarg2= (CDoubleInt *) arg2;
	if (diarg1->dvalue < diarg2->dvalue) {
		return -1;
	} else if (diarg1->dvalue == diarg2->dvalue) {
		return 0;
	} else {
		return 1;
	};
};

class CIntDouble
{
public:
// Data
	int intvalue; // intvalue contains integer value to be sorted
	double dvalue; // dvalue contains double value to be sorted
// Functions
// Constructors and destructor
	CIntDouble () {intvalue = 0; dvalue = 0.0e0;}; // Default constructor
	CIntDouble (int _intvalue, double _dvalue) { // Memory allocation constructor with zero data
		intvalue = _intvalue;
		dvalue = _dvalue;
	};
	~CIntDouble () {intvalue = 0; dvalue = 0.0e0;}; // Destructor
// Operator functions
	CIntDouble &operator= (const CIntDouble &_id2) { // Equality operator
		intvalue = _id2.intvalue;
		dvalue = _id2.dvalue;
		return *this;
	};
// Comparison function
	static int CompareIntDouble (const void *_arg1, const void *_arg2); // Compare integer-double values
};

//========================================================================
// Compare two integer double values
inline int CIntDouble::CompareIntDouble (const void *arg1, const void *arg2) { // Compare two integer double values
	CIntDouble *idarg1= (CIntDouble *) arg1;
	CIntDouble *idarg2= (CIntDouble *) arg2;
	if (idarg1->intvalue < idarg2->intvalue) {
		return -1;
	} else if (idarg1->intvalue == idarg2->intvalue) {
		return 0;
	} else {
		return 1;
	};
};

class CIntInt
{
public:
// Data
	int intvalue; // intvalue contains integer value to be sorted
	int int2value; // int2value contains secondary integer value to be sorted
// Functions
// Constructors and destructor
	CIntInt () {intvalue = 0; int2value = 0;}; // Default constructor
	CIntInt (int _intvalue, int _int2value) { // Memory allocation constructor with zero data
		intvalue = _intvalue;
		int2value = _int2value;
	};
	~CIntInt () {intvalue = 0; int2value = 0;}; // Destructor
// Operator functions
	CIntInt &operator= (const CIntInt &_ii2) { // Equality operator
		intvalue = _ii2.intvalue;
		int2value = _ii2.int2value;
		return *this;
	};
// Comparison functions
	static int CompareIntInt (const void *_arg1, const void *_arg2); // Compare integer-integer values
	inline bool operator< (const CIntInt& _obj) const {
		return intvalue < _obj.intvalue;
	}
};

//========================================================================
// Compare two integer integer values
inline int CIntInt::CompareIntInt (const void *arg1, const void *arg2) { // Compare two integer integer values
	CIntInt *iiarg1= (CIntInt *) arg1;
	CIntInt *iiarg2= (CIntInt *) arg2;
	if (iiarg1->intvalue < iiarg2->intvalue) {
		return -1;
	} else if (iiarg1->intvalue == iiarg2->intvalue) {
		return 0;
	} else {
		return 1;
	};
};

class CInd2Int
{
public:
// Data
	int indx; // indx contains integer x axis value to be sorted
	int indy; // indy contains integer y axis value to be sorted
	int intvalue; // intvalue contains secondary integer value to be sorted
// Functions
// Constructors and destructor
	CInd2Int () {indx=0; indy=0; intvalue=0;}; // Default constructor
	~CInd2Int () {indx=0; indy=0; intvalue=0;}; // Destructor
// Operator functions
	CInd2Int &operator= (const CInd2Int &_ind2i2) { // Equality operator
		indx = _ind2i2.indx;
		indy = _ind2i2.indy;
		intvalue = _ind2i2.intvalue;
		return *this;
	};
// Comparison functions
	static int CompareInd2Int (const void *_arg1, const void *_arg2); // Compare ind3-integer values
	inline bool operator< (const CInd2Int& _obj) const {
		if(indx < _obj.indx)
			return true;
		else if(indx > _obj.indx)
			return false;
		return indy < _obj.indy;
	}
};

//========================================================================
// Compare ind2 integer values
inline int CInd2Int::CompareInd2Int (const void *arg1, const void *arg2) { // Compare ind2 integer values
	CInd2Int *iiarg1= (CInd2Int *) arg1;
	CInd2Int *iiarg2= (CInd2Int *) arg2;
	if (iiarg1->indx < iiarg2->indx) {
		return -1;
	} else if (iiarg1->indx == iiarg2->indx) {
		if (iiarg1->indy < iiarg2->indy) {
			return -1;
		} else {
			return 1;
		};
	} else {
		return 1;
	};
};

class CInd3Int
{
public:
// Data
	int indx; // indx contains integer x axis value to be sorted
	int indy; // indy contains integer y axis value to be sorted
	int indz; // indz contains integer z axis value to be sorted
	int intvalue; // intvalue contains secondary integer value to be sorted
// Functions
// Constructors and destructor
	CInd3Int () {indx=0; indy=0; indz=0; intvalue=0;}; // Default constructor
	~CInd3Int () {indx=0; indy=0; indz=0; intvalue=0;}; // Destructor
// Operator functions
	CInd3Int &operator= (const CInd3Int &_ind3i2) { // Equality operator
		indx = _ind3i2.indx;
		indy = _ind3i2.indy;
		indz = _ind3i2.indz;
		intvalue = _ind3i2.intvalue;
		return *this;
	};
// Comparison function
	static int CompareInd3Int (const void *_arg1, const void *_arg2); // Compare ind3-integer values
	static int CompareZYXInd3Int (const void *_arg1, const void *_arg2); // Compare ind3-integer values
};

//========================================================================
// Compare ind3 integer values
inline int CInd3Int::CompareInd3Int (const void *arg1, const void *arg2) { // Compare ind3 integer values
	CInd3Int *iiarg1= (CInd3Int *) arg1;
	CInd3Int *iiarg2= (CInd3Int *) arg2;
	if (iiarg1->indx < iiarg2->indx) {
		return -1;
	} else if (iiarg1->indx == iiarg2->indx) {
		if (iiarg1->indy < iiarg2->indy) {
			return -1;
		} else if (iiarg1->indy == iiarg2->indy) {
			if (iiarg1->indz < iiarg2->indz) {
				return -1;
			} else if (iiarg1->indz == iiarg2->indz) {
				return 0;
			} else {
				return 1;
			};
		} else {
			return 1;
		};
	} else {
		return 1;
	};
};

//========================================================================
// Compare ind3 integer values
inline int CInd3Int::CompareZYXInd3Int (const void *arg1, const void *arg2) { // Compare ind3 integer values
	CInd3Int *iiarg1= (CInd3Int *) arg1;
	CInd3Int *iiarg2= (CInd3Int *) arg2;
	if (iiarg1->indz < iiarg2->indz) {
		return -1;
	} else if (iiarg1->indz == iiarg2->indz) {
		if (iiarg1->indy < iiarg2->indy) {
			return -1;
		} else if (iiarg1->indy == iiarg2->indy) {
			if (iiarg1->indx < iiarg2->indx) {
				return -1;
			} else if (iiarg1->indx == iiarg2->indx) {
				return 0;
			} else {
				return 1;
			};
		} else {
			return 1;
		};
	} else {
		return 1;
	};
};

template <class _T>
void ReallocateAndCopy (int _nold, int _nnew, _T *&_arr) {
	char *funcname = "ReallocateAndCopy";
	_T *arrnew;
	arrnew = new _T [_nnew];
	if (!arrnew) MemoryFail (funcname);
	int i;
	for (i=0;i<_nold;i++) arrnew[i] = _arr[i];
	delete [] _arr;
	_arr = arrnew;
};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif

#endif
