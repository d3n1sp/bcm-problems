//------------------------------------------------------------------------------------------------
// File: fct.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"
#include "hyst.h"

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Fct.h: Description of the factorization support data
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Fct
#define __Fct

// Preliminary declarations

class CHyst;
class CSMatrix;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFct
{
protected:
// Data
	int nblks;     // nblks is the total number of blocks
	int nsups;     // nsups is the total number of supernodes
	int icycle;    // icycle is the cycle variable for imask array
	int icyclegl;  // icyclegl is the cycle variable for imaskgl array
	int nlist;     // nlist is the total number of elements in the current list
	int nupd;      // nupd is the total number updates to be performed
	int nzr;       // nzr is the total number of second order elements
	int *ibsdia;   // ibsdia  [nblks+1] is fct working array
	int *ibsimask; // ibsimask[nblks+1] is fct working array
	int *ibsfmask; // ibsfmask[nblks+1] is fct working array
	int *imask;    // imask   [nsups+1] is fct working array
	int *imasklev; // imasklev[nsups+1] is fct working array
	int *lstloc;   // lstloc  [nsups+1] is fct working array
	int *lstloc2;  // lstloc2 [nsups+1] is fct working array
	int *iv;       // iv      [nsups+1] is fct working array
	int *madj;     // madj    [nsups+1] is fct working array
	int *madj2;    // madj2   [nsups+1] is fct working array
	int *ibegm;    // ibegm   [nsups+1] is fct working array
	int *ibegm2;   // ibegm2  [nsups+1] is fct working array
	int *ibsblk;   // ibsblk  [nblks+1] is fct working array
	int *lstupd;   // lstupd  [nsups+1] is fct working array
	int *lstupd2;  // lstupd2 [nsups+1] is fct working array
	int *irw2lst;  // irw2lst [nsups+1] is fct working array
	int *ivgl;     // ivgl    [nsups+1] is fct working array
	int *madjgl;   // madjgl  [nsups+1] is fct working array
	int *madj2gl;  // madj2gl [nsups+1] is fct working array
	int *ibegmgl;  // ibegmgl [nsups+1] is fct working array
	int *ibegm2gl; // ibegm2gl[nsups+1] is fct working array
	int *ibsblkgl; // ibsblkgl[nblks+1] is fct working array
	int *imaskgl;  // ibsblkgl[nblks+1] is fct working array
	int *bsdia;    // bsdia   [nsups+1] is fct working array
	int *adrfmask; // adrfmask[nsups+1] is fct working array
	int *ig;       // ig      [nsups+1] is fct working array
	int **jg;      // jg      [nsups+1] is fct working array
	int **jg2;     // jg2     [nsups+1] is fct working array
	int **jg3;     // jg2     [nsups+1] is fct working array
	int elemlev;   // level of the current update element
	int **addrg;   // addrg   [nsups+1] is fct working array
	CDoubleInt *diarr; // diarr [nsups+1] is fct working array
	CInd2Int *ind2arr; // ind2arr [nsups+1] is fct working array
	CInd2Int *ind2arrgl; // ind2arrgl [nsups+1] is fct working array
	int nmodsv;    // nmodsv is the total number of modified singular values
	double ops;    // ops variable stores the total arithmetic costs during factorization
	CHyst hyst;    // hyst is the hystogram of values
public:
// Functions
	int GetNzr () const {return nzr;}; // Get nzr
	int GetNmodsv () const {return nmodsv;}; // Get nmodsv
	double GetOps () const {return ops;}; // Get ops
// Constructors and destructor
	CFct (); // Memory allocation zero data constructor
	CFct (int nblks, int _nsups); // Memory allocation zero data constructor
	CFct (const CFct &_aa) {throw " Copy constructor for Fct class called";}; // Copy constructor
	virtual ~CFct (); // Destructor
// Operator functions
	CFct &operator= (const CFct &_aa) {throw " Equality operator for Fct class called";}; // Equality operator
// Get/set functions
	CHyst * GetHyst () {return &hyst;};  // Get hyst
// Fct functions
	void AllocateGlobalArrays2Index (int _nblks, int *_blks, int _nlistschur, int *_listschur); // Init arrays that support global update factorization routines
	void FreeGlobalArrays2Index (); // Free arrays that support global update factorization routines
	void InitGlobalArraysByBlockData2Index (CSMatrix &_mtra); // Init global arrays by block data
	void ReInitGlobalArraysByBlockData2Index (CSMatrix &_mtra); // ReInit global arrays by block data
// Input/Output functions
// Friend classes
	friend class CGSMatrixCS;
};

#endif

#ifndef __FctR
#define __FctR

// Preliminary declarations

class CGSMatrix;
class CSMatrixR;
class CSMatrixRS;
class CSlvParam;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctR: public CFct
{
protected:
// Data
	int blamx;       // blamx is the maximal supernode size
	int lwork;       // lwork is the size of the working array
	int iblkalloc;   // iblkalloc is the current allocation block number
	int nzused;      // nzused is the current number of used locations in the allocation block
	int nzalloc;     // nzalloc is the current size of allocated data
	int nzallocmax;  // nzallocmax is the maximal size of allocated data
	int **jgblk;     // jgblk[nsups+1] array contains pointers to the allocated jg blocks
	int **jg2blk;    // jg2blk[nsups+1] array contains pointers to the allocated jg2 blocks
	int **jg3blk;    // jg2blk[nsups+1] array contains pointers to the allocated jg3 blocks
	int **adrgblk;   // adrgblk[nsups+1] array contains pointers to the allocated adrg blocks
	double **glblk;  // glblk[nsups+1] array contains pointers to the allocated gl blocks
	double **gublk;  // gublk[nsups+1] array contains pointers to the allocated gu blocks
	double *dpiv;    // dpiv   [           ] is fct working array
	double *sclpt;   // sclpt  [           ] is fct working array
	double *scll;    // scll   [           ] is fct working array
	double *sclu;    // sclu   [           ] is fct working array
	double *sclinvl; // sclinvl[           ] is fct working array
	double *sclinvu; // sclinvu[           ] is fct working array
	double *diagg;   // diagg  [           ] is fct working array
	double *fmaskl;  // fmaskl [           ] is fct working array
	double *fmasku;  // fmasku [           ] is fct working array
	double **adrupdl;// adrupdl[nsups+1    ] is fct working array
	double **adrupdu;// adrupdu[nsups+1    ] is fct working array
	double **gl;     // gl     [nsups+1    ] is fct working array
	double **gu;     // gu     [nsups+1    ] is fct working array
	double *lelem;   // lelem  [blamx*blamx] is fct working array
	double *uelem;   // uelem  [blamx*blamx] is fct working array
	double *deleml;  // deleml [blamx*blamx] is fct working array
	double *delemu;  // delemu [blamx*blamx] is fct working array
	double *aloc;    // aloc   [blamx*blamx] is fct working array
	double *uloc;    // uloc   [blamx*blamx] is fct working array
	double *vloc;    // vloc   [blamx*blamx] is fct working array
	double *eig;     // eig    [blamx]       is fct working array
	double *work;    // work   [lwork]       is fct working array
	double *rwork;   // rwork  [lwork]       is fct working array
public:
// Functions
// Constructors and destructor
	CFctR (int _nblks, int _nsups, int _blamx); // Memory allocation zero data constructor
	CFctR (const CFctR &_aa) : CFct(_aa) {throw " Copy constructor for FctR class called";}; // Copy constructor
	~CFctR (); // Destructor
// Operator functions
	CFctR &operator= (const CFctR &_aa) {throw " Equality operator for FctR class called";}; // Equality operator
// Up level incomplete factorization functions
	void Ich2InitBlkRow (int _iblkrow, double _dshift, // Init current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau, CSMatrixRS &_mtru);
	void Ilu2InitBlkRow2Index (int _iblkrow, // Init current block row
								double _dshift,
								const CSMatrixR &_mtral, const CSMatrixR &_mtrau, 
								CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ich2InitBlkRow2Index (int _iblkrow, // Init current block row
										double _dshift,
										const CSMatrixR &_mtral,
										CSMatrixR &_mtrl);
	void Ilu2InitBlkRow (int _iblkrow, double _dshift, // Init current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, 
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
	void Ich2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau, const CSMatrixRS &_mtru2,
							CSMatrixRS &_mtru);
	void Ilu2UpdateBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
								const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
								const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2,
								CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ilu2UpdateBlkRowBySetOfBlocks2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows
															const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
															int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr, CSMatrixR *_mtru2arr,
															CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows with the global transpose search
																	const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
																	int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr, CSMatrixR *_mtru2arr,
																	CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ich2UpdateBlkRowBySetOfBlocksGlobal2Index (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row by the set of previous block rows with the global transpose search
																	const CGSMatrix &_gmtra, const CSMatrixR &_mtral, 
																	int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr,
																	CSMatrixR &_mtrl);
	void Ilu2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
	void Ilu2PointUpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2,
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
	void Ich2FctBlkRow (int _iblkrow, const CSlvParam &_param, int _isupbrd, // Perform factorization of the current block row
						const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau, CSMatrixRS &_mtru,
						CSMatrixRS &_mtru2);
	void Ilu2FctBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
								const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau, const CSMatrixR &_mtraini,
								CSMatrixR &_mtrl, CSMatrixR &_mtru,
								CSMatrixR &_mtrl2, CSMatrixR &_mtru2);
	void Ich2FctBlkRow2Index (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
										const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtraini,
										CSMatrixR &_mtrl,
										CSMatrixR &_mtrl2);
	void Ilu2FctBlkRow (int _iblkrow, const CSlvParam &_param, int _isupbrd, // Perform factorization of the current block row
								const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, const CSMatrixRS &_mtraini,
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
								CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2);
	void Ilu2PointFctBlkRow (int _iblkrow, const CSlvParam &_param, int _isupbrd, // Perform factorization of the current block row
										const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, const CSMatrixRS &_mtraini,
										CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
										CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2);
	CSMatrixR AddBlocks2Index ( // Add two block rows/columns assuming that the list of indices coincide
									CSMatrixR &_mtru, CSMatrixR &_mtru2);
	CSMatrixRS AddBlocks (const CGSMatrix &_gmtra, // Add two block rows/columns assuming that the list of indices coincide
							CSMatrixRS &_mtru, CSMatrixRS &_mtru2);
	CSMatrixRS CreateBlock (int _isupbeg, int _isupend, const CGSMatrix &_gmtra); // Init zero data diagonal block
	CSMatrixR CreateBlock2Index (int _iblk, const CGSMatrix &_gmtra); // Init zero data diagonal block
// Low-level incomplete factorization functions
	void Ich2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
						const CGSMatrix &_gmtra, const CSMatrixRS &_mtrau);
	void Ilu2InitRow2Index (char _type, int _iblk, int _irow, // Init current supernode row
								const CGSMatrix &_gmtra, const CSMatrixR &_mtral, const CSMatrixR &_mtrau);
	void Ich2InitRow2Index (char _type, int _iblk, int _irow, // Init current supernode row
											const CGSMatrix &_gmtra, const CSMatrixR &_mtral);
	void Ilu2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
						const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau);
	void Ilu2PointInitRow (char _type, int _ilist, int _irow, // Init current supernode row
							const CGSMatrix &_gmtra, const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau);
	void Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra);
	void Ilu2InitMask2Index (bool _elmisr, int _ilstprv); // Init mask and update arrays
	void Ich2InitMask2Index (bool _elmisr, int _ilstprv); // Init mask and update arrays
	void Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra);
	void Ilu2PointInitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
									const CGSMatrix &_gmtra);
	void Ich2UpdateRow (int _irow, int _irwprv, // Update current row
						const CGSMatrix &_gmtra);
	void Ilu2UpdateRow2Index (); // Update current row
	void Ilu2UpdateRow (int _irow, int _irwprv, // Update current row
						const CGSMatrix &_gmtra);
	void Ilu2PointUpdateRow (int _irow, int _irwprv, // Update current row
							const CGSMatrix &_gmtra);
	void Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend); // Update arrays that support transposed structures search
	void Ilu2UpdateTransp2Index (int _ilist, int _nlist); // Update arrays that support transposed structures search
	void Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						const CGSMatrix &_gmtra);
	void Ilu2FiltrRow2Index (int _iblk, int _irow, // Perform filtering of the current row
								int _fcttyp, double _theta, double _tau2, int _ndiagfct2, int _nlevfct2, 
								const CSMatrix &_mtraini);
	void Ich2FiltrRow2Index (int _iblk, int _irow, // Perform filtering of the current row
								int _fcttyp, double _theta, double _tau2, int _ndiagfct2, int _nlevfct2, 
								const CSMatrix &_mtraini);
	void Ilu2FiltrRow (int _ilist, int _irow, int _fcttyp, // Perform filtering of the current row
								double _theta, double _tau2, int _ndiagfct2, 
								int _isupbrd, double _tau2bord,
								const CGSMatrix &_gmtra, const CSMatrix &_mtraini);
	void Ilu2PointFiltrRow (int _ilist, int _irow, int _fcttyp, // Perform filtering of the current row
									double _theta, double _tau2, int _ndiagfct2,
									int _isupbrd, double _tau2bord,
									const CGSMatrix &_gmtra, const CSMatrix &_mtraini);
	void Ich2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
							const CGSMatrix &_gmtra);
	void Ilu2PivScaleRow2Index (int _iblk, int _irow, double _pivmin); // Compute diagonal pivot and scale current row
	void Ich2PivScaleRow2Index (int _iblk, int _irow, double _pivmin); // Compute diagonal pivot and scale current row
	void Ilu2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
							const CGSMatrix &_gmtra);
	void Ilu2PointPivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
								const CGSMatrix &_gmtra);
	void Ich2StoreRow (int _ilist, int _irow, double _tau1, // Store current row
						int &_nzrtot,
						const CGSMatrix &_gmtra);
	void Ilu2StoreRow2Index (int _iblk, int _irow, // Store current row
								int _fcttyp, double _tau1, int _ndiagfct, int _nlevfct,
								int &_nzrtot);
	void Ich2StoreRow2Index (int _iblk, int _irow, // Store current row
									int _fcttyp, double _tau1, int _ndiagfct, int _nlevfct, 
									int &_nzrtot);
	void Ilu2StoreRow (int _ilist, int _irow, 
								int _fcttyp, double _tau1, int _ndiagfct, 
								int _isupbrd, double _tau1brd, // Store current row
								int &_nzrtot,
								const CGSMatrix &_gmtra);
	void Ilu2PointStoreRow (int _ilist, int _irow, // Store current row
									int _fcttyp, double _tau1, int _ndiagfct, 
									int _isupbrd, double _tau1brd, 
									int &_nzrtot, 
									const CGSMatrix &_gmtra);
	void Ilu2PointStoreRow (int _ilist, int _irow, // Store current row
							const CGSMatrix &_gmtra);
	void Ich2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixRS &_mtru, CSMatrixRS &_mtru2);
	void Ilu2StoreScaleU2Index (int _nlistl, int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
								CSMatrixR &_mtrl, CSMatrixR &_mtru,
								CSMatrixR &_mtrl2, CSMatrixR &_mtru2);
	void Ich2StoreScaleU2Index (int _nlistl, int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
									CSMatrixR &_mtrl,
									CSMatrixR &_mtrl2);
	void Ilu2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
							CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2);
	void Ilu2PointStoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
								const CGSMatrix &_gmtra, 
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru,
								CSMatrixRS &_mtrl2, CSMatrixRS &_mtru2);
	void Ilu2PointStoreU (int _iblkrow, // Store the resulting U objects
								const CGSMatrix &_gmtra, 
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
	void Ilu2InitTransp2Index (int _iblk, int _nlist, // Init arrays that support transposed structures search
							const CSMatrixR &_mtrl2);
	void Ilu2InitTransp2Index (int _iblk, int _nlist, // Init arrays that support transposed structures search
										int _nlistupd, int *_listupd, CSMatrixR *_mtrl2arr);
	void Ilu2InitTransp (int _isupbeg, int _isupend, // Init arrays that support transposed structures search
							const CSMatrixRS &_mtrl2);
	void Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra, const CSMatrixRS &_mtru2);
	void Ilu2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, // Init mask and update arrays
								const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2);
	void Ilu2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, int _irowprv, // Init mask, update arrays and perform computations
										const CSMatrixR &_mtrl2, const CSMatrixR &_mtru2);
	void Ich2InitMask2Index (bool _elmisr, int _ibs, int _ilstprv, int _irowprv, // Init mask, update arrays and perform computations
										const CSMatrixR &_mtrl2);
	void Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra, 
						const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2);
	void Ilu2PointInitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
							const CGSMatrix &_gmtra, 
							const CSMatrixRS &_mtrl2, const CSMatrixRS &_mtru2);
	void Ilu2UpdateTransp2Index (int _iblk, int _ilist, // Update arrays that support transposed structures search
								const CSMatrixR &_mtrl2);
	void Ilu2UpdateTransp2Index (int _iblk, int _ilist, // Update arrays that support transposed structures search
											int *_listupd, CSMatrixR *_mtrl2arr);
	void Ilu2UpdateTranspGlobal2Index (int _iblk, int _irow, // Update arrays that support transposed structures search
													int *_listupd, CSMatrixR *_mtrl2arr);
	void Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend, // Update arrays that support transposed structures search
							const CSMatrixRS &_mtrl2);
	void Ich2StoreRow (int _ilist, int _irow, // Store current row
						const CGSMatrix &_gmtra);
	void Ilu2StoreRow2Index (int _irow, int _nlevshur); // Store current row
	void Ich2StoreRow2Index (int _irow, int _nlevshur); // Store current row
	void Ilu2StoreRow (int _ilist, int _irow, // Store current row
						const CGSMatrix &_gmtra);
	void Ich2StoreU (int _iblkrow, // Store the resulting U objects
						const CGSMatrix &_gmtra, CSMatrixRS &_mtru);
	void Ilu2StoreU2Index (int _nlistl, int _iblkrow, // Store the resulting U objects
								CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ich2StoreU2Index (int _nlistl, int _iblkrow, // Store the resulting U objects
										CSMatrixR &_mtrl);
	void Ilu2StoreU (int _iblkrow, // Store the resulting U objects
						const CGSMatrix &_gmtra, 
						CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
// Input/Output functions
// Friend classes
	friend class CFctDiagR;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
};

#endif

#ifndef __FctC
#define __FctC

// Preliminary declarations

class CGSMatrix;
class CSMatrixCS;
class CSlvParam;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctC: public CFct
{
protected:
// Data
	int blamx;       // blamx is the maximal supernode size
	int lwork;       // lwork is the size of the working array
	int *sp2blk;     // sp2blk [nsups      ] is fct working array
	double *dpiv;    // dpiv   [           ] is fct working array
	dcmplx *sclpt;   // sclpt  [           ] is fct working array
	dcmplx *scll;    // scll   [           ] is fct working array
	dcmplx *sclu;    // sclu   [           ] is fct working array
	dcmplx *sclinvl; // sclinvl[           ] is fct working array
	dcmplx *sclinvu; // sclinvu[           ] is fct working array
	dcmplx *diagg;   // diagg  [           ] is fct working array
	dcmplx *fmaskl;  // fmaskl [           ] is fct working array
	dcmplx *fmasku;  // fmasku [           ] is fct working array
	dcmplx **adrupdl;// adrupdl[nsups+1    ] is fct working array
	dcmplx **adrupdu;// adrupdu[nsups+1    ] is fct working array
	dcmplx **gl;     // gl     [nsups+1    ] is fct working array
	dcmplx **gu;     // gu     [nsups+1    ] is fct working array
	dcmplx *lelem;   // lelem  [blamx*blamx] is fct working array
	dcmplx *uelem;   // uelem  [blamx*blamx] is fct working array
	dcmplx *deleml;  // deleml [blamx*blamx] is fct working array
	dcmplx *delemu;  // delemu [blamx*blamx] is fct working array
	dcmplx *aloc;    // aloc   [blamx*blamx] is fct working array
	dcmplx *uloc;    // uloc   [blamx*blamx] is fct working array
	dcmplx *vloc;    // vloc   [blamx*blamx] is fct working array
	double *eig;     // eig    [blamx]       is fct working array
	dcmplx *work;    // work   [lwork]       is fct working array
	double *rwork;   // rwork  [lwork]       is fct working array
public:
// Functions
// Constructors and destructor
	CFctC (int _nblks, int _nsups, int _blamx); // Memory allocation zero data constructor
	CFctC (const CFctC &_aa) : CFct(_aa)
		{throw " Copy constructor for FctC class called";}; // Copy constructor
	~CFctC (); // Destructor
// Get/set functions
	int * GetSp2blk () {return sp2blk;}; // Get sp2blk
// Operator functions
	CFctC &operator= (const CFctC &_aa) {throw " Equality operator for FctC class called";}; // Equality operator
// Up level incomplete factorization functions
	void Ich2InitBlkRow (int _iblkrow, double _dshift, double _dshiftim, // Init current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, CSMatrixCS &_mtru);
	void Ilu2InitBlkRow (int _iblkrow, double _dshift, double _dshiftim, // Init current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau, 
							CSMatrixCS &_mtrl, CSMatrixCS &_mtru);
	void Ich2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, const CSMatrixCS &_mtru2,
							CSMatrixCS &_mtru);
	void Ilu2UpdateBlkRow (int _iblkrow, const CSlvParam &_param, // Perform updates of the current block row
							const CGSMatrix &_gmtra, const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
							const CSMatrixCS &_mtrl2, const CSMatrixCS &_mtru2,
							CSMatrixCS &_mtrl, CSMatrixCS &_mtru);
	void Ich2FctBlkRow (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
						const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau, CSMatrixCS &_mtru,
						CSMatrixCS &_mtru2);
	void Ilu2FctBlkRow (int _iblkrow, const CSlvParam &_param, // Perform factorization of the current block row
						const CGSMatrix &_gmtra, const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
						CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
						CSMatrixCS &_mtrl2, CSMatrixCS &_mtru2);
	CSMatrixCS AddBlocks (const CGSMatrix &_gmtra, // Add two block rows/columns assuming that the list of indices coincide
							CSMatrixCS &_mtru, CSMatrixCS &_mtru2);
	CSMatrixCS CreateBlock (int _isupbeg, int _isupend, const CGSMatrix &_gmtra); // Init zero data diagonal block
// Low-level incomplete factorization functions
	void Ich2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
						const CGSMatrix &_gmtra, const CSMatrixCS &_mtrau);
	void Ilu2InitRow (char _type, int _ilist, int _irow, // Init current supernode row
						const CGSMatrix &_gmtra, const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau);
	void Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra);
	void Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra);
	void Ich2UpdateRow (int _irow, int _irwprv, // Update current row
						const CGSMatrix &_gmtra);
	void Ilu2UpdateRow (int _irow, int _irwprv, // Update current row
						const CGSMatrix &_gmtra);
	void Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend); // Update arrays that support transposed structures search
	void Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						const CGSMatrix &_gmtra);
	void Ilu2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						const CGSMatrix &_gmtra);
	void Ich2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
							const CGSMatrix &_gmtra);
	void Ilu2PivScaleRow (int _irow, double _pivmin, // Compute diagonal pivot and scale current row
							const CGSMatrix &_gmtra);
	void Ich2StoreRow (int _ilist, int _irow, double _tau1, // Store current row
						int &_nzrtot,
						const CGSMatrix &_gmtra);
	void Ilu2StoreRow (int _ilist, int _irow, double _tau1, // Store current row
						int &_nzrtot,
						const CGSMatrix &_gmtra);
	void Ich2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixCS &_mtru, CSMatrixCS &_mtru2);
	void Ilu2StoreScaleU (int _iblkrow, const CSlvParam &_param, // Store the resulting U objects
							const CGSMatrix &_gmtra, 
							CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
							CSMatrixCS &_mtrl2, CSMatrixCS &_mtru2);
	void Ilu2InitTransp (int _isupbeg, int _isupend, // Init arrays that support transposed structures search
							const CSMatrixCS &_mtrl2);
	void Ich2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra, const CSMatrixCS &_mtru2);
	void Ilu2InitMask (bool _elmisr, int _irow, int _ilstprv, // Init mask and update arrays
						const CGSMatrix &_gmtra, 
						const CSMatrixCS &_mtrl2, const CSMatrixCS &_mtru2);
	void Ilu2UpdateTransp (int _ilist, int _isupbeg, int _isupend, // Update arrays that support transposed structures search
							const CSMatrixCS &_mtrl2);
	void Ich2StoreRow (int _ilist, int _irow, // Store current row
						const CGSMatrix &_gmtra);
	void Ilu2StoreRow (int _ilist, int _irow, // Store current row
						const CGSMatrix &_gmtra);
	void Ich2StoreU (int _iblkrow, // Store the resulting U objects
						const CGSMatrix &_gmtra, CSMatrixCS &_mtru);
	void Ilu2StoreU (int _iblkrow, // Store the resulting U objects
						const CGSMatrix &_gmtra, 
						CSMatrixCS &_mtrl, CSMatrixCS &_mtru);
// Input/Output functions
// Friend classes
	friend class CFctDiagC;
	friend class CGSMatrixCS;
};

#endif

#ifndef __FctSchur
#define __FctSchur

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctSchur
{
public:
// Functions
// Constructors and destructor
	CFctSchur () {}; // Memory allocation zero data constructor
	CFctSchur (const CFctSchur &_aa) {throw " Copy constructor for CFctSchur class called";}; // Copy constructor
	virtual ~CFctSchur () {}; // Destructor
// Operator functions
	CFctSchur &operator= (const CFctSchur &_aa) {throw " Equality operator for CFctSchur class called";}; // Equality operator
// Up level incomplete factorization functions
	void AbstractBlockIlu2 (CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu); // Implement abstract block Fct for ILU2 with Schur complement parallelization
// Low-level incomplete factorization functions
	virtual void PrepareWorkData () = 0; // Prepare fct working data
	virtual void PrepareScaling (int _nlistschur, int *_listschur) = 0; // Prepare the scaling
	virtual void FctAndFilterBlocks (int _nlistfct, int *_listfct) = 0; // Factorize the set of sequential blocks
	virtual void FreeSecondOrderBlocks (int _nlistfct, int *_listfct) = 0; // Free second order data for the set of blocks
	virtual void InitSchurBlocks (int _nlistschur, int *_listschur) = 0; // Init the set of Schur blocks
	virtual void UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur) = 0; // Update Schur data for the set of blocks
	virtual void FilterSchurBlocks (int _nlistschur, int *_listschur) = 0; // Filter the set of Schur blocks
	virtual void SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur) = 0; // Send the set of Schur blocks
	virtual void ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur) = 0; // Receive the set of Schur blocks
	virtual void AddSchurBlocks (int _nlistschur, int *_listschur) = 0; // Add the set of Schur blocks
	virtual void FreeSchurBlocks (int _nlistschur, int *_listschur) = 0; // Free the set of Schur blocks
	virtual void WaitSends () = 0; // Wait for completion of sends
	virtual void DestroyWorkData () = 0; // Destroy fct working data
// Input/Output functions
// Friend classes
};

#endif

#ifndef __FctSchurR
#define __FctSchurR

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctSchurR: public CFctSchur
{
protected:
// Pointers to the external data
	std::ofstream *pfout; // Pointer to the output file
	int nblks; // The number of blocks
	int *pblks; // Pointer to the blocks array
	int *pblk2cpu; // Pointer to the blk2cpu array
	CTree *ptree; // Pointer to the tree
	CSMatrix *pablkstr; // Pointer to the ablkstr
	void *pgmtral; // Pointer to the gmtral
	void *pgmtrau; // Pointer to the gmtrau
	void *pgmtrl; // Pointer to the gmtrl
	void *pgmtru; // Pointer to the gmtru
	CSlvParam *pparam; // Pointer to the params array
// Local data
	int icycleflt; // The block mask index
	int *imaskflt; // The block mask
	void *pgmtrlschr; // Pointer to the second order/Schur L data
	void *pgmtruschr; // Pointer to the second order/Schur U data
	void *pgmtrladd; // Pointer to the received Schur L data
	void *pgmtruadd; // Pointer to the received Schur U data
	int nsends; // The number of sends performed
	int nsendsmax; // The maximal allowable number of sends
	int *imasksend; // The mask of sent data
	char **psenddata; // Pointer to the array of sent data
	CMPIRequest *sendreqvarr; // Array of send requests
	void *pfct; // Pointer to the fct
	void *pfctdiag; // Pointer to the fctdiag
public:
// Functions
// Constructors and destructor
	CFctSchurR (); // Memory allocation zero data constructor
	CFctSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
					int _nblks, int *_blks, int *_blk2cpu, 
					CTree *_ptree, CSMatrix *_pablkstr, 
					void *_pgmtral, void *_pgmtrau, 
					void *_pgmtrl, void *_pgmtru,
					CSlvParam *_pparam);
	CFctSchurR (const CFctSchurR &_aa) {throw " Copy constructor for CFctSchurR class called";}; // Copy constructor
	~CFctSchurR (); // Destructor
// Operator functions
	CFctSchurR &operator= (const CFctSchurR &_aa) {throw " Equality operator for CFctSchurR class called";}; // Equality operator
// Low-level incomplete factorization functions
	void PrepareWorkData (); // Prepare fct working data
	void PrepareScaling (int _nlistschur, int *_listschur); // Prepare the scaling
	void FctAndFilterBlocks (int _nlistfct, int *_listfct); // Factorize the set of sequential blocks
	void FreeSecondOrderBlocks (int _nlistfct, int *_listfct); // Free second order data for the set of blocks
	void InitSchurBlocks (int _nlistschur, int *_listschur); // Init the set of Schur blocks
	void UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur); // Update Schur data for the set of blocks
	void FilterSchurBlocks (int _nlistschur, int *_listschur); // Filter the set of Schur blocks
	void SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Send the set of Schur blocks
	void ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Receive the set of Schur blocks
	void AddSchurBlocks (int _nlistschur, int *_listschur); // Add the set of Schur blocks
	void FreeSchurBlocks (int _nlistschur, int *_listschur); // Free the set of Schur blocks
	void WaitSends (); // Wait for completion of sends
	void DestroyWorkData (); // Destroy fct working data
// Input/Output functions
// Friend classes
};

#endif

#ifndef __FctSymmSchurR
#define __FctSymmSchurR

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctSymmSchurR: public CFctSchurR
{
public:
// Functions
// Constructors and destructor
	CFctSymmSchurR () {}; // Memory allocation zero data constructor
	CFctSymmSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
					int _nblks, int *_blks, int *_blk2cpu, 
					CTree *_ptree, CSMatrix *_pablkstr, 
					void *_pgmtral, 
					void *_pgmtrl, 
					CSlvParam *_pparam);
	CFctSymmSchurR (const CFctSymmSchurR &_aa) {throw " Copy constructor for CFctSymmSchurR class called";}; // Copy constructor
	~CFctSymmSchurR () {}; // Destructor
// Operator functions
	CFctSymmSchurR &operator= (const CFctSymmSchurR &_aa) {throw " Equality operator for CFctSymmSchurR class called";}; // Equality operator
// Low-level incomplete factorization functions
	void FctAndFilterBlocks (int _nlistfct, int *_listfct); // Factorize the set of sequential blocks
	void InitSchurBlocks (int _nlistschur, int *_listschur); // Init the set of Schur blocks
	void UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur); // Update Schur data for the set of blocks
	void SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Send the set of Schur blocks
	void ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Receive the set of Schur blocks
	void AddSchurBlocks (int _nlistschur, int *_listschur); // Add the set of Schur blocks
// Input/Output functions
// Friend classes
};

#endif

#ifndef __FctSchurC
#define __FctSchurC

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFctSchurC: public CFctSchur
{
protected:
// Pointers to the external data
	std::ofstream *pfout; // Pointer to the output file
	int nblks; // The number of blocks
	int *pblks; // Pointer to the blocks array
	int *pblk2cpu; // Pointer to the blk2cpu array
	CTree *ptree; // Pointer to the tree
	CSMatrix *pablkstr; // Pointer to the ablkstr
	void *pgmtral; // Pointer to the gmtral
	void *pgmtrau; // Pointer to the gmtrau
	void *pgmtrl; // Pointer to the gmtrl
	void *pgmtru; // Pointer to the gmtru
	CSlvParam *pparam; // Pointer to the params array
// Local data
	int icycleflt; // The block mask index
	int *imaskflt; // The block mask
	void *pgmtrlschr; // Pointer to the second order/Schur L data
	void *pgmtruschr; // Pointer to the second order/Schur U data
	void *pgmtrladd; // Pointer to the received Schur L data
	void *pgmtruadd; // Pointer to the received Schur U data
	int nsends; // The number of sends performed
	int nsendsmax; // The maximal allowable number of sends
	int *imasksend; // The mask of sent data
	char **psenddata; // Pointer to the array of sent data
	CMPIRequest *sendreqvarr; // Array of send requests
	void *pfct; // Pointer to the fct
	void *pfctdiag; // Pointer to the fctdiag
public:
// Functions
// Constructors and destructor
	CFctSchurC (); // Memory allocation zero data constructor
	CFctSchurC (std::ofstream *_pfout, // Memory allocation zero data constructor
					int _nblks, int *_blks, int *_blk2cpu, 
					CTree *_ptree, CSMatrix *_pablkstr, 
					void *_pgmtral, void *_pgmtrau, 
					void *_pgmtrl, void *_pgmtru,
					CSlvParam *_pparam);
	CFctSchurC (const CFctSchurC &_aa) {throw " Copy constructor for CFctSchurC class called";}; // Copy constructor
	~CFctSchurC (); // Destructor
// Operator functions
	CFctSchurC &operator= (const CFctSchurC &_aa) {throw " Equality operator for CFctSchurC class called";}; // Equality operator
// Low-level incomplete factorization functions
	void PrepareWorkData (); // Prepare fct working data
	void PrepareScaling (int _nlistschur, int *_listschur); // Prepare the scaling
	void FctAndFilterBlocks (int _nlistfct, int *_listfct); // Factorize the set of sequential blocks
	void FreeSecondOrderBlocks (int _nlistfct, int *_listfct); // Free second order data for the set of blocks
	void InitSchurBlocks (int _nlistschur, int *_listschur); // Init the set of Schur blocks
	void UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur); // Update Schur data for the set of blocks
	void FilterSchurBlocks (int _nlistschur, int *_listschur); // Filter the set of Schur blocks
	void SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Send the set of Schur blocks
	void ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur); // Receive the set of Schur blocks
	void AddSchurBlocks (int _nlistschur, int *_listschur); // Add the set of Schur blocks
	void FreeSchurBlocks (int _nlistschur, int *_listschur); // Free the set of Schur blocks
	void WaitSends (); // Wait for completion of sends
	void DestroyWorkData (); // Destroy fct working data
// Input/Output functions
// Friend classes
};

#endif

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif// 
