//------------------------------------------------------------------------------------------------
// File: mvm.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "svector.h"
#include "smatrix.h"
#include "globals.h"
#include "ExchangeMPI.h"

// Mvm.h: Description of the multiplication support data
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Mvm
#define __Mvm

// Preliminary declarations

class CTree;
class CGSMatrix;
class CGSMatrixCS;
class CSVectorC;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvm
{
protected:
// Data
	int nlev;            // nlev    is the total number of levels
	int nsups;           // nsups   is the total number of supernodes
	int nloc;            // nloc    is the local          length of the vector data
	int nlocext;         // nlocext is the local extended length of the vector data
	int nrhs;            // nrhs    is the total number of vectors
	int *bsx;            // bsx [nsups||nblks] is the array that contains base addresses for the necessary supernode data in the local vector
	int *bspx;           // bspx[nsups||nblks] is the array that contains base addresses for the necessary supernode data in the local extended vector
	int *bsx2;           // bsx2 [nblks] is the array that contains secondary base addresses for the necessary supernode data in the local vector
	int *bspx2;          // bspx2[nblks] is the array that contains secondary base addresses for the necessary supernode data in the local extended vector
	int *lev2node;       // lev2node      [nlev]  is mvm working array
	int *blk2cpu;        // blk2cpu       [nblks] is mvm working array
	int *nchildsarr;     // nchildsarr    [nlev]  is mvm working array
	int *nlistupsndarr;  // nlistupsndarr [nlev]  is mvm working array
	int *nlistuprcvarr;  // nlistuprcvarr [nlev]  is mvm working array
	int *nlistdnsndarr;  // nlistdnsndarr [nlev]  is mvm working array
	int *nlistdnrcvarr;  // nlistdnrcvarr [nlev]  is mvm working array
	int *nlistfctarr;    // nlistfctarr   [nlev]  is mvm working array
	int **listupsndarr;  // listupsndarr  [nlev]  is mvm working array
	int **listuprcvarr;  // listuprcvarr  [nlev]  is mvm working array
	int **listdnsndarr;  // listdnsndarr  [nlev]  is mvm working array
	int **listdnrcvarr;  // listdnrcvarr  [nlev]  is mvm working array
	int **list2upsndarr; // list2upsndarr [nlev]  is mvm working array
	int **list2uprcvarr; // list2uprcvarr [nlev]  is mvm working array
	int **list2dnsndarr; // list2dnsndarr [nlev]  is mvm working array
	int **list2dnrcvarr; // list2dnrcvarr [nlev]  is mvm working array
	int **listfctarr;    // listfctarr    [nlev]  is mvm working array
public:
// Functions
// Constructors and destructor
	CMvm () {throw " Default constructor for Mvm class called";}; // Memory allocation zero data constructor
	CMvm (int _nblks, int _nsups, int _nlev); // Constructor
	CMvm (bool _is2index, const CTree &_tree, const CGSMatrix &_gmtra); // Memory allocation zero data constructor
	CMvm (const CMvm &_mvm) {throw " Copy constructor for Mvm class called";}; // Copy constructor
	virtual ~CMvm (); // Destructor
// Operator functions
	CMvm &operator= (const CMvm &_mvm) {throw " Equality operator for Mvm class called";}; // Equality operator
// Input/Output functions
// Friend classes
	friend class CFctDiagR;
	friend class CFctDiagC;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CGSMatrixCS;
};

#endif

#ifndef __MvmR
#define __MvmR

// Preliminary declarations

class CSVectorR;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmR: public CMvm
{
protected:
// Data
	CSVector px;        // px contains main extended vector data
	CSVector pax;       // pax contains additional extended vector data
	double *vectsndrcv;  // vectsndrcv is mvm working array
public:
// Functions
// Constructors and destructor
	CMvmR (const CTree &_tree, const CGSMatrixRS &_gmtra, int _nrhs); // Memory allocation zero data constructor
	CMvmR (const CTree &_tree, const CGSMatrixR &_gmtra, int _nrhs); // Memory allocation zero data constructor
	CMvmR (bool _bsecond, const CTree &_tree, const CGSMatrixRS &_gmtra, int _nrhs); // Memory allocation zero data constructor
	CMvmR (bool _bsecond2ind, const CTree &_tree, const CGSMatrixR &_gmtra, int _nrhs); // Constructor
	CMvmR (const CMvmR &_mvm) : CMvm(_mvm)
		{throw " Copy constructor for MvmR class called";}; // Copy constructor
	~CMvmR (); // Destructor
// Operator functions
	CMvmR &operator= (const CMvmR &_mvm) {throw " Equality operator for MvmR class called";}; // Equality operator
// Send/Receive functions
	void SendVectData    (bool _is2index, const CMPIComm &_comm, char _sendtype, // Send vector data to the prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVector &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr);
	void ReceiveVectData (bool _is2index, const CMPIComm &_comm, char _recvtype, char _optype, // Receive vector data from prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVector &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr);
// Input/Output functions
// Friend classes
	friend class CGSMatrixRS;
	friend class CGSMatrixR;
};

#endif

#ifndef __MvmC
#define __MvmC

// Preliminary declarations

class CSVectorC;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmC: public CMvm
{
protected:
// Data
	CSVectorC px;        // px contains main extended vector data
	CSVectorC pax;       // pax contains additional extended vector data
	dcmplx *vectsndrcv;  // vectsndrcv is mvm working array
public:
// Functions
// Constructors and destructor
	CMvmC (bool _is2index, const CTree &_tree, const CGSMatrixCS &_gmtra, int _nrhs); // Memory allocation zero data constructor
	CMvmC (const CMvmC &_mvm) : CMvm(_mvm)
		{throw " Copy constructor for MvmC class called";}; // Copy constructor
	~CMvmC (); // Destructor
// Operator functions
	CMvmC &operator= (const CMvmC &_mvm) {throw " Equality operator for MvmC class called";}; // Equality operator
// Send/Receive functions
	void SendVectData    (const CMPIComm &_comm, char _sendtype, // Send vector data to the prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVectorC &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr);
	void ReceiveVectData (const CMPIComm &_comm, char _recvtype, char _optype, // Receive vector data from prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVectorC &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr);
// Input/Output functions
// Friend classes
	friend class CGSMatrixCS;
};

#endif

#ifndef __MvmCLUOrd
#define __MvmCLUOrd

// Preliminary declarations

class CSVectorC;
class CSMatrixCS;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmCLUOrd
{
public:
// Data
	CSMatrixCS *pl; // pl contains pointer to L part
	CSMatrixCS *pu; // pu contains pointer to U part
	int *porder;     // porder
public:
// Functions
// Constructors and destructor
	CMvmCLUOrd (CSMatrixCS *_pl, CSMatrixCS *_pu, int *_porder) { // Memory allocation zero data constructor
		pl = _pl;
		pu = _pu;
		porder = _porder;
	};
	CMvmCLUOrd (const CMvmCLUOrd &_mvm) {throw " Copy constructor for MvmCLUOrd class called";}; // Copy constructor
	~CMvmCLUOrd () {}; // Destructor
// Operator functions
	CMvmCLUOrd &operator= (const CMvmCLUOrd &_mvm) {throw " Equality operator for CMvmCLUOrd class called";}; // Equality operator
// Solve function
	static void SolveLUOrd (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CMvmCLUOrd *pobjloc = (CMvmCLUOrd *)_objlu;
		CSMatrixCS *plloc = pobjloc->pl;
		CSMatrixCS *puloc = pobjloc->pu;
		int *porderloc = pobjloc->porder;
		plloc->SolveLUOrd (porderloc, *plloc, *puloc, _x, _px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __MvmSchur
#define __MvmSchur

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmSchur
{
public:
// Functions
// Constructors and destructor
	CMvmSchur () {}; // Memory allocation zero data constructor
	CMvmSchur (const CMvmSchur &_aa) {throw " Copy constructor for CMvmSchur class called";}; // Copy constructor
	virtual ~CMvmSchur () {}; // Destructor
// Operator functions
	CMvmSchur &operator= (const CMvmSchur &_aa) {throw " Equality operator for CMvmSchur class called";}; // Equality operator
// Up level triangular solve and multiplication functions
// Double version
	void AbstractBlockMvmA (char _transp, // Implement abstract block Mvm with Schur complement parallelization
									CSVector &_x, CSVector &_ax);
	void AbstractBlockSolveL (char _lutype, // Implement abstract block L solve with Schur complement parallelization
										CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
										CSVector &_x, CSVector &_px);
	void AbstractBlockSolveU (char _lutype, // Implement abstract block U solve with Schur complement parallelization
										CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
										CSVector &_x, CSVector &_px);
// Double complex version
	void AbstractBlockMvmA (char _transp, // Implement abstract block Mvm with Schur complement parallelization
									CSVectorC &_x, CSVectorC &_ax);
	void AbstractBlockSolveL (char _lutype, // Implement abstract block L solve with Schur complement parallelization
										CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
										CSVectorC &_x, CSVectorC &_px);
	void AbstractBlockSolveU (char _lutype, // Implement abstract block U solve with Schur complement parallelization
										CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
										CSVectorC &_x, CSVectorC &_px);
// Low-level triangular solves functions
// Init/destroy
	virtual void PrepareWorkData () = 0; // Prepare solve working data
	virtual void DestroyWorkData () = 0; // Destroy solve working data
// MvmA
	virtual void ExchangeX (CSVector  &_x) = 0; // Exchange X data
	virtual void ExchangeX (CSVectorC &_x) = 0; // Exchange X data
	virtual void InitAxByZeroes () = 0; // Init necessary Ax by zeroes
	virtual void PerformLocalMultiplications (char _transp) = 0; // Multiply by local block rows/columns of A
	virtual void ExchangeAX (CSVector  &_ax) = 0; // Exchange AX data
	virtual void ExchangeAX (CSVectorC &_ax) = 0; // Exchange AX data
// SlvL/SlvU
	virtual void InitPxBlocks (CSVector  &_x) = 0; // Init the set of Px blocks
	virtual void InitPxBlocks (CSVectorC &_x) = 0; // Init the set of Px blocks
	virtual void SolveLBlocks (char _lutype, int _nlistfct, int *_listfct) = 0; // Solve L with the set of blocks
	virtual void SolveUBlocks (char _lutype, int _nlistfct, int *_listfct) = 0; // Solve U with the set of blocks
	virtual void SendBlocks (int _nlistcpu, int *_listcpu, int _nlistschur, int *_listschur) = 0; // Send the set of X blocks
	virtual void ReceiveBlocks (bool _add, int _jproc, int _nlistschur, int *_listschur) = 0; // Receive the set of Px blocks
	virtual void WaitSends () = 0; // Wait for completion of sends
	virtual void StorePx (CSVector  &_px) = 0; // Store the result
	virtual void StorePx (CSVectorC &_px) = 0; // Store the result
// Input/Output functions
// Friend classes
};

#endif

#ifndef __MvmSchurR
#define __MvmSchurR

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmSchurR: public CMvmSchur
{
protected:
// Prototypes of the called internal functions
	typedef void (FUNC_MULT)  (void *_obj, int _iblk, int *_ibsx, const CSVector &_x, int *_ibspx, CSVector &_px); // Prototype of called mvm functions
	typedef void (FUNC_SOLVE) (void *_obj, int _iblk, int *_ibsx, const CSVector &_x, int *_ibspx, CSVector &_px); // Prototype of called solve functions
protected:
// Pointers to the external data
	std::ofstream *pfout; // Pointer to the output stream
	int nblks; // The number of blocks
	int *pblks; // Pointer to the blocks array
	int *pblk2cpu; // Pointer to the blk2cpu array
	CTree *ptree; // Pointer to the tree
	CSMatrix *pablkstr; // Pointer to the ablkstr
	void *pgmtral; // Pointer to the gmtral
	void *pgmtrau; // Pointer to the gmtrau
	void *pgmtrl; // Pointer to the gmtrl
	void *pgmtru; // Pointer to the gmtru
	FUNC_MULT *pmvmal; // Pointer to the mvmal function
	FUNC_MULT *pmvmau; // Pointer to the mvmau function
	FUNC_SOLVE *psolvel; // Pointer to the solve L function
	FUNC_SOLVE *psolveu; // Pointer to the solve U function
	CSlvParam *pparam; // Pointer to the params array
// Local data
	int icycleflt; // The block mask index
	int *imaskflt; // The block mask
	int *ibsvect; // The array of block base addresses
	int nlistblk; // The number of local blocks on cpu
	int *listblk; // The list of local blocks on cpu
	int nlistblkext; // The extended number of local blocks on cpu
	int *listblkext; // The extended list of local blocks on cpu
	int *iasndblk; // The IA type array that describes the blocks to be sent to each cpu
	int *jasndblk; // The JA type array that describes the blocks to be sent to each cpu
	int *iarcvblk; // The IA type array that describes the blocks to be received from each cpu
	int *jarcvblk; // The JA type array that describes the blocks to be received from each cpu
	int nsends; // The number of sends performed
	int nsendsmax; // The maximal allowable number of sends
	int *imasksend; // The mask of sent data
	char **psenddata; // Pointer to the array of sent data
	CMPIRequest *sendreqvarr; // Array of send requests
	CSVector *pvect; // Pointer to the work vector array
	CSVector *pvect1; // Pointer to the work vector array
public:
// Functions
// Constructors and destructor
	CMvmSchurR (); // Memory allocation zero data constructor
	CMvmSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
					int _nblks, int *_blks, int *_blk2cpu,
					CTree *_ptree, CSMatrix *_pablkstr, 
					void *_pgmtral, void *_pgmtrau, 
					void *_pgmtrl, void *_pgmtru,
					FUNC_MULT *_pmvmal, FUNC_MULT *_pmvmau,
					FUNC_SOLVE *_psolvel, FUNC_MULT *_psolveu,
					CSlvParam *_pparam);
	CMvmSchurR (const CMvmSchurR &_aa) {throw " Copy constructor for CMvmSchurR class called";}; // Copy constructor
	virtual ~CMvmSchurR (); // Destructor
// Operator functions
	CMvmSchurR &operator= (const CMvmSchurR &_aa) {throw " Equality operator for CMvmSchurR class called";}; // Equality operator
// Get/set functions
	int GetNblks () const {return nblks;}; // Get nblks
	int * GetPblks () const {return pblks;}; // Get pblks
	int * GetPblk2cpu () const {return pblk2cpu;}; // Get pblk2cpu
	CSMatrix * GetPablkstr () const {return pablkstr;}; // Get pablkstr
	int * GetIbsvect () const {return ibsvect;}; // Get ibsvect
// Low-level triangular solves functions
// Init/destroy
	void PrepareWorkData (); // Prepare solve working data
	void DestroyWorkData (); // Destroy solve working data
// MvmA
	void ExchangeX (CSVector  &_x); // Exchange X data
	void ExchangeX (CSVectorC &_x) {}; // Exchange X data
	void InitAxByZeroes (); // Init necessary Ax by zeroes
	void PerformLocalMultiplications (char _transp); // Multiply by local block rows/columns of A
	void ExchangeAX (CSVector  &_ax); // Exchange AX data
	void ExchangeAX (CSVectorC &_ax) {}; // Exchange AX data
// SlvL/SlvU
	void InitPxBlocks (CSVector  &_x); // Init the set of Px blocks
	void InitPxBlocks (CSVectorC &_x) {}; // Init the set of Px blocks
	void SolveLBlocks (char _lutype, int _nlistfct, int *_listfct); // Solve L with the set of blocks
	void SolveUBlocks (char _lutype, int _nlistfct, int *_listfct); // Solve U with the set of blocks
	void SendBlocks (int _nlistcpu, int *_listcpu, int _nlistschur, int *_listschur); // Send the set of X blocks
	void ReceiveBlocks (bool _add, int _jproc, int _nlistschur, int *_listschur); // Receive the set of Px blocks
	void WaitSends (); // Wait for completion of sends
	void StorePx (CSVector  &_px); // Store the result
	void StorePx (CSVectorC &_px) {}; // Store the result
// Input/Output functions
// Friend classes
};

#endif

#ifndef __MvmSchurC
#define __MvmSchurC

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CMvmSchurC: public CMvmSchur
{
protected:
// Prototypes of the called internal functions
	typedef void (FUNC_MULTC)  (void *_obj, void *_mtr, int _iblk, int *_ibsx, const CSVectorC &_x, int *_ibspx, CSVectorC &_px); // Prototype of called mvm functions
	typedef void (FUNC_SOLVEC) (void *_obj, void *_mtr, int _iblk, int *_ibsx, const CSVectorC &_x, int *_ibspx, CSVectorC &_px); // Prototype of called solve functions
protected:
// Pointers to the external data
	std::ofstream *pfout; // Pointer to the output stream
	int nblks; // The number of blocks
	int *pblks; // Pointer to the blocks array
	int *pbl2nd; // Pointer to the bl2nd array
	int *pblk2cpu; // Pointer to the blk2cpu array
	CTree *ptree; // Pointer to the tree
	CSMatrix *pablkstr; // Pointer to the ablkstr
	void *pschurmvm; // Pointer to the schurmvm
	void *pgmtral; // Pointer to the gmtral
	void *pgmtrau; // Pointer to the gmtrau
	void *pgmtrl; // Pointer to the gmtrl
	void *pgmtru; // Pointer to the gmtru
	FUNC_MULTC *pmvmal; // Pointer to the mvmal function
	FUNC_MULTC *pmvmau; // Pointer to the mvmau function
	FUNC_SOLVEC *psolvel; // Pointer to the solve L function
	FUNC_SOLVEC *psolveu; // Pointer to the solve U function
	CSlvParam *pparam; // Pointer to the params array
// Local data
	int icycleflt; // The block mask index
	int *imaskflt; // The block mask
	int *ibsvect; // The array of block base addresses
	int nlistblk; // The number of local blocks on cpu
	int *listblk; // The list of local blocks on cpu
	int nlistblkext; // The extended number of local blocks on cpu
	int *listblkext; // The extended list of local blocks on cpu
	int *iasndblk; // The IA type array that describes the blocks to be sent to each cpu
	int *jasndblk; // The JA type array that describes the blocks to be sent to each cpu
	int *iarcvblk; // The IA type array that describes the blocks to be received from each cpu
	int *jarcvblk; // The JA type array that describes the blocks to be received from each cpu
	int nsends; // The number of sends performed
	int nsendsmax; // The maximal allowable number of sends
	int *imasksend; // The mask of sent data
	char **psenddata; // Pointer to the array of sent data
	CMPIRequest *sendreqvarr; // Array of send requests
	CSVectorC *pvect; // Pointer to the work vector array
	CSVectorC *pvect1; // Pointer to the work vector array
public:
// Functions
// Constructors and destructor
	CMvmSchurC (); // Memory allocation zero data constructor
	CMvmSchurC (std::ofstream *_pfout, // Memory allocation zero data constructor
					int _nblks, int *_blks, int *_bl2nd, int *_blk2cpu,
					CTree *_ptree, CSMatrix *_pablkstr, 
					void *_pgmtral, void *_pgmtrau, 
					void *_pgmtrl, void *_pgmtru,
					FUNC_MULTC *_pmvmal, FUNC_MULTC *_pmvmau,
					FUNC_SOLVEC *_psolvel, FUNC_SOLVEC *_psolveu,
					CSlvParam *_pparam);
	CMvmSchurC (const CMvmSchurC &_aa) {throw " Copy constructor for CMvmSchurC class called";}; // Copy constructor
	virtual ~CMvmSchurC (); // Destructor
// Operator functions
	CMvmSchurC &operator= (const CMvmSchurC &_aa) {throw " Equality operator for CMvmSchurC class called";}; // Equality operator
// Get/set functions
	void SetSchurmvm (void *_pschurmvm) {pschurmvm = _pschurmvm;}; // Set schurmvm
	void * GetSchurmvm () const {return pschurmvm;}; // Get schurmvm
	int GetNblks () const {return nblks;}; // Get nblks
	int * GetPblks () const {return pblks;}; // Get pblks
	int * GetPblk2cpu () const {return pblk2cpu;}; // Get pblk2cpu
	CSMatrix * GetPablkstr () const {return pablkstr;}; // Get pablkstr
	int * GetIbsvect () const {return ibsvect;}; // Get ibsvect
// Low-level triangular solves functions
// Init/destroy
	void PrepareWorkData (); // Prepare solve working data
	void DestroyWorkData (); // Destroy solve working data
// MvmA
	void ExchangeX (CSVector  &_x) {}; // Exchange X data
	void ExchangeX (CSVectorC &_x); // Exchange X data
	void InitAxByZeroes (); // Init necessary Ax by zeroes
	void PerformLocalMultiplications (char _transp); // Multiply by local block rows/columns of A
	void ExchangeAX (CSVector  &_ax) {}; // Exchange AX data
	void ExchangeAX (CSVectorC &_ax); // Exchange AX data
// SlvL/SlvU
	void InitPxBlocks (CSVector  &_x) {}; // Init the set of Px blocks
	void InitPxBlocks (CSVectorC &_x); // Init the set of Px blocks
	void SolveLBlocks (char _lutype, int _nlistfct, int *_listfct); // Solve L with the set of blocks
	void SolveUBlocks (char _lutype, int _nlistfct, int *_listfct); // Solve U with the set of blocks
	void SendBlocks (int _nlistcpu, int *_listcpu, int _nlistschur, int *_listschur); // Send the set of X blocks
	void ReceiveBlocks (bool _add, int _jproc, int _nlistschur, int *_listschur); // Receive the set of Px blocks
	void WaitSends (); // Wait for completion of sends
	void StorePx (CSVector  &_px) {}; // Store the result
	void StorePx (CSVectorC &_px); // Store the result
// Input/Output functions
// Friend classes
};

#endif
