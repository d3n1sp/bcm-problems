//------------------------------------------------------------------------------------------------
// File: gsmatrix.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "hyst.h"
#include "smatrix.h"
#include "svector.h"
#include "globals.h"
#include "slvparam.h"
#include "mvm.h"
#include "fct.h"
#include "fctdiag.h"

// GSMatrix.h: Description of the block structure of the matrix
//
//////////////////////////////////////////////////////////////////////////////
// Preliminary declarations

class CSlvParam;
class CSMatrixRS;
class CSMatrixCS;
class CMvmR;
class CMvmC;
class CFctR;
class CFctC;
class CFctDiagR;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

#ifndef __GSMatrix
#define __GSMatrix

class CGSMatrix
{
protected:
// Data
	int m;       // m is the total number of rows    in the matrix
	int n;       // n is the total number of columns in the matrix
	int nblksr;  // nblksr is the number of supernode block rows    in the matrix
	int nblksc;  // nblksc is the number of supernode block columns in the matrix
	int nsupr;   // nsupr is the number of supernode rows    in the matrix
	int nsupc;   // nsupc is the number of supernode columns in the matrix
	int nlist;   // nlist is the local number of block rows/columns in the matrix
	int nlistcnd;// nlistcnd is the local number condensed rows/columns in the matrix
	int nfiles;  // nfiles is the total number of files used to store the matrix
	int ifile;   // ifile is the current file number
	int *blksr;  // blksr[nblksr+1]  array describes the block supernode rows    partitioning of the matrix
	int *blksc;  // blksc[nblksc+1]  array describes the block supernode columns partitioning of the matrix
	int *sprndr; // sprndr[nsupr+1]  array describes the rows    supernodes partitioning of the matrix
	int *sprndc; // sprndc[nsupc+1]  array describes the columns supernodes partitioning of the matrix
	int *bl2ndr; // bl2ndr[nblksr+1] array describes the block point rows    partitioning of the matrix
	int *bl2ndc; // bl2ndc[nblksc+1] array describes the block point columns partitioning of the matrix
	int *listb;  // listb[nlist]     array describes the local block column/row numbers of the matrix
	int *blkscnd;// blkscnd[nblksr+1] array describes the condensed block partitioning of the matrix
	int *listcnd;// listcnd[nlistcnd]  array describes the indices of the condensed matrix for each original block
	int *iblkcopymask; // iblkcopymask[nblksr] array describes the mask of the sparsity copy of the blocks
	int **ptrsparr; // ptrsparr[nblksr*6] array contains if necessary the sparsity pointers
	int *bl2file;// bl2file[nblksc]  array describes the block to file correspondence
	int *bl2ibsfile;// bl2ibsfile[nblksc] array describes the base address of the block in file
	int *ibsfile;// ibsfile[nfiles]  array describes the current free base addresses in file
	FILE **files;// files[nfiles]    array describes the files used to store the matrix
	char **fnames;//fnames[nfiles]   array describes the file names used to store the matrix
public:
// Functions
// Constructors and destructor
	CGSMatrix (); // Memory allocation zero data constructor
	CGSMatrix (int _nblks, int _nsups, int _nlist); // Memory allocation zero data constructor
	CGSMatrix (const CGSMatrix &_aa); // Copy constructor
	virtual ~CGSMatrix (); // Destructor
// Get/set functions
	int GetM       () const {return m;};      // Get m
	int GetN       () const {return n;};      // Get n
	int GetNblksc  () const {return nblksc;}; // Get nblksc
	int GetNblksr  () const {return nblksr;}; // Get nblksr
	int GetNsupc   () const {return nsupc;};  // Get nsupc
	int GetNsupr   () const {return nsupr;};  // Get nsupr
	int GetNlist   () const {return nlist;};  // Get nlist
	int * GetListb () const {return listb;};  // Get listb
	int * GetBlksc () const {return blksc;};  // Get blksc
	int * GetBl2ndc() const {return bl2ndc;};  // Get bl2ndc
	void GetBlks (int *&_blksc, int *&_blksr) const { // Get the blocks partitionings
		_blksc = blksc;
		_blksr = blksr;
	};
	void GetBlsiz (int *&_bl2ndc, int *&_bl2ndr) const { // Get the blocks partitionings in terms of nodes
		_bl2ndc = bl2ndc;
		_bl2ndr = bl2ndr;
	};
	void GetSups (int *&_sprndc, int *&_sprndr) const { // Get the supernodes partitionings
		_sprndc = sprndc;
		_sprndr = sprndr;
	};
	int * GetSprndr () const {return sprndr;}; // Get the supernodes rows partitionings
	void GetSparsity (int *&_listb) const { // Get the sparsity arrays
		_listb = listb;
	};
	void SetM      (int _m)      {m=_m;};             // Set m
	void SetN      (int _n)      {n=_n;};             // Set n
	void SetNblksr (int _nblksr) {nblksr = _nblksr;}; // Set nblksr
	void SetNblksc (int _nblksc) {nblksc = _nblksc;}; // Set nblksc
	void SetNsupr  (int _nsupr)  {nsupr  = _nsupr;};  // Set nsupr
	void SetNsupc  (int _nsupc)  {nsupc  = _nsupc;};  // Set nsupc
	void SetNlist  (int _nlist)  {nlist  = _nlist;};  // Set nlist
	void SetList   (int *_listb) {listb  = _listb;};  // Set listb
	void SetSprndc (int *_sprndc){sprndc = _sprndc;}; // Set sprndc
	void SetSprndr (int *_sprndr){sprndr = _sprndr;}; // Set sprndr
// Operator functions
	CGSMatrix &operator= (const CGSMatrix &_aa); // Equality operator
// Transformation functions
	void BlockSparsity2Index (const CSMatrix &_mtra, // Compute block sparsity of the block row
								int &_nlstblk, int *_lstblk,
								int &_icycle, int *_imask) const;
	void BlockSparsity (const CSMatrix &_mtra, // Compute block sparsity of the block row
						int &_nlstblk, int *_lstblk,
						int &_icycle, int *_imask, 
						const int *_sp2blk) const;
	void MtrBlockSparsity2Index (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
							int *&_iab, int *&_jab) const;
	void MtrBlockSparsity (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
							int *&_iab, int *&_jab) const;
	void MtrBlockSparsity2Index (CSMatrix **_mtrarr, // Compute block sparsity of the matrix
										int _iblkbeg, int _iblkend, 
										int *&_iab, int *&_jab) const;
	void MtrBlockSparsity (CSMatrix **_mtrarr, int _iblkbeg, int _iblkend, // Compute block sparsity of the matrix
									int *&_iab, int *&_jab) const;
	void MtrBlockSparsity2Index (bool _check, const CTree &_tree, int *_bl2cpu, // Compute block sparsity of the matrix in the parallel mode
									CSMatrix **_mtrarr, int *&_iab, int *&_jab) const;
	void MtrBlockSparsity (const CTree &_tree, int *_bl2cpu, // Compute block sparsity of the matrix in the parallel mode
							CSMatrix **_mtrarr,
							int *&_iab, int *&_jab) const;
	void ExtendBlockSparsity (int _nblks, // Compute extended block sparsity with respect to triangular factorization
								int *_iab, int *_jab, int *&_iabext, int *&_jabext) const;
	void ExtendBlockSparsityTri (int _nblks, // Compute extended block sparsity with respect to triangular factorization
									int *_iab, int *_jab, int *&_iabext, int *&_jabext) const;
	void ExtendBlockSparsityOpt (int _nblks, // Compute extended block sparsity with respect to triangular factorization
										int *_iab, int *_jab, int *&_iabext, int *&_jabext) const;
	static void CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												const CSMatrix *_mtrarr,
												CSMatrix &_mtrnew);
	static void CombineMatricesSort2Index (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result (2index version)
															const CSMatrix *_mtrarr,
															CSMatrix &_mtrnew);
	static void SymmetrizeMatrix2Index (int _nblks, int *_blks, // Symmetrize the matrix
														int _nlist, int *_list, const CSMatrix *_mtrarr, 
														int &_nlistnew, int *_listnew, CSMatrix *_mtrarrnew);
// Multiplication functions
	virtual void SolveL   (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.GetNv ();
		i = _lx.GetNv ();
		std::cout << " In SolveL   CSMatrix " << std::endl; 
	};
	virtual void SolveU   (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.GetNv ();
		i = _ux.GetNv ();
		std::cout << " In SolveU   CSMatrix " << std::endl; 
	};
	virtual void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.GetNv ();
		i = _lx.GetNv ();
		std::cout << " In SolveL   CSMatrix " << std::endl; 
	};
	virtual void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.GetNv ();
		i = _ux.GetNv ();
		std::cout << " In SolveU   CSMatrix " << std::endl; 
	};
// Compute ordering functions
// Input/Output functions
	void SetupFiles (int _nfiles, char *_name); // Setup files
	void CloseFiles (); // Close files
	void FreeFiles (); // Free files
	friend std::ostream &operator<< (std::ostream &_stream, const CGSMatrix &_a); // Output matrix
// Friend classes
	friend class CFctR;
	friend class CFctDiagR;
	friend class CFctC;
	friend class CFctDiagC;
	friend class CMvm;
	friend class CMvmR;
	friend class CMvmC;
#ifndef __FV_2004__
	friend class CModel;
#endif
	friend class CSVector;
};

#endif

#ifndef __GSMatrixR
#define __GSMatrixR

class CGSMatrixR: public CGSMatrix
{
	CSMatrixR *mtrarr; // mtrarr[nblksr] array contains the set of CSMatrixR structures that specify each block row of the matrix
	int ivect;          // ivect describes the type of additional vector data
	CSVector vector1;   // vector1 contains some additional vector data
	CSVector vector2;   // vector2 contains some additional vector data
public:
// Constructors and destructor
	CGSMatrixR (); // Memory allocation zero data constructor
	CGSMatrixR (int _nblks, int _nsups, int _nlist); // Memory allocation zero data constructor
	~CGSMatrixR (); // Destructor
// Operator functions
	CGSMatrixR &operator= (const CGSMatrixR &_aa); // Equality operator
// Get/set functions
	CSMatrixR * GetMtrarr () const {return mtrarr;}; // Get mtrarr
	void SetMtrarr (CSMatrixR *_mtrarr) {mtrarr = _mtrarr;};   // Set mtrarr
// Functions
	void TransformIntoHyperSparseFormat (CTree &_tree); // Transform secondary matrix into the hypersparse format
	void R2GR2Ind (const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only
					int _myid, int *_bl2cpu);
	void R2GRSubmatrix2Ind (const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only, matrix contains only necessary rows
									int _myid, int *_bl2cpu);
	void R2GRSubmatrix2IndFastSearch (bool _is_all_rows, const CSMatrixR &_mtra, int _nblks, int *_blks, // Create GSMatrixR structure via SMatrixR for specified indices only, matrix contains only necessary rows
												int _myid, int *_bl2cpu);
	CSMatrixR CombineMatrices2Index (); // Combine matrices into one
	CSMatrixR CombineMatrices2IndexInto2Index (); // Combine matrices into one
	static void CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												const CSMatrixR *_mtrarr,
												CSMatrixR &_mtrnew);
	static void Ich2StoreScaleUDynamic (int _nblks, int *_blks, double *_scl, // Store the resulting U and scale it
								int *_ig, int **_jgblk, double **_gblk,
								CGSMatrixR &_gich);
	void Ilu2Index (std::ofstream &_fout, const CSlvParam &_param, // Second order Incomplete LU decomposition
						const CGSMatrixR &_mtral, const CGSMatrixR &_mtrau,
						CGSMatrixR &_mtrl, CGSMatrixR &_mtru);
	void Ilu2Schur2Index (std::ofstream &_fout, const CTree &_tree, CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
								CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
								CGSMatrixR &_mtrl, CGSMatrixR &_mtru);
	void Ilu2SchurDynamic2Index (std::ofstream &_fout, const CTree &_tree, CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode (dynamic version)
											CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
											CGSMatrixR &_mtrl, CGSMatrixR &_mtru);
	void SolveL (const CSVector &_x, CSVector &_lx) const; // Solve triangular system with L stored by block columns
	void SolveU (const CSVector &_x, CSVector &_ux) const; // Solve triangular system with U stored by block rows
	void SolveL (const CSVectorC &_x, CSVectorC &_lx) const {}; // Solve triangular system with L stored by block columns
	void SolveU (const CSVectorC &_x, CSVectorC &_ux) const {}; // Solve triangular system with U stored by block rows
	void MvmA2Index (const CTree &_tree, CMvmR &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
							const CSVector &_x, CSVector &_ax) const;
	void MvmA2Index2 (const CTree &_tree, CMvmR &_mvm, // General secondary matrix by vector multiplication in parallel mode
							const CSVector &_x, CSVector &_ax);
	void MvmAh2Index2 (const CTree &_tree, CMvmR &_mvm, // General transposed secondary matrix by vector multiplication in parallel mode
									const CSVector &_x, CSVector &_ax);
	void SolveL2Index (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with L stored by block columns in parallel mode
							const CSVector &_x, CSVector &_lx);
	void SolveU2Index (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with U stored by block rows in parallel mode
							const CSVector &_x, CSVector &_ux);
	void SplitMatrix2Index (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
							const CSMatrixR &_mtr, 
							CSMatrixR *_mtrarr);
	void SplitMatrixList2Index (int &_nlistblk, int *_listblk, // Split matrix into the set of matrices
											const CSMatrixR &_mtr, 
											CSMatrixR *_mtrarr);
	void SendMatrix2Index (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixR &_mtra) const; // Send matrix
	void ReceiveMatrix2Index (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixR &_mtra) const; // Send matrix
	void FilterSchurBlock2Index (int _iblk, // Filter Schur block row/column according to the indices
									int _iblkbeg, int _iblkend, int _iblkbegrest, double _theta,
									CFctR &_fct);
	void CombineMatrices2Index (bool _cpyindex3, int _nlist, int *_list, // Combine matrices into one
											const CSMatrixR *_mtrarr,
											CSMatrixR &_mtrnew) const;
	static void ReorderDiagonalPart2Index (int _nblks, int *_blks, // Reorder diagonal part of the matrix
													int _nlist, int *_list, const CSMatrixR *_mtrarr, 
													int _iblkbeg, int _iblkend, int *_order,
													CSMatrixR *_mtrarrnew);
	static void SplitLU2Index (int _nblks, int *_blks, // Split matrix into L and U parts
										int _nlist, int *_list, const CSMatrixR *_mtrarr, 
										CSMatrixR *_mtrlarr, CSMatrixR *_mtruarr);
	void AddReplaceSetOfBlocks2Index (CFctR &_fct, int _nlistblk, int *_listblk, // Add and replace set of block rows/columns
													CSMatrixR *_mtrlarr, CSMatrixR *_mtrl2arr, CSMatrixR *_mtrarru, CSMatrixR *_mtrarru2);
	void FilterListOfSchurBlocks2Index (int _nlistblk, int *_listblk, // Filter list of Schur blocks according to the block indices
													int *_iabextflt, int *_jabextflt, double _theta,
													CFctR &_fct);
	CSMatrixR CreateDiagonalUpdate2Index (int _nlistblk, int *_listblk, // Create diagonal update matrix
														CFctDiagR &_fctdiag);
	void UpdateDiagonal2Index (CSMatrixR &_mtrd, // Update diagonal data
										CFctDiagR &_fctdiag);
	void FctSubTree2Index (const CTree &_treeini, const CTree &_tree, // Perform parallel computations with the current subtree
									CFctR &_fct, CFctDiagR &_fctdiag,
									int *_bl2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
									int _nlistschur, int *_listschur, int *_lev2nodes,
									CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
									CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
									CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
									CSlvParam &_param);
	void SchurUpdates2Index (bool _last_work_with_Schur, // Perform Schur updates via exchanges from all involved processors
										int _ilev, const CTree &_treeini, const CTree &_tree,
										CFctR &_fct, CFctDiagR &_fctdiag,
										int *_blk2cpu,
										int *_iabextflt, int *_jabextflt,
										int _nlistschur, int *_listschur,
										int *_lev2nodearr,
										CGSMatrixR &_mtralsch, CGSMatrixR &_mtrausch,
										CSlvParam &_param);
	void FctNode2Index (bool _update_by_Schur, int _inode, // Perform computations with the current tree node including local Schur updates
								const CTree &_treeini, const CTree &_tree,
								CFctR &_fct, CFctDiagR &_fctdiag,
								int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
								int _nlistschur, int *_listschur,
								CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
								CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
								CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
								CSlvParam &_param);
	void PrepareDynamicArrays2Index (const CTree &_tree, // Prepare all arrays for dynamic parallel preconditioner computation
												int &_n_lev_work, int *&_lev_work, 
												int *&_i_beg_work, int *&_i_end_work, int *&_i_beg_schur,
												int _nblks, int *_iabextt, int _nlistschur, int *_listschur,
												int *&_iptrextt, int *&_updblkschur,
												int &_nsendmax, int *&_imasksend, char **&_psenddata, CMPIRequest *&_sendreqvarr, 
												int *&_imaskrecv, char **&_precvdata, CMPIRequest *&_recvreqvarr, 
												CSMatrixR *&_mtrarrwork, CSMatrixR *&_mtrarr2work);
	void FctNode2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
								int *_i_beg_work, int *_i_end_work,
								const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag,
								int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
								CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
								CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
								CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
								CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd, 
								CSlvParam &_param);
	void PriorityUpdates2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
											int *_i_beg_work, int *_i_end_work, int _iblkbeg,
											int &_remained_updates, 
											const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
											int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
											CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
											CSlvParam &_param);
	void NonPriorityUpdates2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
												int *_i_beg_work, int *_i_end_work,
												int &_ip_work, int _ip_work_end, int *_listschur, int *_updblkschur,
												const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
												int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
												CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
												CSlvParam &_param);
	void UpdateIblk2Index (int _iblk, int _iblkbeg, // Perform computations with the current iblk
									const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
									int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
									CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
									CSlvParam &_param);
	void AsyncSendBlocks2Index (int _i_lev_work, int *_lev_work, // Send data to the other cpu's
											int *_i_beg_work, int *_i_end_work, int _nlistschur, int *_listschur,
											const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag,
											int _id_shift, int _nsendsmax, int &_nsends, int *_imasksend, char **_psenddata, CMPIRequest *_sendreqvarr,
											int *_blk2cpu,
											int *_iabextflt, int *_jabextflt,
											int *_iabextt, int *_iptrextt, int *_jabextt,
											CGSMatrixR &_mtralsch, CGSMatrixR &_mtrausch,
											CSlvParam &_param);
	void InitReceiveBlocks2Index (int _i_lev_work, int *_lev_work, const CTree &_tree, // Init to receive some data
											int &_nrecvs, int *_imaskrecv, char **_precvdata, CMPIRequest *_recvreqvarr);
	void TestSendsRecvsAndUpdate2Index (const CTree &_tree, // Test/wait sends and receives for completion and update if necessary
													bool _remained_work,
													int _id_shift, int _nsends, int &_isend, 
													int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr,
													int &_nrecvs, int &_irecv,
													int *_imaskrecv, char** _precvdata, CMPIRequest *_recvreqvarr,
													CFctR &_fct, CFctDiagR &_fctdiag,
													CGSMatrixR &_mtrlsch, CGSMatrixR &_mtrusch, CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd, 
													CSMatrixR *_mtrarrwork, CSMatrixR *_mtrarr2work);
	void UseReceivedData2Index (const CTree &_tree, char *_buffer, CFctR &_fct, CFctDiagR &_fctdiag, // Use receive data for updates
											CGSMatrixR &_mtrlsch, CGSMatrixR &_mtrusch, CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd, 
											CSMatrixR *_mtrarrwork, CSMatrixR *_mtrarr2work);
	void CheckSends2Index (int _nsends, int &_isend, // Check sends for completion
									int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr);
	void WaitSends2Index (int _nsends, int _isend, // Wait sends for completion
									int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr);
	void FreeAndCopySparsityPointers (int _iblk, CSMatrix &_mtr); // Free pointers of the current block and copy them from given matrix
	void RestoreSparsityPointers (int _iblk); // Restore sparsity pointers of the current block
// Input/Output functions
	void FilterBlock2Index (int _iblk, int _iblkend); // Filter block row/column according to the index
	void WriteBlock    (int _iblk); // Write block row/column into the disk
	void ReWriteBlock  (int _iblk); // Rewrite block row/column into the disk
	void ReadBlock     (int _iblk); // Read block row/column from the disk into the main memory if necessary
	void AllocateBlock (int _iblk); // Allocate block row/column
	void FreeBlock     (int _iblk); // Free block row/column from main memory
	friend std::ostream &operator<< (std::ostream &_stream, const CGSMatrixR &_a); // Output matrix
// Friend classes
	friend class CFctR;
	friend class CFctDiagR;
	friend class CSVector;
};

#endif

#ifndef __GSMatrixRS
#define __GSMatrixRS

class CGSMatrixRS: public CGSMatrix
{
protected:
// Data
	CSMatrixRS *mtrarr; // mtrarr[nblksr] array contains the set of CSMatrixRS structures that specify each block row of the matrix
	int ivect;          // ivect describes the type of additional vector data
	CSVector vector1;   // vector1 contains some additional vector data
	CSVector vector2;   // vector2 contains some additional vector data
public:
// Functions
// Constructors and destructor
	CGSMatrixRS (); // Memory allocation zero data constructor
	CGSMatrixRS (int _nblks, int _nsups, int _nlist); // Memory allocation zero data constructor
	CGSMatrixRS (const CGSMatrixRS &_aa) : CGSMatrix(_aa)
		{throw " Copy constructor for CGSMatrixRS class called";}; // Copy constructor
	~CGSMatrixRS (); // Destructor
// Get/set functions
	CSMatrixRS * GetMtrarr () const {return mtrarr;}; // Get mtrarr
	void SetMtrarr (CSMatrixRS *_mtrarr) {mtrarr = _mtrarr;};   // Set mtrarr
// Operator functions
	CGSMatrixRS &operator= (const CGSMatrixRS &_aa); // Equality operator
// Transformation functions
	void SupersparsifyBlockRows (); // Supersparsify block rows
	void RS2GRS (const CSMatrixRS &_mtra, int _nblks, int *_blks); // Create GSMatrixRS structure via SMatrixRS
	void RS2GRS (bool _is_all_rows, const CSMatrixRS &_mtra, int _nblks, int *_blks, // Create GSMatrixRS structure via SMatrixRS for specified indices only
					int _myid, int *_bl2cpu);
	void FullFormat (int _iblk, int _ntotal, double *_a); // Write current block into full format
	CSMatrix GlobalSparsity () const; // Compute the sparsity structure of all local block rows
	void CopyBlk2Blk  (char _diatype, const CSMatrixRS &_a, CSMatrixRS &_anew) const; // Copy a part of data of the block into the new one
	void CopyBlkT2Blk (char _diatype, int &_icycle, int *_imask, 
						const CSMatrixRS &_a, CSMatrixRS &_anew) const; // Copy a part of data of the block into the new one
// Ordering and scaling functions
	void PointScaling (const CSVector &_sclpt); // Perform explicit point scaling of the matrix
	void BlockScaling (const CSVector &_sclblkl, const CSVector &_sclblku); // Perform explicit block scaling of the matrix
// Compute ordering functions
// Send/receive functions
	void CombineMatrices (int _nlist, int *_list, // Combine matrices into one
								const CSMatrixRS *_mtrarr,
								CSMatrixRS &_mtrnew) const;
	static void CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												CSMatrixR *_mtrarr,
												CSMatrixR &_mtrnew);
	static void CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
												int _nsup, int *_sprnds,
												CSMatrixRS *_mtrarr,
												CSMatrixRS &_mtrnew);
	void SplitMatrix (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
						const CSMatrixRS &_mtr, 
						CSMatrixRS *_mtrarr);
	void PrepareSendMatrix (int _inode, const CTree &_tree, // Prepare the matrix for exchange
							CGSMatrixRS &_gmtru2, const CSMatrixRS *_mtru2charr,
							CSMatrixRS &_mtru2send) const;
	void SendMatrix    (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixRS &_mtra) const; // Send matrix
	void ReceiveMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixRS &_mtra) const; // Receive matrix
// Incomplete factorization functions
	void AtAMatr (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute complex A^hA matrix in parallel mode
					CGSMatrixRS &_mtral, CGSMatrixRS &_mtrau,
					CGSMatrixRS &_mtratau);
	void FilterSchurBlock (int _iblk, // Filter Schur block row/column according to the indices
							int _isupbeg, int _isupend, int _isupbegrest, double _theta,
							CFctR &_fct);
	void Ilu2 (std::ofstream &_fout, const CSlvParam &_param, // Second order Incomplete LU decomposition
				const CGSMatrixRS &_mtral, const CGSMatrixRS &_mtrau,
				CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru);
	void Ilu2 (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
				CGSMatrixRS &_mtral, CGSMatrixRS &_mtrau,
				CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru);
	void Ich2Schur (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete Cholessky decomposition in parallel mode
					CGSMatrixRS &_mtrau,
					CGSMatrixRS &_mtru);
	void Ilu2Schur (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
					CGSMatrixRS &_mtral, CGSMatrixRS &_mtrau,
					CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru);
	void FindSups (std::ofstream &_fout, int _nthresh, double *_thresh, int *_iorder); // Find supernodes with small singular values
// Matrix by vector multiplication and triangular solve functions
	void MvmA (const CGSMatrixRS &_gmtral, const CGSMatrixRS &_gmtrau, // General matrix by vector multiplication
				const CSVector &_x, CSVector &_ax) const;
	void MvmA (const CTree &_tree, CMvmR &_mvm, // General matrix by vector multiplication in parallel mode
				CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau, 
				const CSVector &_x, CSVector &_ax) const;
	void MvmA2 (const CTree &_tree, CMvmR &_mvm, // General second order matrix by vector multiplication in parallel mode
				const CSVector &_x, CSVector &_ax);
	void MvmAh (const CTree &_tree, CMvmR &_mvm, // General hermitian transposed matrix by vector multiplication in parallel mode
				CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau, 
				const CSVector &_x, CSVector &_ax) const;
	void MvmAh2 (const CTree &_tree, CMvmR &_mvm, // General second order matrix by vector multiplication in parallel mode
				const CSVector &_x, CSVector &_ax);
	void SolveL (const CSVector &_x, CSVector &_lx) const; // Solve triangular system with L stored by block columns
	void SolveL (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with L stored by block columns in parallel mode
					const CSVector &_x, CSVector &_lx);
	void SolveLConj (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with conj(L) stored by block columns in parallel mode
						const CSVector &_x, CSVector &_lx);
	void SolveU (const CSVector &_x, CSVector &_ux) const; // Solve triangular system with U stored by block rows
	void SolveU (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with U stored by block rows in parallel mode
					const CSVector &_x, CSVector &_ux);
	void SolveUConj (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with conj(U) stored by block rows in parallel mode
						const CSVector &_x, CSVector &_ux);
	void SolveL (const CSVectorC &_x, CSVectorC &_lx) const {}; // Solve triangular system with L stored by block columns
	void SolveU (const CSVectorC &_x, CSVectorC &_ux) const {}; // Solve triangular system with U stored by block rows
// Input/Output functions
	void FilterBlock   (int _iblk, int _isupend); // Filter block row/column according to the index
	void WriteBlock    (int _iblk); // Write block row/column into the disk
	void ReWriteBlock  (int _iblk); // Rewrite block row/column into the disk
	void ReadBlock     (int _iblk); // Read block row/column from the disk into the main memory if necessary
	void AllocateBlock (int _iblk); // Allocate block row/column
	void FreeBlock     (int _iblk); // Free block row/column from main memory
	friend std::ostream &operator<< (std::ostream &_stream, const CGSMatrixRS &_a); // Output matrix
// Friend classes
	friend class CFctR;
	friend class CFctDiagR;
	friend class CSVector;
#ifndef __FV_2004__
	friend class CModel;
#endif
};

#endif

#ifndef __GSMatrixCS
#define __GSMatrixCS

class CGSMatrixCS: public CGSMatrix
{
protected:
// Data
	CSMatrixCS *mtrarr; // mtrarr[nblksr] array contains the set of CSMatrixCS structures that specify each block row of the matrix
	int ivect;          // ivect describes the type of additional vector data
	CSVectorC vector1;  // vector1 contains some additional vector data
	CSVectorC vector2;  // vector2 contains some additional vector data
public:
// Functions
// Constructors and destructor
	CGSMatrixCS (); // Memory allocation zero data constructor
	CGSMatrixCS (int _nblks, int _nsups, int _nlist); // Memory allocation zero data constructor
	CGSMatrixCS (const CGSMatrixCS &_aa) : CGSMatrix(_aa)
	{throw " Copy constructor for CGSMatrixCS class called";}; // Copy constructor
	~CGSMatrixCS (); // Destructor
// Get/set functions
	CSMatrixCS *GetMtrarr () const {return mtrarr;}; // Get mtrarr
// Operator functions
	CGSMatrixCS &operator= (const CGSMatrixCS &_aa); // Equality operator
// Transformation functions
	void CS2GCS (const CSMatrixCS &_mtra, int _nblks, int *_blks); // Create GSMatrixCS structure via SMatrixCS
	void CS2GCS (const CSMatrixCS &_mtra, int _nblks, int *_blks, // Create GSMatrixCS structure via SMatrixCS for specified indices only
					int _myid, int *_bl2cpu);
	void CS2GCS (bool _is_all_rows, const CSMatrixCS &_mtra, int _nblks, int *_blks, // Create GSMatrixCS structure via SMatrixCS for specified indices only
					int _myid, int *_bl2cpu);
	void FullFormat (int _iblk, int _ntotal, dcmplx *_a); // Write current block into full format
	CSMatrix GlobalSparsity () const; // Compute the sparsity structure of all local block rows
	void CopyBlk2Blk  (char _diatype, const CSMatrixCS &_a, CSMatrixCS &_anew) const; // Copy a part of data of the block into the new one
	void CopyBlkT2Blk (char _diatype, int &_icycle, int *_imask, 
						const CSMatrixCS &_a, CSMatrixCS &_anew) const; // Copy a part of data of the block into the new one
	CSMatrixCS Sparsify (int _nitermax, int _nza); // Compute the point sparse matrix via very dense rectangular matrix
// Ordering and scaling functions
	void PointScaling (const CSVectorC &_sclpt); // Perform explicit point scaling of the matrix
	void BlockScaling (const CSVectorC &_sclblkl, const CSVectorC &_sclblku); // Perform explicit block scaling of the matrix
	void ExplicitPointScaling (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute point scaling of the rectangular matrix A and scale explicitely
										CGSMatrixCS &_mtra, int *_iamask, int *_jamask,
										CCorrectorC &_scla);
	void ExplicitSymmetricPointScaling (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute point scaling of the rectangular matrix A and scale explicitely
													CGSMatrixCS &_mtra, CCorrectorC &_scla);
	void ExplicitSymmetricBlockScaling (std::ofstream &_fout, const CSlvParam &_param, // Perform in-place symmetric block scaling of the square matrix stored by block columns
													CGSMatrixCS &_mtra, CSVectorC &_rhs,
													int *_sprndsflt, CCorrectorC &_scla);
	void ExplicitBlockScaling (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute block scaling of the rectangular matrix A and scale explicitely
											CGSMatrixCS &_mtra, int *_iamask, int *_jamask,
											CCorrectorC &_scl);
	void FilterLinearSystem (// Perform in-place filtering of the system of equations
										CGSMatrixCS &_mtra, CSVectorC &_rhs,
										int *_sprndsflt);
// Compute ordering functions
// Send/receive functions
	void CombineMatrices (int _nlist, int *_list, // Combine matrices into one
							const CSMatrixCS *_mtrarr,
							CSMatrixCS &_mtrnew) const;
	static void CombineMatricesSort (int _n, int _nlist, int *_list, // Combine matrices into one and sort the result
											int _nsup, int *_sprnds,
											CSMatrixCS *_mtrarr,
											CSMatrixCS &_mtrnew);
	void SplitMatrix (int _nlist, int *_list, // Split matrix into the set of matrices according to the list
						const CSMatrixCS &_mtr, 
						CSMatrixCS *_mtrarr);
	void SplitMatrixList (int &_nlistblk, int *_listblk, int *_sp2blk, // Split matrix into the set of matrices
											const CSMatrixCS &_mtr, 
											CSMatrixCS *_mtrarr);
	CSMatrixC CreateDiagonalUpdate2Index (int _nlistblk, int *_listblk, // Create diagonal update matrix
														CFctDiagC &_fctdiag);
	void UpdateDiagonal2Index (CSMatrixC &_mtrd, // Update diagonal data
										CFctDiagC &_fctdiag);
	void PrepareSendMatrix (int _inode, const CTree &_tree, // Prepare the matrix for exchange
							CGSMatrixCS &_gmtru2, const CSMatrixCS *_mtru2charr,
							CSMatrixCS &_mtru2send) const;
	void SendMatrix    (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixCS &_mtra) const; // Send matrix
	void ReceiveMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixCS &_mtra) const; // Receive matrix
// Incomplete factorization functions
	void ComputeSymmetricBlockScaling (int _iblk, double _threshflt, // Compute scaling of the current block
													CHyst &_hyst, CGSMatrixCS &_mtra, int *_sprndsflt, CCorrectorC &_scla);
	void SymmetricScaleBlockColumn (int _iblk, // Perform symmetric scaling of the current block
												CGSMatrixCS &_mtra, int *_sp2blk, CSVectorC &_rhs,
												int *_ibsscl, CCorrectorC &_scla);
	void ScaleBlockRows (std::ofstream &_fout, const CSlvParam &_param, // Perform in-place scaling of the block rows of the square matrix stored by block columns
								CSMatrix &_mtrasp, CGSMatrixCS &_mtra, CSVectorC &_rhs);
	void FilterBlockRows (std::ofstream &_fout, const CSlvParam &_param, // Perform in-place filtering of the block rows of the rectangular matrix
									CSMatrix &_mtrasp, CGSMatrixCS &_mtra, CSVectorC &_rhs,
									int *_iamask, int *_jamask);
	void IndicesForFilteringOrScaling (int _iblk, // Compute indices of the nonlocal data for block rows filtering or scaling
										CSMatrix &_mtrasp, CGSMatrixCS &_mtra, 
										int *_jamask, int *_sp2blk,
										int &_icycle, int *_imask,
										CSMatrix &_mtrasptran, int *_indarr,
										CSMatrixCS &_mtraoff, 
										int *&_ja2isup, int *&_ja2blk, int *&_ja2ind, int *&_ja2bs);
	void IndicesForScaling (int _iblk, // Compute indices of the nonlocal data for block rows scaling
									CSMatrix &_mtrasp, CGSMatrixCS &_mtra, 
									int *_jamask, int *_sp2blk,
									int &_icycle, int *_imask,
									CSMatrix &_mtrasptran, int *_indarr,
									CSMatrixCS &_mtraoff, 
									int *&_ja2isup, int *&_ja2blk, int *&_ja2ind, int *&_ja2bs);
	void ReadBlockRows ( // Read some block rows
								CGSMatrixCS &_mtra, 
								CSMatrixCS &_mtraoff, 
								int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs,
								int &_icycle, int *_imask);
	void WriteBlockRows ( // Write some block rows
								CGSMatrixCS &_mtra, 
								CSMatrixCS &_mtraoff, 
								int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs,
								int &_icycle, int *_imask);
	void FilterBlockRows (int _iblk, // Filter some block rows
									CSMatrix &_mtrasp, CSMatrix &_mtrasptran, int *_indarr,
									CGSMatrixCS &_mtra, int *_jamask, int *_sp2blk,
									CSMatrixCS &_mtraoff, CSVectorC &_rhs,
									int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs);
	void ScaleBlockRows (int _iblk, // Scale some block rows
												CSMatrix &_mtrasp, CSMatrix &_mtrasptran, int *_indarr,
												CGSMatrixCS &_mtra, int *_jamask, int *_sp2blk,
												CSMatrixCS &_mtraoff, CSVectorC &_rhs,
												int *_ja2isup, int *_ja2blk, int *_ja2ind, int *_ja2bs);
	void AtAMatr (std::ofstream &_fout, const CSlvParam &_param, // Compute complex A^hA matrix for rectangular A in serial mode
						CGSMatrixCS &_mtra,
						CSMatrix &_mtratasp,
						CGSMatrixCS &_mtratau);
	void AtAMatr (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute complex A^hA matrix for rectangular A in parallel mode
						CSMatrix &_mtrasp, CGSMatrixCS &_mtra,
						CSMatrix &_mtratasp,
						CGSMatrixCS &_mtratau);
	void AtAMatr (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute complex A^hA matrix in parallel mode
					CGSMatrixCS &_mtral, CGSMatrixCS &_mtrau,
					CGSMatrixCS &_mtratau);
	void UpdateAtAMatrixBlock (int _iblk, int _jblk, // Compute current block of the complex A^hA matrix for rectangular A
										 CGSMatrixCS &_mtra, int *_sp2blk);
	void FilterSchurBlock (int _iblk, // Filter Schur block row/column according to the indices
							int _isupbeg, int _isupend, int _isupbegrest, double _theta,
							CFctC &_fct);
	static void Ich2 (std::ofstream &_fout, const CSlvParam &_param, // Second order ICH decomposition
						const CGSMatrixCS &_mtrau,
						CGSMatrixCS &_mtru);
	void Ilu2 (std::ofstream &_fout, const CSlvParam &_param, // Second order Incomplete LU decomposition
				const CGSMatrixCS &_mtral, const CGSMatrixCS &_mtrau,
				CGSMatrixCS &_mtrl, CGSMatrixCS &_mtru);
	void Ilu2 (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
				CGSMatrixCS &_mtral, CGSMatrixCS &_mtrau,
				CGSMatrixCS &_mtrl, CGSMatrixCS &_mtru);
	void Ich2Schur (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete Cholessky decomposition in parallel mode
					CGSMatrixCS &_mtrau,
					CGSMatrixCS &_mtru);
	void Ilu2Schur (std::ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
					CGSMatrixCS &_mtral, CGSMatrixCS &_mtrau,
					CGSMatrixCS &_mtrl, CGSMatrixCS &_mtru);
	void FindSups (std::ofstream &_fout, int _nthresh, double *_thresh, int *_iorder); // Find supernodes with small singular values
// Matrix by vector multiplication and triangular solve functions
	void VectorPointScale (int _ivar,  // Perform explicit point scaling of the solution
									CGSMatrixCS &_mtra,
									CCorrectorC &_scl, CSVectorC &_sol);
	void VectorBlockScale (int _ivar,  // Perform explicit block scaling of the solution
									CGSMatrixCS &_mtra,
									CCorrectorC &_scl, CSVectorC &_sol);
	void ExtendedVectorBlockScale (int _ivar,  // Perform explicit block scaling of the solution with extension
												CGSMatrixCS &_mtraflt,
												int *_sprndsini, CCorrectorC &_scl, 
												CSVectorC &_sol, CSVectorC &_solext);
	void MvmA (const CGSMatrixCS &_gmtral, const CGSMatrixCS &_gmtrau, // General matrix by vector multiplication
				const CSVectorC &_x, CSVectorC &_ax) const;
	void MvmASymm ( // Complex hermitian matrix by vector multiplication
						const CSVectorC &_x, CSVectorC &_ax);
	void MvmASymm (const CTree &_tree, CMvmC &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixCS &_gmtrau, 
							const CSVectorC &_x, CSVectorC &_ax) const;
	void MvmA (const CTree &_tree, CMvmC &_mvm, // General matrix by vector multiplication in parallel mode
				CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
				const CSVectorC &_x, CSVectorC &_ax) const;
	void MvmAh (const CTree &_tree, CMvmC &_mvm, // General hermitian transposed matrix by vector multiplication in parallel mode
				CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
				const CSVectorC &_x, CSVectorC &_ax) const;
	void MvmARect ( // General rectangular matrix by vector multiplication
							const CSVectorC &_x, CSVectorC &_ax);
	void MvmAhRect ( // General complex conjugate transposed rectangular matrix by vector multiplication
							const CSVectorC &_x, CSVectorC &_ax);
	void MvmAhRect (const CTree &_tree, // General complex conjugate transposed rectangular matrix by vector multiplication
							const CSVectorC &_x, CSVectorC &_ax);
	void MvmAhAOrdShiftSymm (double _dshift, int *_order,  // Complex hermitian matrix by vector multiplication
								const CSVectorC &_x, CSVectorC &_ax);
	void SolveL (const CSVectorC &_x, CSVectorC &_lx) const; // Solve triangular system with L stored by block columns
	void SolveL (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with L stored by block columns in parallel mode
					const CSVectorC &_x, CSVectorC &_lx);
	void SolveLConj (const CSVectorC &_x, CSVectorC &_lx); // Solve triangular system with L stored by block columns
	void SolveLConj (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with conj(L) stored by block columns in parallel mode
						const CSVectorC &_x, CSVectorC &_lx);
	void SolveU (const CSVectorC &_x, CSVectorC &_ux) const; // Solve triangular system with U stored by block rows
	void SolveU (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with U stored by block rows in parallel mode
					const CSVectorC &_x, CSVectorC &_ux);
	void SolveUConj (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with conj(U) stored by block rows in parallel mode
						const CSVectorC &_x, CSVectorC &_ux);
	void SolveL (const CSVector &_x, CSVector &_lx) const {}; // Solve triangular system with L stored by block columns
	void SolveU (const CSVector &_x, CSVector &_ux) const {}; // Solve triangular system with U stored by block rows
// Implementation of the abstract interface functions
	static int GETM (void *_obj) { // Define interface of the columns object size
		CGSMatrixCS *objloc = (CGSMatrixCS *) _obj;
		int mloc = objloc->GetM ();
		return mloc;
	};
	static int GETN (void *_obj) { // Define interface of the columns object size
		CGSMatrixCS *objloc = (CGSMatrixCS *) _obj;
		int nloc = objloc->GetN ();
		return nloc;
	};
	static void MVMC (void *_obj, CSVectorC &_x, CSVectorC &_ax) { // Define interface of the matrix by vector multiplication
		CGSMatrixCS *objloc = (CGSMatrixCS *) _obj;
		objloc->MvmARect (_x, _ax);
	};
	static void MVMHC (void *_obj, CSVectorC &_x, CSVectorC &_ax) { // Define interface of the matrix by vector multiplication
		CGSMatrixCS *objloc = (CGSMatrixCS *) _obj;
		objloc->MvmAhRect (_x, _ax);
	};
	static CSMatrixCS GETSPARSITYCS (void *_obj, CSlvParam &_params) { // Define interface of the sparse matrix generation
		CGSMatrixCS *objloc = (CGSMatrixCS *) _obj;
		int nloc = objloc->GetN ();
		int nitermax = 10;
		int ndiag = _params.ndiag;
		return objloc->Sparsify (nitermax, nloc*ndiag);
	};
// Input/Output functions
	void FilterBlock   (int _iblk, int _isupend); // Filter block row/column according to the index
	void WriteBlock    (int _iblk); // Write block row/column into the disk
	void ReWriteBlock  (int _iblk); // Rewrite block row/column into the disk
	void ReadBlock     (int _iblk); // Read block row/column from the disk into the main memory if necessary
	void AllocateBlock (int _iblk); // Allocate block row/column
	void FreeBlock     (int _iblk); // Free block row/column from main memory
	friend std::ostream &operator<< (std::ostream &_stream, const CGSMatrixCS &_a); // Output matrix
// Friend classes
	friend class CFctC;
	friend class CFctDiagC;
	friend class CSVectorC;
#ifndef __FV_2004__
	friend class CModel;
#endif
};

#endif

#ifndef __ParMvmA
#define __ParMvmA

// Preliminary declarations

class CTree;
class CMvmR;

class CParMvmA
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CMvmR *pmvm;     // pmvm contains pointer to mvm structure
	CGSMatrixRS *pal; // pl contains pointer to L part
	CGSMatrixRS *pau; // pu contains pointer to U part
	bool second;      // second specifies second order data
	CMvmR *pmvm2;     // pmvm2 contains pointer to mvm2 structure
	CGSMatrixRS *pa2; // pa contains pointer to A2 part
public:
// Functions
// Constructors and destructor
	CParMvmA (CTree *_ptree, CMvmR *_pmvm, CGSMatrixRS *_pal, CGSMatrixRS *_pau, // Memory allocation zero data constructor
					bool _second, CMvmR *_pmvm2, CGSMatrixRS *_pa2) { 
		ptree = _ptree;
		pmvm = _pmvm;
		pal = _pal;
		pau = _pau;
		second = _second;
		pmvm2 = _pmvm2;
		pa2 = _pa2;
	};
	CParMvmA (const CParMvmA &_mvm) {throw " Copy constructor for CParMvmA class called";}; // Copy constructor
	~CParMvmA () {}; // Destructor
// Operator functions
	CParMvmA &operator= (const CParMvmA &_mvm) {throw " Equality operator for CParMvmA class called";}; // Equality operator
// Solve function
	static void MvmA (void *_objlu, CSVector &_x, CSVector &_px) {
		CParMvmA *pobjloc = (CParMvmA *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *palloc = pobjloc->pal;
		CGSMatrixRS *pauloc = pobjloc->pau;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixRS *pa2loc = pobjloc->pa2;
		CGSMatrixRS gmtrdummy;
		gmtrdummy.MvmA (*ptreeloc, *pmvmloc, *palloc, *pauloc, _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmA2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
	static void MvmAt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParMvmA *pobjloc = (CParMvmA *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *palloc = pobjloc->pal;
		CGSMatrixRS *pauloc = pobjloc->pau;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixRS *pa2loc = pobjloc->pa2;
		CGSMatrixRS gmtrdummy;
		gmtrdummy.MvmA (*ptreeloc, *pmvmloc, *pauloc, *palloc, _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmAh2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParMvmA2Index
#define __ParMvmA2Index

// Preliminary declarations

class CTree;
class CMvmR;

class CParMvmA2Index
{
public:
// Data
	CTree * ptree;    // ptree contains pointer to the tree
	CMvmR *pmvm;      // pmvm contains pointer to mvm structure
	CGSMatrixR *pal;  // pl contains pointer to L part
	CGSMatrixR *pau;  // pu contains pointer to U part
	bool second;      // second specifies second order data
	CMvmR *pmvm2;     // pmvm2 contains pointer to mvm2 structure
	CGSMatrixR *pa2;  // pa contains pointer to A2 part
public:
// Functions
// Constructors and destructor
	CParMvmA2Index (CTree *_ptree, CMvmR *_pmvm, CGSMatrixR *_pal, CGSMatrixR *_pau, // Memory allocation zero data constructor
							bool _second, CMvmR *_pmvm2, CGSMatrixR *_pa2) { 
		ptree = _ptree;
		pmvm = _pmvm;
		pal = _pal;
		pau = _pau;
		second = _second;
		pmvm2 = _pmvm2;
		pa2 = _pa2;
	};
	CParMvmA2Index (const CParMvmA2Index &_mvm) {throw " Copy constructor for CParMvmA2Index class called";}; // Copy constructor
	~CParMvmA2Index () {}; // Destructor
// Operator functions
	CParMvmA2Index &operator= (const CParMvmA2Index &_mvm) {throw " Equality operator for CParMvmA2Index class called";}; // Equality operator
// Solve function
	static void MvmA (void *_objlu, CSVector &_x, CSVector &_px) {
		CParMvmA2Index *pobjloc = (CParMvmA2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *palloc = pobjloc->pal;
		CGSMatrixR *pauloc = pobjloc->pau;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixR *pa2loc = pobjloc->pa2;
		CGSMatrixR gmtrdummy;
		gmtrdummy.MvmA2Index (*ptreeloc, *pmvmloc, *palloc, *pauloc, _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmA2Index2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
	static void MvmAt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParMvmA2Index *pobjloc = (CParMvmA2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *palloc = pobjloc->pal;
		CGSMatrixR *pauloc = pobjloc->pau;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixR *pa2loc = pobjloc->pa2;
		CGSMatrixR gmtrdummy;
		gmtrdummy.MvmA2Index (*ptreeloc, *pmvmloc, *pauloc, *palloc, _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmAh2Index2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParMvmAC
#define __ParMvmAC

// Preliminary declarations

class CTree;
class CMvmC;

class CParMvmAC
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CGSMatrixCS *pa; // pa contains pointer to A matrix
public:
// Functions
// Constructors and destructor
	CParMvmAC (CTree *_ptree, CGSMatrixCS *_pa) { // Memory allocation zero data constructor
		ptree = _ptree;
		pa = _pa;
	};
	CParMvmAC (const CParMvmAC &_mvm) {throw " Copy constructor for CParMvmAC class called";}; // Copy constructor
	~CParMvmAC () {}; // Destructor
// Operator functions
	CParMvmAC &operator= (const CParMvmAC &_mvm) {throw " Equality operator for CParMvmAC class called";}; // Equality operator
// Solve function
	static void MvmAC (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParMvmAC *pobjloc = (CParMvmAC *)_objlu;
		CGSMatrixCS *paloc = pobjloc->pa;
		paloc->MvmARect (_x, _px);
	};
	static void MvmACH (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParMvmAC *pobjloc = (CParMvmAC *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CGSMatrixCS *paloc = pobjloc->pa;
		paloc->MvmAhRect (*ptreeloc, _x, _px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParMvmALsqC
#define __ParMvmALsqC

// Preliminary declarations

class CTree;
class CMvmR;

class CParMvmALsqDenseC
{
public:
// Data
	int m;        // the number of columns
	int n;        // the number of columns (without regularization)
	double alpha; // Tikhonov regularization parameter
	dcmplx *pa;   // pa contains pointer to the complex rectangular dense coefficient matrix
public:
// Functions
// Constructors and destructor
	CParMvmALsqDenseC (int _m, int _n, double _alpha, dcmplx *_pa) { // Memory allocation zero data constructor
		m = _m;
		n = _n;
		alpha = _alpha;
		pa = _pa;
	};
	CParMvmALsqDenseC (const CParMvmALsqDenseC &_mvm) {throw " Copy constructor for CParMvmALsqDenseC class called";}; // Copy constructor
	~CParMvmALsqDenseC () {}; // Destructor
// Operator functions
	CParMvmALsqDenseC &operator= (const CParMvmALsqDenseC &_mvm) {throw " Equality operator for CParMvmALsqDenseC class called";}; // Equality operator
// Solve function
	static void MvmAC (void *_objlu, CSVectorC &_x, CSVectorC &_ax) {
		dcmplx czero (0.0e0,0.0e0);
		CParMvmALsqDenseC *pobjloc = (CParMvmALsqDenseC *)_objlu;
		int mloc = pobjloc->m;
		int nloc = pobjloc->n;
		double alphaloc  = pobjloc->alpha;
		dcmplx *paloc  = pobjloc->pa;
		_ax.SetSVect (czero);
		dcmplx *pxloc = _x.GetVect ();
		dcmplx *paxloc = _ax.GetVect ();
		CSVectorC::MvmBlock (mloc, nloc, paloc, pxloc, paxloc);
		int i;
		for (i=0;i<nloc;i++) paxloc[mloc+i] += pxloc[i] * alphaloc;
	};
	static void MvmAHC (void *_objlu, CSVectorC &_x, CSVectorC &_ax) {
		dcmplx czero (0.0e0,0.0e0);
		CParMvmALsqDenseC *pobjloc = (CParMvmALsqDenseC *)_objlu;
		int mloc = pobjloc->m;
		int nloc = pobjloc->n;
		double alphaloc  = pobjloc->alpha;
		dcmplx *paloc  = pobjloc->pa;
		_ax.SetSVect (czero);
		dcmplx *pxloc = _x.GetVect ();
		dcmplx *paxloc = _ax.GetVect ();
		CSVectorC::MvmBlockH (mloc, nloc, paloc, pxloc, paxloc);
		int i;
		for (i=0;i<nloc;i++) paxloc[i] += pxloc[mloc+i] * alphaloc;
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSolveLU
#define __ParSolveLU

// Preliminary declarations

class CTree;
class CMvmR;
class CGVector;

class CParSolveLU
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CMvmR *pmvm;     // pmvm contains pointer to mvm structure
	CGSMatrixRS *pl; // pl contains pointer to L part
	CGSMatrixRS *pu; // pu contains pointer to U part
public:
// Functions
// Constructors and destructor
	CParSolveLU (CTree *_ptree, CMvmR *_pmvm, CGSMatrixRS *_pl, CGSMatrixRS *_pu) { // Memory allocation zero data constructor
		ptree = _ptree;
		pmvm = _pmvm;
		pl = _pl;
		pu = _pu;
	};
	CParSolveLU (const CParSolveLU &_mvm) {throw " Copy constructor for CParSolveLU class called";}; // Copy constructor
	~CParSolveLU () {}; // Destructor
// Operator functions
	CParSolveLU &operator= (const CParSolveLU &_mvm) {throw " Equality operator for CParSolveLU class called";}; // Equality operator
// Solve function
	static void SolveLU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU *pobjloc = (CParSolveLU *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *plloc = pobjloc->pl;
		CGSMatrixRS *puloc = pobjloc->pu;
		CSVector p;
		p = _x;
		plloc->SolveL (*ptreeloc,*pmvmloc,_x,p);
		puloc->SolveU (*ptreeloc,*pmvmloc,p,_px);
	};
	static void SolveL (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU *pobjloc = (CParSolveLU *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *plloc = pobjloc->pl;
		plloc->SolveL (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveLt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU *pobjloc = (CParSolveLU *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *plloc = pobjloc->pl;
		plloc->SolveU (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU *pobjloc = (CParSolveLU *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *puloc = pobjloc->pu;
		puloc->SolveU (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveUt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU *pobjloc = (CParSolveLU *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixRS *puloc = pobjloc->pu;
		puloc->SolveL (*ptreeloc,*pmvmloc,_x,_px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSolveLU2Ind
#define __ParSolveLU2Ind

// Preliminary declarations

class CTree;
class CMvmR;
class CSVector;

class CParSolveLU2Index
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CMvmR *pmvm;     // pmvm contains pointer to mvm structure
	CGSMatrixR *pl;  // pl contains pointer to L part
	CGSMatrixR *pu;  // pu contains pointer to U part
public:
// Functions
// Constructors and destructor
	CParSolveLU2Index (CTree *_ptree, CMvmR *_pmvm, CGSMatrixR *_pl, CGSMatrixR *_pu) { // Memory allocation zero data constructor
		ptree = _ptree;
		pmvm = _pmvm;
		pl = _pl;
		pu = _pu;
	};
	CParSolveLU2Index (const CParSolveLU2Index &_mvm) {throw " Copy constructor for CParSolveLU2Ind class called";}; // Copy constructor
	~CParSolveLU2Index () {}; // Destructor
// Operator functions
	CParSolveLU2Index &operator= (const CParSolveLU2Index &_mvm) {throw " Equality operator for CParSolveLU2Index class called";}; // Equality operator
// Solve function
	static void SolveLU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU2Index *pobjloc = (CParSolveLU2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *plloc = pobjloc->pl;
		CGSMatrixR *puloc = pobjloc->pu;
		CSVector p;
		p = _x;
		plloc->SolveL2Index (*ptreeloc,*pmvmloc,_x,p);
		puloc->SolveU2Index (*ptreeloc,*pmvmloc,p,_px);
	};
	static void SolveL (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU2Index *pobjloc = (CParSolveLU2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *plloc = pobjloc->pl;
		plloc->SolveL2Index (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveLt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU2Index *pobjloc = (CParSolveLU2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *plloc = pobjloc->pl;
		plloc->SolveU2Index (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU2Index *pobjloc = (CParSolveLU2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *puloc = pobjloc->pu;
		puloc->SolveU2Index (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveUt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSolveLU2Index *pobjloc = (CParSolveLU2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmR *pmvmloc = pobjloc->pmvm;
		CGSMatrixR *puloc = pobjloc->pu;
		puloc->SolveL2Index (*ptreeloc,*pmvmloc,_x,_px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSchurMvmSlv2Index
#define __ParSchurMvmSlv2Index

// Preliminary declarations

class CTree;
class CMvmR;
class CMvmSchurR;
class CSVector;

class CParSchurMvmSlv2Index
{
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CMvmSchurR *pmvm;     // pmvm contains pointer to mvm structure
	bool second;      // second specifies second order data
	CMvmR *pmvm2;     // pmvm2 contains pointer to mvm2 structure
	CGSMatrixR *pa2;  // pa contains pointer to A2 part
public:
// Functions
// Constructors and destructor
	CParSchurMvmSlv2Index (CTree *_ptree, CMvmSchurR *_pmvm, bool _second, CMvmR *_pmvm2, CGSMatrixR *_pa2) { // Memory allocation zero data constructor
		ptree = _ptree;
		pmvm = _pmvm;
		second = _second;
		pmvm2 = _pmvm2;
		pa2 = _pa2;
	};
	CParSchurMvmSlv2Index (const CParSchurMvmSlv2Index &_mvm) {throw " Copy constructor for CParSchurMvmSlv2Index class called";}; // Copy constructor
	~CParSchurMvmSlv2Index () {}; // Destructor
// Operator functions
	CParSchurMvmSlv2Index &operator= (const CParSchurMvmSlv2Index &_mvm) {throw " Equality operator for CParSchurMvmSlv2Index class called";}; // Equality operator
// Mvm functions
	static void MvmA (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixR *pa2loc = pobjloc->pa2;
		pmvmloc->AbstractBlockMvmA ('N', _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmA2Index2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
	static void MvmAt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		CMvmR *pmvm2loc = pobjloc->pmvm2;
		CGSMatrixR *pa2loc = pobjloc->pa2;
		pmvmloc->AbstractBlockMvmA ('T', _x, _px);
		if (pobjloc->second) {
			pa2loc->MvmAh2Index2 (*ptreeloc, *pmvm2loc, _x, _px);
		};
	};
// Solve functions
	static void SolveLU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		CSVector p;
		p = _x;
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, p);
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												p, _px);
	};
	static void SolveL (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveAL (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		pmvmloc->AbstractBlockSolveL ('A', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveLt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		pmvmloc->AbstractBlockSolveU ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveUt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlv2Index *pobjloc = (CParSchurMvmSlv2Index *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		pmvmloc->AbstractBlockSolveL ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSchurMvmSlvRS
#define __ParSchurMvmSlvRS

// Preliminary declarations

class CTree;
class CMvmR;
class CMvmSchurR;
class CSVector;

class CParSchurMvmSlvRS
{
// Pointers to the existing data
	CTree *ptree;         // ptree contains pointer to the tree
	CGSMatrixRS *pgmtral; // pgmtral contains pointer to gmtral part
	CGSMatrixRS *pgmtrau; // pgmtrau contains pointer to gmtrau part
	CGSMatrixRS *pgmtrl;  // pgmtral contains pointer to gmtrl  part
	CGSMatrixRS *pgmtru;  // pgmtrau contains pointer to gmtru  part
	CMvmSchurR *pmvm;     // pmvm contains pointer to mvm structure
// Local data
	char optype; // Type of the operation to be performed
	int *ibssup; // Bases of supernodes in vector structures
public:
// Functions
// Constructors and destructor
	CParSchurMvmSlvRS (CTree *_ptree, // Memory allocation zero data constructor
								CGSMatrixRS *_pgmtral, CGSMatrixRS *_pgmtrau, CGSMatrixRS *_pgmtrl, CGSMatrixRS *_pgmtru, 
								CMvmSchurR *_pmvm) {
		ptree = _ptree;
		pgmtral = _pgmtral;
		pgmtrau = _pgmtrau;
		pgmtrl = _pgmtrl;
		pgmtru = _pgmtru;
		pmvm = _pmvm;
		optype = ' ';
		ibssup = 0;
	};
	CParSchurMvmSlvRS (const CParSchurMvmSlvRS &_mvm) {throw " Copy constructor for CParSchurMvmSlvRS class called";}; // Copy constructor
	~CParSchurMvmSlvRS () {if (ibssup != 0) delete [] ibssup;}; // Destructor
// Operator functions
	CParSchurMvmSlvRS &operator= (const CParSchurMvmSlvRS &_mvm) {throw " Equality operator for CParSchurMvmSlvRS class called";}; // Equality operator
// Get/set functions
	char GetOptype () {return optype;}; // Get optype
	int * GetIbssup () {return ibssup;}; // Get ibssup
	void SetOptype (char _optype) {optype = _optype;}; // Set optype
// Global mvm functions
	static void MvmA (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('N');
		pmvmloc->AbstractBlockMvmA ('N', _x, _px);
	};
	static void MvmAt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('T');
		pmvmloc->AbstractBlockMvmA ('T', _x, _px);
	};
// Global solve functions
	static void SolveLU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		CSVector p;
		p = _x;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, p);
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												p, _px);
	};
	static void SolveL (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveLt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveU ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveU (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveUt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvRS *pobjloc = (CParSchurMvmSlvRS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurR *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveL ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
// Local mvm functions
	static void AbsMvmBlockRowAu (void *_obj, int _iblk, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax);
	static void AbsMvmBlockColumnAl (void *_obj, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax);
// Local solve functions
	static void AbsSolveBlockColumnL (void *_obj, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVector &_x, 
								int *_bslx, CSVector &_lx);
	static void AbsSolveBlockRowU (void *_obj, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVector &_x, 
								int *_bsux, CSVector &_ux);
// Init functions
	void CreateIbssup (); // Create ibssup data for local operations
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSchurMvmSlvCS
#define __ParSchurMvmSlvCS

// Preliminary declarations

class CTree;
class CMvmC;
class CMvmSchurC;
class CSVectorC;

class CParSchurMvmSlvCS
{
// Pointers to the existing data
	CTree *ptree;         // ptree contains pointer to the tree
	CGSMatrixCS *pgmtral; // pgmtral contains pointer to gmtral part
	CGSMatrixCS *pgmtrau; // pgmtrau contains pointer to gmtrau part
	CGSMatrixCS *pgmtrl;  // pgmtral contains pointer to gmtrl  part
	CGSMatrixCS *pgmtru;  // pgmtrau contains pointer to gmtru  part
	CMvmSchurC *pmvm;     // pmvm contains pointer to mvm structure
// Local data
	char optype; // Type of the operation to be performed
	int *ibssup; // Bases of supernodes in vector structures
public:
// Functions
// Constructors and destructor
	CParSchurMvmSlvCS (CTree *_ptree, // Memory allocation zero data constructor
								CGSMatrixCS *_pgmtral, CGSMatrixCS *_pgmtrau, CGSMatrixCS *_pgmtrl, CGSMatrixCS *_pgmtru, 
								CMvmSchurC *_pmvm) {
		ptree = _ptree;
		pgmtral = _pgmtral;
		pgmtrau = _pgmtrau;
		pgmtrl = _pgmtrl;
		pgmtru = _pgmtru;
		pmvm = _pmvm;
		optype = ' ';
		ibssup = 0;
	};
	CParSchurMvmSlvCS (const CParSchurMvmSlvCS &_mvm) {throw " Copy constructor for CParSchurMvmSlvCS class called";}; // Copy constructor
	~CParSchurMvmSlvCS () {if (ibssup != 0) delete [] ibssup;}; // Destructor
// Operator functions
	CParSchurMvmSlvCS &operator= (const CParSchurMvmSlvCS &_mvm) {throw " Equality operator for CParSchurMvmSlvCS class called";}; // Equality operator
// Get/set functions
	char GetOptype () {return optype;}; // Get optype
	int * GetIbssup () {return ibssup;}; // Get ibssup
	void SetOptype (char _optype) {optype = _optype;}; // Set optype
// Global mvm functions
	static void MvmA (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('N');
		pmvmloc->AbstractBlockMvmA ('N', _x, _px);
	};
	static void MvmAt (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('T');
		pmvmloc->AbstractBlockMvmA ('T', _x, _px);
	};
// Global solve functions
	static void SolveLU (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		CSVectorC p;
		p = _x;
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, p);
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												p, _px);
	};
	static void SolveL (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveL ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveLt (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('L');
		pmvmloc->AbstractBlockSolveU ('L', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveU (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveU ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
	static void SolveUt (void *_objlu, CSVector &_x, CSVector &_px) {
		CParSchurMvmSlvCS *pobjloc = (CParSchurMvmSlvCS *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmSchurC *pmvmloc = pobjloc->pmvm;
		int nblksloc = pmvmloc->GetNblks ();
//		int *pblksloc = pmvmloc->GetPblks ();
		int *pblk2cpuloc = pmvmloc->GetPblk2cpu ();
		CSMatrix *pablkstrloc = pmvmloc->GetPablkstr ();
		if (pobjloc->GetIbssup ()==0) pobjloc->CreateIbssup ();
		pobjloc->SetOptype ('U');
		pmvmloc->AbstractBlockSolveL ('U', 
												*ptreeloc, *pablkstrloc, nblksloc, pblk2cpuloc,
												_x, _px);
	};
// Local mvm functions
	static void AbsMvmBlockRowAu (void *_obj, void *_mtr, int _iblk, // Multiply by the current block row of Au
								int *_bsx, const CSVectorC &_x, 
								int *_bsax, CSVectorC &_ax);
	static void AbsMvmBlockColumnAl (void *_obj, void *_mtr, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVectorC &_x, 
								int *_bsax, CSVectorC &_ax);
// Local solve functions
	static void AbsSolveBlockColumnL (void *_obj, void *_mtr, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVectorC &_x, 
								int *_bslx, CSVectorC &_lx);
	static void AbsSolveBlockRowU (void *_obj, void *_mtr, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVectorC &_x, 
								int *_bsux, CSVectorC &_ux);
// Init functions
	void CreateIbssup (); // Create ibssup data for local operations
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSolveUhUC
#define __ParSolveUhUC

// Preliminary declarations

class CTree;
class CMvmR;
class CGVector;

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CParSolveUhUC
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	CMvmC *pmvm;     // pmvm contains pointer to mvm structure
	CGSMatrixCS *pu; // pu contains pointer to U part
public:
// Functions
// Constructors and destructor
	CParSolveUhUC (CTree *_ptree, CMvmC *_pmvm, CGSMatrixCS *_pu) { // Memory allocation zero data constructor
		ptree = _ptree;
		pmvm = _pmvm;
		pu = _pu;
	};
	CParSolveUhUC (const CParSolveUhUC &_mvm) {throw " Copy constructor for CParSolveUhUC class called";}; // Copy constructor
	~CParSolveUhUC () {}; // Destructor
// Operator functions
	CParSolveUhUC &operator= (const CParSolveUhUC &_mvm) {throw " Equality operator for CParSolveUhUC class called";}; // Equality operator
// Solve functions
	static void SolveUh (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSolveUhUC *pobjloc = (CParSolveUhUC *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmC *pmvmloc = pobjloc->pmvm;
		CGSMatrixCS *puloc = pobjloc->pu;
		puloc->SolveLConj (*ptreeloc,*pmvmloc,_x,_px);
	};
	static void SolveU (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSolveUhUC *pobjloc = (CParSolveUhUC *)_objlu;
		CTree *ptreeloc = pobjloc->ptree;
		CMvmC *pmvmloc = pobjloc->pmvm;
		CGSMatrixCS *puloc = pobjloc->pu;
		puloc->SolveU (*ptreeloc,*pmvmloc,_x,_px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __ParSolveIchC
#define __ParSolveIchC

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CParSolveIchC
{
public:
// Data
	CTree * ptree;   // ptree contains pointer to the tree
	int *porder;     // porder contains pointer to the order array
	CGSMatrixCS *pu; // pu contains pointer to U part
public:
// Functions
// Constructors and destructor
	CParSolveIchC (CTree *_ptree, int *_porder, CGSMatrixCS *_pu) { // Memory allocation zero data constructor
		ptree = _ptree;
		porder = _porder;
		pu = _pu;
	};
	CParSolveIchC (const CParSolveIchC &_mvm) {throw " Copy constructor for CParSolveIchC class called";}; // Copy constructor
	~CParSolveIchC () {}; // Destructor
// Operator functions
	CParSolveIchC &operator= (const CParSolveIchC &_mvm) {throw " Equality operator for CParSolveIchC class called";}; // Equality operator
// Solve function
	static void SolveU (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSolveIchC *pobjloc = (CParSolveIchC *)_objlu;
//		CTree *ptreeloc = pobjloc->ptree;
		int *porderloc = pobjloc->porder;
		CGSMatrixCS *puloc = pobjloc->pu;
		int nloc = puloc->GetN ();
		CSVectorC temp (nloc);
		puloc->SolveU (_x, temp);
		dcmplx *ptemp = temp.GetVect ();
		dcmplx *ppx = _px.GetVect ();
		CSVectorC::OrderVector (-1, nloc, porderloc, ptemp, ppx);
	};
	static void SolveUH (void *_objlu, CSVectorC &_x, CSVectorC &_px) {
		CParSolveIchC *pobjloc = (CParSolveIchC *)_objlu;
//		CTree *ptreeloc = pobjloc->ptree;
		int *porderloc = pobjloc->porder;
		CGSMatrixCS *puloc = pobjloc->pu;
		int nloc = puloc->GetN ();
		CSVectorC temp (nloc);
		dcmplx *ptemp = temp.GetVect ();
		dcmplx *px = _x.GetVect ();
		CSVectorC::OrderVector (1, nloc, porderloc, px, ptemp);
		puloc->SolveLConj (temp, _px);
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __AbsSymmSolve
#define __AbsSymmSolve

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CAbsSymmSolve
{
public:
// Data
	CSMatrixRS *pmatrix; // pmatrix contains pointer to the matrix
	CSVector *prhs;      // prhs contains pointer to rhs
	CSVector *px;        // px contains pointer to the initial guess
	CSVector *psol;      // psol contains pointer to the solution
public:
// Functions
// Constructors and destructor
	CAbsSymmSolve (CSMatrixRS *_pmatrix, CSVector *_prhs, CSVector *_px, CSVector *_psol) { // Memory allocation zero data constructor
		pmatrix = _pmatrix;
		prhs = _prhs;
		px = _px;
		psol = _psol;
	};
	CAbsSymmSolve (const CAbsSymmSolve &_mvm) {throw " Copy constructor for CAbsSymmSolve class called";}; // Copy constructor
	~CAbsSymmSolve () {}; // Destructor
// Operator functions
	CAbsSymmSolve &operator= (const CAbsSymmSolve &_mvm) {throw " Equality operator for CAbsSymmSolve class called";}; // Equality operator
// Solve function
	static void Nsup (void *_pobj, int &_nsup) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		_nsup = pmatrixloc->GetNsupr ();
	};
	static void Sprnds (void *_pobj, int *_sprnds) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		int nsuploc = pmatrixloc->GetNsupr ();
		int *psprnds = pmatrixloc->GetSprndr ();
		int i;
		for (i=0;i<=nsuploc;i++) _sprnds[i] = psprnds[i];
	};
	static void IaJa (void *_pobj, int *&_ia, int *&_ja) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		int nsuploc = pmatrixloc->GetNsupr ();
		int nzjaloc = pmatrixloc->GetNzja ();
		int *pia = pmatrixloc->GetIa ();
		int *pja = pmatrixloc->GetJa ();
		_ia = new int [nsuploc+1];
		_ja = new int [nzjaloc];
		int i;
		for (i=0;i<=nsuploc;i++) _ia[i] = pia[i];
		for (i=0;i<nzjaloc;i++) _ja[i] = pja[i];
	};
	static void Blockrow (void *_pobj, int _isup, double *_arowdata) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		int *psprnds = pmatrixloc->GetSprndr ();
		int *pia = pmatrixloc->GetIa ();
		int *pja = pmatrixloc->GetJa ();
		int *pbsa = pmatrixloc->GetBsa ();
		double *pa = pmatrixloc->GetA ();
		int blai = psprnds[_isup+1]-psprnds[_isup];
		int j, jj, ibs, blaj, blaij, k;
		int nz = 0;
		for (j=pia[_isup];j<pia[_isup+1];j++) {
			jj = pja[j];
			blaj = psprnds[jj+1]-psprnds[jj];
			ibs = pbsa[j];
			blaij = blai*blaj;
			for (k=0;k<blaij;k++) _arowdata[nz+k] = pa[ibs+k];
			nz += blaij;
		};
	};
	static void Rhs (void *_pobj, int _isup, double *_rhs) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		CSVector *prhs = pobjloc->prhs;
		int *psprnds = pmatrixloc->GetSprndr ();
		double *pv = prhs->GetVect();
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) _rhs[i-psprnds[_isup]] = pv[i];
	};
	static void Iniguess (void *_pobj, int _isup, double *_x) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		CSVector *px = pobjloc->px;
		int *psprnds = pmatrixloc->GetSprndr ();
		double *pv = px->GetVect();
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) _x[i-psprnds[_isup]] = pv[i];
	};
	static void Storesol (void *_pobj, int _isup, double *_sol) {
		CAbsSymmSolve *pobjloc = (CAbsSymmSolve *)_pobj;
		CSMatrixRS *pmatrixloc = pobjloc->pmatrix;
		CSVector *psol = pobjloc->psol;
		int *psprnds = pmatrixloc->GetSprndr ();
		double *pv = psol->GetVect();
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) pv[i] = _sol[i-psprnds[_isup]];
	};
// Input/Output functions
// Friend classes
};

#endif

#ifndef __AbsCSolve
#define __AbsCSolve

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CAbsCSolve
{
public:
// Data
	CSMatrixCS *pmatrix; // pmatrix contains pointer to the matrix
	CSVectorC *prhs;      // prhs contains pointer to rhs
	CSVectorC *px;        // px contains pointer to the initial guess
	CSVectorC *psol;      // psol contains pointer to the solution
public:
// Functions
// Constructors and destructor
	CAbsCSolve (CSMatrixCS *_pmatrix, CSVectorC *_prhs, CSVectorC *_px, CSVectorC *_psol) { // Memory allocation zero data constructor
		pmatrix = _pmatrix;
		prhs = _prhs;
		px = _px;
		psol = _psol;
	};
	CAbsCSolve (const CAbsCSolve &_mvm) {throw " Copy constructor for CAbsCSolve class called";}; // Copy constructor
	~CAbsCSolve () {}; // Destructor
// Operator functions
	CAbsCSolve &operator= (const CAbsCSolve &_mvm) {throw " Equality operator for CAbsCSolve class called";}; // Equality operator
// Solve function
	static void Nsup (void *_pobj, int &_nsup) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		_nsup = pmatrixloc->GetNsupr ();
	};
	static void Sprnds (void *_pobj, int *_sprnds) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		int nsuploc = pmatrixloc->GetNsupr ();
		int *psprnds = pmatrixloc->GetSprndr ();
		int i;
		for (i=0;i<=nsuploc;i++) _sprnds[i] = psprnds[i];
	};
	static void IaJa (void *_pobj, int *&_ia, int *&_ja) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		int nsuploc = pmatrixloc->GetNsupr ();
		int nzjaloc = pmatrixloc->GetNzja ();
		int *pia = pmatrixloc->GetIa ();
		int *pja = pmatrixloc->GetJa ();
		_ia = new int [nsuploc+1];
		_ja = new int [nzjaloc];
		int i;
		for (i=0;i<=nsuploc;i++) _ia[i] = pia[i];
		for (i=0;i<nzjaloc;i++) _ja[i] = pja[i];
	};
	static void Blockrow (void *_pobj, int _isup, void *_arowdata) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		int *psprnds = pmatrixloc->GetSprndr ();
		int *pia = pmatrixloc->GetIa ();
		int *pja = pmatrixloc->GetJa ();
		int *pbsa = pmatrixloc->GetBsa ();
		dcmplx *pa = pmatrixloc->GetA ();
		dcmplx *arowdata = (dcmplx *)_arowdata;
		int blai = psprnds[_isup+1]-psprnds[_isup];
		int j, jj, ibs, blaj, blaij, k;
		int nz = 0;
		for (j=pia[_isup];j<pia[_isup+1];j++) {
			jj = pja[j];
			blaj = psprnds[jj+1]-psprnds[jj];
			ibs = pbsa[j];
			blaij = blai*blaj;
			for (k=0;k<blaij;k++) arowdata[nz+k] = pa[ibs+k];
			nz += blaij;
		};
	};
	static void Rhs (void *_pobj, int _isup, void *_rhs) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		CSVectorC *prhs = pobjloc->prhs;
		int *psprnds = pmatrixloc->GetSprndr ();
		dcmplx *pv = prhs->GetVect();
		dcmplx *rhsloc = (dcmplx *)_rhs;
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) rhsloc[i-psprnds[_isup]] = pv[i];
	};
	static void Iniguess (void *_pobj, int _isup, void *_x) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		CSVectorC *px = pobjloc->px;
		int *psprnds = pmatrixloc->GetSprndr ();
		dcmplx *pv = px->GetVect();
		dcmplx *xloc = (dcmplx *)_x;
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) xloc[i-psprnds[_isup]] = pv[i];
	};
	static void Storesol (void *_pobj, int _isup, void *_sol) {
		CAbsCSolve *pobjloc = (CAbsCSolve *)_pobj;
		CSMatrixCS *pmatrixloc = pobjloc->pmatrix;
		CSVectorC *psol = pobjloc->psol;
		int *psprnds = pmatrixloc->GetSprndr ();
		dcmplx *pv = psol->GetVect();
		dcmplx *sol = (dcmplx *)_sol;
		int i;
		for (i=psprnds[_isup];i<psprnds[_isup+1];i++) pv[i] = sol[i-psprnds[_isup]];
	};
// Input/Output functions
// Friend classes
};

#endif
