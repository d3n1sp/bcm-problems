//------------------------------------------------------------------------------------------------
// File: smatrix.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"
#include "svector.h"
#include "Decomp.h"

typedef void (*CNSUP)      (void *_obj, int &_nsup);                   // Define interface of the slae generation, number of supernodes
typedef void (*CSPRNDS)    (void *_obj, int *_sprnds);                 // Define interface of the slae generation, supernodes description
typedef void (*CIAJA)      (void *_obj, int *&_ia, int *&_ja);         // Define interface of the slae generation, supernodes sparsity description
typedef void (*CBLOCKROW)  (void *_obj, int _isup, double *_arowdata); // Define interface of the slae generation, current supernode row data
typedef void (*CRHS)       (void *_obj, int _isup, double *_rhs);      // Define interface of the slae generation, current part of right hand side data
typedef void (*CINIGUESS)  (void *_obj, int _isup, double *_xini);     // Define interface of the slae generation, current part of initial guess
typedef void (*CSTORESOL)  (void *_obj, int _isup, double *_sol);      // Define interface of the slae generation, store current part of the solution
typedef void (*CBLOCKROWC) (void *_obj, int _isup, void   *_arowdata); // Define interface of the slae generation, current supernode row data
typedef void (*CRHSC)      (void *_obj, int _isup, void   *_rhs);      // Define interface of the slae generation, current part of right hand side data
typedef void (*CINIGUESSC) (void *_obj, int _isup, void   *_xini);     // Define interface of the slae generation, current part of initial guess
typedef void (*CSTORESOLC) (void *_obj, int _isup, void   *_sol);      // Define interface of the slae generation, store current part of the solution

void AbstractSymmSolver (void *_pgener, // Solve real block symmetric positive system of equations in abstract interface
									std::ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol);
void AbstractSolver (void *_pgener, // Solve real block system of equations in abstract interface
									std::ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol);
void AbstractParSolver (void *_pgener, // Solve real block system of equations in parallel mode in abstract interface
									void *_comm, std::ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROW _fblockrow, 
									CRHS _frhs, CINIGUESS _finiguess, CSTORESOL _fstoresol);
void AbstractParSolver (void *_pgener, // Solve real block system of equations in parallel mode in abstract interface
									void *_comm, std::ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROWC _fblockrow, 
									CRHSC _frhs, CINIGUESSC _finiguess, CSTORESOLC _fstoresol);
void AbstractSolver (void *_pgener, // Solve complex block system of equations in abstract interface
									std::ofstream &_fout, CSlvParam &_param,
									CNSUP _fnsup, CSPRNDS _fsprnds, CIAJA _fiaja, CBLOCKROWC _fblockrow, 
									CRHSC _frhs, CINIGUESSC _finiguess, CSTORESOLC _fstoresol);

class CSlvParam;
class CTree;
class CSVector;
class CSVectorC;

// SMatrix.h: Description of the sparsity structure of the matrix stored in super sparse format
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrix
#define __SMatrix

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CSMatrix
{
protected:
// Data
	int m;     // m is the total number of rows    in the matrix
	int n;     // n is the total number of columns in the matrix
	int nsupr; // nsupr is the number of supernode rows    in the matrix
	int nsupc; // nsupc is the number of supernode columns in the matrix
	int nlist; // nlist is the number of nonzero rows or columns in the matrix
	int nlist2;// nlist2 is the number of secondary nonzero rows or columns in the matrix
	int nlist3;// nlist3 is the number of secondary nonzero rows or columns in the matrix
	int nzja;  // nzja is the number of nonzero elements in the matrix
	int nzja2; // nzja2 is the number of secondary nonzero elements in the matrix
	int nzja3; // nzja3 is the number of thirdly nonzero elements in the matrix
	int *list; // list[nlist] array describes the list of nonzero supernode rows / colunms
	int *list2;// list2[nlist2] array describes the list of secondary nonzero supernode rows / colunms
	int *list3;// list3[nlist3] array describes the list of secondary nonzero supernode rows / colunms
	int *ia;   // ia[nlist+1] array describes the number of supernode elements in each nonzero supernode row / colunm
	int *ja;   // ja[nzja] array describes the consecutive supernode column / row numbers of each supernode row /column
	int *ja2;  // ja2[nzja2] array describes the secondary consecutive supernode column / row numbers of each supernode row /column
	int *ja3;  // ja3[nzja3] array describes the secondary consecutive supernode column / row numbers of each supernode row /column
protected:
// Factorization support functions
// Ich2:
	void Ich2AllocInt (int &_icycle, int *&_imask, int *&_lstloc, // Allocate and init integer data for Ich2/Ilu2
						int *&_iv, int *&_madj, int *&_ibegm, int *&_lstupd, int *&_indlst) const; 
	void Ich2UpdateTransp (int _irow, int *_ibegm, int *_madj, // Update arrays that support transposed structures search
							int *_iv, int *_ig, int *_jg) const;
	void Ich2UpdateTranspDynamic (int _irow, int *_ibegm, int *_madj, // Update arrays that support transposed structures search
									int *_iv, int *_ig, int **_jg) const;
public:
// Functions
// Constructors and destructor
	CSMatrix (); // Memory allocation zero data constructor
	CSMatrix (int _nlist, int _nzja); // Memory allocation zero data constructor
	CSMatrix (int _nlist, int _nlist2, int _nzja, int _nzja2); // Memory allocation zero data constructor
	CSMatrix (const CSMatrix &_aa); // Copy constructor
	virtual ~CSMatrix (); // Destructor
// Get/set functions
	int GetM      () const {return m;};     // Get m
	int GetN      () const {return n;};     // Get n
	int GetNsupr  () const {return nsupr;}; // Get nsupr
	int GetNsupc  () const {return nsupc;}; // Get nsupc
	int GetNlist  () const {return nlist;}; // Get nlist
	int GetNlist2 () const {return nlist2;};// Get nlist2
	int GetNzja   () const {return nzja;};  // Get nzja
	int GetNzja2  () const {return nzja2;}; // Get nzja2
	int GetNzja3  () const {return nzja3;}; // Get nzja3
	int GetNza    () const {return nzja;};  // Get nza
	int * GetList () const {return list;};  // Get list
	int * GetList2() const {return list2;}; // Get list2
	int * GetList3() const {return list3;}; // Get list3
	int * GetIa   () const {return ia;};    // Get ia
	int * GetJa   () const {return ja;};    // Get ja
	int * GetJa2  () const {return ja2;};    // Get ja2
	int * GetJa3  () const {return ja3;};    // Get ja3
	void GetSparsity (int *&_list, int *&_ia, int *&_ja) const { // Get the sparsity arrays
		_list = list;
		_ia = ia;
		_ja = ja;
	};
	void SetM     (int _m)     {m=_m;};           // Set m
	void SetN     (int _n)     {n=_n;};           // Set n
	void SetNsupr (int _nsupr) {nsupr = _nsupr;}; // Set nsupr
	void SetNsupc (int _nsupc) {nsupc = _nsupc;}; // Set nsupc
	void SetNsupcBase (int _nsupc) {nsupc = _nsupc;}; // Set nsupc
	void SetNlist (int _nlist) {nlist = _nlist;}; // Set nlist
	void SetNlist2(int _nlist2) {nlist2 = _nlist2;}; // Set nlist2
	void SetNlist3(int _nlist3) {nlist3 = _nlist3;}; // Set nlist3
	void SetNzja  (int _nzja ) {nzja  = _nzja; }; // Set nzja
	void SetNzja2 (int _nzja2) {nzja2 = _nzja2; }; // Set nzja2
	void SetNzja3 (int _nzja3) {nzja3 = _nzja3; }; // Set nzja3
	void SetList3 (int *_list3) {list3 = _list3;}; // Set list3
	void SetIa    (int *_ia) {ia = _ia;};         // Set ia
	void SetJa    (int *_ja) {ja = _ja;};         // Set ja
	void SetJa2   (int *_ja2) {ja2 = _ja2;};      // Set ja2
	void SetJa3   (int *_ja3) {ja3 = _ja3;};      // Set ja3
// Operator functions
	CSMatrix &operator= (const CSMatrix &_aa); // Equality operator
	CSMatrix  operator+ (const CSMatrix &_mtr2) const; // Add two matrices
	void AllocateJa3 (int _nzja3) {delete [] ja3; ja3 = new int [_nzja3];}; // Allocate ja3
// Transformation functions
	void ComputeBlocks2Index (int &_nblks, int *&_blks) const; // Compute minimal blocks partitioning for 2 index matrix
	static CSMatrix Sparsity3D (int _nx, int _ny, int _nz); // Compute sparsity of the 3D matrix
	static CSMatrix Sparsity3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz); // Compute sparsity of the 3D matrix
	static CSMatrix TridiagonalMtr (int _n); // Compute tridiagonal matrix 
	static CSMatrix DiagonalMtr (int _n); // Compute diagonal matrix 
	void SuperSparsify (); // Super sparsify the matrix
	void SuperSparsify2Index (); // Super sparsify the matrix
	CSMatrix TranspMtr     () const; // Compute transposed matrix 
	CSMatrix TranspMtrRect () const; // Compute transposed to the rectangular matrix stored by columns
	CSMatrix TranspMtrRectRows () const; // Compute transposed to the rectangular matrix stored by rows
	CSMatrix TranspMtrRect (int *_indarr) const; // Compute transposed to the rectangular matrix stored by columns
	CSMatrix TranspMtrMask () const; // Compute transposed matrix via the mask (list is incomplete)
	CSMatrix TranspMtrMask (int &_icycle, int *_imask) const; // Compute transposed matrix via the mask (list is incomplete)
	CSMatrix TranspMtrMask2Index () const; // Compute transposed 2 index matrix via the mask (list is incomplete)
	CSMatrix SymmMtr     () const; // Compute symmetrized matrix 
	CSMatrix SymmMtr     (CMPIComm &_comm, // Compute symmetrized matrix 
									int _nblks, int *_blks, int *_blk2cpu) const; 
	CSMatrix SymmMtrRect () const; // Compute general symmetrized rectangular matrix
	CSMatrix SymmMtrRect2Index () const; // Compute general symmetrized rectangular matrix
	static CSMatrix CombineStructures (const CSMatrix &_a1, const CSMatrix &_a2); // Combine sparsity structures of the two matrices
	CSMatrix AtAMatrSpars  () const; // Compute the upper triangular part of the sparsity of the AtA matrix 
	CSMatrix SuperSparseMultiplyByRows (int &_icycle, int *_imask, // Multiply two matrices in supersparse format stored by rows
																int &_icycle1, int *_imask1, 
																CSMatrix &_amatr2) const;
	void SuperSparseMultiplyByRowsWithBlockList (int &_icycle, int *_imask, // Multiply two matrices in supersparse format stored by rows with block list information
																int &_icycle1, int *_imask1, 
																CSMatrix &_amatr2,
																int _iblkstr, CSMatrix &_amatr2str,
																CSMatrix &_mvm, CSMatrix &_mvmstr) const;
	void SuperSparseAddWithBlockList (int &_icycle, int *_imask, // Add the sparsity in supersparse format with block lists
																int &_icycle1, int *_imask1, 
																CSMatrix &_ablkstr,
																CSMatrix &_amtradd, CSMatrix &_ablkstradd) const;
	int UpdateAtAMatrixStruct (CSMatrix &_mtratasp, // Update mask of the current jblk block elements of the complex A^hA matrix for rectangular A
											int _iblk, int _jblk, 
											int *_blksc, int *_sp2blk, int *_imaskja);
	CSMatrix FilterMatrix (int _isupbeg, int _isupend) const; // Filter matrix according to the indices
	CSMatrix BlockSparsity (int _nblks, int *_sp2blk) const; // Compute block sparsity of a sparse matrix (incomplete lists version)
	CSMatrix BlockSparsityDiag (int _nblks, int *_sp2blk) const; // Compute block sparsity of a sparse matrix (incomplete lists version) with added diagonal
	void Transform2IndexIntoContinuous (int *_blks); // Transform 2Index matrix sparsity into one index
	void Transform2IndexIntoContinuous (CMPIComm &_comm); // Transform 2Index matrix sparsity into one index
	void ShiftSparsity (bool _is_ja_shift, int _ishift); // Shift the sparsity
	static void SparsifyExtendedBinaryTree (const CTree &_tree, // Sparsify Schur part of the extended block sparsity assuming binary tree partitioning
															int _nblks, int *_iabext, int *_jabext,
															int *_iabextflt, int *_jabextflt);
	static void SparsifyExtendedSchurBinaryTree (const CTree &_tree, // Sparsify Schur part of the extended block sparsity assuming binary tree partitioning
														int _nblks, int *_iabext, int *_jabext,
														int *_iabextflt, int *_jabextflt);
	static void ModifyPartitioning (int *_order, int _nblks, int *_blks, // Modify partitioning and ordering according to the list
												int _nlist, int *_list);
	void SortListAndColumns (); // Sort list indices and column indices
	void SortListAndColumns2Index (); // Sort list indices and column indices (2index version)
	void OptimalBoundary (int *_imask, // Compute optimal boundary
									int &_nlistbnd, int *_listbnd) const;
	bool BinarySearch (int *_imask, int *_imasknd, // For the list of nodes perform binary search of the correct combination
										int _nlist, int *_list, int _nlistchk, int *_listchk,
										int *_bitmask, int *_listfailed) const;
	int BinarySearchByPairs (int *_imask, int *_imasknd, // For the list of nodes perform binary search of the correct combination
										int _nlist, int *_list, int _nlistchk, int *_listchk,
										int *_bitmask, int *_listfailed) const;
	bool CheckNeibours (int *_imask, int *_imasknd, // For the list of nodes check their neibours for correctness
								int _nlist, int *_list) const;
	bool CheckList (int *_imask, int *_imasknd, // For the list of nodes check their correctness
							int _nlist, int *_list) const;
	void CheckList (int *_imask, int *_imasknd, // For the list of nodes check their correctness
							int _nlist, int *_list, int &_nfailed, int *_listfailed) const;
	void CheckWholes (int *_imask, int *_imasknd, // For the list of nodes check that there are no wholes
										int _nlist, int *_list, int &_nfailed, int *_listfailed) const;
	static void ComputeJa2Blk (int _nblks, int *_blks, // Compute ja2blk array with the previous row block history
										int _nlist, int *_ia, int *_ja,
										int *_ja2blk);
// Ordering functions
// Compute ordering functions
	void FindSeparator (int &_ind1, int &_ind2) const; // Find separators of the symmetric square matrix if any
	void IndependentSubsets (int &_nblks, int *&_blks, int *_order) const; // Compute the independent subsets of the matrix nodes
	void ReordMtrWeights (int *_ndweights, int *_adjweights, int *_order); // Reorder matrix and weights data in place
	void Submatrix (int _isupbeg, int _isupend, // Take current submatrix
					CSMatrix &_asub) const;
	void Submatrix (int _ilistbeg, int _ilistend, int _isupbeg, int _isupend, // Take current submatrix
							CSMatrix &_asub) const;
	void SubmatrixList (int _nlist, int *_list, // Take current submatrix according to the list
								int &_icycle, int *_imask, int *_irow2ind,
								CSMatrix &_asub) const;
	void SubmatrixByRows (int _isupbeg, int _isupend, // Take current submatrix by rows
					CSMatrix &_asub) const;
	void Submatrix2Index (int _iblkbeg, int _iblkend, // Take current submatrix
									CSMatrix &_asub) const;
	void ComputeSchur (int _nbord, int _ncycle, CSMatrix &_ashur) const; // Compute estimate of the Schur complement matrix
	void ComputeSchur (int _nblks, int *_blks, // Compute estimate of the Schur complement matrix
								int _ibegbrd, int _iendbrd, int _ncycle,
								CSMatrix &_ashur, CSMatrix &_ashur_blk) const;
	void Compute2LevelOrderingWithBlockSurfaceBoundary (int _nblks, int *_blks, // Compute ordering of the nodes into two parts with the boundary on blocks surfaces
																			int *_orderblk, int *_blkso_nobnd, int *_order, 
																			int &_ind1blk, int &_ind2blk,
																			int &_ind1, int &_ind2) const;
	void SubmatrixWeights (int _isupbeg, int _isupend, // Take current submatrix with weights
							int *_ndweights, int *_adjweights, 
							CSMatrix &_asub, int *&_ndweightssub, int *&_adjweightssub) const;
	void A2MatrWeights (int _nrows, int *_listrows, // Compute the sparsity of the AtA matrix in prescribed rows and the resulting weights assuming that matrix is symmetric
						int *_ndweights, int *_adjweights, 
						CSMatrix &_asub, int *&_ndweightssub, int *&_adjweightssub) const;
	void OrderPrfMtr     (int _ordtyp, int *_order) const; // Compute Nested Dissection or profile optimization ordering of the matrix
	void OrderNDParMetis (CMPIComm &_comm, int _nparts, // Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS
									CTree &_tree, int *&_order, int &_nblks, int *&_blks, int *&_blk2cpu) const;
	void OrderNDParMetisSchur (CMPIComm &_comm, int _nparts, // Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS with Schur partitioning
										CTree &_tree, int *&_order, int &_nblks, int *&_blks, int *&_blk2cpu) const;
	void OrderNodeNDParMetis (CMPIComm &_comm, // Compute new ND-type ordering and partioning via ParMeTiS
										int *&_order, int &_nblks, int *&_blks) const;
	void SplitNDParMetis (CMPIComm &_comm, // Split the matrix into two parts and create the corresponding ordering via ParMeTiS
									int *&_order, int &_nblks, int *_blks);
	void SplitBlockMatrixRecursive (CMPIComm &_comm, int _nparts, // Compute in parallel partitioning of the matrix into many parts taking into account blocks partitioning
												int _nblks, int *_blks, 
												int *_blk2cpu,
												int *_ordernd, int *_orderblk) const;
	void SplitBlockMatrixRecursive (int _nparts, // Compute partitioning of the matrix into many parts taking into account blocks partitioning
												int _nblks, int *_blks, 
												int *_blk2cpu) const;
	void SplitBlockMatrix (CMPIComm &_comm, // Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
									int _nblks, int *_blks,
									int *&_orderblk, int *_blkspartblk,
									int *&_ordernd, int *_blkspartnd,
									int _nproc1, int _igroup, int _myid1,
									int &_nblkslev, int *&_blkslev,
									int &_nblkslevblk, int *&_blkslevblk,
									int &_nblks1, int *&_blks1) const;
	void SplitBlockMatrix (int _nblks, int *_blks, // Compute partitioning of the matrix into two parts taking into account blocks partitioning
									int *_blk2cpu) const;
	static void ModifyDisjointParts (CMPIComm &_comm, // Modify in parallel the disjoint sets
												int *_blkscpu, 
												int *_ia, int *_ja, int *_wvertices, int *_partition);
	static void ModifyDisjointParts ( // Modify the disjoint sets
												int _ni, int *_ia, int *_ja, int *_wvertices, int *_partition);
	static void ComputeDisjointSets (CMPIComm &_comm, // Compute in parallel the set of disjoint sets for the interval and ipart
												int *_blkscpu, int _ipart,
												int *_ia, int *_ja, int *_partition,
												int &_nsets, int *&_iasets, int *&_jasets);
	static void ComputeDisjointSets (int _ibeg, int _iend, int _ipart, // Compute the set of disjoint sets for the interval and ipart
												int *_ia, int *_ja, int *_partition,
												int &_nsets, int *&_iasets, int *&_jasets);
	void ComputeTwoSidedBoundary (CMPIComm &_comm, // Compute two sided boundary in parallel mode
														int *_blkscpu, int *_partition,
														int &_nlistbnd, int *&_listbnd) const;
	CSMatrix ComputeBlockSparsity (CMPIComm &_comm, int _nblks, int *_blks) const; // Compute block sparsity
	CSMatrix ComputeBlockSparsityLocal (CMPIComm &_comm, int _nblks, int *_blks) const; // Compute block sparsity with local blocks partitioning
	CSMatrix ComputeBlockSparsityWeights (CMPIComm &_comm, int _nblks, int *_blks, // Compute block sparsity with weights and partitioning indices
														int *_blkscpu, int *_gblkscpu, int *_partition_blk,
														int *&_wvertices, int *&_wedges) const;
	CSMatrix ComputeBlockSparsityWeightsSets (CMPIComm &_comm, int _nblks, int *_blks, int _nsets, // Compute block sparsity with weights and partitioning indices
															int *_blkscpu, int *_gblkscpu, int *_partition_blk,
															int *&_wvertices, int *&_wedges) const;
	void DecompSchurBndRecursive (CMPIComm &_comm, int _decomptype, int _ncycle, int _ncolors, // Compute in parallel decomposition into the set of volume and boundary blocks
											int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
											int *_order,
											int &_nblksnew, int *&_blknew2ext, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const;
	void DecompSurfaceRecursive (CMPIComm &_comm, int _decomptype, int _ncolors, // Compute in parallel the surfaces that separate processors data
														int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
														int *_order,
														int &_nblksnew, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const;
	void DecompSurfaceRecursive1Cpu (CMPIComm &_comm1cpu, int _decomptype, // Compute on one cpu partitioning of the matrix into many parts taking into account blocks partitioning
														int _nblks, int *_blks,
														int *&_order,
														int &_nblksnew, int *&_blksnew, int *&_blk2colornew) const;
	void DecompSchurSurface (CMPIComm &_comm, int _decomptype, int _ncycle, int _ncolors, // Compute in parallel decomposition into the set of boundary blocks
										int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
										int *_order,
										int &_nblksnew, int *&_blknew2ext, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const;
	void PartSchurSubmatrix (int _decomptype, int _icolor, // Partition Schur complement submatrix by ND METIS with boundary local assignments
										int *_blk2color, int *_imaskblk, CSMatrix &_ablkref,
										int *_order,
										int &_nblksnew, int *&_blksnew, int *&_blk2colornew);
	void ComputeGroupsBoundary (CMPIComm &_comm, // Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
											int _nblks, int *_blks, int *_blk2cpu, int *_partitionblk,
											int *&_order, 
											int &_nblksnew, int *&_blksnew, int *&_blk2cpunew, int *&_partitionblknew) const;
	CSMatrix ExtendGroupsBoundary (CMPIComm &_comm, int _ncycle, // Compute in parallel the sparsity of the extended boundary
												int _nblks, int *_blks, int *_blk2cpu, int *_partitionblk) const;
	void ComputeSetOfSchurComplements (CMPIComm &_comm, int _ncycle, // Compute in parallel the sparsity of the extended boundary
													int _nblks, int *_blks, int *_blk2cpu, int *_imaskblk,
													CSMatrix *_mtrarr, CSMatrix *_mtrblkarr) const;
	static void Complete3DReduction (int _decomp_type, // Compute 3D partitioning of the surface boundary according to the partitioning type
												CSMatrix &_amtrschur, CSMatrix &_amtrblk,
												int _nblks, int *_blktype, int *_blk2cpu,
												int *_order, 
												int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew);
	static void Complete3DReduction_3Dtype0 ( // 3D partitioning of the surface boundary according to the partitioning type==0
												CSMatrix &_amtrschur, CSMatrix &_amtrblk,
												int _nblks, int *_blktype, int *_blk2cpu,
												int *_order, 
												int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew);
	static void Complete3DReduction_3Dtype1_2 (int _decomp_type, // 3D partitioning of the surface boundary according to the partitioning type==1 || type==2
												CSMatrix &_amtrschur, CSMatrix &_amtrblk,
												int _nblks, int *_blktype, int *_blk2cpu,
												int *_order, 
												int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew);
	static void Complete3DReduction_FindAndPartitionSurfacePairs ( // Compute the set of surface pairs and partition them
												CSMatrix &_amtrschur, CSMatrix &_amtrblk,
												int _nblks, int *_blktype, int *_blk2cpu,
												int &_npairs, int *&_iapairs, int *&_listpairs, int *&_cpupairs,
												int *&_order, int *&_decomp_pairs);
	void ExtendNodes (CMPIComm &_comm, // Extend the nodes in parallel mode
								int *_blkscpu, 
								int &_nlistbnd, int *&_listbnd) const;
	CSMatrix GetSubmatrixViaIntervals (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals
													int _nparts, int *_ibegarr, int *_iendarr);
	CSMatrix GetSubmatrixViaIntervals2Index (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals (2index version)
															int _nparts, int *_ibegarr, int *_iendarr);
	void ChangeIndices (CMPIComm &_comm, // Change rows/column indices
								int *_listord);
	void ChangeIndicesViaIntervals (CMPIComm &_comm, // Change rows/column indices
												int *_listord,
												int _nlisttree, int *_listtree, int *_invordertree, 
												int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart);
	static void RestoreOrderNd (int _ilev, int _nlev, // Restore local ordering at the previous level
											CMPIComm *_commlev, 
											int *_blks, int **_pporderlev, int *_nblkslev, int **_ppblkslev);
	void OrderBrdMtr     (int _nbord,  int *_order) const; // Reorder the matrix according to the bordering parameter and compute profile optimization ordering for the remaining subblock
	void PartitionMtr    (int _nparts, // Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
							int *_ndweights, int *_adjweights, int *_blks, int *_partition, int *_order) const;
	void PartitionMtr    (int _nparts, // Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
							int *_blks, int *_partition, int *_order) const;
	void PartMtr         (int _nparts, float *_weights, // Partition matrix by METIS into not equal parts with weights of each node/edge and compute its ordering
							int *_ndweights, int *_adjweights, int *_blks, int *_partition, int *_order) const;
	void PartMtr         (int _nparts, float *_weights, int *_blks, int *_partition, int *_order) const; // Partition matrix by METIS into not equal parts and compute its ordering
	void PartMtr         (int _nparts, int *_blks, int *_partition, int *_order) const; // Partition matrix by METIS and compute its ordering
	void PartOrdMtr      (int _nparts, float *_weights, int *_ndweights, int *_adjweights, int *_blks, int *_order) const; // Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks
	void PartOrdMtr      (int _nparts, float *_weights, int *_blks, int *_order) const; // Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks
	void PartOrdMtr      (int _nparts, int *_blks, int *_order) const; // Partition matrix by METIS and compute profile optimization ordering for subblocks
	void PartBinaryTreeND(CTree &_tree) const; // Partition matrix according to the prescribed binary tree, matrix ordering is assumed to be nested dessection
	void PartBinaryTreeND(CTree &_tree, int *_order, const CSlvParam &_param, // Partition matrix by ND METIS
									int &_nblks, int *&_blks, int *&_blk2cpu) const;
	void PartBinaryTreeNDSups (CTree &_tree, CTree &_treend, // Partition matrix according to the prescribed binary tree and local supernodes partitioning, internal matrix ordering is assumed to be ND
												int _nsups, int *_sprnds, 
												int *_ordersp, int *_ordernd,
												int &_nblkssp, int *&_blkssp, int *&_blk2cpusp, 
												int &_nblksnd, int *&_blksnd, int *&_blk2cpund) const;
	void PartBinaryTreeNDSchur (CTree &_tree, int *_order, const CSlvParam &_param) const; // Partition matrix by METIS including Schur complement
	void PartBinaryTreeNDSchurBnd (CTree &_tree, int *_order, const CSlvParam &_param, // Partition matrix by ND METIS including Schur complement with boundary local assignments
												int &_nblks, int *&_blks, int *&_blk2cpu) const;
	void PartBinaryTreeNDSchurBnd_FixedInitialOrdering ( // Partition matrix by ND METIS including Schur complement with boundary local assignments (fixed initial ordering on entry)
																			CTree &_tree, int *_order, const CSlvParam &_param,
																			int &_nblks, int *&_blks, int *&_blk2cpu) const;
	void PartOrdMtr      (CTree &_tree, int *_ndweights, int *_adjweights, int *_order) const; // Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks
	void PartOrdMtr      (CTree &_tree, int *_order) const; // Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks
	void PartOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_ndweights, int *_adjweights, int *_order, const CSlvParam &_param) const; // Partition matrix by METIS including Schur complement
	void PartOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_order, const CSlvParam &_param) const; // Partition matrix by METIS including Schur complement
	void PartAtAOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_order, const CSlvParam &_param) const; // Partition AtA matrix by METIS including Schur complement
	void PartitionOrdMtr (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtr (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrForSortedTree (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrForSortedTree (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrSchur (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrSchur (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrSchurForSortedTree (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	void PartitionOrdMtrSchurForSortedTree (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const;
	static void OrderViaIntervals (CMPIComm &_comm, // Compute new reordered list for the ordering presented for intervals
												int _nlist, int *_listini,
												int _nparts, int *_ibegarr, int *_iendarr, 
												int *_order,
												int *_listnew);
	static void OrderViaIntervalsAndInverseOrder (CMPIComm &_comm, // Compute new reordered list for the inverse ordering presented for intervals
																	int _nlistini, int *_listini,
																	int _nlisttree, int *_listtree, int *_invordertree, 
																	int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart,
																	int *_listord);
	static void DirOrderViaIntervals (CMPIComm &_comm, // Compute inverse order for the list described by the intervals
													int _nlist, int *_listini, int *_listord,
													int _nparts, int *_ibegarr, int *_iendarr, 
													int &_nlistpart, int *&_orderpart);
	static void InvOrderViaIntervals (CMPIComm &_comm, // Compute inverse order for the list described by the intervals
													int _nlist, int *_listini, int *_listord,
													int _nparts, int *_ibegarr, int *_iendarr, 
													int &_nlistpart, int *&_invorderpart);
	static void ChangeVectorIndices (int _nrhs, // Compute new reordered local vector data according to the new global ordering of the local data
												int _nlist, int *_list, double *_vect, int *_order);
	static void ComputeIntervals (CMPIComm &_comm, // Compute global intervals
															int _nparts, int *_ibegarr, int *_iendarr, 
															int &_nptot, int *&_ibegtot, int *&_iendtot, int *&_icputot);
// Reorder the matrix functions
	CSMatrix OrdMtr     (const int *_order) const; // Matrix ordering routine
	CSMatrix OrdMtrSymm (const int *_order) const; // Symmetric matrix ordering routine
	CSMatrix SortDataViaList (const int *_orderr) const; // Sort matrix data by ordering the list (not indices)
	CSMatrix Sort2IndexDataViaList (const int *_orderr) const; // Sort matrix data by ordering the list (not indices)
// Matrix by vector multiplication and triangular solve functions
	virtual void MvmA     (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrix " << std::endl; 
	}; 
	virtual void MvmAt    (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrix " << std::endl; 
	}; 
	virtual void MvmACol  (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrix " << std::endl; 
	};
	virtual void MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrix " << std::endl; 
	};
	virtual void MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrix " << std::endl; 
	};
	virtual void SolveL   (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrix " << std::endl; 
	};
	virtual void SolveU   (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrix " << std::endl; 
	};
	virtual void MvmA     (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrix " << std::endl; 
	}; 
	virtual void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrix " << std::endl; 
	}; 
	virtual void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrix " << std::endl; 
	};
	virtual void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrix " << std::endl; 
	};
	virtual void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrix " << std::endl; 
	};
	virtual void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrix " << std::endl; 
	};
	virtual void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrix " << std::endl; 
	};
// Input/Output functions
	virtual void PackMatrix (int &_length, char *&_obj); // Pack matrix data
	virtual void UnPackMatrix (int _length, char *_obj); // UnPack matrix data
	void A2Ps (int _collap, const char *_fname, int _Nbl, const int *_Blks) const; // PostScript collapsed sparsity structure
	void A2Ps (const char *_fname, int _Nbl, const int *_Blks) const; // PostScript sparsity structure
	static void PrintSparsity (int _collap, const char *_fname, int _n, int *_ia, int *_ja); // PostScript print of collapsed sparsity structure
	CSMatrix ReadSparsity (std::istream &_fin); // Read sparsity from disk
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrix &_a); // Output matrix
// Friend classes
	friend class CGSMatrix;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CGSMatrixCS;
	friend class CSMatrixR;
	friend class CSMatrixRS;
	friend class CSMatrixCS;
	friend class CMvmR;
	friend class CFctR;
	friend class CFctC;
	friend class CFctDiagR;
	friend class CFctDiagC;
	friend class CModel;
	friend class CDecomp;
};

#endif

// SMatrixR.h: Description of the real point matrix stored in super sparse format
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrixR
#define __SMatrixR

class CSMatrixR: public CSMatrix
{
protected:
// Data
	int nzatot; // nzatot is the total number of elements in the matrix
	int nza;    // nza is the number of elements in the matrix (array a)
	double *a;  // a[nza] array describes the elements of the matrix
public:
// Functions
// Constructors and destructor
	CSMatrixR (); // Memory allocation zero data constructor
	CSMatrixR (int _nsupa, int _nzja); // Memory allocation zero data constructor
	CSMatrixR (int _nsupa, int _nzja, int _nza); // Memory allocation zero data constructor
	CSMatrixR (int _nsupa, int _nsupa2, int _nzja, int _nzja2, int _nza); // Memory allocation zero data constructor
	CSMatrixR (int _nn); // Constructor for Laplace
	CSMatrixR (int _nn, double _dshift); // Constructor for Laplace^2 with diagonal shift
	CSMatrixR (int _nn, double _dshft, double _dunsy); // Laplace^2 structure unsymmetric matrix constructor with diagonal shift
	CSMatrixR (int _nx, int _ny, int _nz, double _dshft, double _dunsy); // Laplace^3 structure unsymmetric matrix constructor with diagonal shift
	CSMatrixR (int _nx, int _ny, int _nz, int _ibegz, int _iendz, double _dshft, double _dunsy); // Laplace structure 3D unsymmetric matrix constructor with diagonal shift
	CSMatrixR (bool _ivar, int _nx, int _ny, int _nz, int _ibegz, int _iendz, double _eps, double _dunsy); // Highly ill-conditioned 3D unsymmetric matrix constructor
	CSMatrixR (char _ch, int _n, int _lb, int _rb); // Unsymmetric banded matrix
	CSMatrixR (int _nlist, int _nlist2, int _nzja, int _nzja2); // Memory allocation zero data constructor
	CSMatrixR (const CSMatrixR &_aa); // Copy constructor
	~CSMatrixR () { // Destructor
//		cout << " On entry to CSMatrixR destructor " << endl;
		nza = 0;
		nzatot = 0;
		delete [] a;
//		cout << " On return from CSMatrixR destructor " << endl;
	};
// Get/set functions
	int   GetNza    () const {return nza;};    // Get nza
	int   GetNzatot () const {return nzatot;}; // Get nzatot
	double* GetA    () const {return a;};      // Get a
	void SetNza (int _nza) {nza = _nza;};      // Set nza
	void SetNzatot (int _nzatot) {nzatot = _nzatot;}; // Set nzatot
// Operator functions
	CSMatrixR &operator= (const CSMatrixR &_aa); // Equality operator
	CSMatrixR operator+ (const CSMatrixR &_mtr2) const; // Add two matrices
// Transformation functions
	CSMatrixR TranspMtr () const; // Compute transposed matrix 
	CSMatrixR SymmMtr   () const; // Compute symmetrized matrix 
	CSMatrixR SymmMtrRect () const; // Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes
	CSMatrixR SymmMtrRect2Index () const; // Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes
	void CombineLU (CSMatrixR &_mtrl, CSMatrixR &_mtru) const; // Compute L and U data for given matrix
	void CombineLU (CMPIComm &_comm, int _nblks, int *_blks, int *_blks2cpu, // Compute L and U data for given matrix in parallel
							CSMatrixR &_mtrl, CSMatrixR &_mtru) const;
	void SplitAccordingBlockSparsity (int _nblks, int *_blks, // Split the matrix according to the specified block sparsity
													int *_iabmax, int *_jabmax, 
													CSMatrixR &_mtr2);
	void FilterBlockRow2Index (bool _is_3index, int _nlistblk, int *_listblk, int &_icycleblk, int *_imaskblk); // Filter block row according to the list of indices
	void RemoveEquivalentRowsNumbers (); // Remove from the list of rows equivalent rows numbers
	void AddLU (CSMatrixR &_mtrl, CSMatrixR &_mtru); // Compute new matrix as a sum of L and U data for given matrix
	CSMatrixR AddMatrSS2Index (bool _is_index3, const CSMatrixR &_mtr2) const; // Add two matrices in super sparse two index format
	static void DoubleAddMatrSS2Index (bool _is_index3, // Double add two matrices in super sparse two index format
														const CSMatrixR &_mtrl, const CSMatrixR &_mtrl2, CSMatrixR &_mtrlsum,
														const CSMatrixR &_mtru, const CSMatrixR &_mtru2, CSMatrixR &_mtrusum);
	void ExtendSparsity (const CSMatrixR &_mtr2); // Extend sparsity of the matrix according to the second matrix
	void ExtendedCopy2Index (CSMatrix &_mtrsp, const CSMatrixR &_mtr2); // Copy matrix into extended sparsity
	CSMatrixR CutSparsity (const CSMatrix &_mtr2); // Cut the sparsity of the matrix according to the prescribed one
	CSMatrixR ExcludeList (int _nlist, int *_jalist); // Cut the sparsity of the matrix according to the list of index numbers
	CSMatrixR IncludeList (int _nlist, int *_jalist); // Include into the sparsity of the new martrix only elements from the list of index numbers
	CSMatrixR IncludeListNoDiag (int _nlist, int *_jalist); // Include into the sparsity of the new matrix only elements from the list of index numbers excluding diagonal ones
	void ConvertData (const CSMatrixR &_mtra, const CSVector &_x, // Perform transformation of the data for nonlinear solver
						const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, const bool *_wall_E,
						const CSVector &_norm,
						CSMatrixR &_mtranew, CSVector &_xnew, CSVector &_poly0new, CSVector &_poly2new, CSVector &_normnew) const;
	static void SuperSparsifySymm2Index (bool _cpy3index, CSMatrixR &_mtrl, CSMatrixR &_mtru); // Super sparsify symmetrically the matrices L and U by numerically zero values
	static void SuperSparsifyL2Index (bool _cpy3index, CSMatrixR &_mtrl); // Super sparsify the matrix L by numerically zero values
	void ListOfSmallDiagonalValues (double _thresh, int &_nlist, int *_list) const; // Compute the list of small diagonal elements
	void ChangeIndices (CMPIComm &_comm, int *_order); // Change rows/column indices
	void ChangeIndicesViaIntervals (CMPIComm &_comm, // Change rows/column indices
												int *_listord,
												int _nlisttree, int *_listtree, int *_invordertree, 
												int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart);
	CSMatrixR GetSubmatrixViaIntervals (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals
													int _nparts, int *_ibegarr, int *_iendarr);
	static void GetSubvectorsViaIntervals (CMPIComm &_comm, // Get the subvectors from other processors described by the list of intervals
														int _nparts, int *_ibegarr, int *_iendarr, 
														int _nrhs, int _nlist, int *_list, double *_vectini,
														int &_nlistnew, int *&_listnew, double *&_vectnew);
	static void GetSubvectorsViaEntryIntervals (CMPIComm &_comm, // Get the subvectors from other processors described by the list of entry intervals and list of final indices
														int _nparts, int *_ibegarr, int *_iendarr, 
														int _nrhs, double *_vectini,
														int _nlistnew, int *_listnew, double *_vectnew);
// Reorder the matrix functions
	CSMatrixR OrdMtr     (const int *_order) const; // Matrix ordering routine
	CSMatrixR OrdMtrSymm (const int *_order) const; // Symmetric matrix ordering routine
	CSMatrixR SortDataViaList (const int *_orderr) const; // Sort matrix data by ordering the list (not indices)
	CSMatrixR Sort2IndexDataViaList (const int *_orderr) const; // Sort matrix data by ordering the list (not indices)
	void SortListAndColumns (); // Sort list indices and column indices
// Incomplete factorization functions
	CSMatrixR Ich () const; // Incomplete Cholessky decomposition
	CSMatrixR Ich2 (CSlvParam &_param) const; // Second order Incomplete Cholessky decomposition
	CSMatrixR Ich2Dynamic (std::ofstream &_fout, CSlvParam &_param) const; // Second order Incomplete Cholessky decomposition with dynamic memory allocation
	CSMatrixR Ich2DynamicByBlocks (std::ofstream &_fout, CSlvParam &_param) const; // Second order Incomplete Cholessky decomposition with dynamic memory allocation by blocks
	void Ich2DynamicByBlocks (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete Cholessky decomposition with dynamic memory allocation by blocks
										int _istore, CGSMatrixR &_gich, CSMatrixR &_ich) const;
	void Ilu2 (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition
						const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
						CSMatrixR &_mtrl, CSMatrixR &_mtru); 
	void Ilu2Dynamic (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation
						const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
						CSMatrixR &_mtrl, CSMatrixR &_mtru);
	void Ilu2DynamicByBlocks (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation by blocks
								const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
								CSMatrixR &_mtrl, CSMatrixR &_mtru);
// Matrix by vector multiplication and triangular solve functions
	void MvmA     (const CSVector &_x, CSVector &_ax) const; // General matrix by vector multiplication
	void MvmA     (const CSVector &_x, CSVector &_ax, const bool *_barr) const; // General matrix by vector multiplication for chosen nodes
	void MvmAt    (const CSVector &_x, CSVector &_ax) const; // General transposed matrix by vector multiplication
	void MvmACol  (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixR " << std::endl; 
	};
	void MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixR " << std::endl; 
	};
	void MvmASymm (const CSVector &_x, CSVector &_ax) const; // Symmetric matrix by vector multiplication
	void SolveL   (const CSVector &_x, CSVector &_lx) const; // Solve triangular system with L stored by columns
	void SolveBlockColumnL (// Solve triangular system with block column of L
									const CSVector &_x, CSVector &_lx) const;
	void SolveU   (const CSVector &_x, CSVector &_ux) const; // Solve triangular system with U stored by rows
	void SolveBlockRowU (// Solve triangular system with the block row of U
								const CSVector &_x, 
								CSVector &_ux) const;
	void SolveLA  (const CSVector &_diagl, const CSVector &_x, CSVector &_lx, const bool *_barr) const; // Solve for chosen components triangular system with L stored by rows in the whole matrix stored by rows
	void MvmA     (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixR " << std::endl; 
	}; 
	void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixR " << std::endl; 
	}; 
	void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixR " << std::endl; 
	};
	void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixR " << std::endl; 
	};
	void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixR " << std::endl; 
	};
	void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrixR " << std::endl; 
	};
	void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrixR " << std::endl; 
	};
	void MvmBlockRowAu2Index (// Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	static void AbsMvmBlockRowAu2Index (void *_obj, int _iblk, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax);
	void MvmBlockRowAu2Index2 (// Multiply by the current block row of Au
								int *_bsx, int *_bsx2, const CSVector &_x, int *_bsax, int *_bsax2, CSVector &_ax) const;
	void MvmBlockColumnAl2Index (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	static void AbsMvmBlockColumnAl2Index (void *_obj, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax);
	void MvmBlockColumnAl2Index2 (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
								int *_bsx, int *_bsx2, const CSVector &_x, int *_bsax, int *_bsax2, CSVector &_ax) const;
	void SolveBlockColumnL2Index (// Solve triangular system with block column of L
									int *_bsx, const CSVector &_x, 
									int *_bslx, CSVector &_lx) const;
	static void AbsSolveBlockColumnL2Index (void *_obj, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVector &_x, 
								int *_bslx, CSVector &_lx);
	void SolveBlockRowU2Index (// Solve triangular system with the block row of U
								int *_bsx, const CSVector &_x, 
								int *_bsux, CSVector &_ux) const;
	static void AbsSolveBlockRowU2Index (void *_obj, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVector &_x, 
								int *_bsux, CSVector &_ux);
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrixR &_a); // Output matrix
	void WriteMtr (std::ostream &_stream); // Output matrix to disk
	void WriteMtrBin (std::ostream &_stream); // Output matrix to disk in binary form
	static CSMatrixR ReadMtr    (std::istream &_fin); // Read matrix from disk
	static CSMatrixR ReadMtrBin (std::istream &_stream); // Read matrix from disk in binary form
	static void ConvertToBinary (std::ifstream &_fin, std::ofstream &_fout); // Read matrix, right hand side and norm vector from disk and write in binary form
	static void ReadFromBinary (std::ifstream &_fin, CSMatrix &_amatr, CSVector &_rhs, CSVector &_norm); // Read matrix, right hand side and norm vector from disk
	void PackMatrix (int &_length, char *&_obj); // Pack matrix data
	void UnPackMatrix (int _length, char *_obj); // UnPack matrix data
// Friend classes
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CSMatrixRB;
	friend class CSMatrixRS;
	friend class CSMatrixC;
	friend class CSMatrixCS;
	friend class CSVector;
	friend class CMvmR;
	friend class CFctR;
	friend class CFctDiagR;
private:
// Factorization support functions
// Ich2:
	void Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
					int *&_ig, int *&_jg, int *&_addrg, double *&_g) const;
	void Ich2AllocGDynamic (int *&_ig, int **&_jg, double **&_g) const; // Allocate G data for Ich2
	void Ich2AllocGDynamicByBlocks (int *&_ig, int **&_jg, double **&_g, // Allocate G data for Ich2
									int *&_blks, int *&_nd2blk, 
									int **&_jgblk, double **&_gblk) const;
	void Ich2AllocScale (int _msglev, double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
						double *&_diagg, double *&_scl) const;
	void Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
						int *_lstloc, int _icycle, int *_imask, 
						double *_diagg, double *_fmask, double *_scl, int &_iops) const;
	void Ich2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
						int *_iv, int *_ig, int *_jg, 
						int &_nupd, int *_lstupd, int *_indlst,
						int &_nlist1, int *_lstloc, 
						int _icycle, int *_imask, double *_fmask) const;
	void Ich2InitMaskDynamic (bool _elmisr, int _irwprv, // Init mask and update arrays
								int *_iv, int *_ig, int **_jg, 
								int &_nupd, int *_lstupd, int *_indlst,
								int &_nlist1, int *_lstloc, 
								int _icycle, int *_imask, double *_fmask) const;
	void Ich2UpdateRow (int _nupd, int *_lstupd, int *_indlst, // Update current row
						double _lelem, double *_fmask, double *_g, int &_iops) const;
	void Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						int &_nlist1, int *_lstloc, int *_lstupd,
						double *_fmask, double *_diagg, int &_iops) const;
	void Ich2PivScaleRow (int _msglev, int _irow, double &_delem, // Compute diagonal pivot and scale current row
							int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, int &_iops) const;
	void Ich2FreeR (int _msglev, int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, double *_g) const;
	void Ich2FreeRDynamic (int _msglev, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int **_jg, double **_g,
							int *_list, double *_fmask) const;
	void Ich2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
									int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
									int *_iv, int *_ig, int **_jg, double **_g,
									int _iblk, int _nzalloc, int &_nzused,
									int *_blks, int *_nd2blk, int **_jgblk, double **_gblk,
									int *_list, double *_fmask) const;
	void Ich2StoreRow (int _irow, double _tau1, // Store current row
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
						double _delem, int _nlist1, int *_lstupd, double *_fmask, 
						int *_ig, int *_jg, double *_g) const;
	void Ich2StoreRowDynamic (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _delem, 
								int _nlist1, int *_lstupd, double *_fmask,
								int *_ig, int **_jg, double **_g) const;
	void Ich2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
										int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
										double _delem,
										int _nlist1, int *_lstupd, double *_fmask,
										int *_ig, int **_jg, double **_g,
										int &_iblk, int _nzallocmax, int &_nzalloc, int &_nzused,
										int *_blks, int *_nd2blk, int **_jgblk, 
										double **_gblk) const;
	void Ich2StoreScaleU (double *_scl, // Store the resulting U and scale it
							int *_ig, int *_jg, double *_g,
							CSMatrixR &_temp, double &_ops) const;
	void Ich2StoreScaleUDynamic (double *_scl, // Store the resulting U and scale it
								int *_ig, int **_jg, double **_g,
								CSMatrixR &_temp, double &_ops) const;
// Ilu2:
	void Ilu2AllocScale (double _sclmin, // Allocate double data and init scaling for Ilu2
							double *&_fmaskl, double *&_fmasku, double *&_dpiv,
							double *&_diagg, double *&_scll, double *&_sclu) const;
	void Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
						int *&_ig, int *&_jg, int *&_addrg, double *&_gl, double *&_gu) const;
	void Ilu2AllocGDynamic ( // Allocate G data for Ilu2
							int *&_ig, int **&_jg, double **&_gl, double **&_gu) const;
	void Ilu2AllocGDynamicByBlocks ( // Allocate G data for Ilu2
									int *&_ig, int **&_jg, double **&_gl, double **&_gu,
									int *&_blks, int *&_nd2blk, int **&_jgblk,
									double **&_glblk, double **&_gublk) const;
	void Ilu2InitRow (const CSMatrixR &_mtral, const CSMatrixR &_mtrau, // Init current scaled row
						int _irow, int &_nlist1,
						int *_lstloc, int _icycle, int *_imask, 
						double *_diagg, double *_fmaskl, double *_fmasku, 
						double *_scll, double *_sclu, int &_iops) const;
	void Ilu2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
						int *_iv, int *_ig, int *_jg, 
						int &_nupd, int *_lstupd, int *_indlst,
						int &_nlist1, int *_lstloc, 
						int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const;
	void Ilu2InitMaskDynamic (bool _elmisr, int _irwprv, // Init mask and update arrays
								int *_iv, int *_ig, int **_jg, 
								int &_nupd, int *_lstupd, int *_indlst,
								int &_nlist1, int *_lstloc, 
								int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const;
	void Ilu2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						int &_nlist1, int *_lstloc, int *_lstupd,
						double *_fmaskl, double *_fmasku, double *_diagg, int &_iops) const;
	void Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
							int _irow, double &_deleml, double &_delemu,
							int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
							double *_dpiv, int &_iops) const;
	void Ilu2FreeR (int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, double *_gl, double *_gu) const;
	void Ilu2FreeRDynamic (int _msglev, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int **_jg, double **_gl, double **_gu,
							int *_list, double *_fmaskl, double *_fmasku) const;
	void Ilu2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
									int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
									int *_iv, int *_ig, int **_jg, double **_gl, double **_gu,
									int _iblk, int _nzalloc, int &_nzused,
									int *_blks, int *_nd2blk, int **_jgblk, double **_glblk, double **_gublk,
									int *_list, double *_fmaskl, double *_fmasku) const;
	void Ilu2StoreRow (int _irow, double _tau1, // Store current row
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
						double _deleml, double _delemu, 
						int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
						int *_ig, int *_jg, double *_gl, double *_gu) const;
	void Ilu2StoreRowDynamic (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _deleml, double _delemu, 
								int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
								int *_ig, int **_jg, double **_gl, double **_gu) const;
	void Ilu2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
										int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
										double _deleml, double _delemu, 
										int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
										int *_ig, int **_jg, double **_gl, double **_gu,
										int &_iblk, int _nzallocmax, int &_nzalloc, int &_nzused,
										int *_blks, int *_nd2blk, int **_jgblk, 
										double **_glblk, double **_gublk) const;
};

#endif

// SMatrixRB.h: Description of the real block matrix stored in super sparse format with constant supernode size
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrixRB
#define __SMatrixRB

class CSMatrixRB: public CSMatrixR
{
protected:
// Data
	int bla;   // bla is the size of each supernode
	int bla_2; // bla_2 is the squared size of each supernode
	int bla_3; // bla_3 is the size of each supernode in power 3
public:
// Functions
// Constructors and destructor
	CSMatrixRB (); // Memory allocation zero data constructor
	CSMatrixRB (int _nsupa, int _nzja, int _bla); // Memory allocation zero data constructor
	CSMatrixRB (int _nn, int _bla); // Constructor for Laplace
	CSMatrixRB (int _nn, double _dshft, int _bla); // Constructor for Laplace^2 with diagonal shift
	CSMatrixRB (const CSMatrixRB &_aa); // Copy constructor
	~CSMatrixRB () { // Destructor
//		cout << " On entry to CSMatrixRB destructor " << endl;
		bla = 0;
		bla_2 = 0;
		bla_3 = 0;
//		cout << " On return from CSMatrixRB destructor " << endl;
	};
// Operator functions
	CSMatrixRB &operator= (const CSMatrixRB &_aa); // Equality operator
// Transformation functions
	CSMatrixRB R2RBSymm (const CSMatrixR &_ar, int _bla); // Transform matrix from point to bla format
// Incomplete factorization functions
	CSMatrixRB Ich2 (const CSlvParam _param) const; // Second order Incomplete Cholessky decomposition
// Matrix by vector multiplication and triangular solve functions
	void MvmA     (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixRB " << std::endl; 
	};
	void MvmAt    (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixRB " << std::endl; 
	};
	void MvmACol  (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixRB " << std::endl; 
	};
	void MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixRB " << std::endl; 
	};
	void MvmASymm (const CSVector &_x, CSVector &_ax) const; // Symmetric matrix by vector multiplication
	void SolveL   (const CSVector &_x, CSVector &_lx) const; // Solve triangular system with L stored by columns
	void SolveU   (const CSVector &_x, CSVector &_ux) const; // Solve triangular system with U stored by rows
	void MvmA     (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixRB " << std::endl; 
	}; 
	void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixRB " << std::endl; 
	}; 
	void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixRB " << std::endl; 
	};
	void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixRB " << std::endl; 
	};
	void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixRB " << std::endl; 
	};
	void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrixRB " << std::endl; 
	};
	void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrixRB " << std::endl; 
	};
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrixRB &_a); // Output matrix
// Friend classes
	friend class CSVector;
private:
// Factorization support functions
// Ich2:
	void Ich2AllocScale (double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
							double *&_diagg, double *&_scl, double *&_sclinv,
							double *&_lelem, double *&_delem,
							double *&_aloc, double *&_vloc, 
							double *&_eig, int &_lwork, double *&_work,
							double &_ops) const;
	void Ich2AllocScale (double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
							double *&_diagg, double *&_scl, double *&_sclinv,
							double &_ops) const;
	void Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
						int *&_ig, int *&_jg, int *&_addrg, double *&_g) const;
	void Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
						int *_lstloc, int _icycle, int *_imask, 
						double *_diagg, double *_fmask, double *_scl, 
						double *_aloc, double *_vloc,
						int &_iops) const;
	void Ich2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmask) const;
	void Ich2UpdateRow (int _nupd, int *_lstupd, int *_indlst, // Update current row
						double *_lelem, double *_fmask, double *_g, int &_iops) const;
	void Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						int &_nlist1, int *_lstloc, int *_lstupd,
						double *_fmask, double *_diagg, int &_iops) const;
	void Ich2PivScaleRow (int _irow, double *_delem, // Compute diagonal pivot and scale current row
							int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, 
							double *_aloc, double *_eig, int _lwork, double *_work,
							int &_iops) const;
	void Ich2FreeR (int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, double *_g) const;
	void Ich2StoreRow (int _irow, double _tau1, // Store current row
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
						double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
						int *_ig, int *_jg, double *_g) const;
	void Ich2StoreScaleU (double *_scl, double *_sclinv, // Store the resulting U and scale it
							int *_ig, int *_jg, double *_g,
							CSMatrixRB &_temp, 
							double *_aloc,
							double &_ops) const;
};

#endif

// SMatrixRS.h: Description of the real block matrix stored in super sparse format with variable supernode sizes
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrixRS
#define __SMatrixRS

class CSMatrixRS: public CSMatrixR
{
protected:
// Data
	bool is_bsa; // is_bsa specifies if arrays sprndr, sprndc and bsa contain some data
	int blamx;   // blamx is the maximal supernode size
	int *sprndr; // sprndr[nsupr+1] array describes the rows    supernodes partitioning of the matrix
	int *sprndc; // sprndc[nsupc+1] array describes the columns supernodes partitioning of the matrix
	int *bsa;    // bsa[nzja] array describes the base addresses of each supernodes block in the array a
public:
// Functions
// Constructors and destructor
	CSMatrixRS (); // Memory allocation zero data constructor
	CSMatrixRS (int _nsupa, int _nzja, int _nza); // Memory allocation zero data constructor
	CSMatrixRS (int _nn, int _nsupamx, const int *_sprnda); // Constructor for Laplace
	CSMatrixRS (int _nn, double _dshft, int _nsupamx, const int *_sprnda); // Constructor for Laplace^2 with diagonal shift
	CSMatrixRS (int _nn, double _dshft, double _dunsy, int _nsupamx, const int *_sprnda); // Laplace^2 structure unsymmetric matrix constructor with diagonal shift
	CSMatrixRS (const CSMatrixRS &_aa); // Copy constructor
	~CSMatrixRS () { // Destructor
//		cout << " On entry to CSMatrixRS destructor " << endl;
		delete [] sprndr;
		delete [] sprndc;
		delete [] bsa;
		blamx = 0;
//		cout << " On return from CSMatrixRS destructor " << endl;
	};
// Get/set functions
	int GetBlamx () const {return blamx;}; // Get blamx
	void GetPartition (int *&_sprndr, int *&_sprndc, int *&_bsa) const { // Get partitioning
		_sprndr = sprndr;
		_sprndc = sprndc;
		_bsa    = bsa;
	};
	int * GetSprndr () {return sprndr;}; // Get partitioning
	int * GetBsa () {return bsa;}; // Get partitioning
	void SetIsbsa (bool _is_bsa) {is_bsa = _is_bsa;}; // Set _is_bsa
	void SetBlamx (int _blamx) {blamx = _blamx;}; // Set blamx
	void SetNsupc (int _nsupc) { // Set nsupc
		nsupc = _nsupc;
		nlist = _nsupc;
		n = sprndc[nsupc];
		nzja = ia[nsupc];
	};
	void SetSprndr (int *_sprndr) {sprndr = _sprndr;}; // Set sprndr
	void SetSprndc (int *_sprndc) {sprndc = _sprndc;}; // Set sprndc
	void SetPartition (int *_sprndr, int *_sprndc, int *_bsa) { // Set partitioning
		sprndr = _sprndr;
		sprndc = _sprndc;
		bsa    = _bsa;
	};
// Operator functions
	CSMatrixRS &operator= (const CSMatrixRS &_aa); // Equality operator
// Transformation functions
	void CleanBsa (); // Clean bsa arrays
	CSMatrixRS TranspMtr () const; // Compute transposed matrix 
	CSMatrixRS SymmMtr   () const; // Compute symmetrized matrix 
	void CombineLU (CSMatrixRS &_mtrl, CSMatrixRS &_mtru) const; // Compute L and U data for given matrix
	void CombineLU (CMPIComm &_comm, int _nblks, int *_blks, int *_blks2cpu, // Compute L and U data for given matrix in parallel
									CSMatrixRS &_mtrl, CSMatrixRS &_mtru) const;
	void SortListAndColumns (); // Sort list indices and column indices
	void RemoveEquivalentRowsNumbers (); // Remove from the list of rows equivalent rows numbers
	void ExtendSparsity (const CSMatrixRS &_mtr2); // Extend sparsity of the matrix according to the second matrix
	CSMatrixRS AtAMatrElem (std::ofstream &_fout) const; // Compute the elements of the triangular part of the AtA matrix
	CSMatrixRS AtAMatrElem (std::ofstream &_fout, // Compute the elements of the triangular part of the AtA matrix
							int _nblks, int *_blks, FILE **_mtrfiles, 
							FILE **_mtratafiles) const; 
	CSMatrixRS R2RS     (const CSMatrixR &_ar, int _nsupamx, const int *_sprnda); // Transform matrix from point to sprnda format
	CSMatrixRS R2RSSymm (const CSMatrixR &_ar, int _nsupamx, const int *_sprnda); // Transform matrix from point to sprnda format
	CSMatrixRS RS2Disk  (int _nblks, int *_blks, FILE **_mtrfiles) const; // Copy matrix into disk
	void PointScaling (int *_sprndr, int *_ibsscl, double *_scl); // Perform point scaling of the block row
	void BlockScaling (int *_sprndr, int *_ibsscl, double *_scll, double *_sclu); // Perform supernode scaling of the block row
// Reordering functions
	CSMatrixRS OrdMtr        (const int *_order) const; // Matrix ordering routine
	CSMatrixRS OrdMtrSymm    (const int *_order) const; // Symmetric matrix ordering routine
	CSMatrixRS OrdMtrRectCol (const int *_order) const; // Ordering of the matrix columns
	CSMatrixRS OrdMtrRectCol (int _nblks, int *_blks, FILE **_mtrfiles, // Ordering of the matrix columns
								const int *_order, 
								FILE **_mtrofiles) const; 
// Incomplete factorization functions
	CSMatrixRS Ich2 (std::ofstream &_fout, CSlvParam &_param) const; // Second order Incomplete Cholessky decomposition
	CSMatrixRS Ich2 (std::ofstream &_fout, // Second order Incomplete Cholessky decomposition in the out-of-core
						int _nblks, int *_blks, 
						FILE **_mtrafiles, FILE **_mtrlufiles, 
						CSlvParam &_param) const;
	void Ilu2 (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition
				const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
				CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
	void Ilu2DynamicByBlocks (std::ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation by blocks
								const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
								CSMatrixRS &_mtrl, CSMatrixRS &_mtru);
// Matrix by vector multiplication and triangular solve functions
	void MvmA      (const CSVector &_x, CSVector &_ax) const;  // General matrix by vector multiplication
	void MvmAt     (const CSVector &_x, CSVector &_ax) const;  // General transposed matrix by vector multiplication
	void MvmACol   (const CSVector &_x, CSVector &_ax) const;  // General matrix by vector multiplication for matrix stored by columns
	void MvmACol   (int _nblks, int *_blks, FILE **_mtrafiles, // General matrix by vector multiplication for matrix stored by columns in out-of-core mode
					const CSVector &_x, CSVector &_ax) const; 
	void MvmAtCol  (const CSVector &_x, CSVector &_ax) const;  // General transposed matrix by vector multiplication for matrix stored by columns
	void MvmAtCol  (int _nblks, int *_blks, FILE **_mtrafiles, // General transposed matrix by vector multiplication for matrix stored by columns in out-of-core mode
					const CSVector &_x, CSVector &_ax) const; 
	void MvmASymm  (const CSVector &_x, CSVector &_ax) const;  // Symmetric matrix by vector multiplication
	void MvmBlockRowAu (int *_sprnds, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	void MvmPointBlockRowAu (int *_sprnds, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	void MvmBlockColumnAl (bool _bdiag, int *_sprnds, // Multiply by the current block column of Al with/without diagonal supernode
									int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	void MvmPointBlockColumnAl (bool _bdiag, int *_sprnds, // Multiply by the current block column of Al with/without diagonal supernode
									int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const;
	void SolveL    (const CSVector &_x, CSVector &_lx) const;  // Solve triangular system with L stored by columns
	void SolveL    (int _nblks, int *_blks, FILE **_mtrlfiles, // Solve triangular system with L stored by columns in out-of-core mode
					const CSVector &_x, CSVector &_lx) const;
	void SolveBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
									int *_bsx, const CSVector &_x, 
									int *_bslx, CSVector &_lx) const;
	void SolvePointBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
									int *_bsx, const CSVector &_x, 
									int *_bslx, CSVector &_lx) const;
	void SolveU    (const CSVector &_x, CSVector &_ux) const;  // Solve triangular system with U stored by rows
	void SolveU    (int _nblks, int *_blks, FILE **_mtrufiles, // Solve triangular system with U stored by rows in out-of-core mode
					const CSVector &_x, CSVector &_ux) const;
	void SolveBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
									int *_bsx, const CSVector &_x, 
									int *_bsux, CSVector &_ux) const;
	void SolvePointBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
									int *_bsx, const CSVector &_x, 
									int *_bsux, CSVector &_ux) const;
	void MvmA     (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixRS " << std::endl; 
	}; 
	void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixRS " << std::endl; 
	}; 
	void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixRS " << std::endl; 
	};
	void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixRS " << std::endl; 
	};
	void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixRS " << std::endl; 
	};
	void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrixRS " << std::endl; 
	};
	void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrixRS " << std::endl; 
	};
// Input/Output functions
	void PackMatrix (int &_length, char *&_obj); // Pack matrix data
	void UnPackMatrix (int _length, char *_obj); // UnPack matrix data
	void OutputMatrix (std::ostream _stream, int _nblks, int *_blks, std::fstream *_mtrfiles) const; // Output matrix
	void WriteMtrBin (std::ostream &_stream); // Output matrix to disk in binary form
	void WriteMtrBin (FILE * _file, int &_offset); // Output matrix to disk in binary form
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrixRS &_a); // Output matrix
// Friend classes
	friend class CSMatrixCS;
	friend class CGSMatrixRS;
	friend class CSVector;
	friend class CFctR;
	friend class CFctDiagR;
	friend class CModel;
	friend class CMvmR;
private:
// Factorization support functions
// Ich2/Ilu2:
	void Ich2AllocScale (std::ofstream &_fout, double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
							int *&_bsdia, double *&_diagg, double *&_scl, double *&_sclinv,
							double *&_lelem, double *&_delem,
							double *&_aloc, double *&_vloc, 
							double *&_eig, int &_lwork, double *&_work,
							double &_ops) const;
	void Ich2AllocScale (std::ofstream &_fout, // Allocate double data and init scaling for Ich2
							int _nblks, int *_sp2blk, FILE **_mtrafiles,
							double *&_gread, int &_nzgread,
							double *&_fmask, double *&_dpiv, 
							int *&_bsdia, double *&_diagg, 
							int *&_bsblk, FILE **_mtrlufiles,
							int *&_addrscl, int *&_addrsclinv,
							double *&_lelem, double *&_delem,
							double *&_aloc, double *&_vloc, 
							double *&_eig, int &_lwork, double *&_work,
							double &_ops) const;
	void Ilu2AllocScale (std::ofstream &_fout, int _msglev, double _sclmin, // Allocate double data and init scaling for Ilu2
							double *&_fmaskl, double *&_fmasku, double *&_dpiv, 
							int *&_bsdia, double *&_diagg, 
							double *&_scll, double *&_sclu, double *&_sclinvl, double *&_sclinvu,
							double *&_lelem, double *&_uelem, double *&_deleml, double *&_delemu,
							double *&_aloc, double *&_uloc, double *&_vloc, 
							double *&_eig, int &_lwork, double *&_work,
							double &_ops) const;
	void Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
						int *&_ig, int *&_jg, int *&_addrg, double *&_g) const;
	void Ich2AllocG (double _memory, int &_nzjgmx, // Allocate G data for Ich2
						int *&_ig, int *&_jg, int *&_addrg) const;
	void Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
						int *&_ig, int *&_jg, int *&_addrg, double *&_gl, double *&_gu) const;
	void Ilu2AllocGDynamicByBlocks (int *&_ig, int *&_blks, int *&_sp2blk, // Allocate G data for Ilu2
									int **&_jg, int **&_addrg, double **&_gl, double **&_gu,
									int **&_jgblk, int **&_addrgblk, double **&_glblk, double **&_gublk) const;
	void Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
						int *_lstloc, int _icycle, int *_imask, 
						int *_bsdia, double *_diagg, double *_fmask, double *_scl, 
						double *_aloc, double *_vloc,
						int &_iops) const;
	void Ich2InitRow (int _irow, // Init current scaled row
						int *_sp2blk, FILE *_mtrafile, FILE **_mtrlufiles,
						int &_nlist1, 
						int *_lstloc, int _icycle, int *_imask, 
						int *_bsdia, double *_diagg, double *_fmask, int *_addrscl, double *_delem,
						double *_aloc, double *_vloc,
						int &_iops) const;
	void Ilu2InitRow (const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, // Init current scaled row
						int _irow, int &_nlist1, 
						int *_lstloc, int _icycle, int *_imask, 
						int *_bsdia, double *_diagg, double *_fmaskl, double *_fmasku, 
						double *_scll, double *_sclu,
						double *_aloc, double *_vloc,
						int &_iops) const;
	void Ich2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
						int *_iv, int *_ig, int *_jg, 
						int &_nupd, int *_lstupd, int *_indlst,
						int &_nlist1, int *_lstloc, 
						int _icycle, int *_imask, double *_fmask) const;
	void Ilu2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
						int *_iv, int *_ig, int *_jg, 
						int &_nupd, int *_lstupd, int *_indlst,
						int &_nlist1, int *_lstloc, 
						int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const;
	void Ilu2InitMaskDynamicByBlocks (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const;
	void Ich2UpdateRow (int _irow, int _irwprv, int _nupd, int *_lstupd, int *_indlst, // Update current row
						double *_lelem, double *_fmask, int *_addrg, double *_g, int &_iops) const;
	void Ich2UpdateRow (int _irow, // Update current row
						FILE *_mtrlufile, double *_gread,
						int _irwprv, int _nupd, int *_lstupd, int *_indlst, 
						double *_lelem, double *_fmask, int *_addrg, int &_iops) const;
	void Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
								int &_nlist1, int *_lstloc, int *_lstupd,
								double *_fmask, int *_bsdia, double *_diagg, int &_iops) const;
	void Ilu2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
						int &_nlist1, int *_lstloc, int *_lstupd,
						double *_fmaskl, double *_fmasku, int *_bsdia, double *_diagg, int &_iops) const;
	void Ich2PivScaleRow (int _irow, double *_delem, // Compute diagonal pivot and scale current row
							int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, 
							double *_aloc, double *_eig, int _lwork, double *_work,
							int &_iops) const;
	void Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
							int _irow, double *_deleml, double *_delemu, 
							int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, double *_dpiv, 
							double *_aloc, double *_eig, double *_uloc, double *_vloc, int _lwork, double *_work,
							int &_iops) const;
	void Ich2FreeR (std::ofstream &_fout, int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, int &_nznew, int *_addrg, double *_g) const;
	void Ich2FreeR (int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, int &_nznew, int *_addrg) const;
	void Ilu2FreeR (std::ofstream &_fout, int _msglev, int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, int &_nznew, 
					int *_addrg, double *_gl, double *_gu) const;
	void Ilu2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
									int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
									int *_iv, int *_ig, int **_jg, int **_addrg, double **_gl, double **_gu,
									int _iblk, 
									int *_blks, int *_sp2blk, int **_jgblk, int **_addrgblk, double **_glblk, double **_gublk) const;
	void Ich2StoreRow (int _irow, double _tau1, // Store current row
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
						double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
						int *_ig, int *_jg, int &_nznew, int *_addrg, double *_g) const;
	void Ich2StoreRow (int _irow, // Store current row
						FILE *_mtrlufile, double *_aloc,
						double _tau1,
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
						double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
						int *_ig, int *_jg, int &_nznew, int &_nzgblk, int *_addrg) const;
	void Ilu2StoreRow (int _irow, double _tau1, // Store current row
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
						double *_deleml, double *_delemu, int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
						int *_ig, int *_jg, int &_nznew, int *_addrg, double *_gl, double *_gu) const;
	void Ilu2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
										int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
										double *_deleml, double *_delemu, int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
										int *_ig, int **_jg, int **_addrg, double **_gl, double **_gu,
										int &_iblk, int _nzjallocmax, int &_nzjalloc, int &_nzjused, int _nzallocmax, int &_nzalloc, int &_nzused,
										int *_blks, int *_sp2blk, int **_addrgblk, int **_jgblk, 
										double **_glblk, double **_gublk) const;
	void Ich2StoreScaleU (int *_bsdia, double *_scl, double *_sclinv, // Store the resulting U and scale it
							int *_ig, int *_jg, int *_addrg, double *_g,
							CSMatrixRS &_temp, 
							double *_aloc,
							double &_ops) const;
	void Ich2StoreScaleU (int _nblks, int *_blks, FILE **_mtrlufiles, // Store the resulting U and scale it
							int *_sp2blk, int *_bsdia, int *_addrscl, int *_addrsclinv, int *_bsblk,
							int *_ig, int *_jg, int *_addrg,
							CSMatrixRS &_temp, 
							double *_aloc, double *_vloc, double *_delem, 
							double &_ops) const;
	void Ich2StoreScaleUDynamic (int *_bsdia, double *_scl, double *_sclinv, // Store the resulting U and scale it
								int *_ig, int **_jg, int **_addrg, double **_g,
								CSMatrixRS &_temp, 
								double *_aloc,
								double &_ops) const;
};

#endif

// SMatrixC.h: Description of the complex point matrix stored in super sparse format
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrixC
#define __SMatrixC

class CSMatrixC: public CSMatrix
{
protected:
// Data
	int nzatot; // nzatot is the total number of elements in the matrix
	int nza;    // nza is the number of elements in the matrix (array a)
	dcmplx *a;  // a[nza] array describes the elements of the matrix
private:
public:
// Functions
// Constructors and destructor
	CSMatrixC (); // Memory allocation zero data constructor
	CSMatrixC (int _nsupa, int _nzja); // Memory allocation zero data constructor
	CSMatrixC (int _nsupa, int _nzja, int _nza); // Memory allocation zero data constructor
	CSMatrixC (int _nsupa, int _nsupa2, int _nzja, int _nzja2, int _nza); // Memory allocation zero data constructor
	CSMatrixC (const CSMatrixC &_aa); // Copy constructor
	~CSMatrixC () { // Destructor
//		std::cout << " On entry to CSMatrixC destructor " << std::endl;
		nza = 0;
		nzatot = 0;
		delete [] a;
//		std::cout << " On return from CSMatrixC destructor " << std::endl;
	};
// Get/set functions
	int GetNza    () const {return nza;};    // Get nza
	int GetNzatot () const {return nzatot;}; // Get nzatot
	dcmplx * GetA    () const {return a;};   // Get a
	void SetNza (int _nza) {nza = _nza;};             // Set nza
	void SetNzatot (int _nzatot) {nzatot = _nzatot;}; // Set nzatot
// Operator functions
	CSMatrixC &operator= (const CSMatrixC &_aa); // Equality operator
// Transformation functions
	CSMatrixC TranspMtr () const; // Compute transposed matrix 
	CSMatrixC TranspMtrNoConj () const; // Compute transposed matrix without complex conjugation
	static CSMatrixC SparsifyRectangularMatrix (double _alpha, double _thresh_col, // Perform matrix sparsification according to the column theshold (Tikhonov regularization block is included)
																int _m, int _n, dcmplx *_amatr);
// Reorder the matrix functions
	CSMatrixC OrdMtrRectCol (const int *_order) const; // Ordering of the matrix columns
// Incomplete factorization functions
// Matrix by vector multiplication and triangular solve functions
	void MvmA     (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixC " << std::endl; 
	}; 
	void MvmAt    (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixC " << std::endl; 
	}; 
	void MvmACol  (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixC " << std::endl; 
	};
	void MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixC " << std::endl; 
	};
	void SolveL   (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrixC " << std::endl; 
	};
	void SolveU   (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrixC " << std::endl; 
	};
	void MvmA     (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixC " << std::endl; 
	};
	void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void SolveL   (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
	void SolveU   (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In MvmACol  CSMatrixC " << std::endl; 
	};
// Input/Output functions
	void PackMatrix (int &_length, char *&_obj); // Pack matrix data
	void UnPackMatrix (int _length, char *_obj); // UnPack matrix data
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrixC &_a); // Output matrix
// Friend classes
	friend class CSMatrixCS;
	friend class CGSMatrixCS;
	friend class CSVector;
	friend class CModel;
};

#endif

// SMatrixCS.h: Description of the complex block matrix stored in super sparse format with variable supernode sizes
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SMatrixCS
#define __SMatrixCS

class CSMatrixCS: public CSMatrixC
{
protected:
// Data
	int blamx;   // blamx is the maximal supernode size
	int *sprndr; // sprndr[nsupr+1] array describes the rows    supernodes partitioning of the matrix
	int *sprndc; // sprndc[nsupc+1] array describes the columns supernodes partitioning of the matrix
	int *bsa;    // bsa[nzja] array describes the base addresses of each supernodes block in the array a
private:
// Factorization support functions
// Ilu2:
	void Ilu2AllocScale (std::ofstream &_fout, double _sclmin, // Allocate double complex data and init scaling for Ilu2
							dcmplx *&_fmaskl, dcmplx *&_fmasku, double *&_dpiv, 
							int *&_bsdia, dcmplx *&_diagg, 
							dcmplx *&_scll, dcmplx *&_sclu, dcmplx *&_sclinvl, dcmplx *&_sclinvu,
							dcmplx *&_lelem, dcmplx *&_uelem, dcmplx *&_deleml, dcmplx *&_delemu,
							dcmplx *&_aloc, dcmplx *&_uloc, dcmplx *&_vloc, 
							double *&_eig, int &_lwork, dcmplx *&_work, double *&_rwork, 
							double &_ops) const;
	void Ilu2AllocScale (std::ofstream &_fout, double _sclmin, // Allocate double complex data and init scaling for Ilu2
							int _nblks, int *_sp2blk, FILE **_mtraufiles,
							dcmplx *&_fmaskl, dcmplx *&_fmasku, double *&_dpiv, 
							int *&_bsdia, dcmplx *&_diagg, 
							int *&_bsblk, FILE **_mtrlfiles, FILE **_mtrufiles,
							int *&_addrscll, int *&_addrsclu, int *&_addrsclinvl, int *&_addrsclinvu,
							dcmplx *&_gread, dcmplx *&_lelem, dcmplx *&_uelem, dcmplx *&_deleml, dcmplx *&_delemu,
							dcmplx *&_aloc, dcmplx *&_uloc, dcmplx *&_vloc, 
							double *&_eig, int &_lwork, dcmplx *&_work, double *&_rwork, 
							double &_ops) const;
	void Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
						int *&_ig, int *&_jg, int *&_addrg, dcmplx *&_gl, dcmplx *&_gu) const;
	void Ilu2AllocG (double _memory, int &_nzjgmx, // Allocate G data for Ilu2
						int *&_ig, int *&_jg, int *&_addrg) const;
	void Ilu2InitRow (const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau, // Init current scaled row
						int _irow, int &_nlist1, 
						int *_lstloc, int _icycle, int *_imask, 
						int *_bsdia, dcmplx *_diagg, dcmplx *_fmaskl, dcmplx *_fmasku, 
						dcmplx *_scll, dcmplx *_sclu,
						dcmplx *_aloc, dcmplx *_vloc,
						int &_iops) const;
	void Ilu2InitRow (int *_sp2blk, FILE *_mtralfile, FILE *_mtraufile, // Init current scaled row
						const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau, 
						FILE **_mtrlfiles, FILE **_mtrufiles,
						int _irow, int &_nlist1, 
						int *_lstloc, int _icycle, int *_imask, 
						int *_bsdia, dcmplx *_diagg, dcmplx *_fmaskl, dcmplx *_fmasku, dcmplx *_gread,
						int *_addrscll, int *_addrsclu,
						dcmplx *_aloc, dcmplx *_vloc,
						int &_iops) const;
	void Ilu2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
						int *_iv, int *_ig, int *_jg, 
						int &_nupd, int *_lstupd, int *_indlst,
						int &_nlist1, int *_lstloc, 
						int _icycle, int *_imask, dcmplx *_fmaskl, dcmplx *_fmasku) const;
	void Ich2UpdateRow (int _irow, int _irwprv, int _nupd, int *_lstupd, int *_indlst, // Update current row
						dcmplx *_lelem, dcmplx *_fmask, int *_addrg, dcmplx *_g, int &_iops) const;
	void Ich2UpdateRow (int _irow, // Update current row
						FILE *_mtrlufile, dcmplx *_gread,
						int _irwprv, int _nupd, int *_lstupd, int *_indlst,
						dcmplx *_lelem, dcmplx *_fmask, int *_addrg, int &_iops) const;
	void Ilu2FiltrRow (int _irow, int _fcttyp, double _theta, double _tau2, // Perform filtering of the current row
						int *_iaini, int *_jaini, int &_icycleini, int *_imaskini,
						int &_nlist1, int *_lstloc, int *_lstupd,
						dcmplx *_fmaskl, dcmplx* _fmasku, int *_bsdia, dcmplx *_diagg, int &_iops) const;
	void Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
							int _irow, dcmplx *_deleml, dcmplx *_delemu, 
							int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, double *_dpiv, 
							dcmplx *_aloc, double *_eig, dcmplx *_uloc, dcmplx *_vloc, 
							int _lwork, dcmplx *_work, double *_rwork, int &_iops) const;
	void Ilu2FreeR (std::ofstream &_fout, int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, int &_nznew, 
					int *_addrg, dcmplx *_gl, dcmplx *_gu) const;
	void Ilu2FreeR (std::ofstream &_fout, int _irow, // Free R part of the data
					int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
					int *_iv, int *_ig, int *_jg, int &_nznew, 
					int *_addrg) const;
	void Ilu2StoreRow (int _irow, int _fcttyp, double _tau1, // Store current row
						int *_iaini, int *_jaini, int &_icycleini, int *_imaskini,
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
						dcmplx *_deleml, dcmplx *_delemu, int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, 
						int *_ig, int *_jg, int &_nznew, int *_addrg, dcmplx *_gl, dcmplx *_gu) const;
	void Ilu2StoreRow (int _irow, double _tau1, // Store current row
						FILE *_mtrlfile, FILE *_mtrufile, dcmplx *_gread,
						int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
						dcmplx *_deleml, dcmplx *_delemu, int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, 
						int *_ig, int *_jg, int &_nznew, int &_nzgblk, int *_addrg) const;
	void Ich2StoreScaleU (int *_bsdia, dcmplx *_scl, dcmplx *_sclinv, // Store the resulting U and scale it
							int *_ig, int *_jg, int *_addrg, dcmplx *_g,
							CSMatrixCS &_temp, 
							dcmplx *_aloc,
							double &_ops) const;
	void Ich2StoreScaleU (int _nblks, int *_blks, FILE **_mtrlufiles, // Store the resulting U and scale it
							int *_sp2blk, int *_bsdia, int *_addrscl, int *_addrsclinv, int *_bsblk,
							int *_ig, int *_jg, int *_addrg,
							CSMatrixCS &_temp, 
							dcmplx *_aloc, dcmplx *_gread, dcmplx *_delem,
							double &_ops) const;
public:
// Functions
// Constructors and destructor
	CSMatrixCS (); // Memory allocation zero data constructor
	CSMatrixCS (int _nsupa, int _nzja, int _nza); // Memory allocation zero data constructor
	CSMatrixCS (int _nn, double _dshft, double _dunsyr, double _dunsyi, // Laplace^2 structure unsymmetric complex matrix constructor with diagonal shift
				double _alphar, double _alphai, 
				int _nsupamx, const int *_sprnda); 
	CSMatrixCS (const CSMatrixCS &_aa); // Copy constructor
	~CSMatrixCS () { // Destructor
//		std::cout << " On entry to CSMatrixCS destructor " << std::endl;
		delete [] sprndr;
		delete [] sprndc;
		delete [] bsa;
		blamx = 0;
//		std::cout << " On return from CSMatrixCS destructor " << std::endl;
	};
// Get/set functions
	int GetBlamx () const {return blamx;}; // Get blamx
	void GetPartition (int *&_sprndr, int *&_sprndc, int *&_bsa) const { // Get partitioning
		_sprndr = sprndr;
		_sprndc = sprndc;
		_bsa    = bsa;
	};
	int * GetSprndr () {return sprndr;}; // Get sprndr
	int * GetSprndc () {return sprndc;}; // Get sprndc
	int * GetBsa () {return bsa;}; // Get bsa
	void SetBlamx (int _blamx) {blamx = _blamx;}; // Set blamx
	void SetNsupc (int _nsupc) { // Set nsupc
		nsupc = _nsupc;
		nlist = _nsupc;
		n = sprndc[nsupc];
		nzja = ia[nsupc];
	};
	void SetSprndr (int *_sprndr) {sprndr = _sprndr;}; // Set sprndr
	void SetSprndc (int *_sprndc) {sprndc = _sprndc;}; // Set sprndc
// Operator functions
	CSMatrixCS &operator= (const CSMatrixCS &_aa); // Equality operator
// Transformation functions
	CSMatrixCS TranspMtr () const; // Compute transposed matrix 
	CSMatrixCS TranspMtrNoConj () const; // Compute transposed matrix without conjugation
	CSMatrixCS TranspMtrRowsNoConj () const; // Compute transposed to the matrix stored by rows without conjugation
	CSMatrixCS SymmMtr   () const; // Compute symmetrized matrix 
	void CombineLU (CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const; // Compute L and U data for given matrix
	void CombineLU (int _nblks, int *_blks, FILE **_mtrfiles, // Compute L and U data for given matrix
					FILE **_mtrlfiles, FILE **_mtrufiles,
					CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const;
	void CombineLU (CMPIComm &_comm, int _nblks, int *_blks, int *_blks2cpu, // Compute L and U data for given matrix in parallel
									CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const;
	void FilterBlockRow (const CGSMatrix &_gmtra, // Filter block row according to the list of indices
								int _nlistblk, int *_listblk, int *_sp2blk,
								int &_icycleblk, int *_imaskblk);
	void SortListAndColumns (); // Sort list indices and column indices
	void RemoveEquivalentRowsNumbers (); // Remove from the list of rows equivalent rows numbers
	void ExtendSparsity (const CSMatrixCS &_mtr2); // Extend sparsity of the matrix according to the second matrix
	CSMatrixCS AtAMatrElem (std::ofstream &_fout, CSlvParam &_param) const; // Compute the elements of the triangular part of the AtA matrix
	CSMatrixCS AtAMatrElem (std::ofstream &_fout, // Compute the elements of the triangular part of the AtA matrix
							int _nblks, int *_blks, FILE **_mtrfiles, 
							FILE **_mtratafiles) const; 
	static CSMatrixCS C2CS     (const CSMatrixC &_ac, int _nsupamx, const int *_sprnda); // Transform matrix from point to sprnda format
	static CSMatrixCS C2CSRectCol (const CSMatrixC &_ac, int _nsupamx, const int *_sprnda); // Transform matrix from point to sprnda format for rectangular matrix stored by columns
	CSMatrixCS RS2CS (CSMatrixRS &_ar, CSMatrixRS &_ai, // Construct complex matrix from two real matrices with the same sparsity
						double _alphar, double _alphai);
	CSMatrixCS CS2Disk (int _nblks, int *_blks, FILE **_mtrfiles) const; // Copy matrix into disk
// Reordering and scaling functions
	CSMatrixCS OrdMtr (const int *_order) const; // Matrix ordering routine
	CSMatrixCS OrdMtr (int _nblks, int *_blks, FILE **_mtrfiles, // Matrix ordering routine
						const int *_order, 
						FILE **_mtrofiles) const;
	CSMatrixCS OrdMtrSymm    (const int *_order) const; // Symmetric matrix ordering routine
	CSMatrixCS OrdMtrRectCol (const int *_order) const; // Ordering of the matrix columns
	CSMatrixCS OrdMtrRectCol (int _nblks, int *_blks, FILE **_mtrfiles, // Ordering of the matrix columns
								const int *_order, 
								FILE **_mtrofiles) const; 
	void PointScaling (int *_sprndr, int *_ibsscl, dcmplx *_scl); // Perform point scaling of the block row
	void BlockScaling (int *_sprndr, int *_ibsscl, dcmplx *_scll, dcmplx *_sclu); // Perform supernode scaling of the block row
// Incomplete factorization functions
	void Ilu2 (std::ofstream &_fout, const CSlvParam _param, // Second order Incomplete LU decomposition
				const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
				CSMatrixCS &_mtrl, CSMatrixCS &_mtru);
	void Ilu2 (std::ofstream &_fout, const CSlvParam _param, // Second order Incomplete LU decomposition
				int _nblks, int *_blks, 
				FILE **_mtralfiles, FILE **_mtraufiles, 
				const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
				FILE **_mtrlfiles, FILE **_mtrufiles, 
				CSMatrixCS &_mtrl, CSMatrixCS &_mtru);
// Matrix by vector multiplication and triangular solve functions
	void MvmA     (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmA     CSMatrixCS " << std::endl; 
	}; 
	void MvmAt    (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixCS " << std::endl; 
	}; 
	void MvmACol  (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixCS " << std::endl; 
	};
	void MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixCS " << std::endl; 
	};
	void MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixCS " << std::endl; 
	};
	void SolveL   (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns
		int i = _x.nv;
		i = _lx.nv;
		std::cout << " In SolveL   CSMatrixCS " << std::endl; 
	};
	void SolveU   (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows
		int i = _x.nv;
		i = _ux.nv;
		std::cout << " In SolveU   CSMatrixCS " << std::endl; 
	};
	void MvmA      (const CSVectorC &_x, CSVectorC &_ax) const; // General matrix by vector multiplication
	void MvmA      (int _nblks, int *_blks, FILE **_mtrafiles,  // General matrix by vector multiplication in out-of-core mode
					const CSVectorC &_x, CSVectorC &_ax) const;
	void MvmBlockRowAu    (int *_sprnds, // Mupliply by the current block row    of Au
							int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmBlockRowAuConj (int *_sprnds, // Mupliply by the complex conjugate current block row of Au
							int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmBlockColumnAl (int *_sprnds, // Mupliply by the current block column of Al without diagonal supernode
							int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmBlockColumnAlConj (int *_sprnds, // Mupliply by the complex conjugate current block column of Al without diagonal supernode
								int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmBlockColumnA (int *_sprndr, int *_sprndc, // Mupliply by the current block column of A
									int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmBlockColumnAConj (int *_sprndr, int *_sprndc, // Mupliply by the complex conjugate current block column of rectangular A
										int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const;
	void MvmAh     (const CSVectorC &_x, CSVectorC &_ax) const; // General hermittian transposed matrix by vector multiplication
	void SolveL    (const CSVectorC &_x, CSVectorC &_lx) const; // Solve triangular system with L stored by columns
	void SolveL    (int _nblks, int *_blks, FILE **_mtrlfiles,  // Solve triangular system with L stored by columns in out-of-core mode
					const CSVectorC &_x, CSVectorC &_lx) const;
	void SolveBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
							int *_bsx, const CSVectorC &_x, 
							int *_bslx, CSVectorC &_lx) const;
	void SolveBlockColumnLConj (int *_sprnds, // Solve triangular system with the complex conjugate block column of L
										int *_bsx, const CSVectorC &_x, 
										int *_bslx, CSVectorC &_lx) const;
	void SolveU    (const CSVectorC &_x, CSVectorC &_ux) const; // Solve triangular system with U stored by rows
	void SolveU    (int _nblks, int *_blks, FILE **_mtrufiles,  // Solve triangular system with U stored by rows in out-of-core mode
					const CSVectorC &_x, CSVectorC &_ux) const;
	void SolveBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
							int *_bsx, const CSVectorC &_x, 
							int *_bsux, CSVectorC &_ux) const;
	void SolveBlockRowUConj (int *_sprnds, // Solve triangular system with the complex conjugate block row of U
								int *_bsx, const CSVectorC &_x, 
								int *_bsux, CSVectorC &_ux) const;
	void SolveLUOrd (int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru, // Solve system with ordered preconditioner
							CSVectorC &_x, CSVectorC &_px);
	void MvmAt    (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAt    CSMatrixCS " << std::endl; 
	}; 
	void MvmACol  (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmACol  CSMatrixCS " << std::endl; 
	};
	void MvmAtCol (const CSVectorC &_x, CSVectorC &_ax) const { // General transposed matrix by vector multiplication for matrix stored by columns
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmAtCol CSMatrixCS " << std::endl; 
	};
	void MvmASymm (const CSVectorC &_x, CSVectorC &_ax) const { // Symmetric matrix by vector multiplication
		int i = _x.nv;
		i = _ax.nv;
		std::cout << " In MvmASymm CSMatrixCS " << std::endl; 
	};
// Input/Output functions
	void PackMatrix (int &_length, char *&_obj); // Pack matrix data
	void UnPackMatrix (int _length, char *_obj); // UnPack matrix data
	void OutputMatrix (std::ostream _stream, // Output matrix
						int _nblks, int *_blks, FILE **_mtrfiles) const;
	friend std::ostream &operator<< (std::ostream &_stream, const CSMatrixCS &_a); // Output matrix
// Friend classes
	friend class CGSMatrixCS;
	friend class CSVector;
	friend class CSVectorC;
	friend class CFctC;
	friend class CFctDiagC;
	friend class CModel;
};

#endif
