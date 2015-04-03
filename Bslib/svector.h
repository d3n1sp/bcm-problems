//------------------------------------------------------------------------------------------------
// File: svector.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// svector.h: Description of real sparse vector
//
//////////////////////////////////////////////////////////////////////////////

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

#ifndef __SVector
#define __SVector

// Preliminary declarations

typedef double (*FUNCVALUE) (void *_obj, int _itype, int _param, double _value);

class CTree;
class CMvmR;
class CSMatrix;
class CSMatrixR;
class CSMatrixRS;
class CSMatrixRB;
class CGSMatrix;
class CGSMatrixR;
class CGSMatrixRS;
class CSVector;
class CSlvParam;

typedef int  (*CGETM)     (void *_obj); // Define interface of the rows    object size
typedef int  (*CGETN)     (void *_obj); // Define interface of the columns object size
typedef void (*CMVM)      (void *_obj, CSVector &_x, CSVector &_ax);  // Define interface of the matrix by vector multiplication
typedef void (*CMVMT)     (void *_obj, CSVector &_x, CSVector &_ax);  // Define interface of the complex conjugate matrix by vector multiplication
typedef void (*CSOLVEL)   (void *_obj, CSVector &_x, CSVector &_lx);  // Define interface of the solve with the lower triangular part of the preconditioner
typedef void (*CSOLVEU)   (void *_obj, CSVector &_x, CSVector &_ux);  // Define interface of the solve with the upper triangular part of the preconditioner
typedef void (*CSOLVELU)  (void *_obj, CSVector &_x, CSVector &_lux); // Define interface of the solve with the preconditioner (including l and u parts)
typedef void (*CSOLVEAL)  (void *_obj, CSVector &_x, CSVector &_alx); // Define interface of the solve with the lower triangular part of the coefficient matrix (SOR)

class CSVector 
{
// Data
	int nv;       // nv is the vector length
	int nrhs;     // nrhs is the total number of right hand sides
	double *vect; // vect [nv*nrhs] array contains the elements of the vector
	int nparts;   // nparts is the number of sparse parts in the vector
	int *blksv;   // blksv [nparts+1] array describes partitioning of the sparse vector
	int *inda;    // inda [nparts] array describes first positions of each dense part of the vector in the original matrix
public:
// Functions
// Constructors and destructor
	CSVector (); // Zero length vector constructor
	CSVector (int _n); // Memory allocation constructor with zero data
	CSVector (int _n, int _nrhs); // Memory allocation constructor with zero data
	CSVector (int _n, int _np, int _nrhs); // One more memory allocation constructor with zero data
	CSVector (const CGSMatrix &_mtra, const CSVector &_x); // Constructor for parallel mode
	CSVector (const CSVector &_vv); // Copy constructor
	~CSVector (); // Destructor
// Get/set functions
	int GetNv       () const {return nv;};     // Get nv
	int GetNrhs     () const {return nrhs;};   // Get nrhs
	int GetNparts   () const {return nparts;}; // Get nparts
	double* GetVect () const {return vect;};   // Get vect
	void SetNv   (int _nv)   {nv = _nv;};      // Set nv
	void SetNrhs (int _nrhs) {nrhs = _nrhs;};  // Set nrhs
// Set functions
	void SetSVect (double _value); // Set a float value for all components
	void SetComponentSVect (int _index, double _value); // Set a float value for component
// Operator functions
	CSVector operator+ (const CSVector &_v2); // Add two vectors
	CSVector operator- (const CSVector &_v2); // Substract two vectors
	CSVector operator* (double _alpha); // Multiply vector by a scalar
	CSVector &operator= (const CSVector &_v2); // Equality operator
	void DaxpySVect (double _alpha, const CSVector &_x); // Compute y += alpha*x
	void DaxpySVectBlk (double _alpha, const double *_a, const CSVector &_x); // Compute y += alpha*x*a
	void ScaleSVect (double _alpha); // Compute x *= alpha
	double CompProdSVect (int _i, const CSVector &_v2); // Component product
	double ScProdSVect (const CSVector &_v2); // Scalar product
	double MaxProdSVect (const CSVector &_v2) const; // Maximum of component products
	double MaxProdSVect (const CSVector &_v2, const bool *_barr) const; // Maximum of chosen component products
// Ordering functions
	void OrdVect (char _ordtyp, const int *_order, const CSVector &_x); // Reorder the vector
	void OrdVect (char _ortype, // Order the vector according to the supernodes partitioning
					int _nsupc, int *_sprndc, int *_sprndcnew, int *_order);
// Iterative schemes
	void Cg (std::ofstream &_fout, // Perform preconditioned CG iterations
				const CSMatrix *_pa, int _istore, CGSMatrix *_pglu, const CSMatrix *_plu, 
				const CSVector &_x, const CSVector &_b, 
				CSlvParam &_param, const CSVector &_norm);
	void Cg (std::ofstream &_fout, // Perform preconditioned CG iterations
				const CTree &_tree,
				void *_obja, CMVM _mvm,
				void *_objlu, CSOLVELU _solvelu,
				const CSVector &_x, const CSVector &_b, 
				CSlvParam &_param);
	void Lanczos (std::ofstream &_fout, // Perform preconditioned Lanczos iterations
					const CSMatrix *_mtra, const CSMatrix *_mtrl, const CSMatrix *_mtru, 
					const CSVector &_x, const CSVector &_b, 
					CSlvParam &_param, const CSVector &_norm);
//	void Lanczos (std::ofstream &_fout, // Perform preconditioned Lanczos iterations
//					const CTree &_tree, CMvmR &_mvm,
//					CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau, 
//					bool _second, CMvmR &_mvm2, CGSMatrixRS &_gmtra2,
//					CGSMatrixRS &_gmtrl, CGSMatrixRS &_gmtru,
//					const CSVector &_x, const CSVector &_b, const CSVector &_norm,
//					CSlvParam &_param); 
	void Lanczos (std::ofstream &_fout, // Perform preconditioned Lanczos iterations
					const CTree &_tree, 
					void *_obja, CMVM _mvm, CMVM _mvmt,
					void *_objl, CSOLVEL _solvel, CSOLVEL _solvelt,
					void *_obju, CSOLVEU _solveu, CSOLVEU _solveut,
					const CSVector &_x, const CSVector &_b, const CSVector &_norm,
					CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
										const CTree &_tree, 
										void *_obja, CMVM _mvm,
										void *_objlu, CSOLVELU _solvelu,
										const CSVector &_x, const CSVector &_b);
	void BlockLanczosCentral (std::ofstream &_fout, // Perform central preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
								const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, 
								const CSVector &_x, const CSVector &_b, 
								const CSlvParam &_param);
	void BlockLanczosLeft (std::ofstream &_fout, // Perform left preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
							const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, 
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param);
	void BlockLanczosRight (std::ofstream &_fout, // Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
							const CSMatrixRS *_mtra, const CSMatrixRS *_mtrl, const CSMatrixRS *_mtru, 
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param);
	void LanczosLsq (std::ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrlu, // Perform preconditioned Lanczos iterations with rectangular matrix
						const CSVector &_x, const CSVector &_b, 
						const CSlvParam &_param);
	void BlockLanczosLsq (std::ofstream &_fout, const CSMatrixRS *_mtra, const CSMatrixRS *_mtrlu, // Perform preconditioned block Lanczos iterations with rectangular matrix
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param);
	void LanczosLsq (std::ofstream &_fout, // Perform preconditioned Lanczos iterations with rectangular matrix in out-of-core mode
						int _nblks, int *_blks, 
						const CSMatrixRS *_mtra, FILE **_mtrafiles, 
						const CSMatrixRS *_mtrlu, FILE **_mtrlufiles, 
						const CSVector &_x, const CSVector &_b, 
						const CSlvParam &_param);
	void BlockLanczosLsq (std::ofstream &_fout, // Perform preconditioned block Lanczos iterations with rectangular matrix in the out-of-core mode
							int _nblks, int *_blks, 
							const CSMatrixRS *_mtra, FILE **_mtrafiles, 
							const CSMatrixRS *_mtrlu, FILE **_mtrlufiles, 
							const CSVector &_x, const CSVector &_b, 
							const CSlvParam &_param);
	void ComputeRhs (int _typePol, const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, // Compute right hand side for nonlinear solver
						const bool *_cmprhs, const CSVector &_sol);
	void NonlinearSor (const CSMatrixR &_mtra, const CSVector &_x, // Perform iterations of the nonlinear SOR algorithm
								const CSVector &_poly0, const CSVector &_poly1, const CSVector &_poly2, const bool *_wall_E,
								CSlvParam &_param, const CSVector &_norm);
	void NonlinearSor (std::ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
								CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau,
								bool _second, CMvmR &_mvm2, CGSMatrixRS &_gmtra2,
								int _itype, int _iparam, void **_objs, FUNCVALUE *_funcs, 
								CSVector &_poly0, CSVector &_poly2, 
								CSVector &_x, CSVector &_norm, 
								CSlvParam &_param);
	void NonlinearSor2Index (std::ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
										CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
										bool _second, CMvmR &_mvm2, CGSMatrixR &_gmtra2,
										int _itype, 
										int _iparam, void **_objs, FUNCVALUE *_funcs, 
										CSVector &_poly0, CSVector &_poly2, 
										CSVector &_x, CSVector &_norm, 
										CSlvParam &_param);
	void NonlinearSorSchur2Index (std::ofstream &_fout, const CTree &_tree, // Perform in parallel sor iterations for general or quadratic nonlinearity in the rhs (diagonal one)
											CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
											void *_obja, CMVM _mvm,
											void *_objal, CSOLVEAL _solveal,
											int _itype, 
											int _iparam, void **_objs, FUNCVALUE *_funcs, 
											CSVector &_poly0, CSVector &_poly2, 
											CSVector &_x, CSVector &_norm, 
											CSlvParam &_param);
	void ConvertDataBack (const CSVector &_x, const CSVector &_sol, const bool *_wall_E); // Perform transformation of the data back from the nonlinear solver
// Iterative solvers
	void Solver (std::ofstream &_fout, CSMatrixR &_mtra, // Solve unsymmetric system
					const CSVector &_b, 
					CSlvParam &_param, const CSVector &_norm);
	void Solver (std::ofstream &_fout, CSMatrixRS &_mtra, // Solve unsymmetric system
					const CSVector &_b, 
					CSlvParam &_param, const CSVector &_norm);
	void Solver (std::ofstream &_fout, const CTree &_tree, // Solve unsymmetric system in parallel mode
						CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau,
						CGSMatrixRS &_gmtralfct, CGSMatrixRS &_gmtraufct,
						bool _second, CGSMatrixRS &_gmtra2,
						CGSMatrixRS &_gmtrl, CGSMatrixRS &_gmtru,
						const CSVector &_b, const CSVector &_x, CSVector &_norm, 
						CSlvParam &_param);
	void Solver2Ind (std::ofstream &_fout, const CTree &_tree, // Solve unsymmetric system in parallel mode
							CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
							CGSMatrixR &_gmtralfct, CGSMatrixR &_gmtraufct,
							bool _second, CGSMatrixR &_gmtra2,
							CGSMatrixR &_gmtrl, CGSMatrixR &_gmtru,
							const CSVector &_b, const CSVector &_x, CSVector &_norm, 
							CSlvParam &_param);
	void SolverSymm2IndSchur (std::ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
										CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
										CGSMatrixR &_gmtralfct,
										bool _second, CGSMatrixR &_gmtra2,
										CGSMatrixR &_gmtrl,
										const CSVector &_b, const CSVector &_x, CSVector &_norm, 
										CSlvParam &_param);
	void Solver2IndSchur (std::ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
									CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
									CGSMatrixR &_gmtralfct, CGSMatrixR &_gmtraufct,
									bool _second, CGSMatrixR &_gmtra2,
									CGSMatrixR &_gmtrl, CGSMatrixR &_gmtru,
									const CSVector &_b, const CSVector &_x, CSVector &_norm, 
									CSlvParam &_param);
	void SolveSymm (std::ofstream &_fout, CSMatrixR &_mtra, // Solve symmetric system
					const CSVector &_b, 
					CSlvParam &_param, const CSVector &_norm);
	void SolveSymm (std::ofstream &_fout, CSMatrixRS &_mtra, // Solve symmetric system
						const CSVector &_b, 
						CSlvParam &_param, const CSVector &_norm);
	void SolveLsq (std::ofstream &_fout, CSMatrixRS &_mtra, // Solve Lsq problem
					const CSVector &_b, 
					CSlvParam &_param);
	void SolveLsq (std::ofstream &_fout, // Solve Lsq problem in the out-of-core mode
					int _nblks, int *_blks, FILE **_mtrfiles,
					const CSMatrixRS &_mtra,
					char *_path,
					const CSVector &_b, CSlvParam &_param);
	void NewtonSolve (std::ofstream &_fout, CSMatrixR &_mtra, // Perform Newton iterations for quadratic nonlinearity in the rhs
						CSVector &_x, CSVector &_poly0, CSVector &_poly2, 
						CSVector &_norm, CSlvParam &_param);
	void NewtonSolve (std::ofstream &_fout, CSMatrixR &_mtra, // Perform Newton iterations for general nonlinearity in the rhs (with one discontinuity)
						int _iparam, void **_objs, FUNCVALUE *_funcs,
						CSVector &_x, CSVector &_poly0, 
						CSVector &_norm, CSlvParam &_param);
	void NewtonSolve (std::ofstream &_fout, const CTree &_tree,  CMvmR &_mvm, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
							CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau,
							bool _second, CMvmR &_mvm2, CGSMatrixRS &_gmtra2,
							int _itype, 
							int _iparam, void **_objs, FUNCVALUE *_funcs, 
							CSVector &_poly0, CSVector &_poly2, 
							CSVector &_x, CSVector &_norm, 
							CSlvParam &_param);
	void NewtonSolve2Index (std::ofstream &_fout, const CTree &_tree, CMvmR &_mvm, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
									CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
									bool _second, CMvmR &_mvm2, CGSMatrixR &_gmtra2,
									int _itype, 
									int _iparam, void **_objs, FUNCVALUE *_funcs, 
									CSVector &_poly0, CSVector &_poly2, 
									CSVector &_x, CSVector &_norm, 
									CSlvParam &_param);
	void NewtonSolveSchur2Index (std::ofstream &_fout, CTree &_tree, // Perform in parallel Newton iterations for general or quadratic nonlinearity in the rhs (diagonal one)
											CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau,
											bool _second, CGSMatrixR &_gmtra2,
											int _itype, 
											int _iparam, void **_objs, FUNCVALUE *_funcs, 
											CSVector &_poly0, CSVector &_poly2, 
											CSVector &_x, CSVector &_norm, 
											CSlvParam &_param);
// Input and output
	static CSVector ReadVect    (std::istream &_fin, int _nv); // Read vector from disk
	static CSVector ReadVectBin (std::istream &_stream); // Read vector from disk in binary form
	friend std::ostream &operator<< (std::ostream &_stream, const CSVector &_v); // Output vector
	void WriteVect (std::ostream &_stream); // Output vector to the disk
	void WriteVectBin (std::ostream &_stream); // Output vector to the disk in a binary form
	void WriteVectBin (FILE * _file, int &_offset); // Output vector to the disk in a binary form
// Friend classes
	friend class CMvmR;
	friend class CFctDiagR;
	friend class CSMatrix;
	friend class CSMatrixR;
	friend class CSMatrixRB;
	friend class CSMatrixRS;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CSMatrixC;
	friend class CSMatrixCS;
	friend class CSlvParam;
};

#endif

// svectorc.h: Description of complex sparse vector
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __SVectorC
#define __SVectorC

// Preliminary declarations

class CSMatrixC;
class CSMatrixCS;
class CGSMatrixCS;
class CSlvParam;
class CMvmC;
class CQrdC;
class CCorrectorC;
class CTree;

class CSVectorC;

typedef void (*CMVMC)     (void *_obj, CSVectorC &_x, CSVectorC &_ax);  // Define interface of the matrix by vector multiplication
typedef void (*CMVMHC)    (void *_obj, CSVectorC &_x, CSVectorC &_ax);  // Define interface of the complex conjugate matrix by vector multiplication
typedef void (*CSOLVELC)  (void *_obj, CSVectorC &_x, CSVectorC &_lx);  // Define interface of the solve with the lower triangular part of the preconditioner
typedef void (*CSOLVEUC)  (void *_obj, CSVectorC &_x, CSVectorC &_ux);  // Define interface of the solve with the upper triangular part of the preconditioner
typedef void (*CSOLVELUC) (void *_obj, CSVectorC &_x, CSVectorC &_lux); // Define interface of the solve with the preconditioner (including l and u parts)
typedef CSMatrixCS (*CGETSPARSITYCS) (void *_obj, CSlvParam &_params);  // Define interface of the sparse matrix generation

class CSVectorC 
{
// Data
	int nv;       // nv is the vector length
	int nrhs;     // nrhs is the total number of right hand sides
	dcmplx *vect; // vect [nv*nrhs] array contains the elements of the vector
	int nparts;   // nparts is the number of sparse parts in the vector
	int *blksv;   // blksv [nparts+1] array describes partitioning of the sparse vector
	int *inda;    // inda [nparts] array describes first positions of each dense part of the vector in the original matrix
public:
// Functions
// Constructors and destructor
	CSVectorC (); // Zero length vector constructor
	CSVectorC (int _n); // Memory allocation constructor with zero data
	CSVectorC (int _n, int _nrhs); // Memory allocation constructor with zero data
	CSVectorC (int _n, int _np, int _nrhs); // One more memory allocation constructor with zero data
	CSVectorC (const CGSMatrixCS &_mtra, const CSVectorC &_x); // Constructor for parallel mode
	CSVectorC (const CSVectorC &_vv); // Copy constructor
	~CSVectorC (); // Destructor
// Get/set functions
	int GetNv        () const {return nv;};     // Get nv
	int GetNrhs      () const {return nrhs;};   // Get nrhs
	int GetNparts    () const {return nparts;}; // Get nparts
	int * GetBlksV   () const {return blksv;};  // Get blksv
	int * GetIndA    () const {return inda;};   // Get inda
	dcmplx * GetVect () const {return vect;};   // Get vect
// Set functions
	void SetNv   (int _nv)   {nv = _nv;};       // Set nv
	void SetNrhs (int _nrhs) {nrhs = _nrhs;};   // Set nrhs
	void SetSVect (dcmplx _value); // Set a complex value for all components
	void SetComponentSVect (int _index, dcmplx _value); // Set a complex value for component
// Operator functions
	CSVectorC operator+ (const CSVectorC &_v2); // Add two vectors
	CSVectorC operator- (const CSVectorC &_v2); // Substract two vectors
	CSVectorC operator* (dcmplx _alpha); // Multiply vector by a scalar
	CSVectorC &operator= (const CSVectorC &_v2); // Equality operator
	void DaxpySVect (dcmplx _alpha, const CSVectorC &_x); // Compute y += alpha*x
	void DaxpySVectBlk (dcmplx _alpha, const dcmplx *_a, const CSVectorC &_x); // Compute y += alpha*x*a
	void ScaleSVect (dcmplx _alpha); // Compute x *= alpha
	double CompProdSVect (int _i, const CSVectorC &_v2); // Component product
	dcmplx ScProdSVect (const CSVectorC &_v2); // Scalar product
	double MaxProdSVect (const CSVectorC &_v2) const; // Maximum of component products
	double MaxProdSVect (const CSVectorC &_v2, const bool *_barr) const; // Maximum of chosen component products
	void BlockDaxpy (const CSVectorC &_x, int _lda, dcmplx *_a); // Perform block Daxpy y += x*a
	void BlockDaxpy (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a); // Perform block Daxpy y += x*a
	void BlockDot   (const CSVectorC &_x, int _lda, dcmplx *_a) const; // Perform block Dot a = y^h*x
	void BlockDot   (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a) const; // Perform block Dot a = y^h*x
	void BlockDotH  (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a) const; // Perform block Dot a = x^h*y
	void MvCorr (const CTree &_tree, // Perform multiplication by the corrector matrix
					const CSVectorC &_ycorr, const CSVectorC &_wcorr,
					int _lddcorr, dcmplx *_dcorr, 
					int _ldcoef, dcmplx *_coef1, dcmplx *_coef2);
	void MvCorr (const CTree &_tree, // Perform multiplication by the corrector matrix
						const CCorrectorC &_corr,
						int _ldcoef, dcmplx *_coef1, dcmplx *_coef2);
	void MvmAhAOrdShiftSymm (void *_obj, CGETM _getm, CGETN _getn, CMVMC _mvmc, CMVMHC _mvmhc, // Complex hermitian matrix by vector multiplication
										double _dshift, int *_order,
										const CSVectorC &_x, CSVectorC &_ax) const;
	static void OrderVector (int _ordtype, int _n, int *_order, dcmplx *_x, dcmplx *_ax); // Reorder vector components
	static void MvmBlock  (int _m, int _n, dcmplx *_ablk, dcmplx *_x, dcmplx *_ax); // Multiply by the rectangular block and add
	static void MvmBlockH (int _m, int _n, dcmplx *_ablk, dcmplx *_x, dcmplx *_ax); // Multiply by the hermitian transposed rectangular block and add
// Ordering and scaling functions
	void OrdVect (char _ordtyp, const int *_order, const CSVectorC &_x); // Reorder the vector
	void OrdVect (char _ortype, // Order the vector according to the supernodes partitioning
					int _nsupc, int *_sprndc, int *_sprndcnew, int *_order);
	void ScaleVect (char _sctype, // Perform point scaling of the vector according to the blocks partitioning
					const CSVectorC &_sclpt);
	void BlockScaleVect (char _transp, const CGSMatrixCS &_mtra, // Perform block scaling of the vector according to the blocks partitioning
							const CSVectorC &_sclblk);
// Iterative schemes
	void Cg (std::ofstream &_fout, // Perform preconditioned CG iterations
					const CTree &_tree, CMvmC &_mvm,
					CGSMatrixCS &_gmtrau, CGSMatrixCS &_gmtru,
					const CSVectorC &_x, const CSVectorC &_b, 
					const CSlvParam &_param);
	void BlockLanczosRight (std::ofstream &_fout, // Perform right preconditioned block Lanczos iterations with square ILU2 preconditioned matrix
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
							CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param);
	void BlockLanczosRight (std::ofstream &_fout, // Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
									const CTree &_tree, 
									void *_obja, CMVMC _mvm, CMVMHC _mvmh,
									void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
									const CSVectorC &_x, const CSVectorC &_b, 
									const CSlvParam &_param);
	void BlockLanczosReorth (std::ofstream &_fout, // Perform right preconditioned block Lanczos iterations with rectangular ICH2 preconditioned matrix
									const CTree &_tree, CMvmC &_mvm,
									CGSMatrixCS &_gmtra,
									CGSMatrixCS &_gmtru,
									const CSVectorC &_x, const CSVectorC &_b, 
									const CSlvParam &_param);
	void SOFLanczosNE_old (std::ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
								const CTree &_tree, 
								void *_obja, CMVMC _mvm, CMVMHC _mvmh,
								void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
								const CSVectorC &_x, const CSVectorC &_b);
	void SOFLanNEFilterRelations_old (std::ofstream &_fout, // Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
											const CTree &_tree, const CTree &_treeclp, 
											void *_obja, CMVMC _mvm, CMVMHC _mvmh,
											void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
											const CSlvParam &_param,
											int &_ibegblk, int &_iendblk, int &_inextblk, 
											int _m, int _n, 
											CQrdC &_qrdu, CQrdC &_qrdv, dcmplx *_ppu, dcmplx *_ppv, 
											dcmplx *_tauu, dcmplx *_tauv,
											int _ldr, dcmplx *_rmatr, dcmplx *_rcolmn,
											CSVectorC &_u, CSVectorC &_v, CSVectorC &_pwork,
											dcmplx *_givarr, dcmplx *_giv, dcmplx *_giv1, 
											double &_svmin_r, double &_svmax_r);
	void SOFLanczosNE (std::ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned SOFLanNE iterations with rectangular ICH2 preconditioned matrix
									const CTree &_tree, CCorrectorC &_corr,
									void *_obja, CMVMC _mvm, CMVMHC _mvmh,
									void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
									const CSVectorC &_x, const CSVectorC &_b);
	void SOFLanNEFilterRelations (std::ofstream &_fout, // Perform filtering of the SOFLanNE matrix relations according to solution and singular values/vectors
													const CTree &_tree, const CTree &_treeclp, 
													void *_obja, CMVMC _mvm, CMVMHC _mvmh,
													void *_objlu, CSOLVELC _solvel, CSOLVEUC _solveu,
													const CSlvParam &_param,
													int &_ibegblk, int &_iendblk, int &_inextblk, 
													int _m, int _n, 
													CQrdC &_qrdu, CQrdC &_qrdv, dcmplx *_ppu, dcmplx *_ppv, 
													dcmplx *_tauu, dcmplx *_tauv,
													int _ldr, dcmplx *_rmatr, dcmplx *_rcolmn,
													CSVectorC &_r0, CSVectorC &_u, CSVectorC &_v, CSVectorC &_pwork,
													dcmplx *_givarr, dcmplx *_giv, dcmplx *_giv1,
													double &_svmin_r, double &_svmax_r);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
						const CSMatrixCS *_mtra, const CSMatrixCS *_mtrl, const CSMatrixCS *_mtru, 
						const CSVectorC &_x, const CSVectorC &_b, 
						const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix in the out-of-core mode
						int _nblks, int *_blks, FILE **_mtrafiles, const CSMatrixCS *_mtra, 
						FILE **_mtrlfiles, FILE **_mtrufiles, const CSMatrixCS *_mtrl, const CSMatrixCS *_mtru, 
						const CSVectorC &_x, const CSVectorC &_b, 
						const CSlvParam &_param, FILE **_qfiles);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
						const CGSMatrixCS &_gmtral, const CGSMatrixCS &_gmtrau, 
						const CGSMatrixCS &_gmtrl, const CGSMatrixCS &_gmtru,
						const CSVectorC &_x, const CSVectorC &_b, 
						const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							int _ivar, 
							CGSMatrixCS &_gmtrau, 
							double _dshift, CGSMatrixCS &_gmtrarect, int *_order, 
							CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							void *_obj, CGETM _getm, CGETN _getn, CMVMC _mvmc, CMVMHC _mvmhc,
							double _dshift, int *_order, 
							CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							CGSMatrixCS &_gmtra, 
							int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
								const CTree &_tree, 
								CGSMatrixCS &_gmtra, 
								int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
								const CSVectorC &_x, const CSVectorC &_b, 
								const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, const CSlvParam &_param, // Perform right preconditioned OOC block GMRES iterations with square ILU2 preconditioned matrix
							const CTree &_tree, 
							void *_obja, CMVMC _mvmc,
							void *_objlu, CSOLVELUC _solvelu,
							const CSVectorC &_x, const CSVectorC &_b);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
						const CTree &_tree, CMvmC &_mvm,
						CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
						CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
						const CSVectorC &_x, const CSVectorC &_b, 
						const CSlvParam &_param);
	void BlockGMRES (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2 preconditioned matrix
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtrau, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, 
							const CSlvParam &_param);
	void BlockGMRESCorr (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix
							const CTree &_tree, CMvmC &_mvm,
							CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
							CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
							const CSVectorC &_x, const CSVectorC &_b, CSVectorC &_ycorr,
							const CSlvParam &_param);
	int BlockGMRESCorr (std::ofstream &_fout, // Perform right preconditioned block GMRES iterations with square ILU2+Corr preconditioned matrix
								const CTree &_tree, CMvmC &_mvm,
								CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
								CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
								const CSVectorC &_x, const CSVectorC &_b, 
								CCorrectorC &_corr,
								const CSlvParam &_param);
	void BlockCorrector (std::ofstream &_fout, // Compute corrector subspace directions
							const CTree &_tree,
							int _n, int _nlocext, int _nrhs, int _iter,
							CQrdC &_qrd, dcmplx *_p, int _ldp, dcmplx *_taup,
							int _ldrmatr, dcmplx *_rmatr, dcmplx *_givarr,
							CSVectorC &_ycorr,
							const CSlvParam &_param, double &_ops);
	void BlockCorrector (std::ofstream &_fout, // Compute corrector
							const CTree &_tree, 
							int _n, int _nlocext, int _nrhs, int _iter, 
							CQrdC &_qrd, dcmplx *_p, int _ldp, dcmplx *_taup,
							int _ldrmatr, dcmplx *_rmatr, dcmplx *_givarr,
							CCorrectorC &_corr,
							const CSlvParam &_param, double &_ops);
// Iterative solvers
	void Solver (std::ofstream &_fout, CSMatrixCS &_mtra, // Solve unsymmetric system
					int _job, int *&_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
					const CSVectorC &_b, 
					const CSlvParam &_param);
	void Solver (std::ofstream &_fout, // Solve unsymmetric system in the out-of-core mode
					int _nblks, int *_blks, FILE **_mtrafiles,
					const CSMatrixCS &_mtra,
					char *_path,
					int _job, int *&_order, 
					FILE **_mtrlfiles, FILE **_mtrufiles,
					CSMatrixCS &_mtrl, CSMatrixCS &_mtru,
					const CSVectorC &_b, 
					const CSlvParam &_param);
	void Solver (std::ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
					CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau,
					int _job, CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
					const CSVectorC &_b, CSVectorC &_ycorr, 
					CSlvParam &_param);
	int Solver (std::ofstream &_fout, CTree &_tree, // Solve unsymmetric system in parallel mode
					CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau,
					int _job, CGSMatrixCS &_gmtrl, CGSMatrixCS &_gmtru,
					const CSVectorC &_b, CCorrectorC &_corr, 
					const CSlvParam &_param);
	void SolveA (std::ofstream &_fout, // Solve linear system
									CGSMatrixCS &_amatr,
									CSVectorC &_b, 
									CSlvParam &_param);
	void SolveA (std::ofstream &_fout, // Solve linear system
						void *_obj, CGETM _getm, CGETSPARSITYCS _getsparsitycs, CMVMC _mvmc,
						CSVectorC &_b, 
						CSlvParam &_param);
	void SolveAhA (std::ofstream &_fout, double _dshift, // Solve rectangular least squares problem with diagonal shift
								CGSMatrixCS &_amatr,
								CSVectorC &_b, 
								CSlvParam &_param);
	void SolveAhA (std::ofstream &_fout, double _dshift, // Solve rectangular least squares problem with diagonal shift
								void *_obj, CGETM _getm, CGETN _getn, CGETSPARSITYCS _getsparsitycs, CMVMC _mvmc, CMVMHC _mvmhc,
								CSVectorC &_b, 
								CSlvParam &_param);
	void CollectSolutionOnRoot (const CTree &_tree, const CGSMatrixCS &_gmtra, // Collect solution vector on root processor
								int *_bl2cpu, int _nrhs,
								dcmplx *_solloc, dcmplx *_soltot);
// Input and output
	friend std::ostream &operator<< (std::ostream &_stream, const CSVectorC &_v); // Output vector
// Friend classes
	friend class CMvmC;
	friend class CFctDiagC;
	friend class CGSMatrixCS;
	friend class CSMatrix;
	friend class CSMatrixR;
	friend class CSMatrixRB;
	friend class CSMatrixRS;
	friend class CSMatrixC;
	friend class CSMatrixCS;
	friend class CSlvParam;
#ifndef __FV_2004__
	friend class CModel;
#endif
};

#endif

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
