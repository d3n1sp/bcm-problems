//------------------------------------------------------------------------------------------------
// File: parfbss_data.h
// Creator: Kharchenko S.A. 
// Goal: Stores system of equations data and solves system of equations in parallel mode
//------------------------------------------------------------------------------------------------

#include "ExchangeMPI.h"
#include "tree.h"
#include "smatrix.h"
#include "gsmatrix.h"
#include "slvparam.h"

#ifndef __PARFBSS_DATA_H__
#define __PARFBSS_DATA_H__

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

class CParFbssData;

class ParFBSS_Array
{
public:
//const static int NumberOfObjects = 20;
	enum {			 NumberOfObjects = 20};
	int IRegistration [NumberOfObjects];
	CParFbssData * PRegistration [NumberOfObjects];
	ParFBSS_Array () {
		int i;
		for (i=0;i<NumberOfObjects;i++) {
			IRegistration[i] = 0;
		};
	};
// Functions
	static int ParFBSS_Create (void *_Comm, int _NI, int _NZAI, CParFbssData *&_PObj);
	static int ParFBSS_ReCreate (void *_PObj, int _NI, int _NZAI);
	static int ParFBSS_LinearEquationsData (void *_PObj, int *&_LIST, int *&_IA, int *&_JA, double *&_A, 
															double *&_RightHandSide, double *&_InitialGuess);
	static int ParFBSS_NormalizationData (void *_PObj, double *&_Norm);
	static int ParFBSS_NonLinearEquationsData (void *_PObj, double *&_Poly2);
	static int ParFBSS_GetParameters (void *_PObj, void *&_Params);
	static int ParFBSS_SetSolverParameters (void *_PObj, int _NpDistr, 
															double _Tau1, double _Tau2, double _Theta, 
															int _NIterMax, double _Eps);
	static int ParFBSS_SetMsglevAndFilename (void *_PObj, int _Msglev, char *_FileName);
	static int ParFBSS_SetSchurtype (void *_PObj, int _Schurtype, int _Ncycle);
	static int ParFBSS_SetNcpuext (void *_PObj, int _Ncpuext);
	static int ParFBSS_SolveSymm (void *_PObj, int _Job, int *_NIterPerf);
	static int ParFBSS_Solve (void *_PObj, int _Job, int *_NIterPerf);
	static int ParFBSS_RestartSolve_Ini (void *_PObj, int _Job);
	static int ParFBSS_RestartSolve_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_RestartSolveSymm_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_RestartSolve_NonLinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_RestartSolve_Fin (void *_PObj);
	static int ParFBSS_Solve_Ini (void *_PObj);
	static int ParFBSS_Solve_Distrib (void *_PObj);
	static int ParFBSS_Solve_GetDistrib (void *_PObj);
	static int ParFBSS_Solve_Transform (void *_PObj);
	static int ParFBSS_Solve_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_SolveSymm_LinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_Solve_NonLinearIter (void *_PObj, int *_NIterPerf, double *_ResFin);
	static int ParFBSS_Solve_ExchangeSolution (void *_PObj);
	static int ParFBSS_Solution (void *_PObj, double *&_Solution);
	static int ParFBSS_Distribution (void *_PObj, int &_NIOpt, int *&_LISTOpt);
	static int ParFBSS_Destroy (CParFbssData *_PObj);
};

class CSMatrixR;
class CSVector;
class CTree;
class CGSMatrixR;
class CSlvParam;

class CParFbssData
{
public:
// Functions
	CParFbssData();
	CParFbssData (void *_comm, int _ni, int _nzai);
	~CParFbssData(); // Destructor
// Interface functions
// Get functions
	CMPIComm * GetComm () {return comm;}; // Get comm
	int GetNprocsolve () {return nprocsolve;}; // Get nprocsolve
	int GetNtot () {return ntot;}; // Get ntot
	int GetNlistini () {return nlistini;}; // Get nlistini
	int * GetIorderini () {return iorderini;}; // Get iorderini
	int * GetListini () {return listini;}; // Get listini
	int * GetListord () {return listord;}; // Get listord
	int GetNppart () {return nppart;}; // Get nppart
	int *GetIbegpartarr () {return ibegpartarr;}; // Get ibegpartarr
	int *GetIendpartarr () {return iendpartarr;}; // Get iendpartarr
	int GetNhcells () {return nhcells;}; // Get nhcells
	int * GetHcells () {return hcells;}; // Get hcells
	int * GetHcells2cpu () {return hcells2cpu;}; // Get hcells2cpu
	int * GetOrderpart () {return orderpart;}; // Get orderpart
	int GetNlisttree () {return nlisttree;}; // Get nlisttree
	int * GetInvordertree () {return invordertree;}; // Get invordertree
	int * GetListtree () {return listtree;}; // Get listtree
	CSMatrixR * GetAmtr () {return amtr;}; // Get amtr
	CSVector * GetRhs () {return rhs;}; // Get rhs
	CSVector * GetX () {return x;}; // Get x
	CSVector * GetNorm () {return norm;}; // Get norm
	CSVector * GetPoly2 () {return poly2;}; // Get poly2
	CSVector * GetSol () {return sol;}; // Get sol
	CTree * GetTree () {return tree;}; // Get tree
	CSlvParam * GetParams () {return param;}; // Get param
	CGSMatrixR * GetGmtral () {return pgmtral;}; // Get gmtral
	CGSMatrixR * GetGmtrau () {return pgmtrau;}; // Get gmtrau
	CGSMatrixR * GetGmtra2 () {return pgmtra2;}; // Get gmtra2
	CSVector * GetRhsnew () {return prhsnew;}; // Get rhsnew
	CSVector * GetXnew () {return pxnew;}; // Get xnew
	CSVector * GetNormnew () {return pnormnew;}; // Get normnew
	CSVector * GetPoly2new () {return ppoly2new;}; // Get poly2new
	CSVector * GetSolnew () {return psolnew;}; // Get solnew
// Set functions
	void SetNprocsolve (int _nprocsolve) {nprocsolve = _nprocsolve;}; // Set nprocsolve
	void SetNtot (int _ntot) {ntot = _ntot;}; // Set ntot
	void SetNlistini (int _nlistini) {nlistini = _nlistini;}; // Set nlistini
	void SetIorderini (int *_iorderini) {iorderini = _iorderini;}; // Set iorderini
	void SetListini (int *_listini) {listini = _listini;}; // Set listini
	void SetListord (int *_listord) {listord = _listord;}; // Set listord
	void SetNppart (int _nppart) {nppart = _nppart;}; // Set nppart
	void SetIbegpartarr (int *_ibegpartarr) {ibegpartarr = _ibegpartarr;}; // Set ibegpartarr
	void SetIendpartarr (int *_iendpartarr) {iendpartarr = _iendpartarr;}; // Set iendpartarr
	void SetNhcells (int _nhcells) {nhcells = _nhcells;}; // Set nhcells
	void SetHcells (int *_hcells) {hcells = _hcells;}; // Set hcells
	void SetHcells2cpu (int *_hcells2cpu) {hcells2cpu = _hcells2cpu;}; // Set hcells2cpu
	void SetOrderpart (int *_orderpart) {orderpart = _orderpart;}; // Set orderpart
	void SetNlisttree (int _nlisttree) {nlisttree = _nlisttree;}; // Set nlisttree
	void SetListtree (int *_listtree) {listtree = _listtree;}; // Set listtree
	void SetInvordertree (int *_invordertree) {invordertree = _invordertree;}; // Set invordertree
	void SetAmtr (CSMatrixR *_amtr) {amtr = _amtr;}; // Set amtr
	void SetRhs (CSVector *_rhs) {rhs = _rhs;}; // Set rhs
	void SetX (CSVector *_x) {x = _x;}; // Set x
	void SetNorm (CSVector *_norm) {norm = _norm;}; // Set norm
	void SetPoly2 (CSVector *_poly2) {poly2 = _poly2;}; // Set poly2
	void SetSol (CSVector *_sol) {sol = _sol;}; // Set sol
	void SetGmtral (CGSMatrixR *_gmtral) {pgmtral = _gmtral;}; // Set gmtral
	void SetGmtrau (CGSMatrixR *_gmtrau) {pgmtrau = _gmtrau;}; // Set gmtrau
	void SetGmtra2 (CGSMatrixR *_gmtra2) {pgmtra2 = _gmtra2;}; // Set gmtra2
	void SetRhsnew (CSVector *_rhsnew) {prhsnew = _rhsnew;}; // Set rhsnew
	void SetXnew (CSVector *_xnew) {pxnew = _xnew;}; // Set xnew
	void SetNormnew (CSVector *_normnew) {pnormnew = _normnew;}; // Set normnew
	void SetPoly2new (CSVector *_poly2new) {ppoly2new = _poly2new;}; // Set poly2new
	void SetSolnew (CSVector *_solnew) {psolnew = _solnew;}; // Set solnew
private:
// Data
	CMPIComm *comm;   // comm is the MPI communicator
	int nprocsolve;   // nprocsolve is the size number of processors assumed to be used to solve the linear system
	int ntot;         // ntot is the size of the system of linear equations
	int nlistini;     // nlistini is the size of the system of linear equations
	int *iorderini;   // iorderini[nlistini] is the initial local sorting ordering of the rows/columns of the matrix for monotonicity
	int *listini;     // listini[nlistini] is the sorted initial list of rows of the matrix provided by user
	int *listord;     // listord[nlistini] is the globally ordered and partitioned initial list of rows of the matrix provided by user
	int nppart;       // nppart is the number of partitioning intervals
	int *ibegpartarr; // ibegpartarr[nppart] is the array of the first partitioning indices
	int *iendpartarr; // iendpartarr[nppart] is the array of the last  partitioning indices
	int nhcells;      // nhcells is the total number of hypercells
	int *hcells;      // hcells[nhcells+1] is the description of hypercells
	int *hcells2cpu;  // hcells2cpu[nhcells] is the CPU's partitioning of hypercells
	int *orderpart;   // orderpart[nlistpart] is the global ordering of the indices specified by the interval
	int nlisttree;    // nlisttree is the number of local rows on cpu according to the tree
	int *listtree;    // listtree[nlisttree] is the list of indices specified by the tree
	int *invordertree;// invordertree[nlisttree] is the global inverse ordering of the local indices specified by the tree
	CSMatrixR *amtr;  // amtr contains the matrix to be solved
	CSVector *rhs;    // rhs contains the right hand side
	CSVector *x;      // s contains the initial guess to the solution
	CSVector *norm;   // norm contains the normalization
	CSVector *poly2;  // poly2 contains the nonlinear coefficients
	CSVector *sol;    // sol contains the solution
	CTree *tree;      // tree contains the tree of processors
	CSlvParam *param; // param contains parameters of the solver
	CGSMatrixR *pgmtral; // pgmtral contains pointer to the transformed matrix data (L part)
	CGSMatrixR *pgmtrau; // pgmtrau contains pointer to the transformed matrix data (U part)
	CGSMatrixR *pgmtra2; // pgmtrau contains pointer to the transformed matrix data (second order)
	CSVector *prhsnew;   // rhsnew contains pointer to the transformed right hand side
	CSVector *pxnew;     // xnew contains pointer to the transformed initial guess
	CSVector *pnormnew;  // normnew contains pointer to the transformed norm data
	CSVector *ppoly2new; // poly2new contains pointer to the transformed poly2 data
	CSVector *psolnew;   // solnew contains pointer to the transformed solution
// Internal functions
};

#endif

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
