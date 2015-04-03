//------------------------------------------------------------------------------------------------
// File: ilu2grs2ind.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#ifdef __PMPITRACE__
#include "mpi.h"
#include "mpe.h"
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <ctime>

#include "globals.h"
#include "tree.h"
#include "fct.h"
#include "fctdiag.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "slvparam.h"
#include "mvm.h"
#include "ExchangeMPI.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CGSMatrixR: Second order Incomplete LU decomposition
//========================================================================================
void CGSMatrixR::Ilu2Index (ofstream &_fout, const CSlvParam &_param, // Second order Incomplete LU decomposition
							const CGSMatrixR &_mtral, const CGSMatrixR &_mtrau,
							CGSMatrixR &_mtrl, CGSMatrixR &_mtru) { 

	const char *funcname = "Ilu2Index";

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;

// Init time measurement

	time0 = clock ();

// Compute block sparsity structure of the matrix

	CSMatrix **pmtrarr;

	pmtrarr = new CSMatrix * [_mtral.nblksr];
	if (!pmtrarr) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		pmtrarr[iblk] = &_mtral.mtrarr[iblk];
	};

	int *iab, *jab;

	const CGSMatrix *pmtral;

	pmtral = &_mtral;

	pmtral->MtrBlockSparsity2Index (pmtrarr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsityOpt (_mtral.nblksr, iab, jab, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	int nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
						iabextt, jabextt);

// Allocate Fct data

	CFctR fct (_mtral.nblksr, 0, 1); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr); 

// Compute the scaling of the matrix

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		fctdiag.Ilu2ScaleBlkRow2Index (iblk, _param.sclpttype, _param.sclmin,
									_mtral, _mtrau.mtrarr[iblk], fct);
	};

//	csvhyst (cout, "ScDia",nloc, _dpiv);
//	csvhyst (_fout,"ScDia",nloc, _dpiv);

// Create dummy matrix

	CSMatrixR mtrdummy;

// Allocate array containing second order matrices

	CSMatrixR *mtrl2arr, *mtru2arr;

	mtrl2arr = new CSMatrixR [_mtral.nblksr];
	if (!mtrl2arr) MemoryFail (funcname);
	mtru2arr = new CSMatrixR [_mtral.nblksr];
	if (!mtru2arr) MemoryFail (funcname);

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

// Main cycle over the block rows

	CSMatrixR mtrltemp, mtrl2temp;
	CSMatrixR mtrutemp, mtru2temp;

	int j, jblk;

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {

// Load diagonal data into fct structure

		int nlstblk = iabext[iblk+1]-iabext[iblk];
		int ibs = iabext[iblk];
		int nlstblkcol = iabextt[iblk+1]-iabextt[iblk];
		int ibst = iabextt[iblk];

		fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, jabext+ibs, nlstblkcol, jabextt+ibst,
										_mtral, fct);

// Init current block row of the matrix

		fct.Ilu2InitBlkRow2Index (iblk, _param.dshift, 
									_mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
									mtrltemp, mtrutemp);

// Update current block row by the previous ones

		for (j=iabextt[iblk];j<iabextt[iblk+1]-1;j++) {

			jblk = jabextt[j];

//			_fout << " update Jblk = " << jblk << endl;

			fct.Ilu2UpdateBlkRow2Index (iblk, _param,
									_mtral, mtrltemp, mtrutemp,
									mtrl2arr[jblk], mtru2arr[jblk],
									mtrl2temp, mtru2temp);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

		};

// Factorize current block row

		fct.Ilu2FctBlkRow2Index (iblk, _param, 
								_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
								_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
								mtrl2arr[iblk], mtru2arr[iblk]);

		mtrltemp = mtrdummy;
		mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

		nlstblk = iabext[iblk+1]-iabext[iblk];
		ibs = iabext[iblk];

		fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, jabext+ibs,
												_mtral, fct);

	};

// Free second order matrices

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		mtrl2arr[iblk] = mtrdummy;
		mtru2arr[iblk] = mtrdummy;
	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Ilu2 factorization statistics

	int nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	int nzmema = 2 * nzatotal;

	nzmema -= _mtral.blksr[_mtral.nblksr];

	int nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	nzutot *= 2;

	int nzrtot = fct.GetNzr ();

	nzrtot *= 2;

	densu = 1.0e2 * (double) nzutot / (double) nzatotal;
	densr = 1.0e2 * (double) nzrtot / (double) nzatotal;

	ops = fct.GetOps ();

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	if (_param.msglev > 1) {
		cout  << " Ilu2 preconditioner generation statistics: " << endl;
		cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev > 0) {
		_fout << " Ilu2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Free working arrays

	delete [] pmtrarr;
	delete [] iab;
	delete [] jab;
	delete [] iabext;
	delete [] jabext;
	delete [] iabextt;
	delete [] jabextt;
	delete [] mtrl2arr;
	delete [] mtru2arr;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Compute second order Incomplete LU decomposition in parallel mode
//========================================================================================
void CGSMatrixR::Ilu2Schur2Index (ofstream &_fout, const CTree &_tree, CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
												CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
												CGSMatrixR &_mtrl, CGSMatrixR &_mtru) { 

	const char *funcname = "Ilu2Schur2Index";

	int nproc = _tree.GetNproc();
	int myid = _tree.GetMyid();

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fct_beg;
	static int fct_end;
	static int fctscl_beg;
	static int fctscl_end;

	if (i_set == 0) {

		fct_beg = MPE_Log_get_event_number ();
		fct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fct_beg, fct_end, "Fct", "grey");

		fctscl_beg = MPE_Log_get_event_number ();
		fctscl_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctscl_beg, fctscl_end, "FctScl", "magenta");

		i_set = 1;

	};

	MPE_Log_event (fct_beg, 0, NULL);
#endif

// Check the data on entry

	int numprocs;

	CMPIComm pcomm = _tree.GetComm ();
	numprocs = pcomm.GetNproc ();

	if (nproc != numprocs) {
		assert (false); throw " CGSMatrixR::Ilu2Schur2Index: Parallel Ilu2 function called with wrong number of processors ";
	};

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkIlu2Schur_",_tree.myid,".dat");
//	std::ofstream fout (strbuff);

// Statistics data

	double ops, densu, densr, denss, perf, tottim;
	clock_t time0, time1;
	int nzstot=0;

// Init time measurement

	time0 = clock ();

// Compute bl2cpu

	int *bl2cpu;
	int *listsub;

	bl2cpu = new int [_mtral.nblksr];
	if (!bl2cpu) MemoryFail (funcname);
	listsub = new int [_mtral.nblksr];
	if (!listsub) MemoryFail (funcname);

	_tree.Block2CpuSchur (bl2cpu);

// Compute block sparsity structure of the matrix

	CSMatrix **pmtrarr;

	pmtrarr = new CSMatrix * [_mtral.nblksr];
	if (!pmtrarr) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		pmtrarr[iblk] = &_mtral.mtrarr[iblk];
	};

	int *iab, *jab;

	const CGSMatrix *pmtral;

	pmtral = &_mtral;

	pmtral->MtrBlockSparsity2Index (true, _tree, bl2cpu, pmtrarr, iab, jab);

//	if (_tree.GetMyid () == 0) CSMatrix::PrintSparsity (1, "BlkSpMtr.ps", _mtral.nblksr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsityOpt (_mtral.nblksr, 
												iab, jab, iabext, jabext);

// Perform filtering of the extended sparsity according to the tree

	int *iabextflt, *jabextflt;

	int nzjext = iabext[_mtral.nblksr];

	iabextflt = new int [_mtral.nblksr+1];
	if (!iabextflt) MemoryFail (funcname);
	jabextflt = new int [nzjext];
	if (!jabextflt) MemoryFail (funcname);

	CSMatrix::SparsifyExtendedSchurBinaryTree (_tree,
																_mtral.nblksr, iabext, jabext,
																iabextflt, jabextflt);

	delete [] iabext;
	delete [] jabext;

	iabext = iabextflt;
	jabext = jabextflt;

//	if (_tree.GetMyid () == 0) CSMatrix::PrintSparsity (1, "ExtBlkSpMtr.ps", _mtral.nblksr, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
							iabextt, jabextt);

// Compute the list of blocks for which Schur data should be supported

	int nlevloc = _tree.GetNlev ();

	int nlistschur;
	int *listschur;
	int *lev2nodes;

	listschur = new int [_mtral.nblksr+1];
	if (!listschur) MemoryFail (funcname);
	lev2nodes = new int [nlevloc*nproc];
	if (!lev2nodes) MemoryFail (funcname);

	_tree.SchurList  (myid, nlistschur, listschur);
	_tree.LevelNodes (lev2nodes);

//	OutArr (_fout," List of Schur blocks ",nlistschur,listschur);
//	OutArr (_fout," LevNodes ",nlevloc*nproc,lev2nodes);
//	OutArr (_fout," Blk2Cpu ",_mtral.nblksr,bl2cpu);

//	_fout << " Tree = " << _tree << endl;

// Determine blamx

	int blamxl=1;

// Allocate Fct data

	CFctR fct (_mtral.nblksr, 0, blamxl); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr);

// Compute the local scaling of the matrix

#ifdef __PMPITRACE__
	MPE_Log_event (fctscl_beg, 0, NULL);
#endif

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {

		if (bl2cpu[iblk] == myid) {

			fctdiag.Ilu2ScaleBlkRow2Index (iblk, _param.sclpttype, _param.sclmin,
										_mtral, _mtrau.mtrarr[iblk], fct);

		};

	};

// Update and output hystogram

	pcomm = _tree.GetComm ();

	fctdiag.hyst.MakeGlobalHystogram (pcomm);
	fctdiag.hyst2.MakeGlobalHystogram (pcomm);

	if (myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " PointScaling" << fctdiag.hyst2;
			cout  << " Scaling" << fctdiag.hyst;
		};
		if (_param.msglev > 0) {
			_fout << " PointScaling" << fctdiag.hyst2;
			_fout << " Scaling" << fctdiag.hyst;
		};
	};

// Exchange the scaling data from the bordering

	fctdiag.ExchangeScalingGeneral (true, _tree, bl2cpu, nlistschur, listschur,
												iab, jab, _mtral);

#ifdef __PMPITRACE__
	MPE_Log_event (fctscl_end, 0, NULL);
#endif

// Create dummy matrix

	CSMatrixR mtrdummy;

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

	_mtrl.ivect = _param.sclpttype;
	_mtru.ivect = _param.sclpttype;

// Init gmtrl2 and gmtru2

	CGSMatrixR gmtrl2, gmtru2;
	CGSMatrixR gmtrl3, gmtru3;

	gmtrl2 = _mtral;
	gmtru2 = _mtrau;

	gmtrl3 = _mtral;
	gmtru3 = _mtrau;

// Init Schur data blocks

	int ilist;

	for (ilist=0;ilist<nlistschur;ilist++) {

		iblk = listschur[ilist];

		gmtrl2.mtrarr[iblk] = fct.CreateBlock2Index (iblk, _mtral);
		gmtru2.mtrarr[iblk] = fct.CreateBlock2Index (iblk, _mtral);

	};

// Main cycle over the block rows according to the tree

	CSMatrixR mtrltemp, mtrl2temp, mtrlrecv, mtrlsend;
	CSMatrixR mtrutemp, mtru2temp, mtrurecv, mtrusend;

	int ilev, inodecurr;

	for (ilev=nlevloc-1;ilev>=0;ilev--) {

		inodecurr = lev2nodes[myid*nlevloc+ilev];

// Exchange Schur complement and diagonal update data

		bool last_work_with_Schur = true;

		if (_tree.nodes[inodecurr].indtree >= 0) last_work_with_Schur = false;

		if (ilev < nlevloc-1) {

//			cout << " Before SchurUpdates2Index " << endl;

			SchurUpdates2Index (last_work_with_Schur,
										ilev, _tree, _tree, fct, fctdiag,
										bl2cpu,
										iabext, jabext,
										nlistschur, listschur,
										lev2nodes,
										gmtrl2, gmtru2,
										_param);

//			cout << " After  SchurUpdates2Index " << endl;

		};

// Perform computations with blocks of the current node taking into account Schur data

		bool update_by_Schur = true;

		if (ilev == nlevloc-1) update_by_Schur = false;

		if (_tree.nodes[inodecurr].indtree >= 0) {

			int indloc = _tree.nodes[inodecurr].indtree;

// Compute local array lev2nodes

			int nlevloc1  = _tree.subtreearr[indloc].GetNlev  ();
			int nprocloc1 = _tree.subtreearr[indloc].GetNproc ();

			int *lev2nodesloc;

			lev2nodesloc = new int [nlevloc1*nprocloc1];
			if (!lev2nodesloc) MemoryFail (funcname);

			_tree.subtreearr[indloc].LevelNodes (lev2nodesloc);

// Factorize the subtree

//			cout << " Before FctSubTree2Index " << endl;

			FctSubTree2Index (_tree, _tree.subtreearr[indloc],
									fct, fctdiag,
									bl2cpu, iabext, jabext, iabextt, jabextt,
									nlistschur, listschur, lev2nodesloc,
									_mtral, _mtrau,
									_mtrl, _mtru, gmtrl2, gmtru2,
									_param);

//			cout << " After  FctSubTree2Index " << endl;

			delete [] lev2nodesloc;

		} else {

			if (_tree.nodes[inodecurr].nodecpu == myid) {

//				cout << " Before FctNode2Index " << endl;

				FctNode2Index (update_by_Schur, inodecurr, _tree, _tree, 
									fct, fctdiag, 
									iabext, jabext, iabextt, jabextt,
									nlistschur, listschur,
									_mtral, _mtrau,
									_mtrl, _mtru, gmtrl2, gmtru2,
									_param);

//				cout << " After  FctNode2Index " << endl;

			};

		};

	};

// Store point/block scaling data in L and U structures

	if (_mtrl.ivect < 2) {
		fctdiag.StorePointScaling (_mtrl);
		fctdiag.StorePointScaling (_mtru);
	} else {
		fctdiag.StoreBlockScaling ('L', _mtrl);
		fctdiag.StoreBlockScaling ('U', _mtru);
	};

// Finalize time measurement

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute the global hystogram of diagonal pivots and output it

	pcomm = _tree.GetComm ();

	fct.hyst.MakeGlobalHystogram (pcomm);

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " Factorization" << fct.hyst;
		};
		if (_param.msglev > 0) {
			_fout << " Factorization" << fct.hyst;
		};
	};

// Output Ilu2 factorization statistics

	double nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	double nzmema = 2 * nzatotal;

	double nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	nzutot *= 2.0e0;

	double nzrtot = fct.nzr;

	nzrtot *= 2.0e0;

// Compute global data

	double arr0[7], arr1[7];

	arr0[0] = (double) nzutot;
	arr0[1] = (double) nzrtot;
	arr0[2] = (double) nzstot;
	arr0[3] = (double) nzatotal;
	arr0[4] = (double) nzmema;
	arr0[5] = fct.ops;
	arr0[6] = (double) fct.nmodsv;

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													7, arr0, arr1);
	} else {
		for (int i=0;i<=6;i++) arr1[i] = arr0[i];
	};

	double dnzu = arr1[0];
	double dnzr = arr1[1];
	double dnzs = arr1[2];
	double dnzamem = arr1[4];
	ops = arr1[5];
	int nmodtot = (int) arr1[6];

// Output resulting statistics

	densu = 1.0e2 * dnzu / dnzamem;
	densr = 1.0e2 * dnzr / dnzamem;
	denss = 1.0e2 * dnzs / dnzamem;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnzamem;

	if (myid == _tree.rootcpu) {

		if (_param.msglev > 1) {
			cout  << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				cout  << "     NmodS = " << nmodtot << endl;
			};
			cout  << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
			cout  << "     Costs = " << ops << " MvmA flops. " << endl;
			cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

		if (_param.msglev > 0) {
			_fout << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				_fout << "     NmodS = " << nmodtot << endl;
			};
			_fout << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
			_fout << "     Costs = " << ops << " MvmA flops. " << endl;
			_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

	};

// Free working arrays

	delete [] bl2cpu;
	delete [] listsub;
	delete [] pmtrarr;
	delete [] iab;
	delete [] jab;
	delete [] iabext;
	delete [] jabext;
	delete [] iabextt;
	delete [] jabextt;

#ifdef __PMPITRACE__
	MPE_Log_event (fct_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform parallel computations with the current subtree
//========================================================================================
void CGSMatrixR::FctSubTree2Index (const CTree &_treeini, const CTree &_tree, // Perform parallel computations with the current subtree
												CFctR &_fct, CFctDiagR &_fctdiag,
												int *_bl2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
												int _nlistschur, int *_listschur, int *_lev2nodes,
												CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
												CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
												CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
												CSlvParam &_param) {

//	const char *funcname = "FctSubTree2Index";

// Open the tree

	int nlevloc = _tree.GetNlev ();
	int myidloc = _tree.GetMyid ();

// Main cycle over the block rows according to the tree

	CSMatrixR mtrltemp, mtrl2temp, mtrlrecv, mtrlsend;
	CSMatrixR mtrutemp, mtru2temp, mtrurecv, mtrusend;

	int ilev, inodecurr;

	for (ilev=nlevloc-1;ilev>=0;ilev--) {

		inodecurr = _lev2nodes[myidloc*nlevloc+ilev];

// Exchange Schur complement and diagonal update data

		if (ilev < nlevloc-1) {

//			cout << " In subtree: Before SchurUpdates2Index " << endl;

			SchurUpdates2Index (true, 
										ilev, _treeini, _tree, _fct, _fctdiag,
										_bl2cpu,
										_iabext, _jabext,
										_nlistschur, _listschur,
										_lev2nodes,
										_mtrlschr, _mtruschr,
										_param);

//			cout << " In subtree: After  SchurUpdates2Index " << endl;

		};

// Perform computations with blocks of the current node taking into account Schur data

		bool update_by_Schur = true;

//		if (ilev == nlevloc-1) update_by_Schur = false;

		if (_tree.nodes[inodecurr].nodecpu == myidloc) {

//			cout << " In subtree: Before FctNode2Index " << endl;

			FctNode2Index (update_by_Schur, inodecurr, _treeini, _tree, 
								_fct, _fctdiag, 
								_iabext, _jabext, _iabextt, _jabextt,
								_nlistschur, _listschur,
								_mtral, _mtrau,
								_mtrl, _mtru, _mtrlschr, _mtruschr,
								_param);

//			cout << " In subtree: After  FctNode2Index " << endl;

		};

	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform Schur updates via exchanges from all involved processors
//========================================================================================
void CGSMatrixR::SchurUpdates2Index (bool _last_work_with_Schur, // Perform Schur updates via exchanges from all involved processors
													int _ilev, const CTree &_treeini, const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag,
													int *_blk2cpu,
													int *_iabextflt, int *_jabextflt,
													int _nlistschur, int *_listschur,
													int *_lev2nodearr,
													CGSMatrixR &_mtralsch, CGSMatrixR &_mtrausch,
													CSlvParam &_param) {

	const char *funcname = "SchurUpdates2Index";

// Create the mask of schur blocks

	int nblksloc = _mtralsch.GetNblksc ();

	int *imaskschur;

	imaskschur = new int [nblksloc];
	if (!imaskschur) MemoryFail (funcname);

	int i, iblk;

	for (i=0;i<nblksloc;i++) imaskschur[i] = -1;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		imaskschur[iblk] = i;
	};

// Create the list of blocks for exchanges

	int myid = _treeini.GetMyid ();
	int myidloc = _tree.GetMyid ();
	int nlevloc = _tree.GetNlev ();

	int inode = _lev2nodearr[nlevloc*myidloc+_ilev];

	int *imask, *listloc;

	imask = new int [nblksloc];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nblksloc];
	if (!listloc) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nblksloc;i++) imask[i] = -1;

	int nprocloc = _tree.nodes[inode].nprocnode;
	int *pcpuidarrnode = _tree.nodes[inode].cpuidarrnode;
	int *pcpuidarrnodeloc = _tree.nodes[inode].cpuidarrnodeloc;

	icycle++;

	int nlistloc = 0;

	int inodeloc, iproc, iprocloc, j;

	for (i=0;i<nprocloc;i++) {
		iproc = pcpuidarrnode[i];
		iprocloc = pcpuidarrnodeloc[i];
		inodeloc = _lev2nodearr[nlevloc*iprocloc+_ilev];
		for (j=_tree.nodes[inodeloc].indbeg;j<=_tree.nodes[inodeloc].indend;j++) {
			if (imask[j] != icycle) {
				listloc[nlistloc] = j;
				nlistloc++;
				imask[j] = icycle;
			};
		};
	};

	qsort (listloc, nlistloc, sizeof(int), compint);

// Prepare send lists

	int *iasend, *iptr;

	iasend = new int [nprocloc+1];
	if (!iasend) MemoryFail (funcname);
	iptr = new int [nprocloc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<=nprocloc;i++) iasend[i] = 0;

	int jblk;

	for (i=0;i<nprocloc;i++) {
		iproc = pcpuidarrnode[i];
		for (j=0;j<nlistloc;j++) {
			jblk = listloc[j];
			if (imaskschur[jblk] >= 0 && _blk2cpu[jblk] == iproc && iproc != myid) iasend[i+1]++;
		};
	};

	for (i=0;i<nprocloc;i++) iasend[i+1] = iasend[i]+iasend[i+1];
	for (i=0;i<nprocloc;i++) iptr[i] = iasend[i];

	int nzjasend = iasend[nprocloc];

	int *jasend;

	jasend = new int [nzjasend];
	if (!jasend) MemoryFail (funcname);

	int k;

	for (i=0;i<nprocloc;i++) {
		iproc = pcpuidarrnode[i];
		for (j=0;j<nlistloc;j++) {
			jblk = listloc[j];
			if (imaskschur[jblk] >= 0 && _blk2cpu[jblk] == iproc && iproc != myid) {
				k = iptr[i];
				jasend[k] = jblk;
				iptr[i]++;
			};
		};
	};

//	OutArr (cout," Iasend =",nprocloc+1,iasend);
//	OutArr (cout," Jasend =",nzjasend,jasend);

// Create send structures and free Schur blocks to be sent

	CSMatrixR mtrdummy;

	CSMatrixR mtrlsend, mtrusend;

	char **charlusend;
	int *charlen;

	charlusend = new char * [3*nprocloc];
	if (!charlusend) MemoryFail (funcname);
	charlen = new int [3*nprocloc];
	if (!charlen) MemoryFail (funcname);

	double thetaloc = _param.theta / 2.0e0;

	int NObjSend = 0;

	for (i=0;i<nprocloc;i++) {

//		char strbufff11[256];
//		sprintf (strbufff11, "%s%i%s%i%s%i%s","BefFiltrSnd_",_ilev,"_ind_",i,"_cpu_",_treeini.myid,".dat");
//		ofstream fffout (strbufff11,ios::app);
//		for (j=iasend[i];j<iasend[i+1];j++) {
//			jblk = jasend[j];
//			fffout << " Jblk = " << jblk << " L = " << _mtralsch.mtrarr[jblk] << " U = " << _mtrausch.mtrarr[jblk] << endl;
//		};

		_mtralsch.FilterListOfSchurBlocks2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																_iabextflt, _jabextflt, thetaloc,
																_fct);
		_mtrausch.FilterListOfSchurBlocks2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																_iabextflt, _jabextflt, thetaloc,
																_fct);

//		for (j=iasend[i];j<iasend[i+1];j++) {
//			jblk = jasend[j];
//			fffout << " Filtered Jblk = " << jblk << " L = " << _mtralsch.mtrarr[jblk] << " U = " << _mtrausch.mtrarr[jblk] << endl;
//		};
		iproc = pcpuidarrnode[i];
		if (iproc != myid && iasend[i+1] > iasend[i]) {
			_mtralsch.CombineMatrices2Index (true, iasend[i+1]-iasend[i], jasend+iasend[i], 
														_mtralsch.mtrarr, mtrlsend);
//			fffout << " LSend = " << mtrlsend << endl;
			for (j=iasend[i];j<iasend[i+1];j++) {
				jblk = jasend[j];
				_mtralsch.mtrarr[jblk] = mtrdummy;
			};
			if (!_last_work_with_Schur) {
				for (j=iasend[i];j<iasend[i+1];j++) {
					jblk = jasend[j];
					_mtralsch.mtrarr[jblk] = _fct.CreateBlock2Index (jblk, _mtralsch);
				};
			};
			_mtrausch.CombineMatrices2Index (true, iasend[i+1]-iasend[i], jasend+iasend[i], 
														_mtrausch.mtrarr, mtrusend);
//			fffout << " USend = " << mtrusend << endl;
			for (j=iasend[i];j<iasend[i+1];j++) {
				jblk = jasend[j];
				_mtrausch.mtrarr[jblk] = mtrdummy;
			};
			if (!_last_work_with_Schur) {
				for (j=iasend[i];j<iasend[i+1];j++) {
					jblk = jasend[j];
					_mtrausch.mtrarr[jblk] = _fct.CreateBlock2Index (jblk, _mtralsch);
				};
			};
			CSMatrixR::SuperSparsifySymm2Index (true, mtrlsend, mtrusend);
//			fffout << " After sparsify LSend = " << mtrlsend << " USend = " << mtrusend << endl;

//			char strbufff[256];
//			sprintf (strbufff, "%s%i%s%i%s","UpdateSchurSndL_",i,"_cpu_",myid,".dat");
//			ofstream ffout (strbufff);
//			ffout << " L send = " << mtrlsend << endl;
//			char strbufff1[256];
//			sprintf (strbufff1, "%s%i%s%i%s","UpdateSchurSndU_",i,"_cpu_",myid,".dat");
//			ofstream ffout1 (strbufff1);
//			ffout1 << " U send = " << mtrusend << endl;

			mtrlsend.PackMatrix (charlen[3*i],charlusend[3*i]);
			mtrlsend = mtrdummy;
			mtrusend.PackMatrix (charlen[3*i+1],charlusend[3*i+1]);
			mtrusend = mtrdummy;
			NObjSend += 2;
		} else {
			charlusend[3*i] = new char [0];
			charlusend[3*i+1] = new char [0];
			charlen[3*i  ] = 0;
			charlen[3*i+1] = 0;
		};
	};
	for (i=0;i<nprocloc;i++) {
		iproc = pcpuidarrnode[i];
		if (iproc != myid && iasend[i+1] > iasend[i]) {
			mtrlsend = _mtralsch.CreateDiagonalUpdate2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																				_fctdiag);
//			char strbufff1[256];
//			sprintf (strbufff1, "%s%i%s%i%s%i%s","UpdateSchurSndDia_",i,"_lev_",_ilev,"_cpu_",myid,".dat");
//			ofstream ffout1 (strbufff1);
//			ffout1 << " D send = " << mtrlsend << endl;
			mtrlsend.PackMatrix (charlen[3*i+2],charlusend[3*i+2]);
			mtrlsend = mtrdummy;
			NObjSend++;
		} else {
			charlusend[3*i+2] = new char [0];
			charlen[3*i+2] = 0;
		};
	};

// Send/receive

	bool use_local_comm = true;

	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	NObjSend = 0;

	for (i=0;i<nprocloc;i++) {
		iproc = pcpuidarrnode[i];
		if (iproc != myid && iasend[i+1] > iasend[i]) {
			ObjTypeSend[NObjSend] = 1;
			ObjIDSend[NObjSend] = i;
			if (!use_local_comm) {
				CpuIDSend[NObjSend] = iproc;
			} else {
				CpuIDSend[NObjSend] = i;
			};
			ObjSizeSend[NObjSend] = charlen[3*i];
			ObjSend[NObjSend] = charlusend[3*i];
			NObjSend++;
			ObjTypeSend[NObjSend] = 2;
			ObjIDSend[NObjSend] = i;
			if (!use_local_comm) {
				CpuIDSend[NObjSend] = iproc;
			} else {
				CpuIDSend[NObjSend] = i;
			};
			ObjSizeSend[NObjSend] = charlen[3*i+1];
			ObjSend[NObjSend] = charlusend[3*i+1];
			NObjSend++;
			ObjTypeSend[NObjSend] = 3;
			ObjIDSend[NObjSend] = i;
			if (!use_local_comm) {
				CpuIDSend[NObjSend] = iproc;
			} else {
				CpuIDSend[NObjSend] = i;
			};
			ObjSizeSend[NObjSend] = charlen[3*i+2];
			ObjSend[NObjSend] = charlusend[3*i+2];
			NObjSend++;
		};
	};

// Exchange the data

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	CMPIComm commloc;

	if (!use_local_comm) {
		commloc = _tree.GetComm ();
	} else {
		commloc = _tree.nodes[inode].GetCommNode ();
	};

//	OutArr (cout," CpuIDSend = ",NObjSend,CpuIDSend);
//	OutArr (cout," ObjSizeSend = ",NObjSend,ObjSizeSend);

	info = CMPIExchange::DataExchangeMPI (commloc,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixR::SchurUpdates2Index: Error in DataExchangeMPI";

//	OutArr (cout," CpuIDRecv = ",NObjRecv,CpuIDRecv);
//	OutArr (cout," ObjSizeRecv = ",NObjRecv,ObjSizeRecv);

// Free send structures

	for (i=0;i<nprocloc*3;i++) delete [] charlusend[i];

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Update own Schur blocks and free parts of received data

	CSMatrixR *mtrarrloc;

	mtrarrloc = new CSMatrixR [nblksloc];
	if(!mtrarrloc) MemoryFail(funcname);

	CSMatrixR mtrrecv, mtrl2temp, mtru2temp;

	int itype;

	for (i=0;i<NObjRecv;i++) {
		mtrrecv.UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);
		delete [] ObjRecv[i];
		itype = ObjTypeRecv[i];

//		if (itype == 1) {
//		char strbufff[256];
//		sprintf (strbufff, "%s%i%s%i%s%i%s","UpdateSchurRcv_",i,"_lev_",_ilev,"_cpu_",myid,".dat");
//		ofstream ffout (strbufff,ios::app);
//		ffout << " Itype = " << itype << endl;
//		ffout << " Matr = " << mtrrecv << endl;
//		};

		if (itype == 1 || itype == 2) {
			_mtralsch.SplitMatrixList2Index (nlistloc, listloc,
												mtrrecv, 
												mtrarrloc);
//			OutArr (cout," Received blocks ",nlistloc,listloc);
			mtrrecv = mtrdummy;
			if (itype == 1) {
//				_mtralsch.AddReplaceSetOfBlocks2Index (_fct, nlistloc, listloc, _mtralsch.mtrarr, mtrarrloc);
				throw " Replaced interface ";
			} else {
//				_mtralsch.AddReplaceSetOfBlocks2Index (_fct, nlistloc, listloc, _mtrausch.mtrarr, mtrarrloc);
				throw " Replaced interface ";
			};
			for (j=0;j<nlistloc;j++) {
				iblk = listloc[j];
				mtrarrloc[iblk] = mtrdummy;
//				if (itype == 1) {
//					ffout << " Resulted Schur iblk = " << iblk << " L = " << _mtralsch.mtrarr[iblk] << endl;
//				} else {
//					ffout << " Resulted Schur iblk = " << iblk << " U = " << _mtrausch.mtrarr[iblk] << endl;
//				};
			};
		} else {
			_mtralsch.UpdateDiagonal2Index (mtrrecv, _fctdiag);
			mtrrecv = mtrdummy;
		};
	};

// Free receive structures

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Free work arrays

	delete [] imaskschur;
	delete [] imask;
	delete [] listloc;
	delete [] iasend;
	delete [] iptr;
	delete [] jasend;
	delete [] charlusend;
	delete [] charlen;
	delete [] mtrarrloc;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform computations with the current tree node including local Schur updates
//========================================================================================
void CGSMatrixR::FctNode2Index (bool _update_by_Schur, int _inode, const CTree &_treeini, const CTree &_tree, 
											CFctR &_fct, CFctDiagR &_fctdiag, // Perform computations with the current tree node including local Schur updates
											int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
											int _nlistschur, int *_listschur,
											CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
											CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
											CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
											CSlvParam &_param) {

	const char *funcname = "FctNode2Index";

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctdiag_beg;
	static int fctdiag_end;
	static int fctdiagupd_beg;
	static int fctdiagupd_end;
	static int fctdiagfct_beg;
	static int fctdiagfct_end;
	static int fctschr_beg;
	static int fctschr_end;
	static int fctschrblk_beg;
	static int fctschrblk_end;

	if (i_set == 0) {

		fctdiag_beg = MPE_Log_get_event_number ();
		fctdiag_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiag_beg, fctdiag_end, "FctDiag", "blue");

		fctdiagupd_beg = MPE_Log_get_event_number ();
		fctdiagupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagupd_beg, fctdiagupd_end, "FctDiagUpd", "gold");

		fctdiagfct_beg = MPE_Log_get_event_number ();
		fctdiagfct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagfct_beg, fctdiagfct_end, "FctDiagFct", "cyan");

		fctschr_beg = MPE_Log_get_event_number ();
		fctschr_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschr_beg, fctschr_end, "FctSchur", "white");

		fctschrblk_beg = MPE_Log_get_event_number ();
		fctschrblk_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblk_beg, fctschrblk_end, "FctSchurBlk", "purple");

		i_set = 1;

	};

#endif

// Compute block sparsity structure of the matrix

	CSMatrix **pmtrarr;

	pmtrarr = new CSMatrix * [_mtral.nblksr];
	if (!pmtrarr) MemoryFail (funcname);

// Perform computations with blocks of the current node taking into account Schur data

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_beg, 0, NULL);
#endif

	int iblkbeg = _tree.nodes[_inode].indbeg;
	int iblkend = _tree.nodes[_inode].indend;

	int nfrees_per_node = 5;
	int istepfree = (_tree.nodes[_inode].indend-_tree.nodes[_inode].indbeg+1) / nfrees_per_node;

	if (istepfree <= 0) istepfree = 1;

	int icount = -1;

	int iblk, j, jblk;

	CSMatrixR mtrltemp, mtrl2temp;
	CSMatrixR mtrutemp, mtru2temp;
	CSMatrixR mtrdummy;

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		icount++;

//		if (_param.msglev > 1) {
//			cout << " Myid = " << _treeini.myid << " Iblk = " << iblk << endl;
//		};
//		_fout << " Myid = " << _treeini.myid << " Iblk = " << iblk << endl;

// Load diagonal data into fct structure

		int nlstblk = _iabext[iblk+1]-_iabext[iblk];
		int ibs = _iabext[iblk];
		int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
		int ibst = _iabextt[iblk];

		_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
														_mtral, _fct);

// Init current block row of the matrix

		_fct.Ilu2InitBlkRow2Index (iblk, _param.dshift,
											_mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
											mtrltemp, mtrutemp);

// Update current block row by the Schur data

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagupd_beg, 0, NULL);
#endif

		if (_update_by_Schur) {

			mtrl2temp = _fct.AddBlocks2Index (mtrltemp, _mtrlschr.mtrarr[iblk]);
			mtru2temp = _fct.AddBlocks2Index (mtrutemp, _mtruschr.mtrarr[iblk]);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

			_mtrlschr.mtrarr[iblk] = mtrdummy;
			_mtruschr.mtrarr[iblk] = mtrdummy;

		};

// Update current block row by the previous ones

		for (j=_iabextt[iblk];j<_iabextt[iblk+1]-1;j++) {

			jblk = _jabextt[j];

			if (jblk >= _tree.nodes[_inode].indbeg) {

//				char strbufff[256];
//				sprintf (strbufff, "%s%i%s%i%s","UpdateLU_",iblk,"_",jblk,".dat");
//				ofstream ffout (strbufff);

//				ffout << " Iblk = " << iblk << endl;
//				ffout << " L bef fct = " << mtrltemp << " U bef fct = " << mtrutemp << endl;
//				ffout << " Jblk = " << jblk << endl;
//				ffout << " L upd = " << _mtrlschr.mtrarr[jblk] << " U upd = " << _mtruschr.mtrarr[jblk] << endl;

				_fct.Ilu2UpdateBlkRow2Index (iblk, _param,
														_mtral, mtrltemp, mtrutemp,
														_mtrlschr.mtrarr[jblk], _mtruschr.mtrarr[jblk],
														mtrl2temp, mtru2temp);

//				ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

				mtrltemp = mtrl2temp;
				mtrutemp = mtru2temp;

				mtrl2temp = mtrdummy;
				mtru2temp = mtrdummy;

			};

		};

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagupd_end, 0, NULL);
#endif

// Factorize current block row

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_beg, 0, NULL);
#endif

//		char strbuff[256];
//		sprintf (strbuff, "%s%i%s","LU_",iblk,".dat");
//		ofstream fout (strbuff);
//		fout << " Iblk = " << iblk << endl;
//		fout << " L bef fct = " << mtrltemp << " U bef fct = " << mtrutemp << endl;

		_fct.Ilu2FctBlkRow2Index (iblk, _param, 
											_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
											_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
											_mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk]);

//		fout << " L aft fct = " << _mtrl.mtrarr[iblk] << " U aft fct = " << _mtru.mtrarr[iblk] << endl;

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_end, 0, NULL);
#endif

		mtrltemp = mtrdummy;
		mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

		nlstblk = _iabext[iblk+1]-_iabext[iblk];
		ibs = _iabext[iblk];

		_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
															_mtral, _fct);

// Free second order matrices if necessary

		if (icount%istepfree == 0) {

			for (jblk=iblkbeg;jblk<=iblk;jblk++) {

				_mtrlschr.FilterBlock2Index  (jblk,iblk);
				_mtruschr.FilterBlock2Index  (jblk,iblk);

			};

		};

	};

// Scan second order matrices and filter unnecessary data

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		_mtrlschr.FilterBlock2Index  (iblk,iblkend);
		_mtruschr.FilterBlock2Index  (iblk,iblkend);

	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_end, 0, NULL);
#endif

// Compute the block sparsity structure for the whole bordering

#ifdef __PMPITRACE__
	MPE_Log_event (fctschr_beg, 0, NULL);
#endif

	int *iabbord, *jabbord;

	CGSMatrix *pmtral = &_mtral;

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		pmtrarr[iblk] = &_mtrlschr.mtrarr[iblk];
	};

	pmtral->MtrBlockSparsity2Index (pmtrarr, iblkbeg, iblkend, // Compute block sparsity of the matrix
												iabbord, jabbord);

// Compute transposed block sparsity structure

	int *iabbordtr, *jabbordtr;

	int nzjbord = iabbord[_mtral.nblksr];

	iabbordtr = new int [_mtral.nblksr+1];
	if (!iabbordtr) MemoryFail (funcname);
	jabbordtr = new int [nzjbord];
	if (!jabbordtr) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabbord, jabbord,
							iabbordtr, jabbordtr);

// Update Schur blocks according to the transposed block sparsity

	int ilist;

	for (ilist=0;ilist<_nlistschur;ilist++) {

		iblk = _listschur[ilist];

#ifdef __PMPITRACE__
MPE_Log_event (fctschrblk_beg, 0, NULL);
#endif

		if (iabbordtr[iblk+1] > iabbordtr[iblk]) {

// Load diagonal data

			int nlstblk = _iabext[iblk+1]-_iabext[iblk];
			int ibs = _iabext[iblk];
			int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
			int ibst = _iabextt[iblk];

			_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
															_mtral, _fct);

// Update

			for (j=iabbordtr[iblk];j<iabbordtr[iblk+1];j++) {

				jblk = jabbordtr[j];

//				char strbufff[256];
//				sprintf (strbufff, "%s%i%s%i%s%i%s","UpdateSchur_",iblk,"_",jblk,"_cpu_",_treeini.myid,".dat");
//				ofstream ffout (strbufff);
//				ffout << " Iblk = " << iblk << endl;
//				ffout << " L bef fct = " << _mtrlschr.mtrarr[iblk] << " U bef fct = " << _mtruschr.mtrarr[iblk] << endl;
//				ffout << " Jblk = " << jblk << endl;
//				ffout << " L upd = " << _mtrlschr.mtrarr[jblk] << " U upd = " <<_mtruschr.mtrarr[jblk] << endl;

				_fct.Ilu2UpdateBlkRow2Index (iblk, _param,
														_mtral, _mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk],
														_mtrlschr.mtrarr[jblk], _mtruschr.mtrarr[jblk],
														mtrl2temp, mtru2temp);

//				ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

				_mtrlschr.mtrarr[iblk] = mtrl2temp;
				_mtruschr.mtrarr[iblk] = mtru2temp;

				mtrl2temp = mtrdummy;
				mtru2temp = mtrdummy;

			};

// Unload diagonal data

			nlstblk = _iabext[iblk+1]-_iabext[iblk];
			ibs = _iabext[iblk];

			_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
																_mtral, _fct);

		};

#ifdef __PMPITRACE__
MPE_Log_event (fctschrblk_end, 0, NULL);
#endif

	};

// Free work arrays

	delete [] pmtrarr;
	delete [] iabbord;
	delete [] jabbord;
	delete [] iabbordtr;
	delete [] jabbordtr;

// Free second order matrices of the current node

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {
		_mtrlschr.mtrarr[iblk] = mtrdummy;
		_mtruschr.mtrarr[iblk] = mtrdummy;
	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctschr_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixR: Compute second order Incomplete LU decomposition in parallel mode (dynamic version)
//========================================================================================
void CGSMatrixR::Ilu2SchurDynamic2Index (ofstream &_fout, const CTree &_tree, CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode (dynamic version)
														CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
														CGSMatrixR &_mtrl, CGSMatrixR &_mtru) { 

	const char *funcname = "Ilu2SchurDynamic2Index";

	int nproc = _tree.GetNproc();
	int myid = _tree.GetMyid();

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fct_beg;
	static int fct_end;
	static int fctscl_beg;
	static int fctscl_end;

	if (i_set == 0) {

		fct_beg = MPE_Log_get_event_number ();
		fct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fct_beg, fct_end, "Fct", "grey");

		fctscl_beg = MPE_Log_get_event_number ();
		fctscl_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctscl_beg, fctscl_end, "FctScl", "magenta");

		i_set = 1;

	};

	MPE_Log_event (fct_beg, 0, NULL);
#endif

// Check the data on entry

	int numprocs;

	CMPIComm pcomm = _tree.GetComm ();
	numprocs = pcomm.GetNproc ();

	if (nproc != numprocs) {
		assert (false); throw " CGSMatrixR::Ilu2Schur2Index: Parallel Ilu2 function called with wrong number of processors ";
	};

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkIlu2Schur_",_tree.myid,".dat");
//	std::ofstream fout (strbuff);

// Statistics data

	double ops, densu, densr, denss, perf, tottim;
	clock_t time0, time1;
	int nzstot=0;

// Init time measurement

	time0 = clock ();

// Compute bl2cpu

	int *bl2cpu;
	int *listsub;

	bl2cpu = new int [_mtral.nblksr];
	if (!bl2cpu) MemoryFail (funcname);
	listsub = new int [_mtral.nblksr];
	if (!listsub) MemoryFail (funcname);

	_tree.Block2CpuSchur (bl2cpu);

// Compute block sparsity structure of the matrix

	CSMatrix **pmtrarr;

	pmtrarr = new CSMatrix * [_mtral.nblksr];
	if (!pmtrarr) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		pmtrarr[iblk] = &_mtral.mtrarr[iblk];
	};

	int *iab, *jab;

	const CGSMatrix *pmtral;

	pmtral = &_mtral;

	pmtral->MtrBlockSparsity2Index (true, _tree, bl2cpu, pmtrarr, iab, jab);

//	if (_tree.GetMyid () == 0) CSMatrix::PrintSparsity (1, "BlkSpMtr.ps", _mtral.nblksr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsityOpt (_mtral.nblksr, 
												iab, jab, iabext, jabext);

// Perform filtering of the extended sparsity according to the tree

	int *iabextflt, *jabextflt;

	int nzjext = iabext[_mtral.nblksr];

	iabextflt = new int [_mtral.nblksr+1];
	if (!iabextflt) MemoryFail (funcname);
	jabextflt = new int [nzjext];
	if (!jabextflt) MemoryFail (funcname);

	CSMatrix::SparsifyExtendedSchurBinaryTree (_tree,
																_mtral.nblksr, iabext, jabext,
																iabextflt, jabextflt);

	delete [] iabext;
	delete [] jabext;

	iabext = iabextflt;
	jabext = jabextflt;

//	if (_tree.GetMyid () == 0) CSMatrix::PrintSparsity (1, "ExtBlkSpMtr.ps", _mtral.nblksr, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
							iabextt, jabextt);

// Compute the list of blocks for which Schur data should be supported

	int nlistschur;
	int *listschur;

	listschur = new int [_mtral.nblksr+1];
	if (!listschur) MemoryFail (funcname);

	_tree.SchurList  (myid, nlistschur, listschur);

//	OutArr (_fout," Blks ",_mtral.nblksr,_mtral.blksr);
//	OutArr (_fout," Blk2Cpu ",_mtral.nblksr,bl2cpu);
//	OutArr (_fout," List of Schur blocks ",nlistschur,listschur);

//	_fout << " Tree = " << _tree << endl;

// Determine blamx

	int blamxl=1;

// Allocate Fct data

	CFctR fct (_mtral.nblksr, 0, blamxl); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr);

// Compute the local scaling of the matrix

#ifdef __PMPITRACE__
	MPE_Log_event (fctscl_beg, 0, NULL);
#endif

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {

		if (bl2cpu[iblk] == myid) {

			fctdiag.Ilu2ScaleBlkRow2Index (iblk, _param.sclpttype, _param.sclmin,
										_mtral, _mtrau.mtrarr[iblk], fct);

		};

	};

// Update and output hystogram

	pcomm = _tree.GetComm ();

	fctdiag.hyst.MakeGlobalHystogram (pcomm);
	fctdiag.hyst2.MakeGlobalHystogram (pcomm);

	if (myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " PointScaling" << fctdiag.hyst2;
			cout  << " Scaling" << fctdiag.hyst;
		};
		if (_param.msglev > 0) {
			_fout << " PointScaling" << fctdiag.hyst2;
			_fout << " Scaling" << fctdiag.hyst;
		};
	};

// Exchange the scaling data from the bordering

	fctdiag.ExchangeScalingGeneral (true, _tree, bl2cpu, nlistschur, listschur,
												iab, jab, _mtral);

#ifdef __PMPITRACE__
	MPE_Log_event (fctscl_end, 0, NULL);
#endif

// Create dummy matrix

	CSMatrixR mtrdummy;

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

	_mtrl.ivect = _param.sclpttype;
	_mtru.ivect = _param.sclpttype;

// Create structures that support global updates

	fct.AllocateGlobalArrays2Index (_mtral.nblksr, _mtral.blksr, nlistschur, listschur);

// Init gmtrl2 and gmtru2

	CGSMatrixR gmtrl2, gmtru2;
	CGSMatrixR gmtrl3, gmtru3;
	CGSMatrixR gmtrladd, gmtruadd;

	gmtrl2 = _mtral;
	gmtru2 = _mtrau;

	gmtrl3 = _mtral;
	gmtru3 = _mtrau;
	gmtrladd = _mtral;
	gmtruadd = _mtrau;

// Init Schur data blocks

	int ilist;

	for (ilist=0;ilist<nlistschur;ilist++) {

		iblk = listschur[ilist];

		gmtrl2.mtrarr[iblk] = fct.CreateBlock2Index (iblk, _mtral);
		gmtru2.mtrarr[iblk] = fct.CreateBlock2Index (iblk, _mtral);

	};

// Prepare dynamic computations

	int nblksloc = _mtral.nblksr;

	int n_lev_work, nsendmax;
	int *i_beg_work, *i_end_work, *i_beg_schur, *lev_work;
	int *imasksend, *imaskrecv;
	int *iptrextt, *updblkschur;
	char **psenddata, **precvdata;
	CMPIRequest *sendreqvarr, *recvreqvarr;
	CSMatrixR *mtrarrwork, *mtrarr2work;

	PrepareDynamicArrays2Index (_tree, 
											n_lev_work, lev_work,
											i_beg_work, i_end_work, i_beg_schur,
											nblksloc, iabextt, nlistschur, listschur, iptrextt, updblkschur,
											nsendmax, imasksend, psenddata, sendreqvarr, imaskrecv, precvdata, recvreqvarr, 
											mtrarrwork, mtrarr2work);

	int id_shift = n_lev_work + nblksloc;

//	_fout << " Tree = " << _tree << endl;
//	OutArr (_fout," lev_work ",3*n_lev_work,lev_work);
//	OutArr (_fout," i_beg_work ",n_lev_work,i_beg_work);
//	OutArr (_fout," i_end_work ",n_lev_work,i_end_work);
//	OutArr (_fout," i_beg_schur ",n_lev_work,i_beg_schur);
//	OutArr (_fout," listschur ",nlistschur,listschur);

// Main computational scheme

	int i_lev_work;
	int itree, inode, inodef, iproc, iproc1, nsends, isend, nrecvs, irecv;
	int ip_work_loc, ip_work, ip_work_next, ip_work_end;
	int remained_updates, remained_updates_save;
	bool is_completed, remained_work;

	int nsubtree = _tree.nsubtree;

	int j, jblk;

	nsends = 0;
	isend = 0;
	nrecvs = 0;
	irecv = 0;

	for (i_lev_work=0;i_lev_work<n_lev_work;i_lev_work++) {

//		_fout << " i_lev_work = " << i_lev_work << endl;
		itree = lev_work[i_lev_work*3];
		inode = lev_work[i_lev_work*3+1];
		inodef = lev_work[i_lev_work*3+2];
//		_fout << " itree = " << itree << " inode = " << inode << endl;

// Init variables and arrays that control priority work and receives

		if (i_lev_work == n_lev_work-1) {
			remained_updates = 0;
		} else {
			remained_updates = i_end_work[i_lev_work+1]-i_beg_work[i_lev_work+1]+1;

// Mark: log2 scanning the tree

			if (true) {
				if (itree < 0) {
					iproc = _tree.nodes[inode].nodecpu;
					if (iproc != myid) remained_updates = 0;
				} else {
					remained_updates = 0;
				};
			};
		};

//		_fout << " remained_receives = " << remained_receives << " remained_updates = " << remained_updates << endl;

		remained_updates_save = remained_updates;

		ip_work = i_beg_schur[i_lev_work];
		if (i_lev_work == 0) {
			ip_work_next = 0;
		} else if (i_lev_work == n_lev_work-1) {
			ip_work_next = nlistschur-1;
		} else {
			ip_work_next = i_beg_schur[i_lev_work+1];
		};
		ip_work = nlistschur;
		ip_work_next = nlistschur;
		ip_work_end = nlistschur-1;

//		_fout << " ip_work = " << ip_work << " ip_work_next = " << ip_work_next << " ip_work_end = " << ip_work_end << endl;

// Init async receives of the data from other cpu's

		InitReceiveBlocks2Index (i_lev_work, lev_work, _tree,
											nrecvs, imaskrecv, precvdata, recvreqvarr);

// Check previous sends for completion and free the memory

		CheckSends2Index (nsends, isend, imasksend, psenddata, sendreqvarr);

// Preliminary updates of the current block data

		ip_work_loc = ip_work;
/*
		while (ip_work_loc <= ip_work_next-1) {

				NonPriorityUpdates2Index (i_lev_work, 
													i_beg_work, i_end_work, 
													ip_work_loc, ip_work_next-1, listschur, updblkschur,
													_tree, fct, fctdiag, 
													bl2cpu, iabext, jabext, iabextt, iptrextt, jabextt,
													gmtrl2, gmtru2, _param);

		};
*/
		ip_work = ip_work_next;

// Main updates cycle including asynchronous receives and future updates while waiting

		is_completed = false;

		remained_work = true;

		while (!is_completed) {

			if (remained_updates == 0 && nrecvs == irecv) is_completed = true;
			if (remained_updates == 0 && ip_work > ip_work_end) remained_work = false;

			if (is_completed) goto NextCycle;

// Try to complete async sends or receives and perform updates

			if (nsends != isend || nrecvs != irecv) {

				TestSendsRecvsAndUpdate2Index (_tree, remained_work,
															id_shift, nsends, isend, 
															imasksend, psenddata, sendreqvarr,
															nrecvs, irecv,
															imaskrecv, precvdata, recvreqvarr,
															fct, fctdiag,
															gmtrl2, gmtru2, gmtrladd, gmtruadd, mtrarrwork, mtrarr2work);

			};

// Try to perform priority update if possible

			if (remained_updates != 0) {

				PriorityUpdates2Index (i_lev_work, 
												i_beg_work, i_end_work, i_beg_work[i_lev_work], 
												remained_updates, 
												_tree, fct, fctdiag, 
												bl2cpu, iabext, jabext, iabextt, iptrextt, jabextt,
												gmtrl2, gmtru2, _param);

				goto NextCycle;

			};

// Perform not prioritative work if possible

/*
			if (ip_work <= ip_work_end) {

				NonPriorityUpdates2Index (i_lev_work, 
													i_beg_work, i_end_work,
													ip_work, ip_work_end, listschur, updblkschur,
													_tree, fct, fctdiag, 
													bl2cpu, iabext, jabext, iabextt, iptrextt, jabextt,
													gmtrl2, gmtru2, _param);

			};
*/

NextCycle:;

			if (remained_updates == 0 && ip_work > ip_work_end) remained_work = false;

		};

// Check previous sends for completion and free the memory

		CheckSends2Index (nsends, isend, imasksend, psenddata, sendreqvarr);

// Factorize current node data if necessary

		iproc = -1;
		if (itree < 0) {
			if (_tree.nodes[inode].indtree < 0) {
				iproc = _tree.nodes[inode].nodecpu;
			};
		} else {
			iproc1 = _tree.subtreearr[itree].nodes[inode].nodecpu;
			iproc = _tree.nodes[inodef].cpuidarrnode[iproc1];
		};

		if (i_lev_work == 0 || (itree < 0 && iproc == myid && nsubtree == 0) || (itree >= 0 && iproc == myid)) {

			FctNode2Index (i_lev_work,
								i_beg_work, i_end_work,
								_tree, fct, fctdiag,
								iabext, jabext, iabextt, jabextt,
								_mtral, _mtrau,
								_mtrl, _mtru, gmtrl2, gmtru2, gmtrladd, gmtruadd,
								_param);

			if (i_lev_work < n_lev_work-1) {

				remained_updates = remained_updates_save;

				while (remained_updates != 0) {

					PriorityUpdates2Index (i_lev_work, 
													i_beg_work, i_end_work, i_beg_work[i_lev_work+1], 
													remained_updates, 
													_tree, fct, fctdiag, 
													bl2cpu, iabext, jabext, iabextt, iptrextt, jabextt,
													gmtrl2, gmtru2, _param);

				};
			};

		};

// Prepare and send node data if necessary

		if (i_lev_work < n_lev_work-1) {

			AsyncSendBlocks2Index (i_lev_work, lev_work,
											i_beg_work, i_end_work, nlistschur, listschur,
											_tree, fct, fctdiag,
											id_shift, nsendmax, nsends, imasksend, psenddata, sendreqvarr,
											bl2cpu, iabextflt, jabextflt, iabextt, iptrextt, jabextt,
											gmtrl2, gmtru2,
											_param);

		};

// Check previous sends for completion and free the memory

		CheckSends2Index (nsends, isend, imasksend, psenddata, sendreqvarr);

// Free parts of Schur blocks that are not in use any more

		if (i_lev_work != 0 && ((itree < 0 && nsubtree == 0) || (itree >= 0))) {
			int iblkbeg = i_beg_work[i_lev_work]-1;
			for (j=0;j<nlistschur;j++) {
				jblk = listschur[j];
				if (jblk <= iblkbeg) {
//					gmtrl2.FilterBlock2Index  (jblk,iblkbeg);
//					gmtru2.FilterBlock2Index  (jblk,iblkbeg);
//					fct.ReInitGlobalArraysByBlockData2Index (gmtrl2.mtrarr[jblk]);
				};
			};
		};

	};

// Wait sends for completion and free the memory

	WaitSends2Index (nsends, isend, imasksend, psenddata, sendreqvarr);

// Store point/block scaling data in L and U structures

	if (_mtrl.ivect < 2) {
		fctdiag.StorePointScaling (_mtrl);
		fctdiag.StorePointScaling (_mtru);
	} else {
		fctdiag.StoreBlockScaling ('L', _mtrl);
		fctdiag.StoreBlockScaling ('U', _mtru);
	};

// Free structures that support global updates

	fct.FreeGlobalArrays2Index ();

// Finalize time measurement

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Compute the global hystogram of diagonal pivots and output it

	pcomm = _tree.GetComm ();

	fct.hyst.MakeGlobalHystogram (pcomm);

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout  << " Factorization" << fct.hyst;
		};
		if (_param.msglev > 0) {
			_fout << " Factorization" << fct.hyst;
		};
	};

// Output Ilu2 factorization statistics

	double nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	double nzmema = 2 * nzatotal;

	double nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	nzutot *= 2.0e0;

	double nzrtot = fct.nzr;

	nzrtot *= 2.0e0;

// Compute global data

	double arr0[7], arr1[7];

	arr0[0] = (double) nzutot;
	arr0[1] = (double) nzrtot;
	arr0[2] = (double) nzstot;
	arr0[3] = (double) nzatotal;
	arr0[4] = (double) nzmema;
	arr0[5] = fct.ops;
	arr0[6] = (double) fct.nmodsv;

	if (nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													7, arr0, arr1);
	} else {
		for (int i=0;i<=6;i++) arr1[i] = arr0[i];
	};

	double dnzu = arr1[0];
	double dnzr = arr1[1];
	double dnzs = arr1[2];
	double dnzamem = arr1[4];
	ops = arr1[5];
	int nmodtot = (int) arr1[6];

// Output resulting statistics

	densu = 1.0e2 * dnzu / dnzamem;
	densr = 1.0e2 * dnzr / dnzamem;
	denss = 1.0e2 * dnzs / dnzamem;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnzamem;

	if (myid == _tree.rootcpu) {

		if (_param.msglev > 1) {
			cout  << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				cout  << "     NmodS = " << nmodtot << endl;
			};
			cout  << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
			cout  << "     Costs = " << ops << " MvmA flops. " << endl;
			cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

		if (_param.msglev > 0) {
			_fout << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				_fout << "     NmodS = " << nmodtot << endl;
			};
			_fout << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
			_fout << "     Costs = " << ops << " MvmA flops. " << endl;
			_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

	};

// Free working arrays

	delete [] bl2cpu;
	delete [] listsub;
	delete [] pmtrarr;
	delete [] iab;
	delete [] jab;
	delete [] iabext;
	delete [] jabext;
	delete [] iabextt;
	delete [] jabextt;
	delete [] listschur;
	delete [] lev_work;
	delete [] i_beg_work;
	delete [] i_end_work;
	delete [] i_beg_schur;
	delete [] iptrextt;
	delete [] updblkschur;
	delete [] imasksend;
	delete [] imaskrecv;
	delete [] psenddata;
	delete [] sendreqvarr;
	delete [] precvdata;
	delete [] recvreqvarr;
	delete [] mtrarrwork;
	delete [] mtrarr2work;

#ifdef __PMPITRACE__
	MPE_Log_event (fct_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixR: Prepare all arrays for dynamic parallel preconditioner computation
//========================================================================================
void CGSMatrixR::PrepareDynamicArrays2Index (const CTree &_tree, // Prepare all arrays for dynamic parallel preconditioner computation
															int &_n_lev_work, int *&_lev_work, 
															int *&_i_beg_work, int *&_i_end_work, int *&_i_beg_schur,
															int _nblks, int *_iabextt, int _nlistschur, int *_listschur,
															int *&_iptrextt, int *&_updblkschur,
															int &_nsendmax, int *&_imasksend, char **&_psenddata, CMPIRequest *&_sendreqvarr, 
															int *&_imaskrecv, char **&_precvdata, CMPIRequest *&_recvreqvarr, 
															CSMatrixR *&_mtrarrwork, CSMatrixR *&_mtrarr2work) {

	const char *funcname = "PrepareDynamicArrays2Index";

// Scan the tree and count the total number of levels

	_n_lev_work = 0;

	int myid = _tree.myid;
	int inode = _tree.cpuidend[myid];

	int itree, inode1, myidloc, ifather, ifather1;

	while (inode != -1) {
		_n_lev_work++;
		itree = _tree.nodes[inode].indtree;
		if (itree >= 0) {
			myidloc = _tree.subtreearr[itree].myid;
			inode1 = _tree.subtreearr[itree].cpuidend[myidloc];
			while (inode1 != -1) {
				_n_lev_work++;
				ifather1 = _tree.subtreearr[itree].nodes[inode1].fatherid;
				if (ifather1 != inode1) {
					inode1 = ifather1;
				} else {
					inode1 = -1;
				};
			};
		};
		ifather = _tree.nodes[inode].fatherid;
		if (ifather != inode) {
			inode = ifather;
		} else {
			inode = -1;
		};
	};

// Scan the tree once again and fill control array

	_lev_work = new int [_n_lev_work*3];
	if (!_lev_work) MemoryFail (funcname);

	_n_lev_work = 0;

	inode = _tree.cpuidend[myid];

	while (inode != -1) {
		_lev_work[_n_lev_work*3] = -1;
		_lev_work[_n_lev_work*3+1] = inode;
		_lev_work[_n_lev_work*3+2] = -1;
		_n_lev_work++;
		itree = _tree.nodes[inode].indtree;
		if (itree >= 0) {
			myidloc = _tree.subtreearr[itree].myid;
			inode1 = _tree.subtreearr[itree].cpuidend[myidloc];
			while (inode1 != -1) {
				_lev_work[_n_lev_work*3] = itree;
				_lev_work[_n_lev_work*3+1] = inode1;
				_lev_work[_n_lev_work*3+2] = inode;
				_n_lev_work++;
				ifather1 = _tree.subtreearr[itree].nodes[inode1].fatherid;
				if (ifather1 != inode1) {
					inode1 = ifather1;
				} else {
					inode1 = -1;
				};
			};
		};
		ifather = _tree.nodes[inode].fatherid;
		if (ifather != inode) {
			inode = ifather;
		} else {
			inode = -1;
		};
	};

// Fill other arrays

	_i_beg_work = new int [_n_lev_work];
	if (!_i_beg_work) MemoryFail (funcname);
	_i_end_work = new int [_n_lev_work];
	if (!_i_end_work) MemoryFail (funcname);

	_nsendmax = 0;

	int i;

	for (i=0;i<_n_lev_work;i++) {
		itree = _lev_work[i*3];
		inode = _lev_work[i*3+1];
		if (itree < 0) {
			_i_beg_work[i] = _tree.nodes[inode].indbeg;
			_i_end_work[i] = _tree.nodes[inode].indend;
			if (i>0) _nsendmax += (_tree.nodes[inode].nprocnode-1);
		} else {
			_i_beg_work[i] = _tree.subtreearr[itree].nodes[inode].indbeg;
			_i_end_work[i] = _tree.subtreearr[itree].nodes[inode].indend;
			if (i>0) _nsendmax += (_tree.subtreearr[itree].nodes[inode].nprocnode-1);
		};
	};

	_nsendmax *= 2;

	_i_beg_schur = new int [_n_lev_work];
	if (!_i_beg_schur) MemoryFail (funcname);

	int iendblk, j;

	for (i=0;i<_n_lev_work;i++) {
		iendblk = _i_end_work[i];
		_i_beg_schur[i] = _nlistschur;
		for (j=_nlistschur-1;j>=0;j--) {
			if (_listschur[j] > iendblk) _i_beg_schur[i] = j;
		};
	};

// Allocate and init structures that support forward computations

	_iptrextt = new int [_nblks];
	if (!_iptrextt) MemoryFail (funcname);
	_updblkschur = new int [_nlistschur];
	if (!_updblkschur) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) _iptrextt[i] = _iabextt[i];
	for (i=0;i<_nlistschur;i++) _updblkschur[i] = 0;

// Allocate and init send/receive data

	_imasksend = new int [_nsendmax];
	if (!_imasksend) MemoryFail (funcname);
	_psenddata = new char * [_nsendmax];
	if (!_psenddata) MemoryFail (funcname);
	_sendreqvarr = new CMPIRequest [_nsendmax];
	if (!_sendreqvarr) MemoryFail (funcname);
	_imaskrecv = new int [_nsendmax+1];
	if (!_imaskrecv) MemoryFail (funcname);
	_precvdata = new char * [_nsendmax];
	if (!_precvdata) MemoryFail (funcname);
	_recvreqvarr = new CMPIRequest [_nsendmax];
	if (!_recvreqvarr) MemoryFail (funcname);

	for (i=0;i<_nsendmax;i++) _imasksend[i] = 0;
	for (i=0;i<_nsendmax;i++) _imaskrecv[i] = 0;

// Allocate computations support array

	_mtrarrwork = new CSMatrixR [_nblks];
	if (!_mtrarrwork) MemoryFail (funcname);
	_mtrarr2work = new CSMatrixR [_nblks];
	if (!_mtrarr2work) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform computations with the current tree node including local Schur updates
//========================================================================================
void CGSMatrixR::FctNode2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
											int *_i_beg_work, int *_i_end_work,
											const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag,
											int *_iabext, int *_jabext, int *_iabextt, int *_jabextt,
											CGSMatrixR &_mtral, CGSMatrixR &_mtrau,
											CGSMatrixR &_mtrl, CGSMatrixR &_mtru,
											CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr, CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd,
											CSlvParam &_param) {

	const char *funcname = "FctNode2Index";

//	cout << " FctNode2Index: begin " << endl;

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctdiag_beg;
	static int fctdiag_end;
	static int fctdiagupd_beg;
	static int fctdiagupd_end;
	static int fctdiagfct_beg;
	static int fctdiagfct_end;

	if (i_set == 0) {

		fctdiag_beg = MPE_Log_get_event_number ();
		fctdiag_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiag_beg, fctdiag_end, "FctDiag", "blue");

		fctdiagupd_beg = MPE_Log_get_event_number ();
		fctdiagupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagupd_beg, fctdiagupd_end, "FctDiagUpd", "gold");

		fctdiagfct_beg = MPE_Log_get_event_number ();
		fctdiagfct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagfct_beg, fctdiagfct_end, "FctDiagFct", "cyan");

		i_set = 1;

	};

#endif

// Perform computations with blocks of the current node taking into account Schur data if necessary

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_beg, 0, NULL);
#endif

	int iblkbeg = _i_beg_work[_i_lev_work];
	int iblkend = _i_end_work[_i_lev_work];

	int nfrees_per_node = 5;
	int istepfree = (iblkend-iblkbeg+1) / nfrees_per_node;

	if (istepfree <= 0) istepfree = 1;

	int icount = -1;

	int iblk, j, jblk;

	CSMatrixR mtrltemp, mtrl2temp;
	CSMatrixR mtrutemp, mtru2temp;
	CSMatrixR mtrdummy;

	int *listupd;

	listupd = new int [_mtral.nblksr];
	if (!listupd) MemoryFail (funcname);

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		icount++;

		if (_param.msglev > 1) {
			cout << " Myid = " << _tree.myid << " Iblk = " << iblk << endl;
		};
//		_fout << " Myid = " << _treeini.myid << " Iblk = " << iblk << endl;

// Load diagonal data into fct structure

		int nlstblk = _iabext[iblk+1]-_iabext[iblk];
		int ibs = _iabext[iblk];
		int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
		int ibst = _iabextt[iblk];

		_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
														_mtral, _fct);

// Init current block row of the matrix

		_fct.Ilu2InitBlkRow2Index (iblk, _param.dshift,
											_mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
											mtrltemp, mtrutemp);

// Update current block row by the Schur data

#ifdef __PMPITRACE__
		MPE_Log_event (fctdiagupd_beg, 0, NULL);
#endif

		if (_mtrladd.mtrarr[iblk].GetNzja () != 0) {
			if (_mtrlschr.mtrarr[iblk].GetNzja () != 0) {
				CSMatrixR::DoubleAddMatrSS2Index (true,
																_mtrlschr.mtrarr[iblk], _mtrladd.mtrarr[iblk], mtrl2temp,
																_mtruschr.mtrarr[iblk], _mtruadd.mtrarr[iblk], mtru2temp);
				_mtrlschr.mtrarr[iblk] = mtrl2temp;
				_mtruschr.mtrarr[iblk] = mtru2temp;
			} else {
				_mtrlschr.mtrarr[iblk] = _mtrladd.mtrarr[iblk];
				_mtruschr.mtrarr[iblk] = _mtruadd.mtrarr[iblk];
			};
			_mtrladd.mtrarr[iblk] = mtrdummy;
			_mtruadd.mtrarr[iblk] = mtrdummy;
		};

		if (_i_lev_work > 0) {

			mtrl2temp = _fct.AddBlocks2Index (mtrltemp, _mtrlschr.mtrarr[iblk]);
			mtru2temp = _fct.AddBlocks2Index (mtrutemp, _mtruschr.mtrarr[iblk]);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

			_mtrlschr.mtrarr[iblk] = mtrdummy;
			_mtruschr.mtrarr[iblk] = mtrdummy;

		};

// Create the list and update

		int nlistupd = 0;

		for (j=_iabextt[iblk];j<_iabextt[iblk+1]-1;j++) {

			jblk = _jabextt[j];

			if (jblk >= iblkbeg) {

				listupd[nlistupd] = jblk;
				nlistupd++;

			};

		};

		_fct.Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, _param,
												_mtral, mtrltemp, mtrutemp,
												nlistupd, listupd, _mtrlschr.mtrarr, _mtruschr.mtrarr,
												mtrl2temp, mtru2temp);

		mtrltemp = mtrl2temp;
		mtrutemp = mtru2temp;

		mtrl2temp = mtrdummy;
		mtru2temp = mtrdummy;

// Update current block row by the previous ones
/*
		for (j=_iabextt[iblk];j<_iabextt[iblk+1]-1;j++) {

			jblk = _jabextt[j];

			if (jblk >= iblkbeg) {

//				char strbufff[256];
//				sprintf (strbufff, "%s%i%s%i%s","UpdateLU_",iblk,"_",jblk,".dat");
//				ofstream ffout (strbufff);

//				ffout << " Iblk = " << iblk << endl;
//				ffout << " L bef fct = " << mtrltemp << " U bef fct = " << mtrutemp << endl;
//				ffout << " Jblk = " << jblk << endl;
//				ffout << " L upd = " << _mtrlschr.mtrarr[jblk] << " U upd = " << _mtruschr.mtrarr[jblk] << endl;

				_fct.Ilu2UpdateBlkRow2Index (iblk, _param,
														_mtral, mtrltemp, mtrutemp,
														_mtrlschr.mtrarr[jblk], _mtruschr.mtrarr[jblk],
														mtrl2temp, mtru2temp);

//				ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

				mtrltemp = mtrl2temp;
				mtrutemp = mtru2temp;

				mtrl2temp = mtrdummy;
				mtru2temp = mtrdummy;

			};

		};
*/
#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagupd_end, 0, NULL);
#endif

// Factorize current block row

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_beg, 0, NULL);
#endif

//		if (_i_lev_work > 0) {
//			char strbuff[256];
//			sprintf (strbuff, "%s%i%s","LU_",iblk,".dat");
//			ofstream fout (strbuff);
//			fout << " Iblk = " << iblk << endl;
//			fout << " L bef fct = " << mtrltemp << " U bef fct = " << mtrutemp << endl;
//		};

		_fct.Ilu2FctBlkRow2Index (iblk, _param, 
											_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
											_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
											_mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk]);

		_mtru.FreeAndCopySparsityPointers (iblk, _mtrl.mtrarr[iblk]);

		_fct.InitGlobalArraysByBlockData2Index (_mtrlschr.mtrarr[iblk]);

//		fout << " L aft fct = " << _mtrl.mtrarr[iblk] << " U aft fct = " << _mtru.mtrarr[iblk] << endl;

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_end, 0, NULL);
#endif

		mtrltemp = mtrdummy;
		mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

		nlstblk = _iabext[iblk+1]-_iabext[iblk];
		ibs = _iabext[iblk];

		_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
															_mtral, _fct);

// Free second order matrices if necessary

		if (icount%istepfree == 0 && iblk != iblkend) {

			for (jblk=iblkbeg;jblk<=iblk;jblk++) {

				_mtrlschr.FilterBlock2Index  (jblk,iblk);
				_mtruschr.FilterBlock2Index  (jblk,iblk);
				_fct.ReInitGlobalArraysByBlockData2Index (_mtrlschr.mtrarr[jblk]);

			};

		};

	};

// Scan second order matrices and filter unnecessary data

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		_mtrlschr.FilterBlock2Index  (iblk,iblkend);
		_mtruschr.FilterBlock2Index  (iblk,iblkend);
		_fct.ReInitGlobalArraysByBlockData2Index (_mtrlschr.mtrarr[iblk]);

	};

	delete [] listupd;

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_end, 0, NULL);
#endif

//	cout << " FctNode2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform computations with the current tree node including local Schur updates
//========================================================================================
void CGSMatrixR::PriorityUpdates2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
														int *_i_beg_work, int *_i_end_work, int _iblkbeg,
														int &_remained_updates, 
														const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
														int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
														CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
														CSlvParam &_param) {

	const char *funcname = "PriorityUpdates2Index";

// Create MPE descriptions

//	cout << " PriorityUpdates2Index: begin " << endl;

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctschrblk_beg;
	static int fctschrblk_end;
	static int fctschrblkini_beg;
	static int fctschrblkini_end;
	static int fctschrblkupd_beg;
	static int fctschrblkupd_end;
	static int fctschrblkfin_beg;
	static int fctschrblkfin_end;

	if (i_set == 0) {

		fctschrblk_beg = MPE_Log_get_event_number ();
		fctschrblk_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblk_beg, fctschrblk_end, "FctSchurBlk", "purple");
		fctschrblkini_beg = MPE_Log_get_event_number ();
		fctschrblkini_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkini_beg, fctschrblkini_end, "FctSchurBlkIni", "green");
		fctschrblkupd_beg = MPE_Log_get_event_number ();
		fctschrblkupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkupd_beg, fctschrblkupd_end, "FctSchurBlkUpd", "black");
		fctschrblkfin_beg = MPE_Log_get_event_number ();
		fctschrblkfin_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkfin_beg, fctschrblkfin_end, "FctSchurBlkFin", "white");

		i_set = 1;

	};

#endif

// Update Schur blocks according to the transposed block sparsity

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblk_beg, 0, NULL);
#endif

	int *listupd;

	listupd = new int [_mtrlschr.nblksr];
	if (!listupd) MemoryFail (funcname);

	int myid = _tree.myid;

	int iblkend = _i_end_work[_i_lev_work+1];

	int iblk = iblkend-_remained_updates+1;

// Load diagonal data

	int nlstblk = _iabext[iblk+1]-_iabext[iblk];
	int ibs = _iabext[iblk];
	int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
	int ibst = _iabextt[iblk];

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkini_beg, 0, NULL);
#endif

	_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
														_mtrlschr, _fct);

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkini_end, 0, NULL);
#endif

// Update

	int j, jblk;

	CSMatrixR mtrdummy, mtrl2temp, mtru2temp;

	int nlistupd = 0;

	int ibeg = _iptrextt[iblk];

	for (j=ibeg;j<_iabextt[iblk+1];j++) {

		jblk = _jabextt[j];

		if (jblk < _iblkbeg && _blk2cpu[jblk] == myid) {

			_iptrextt[iblk] = j+1;

			listupd[nlistupd] = jblk;
			nlistupd++;

		};

	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkupd_beg, 0, NULL);
#endif

	if (nlistupd > 0) {

		_fct.Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, _param,
																_mtrlschr, _mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk],
																nlistupd, listupd, _mtrlschr.mtrarr, _mtruschr.mtrarr,
																mtrl2temp, mtru2temp);

//			ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

		_mtrlschr.mtrarr[iblk] = mtrl2temp;
		_mtruschr.mtrarr[iblk] = mtru2temp;

		mtrl2temp = mtrdummy;
		mtru2temp = mtrdummy;

	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkupd_end, 0, NULL);
#endif
/*
	for (j=ibeg;j<_iabextt[iblk+1];j++) {

		jblk = _jabextt[j];

		if (jblk < _iblkbeg && _blk2cpu[jblk] == myid) {

			_iptrextt[iblk] = j+1;

//			char strbufff[256];
//			sprintf (strbufff, "%s%i%s%i%s%i%s","UpdateSchur_",iblk,"_",jblk,"_cpu_",_treeini.myid,".dat");
//			ofstream ffout (strbufff);
//			ffout << " Iblk = " << iblk << endl;
//			ffout << " L bef fct = " << _mtrlschr.mtrarr[iblk] << " U bef fct = " << _mtruschr.mtrarr[iblk] << endl;
//			ffout << " Jblk = " << jblk << endl;
//			ffout << " L upd = " << _mtrlschr.mtrarr[jblk] << " U upd = " <<_mtruschr.mtrarr[jblk] << endl;

			_fct.Ilu2UpdateBlkRow2Index (iblk, _param,
													_mtrlschr, _mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk],
													_mtrlschr.mtrarr[jblk], _mtruschr.mtrarr[jblk],
													mtrl2temp, mtru2temp);

//			ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

			_mtrlschr.mtrarr[iblk] = mtrl2temp;
			_mtruschr.mtrarr[iblk] = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

		};

	};
*/
// Unload diagonal data

	nlstblk = _iabext[iblk+1]-_iabext[iblk];
	ibs = _iabext[iblk];

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkfin_beg, 0, NULL);
#endif

	_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
														_mtrlschr, _fct);

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkfin_end, 0, NULL);
#endif

	delete [] listupd;

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblk_end, 0, NULL);
#endif

	_remained_updates--;

//	cout << " PriorityUpdates2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform computations with the current tree node including local Schur updates
//========================================================================================
void CGSMatrixR::NonPriorityUpdates2Index (int _i_lev_work, // Perform computations with the current tree node including local Schur updates
														int *_i_beg_work, int *_i_end_work,
														int &_ip_work, int _ip_work_end, int *_listschur, int *_updblkschur,
														const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
														int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
														CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
														CSlvParam &_param) {

//	const char *funcname = "NonPriorityUpdates2Index";

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctschrsblk_beg;
	static int fctschrsblk_end;

	if (i_set == 0) {

		fctschrsblk_beg = MPE_Log_get_event_number ();
		fctschrsblk_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrsblk_beg, fctschrsblk_end, "FctSchurSBlk", "black");

		i_set = 1;

	};

#endif

	if (_ip_work > _ip_work_end) return;

// Update Schur blocks according to the transposed block sparsity

//	cout << " PriorityUpdates2Index: begin " << endl;

	int myid = _tree.myid;

	int iblkbeg = _i_beg_work[_i_lev_work];

	int iblk = _listschur[_ip_work];
	int ilev = _updblkschur[_ip_work];

	ilev = _i_lev_work;
	_updblkschur[_ip_work] = ilev;

	_updblkschur[_ip_work]++;

	if (ilev > _i_lev_work) {
		_ip_work++;
		return;
	};

	int jblkbeg = _i_beg_work[ilev];

	if (jblkbeg > iblkbeg) jblkbeg = iblkbeg;

// Load diagonal data

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrsblk_beg, 0, NULL);
#endif

	int nlstblk = _iabext[iblk+1]-_iabext[iblk];
	int ibs = _iabext[iblk];
	int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
	int ibst = _iabextt[iblk];

	_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
														_mtrlschr, _fct);

// Update

	int j, jblk;

	CSMatrixR mtrdummy, mtrl2temp, mtru2temp;

	int ibeg = _iptrextt[iblk];

	for (j=ibeg;j<_iabextt[iblk+1];j++) {

		jblk = _jabextt[j];

		if (jblk < jblkbeg && _blk2cpu[jblk] == myid) {

			_iptrextt[iblk] = j+1;

//			char strbufff[256];
//			sprintf (strbufff, "%s%i%s%i%s%i%s","UpdateSchur_",iblk,"_",jblk,"_cpu_",_treeini.myid,".dat");
//			ofstream ffout (strbufff);
//			ffout << " Iblk = " << iblk << endl;
//			ffout << " L bef fct = " << _mtrlschr.mtrarr[iblk] << " U bef fct = " << _mtruschr.mtrarr[iblk] << endl;
//			ffout << " Jblk = " << jblk << endl;
//			ffout << " L upd = " << _mtrlschr.mtrarr[jblk] << " U upd = " <<_mtruschr.mtrarr[jblk] << endl;

			_fct.Ilu2UpdateBlkRow2Index (iblk, _param,
													_mtrlschr, _mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk],
													_mtrlschr.mtrarr[jblk], _mtruschr.mtrarr[jblk],
													mtrl2temp, mtru2temp);

//			ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

			_mtrlschr.mtrarr[iblk] = mtrl2temp;
			_mtruschr.mtrarr[iblk] = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

		};


	};

// Unload diagonal data

	nlstblk = _iabext[iblk+1]-_iabext[iblk];
	ibs = _iabext[iblk];

	_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
														_mtrlschr, _fct);

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrsblk_end, 0, NULL);
#endif

//	cout << " PriorityUpdates2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Perform computations with the current iblk
//========================================================================================
void CGSMatrixR::UpdateIblk2Index (int _iblk, int _iblkbeg, // Perform computations with the current iblk
												const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag, 
												int *_blk2cpu, int *_iabext, int *_jabext, int *_iabextt, int *_iptrextt, int *_jabextt,
												CGSMatrixR &_mtrlschr, CGSMatrixR &_mtruschr,
												CSlvParam &_param) {

	const char *funcname = "UpdateIblk2Index";

// Create MPE descriptions

//	cout << " PriorityUpdates2Index: begin " << endl;

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctschrblk_beg;
	static int fctschrblk_end;
	static int fctschrblkini_beg;
	static int fctschrblkini_end;
	static int fctschrblkupd_beg;
	static int fctschrblkupd_end;
	static int fctschrblkfin_beg;
	static int fctschrblkfin_end;

	if (i_set == 0) {

		fctschrblk_beg = MPE_Log_get_event_number ();
		fctschrblk_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblk_beg, fctschrblk_end, "FctSchurBlk", "purple");
		fctschrblkini_beg = MPE_Log_get_event_number ();
		fctschrblkini_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkini_beg, fctschrblkini_end, "FctSchurBlkIni", "green");
		fctschrblkupd_beg = MPE_Log_get_event_number ();
		fctschrblkupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkupd_beg, fctschrblkupd_end, "FctSchurBlkUpd", "black");
		fctschrblkfin_beg = MPE_Log_get_event_number ();
		fctschrblkfin_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctschrblkfin_beg, fctschrblkfin_end, "FctSchurBlkFin", "white");

		i_set = 1;

	};

#endif

// Update Schur blocks according to the transposed block sparsity

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblk_beg, 0, NULL);
#endif

	int *listupd;

	listupd = new int [_mtrlschr.nblksr];
	if (!listupd) MemoryFail (funcname);

	int myid = _tree.myid;

	int iblk = _iblk;

// Load diagonal data

	int nlstblk = _iabext[iblk+1]-_iabext[iblk];
	int ibs = _iabext[iblk];
	int nlstblkcol = _iabextt[iblk+1]-_iabextt[iblk];
	int ibst = _iabextt[iblk];

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkini_beg, 0, NULL);
#endif

	_fctdiag.Ilu2LoadDiagonalData2Index (iblk, nlstblk, _jabext+ibs, nlstblkcol, _jabextt+ibst,
														_mtrlschr, _fct);

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkini_end, 0, NULL);
#endif

// Update

	int j, jblk;

	CSMatrixR mtrdummy, mtrl2temp, mtru2temp;

	int nlistupd = 0;

	int ibeg = _iptrextt[iblk];

	for (j=ibeg;j<_iabextt[iblk+1];j++) {

		jblk = _jabextt[j];

		if (jblk < _iblkbeg && _blk2cpu[jblk] == myid) {

			_iptrextt[iblk] = j+1;

			listupd[nlistupd] = jblk;
			nlistupd++;

		};

	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkupd_beg, 0, NULL);
#endif

	if (nlistupd > 0) {

		_fct.Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, _param,
																_mtrlschr, _mtrlschr.mtrarr[iblk], _mtruschr.mtrarr[iblk],
																nlistupd, listupd, _mtrlschr.mtrarr, _mtruschr.mtrarr,
																mtrl2temp, mtru2temp);

//			ffout << " Res L = " << mtrl2temp << " Res U = " << mtru2temp << endl;

		_mtrlschr.mtrarr[iblk] = mtrl2temp;
		_mtruschr.mtrarr[iblk] = mtru2temp;

		mtrl2temp = mtrdummy;
		mtru2temp = mtrdummy;

	};

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkupd_end, 0, NULL);
#endif

// Unload diagonal data

	nlstblk = _iabext[iblk+1]-_iabext[iblk];
	ibs = _iabext[iblk];

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkfin_beg, 0, NULL);
#endif

	_fctdiag.Ilu2UnloadDiagonalData2Index (iblk, nlstblk, _jabext+ibs,
														_mtrlschr, _fct);

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblkfin_end, 0, NULL);
#endif

	delete [] listupd;

#ifdef __PMPITRACE__
	MPE_Log_event (fctschrblk_end, 0, NULL);
#endif

//	cout << " PriorityUpdates2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Send data to the other cpu's
//========================================================================================
void CGSMatrixR::AsyncSendBlocks2Index (int _i_lev_work, int *_lev_work, // Send data to the other cpu's
														int *_i_beg_work, int *_i_end_work, int _nlistschur, int *_listschur,
														const CTree &_tree, CFctR &_fct, CFctDiagR &_fctdiag,
														int _id_shift, int _nsendsmax, int &_nsends, 
														int *_imasksend, char **_psenddata, CMPIRequest *_sendreqvarr,
														int *_blk2cpu,
														int *_iabextflt, int *_jabextflt,
														int *_iabextt, int *_iptrextt, int *_jabextt,
														CGSMatrixR &_mtralsch, CGSMatrixR &_mtrausch,
														CSlvParam &_param) {

	const char *funcname = "AsyncSendBlocks2Index";

//	cout << " AsyncSendBlocks2Index: begin " << endl;

// Get communicator

	int i_lev_work = _i_lev_work+1;

	CMPIComm commloc;
	CMPIComm commglob = _tree.comm;

	int itree = _lev_work[i_lev_work*3];
	int inode = _lev_work[i_lev_work*3+1];

	if (itree < 0) {
		commloc = _tree.nodes[inode].commnode;
	} else {
		commloc = _tree.subtreearr[itree].nodes[inode].commnode;
	};

	int myid = _tree.GetMyid ();

// Init processors description

	int nprocloc;
	int *cpuidarrnodeloc;

	int i, j, iproc;
	int *iasend, *iptr;
	int *jasend;

// Prepare send lists

// Mark: log2 scanning the tree

	if (true) {

		itree = _lev_work[_i_lev_work*3];
		inode = _lev_work[_i_lev_work*3+1];

//		int iblkbeg = _i_beg_work[_i_lev_work];
		int iblkend = _i_end_work[_i_lev_work];

		int iprocsnd, iproc1, ifathernode;

		iprocsnd = -1;
		if (itree < 0) {
			iproc = _tree.nodes[inode].nodecpu;
			ifathernode = _tree.nodes[inode].fatherid;
//			cout << " myid = " << myid << " inode = " << inode << " iproc = " << iproc << " ifather = " << ifathernode << endl;
			if (iproc == myid && ifathernode != inode) {
				iproc1 = _tree.nodes[ifathernode].nodecpu;
				if (iproc1 != myid) iprocsnd = iproc1;
			};
		} else {
			iproc = _tree.subtreearr[itree].nodes[inode].nodecpu;
			ifathernode = _tree.subtreearr[itree].nodes[inode].fatherid;
			if (iproc == myid && ifathernode != inode) {
				iproc1 = _tree.subtreearr[itree].nodes[ifathernode].nodecpu;
				if (iproc1 != myid) iprocsnd = iproc1;
			};
		};

//		cout << " Send init: myid = " << myid << " iblkbeg = " << iblkbeg << " iblkend = " << iblkend << " iprocsnd = " << iprocsnd << endl;

		if (iprocsnd >= 0) {
			nprocloc = 1;
		} else {
			nprocloc = 0;
		};

		cpuidarrnodeloc = new int [nprocloc];
		if (!cpuidarrnodeloc) MemoryFail (funcname);

		if (iprocsnd >= 0) {
			cpuidarrnodeloc[0] = iprocsnd;
			iasend = new int [nprocloc+1];
			if (!iasend) MemoryFail (funcname);
			iptr = new int [nprocloc+1];
			if (!iptr) MemoryFail (funcname);
			jasend = new int [_nlistschur];
			if (!jasend) MemoryFail (funcname);
			int iblk;
			int nz = 0;
			for (i=0;i<_nlistschur;i++) {
				iblk = _listschur[i];
				if (iblk > iblkend) {
					jasend[nz] = iblk;
					nz++;
					UpdateIblk2Index (iblk, iblkend+1,
											_tree, _fct, _fctdiag, 
											_blk2cpu, _iabextflt, _jabextflt, _iabextt, _iptrextt, _jabextt,
											_mtralsch, _mtrausch,
											_param);
				};
			};
//			OutArr (cout," list = ",nz,jasend);
			iasend[0] = 0;
			iasend[1] = nz;
		} else {
			iasend = new int [nprocloc+1];
			if (!iasend) MemoryFail (funcname);
			iptr = new int [nprocloc+1];
			if (!iptr) MemoryFail (funcname);
			jasend = new int [nprocloc+1];
			if (!jasend) MemoryFail (funcname);
		};

	} else {

		int *pcpuidarrnode;

		if (itree < 0) {
			nprocloc = _tree.nodes[inode].nprocnode;
			pcpuidarrnode = _tree.nodes[inode].cpuidarrnode;
		} else {
			nprocloc = _tree.subtreearr[itree].nodes[inode].nprocnode;
			pcpuidarrnode = _tree.subtreearr[itree].nodes[inode].cpuidarrnode;
		};

		int myidloc = -1;

		for (j=0;j<nprocloc;j++) {
			if (pcpuidarrnode[j] == myid) myidloc = j;
		};

		if (myidloc < 0) throw " CGSMatrixR::AsyncSendBlocks2Index: myid not found ";

		cpuidarrnodeloc = new int [nprocloc];
		if (!cpuidarrnodeloc) MemoryFail (funcname);

		for (i=0;i<nprocloc;i++) cpuidarrnodeloc[i] = pcpuidarrnode[i];

		iasend = new int [nprocloc+1];
		if (!iasend) MemoryFail (funcname);
		iptr = new int [nprocloc];
		if (!iptr) MemoryFail (funcname);

		for (i=0;i<=nprocloc;i++) iasend[i] = 0;

		int iblkbeg = _i_beg_work[i_lev_work];
		int iblkend = _i_end_work[i_lev_work];

		int jblk, iproc;

		for (i=0;i<nprocloc;i++) {
			iproc = cpuidarrnodeloc[i];
			for (jblk=iblkbeg;jblk<=iblkend;jblk++) {
				if (_blk2cpu[jblk] == iproc && iproc != myid) iasend[i+1]++;
			};
		};

		for (i=0;i<nprocloc;i++) iasend[i+1] = iasend[i]+iasend[i+1];
		for (i=0;i<nprocloc;i++) iptr[i] = iasend[i];

		int nzjasend = iasend[nprocloc];

		jasend = new int [nzjasend];
		if (!jasend) MemoryFail (funcname);

		int k;

		for (i=0;i<nprocloc;i++) {
			iproc = cpuidarrnodeloc[i];
			for (jblk=iblkbeg;jblk<=iblkend;jblk++) {
				if (_blk2cpu[jblk] == iproc && iproc != myid) {
					k = iptr[i];
					jasend[k] = jblk;
					iptr[i]++;
				};
			};
		};

	};

// Create send structures, free Schur blocks and send

	bool last_work_with_Schur = false;
	if (itree >= 0) last_work_with_Schur = true;
	if (_tree.nsubtree == 0) last_work_with_Schur = true;

	CSMatrixR mtrdummy;

	CSMatrixR mtrlsend, mtrusend, mtrdsend;

	double thetaloc = _param.theta / 2.0e0;

	int lenobj, lenobj1, lenobjl, lenobju, lenobjd, jblk;
	char *ptr, *pobj, *pobjl, *pobju, *pobjd;
	int *piarr;

	for (i=0;i<nprocloc;i++) {

//		char strbufff11[256];
//		sprintf (strbufff11, "%s%i%s%i%s%i%s","BefFiltrSnd_",_ilev,"_ind_",i,"_cpu_",_treeini.myid,".dat");
//		ofstream fffout (strbufff11,ios::app);
		for (j=iasend[i];j<iasend[i+1];j++) {
			jblk = jasend[j];
//			fffout << " Jblk = " << jblk << " L = " << _mtralsch.mtrarr[jblk] << " U = " << _mtrausch.mtrarr[jblk] << endl;
		};

		_mtralsch.FilterListOfSchurBlocks2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																_iabextflt, _jabextflt, thetaloc,
																_fct);
		_mtrausch.FilterListOfSchurBlocks2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																_iabextflt, _jabextflt, thetaloc,
																_fct);

		for (j=iasend[i];j<iasend[i+1];j++) {
			jblk = jasend[j];
//			fffout << " Filtered Jblk = " << jblk << " L = " << _mtralsch.mtrarr[jblk] << " U = " << _mtrausch.mtrarr[jblk] << endl;
		};
		iproc = cpuidarrnodeloc[i];
		if (iproc != myid) {
			_mtralsch.CombineMatrices2Index (true, iasend[i+1]-iasend[i], jasend+iasend[i], 
														_mtralsch.mtrarr, mtrlsend);
//			fffout << " LSend = " << mtrlsend << endl;
			for (j=iasend[i];j<iasend[i+1];j++) {
				jblk = jasend[j];
				_mtralsch.mtrarr[jblk] = mtrdummy;
			};
			if (!last_work_with_Schur) {
				for (j=iasend[i];j<iasend[i+1];j++) {
					jblk = jasend[j];
					_mtralsch.mtrarr[jblk] = _fct.CreateBlock2Index (jblk, _mtralsch);
				};
			};
			_mtrausch.CombineMatrices2Index (true, iasend[i+1]-iasend[i], jasend+iasend[i], 
														_mtrausch.mtrarr, mtrusend);
//			fffout << " USend = " << mtrusend << endl;
			for (j=iasend[i];j<iasend[i+1];j++) {
				jblk = jasend[j];
				_mtrausch.mtrarr[jblk] = mtrdummy;
			};
			if (!last_work_with_Schur) {
				for (j=iasend[i];j<iasend[i+1];j++) {
					jblk = jasend[j];
					_mtrausch.mtrarr[jblk] = _fct.CreateBlock2Index (jblk, _mtralsch);
				};
			};
			CSMatrixR::SuperSparsifySymm2Index (true, mtrlsend, mtrusend);
//			fffout << " After sparsify LSend = " << mtrlsend << " USend = " << mtrusend << endl;

//			char strbufff[256];
//			sprintf (strbufff, "%s%i%s%i%s","UpdateSchurSndL_",i,"_cpu_",myid,".dat");
//			ofstream ffout (strbufff);
//			ffout << " L send = " << mtrlsend << endl;
//			char strbufff1[256];
//			sprintf (strbufff1, "%s%i%s%i%s","UpdateSchurSndU_",i,"_cpu_",myid,".dat");
//			ofstream ffout1 (strbufff1);
//			ffout1 << " U send = " << mtrusend << endl;

//			cout << " Unpack size send = " << mtrlsend.GetNzja () << endl;

			mtrlsend.PackMatrix (lenobjl,pobjl);
			mtrlsend = mtrdummy;
			mtrusend.SetNlist (0);
			mtrusend.SetNlist2 (0);
			mtrusend.SetNzja (0);
			mtrusend.SetNzja2 (0);
			mtrusend.SetNzja3 (0);
			mtrusend.PackMatrix (lenobju,pobju);
			mtrusend = mtrdummy;
			mtrdsend = _mtralsch.CreateDiagonalUpdate2Index (iasend[i+1]-iasend[i], jasend+iasend[i], 
																				_fctdiag);
//			char strbufff1[256];
//			sprintf (strbufff1, "%s%i%s%i%s%i%s","UpdateSchurSndDia_",i,"_lev_",_ilev,"_cpu_",myid,".dat");
//			ofstream ffout1 (strbufff1);
//			ffout1 << " D send = " << mtrdsend << endl;
			mtrdsend.SetNzja (0);
			mtrdsend.SetNzja2 (0);
			mtrdsend.SetNzja3 (0);
			mtrdsend.PackMatrix (lenobjd,pobjd);

			lenobj = 3*sizeof(int) + lenobjl + lenobju + lenobjd;

			pobj = new char [lenobj];
			if(!pobj) MemoryFail(funcname);

			ptr = pobj;

			piarr = (int *) pobj;

			piarr[0] = lenobjl;
			piarr[1] = lenobju;
			piarr[2] = lenobjd;

			pobj += 3*sizeof(int);

			memcpy (pobj,pobjl,lenobjl);
			pobj += lenobjl;
			memcpy (pobj,pobju,lenobju);
			pobj += lenobju;
			memcpy (pobj,pobjd,lenobjd);
			pobj += lenobjd;

			delete [] pobjl;
			delete [] pobju;
			delete [] pobjd;

// Two sequential sends

			if (_nsends+1 >= _nsendsmax) throw " CGSMatrixR::AsyncSendBlocks2Index: to small number of sends reserved ";

			lenobj1 = sizeof(int);

			pobj = new char [lenobj1];
			if(!pobj) MemoryFail(funcname);

			piarr = (int *) pobj;
			piarr[0] = lenobj;

			_imasksend[_nsends] = 1;
			_psenddata[_nsends] = pobj;

//			cout << " Send from myid = " << myid << " to iproc = " << iproc << endl;
			CMPIExchange::ISend (commglob, iproc, i_lev_work,
										lenobj1, pobj, _sendreqvarr[_nsends]);

			_nsends++;

			_imasksend[_nsends] = 1;
			_psenddata[_nsends] = ptr;

//			cout << " Send: i_lev_work = " << i_lev_work << " iproc = " << i << " len = " << lenobj << endl;

			CMPIExchange::ISend (commglob, iproc, _id_shift+i_lev_work,
										lenobj, ptr, _sendreqvarr[_nsends]);

			_nsends++;

		};

	};

// Free work arrays

	delete [] cpuidarrnodeloc;
	delete [] iasend;
	delete [] iptr;
	delete [] jasend;

//	cout << " AsyncSendBlocks2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Init to receive some data
//========================================================================================
void CGSMatrixR::InitReceiveBlocks2Index (int _i_lev_work, int *_lev_work, const CTree &_tree, // Init to receive some data
														int &_nrecvs, int *_imaskrecv, char **_precvdata, CMPIRequest *_recvreqvarr) {

	const char *funcname = "InitReceiveBlocks2Index";

	if (_i_lev_work == 0) return;

//	cout << " ReceiveBlocks2Index: begin " << endl;

// Get communicator

	int myid = _tree.myid;

	CMPIComm commloc;
	CMPIComm commglob = _tree.comm;

	int itree = _lev_work[_i_lev_work*3];
	int inode = _lev_work[_i_lev_work*3+1];

	if (itree < 0) {
		commloc = _tree.nodes[inode].commnode;
	} else {
		commloc = _tree.subtreearr[itree].nodes[inode].commnode;
	};

	int nlistloc;
	int *plist;

	if (itree < 0) {
		nlistloc = _tree.nodes[inode].nprocnode;
		plist = _tree.nodes[inode].cpuidarrnode;
	} else {
		nlistloc = _tree.subtreearr[itree].nodes[inode].nprocnode;
		plist = _tree.subtreearr[itree].nodes[inode].cpuidarrnode;
	};

	int i, iproc;

// Mark: log2 scanning the tree

	if (true) {

// Receive from the child only

		int iprocrcv, nchildsloc, iproc1, iproc2, ichild;

		iprocrcv = -1;
		if (itree < 0) {
			iproc = _tree.nodes[inode].nodecpu;
			nchildsloc = _tree.nodes[inode].nchilds;
			if (iproc == myid && nchildsloc > 1) {
				ichild = _tree.nodes[inode].childs[0];
				iproc1 = _tree.nodes[ichild].nodecpu;
				ichild = _tree.nodes[inode].childs[1];
				iproc2 = _tree.nodes[ichild].nodecpu;
				if (iproc1 != myid) iprocrcv = iproc1;
				if (iproc2 != myid) iprocrcv = iproc2;
			};
		} else {
			iproc = _tree.subtreearr[itree].nodes[inode].nodecpu;
			nchildsloc = _tree.subtreearr[itree].nodes[inode].nchilds;
			if (iproc == myid && nchildsloc > 1) {
				ichild = _tree.subtreearr[itree].nodes[inode].childs[0];
				iproc1 = _tree.subtreearr[itree].nodes[ichild].nodecpu;
				ichild = _tree.subtreearr[itree].nodes[inode].childs[1];
				iproc2 = _tree.subtreearr[itree].nodes[ichild].nodecpu;
				if (iproc1 != myid) iprocrcv = iproc1;
				if (iproc2 != myid) iprocrcv = iproc2;
			};
		};

		if (iprocrcv >= 0) {

			int isize = sizeof(int);

			char *buffer;

			buffer = new char [isize];
			if (!buffer) MemoryFail (funcname);

//			cout << " Receive on myid = " << myid << " from iproc = " << iprocrcv << endl;
			CMPIExchange::IRecv (commglob, iprocrcv, _i_lev_work,
										isize, buffer, _recvreqvarr[_nrecvs]);

			_imaskrecv[_nrecvs] = 1;
			_precvdata[_nrecvs] = buffer;

			_nrecvs++;

		};

	} else {

// Receive from all processors

		for (i=0;i<nlistloc;i++) {

			iproc = plist[i];

			if (iproc != myid) {

// Allocate and receive the message

				int isize = sizeof(int);

				char *buffer;

				buffer = new char [isize];
				if (!buffer) MemoryFail (funcname);

				CMPIExchange::IRecv (commglob, iproc, _i_lev_work,
											isize, buffer, _recvreqvarr[_nrecvs]);

				_imaskrecv[_nrecvs] = 1;
				_precvdata[_nrecvs] = buffer;

				_nrecvs++;

			};

		};

	};

//	cout << " ReceiveBlocks2Index: end " << endl;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Test/wait sends and receives for completion and update if necessary
//========================================================================================
void CGSMatrixR::TestSendsRecvsAndUpdate2Index (const CTree &_tree, // Test/wait sends and receives for completion and update if necessary
																bool _remained_work, 
																int _id_shift, int _nsends, int &_isend, 
																int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr,
																int &_nrecvs, int &_irecv,
																int *_imaskrecv, char** _precvdata, CMPIRequest *_recvreqvarr,
																CFctR &_fct, CFctDiagR &_fctdiag,
																CGSMatrixR &_mtrlsch, CGSMatrixR &_mtrusch, CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd, 
																CSMatrixR *_mtrarrwork, CSMatrixR *_mtrarr2work) {

	const char *funcname = "TestSendsRecvsAndUpdate2Index";

// Allocate the memory

	int nloc = (_nsends+_nrecvs) - (_isend+_irecv);

	int *itypeloc, *indexloc, *indarrloc;
	CMPIRequest *reqvarrloc;
	CMPIStatus *statarrloc;

	itypeloc = new int [nloc];
	if (!itypeloc) MemoryFail (funcname);
	indexloc = new int [nloc];
	if (!indexloc) MemoryFail (funcname);
	indarrloc = new int [nloc];
	if (!indarrloc) MemoryFail (funcname);
	reqvarrloc = new CMPIRequest [nloc];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new CMPIStatus [nloc];
	if (!statarrloc) MemoryFail (funcname);

// Prepare data for checks

NextTest:;

	nloc = 0;

	int i;

	for (i=_isend;i<_nsends;i++) {
		if (_imasksend[i] != 0) {
			itypeloc[nloc] = 1;
			indexloc[nloc] = i;
			reqvarrloc[nloc] = _sendreqvarr[i];
			nloc++;
		};
	};
	for (i=_irecv;i<_nrecvs;i++) {
		if (_imaskrecv[i] != 0) {
			itypeloc[nloc] = 2;
			indexloc[nloc] = i;
			reqvarrloc[nloc] = _recvreqvarr[i];
			nloc++;
		};
	};

	CMPIComm commglob = _tree.comm;

	if (nloc == 0) goto EndTest;

// Check for completion

	int ncompl;
	int index;
	bool flag;

	if (!_remained_work) {
		CMPIExchange::WaitSome (nloc, reqvarrloc, ncompl, indarrloc, statarrloc);
	} else {
		CMPIExchange::TestSome (nloc, reqvarrloc, ncompl, indarrloc, statarrloc);
	};

	flag = false;

	if (ncompl > 0) flag = true;

// Update data if completed

	int itype, ind, iproc, isize, msgtag;
	int *pint;
	char *buffer;

	if (flag) {
		for (i=0;i<ncompl;i++) {
			index = indarrloc[i];
			itype = itypeloc[index];
			ind = indexloc[index];
			if (itype == 1) {
				delete [] _psenddata[ind];
				_imasksend[ind] = 0;
			} else {
				CMPIExchange::GetTag (statarrloc[i],msgtag);
				CMPIExchange::GetSource (statarrloc[i],iproc);
				if (msgtag < _id_shift) {
					pint = (int *)_precvdata[ind];
					isize = pint[0];

					buffer = new char [isize];
					if (!buffer) MemoryFail (funcname);

//					cout << " Receive 2 on myid = " << _tree.myid << " from iproc = " << iproc << endl;

					CMPIExchange::IRecv (commglob, iproc, _id_shift+msgtag,
												isize, buffer, _recvreqvarr[_nrecvs]);

					_imaskrecv[_nrecvs] = 1;
					_precvdata[_nrecvs] = buffer;

					_nrecvs++;

				} else {
//					cout << " Receive from iproc = " << iproc << endl;
					UseReceivedData2Index (_tree, _precvdata[ind], _fct, _fctdiag, 
													_mtrlsch, _mtrusch, _mtrladd, _mtruadd, _mtrarrwork, _mtrarr2work);
				};
				delete [] _precvdata[ind];
				_imaskrecv[ind] = 0;
			};
		};

		goto NextTest;

	};

// Free the memory and update isend and irecv

EndTest:;

	delete [] itypeloc;
	delete [] indexloc;
	delete [] indarrloc;
	delete [] reqvarrloc;
	delete [] statarrloc;

	int ibeg = _isend;

	for (i=ibeg;i<_nsends;i++) {
		if (_imasksend[i] == 0) {
			_isend++;
		} else {
			break;
		};
	};

	ibeg = _irecv;

	for (i=ibeg;i<_nrecvs;i++) {
		if (_imaskrecv[i] == 0) {
			_irecv++;
		} else {
			break;
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Use receive data for updates
//========================================================================================
void CGSMatrixR::UseReceivedData2Index (const CTree &_tree, char *_buffer, CFctR &_fct, CFctDiagR &_fctdiag, // Use receive data for updates
														CGSMatrixR &_mtrlsch, CGSMatrixR &_mtrusch, 
														CGSMatrixR &_mtrladd, CGSMatrixR &_mtruadd, 
														CSMatrixR *_mtrarrwork, CSMatrixR *_mtrarr2work) {

	const char *funcname = "UseReceivedData2Index";

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctupdrecv_beg;
	static int fctupdrecv_end;
	static int fctupdrecvunpk_beg;
	static int fctupdrecvunpk_end;
	static int fctupdrecvspl_beg;
	static int fctupdrecvspl_end;
	static int fctupdrecvadd_beg;
	static int fctupdrecvadd_end;

	if (i_set == 0) {

		fctupdrecv_beg = MPE_Log_get_event_number ();
		fctupdrecv_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctupdrecv_beg, fctupdrecv_end, "FctUpdRecv", "blue");
		fctupdrecvunpk_beg = MPE_Log_get_event_number ();
		fctupdrecvunpk_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctupdrecvunpk_beg, fctupdrecvunpk_end, "FctUpdRecvUnpk", "grey");
		fctupdrecvspl_beg = MPE_Log_get_event_number ();
		fctupdrecvspl_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctupdrecvspl_beg, fctupdrecvspl_end, "FctUpdRecvSplit", "red");
		fctupdrecvadd_beg = MPE_Log_get_event_number ();
		fctupdrecvadd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctupdrecvadd_beg, fctupdrecvadd_end, "FctUpdRecvAdd", "white");

		i_set = 1;

	};

#endif

// Unpack the message

#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecv_beg, 0, NULL);
#endif

	CSMatrixR mtrl, mtru, mtrd;

	int *piarr = (int*) _buffer;

	int isizel, isizeu, isized;

	isizel = piarr[0];
	isizeu = piarr[1];
	isized = piarr[2];

	char *pbuff = _buffer + 3*sizeof(int);

#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvunpk_beg, 0, NULL);
#endif

	mtrl.UnPackMatrix(isizel,pbuff);

	pbuff += isizel;

	mtru.UnPackMatrix(isizeu,pbuff);

	pbuff += isizeu;

	mtrd.UnPackMatrix(isized,pbuff);

//	char strbufff[256];
//	sprintf (strbufff, "%s%i%s","UpdateSchurRcvL_",mtrl.GetNzja(),".dat");
//	ofstream ffout (strbufff);
//	ffout << " L recv = " << mtrl << endl;
//	char strbufff1[256];
//	sprintf (strbufff1, "%s%i%s","UpdateSchurRcvU_",mtrl.GetNzja(),".dat");
//	ofstream ffout1 (strbufff1);
//	ffout1 << " U recv = " << mtru << endl;

#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvunpk_end, 0, NULL);
#endif

// Copy the L sparsity into U

	((CSMatrix &)mtru) = ((CSMatrix &)mtrl);

// Update all necessary data

	int nlistloc;
	int *listloc;
	CSMatrixR mtrdummy;

	listloc = new int [_mtrlsch.nblksr];
	if (!listloc) MemoryFail (funcname);

#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvspl_beg, 0, NULL);
#endif
	_mtrlsch.SplitMatrixList2Index (nlistloc, listloc,
												mtrl, 
												_mtrarrwork);
	_mtrusch.SplitMatrixList2Index (nlistloc, listloc,
												mtru, 
												_mtrarr2work);
#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvspl_end, 0, NULL);
#endif
	mtrl = mtrdummy;
	mtru = mtrdummy;
#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvadd_beg, 0, NULL);
#endif
//	_mtrlsch.AddReplaceSetOfBlocks2Index (_fct, nlistloc, listloc, 
//														_mtrlsch.mtrarr, _mtrarrwork, _mtrusch.mtrarr, _mtrarr2work);

	int j, iblk;

	for (j=0;j<nlistloc;j++) {
		iblk = listloc[j];
		if (_mtrladd.mtrarr[iblk].GetNzja () != 0) {
			CSMatrixR::DoubleAddMatrSS2Index (true,
															_mtrarrwork[iblk], _mtrladd.mtrarr[iblk], mtrl,
															_mtrarr2work[iblk], _mtruadd.mtrarr[iblk], mtru);
			_mtrladd.mtrarr[iblk] = mtrl;
			_mtruadd.mtrarr[iblk] = mtru;
			mtrl = mtrdummy;
			mtru = mtrdummy;
		} else {
			_mtrladd.mtrarr[iblk] = _mtrarrwork[iblk];
			_mtruadd.mtrarr[iblk] = _mtrarr2work[iblk];
		};
//		if (_mtrladd.mtrarr[iblk].GetNzja () > _mtrlsch.mtrarr[iblk].GetNzja ()) {
			if (_mtrlsch.mtrarr[iblk].GetNzja () != 0) {
				CSMatrixR::DoubleAddMatrSS2Index (true,
																_mtrladd.mtrarr[iblk], _mtrlsch.mtrarr[iblk], mtrl,
																_mtruadd.mtrarr[iblk], _mtrusch.mtrarr[iblk], mtru);
				_mtrlsch.mtrarr[iblk] = mtrl;
				_mtrusch.mtrarr[iblk] = mtru;
				mtrl = mtrdummy;
				mtru = mtrdummy;
			} else {
				_mtrlsch.mtrarr[iblk] = _mtrladd.mtrarr[iblk];
				_mtrusch.mtrarr[iblk] = _mtruadd.mtrarr[iblk];
			};
			_mtrladd.mtrarr[iblk] = mtrdummy;
			_mtruadd.mtrarr[iblk] = mtrdummy;
//		};
	};
#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecvadd_end, 0, NULL);
#endif

	for (j=0;j<nlistloc;j++) {
		iblk = listloc[j];
		_mtrarrwork[iblk] = mtrdummy;
		_mtrarr2work[iblk] = mtrdummy;
	};

	_mtrlsch.UpdateDiagonal2Index (mtrd, _fctdiag);
	mtrd = mtrdummy;

	delete [] listloc;

#ifdef __PMPITRACE__
	MPE_Log_event (fctupdrecv_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixR: Check sends for completion
//========================================================================================
void CGSMatrixR::CheckSends2Index (int _nsends, int &_isend, // Check sends for completion
												int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr) {

//	const char *funcname = "CheckSends2Index";

	int ibeg = _isend;

	int i;
	bool flag;
	CMPIStatus stat;

	for (i=ibeg;i<_nsends;i++) {
		if (_imasksend[i] == 1) {
			CMPIExchange::Test (_sendreqvarr[i], flag, stat);
			if (flag) {
				delete [] _psenddata[i];
				_imasksend[i] = 0;
				if (_isend == i) _isend++;
			};
		} else if (_isend == i) {
			_isend++;
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Wait sends for completion
//========================================================================================
void CGSMatrixR::WaitSends2Index (int _nsends, int _isend, // Check sends for completion
												int *_imasksend, char** _psenddata, CMPIRequest *_sendreqvarr) {

	const char *funcname = "WaitSends2Index";

	int nloc = _nsends-_isend;

	CMPIRequest *reqvarrloc;
	CMPIStatus *statarrloc;

	reqvarrloc = new CMPIRequest [nloc];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new CMPIStatus [nloc];
	if (!statarrloc) MemoryFail (funcname);

	nloc = 0;

	int i;

	for (i=_isend;i<_nsends;i++) {
		if (_imasksend[i] == 1) {
			reqvarrloc[nloc] = _sendreqvarr[i];
			nloc++;
		};
	};

	CMPIExchange::WaitAll (nloc, reqvarrloc, statarrloc);

	for (i=_isend;i<_nsends;i++) {
		if (_imasksend[i] == 1) {
			delete [] _psenddata[i];
		};
	};

	delete [] reqvarrloc;
	delete [] statarrloc;

};

// Author: Kharchenko S.A.
// CFctSchurR: Prepare fct working data
//========================================================================================
void CFctSchurR::PrepareWorkData () { // Prepare fct working data

	const char *funcname = "PrepareWorkData";

// Init all work structures

	delete [] imaskflt;

	icycleflt = -1;

	imaskflt = new int [nblks];
	if (!imaskflt) MemoryFail (funcname);

	int i;

	for (i=0;i<nblks;i++) imaskflt[i] = icycleflt;

// Init matrix structures

	*((CGSMatrixR *)pgmtrl) = *((CGSMatrixR *)pgmtral);
	*((CGSMatrixR *)pgmtru) = *((CGSMatrixR *)pgmtral);
	*((CGSMatrixR *)pgmtrlschr) = *((CGSMatrixR *)pgmtral);
	*((CGSMatrixR *)pgmtruschr) = *((CGSMatrixR *)pgmtral);
	*((CGSMatrixR *)pgmtrladd) = *((CGSMatrixR *)pgmtral);
	*((CGSMatrixR *)pgmtruadd) = *((CGSMatrixR *)pgmtral);

// Allocate send/receive structures

	int nblkmax = 0;

	int myidloc = ptree->GetMyid ();
	int nprocloc = ptree->GetNproc ();
	int nnodesloc = ptree->GetNnodes ();
	CNode *pnodes = ptree->GetNode ();

	int iblkbeg, iblkend, ni;

	for (i=0;i<nnodesloc;i++) {
		iblkbeg = pnodes[i].GetIndbeg ();
		iblkend = pnodes[i].GetIndend ();
		ni = iblkend-iblkbeg+1;
		if (ni>nblkmax) nblkmax = ni;
	};

	nsendsmax = nprocloc*nblkmax;

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	for (i=0;i<nsendsmax;i++) imasksend[i] = -1;

	nsends = 0;

// Compute the list of Schur blocks used in Fct

	int *piablksp = pablkstr->GetIa ();
	int *pjablksp = pablkstr->GetJa ();

	int j, iblk, jjblk, k, kkblk;

	icycleflt++;

	int nlistschur = 0;

	int *listschur;

	listschur = new int [nblks];
	if (!listschur) MemoryFail (funcname);

	for (iblk=0;iblk<nblks;iblk++) {
		if (pblk2cpu[iblk] == myidloc) {
			if (imaskflt[iblk] != icycleflt) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskflt[iblk] = icycleflt;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskflt[jjblk] != icycleflt) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskflt[jjblk] = icycleflt;
					};
				} else if (jjblk < iblk) {
					for (k=piablksp[jjblk];k<piablksp[jjblk+1];k++) {
						kkblk = pjablksp[k];
						if (kkblk >= iblk) {
							if (imaskflt[kkblk] != icycleflt) {
								listschur[nlistschur] = kkblk;
								nlistschur++;
								imaskflt[kkblk] = icycleflt;
							};
						};
					};
				};
			};
		};
	};

// Add tree data

	int *pcpuidend = ptree->GetCpuidend ();

	int inode = pcpuidend[myidloc];

	int indbeg, indend, ifather;

	while (true) {
		indbeg = pnodes[inode].GetIndbeg ();
		indend = pnodes[inode].GetIndend ();
		for (i=indbeg;i<=indend;i++) {
			if (imaskflt[i] != icycleflt) {
				listschur[nlistschur] = i;
				nlistschur++;
				imaskflt[i] = icycleflt;
			};
		};
		ifather = pnodes[inode].FatherId ();
		if (ifather == inode) break;
		inode = ifather;
	};

// Sort the list

	qsort (listschur, nlistschur, sizeof(int), compint);

// Create fct structures

	pfct = new CFctR (nblks, 0, 1);
	if (!pfct) MemoryFail (funcname);

	pfctdiag = new CFctDiagR (nblks);
	if (!pfctdiag) MemoryFail (funcname);

	((CFctR *)pfct)->AllocateGlobalArrays2Index (nblks, pblks, nlistschur, listschur);

// Free work arrays

	delete [] listschur;

};

// Author: Kharchenko S.A.
// CFctSchurR: Prepare the scaling
//========================================================================================
void CFctSchurR::PrepareScaling (int _nlistschur, int *_listschur) { // Prepare the scaling

//	const char *funcname = "PrepareScaling";

	int myidloc = ptree->GetMyid ();
	int iblk;

	for (iblk=0;iblk<nblks;iblk++) {

		if (pblk2cpu[iblk] == myidloc) {

			((CFctDiagR *)pfctdiag)->Ilu2ScaleBlkRow2Index (iblk, pparam->sclpttype, pparam->sclmin,
										*((CGSMatrixR *)pgmtral), ((CGSMatrixR *)pgmtrau)->GetMtrarr()[iblk], *((CFctR *)pfct));

		};

	};

// Update and output hystogram

	((CFctDiagR *)pfctdiag)->GetHyst()->MakeGlobalHystogram (ptree->GetComm ());
	((CFctDiagR *)pfctdiag)->GetHyst2()->MakeGlobalHystogram (ptree->GetComm ());

	if (myidloc == ptree->GetRootcpu ()) {
		if (pparam->msglev > 1) {
			cout  << " PointScaling" << *(((CFctDiagR *)pfctdiag)->GetHyst2());
			cout  << " Scaling" << *(((CFctDiagR *)pfctdiag)->GetHyst());
		};
		if (pparam->msglev > 0) {
			*pfout << " PointScaling" << *(((CFctDiagR *)pfctdiag)->GetHyst2());
			*pfout << " Scaling" << *(((CFctDiagR *)pfctdiag)->GetHyst());
		};
	};

// Exchange the scaling data from the bordering

	int *piablksp = pablkstr->GetIa ();
	int *pjablksp = pablkstr->GetJa ();

	((CFctDiagR *)pfctdiag)->ExchangeScalingGeneral (true, *ptree, pblk2cpu, _nlistschur, _listschur,
													piablksp, pjablksp, *((CGSMatrixR *)pgmtral));

};

// Author: Kharchenko S.A.
// CFctSchurR: Factorize the set of sequential blocks
//========================================================================================
void CFctSchurR::FctAndFilterBlocks (int _nlistfct, int *_listfct) { // Factorize the set of sequential blocks

	const char *funcname = "FctAndFilterBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

//	cout << " FctBlocks: begin " << endl;

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctdiag_beg;
	static int fctdiag_end;
	static int fctdiagupd_beg;
	static int fctdiagupd_end;
	static int fctdiagfct_beg;
	static int fctdiagfct_end;

	if (i_set == 0) {

		fctdiag_beg = MPE_Log_get_event_number ();
		fctdiag_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiag_beg, fctdiag_end, "FctDiag", "blue");

		fctdiagupd_beg = MPE_Log_get_event_number ();
		fctdiagupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagupd_beg, fctdiagupd_end, "FctDiagUpd", "gold");

		fctdiagfct_beg = MPE_Log_get_event_number ();
		fctdiagfct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagfct_beg, fctdiagfct_end, "FctDiagFct", "cyan");

		i_set = 1;

	};

#endif

// Perform computations with blocks of the current node taking into account Schur data if necessary

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_beg, 0, NULL);
#endif

	int iblkbeg = _listfct[0];
	int iblkend = _listfct[_nlistfct-1];

	int *piablkstr = pablkstr->GetIa ();
	int *pjablkstr = pablkstr->GetJa ();

	int iblk, j, jblk;

	CSMatrixR mtrltemp, mtrl2temp;
	CSMatrixR mtrutemp, mtru2temp;
	CSMatrixR mtrdummy;

	int *listupd;
	int *listrow;
	int *listcol;

	listupd = new int [nblks];
	if (!listupd) MemoryFail (funcname);
	listrow = new int [nblks];
	if (!listrow) MemoryFail (funcname);
	listcol = new int [nblks];
	if (!listcol) MemoryFail (funcname);

	int ilist, nlistrow, nlistcol, nlistupd, jjblk, k, kkblk;

	for (ilist=0;ilist<_nlistfct;ilist++) {

		iblk = _listfct[ilist];

		if (pparam->msglev > 1) {
			cout << " Myid = " << ptree->GetMyid () << " Iblk = " << iblk << endl;
		};

// Load diagonal data into fct structure

		icycleflt++;

		nlistrow = 0;

		if (imaskflt[iblk] != icycleflt) {
			listrow[nlistrow] = iblk;
			nlistrow++;
			imaskflt[iblk] = icycleflt;
		};

		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				if (imaskflt[jjblk] != icycleflt) {
					listrow[nlistrow] = jjblk;
					nlistrow++;
					imaskflt[jjblk] = icycleflt;
				};
			} else if (jjblk < iblk) {
				for (k=piablkstr[jjblk];k<piablkstr[jjblk+1];k++) {
					kkblk = pjablkstr[k];
					if (kkblk >= iblk) {
						if (imaskflt[kkblk] != icycleflt) {
							listrow[nlistrow] = kkblk;
							nlistrow++;
							imaskflt[kkblk] = icycleflt;
						};
					};
				};
			};
		};

		nlistcol = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk <= iblk) {
				listcol[nlistcol] = jjblk;
				nlistcol++;
			};
		};

		((CFctDiagR *)pfctdiag)->Ilu2LoadDiagonalData2Index (iblk, nlistrow, listrow, nlistcol, listcol, 
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

// Init current block row of the matrix

		((CFctR *)pfct)->Ilu2InitBlkRow2Index (iblk, pparam->dshift,
											((CGSMatrixR *)pgmtral)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtrau)->GetMtrarr()[iblk], 
											mtrltemp, mtrutemp);

// Update current block row by the Schur data

#ifdef __PMPITRACE__
		MPE_Log_event (fctdiagupd_beg, 0, NULL);
#endif

		if (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk].GetNzja () != 0) {

			mtrl2temp = ((CFctR *)pfct)->AddBlocks2Index (mtrltemp, ((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);
			mtru2temp = ((CFctR *)pfct)->AddBlocks2Index (mtrutemp, ((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk]);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

			((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
			((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrdummy;

		};

// Create the list and update

		nlistupd = 0;

		for (j=0;j<nlistcol-1;j++) {
			jblk = listcol[j];
			if (jblk >= iblkbeg) {
				listupd[nlistupd] = jblk;
				nlistupd++;
			};
		};

		if (nlistupd > 0) {

			((CFctR *)pfct)->Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, *pparam,
													*((CGSMatrixR *)pgmtral), 
													mtrltemp, mtrutemp,
													nlistupd, listupd, 
													((CGSMatrixR *)pgmtrlschr)->GetMtrarr(), ((CGSMatrixR *)pgmtruschr)->GetMtrarr(),
													mtrl2temp, mtru2temp);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

		};

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagupd_end, 0, NULL);
#endif

// Factorize current block row

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_beg, 0, NULL);
#endif

		nlistupd = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				listupd[nlistupd] = jjblk;
				nlistupd++;
			};
		};

		mtrltemp.FilterBlockRow2Index (true, nlistupd, listupd, icycleflt, imaskflt);
		mtrutemp.FilterBlockRow2Index (true, nlistupd, listupd, icycleflt, imaskflt);

		((CFctR *)pfct)->Ilu2FctBlkRow2Index (iblk, *pparam, 
											*((CGSMatrixR *)pgmtral), 
											mtrltemp, mtrutemp, ((CGSMatrixR *)pgmtral)->GetMtrarr()[iblk],
											((CGSMatrixR *)pgmtrl)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtru)->GetMtrarr()[iblk], 
											((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk]);

		((CGSMatrixR *)pgmtru)->FreeAndCopySparsityPointers (iblk, ((CGSMatrixR *)pgmtrl)->GetMtrarr()[iblk]);

		((CFctR *)pfct)->InitGlobalArraysByBlockData2Index (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_end, 0, NULL);
#endif

		mtrltemp = mtrdummy;
		mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

		((CFctDiagR *)pfctdiag)->Ilu2UnloadDiagonalData2Index (iblk, nlistrow, listrow,
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

	};

// Scan second order matrices and filter unnecessary data

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		((CGSMatrixR *)pgmtrlschr)->FilterBlock2Index  (iblk,iblkend);
		((CGSMatrixR *)pgmtruschr)->FilterBlock2Index  (iblk,iblkend);
		((CFctR *)pfct)->ReInitGlobalArraysByBlockData2Index (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

	};

	delete [] listupd;
	delete [] listrow;
	delete [] listcol;

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_end, 0, NULL);
#endif

//	cout << " FctBlocks: end " << endl;

};

// Author: Kharchenko S.A.
// CFctSymmSchurR: Factorize the set of sequential blocks
//========================================================================================
void CFctSymmSchurR::FctAndFilterBlocks (int _nlistfct, int *_listfct) { // Factorize the set of sequential blocks

	const char *funcname = "FctAndFilterBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

//	cout << " FctBlocks: begin " << endl;

// Create MPE descriptions

#ifdef __PMPITRACE__

	static int i_set = 0;
	static int fctdiag_beg;
	static int fctdiag_end;
	static int fctdiagupd_beg;
	static int fctdiagupd_end;
	static int fctdiagfct_beg;
	static int fctdiagfct_end;

	if (i_set == 0) {

		fctdiag_beg = MPE_Log_get_event_number ();
		fctdiag_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiag_beg, fctdiag_end, "FctDiag", "blue");

		fctdiagupd_beg = MPE_Log_get_event_number ();
		fctdiagupd_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagupd_beg, fctdiagupd_end, "FctDiagUpd", "gold");

		fctdiagfct_beg = MPE_Log_get_event_number ();
		fctdiagfct_end = MPE_Log_get_event_number ();
		MPE_Describe_state (fctdiagfct_beg, fctdiagfct_end, "FctDiagFct", "cyan");

		i_set = 1;

	};

#endif

// Perform computations with blocks of the current node taking into account Schur data if necessary

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_beg, 0, NULL);
#endif

	int iblkbeg = _listfct[0];
	int iblkend = _listfct[_nlistfct-1];

	int *piablkstr = pablkstr->GetIa ();
	int *pjablkstr = pablkstr->GetJa ();

	int iblk, j, jblk;

	CSMatrixR mtrltemp, mtrl2temp;
	CSMatrixR mtrdummy;

	int *listupd;
	int *listrow;
	int *listcol;

	listupd = new int [nblks];
	if (!listupd) MemoryFail (funcname);
	listrow = new int [nblks];
	if (!listrow) MemoryFail (funcname);
	listcol = new int [nblks];
	if (!listcol) MemoryFail (funcname);

	int ilist, nlistrow, nlistcol, nlistupd, jjblk, k, kkblk;

	for (ilist=0;ilist<_nlistfct;ilist++) {

		iblk = _listfct[ilist];

		if (pparam->msglev > 1) {
			cout << " Myid = " << ptree->GetMyid () << " Iblk = " << iblk << endl;
		};

// Load diagonal data into fct structure

		icycleflt++;

		nlistrow = 0;

		if (imaskflt[iblk] != icycleflt) {
			listrow[nlistrow] = iblk;
			nlistrow++;
			imaskflt[iblk] = icycleflt;
		};

		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				if (imaskflt[jjblk] != icycleflt) {
					listrow[nlistrow] = jjblk;
					nlistrow++;
					imaskflt[jjblk] = icycleflt;
				};
			} else if (jjblk < iblk) {
				for (k=piablkstr[jjblk];k<piablkstr[jjblk+1];k++) {
					kkblk = pjablkstr[k];
					if (kkblk >= iblk) {
						if (imaskflt[kkblk] != icycleflt) {
							listrow[nlistrow] = kkblk;
							nlistrow++;
							imaskflt[kkblk] = icycleflt;
						};
					};
				};
			};
		};

		nlistcol = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk <= iblk) {
				listcol[nlistcol] = jjblk;
				nlistcol++;
			};
		};

		((CFctDiagR *)pfctdiag)->Ilu2LoadDiagonalData2Index (iblk, nlistrow, listrow, nlistcol, listcol, 
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

// Init current block row of the matrix

		((CFctR *)pfct)->Ich2InitBlkRow2Index (iblk, pparam->dshift,
											((CGSMatrixR *)pgmtral)->GetMtrarr()[iblk], 
											mtrltemp);

// Update current block row by the Schur data

#ifdef __PMPITRACE__
		MPE_Log_event (fctdiagupd_beg, 0, NULL);
#endif

		if (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk].GetNzja () != 0) {

			mtrl2temp = ((CFctR *)pfct)->AddBlocks2Index (mtrltemp, ((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

			mtrltemp = mtrl2temp;

			mtrl2temp = mtrdummy;

			((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
			((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrdummy;

		};

// Create the list and update

		nlistupd = 0;

		for (j=0;j<nlistcol-1;j++) {
			jblk = listcol[j];
			if (jblk >= iblkbeg) {
				listupd[nlistupd] = jblk;
				nlistupd++;
			};
		};

		if (nlistupd > 0) {

			((CFctR *)pfct)->Ich2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, *pparam,
													*((CGSMatrixR *)pgmtral), 
													mtrltemp,
													nlistupd, listupd, 
													((CGSMatrixR *)pgmtrlschr)->GetMtrarr(),
													mtrl2temp);

			mtrltemp = mtrl2temp;

			mtrl2temp = mtrdummy;

		};

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagupd_end, 0, NULL);
#endif

// Factorize current block row

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_beg, 0, NULL);
#endif

		nlistupd = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				listupd[nlistupd] = jjblk;
				nlistupd++;
			};
		};

		mtrltemp.FilterBlockRow2Index (true, nlistupd, listupd, icycleflt, imaskflt);

		((CFctR *)pfct)->Ich2FctBlkRow2Index (iblk, *pparam, 
											*((CGSMatrixR *)pgmtral), 
											mtrltemp, ((CGSMatrixR *)pgmtral)->GetMtrarr()[iblk],
											((CGSMatrixR *)pgmtrl)->GetMtrarr()[iblk],
											((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

		((CFctR *)pfct)->InitGlobalArraysByBlockData2Index (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiagfct_end, 0, NULL);
#endif

		mtrltemp = mtrdummy;

// Unload diagonal data from fct structure

		((CFctDiagR *)pfctdiag)->Ilu2UnloadDiagonalData2Index (iblk, nlistrow, listrow,
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

	};

// Scan second order matrices and filter unnecessary data

	for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

		((CGSMatrixR *)pgmtrlschr)->FilterBlock2Index  (iblk,iblkend);
		((CFctR *)pfct)->ReInitGlobalArraysByBlockData2Index (((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);

	};

	delete [] listupd;
	delete [] listrow;
	delete [] listcol;

#ifdef __PMPITRACE__
	MPE_Log_event (fctdiag_end, 0, NULL);
#endif

//	cout << " FctBlocks: end " << endl;

};

// Author: Kharchenko S.A.
// CFctSchurR: Free second order data for the set of blocks
//========================================================================================
void CFctSchurR::FreeSecondOrderBlocks (int _nlistfct, int *_listfct) { // Free second order data for the set of blocks

//	const char *funcname = "FreeSecondOrderBlocks";

	CSMatrixR mtrdummy;

	int i, iblk;

	for (i=0;i<_nlistfct;i++) {
		iblk = _listfct[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrdummy;
	};

};

// Author: Kharchenko S.A.
// CFctSchurR: Init the set of Schur blocks
//========================================================================================
void CFctSchurR::InitSchurBlocks (int _nlistschur, int *_listschur) { // Init the set of Schur blocks

//	const char *funcname = "InitSchurBlocks";

	int i, iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = ((CFctR *)pfct)->CreateBlock2Index (iblk, *((CGSMatrixR *)pgmtral));
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = ((CFctR *)pfct)->CreateBlock2Index (iblk, *((CGSMatrixR *)pgmtral));
	};

};

// Author: Kharchenko S.A.
// CFctSymmSchurR: Init the set of Schur blocks
//========================================================================================
void CFctSymmSchurR::InitSchurBlocks (int _nlistschur, int *_listschur) { // Init the set of Schur blocks

//	const char *funcname = "InitSchurBlocks";

	int i, iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = ((CFctR *)pfct)->CreateBlock2Index (iblk, *((CGSMatrixR *)pgmtral));
	};

};

// Author: Kharchenko S.A.
// CFctSchurR: Update Schur data for the set of blocks
//========================================================================================
void CFctSchurR::UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur) { // Update Schur data for the set of blocks

	const char *funcname = "UpdateSchurBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

// Allocate work array

	int *listupd;
	int *listrow;
	int *listcol;

	listupd = new int [nblks];
	if (!listupd) MemoryFail (funcname);
	listrow = new int [nblks];
	if (!listrow) MemoryFail (funcname);
	listcol = new int [nblks];
	if (!listcol) MemoryFail (funcname);

// Main updates cycle

	int *piablkstr = pablkstr->GetIa ();
	int *pjablkstr = pablkstr->GetJa ();

	int i, iblk;
	int j, jjblk, nlistupd, nlistrow, nlistcol, k, kkblk;

	CSMatrixR mtrltemp, mtrutemp;
	CSMatrixR mtrdummy;

	for (i=0;i<_nlistschur;i++) {

		iblk = _listschur[i];

// Register blocks from the list

		icycleflt++;
		for (j=0;j<_nlistfct;j++) {
			jjblk = _listfct[j];
			imaskflt[jjblk] = icycleflt;
		};

// Create the updates list

		nlistupd = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk < iblk) {
				if (imaskflt[jjblk] == icycleflt) {
					listupd[nlistupd] = jjblk;
					nlistupd++;
				};
			};
		};

// Load diagonal data into fct structure

		icycleflt++;

		nlistrow = 0;

		if (imaskflt[iblk] != icycleflt) {
			listrow[nlistrow] = iblk;
			nlistrow++;
			imaskflt[iblk] = icycleflt;
		};

		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				if (imaskflt[jjblk] != icycleflt) {
					listrow[nlistrow] = jjblk;
					nlistrow++;
					imaskflt[jjblk] = icycleflt;
				};
			} else if (jjblk < iblk) {
				for (k=piablkstr[jjblk];k<piablkstr[jjblk+1];k++) {
					kkblk = pjablkstr[k];
					if (kkblk >= iblk) {
						if (imaskflt[kkblk] != icycleflt) {
							listrow[nlistrow] = kkblk;
							nlistrow++;
							imaskflt[kkblk] = icycleflt;
						};
					};
				};
			};
		};

		nlistcol = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk <= iblk) {
				listcol[nlistcol] = jjblk;
				nlistcol++;
			};
		};

		((CFctDiagR *)pfctdiag)->Ilu2LoadDiagonalData2Index (iblk, nlistrow, listrow, nlistcol, listcol, 
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

// Update Schur data

		if (nlistupd > 0) {

			((CFctR *)pfct)->Ilu2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, *pparam,
														*((CGSMatrixR *)pgmtrlschr), 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk],
														nlistupd, listupd, 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr(), ((CGSMatrixR *)pgmtruschr)->GetMtrarr(),
														mtrltemp, mtrutemp);

			((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrltemp;
			((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrutemp;

			mtrltemp = mtrdummy;
			mtrutemp = mtrdummy;

		};

// Unload diagonal data from fct structure

		((CFctDiagR *)pfctdiag)->Ilu2UnloadDiagonalData2Index (iblk, nlistrow, listrow,
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

	};

	delete [] listupd;
	delete [] listrow;
	delete [] listcol;

};

// Author: Kharchenko S.A.
// CFctSymmSchurR: Update Schur data for the set of blocks
//========================================================================================
void CFctSymmSchurR::UpdateSchurBlocks (int _nlistfct, int *_listfct, int _nlistschur, int *_listschur) { // Update Schur data for the set of blocks

	const char *funcname = "UpdateSchurBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

// Allocate work array

	int *listupd;
	int *listrow;
	int *listcol;

	listupd = new int [nblks];
	if (!listupd) MemoryFail (funcname);
	listrow = new int [nblks];
	if (!listrow) MemoryFail (funcname);
	listcol = new int [nblks];
	if (!listcol) MemoryFail (funcname);

// Main updates cycle

	int *piablkstr = pablkstr->GetIa ();
	int *pjablkstr = pablkstr->GetJa ();

	int i, iblk;
	int j, jjblk, nlistupd, nlistrow, nlistcol, k, kkblk;

	CSMatrixR mtrltemp;
	CSMatrixR mtrdummy;

	for (i=0;i<_nlistschur;i++) {

		iblk = _listschur[i];

// Register blocks from the list

		icycleflt++;
		for (j=0;j<_nlistfct;j++) {
			jjblk = _listfct[j];
			imaskflt[jjblk] = icycleflt;
		};

// Create the updates list

		nlistupd = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk < iblk) {
				if (imaskflt[jjblk] == icycleflt) {
					listupd[nlistupd] = jjblk;
					nlistupd++;
				};
			};
		};

// Load diagonal data into fct structure

		icycleflt++;

		nlistrow = 0;

		if (imaskflt[iblk] != icycleflt) {
			listrow[nlistrow] = iblk;
			nlistrow++;
			imaskflt[iblk] = icycleflt;
		};

		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				if (imaskflt[jjblk] != icycleflt) {
					listrow[nlistrow] = jjblk;
					nlistrow++;
					imaskflt[jjblk] = icycleflt;
				};
			} else if (jjblk < iblk) {
				for (k=piablkstr[jjblk];k<piablkstr[jjblk+1];k++) {
					kkblk = pjablkstr[k];
					if (kkblk >= iblk) {
						if (imaskflt[kkblk] != icycleflt) {
							listrow[nlistrow] = kkblk;
							nlistrow++;
							imaskflt[kkblk] = icycleflt;
						};
					};
				};
			};
		};

		nlistcol = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk <= iblk) {
				listcol[nlistcol] = jjblk;
				nlistcol++;
			};
		};

		((CFctDiagR *)pfctdiag)->Ilu2LoadDiagonalData2Index (iblk, nlistrow, listrow, nlistcol, listcol, 
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

// Update Schur data

		if (nlistupd > 0) {

			((CFctR *)pfct)->Ich2UpdateBlkRowBySetOfBlocksGlobal2Index (iblk, *pparam,
														*((CGSMatrixR *)pgmtrlschr), 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk],
														nlistupd, listupd, 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr(),
														mtrltemp);

			((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrltemp;

			mtrltemp = mtrdummy;

		};

// Unload diagonal data from fct structure

		((CFctDiagR *)pfctdiag)->Ilu2UnloadDiagonalData2Index (iblk, nlistrow, listrow,
															*((CGSMatrixR *)pgmtral), *((CFctR *)pfct));

	};

	delete [] listupd;
	delete [] listrow;
	delete [] listcol;

};

// Author: Kharchenko S.A.
// CFctSchurR: Filter the set of Schur blocks
//========================================================================================
void CFctSchurR::FilterSchurBlocks (int _nlistschur, int *_listschur) { // Filter the set of Schur blocks

	const char *funcname = "FilterSchurBlocks";

	int *listrow;

	listrow = new int [nblks];
	if (!listrow) MemoryFail (funcname);

	int *piablkstr = pablkstr->GetIa ();
	int *pjablkstr = pablkstr->GetJa ();

	int i, iblk, nlistrow, j, jjblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		nlistrow = 0;
		for (j=piablkstr[iblk];j<piablkstr[iblk+1];j++) {
			jjblk = pjablkstr[j];
			if (jjblk >= iblk) {
				listrow[nlistrow] = jjblk;
				nlistrow++;
			};
		};
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk].FilterBlockRow2Index (true, nlistrow, listrow, icycleflt, imaskflt);
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk].FilterBlockRow2Index (true, nlistrow, listrow, icycleflt, imaskflt);
	};

	delete [] listrow;

};

// Author: Kharchenko S.A.
// CFctSchurR: Free Schur data for the set of blocks
//========================================================================================
void CFctSchurR::FreeSchurBlocks (int _nlistschur, int *_listschur) { // Free Schur data for the set of blocks

//	const char *funcname = "FreeSchurBlocks";

	CSMatrixR mtrdummy;

	int i, iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrdummy;
	};

};

// Author: Kharchenko S.A.
// CFctSchurR: Free Schur data for the set of blocks
//========================================================================================
void CFctSchurR::AddSchurBlocks (int _nlistschur, int *_listschur) { // Add the set of Schur blocks

//	const char *funcname = "AddSchurBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

	CSMatrixR mtrl, mtru;
	CSMatrixR mtrdummy;

	int i, iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		CSMatrixR::DoubleAddMatrSS2Index (true,
														((CGSMatrixR *)pgmtrladd)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk], mtrl,
														((CGSMatrixR *)pgmtruadd)->GetMtrarr()[iblk], ((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk], mtru);
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrl;
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtru;
		((CGSMatrixR *)pgmtrladd)->GetMtrarr()[iblk] = mtrdummy;
		((CGSMatrixR *)pgmtruadd)->GetMtrarr()[iblk] = mtrdummy;
	};

};

// Author: Kharchenko S.A.
// CFctSymmSchurR: Free Schur data for the set of blocks
//========================================================================================
void CFctSymmSchurR::AddSchurBlocks (int _nlistschur, int *_listschur) { // Add the set of Schur blocks

//	const char *funcname = "AddSchurBlocks";

//	char strbuff[256];
//	sprintf (strbuff, "%s%i%s","ChkBlk_",ptree->GetMyid(),".dat");
//	ofstream fout (strbuff,ios::app);

	CSMatrixR mtrl;
	CSMatrixR mtrdummy;

	int i, iblk;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		mtrl = ((CGSMatrixR *)pgmtrladd)->GetMtrarr()[iblk].AddMatrSS2Index (true, ((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk]);
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrl;
		((CGSMatrixR *)pgmtrladd)->GetMtrarr()[iblk] = mtrdummy;
	};

};

// Author: Kharchenko S.A.
// CFctSchurR: Send the set of Schur blocks
//========================================================================================
void CFctSchurR::SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur) { // Send the set of Schur blocks

	const char *funcname = "SendSchurBlocks";

// Combine all L data into one matrix

	int i, iblk;
	int lenobjl, lenobju, lenobjd, lenobj, lenobj1;
	int *piarr;
	char *pobjl, *pobju, *pobjd, *pobj, *ptr;

	CSMatrixR mtrlsend, mtrusend, mtrdsend;

	CSMatrixR mtrdummy;

	((CGSMatrixR *)pgmtrlschr)->CombineMatrices2Index (true, _nlistschur, _listschur, 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr(), mtrlsend);

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
	};

	((CGSMatrixR *)pgmtruschr)->CombineMatrices2Index (true, _nlistschur, _listschur, 
														((CGSMatrixR *)pgmtruschr)->GetMtrarr(), mtrusend);

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtruschr)->GetMtrarr()[iblk] = mtrdummy;
	};

	CSMatrixR::SuperSparsifySymm2Index (true, mtrlsend, mtrusend);

	mtrlsend.PackMatrix (lenobjl,pobjl);
	mtrlsend = mtrdummy;

	mtrusend.SetNlist (0);
	mtrusend.SetNlist2 (0);
	mtrusend.SetNzja (0);
	mtrusend.SetNzja2 (0);
	mtrusend.SetNzja3 (0);
	mtrusend.PackMatrix (lenobju,pobju);
	mtrusend = mtrdummy;

	mtrdsend = ((CGSMatrixR *)pgmtrlschr)->CreateDiagonalUpdate2Index (_nlistschur, _listschur, 
																		*((CFctDiagR *)pfctdiag));

	mtrdsend.SetNzja (0);
	mtrdsend.SetNzja2 (0);
	mtrdsend.SetNzja3 (0);
	mtrdsend.PackMatrix (lenobjd,pobjd);

	mtrdsend = mtrdummy;

	lenobj = 3*sizeof(int) + lenobjl + lenobju + lenobjd;

	pobj = new char [lenobj];
	if(!pobj) MemoryFail(funcname);

	ptr = pobj;

	piarr = (int *) pobj;

	piarr[0] = lenobjl;
	piarr[1] = lenobju;
	piarr[2] = lenobjd;

	pobj += 3*sizeof(int);

	memcpy (pobj,pobjl,lenobjl);
	pobj += lenobjl;
	memcpy (pobj,pobju,lenobju);
	pobj += lenobju;
	memcpy (pobj,pobjd,lenobjd);
	pobj += lenobjd;

	delete [] pobjl;
	delete [] pobju;
	delete [] pobjd;

// Two sequential sends

	if (nsends+1 >= nsendsmax) throw " CGSMatrixR::AsyncSendBlocks2Index: to small number of sends reserved ";

	lenobj1 = sizeof(int);

	pobj = new char [lenobj1];
	if(!pobj) MemoryFail(funcname);

	piarr = (int *) pobj;
	piarr[0] = lenobj;

	imasksend[nsends] = 1;
	psenddata[nsends] = pobj;

	CMPIExchange::ISend (ptree->GetComm (), _jproc, 0,
								lenobj1, pobj, sendreqvarr[nsends]);

	nsends++;

	imasksend[nsends] = 1;
	psenddata[nsends] = ptr;

	CMPIExchange::ISend (ptree->GetComm (), _jproc, 1,
								lenobj, ptr, sendreqvarr[nsends]);

	nsends++;

};

// Author: Kharchenko S.A.
// CFctSymmSchurR: Send the set of Schur blocks
//========================================================================================
void CFctSymmSchurR::SendSchurBlocks (int _jproc, int _nlistschur, int *_listschur) { // Send the set of Schur blocks

	const char *funcname = "SendSchurBlocks";

// Combine all L data into one matrix

	int i, iblk;
	int lenobjl, lenobjd, lenobj, lenobj1;
	int *piarr;
	char *pobjl, *pobjd, *pobj, *ptr;

	CSMatrixR mtrlsend, mtrdsend;

	CSMatrixR mtrdummy;

	((CGSMatrixR *)pgmtrlschr)->CombineMatrices2Index (true, _nlistschur, _listschur, 
														((CGSMatrixR *)pgmtrlschr)->GetMtrarr(), mtrlsend);

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		((CGSMatrixR *)pgmtrlschr)->GetMtrarr()[iblk] = mtrdummy;
	};

	CSMatrixR::SuperSparsifyL2Index (true, mtrlsend);

	mtrlsend.PackMatrix (lenobjl,pobjl);
	mtrlsend = mtrdummy;

	mtrdsend = ((CGSMatrixR *)pgmtrlschr)->CreateDiagonalUpdate2Index (_nlistschur, _listschur, 
																		*((CFctDiagR *)pfctdiag));

	mtrdsend.SetNzja (0);
	mtrdsend.SetNzja2 (0);
	mtrdsend.SetNzja3 (0);
	mtrdsend.PackMatrix (lenobjd,pobjd);

	mtrdsend = mtrdummy;

	lenobj = 2*sizeof(int) + lenobjl + lenobjd;

	pobj = new char [lenobj];
	if(!pobj) MemoryFail(funcname);

	ptr = pobj;

	piarr = (int *) pobj;

	piarr[0] = lenobjl;
	piarr[1] = lenobjd;

	pobj += 2*sizeof(int);

	memcpy (pobj,pobjl,lenobjl);
	pobj += lenobjl;
	memcpy (pobj,pobjd,lenobjd);
	pobj += lenobjd;

	delete [] pobjl;
	delete [] pobjd;

// Two sequential sends

	if (nsends+1 >= nsendsmax) throw " CFctSymmSchurR::SendSchurBlocks: to small number of sends reserved ";

	lenobj1 = sizeof(int);

	pobj = new char [lenobj1];
	if(!pobj) MemoryFail(funcname);

	piarr = (int *) pobj;
	piarr[0] = lenobj;

	imasksend[nsends] = 1;
	psenddata[nsends] = pobj;

	CMPIExchange::ISend (ptree->GetComm (), _jproc, 0,
								lenobj1, pobj, sendreqvarr[nsends]);

	nsends++;

	imasksend[nsends] = 1;
	psenddata[nsends] = ptr;

	CMPIExchange::ISend (ptree->GetComm (), _jproc, 1,
								lenobj, ptr, sendreqvarr[nsends]);

	nsends++;

};

// Author: Kharchenko S.A.
// CFctSchurR: Receive the set of Schur blocks
//========================================================================================
void CFctSchurR::ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur) { // Receive the set of Schur blocks

	const char *funcname = "ReceiveSchurBlocks";

	int isize = sizeof(int);

	char *buffer;

	buffer = new char [isize];
	if (!buffer) MemoryFail (funcname);

	CMPIStatus status;

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 0,
								isize, buffer, status);

	int *pint;

	pint = (int *)buffer;
	isize = pint[0];

	delete [] buffer;

	buffer = new char [isize];
	if (!buffer) MemoryFail (funcname);

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 1,
								isize, buffer, status);

	CSMatrixR mtrl, mtru, mtrd;

	int *piarr = (int*) buffer;

	int isizel, isizeu, isized;

	isizel = piarr[0];
	isizeu = piarr[1];
	isized = piarr[2];

	char *pbuff = buffer + 3*sizeof(int);

	mtrl.UnPackMatrix(isizel,pbuff);

	pbuff += isizel;

	mtru.UnPackMatrix(isizeu,pbuff);

	pbuff += isizeu;

	mtrd.UnPackMatrix(isized,pbuff);

	delete [] buffer;

// Copy the L sparsity into U

	((CSMatrix &)mtru) = ((CSMatrix &)mtrl);

// Update all necessary data

	int nlistloc;
	int *listloc;
	CSMatrixR mtrdummy;

	listloc = new int [nblks];
	if (!listloc) MemoryFail (funcname);

	((CGSMatrixR *)pgmtrlschr)->SplitMatrixList2Index (nlistloc, listloc,
												mtrl, 
												((CGSMatrixR *)pgmtrladd)->GetMtrarr());
	((CGSMatrixR *)pgmtruschr)->SplitMatrixList2Index (nlistloc, listloc,
												mtru, 
												((CGSMatrixR *)pgmtruadd)->GetMtrarr());

	mtrl = mtrdummy;
	mtru = mtrdummy;

	((CGSMatrixR *)pgmtrlschr)->UpdateDiagonal2Index (mtrd, *((CFctDiagR *)pfctdiag));
	mtrd = mtrdummy;

	delete [] listloc;

};

// Author: Kharchenko S.A.
// CFctSchurR: Receive the set of Schur blocks
//========================================================================================
void CFctSymmSchurR::ReceiveSchurBlocks (int _jproc, int _nlistschur, int *_listschur) { // Receive the set of Schur blocks

	const char *funcname = "ReceiveSymmSchurBlocks";

	int isize = sizeof(int);

	char *buffer;

	buffer = new char [isize];
	if (!buffer) MemoryFail (funcname);

	CMPIStatus status;

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 0,
								isize, buffer, status);

	int *pint;

	pint = (int *)buffer;
	isize = pint[0];

	delete [] buffer;

	buffer = new char [isize];
	if (!buffer) MemoryFail (funcname);

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 1,
								isize, buffer, status);

	CSMatrixR mtrl, mtru, mtrd;

	int *piarr = (int*) buffer;

	int isizel, isized;

	isizel = piarr[0];
	isized = piarr[1];

	char *pbuff = buffer + 2*sizeof(int);

	mtrl.UnPackMatrix(isizel,pbuff);

	pbuff += isizel;

	mtrd.UnPackMatrix(isized,pbuff);

	delete [] buffer;

// Update all necessary data

	int nlistloc;
	int *listloc;
	CSMatrixR mtrdummy;

	listloc = new int [nblks];
	if (!listloc) MemoryFail (funcname);

	((CGSMatrixR *)pgmtrlschr)->SplitMatrixList2Index (nlistloc, listloc,
												mtrl, 
												((CGSMatrixR *)pgmtrladd)->GetMtrarr());

	mtrl = mtrdummy;

	((CGSMatrixR *)pgmtrlschr)->UpdateDiagonal2Index (mtrd, *((CFctDiagR *)pfctdiag));
	mtrd = mtrdummy;

	delete [] listloc;

};

// Author: Kharchenko S.A.
// CFctSchurR: Wait for completion of sends
//========================================================================================
void CFctSchurR::WaitSends () { // Wait for completion of sends

	const char *funcname = "WaitSends";

	int nloc = nsends;

	CMPIRequest *reqvarrloc;
	CMPIStatus *statarrloc;

	reqvarrloc = new CMPIRequest [nloc];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new CMPIStatus [nloc];
	if (!statarrloc) MemoryFail (funcname);

	nloc = 0;

	int i;

	for (i=0;i<nsends;i++) {
		if (imasksend[i] == 1) {
			reqvarrloc[nloc] = sendreqvarr[i];
			nloc++;
		};
	};

	CMPIExchange::WaitAll (nloc, reqvarrloc, statarrloc);

	for (i=0;i<nsends;i++) {
		if (imasksend[i] == 1) {
			delete [] psenddata[i];
		};
	};

	delete [] reqvarrloc;
	delete [] statarrloc;

	nsends = 0;

};

// Author: Kharchenko S.A.
// CFctSchurR: Destroy fct working data
//========================================================================================
void CFctSchurR::DestroyWorkData () { // Destroy fct working data

	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Compute the global hystogram of diagonal pivots and output it

	((CFctR *)pfct)->GetHyst()->MakeGlobalHystogram (ptree->GetComm ());

	if (ptree->GetMyid () == ptree->GetRootcpu ()) {
		if (pparam->msglev > 1) {
			cout  << " Factorization" << *(((CFctR *)pfct)->GetHyst());
		};
		if (pparam->msglev > 0) {
			*pfout  << " Factorization" << *(((CFctR *)pfct)->GetHyst());
		};
	};

// Compute density ILU2 factorization statistics

	int iblk;

	double nzatotal = 0;
	for (iblk=0;iblk<nblks;iblk++) {
		nzatotal += ((CGSMatrixR *)pgmtral)->GetMtrarr()[iblk].GetNzatot ();
	};

	double nzmema = 2 * nzatotal;

	double nzutot = 0;
	for (iblk=0;iblk<nblks;iblk++) {
		nzutot += ((CGSMatrixR *)pgmtrl)->GetMtrarr()[iblk].GetNzatot ();
	};

	nzutot *= 2.0e0;

	double nzrtot = ((CFctR *)pfct)->GetNzr ();

	nzrtot *= 2.0e0;

// Compute global data

	double arr0[7], arr1[7];

	arr0[0] = (double) nzutot;
	arr0[1] = (double) nzrtot;
	arr0[2] = 0;
	arr0[3] = (double) nzatotal;
	arr0[4] = (double) nzmema;
	arr0[5] = ((CFctR *)pfct)->GetOps ();
	arr0[6] = (double) ((CFctR *)pfct)->GetNmodsv ();

	if (nproc != 1) {
		CMPIExchange::ExchangeArrayMPI (ptree->GetComm (),
													DOUBLEVALUE, ADD,
													7, arr0, arr1);
	} else {
		for (int i=0;i<=6;i++) arr1[i] = arr0[i];
	};

	double dnzu = arr1[0];
	double dnzr = arr1[1];
//	double dnzs = arr1[2];
	double dnzamem = arr1[4];
	double ops = arr1[5];
	int nmodtot = (int) arr1[6];

// Output resulting statistics

	double densu = 1.0e2 * dnzu / dnzamem;
	double densr = 1.0e2 * dnzr / dnzamem;
//	double denss = 1.0e2 * dnzs / dnzamem;

	ops = ops / dnzamem;

	if (myid == ptree->GetRootcpu()) {

		if (pparam->msglev > 1) {
			cout  << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				cout  << "     NmodS = " << nmodtot << endl;
			};
			cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
//			cout  << "     Costs = " << ops << " MvmA flops. " << endl;
//			cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

		if (pparam->msglev > 0) {
			*pfout << " Ilu2 preconditioner generation statistics: " << endl;
			if (nmodtot != 0) {
				*pfout << "     NmodS = " << nmodtot << endl;
			};
			*pfout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
//			_fout << "     Costs = " << ops << " MvmA flops. " << endl;
//			_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
		};

	};

	delete ((CFct *)pfct);
	delete ((CFctDiag *)pfctdiag);

};


#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif// 
