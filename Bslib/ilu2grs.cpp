//------------------------------------------------------------------------------------------------
// File: ilu2grs.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
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
// CGSMatrixRS: Perform explicit point scaling of the matrix
//========================================================================================
void CGSMatrixRS::PointScaling (const CSVector &_sclpt) { // Perform explicit point scaling of the matrix

	const char *funcname = "PointScaling";

// Prepare all necessary data

	int *ibsscl;

	ibsscl = new int [nsupr];
	if (!ibsscl) MemoryFail (funcname);

	int i, j, ilist, iblk, isup, ni;

	for (i=0;i<nsupr;i++) ibsscl[i] = -1;

	int iblkscl = 0;
	int nloc = 0;

	for (i=0;i<_sclpt.nparts;i++) {
		int icheck = 0;
		for (j=iblkscl;j<nblksr;j++) {
			if (_sclpt.inda[i] == bl2ndr[j]) {
				icheck = 1;
				iblkscl = j;
				break;
			};
		};
		if (icheck == 0) throw " PointScaling: required block is not found ";
		for (isup=blksr[iblkscl];isup<blksr[iblkscl+1];isup++) {
			ibsscl[isup] = nloc + sprndr[isup]-sprndr[blksr[iblkscl]];
		};
		ni = _sclpt.blksv[i+1]-_sclpt.blksv[i];
		nloc += ni;
	};

// Main cycle over the block rows

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		ReadBlock (iblk);

		mtrarr[iblk].PointScaling (sprndr, ibsscl, _sclpt.vect);

		ReWriteBlock (iblk);
		FreeBlock    (iblk);
	};

// Free work arrays

	delete [] ibsscl;

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Perform explicit block scaling of the matrix
//========================================================================================
void CGSMatrixRS::BlockScaling (const CSVector &_sclblkl, const CSVector &_sclblku) { // Perform explicit block scaling of the matrix

	const char *funcname = "BlockScaling";

// Prepare all necessary data

	int *ibsscl;

	ibsscl = new int [nsupr];
	if (!ibsscl) MemoryFail (funcname);

	int i, j, ilist, iblk, isup, niloc;

	for (i=0;i<nsupr;i++) ibsscl[i] = -1;

	int iblkscl = 0;
	int nloc = 0;

	for (i=0;i<_sclblkl.nparts;i++) {
		int icheck = 0;
		for (j=iblkscl;j<nblksr;j++) {
			if (_sclblkl.inda[i] == bl2ndr[j]) {
				icheck = 1;
				iblkscl = j;
				break;
			};
		};
		if (icheck == 0) throw " BlockScaling: required block is not found ";
		for (isup=blksr[iblkscl];isup<blksr[iblkscl+1];isup++) {
			niloc = sprndr[isup+1]-sprndr[isup];
			ibsscl[isup] = nloc;
			nloc += niloc*niloc;
		};
	};

// Main cycle over the block rows

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = listb[ilist];
		ReadBlock (iblk);

		mtrarr[iblk].BlockScaling (sprndr, ibsscl, _sclblkl.vect, _sclblku.vect);

		ReWriteBlock (iblk);
		FreeBlock    (iblk);
	};

// Free work arrays

	delete [] ibsscl;

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Second order Incomplete LU decomposition
//========================================================================================
void CGSMatrixRS::Ilu2 (ofstream &_fout, const CSlvParam &_param, // Second order Incomplete LU decomposition
						const CGSMatrixRS &_mtral, const CGSMatrixRS &_mtrau,
						CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru) { 

	const char *funcname = "Ilu2";

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

	pmtral->MtrBlockSparsity (pmtrarr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsity (_mtral.nblksr, iab, jab, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	int nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
						iabextt, jabextt);

// Determine blamx

	int blamxl=0;
	int blai;

	int isup;
	for (isup=0;isup<_mtral.nsupr;isup++) {
		blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
		if (blai > blamxl) blamxl = blai;
	};

// Allocate Fct data

	CFctR fct (_mtral.nblksr, _mtral.nsupr, blamxl); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr); 

// Compute the scaling of the matrix

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		fctdiag.Ilu2ScaleBlkRow (iblk, _param.sclpttype, _param.sclmin,
									_mtral, _mtrau.mtrarr[iblk], fct);
	};

//	csvhyst (cout, "ScDia",nloc, _dpiv);
//	csvhyst (_fout,"ScDia",nloc, _dpiv);

// Create dummy matrix

	CSMatrixRS mtrdummy;

// Allocate array containing second order matrices

	CSMatrixRS *mtrl2arr, *mtru2arr;

	mtrl2arr = new CSMatrixRS [_mtral.nblksr];
	if (!mtrl2arr) MemoryFail (funcname);
	mtru2arr = new CSMatrixRS [_mtral.nblksr];
	if (!mtru2arr) MemoryFail (funcname);

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

// Main cycle over the block rows

	CSMatrixRS mtrltemp, mtrl2temp;
	CSMatrixRS mtrutemp, mtru2temp;

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {

// Load diagonal data into fct structure

		int nlstblk = iabext[iblk+1]-iabext[iblk];
		int ibs = iabext[iblk];

		fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
										_mtral, fct);

// Init current block row of the matrix

		fct.Ilu2InitBlkRow (iblk, _param.dshift, 
							_mtral, _mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
							mtrltemp, mtrutemp);

// Update current block row by the previous ones

		for (int j=iabextt[iblk];j<iabextt[iblk+1]-1;j++) {

			int jblk = jabextt[j];

//			_fout << " update Jblk = " << jblk << endl;

			fct.Ilu2UpdateBlkRow (iblk, _param,
									_mtral, mtrltemp, mtrutemp,
									mtrl2arr[jblk], mtru2arr[jblk],
									mtrl2temp, mtru2temp);

			mtrltemp = mtrl2temp;
			mtrutemp = mtru2temp;

			mtrl2temp = mtrdummy;
			mtru2temp = mtrdummy;

		};

// Factorize current block row

		fct.Ilu2FctBlkRow (iblk, _param, _mtral.nsupr, 
							_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
							_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
							mtrl2arr[iblk], mtru2arr[iblk]);

		mtrltemp = mtrdummy;
		mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

		nlstblk = iabext[iblk+1]-iabext[iblk];
		ibs = iabext[iblk];

		fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
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

// Output Ich2 factorization statistics

	int nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	int nzmema = 2 * nzatotal;

	for (isup=0;isup<_mtral.nsupr;isup++) {
		blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
		nzmema -= blai*blai;
	};

	int nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	int nzrtot = fct.nzr;

	densu = 1.0e2 * (double) nzutot / (double) nzatotal;
	densr = 1.0e2 * (double) nzrtot / (double) nzatotal;

	ops = fct.ops;

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
// CGSMatrixRS: Compute second order Incomplete LU decomposition in parallel mode
//========================================================================================
void CGSMatrixRS::Ilu2 (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
						CGSMatrixRS &_mtral, CGSMatrixRS &_mtrau,
						CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru) { 

	const char *funcname = "Ilu2";

	CMPIComm pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);

// Check the data on entry

	int numprocs;

	numprocs = pcomm.GetNproc ();

	if (_tree.nproc != numprocs) {
		throw " Parallel Ilu2 function called with wrong number of processors ";
	};

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timescl=0.0e0, timesclex=0.0e0;
	double timeload=0.0e0, timeunload=0.0e0;
	double timeini=0.0e0, timeupdlocal=0.0e0, timeupdnonlocal=0.0e0, timefct=0.0e0;
	double timeprep=0.0e0, timesend=0.0e0, timerecv=0.0e0;

// Init time measurement

	time0 = clock ();

// Compute bl2cpu

	int *bl2cpu;

	bl2cpu = new int [_mtral.nblksr];
	if (!bl2cpu) MemoryFail (funcname);

	_tree.Block2Cpu (bl2cpu);

	int nchildsmx = 0;

	for (int inode=0;inode<_tree.nnodes;inode++) {
//		int iproc = _tree.nodes[inode].nodecpu;
		int nchilds = _tree.nodes[inode].nchilds;
		if (nchilds > nchildsmx) nchildsmx = nchilds;
	};

	nchildsmx++;

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

	pmtral->MtrBlockSparsity (_tree, bl2cpu, pmtrarr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsity (_mtral.nblksr, 
									iab, jab, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	int nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
						iabextt, jabextt);

// Determine blamx

	int blamxl=0;
	int blai;

	int isup;
	for (isup=0;isup<_mtral.nsupr;isup++) {
		blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
		if (blai > blamxl) blamxl = blai;
	};

// Allocate Fct data

	CFctR fct (_mtral.nblksr, _mtral.nsupr, blamxl); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr);

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"DiagLU_",_tree.myid,"_");

		fctdiag.SetupFiles (_param.ooctyp, strbuff);

	};

// Create mvm structure

	CMvm mvm (false, _tree, _mtral);

// Compute the local scaling of the matrix

	time2 = clock ();

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		if (bl2cpu[iblk] == _tree.myid) {

			_mtrau.ReadBlock (iblk);

			fctdiag.Ilu2ScaleBlkRow (iblk, _param.sclpttype, _param.sclmin,
										_mtral, _mtrau.mtrarr[iblk], fct);

			_mtrau.FreeBlock (iblk);

		};
	};

	time3 = clock ();

	timescl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

//	csvhyst (cout, "ScDia",nloc, _dpiv);
//	csvhyst (_fout,"ScDia",nloc, _dpiv);

// Exchange the scaling data from the bordering

	time2 = clock ();

	fctdiag.ExchangeScaling (false, _tree, mvm, _mtral);

	time3 = clock ();

	timesclex += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Create dummy matrix

	CSMatrixRS mtrdummy;

// Allocate array containing second order matrices

	CSMatrixRS *mtrl2charr, *mtru2charr;

	mtrl2charr = new CSMatrixRS [nchildsmx];
	if (!mtrl2charr) MemoryFail (funcname);
	mtru2charr = new CSMatrixRS [nchildsmx];
	if (!mtru2charr) MemoryFail (funcname);

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"LMatr_",_tree.myid,"_");

		_mtrl.SetupFiles (_param.ooctyp, strbuff);

		sprintf(strbuff, "%s%s%d%s",_param.path,"UMatr_",_tree.myid,"_");

		_mtru.SetupFiles (_param.ooctyp, strbuff);

	};

// Init gmtrl2 and gmtru2

	CGSMatrixRS gmtrl2, gmtru2;

	gmtrl2 = _mtral;
	gmtru2 = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"L2Matr_",_tree.myid,"_");

		gmtrl2.SetupFiles (_param.ooctyp, strbuff);

		sprintf(strbuff, "%s%s%d%s",_param.path,"U2Matr_",_tree.myid,"_");

		gmtru2.SetupFiles (_param.ooctyp, strbuff);

	};

// Main cycle over the block rows according to the tree

	CSMatrixRS mtrltemp, mtrl2temp, mtrl2send;
	CSMatrixRS mtrutemp, mtru2temp, mtru2send;

	int inodecurr, fatheridcurr, fathercpucurr, ichild, childcpucurr, childnodecurr;

	inodecurr = _tree.cpuidend[_tree.myid];

	while (inodecurr >= 0) {

// Load diagonal data 

//		cout << " Current node:" << _tree.nodes[inodecurr] << endl;
//		_fout << " Current node:" << _tree.nodes[inodecurr] << endl;

		time2 = clock ();

		fctdiag.Ilu2LoadDiagonalData (inodecurr, _tree, mvm, iabext, jabext,
										_mtral, fct);

		time3 = clock ();

		timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Perform computations with blocks the current node taking into account child data

		for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {

			if (_param.msglev > 1) {
				cout << " Iblk = " << iblk << endl;
			};
//			_fout << " Iblk = " << iblk << endl;

// Load diagonal data into fct structure

			int nlstblk = iabext[iblk+1]-iabext[iblk];
			int ibs = iabext[iblk];

			time2 = clock ();

			fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
											_mtral, fct);

			time3 = clock ();

			timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Init current block row of the matrix

			time2 = clock ();

			_mtral.ReadBlock (iblk);
			_mtrau.ReadBlock (iblk);

			fct.Ilu2InitBlkRow (iblk, _param.dshift,
								_mtral, _mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
								mtrltemp, mtrutemp);

			_mtral.FreeBlock (iblk);
			_mtrau.FreeBlock (iblk);

			time3 = clock ();

			timeini += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Update current block row by the child data

			if (_tree.nodes[inodecurr].nchilds > 1) {

				for (int j=0;j<_tree.nodes[inodecurr].nchilds;j++) {

//					cout << " update child J = " << j << endl;
//					_fout << " update child J = " << j << endl;

					time2 = clock ();

					fct.Ilu2UpdateBlkRow (iblk, _param,
											_mtral, mtrltemp, mtrutemp,
											mtrl2charr[j], mtru2charr[j],
											mtrl2temp, mtru2temp);

					time3 = clock ();

					timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					mtrltemp = mtrl2temp;
					mtrutemp = mtru2temp;

					mtrl2temp = mtrdummy;
					mtru2temp = mtrdummy;

				};

			};

// Update current block row by the previous ones

			if (_tree.nodes[inodecurr].nchilds == 1) {

				for (int j=iabextt[iblk];j<iabextt[iblk+1]-1;j++) {

					int jblk = jabextt[j];

//					cout << " update Jblk = " << jblk << endl;
//					_fout << " update Jblk = " << jblk << endl;

					gmtrl2.ReadBlock (jblk);
					gmtru2.ReadBlock (jblk);

					time2 = clock ();

					fct.Ilu2UpdateBlkRow (iblk, _param,
											_mtral, mtrltemp, mtrutemp,
											gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
											mtrl2temp, mtru2temp);

					time3 = clock ();

					timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					gmtrl2.FreeBlock (jblk);
					gmtru2.FreeBlock (jblk);

					mtrltemp = mtrl2temp;
					mtrutemp = mtru2temp;

					mtrl2temp = mtrdummy;
					mtru2temp = mtrdummy;

				};

			} else {

				for (int jblk=_tree.nodes[inodecurr].indbeg;jblk<iblk;jblk++) {

					gmtrl2.ReadBlock (jblk);
					gmtru2.ReadBlock (jblk);

					time2 = clock ();

					fct.Ilu2UpdateBlkRow (iblk, _param,
											_mtral, mtrltemp, mtrutemp,
											gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
											mtrl2temp, mtru2temp);

					time3 = clock ();

					timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					gmtrl2.FreeBlock (jblk);
					gmtru2.FreeBlock (jblk);

					mtrltemp = mtrl2temp;
					mtrutemp = mtru2temp;

					mtrl2temp = mtrdummy;
					mtru2temp = mtrdummy;

				};

			};

// Factorize current block row

//			cout << " Before fct " << endl;

			time2 = clock ();

			fct.Ilu2FctBlkRow (iblk, _param, _mtral.nsupr, 
								_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk], 
								_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
								gmtrl2.mtrarr[iblk], gmtru2.mtrarr[iblk]);

			time3 = clock ();

			timefct += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			_mtrl.WriteBlock (iblk);
			_mtrl.FreeBlock  (iblk);

			_mtru.WriteBlock (iblk);
			_mtru.FreeBlock  (iblk);

			gmtrl2.WriteBlock (iblk);
			gmtrl2.FreeBlock  (iblk);

			gmtru2.WriteBlock (iblk);
			gmtru2.FreeBlock  (iblk);

			mtrltemp = mtrdummy;
			mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

			nlstblk = iabext[iblk+1]-iabext[iblk];
			ibs = iabext[iblk];

			time2 = clock ();

			fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
											_mtral, fct);

			time3 = clock ();

			timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		};

// Prepare the data to be exchanged if necessary

		if (_tree.nodes[inodecurr].fatherid != inodecurr) {

			time2 = clock ();

			_mtral.PrepareSendMatrix (inodecurr, _tree,
										gmtrl2, mtrl2charr,
										mtrl2send);
			_mtral.PrepareSendMatrix (inodecurr, _tree,
										gmtru2, mtru2charr,
										mtru2send);

			time3 = clock ();

			timeprep += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

		};

// Free second order matrices of the current node

		for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {
			gmtrl2.mtrarr[iblk] = mtrdummy;
			gmtru2.mtrarr[iblk] = mtrdummy;
		};

// Free child data

		for (int j=0;j<_tree.nodes[inodecurr].nchilds;j++) {
			mtrl2charr[j] = mtrdummy;
			mtru2charr[j] = mtrdummy;
		};

// Exchange the data if necessary

		if (_tree.nodes[inodecurr].fatherid != inodecurr) {

			fatheridcurr = _tree.nodes[inodecurr].fatherid;
			fathercpucurr = _tree.nodes[inodecurr].fathercpu;

// FatherCpu

			if (fathercpucurr == _tree.myid) {

				for (ichild=0;ichild<_tree.nodes[fatheridcurr].nchilds;ichild++) {
					childnodecurr = _tree.nodes[fatheridcurr].childs[ichild];
					childcpucurr = _tree.nodes[childnodecurr].nodecpu;
					if (childcpucurr == _tree.myid) {
						mtrl2charr[ichild] = mtrl2send;
						mtru2charr[ichild] = mtru2send;
					} else {

						time2 = clock ();

//						_mtral.ReceiveMatrix (childcpucurr, mtrl2charr[ichild]);
//						_mtral.ReceiveMatrix (childcpucurr, mtru2charr[ichild]);

						time3 = clock ();

						timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					};
				};

			} else {

// ChildCpu

				time2 = clock ();

//				_mtral.SendMatrix (fathercpucurr, mtrl2send);
//				_mtral.SendMatrix (fathercpucurr, mtru2send);

				time3 = clock ();

				timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			};

		};

// Free send data

		mtrl2send = mtrdummy;
		mtru2send = mtrdummy;

// Unload diagonal data from fct structure

		time2 = clock ();

		fctdiag.Ilu2UnloadDiagonalData (inodecurr, _tree, mvm, iabext, jabext,
										_mtral, fct);

		time3 = clock ();

		timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Go to the previous tree level if necessary

		if (_tree.nodes[inodecurr].fatherid != inodecurr && _tree.nodes[inodecurr].fathercpu == _tree.myid) {
			fatheridcurr = _tree.nodes[inodecurr].fatherid;
			inodecurr = fatheridcurr;
		} else {
			inodecurr = -1;
		};

	};

// Close files

	if (_param.ooctyp > 0) {

		gmtrl2.CloseFiles();
		gmtru2.CloseFiles();

		fctdiag.CloseFiles();

	};

// Finalize time measurement

	pcomm = _tree.GetComm ();
	CMPIExchange::SynchronizeMPI (pcomm);
	
	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output timing results

	if (_param.msglev > 1) {
		cout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
		cout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s." << endl;
		cout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
	};
	if (_param.msglev > 0) {
		_fout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
		_fout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s." << endl;
		_fout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
	};

// Output Ich2 factorization statistics

	int nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	int nzmema = 2 * nzatotal;

	for (isup=0;isup<_mtral.nsupr;isup++) {
		blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
		nzmema -= blai*blai;
	};

	int nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	int nzrtot = fct.nzr;

	// Compute global data

	double arr0[5], arr1[5];

	arr0[0] = (double) nzutot;
	arr0[1] = (double) nzrtot;
	arr0[2] = (double) nzatotal;
	arr0[3] = (double) nzmema;
	arr0[4] = fct.ops;

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													5, arr0, arr1);
	} else {
		for (int i=0;i<5;i++) arr1[i] = arr0[i];
	};

	double dnzu = arr1[0];
	double dnzr = arr1[1];
	double dnza = arr1[2];
	double dnzamem = arr1[3];
	ops = arr1[4];

// Output resulting statistics

	densu = 1.0e2 * dnzu / dnza;
	densr = 1.0e2 * dnzr / dnza;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnzamem;

	if (_tree.myid == _tree.rootcpu) {

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

	};

// Free working arrays

	delete [] bl2cpu;
	delete [] pmtrarr;
	delete [] iab;
	delete [] jab;
	delete [] iabext;
	delete [] jabext;
	delete [] iabextt;
	delete [] jabextt;
	delete [] mtrl2charr;
	delete [] mtru2charr;

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Compute second order Incomplete LU decomposition in parallel mode
//========================================================================================
void CGSMatrixRS::Ilu2Schur (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete LU decomposition in parallel mode
								CGSMatrixRS &_mtral, CGSMatrixRS &_mtrau,
								CGSMatrixRS &_mtrl, CGSMatrixRS &_mtru) { 

	const char *funcname = "Ilu2Schur";

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

// Determine the case

	int icase = 1;
	if (_mtral.n == _mtral.nsupr) icase = 0;

// Check the data on entry

	int numprocs;

	CMPIComm pcomm = _tree.GetComm ();
	numprocs = pcomm.GetNproc ();

	if (_tree.nproc != numprocs) {
		throw " Parallel Ilu2 function called with wrong number of processors ";
	};

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkIlu2Schur_",_tree.myid,".dat");
//	std::ofstream fout (strbuff);

// Statistics data

	double ops, densu, densr, denss, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timescl=0.0e0, timesclex=0.0e0;
	double timeload=0.0e0, timeunload=0.0e0;
	double timeini=0.0e0, timeupdlocal=0.0e0, timeupdnonlocal=0.0e0, timefct=0.0e0, timefctd=0.0e0;
	double timeprep=0.0e0, timesend=0.0e0, timerecv=0.0e0;
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

	_tree.Block2Cpu (bl2cpu);

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

	pmtral->MtrBlockSparsity (_tree, bl2cpu, pmtrarr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtral->ExtendBlockSparsity (_mtral.nblksr, 
									iab, jab, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	int nzjext = iabext[_mtral.nblksr];

	iabextt = new int [_mtral.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtral.nblksr, iabext, jabext,
						iabextt, jabextt);

// Determine blamx

	int blamxl=0;
	int blai;

	int isup;
	for (isup=0;isup<_mtral.nsupr;isup++) {
		blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
		if (blai > blamxl) blamxl = blai;
	};

// Allocate Fct data

	CFctR fct (_mtral.nblksr, _mtral.nsupr, blamxl); 

// Allocate FctDiag data

	CFctDiagR fctdiag (_mtral.nblksr);

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"DiagLU_",_tree.myid,"_");

		fctdiag.SetupFiles (_param.ooctyp, strbuff);

	};

// Create mvm structure

	CMvm mvm (false, _tree, _mtral);

// Compute the local scaling of the matrix

	for (iblk=0;iblk<_mtral.nblksr;iblk++) {

		if (bl2cpu[iblk] == _tree.myid) {

			_mtrau.ReadBlock (iblk);

			time2 = clock ();

			fctdiag.Ilu2ScaleBlkRow (iblk, _param.sclpttype, _param.sclmin,
										_mtral, _mtrau.mtrarr[iblk], fct);

			time3 = clock ();

			timescl += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			_mtrau.FreeBlock (iblk);

		};

	};

// Update and output hystogram

	pcomm = _tree.GetComm ();

	fctdiag.hyst.MakeGlobalHystogram (pcomm);
	fctdiag.hyst2.MakeGlobalHystogram (pcomm);

	if (_tree.myid == _tree.rootcpu) {
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

	time2 = clock ();

	fctdiag.ExchangeScaling (false, _tree, mvm, _mtral);

	time3 = clock ();

	timesclex += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Create dummy matrix

	CSMatrixRS mtrdummy;

// Init mtrl and mtru

	_mtrl = _mtral;
	_mtru = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"LMatr_",_tree.myid,"_");

		_mtrl.SetupFiles (_param.ooctyp, strbuff);

		sprintf(strbuff, "%s%s%d%s",_param.path,"UMatr_",_tree.myid,"_");

		_mtru.SetupFiles (_param.ooctyp, strbuff);

	};

	_mtrl.ivect = _param.sclpttype;
	_mtru.ivect = _param.sclpttype;

// Init gmtrl2 and gmtru2

	CGSMatrixRS gmtrl2, gmtru2;
	CGSMatrixRS gmtrl3, gmtru3;

	gmtrl2 = _mtral;
	gmtru2 = _mtrau;

	gmtrl3 = _mtral;
	gmtru3 = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"L2Matr_",_tree.myid,"_");

		gmtrl2.SetupFiles (_param.ooctyp, strbuff);

		sprintf(strbuff, "%s%s%d%s",_param.path,"U2Matr_",_tree.myid,"_");

		gmtru2.SetupFiles (_param.ooctyp, strbuff);

	};

// Compute the list of blocks for which Schur data should be supported

	int inodecurr, ichild, childcpucurr, childnodecurr;

	int nlistschur;
	int *listschur;

	inodecurr = _tree.cpuidend[_tree.myid];

	int ilev = _tree.nodes[inodecurr].nodelv;

	nlistschur = mvm.nlistfctarr[ilev];
	listschur = mvm.listfctarr[ilev];

// Init Schur data blocks and store them to disk if necessary

	int ilist;
	for (ilist=0;ilist<nlistschur;ilist++) {

		int iblk = listschur[ilist];
		int isupbeg = _mtral.blksc[iblk];
		int isupend = _mtral.blksc[iblk+1]-1;

		gmtrl2.mtrarr[iblk] = fct.CreateBlock (isupbeg, isupend, _mtral);

		if (blamxl == 1) gmtrl2.mtrarr[iblk].CleanBsa ();

		gmtrl2.WriteBlock (iblk);
		gmtrl2.FreeBlock  (iblk);

		gmtru2.mtrarr[iblk] = fct.CreateBlock (isupbeg, isupend, _mtral);

		if (blamxl == 1) gmtru2.mtrarr[iblk].CleanBsa ();

		gmtru2.WriteBlock (iblk);
		gmtru2.FreeBlock  (iblk);

	};

// Main cycle over the block rows according to the tree

	CSMatrixRS mtrltemp, mtrl2temp, mtrlrecv, mtrlsend;
	CSMatrixRS mtrutemp, mtru2temp, mtrurecv, mtrusend;

	inodecurr = _tree.cpuidend[_tree.myid];

	int inodesave = inodecurr;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inodecurr = mvm.lev2node[ilev];
		int inodetype = _tree.nodes[inodecurr].nodetype;

		int fathernode = _tree.nodes[inodecurr].fatherid;
		int fathercpuid = _tree.nodes[inodecurr].fathercpu;
//		int fathernode2 = _tree.nodes[inodecurr].fatherid2;
		int fathercpuid2 = _tree.nodes[inodecurr].fathercpu2;

		if (_tree.nodes[inodecurr].nodecpu != _tree.myid) goto NextNode;

// Receive the data if necessary

		if (inodetype != compnode2) {

// Receive data from childs

			if (_tree.nodes[inodecurr].nchilds != 1) {

				for (ichild=0;ichild<_tree.nodes[inodecurr].nchilds;ichild++) {

					childnodecurr = _tree.nodes[inodecurr].childs[ichild];
					childcpucurr = _tree.nodes[childnodecurr].nodecpu;

					if (childcpucurr != _tree.myid) {

						int nlistloc = mvm.nlistuprcvarr[ilev];
						int *listloc = mvm.listuprcvarr[ilev];
//						cout << " Recv matrices from " << childcpucurr << endl;
//						OutArr(cout," List of recv ",nlistloc,listloc);

						if (_param.ooctyp == 0) {

							if (nlistloc > 0) {

// Receive and split

								time2 = clock ();

								iblk = listloc[0];

								_mtral.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk,   mtrlrecv);
								_mtral.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk+1, mtrurecv);

								time3 = clock ();

								timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								_mtral.SplitMatrix (nlistloc, listloc,
													mtrlrecv, 
													gmtrl3.mtrarr);
								_mtral.SplitMatrix (nlistloc, listloc,
													mtrurecv, 
													gmtru3.mtrarr);

								mtrlrecv = mtrdummy;
								mtrurecv = mtrdummy;

// Update and free

								for (ilist=0;ilist<nlistloc;ilist++) {

									iblk = listloc[ilist];

									gmtrl2.ReadBlock (iblk);
									gmtru2.ReadBlock (iblk);

									time2 = clock ();

									mtrl2temp = fct.AddBlocks (_mtral, gmtrl2.mtrarr[iblk], gmtrl3.mtrarr[iblk]);
									mtru2temp = fct.AddBlocks (_mtral, gmtru2.mtrarr[iblk], gmtru3.mtrarr[iblk]);

									time3 = clock ();

									timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

									gmtrl2.mtrarr[iblk] = mtrl2temp;
									gmtru2.mtrarr[iblk] = mtru2temp;

									gmtrl3.mtrarr[iblk] = mtrdummy;
									gmtru3.mtrarr[iblk] = mtrdummy;

									mtrl2temp = mtrdummy;
									mtru2temp = mtrdummy;

									mtrlrecv = mtrdummy;
									mtrurecv = mtrdummy;

									gmtrl2.WriteBlock (iblk);
									gmtrl2.FreeBlock  (iblk);

									gmtru2.WriteBlock (iblk);
									gmtru2.FreeBlock  (iblk);

								};

							};

						} else {

// Receive and update

							for (ilist=0;ilist<nlistloc;ilist++) {

								iblk = listloc[ilist];

								time2 = clock ();

								_mtral.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk,   mtrlrecv);
								_mtral.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk+1, mtrurecv);

								time3 = clock ();

								timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								gmtrl2.ReadBlock (iblk);
								gmtru2.ReadBlock (iblk);

								time2 = clock ();

								mtrl2temp = fct.AddBlocks (_mtral, gmtrl2.mtrarr[iblk], mtrlrecv);
								mtru2temp = fct.AddBlocks (_mtral, gmtru2.mtrarr[iblk], mtrurecv);

								time3 = clock ();

								timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								gmtrl2.mtrarr[iblk] = mtrl2temp;
								gmtru2.mtrarr[iblk] = mtru2temp;

								mtrl2temp = mtrdummy;
								mtru2temp = mtrdummy;

								mtrlrecv = mtrdummy;
								mtrurecv = mtrdummy;

								gmtrl2.WriteBlock (iblk);
								gmtrl2.FreeBlock  (iblk);

								gmtru2.WriteBlock (iblk);
								gmtru2.FreeBlock  (iblk);

							};

						};

					};
				};

			};

		} else {

			int nlistloc = mvm.nlistuprcvarr[ilev];
			int *listloc = mvm.listuprcvarr[ilev];

//			cout << " Recv matrices from " << fathercpuid2 << endl;
//			OutArr(cout," List of recv ",nlistloc,listloc);

			if (_param.ooctyp == 0) {

				if (nlistloc > 0) {

					time2 = clock ();

					iblk = listloc[0];

					_mtral.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk,   mtrlrecv);
					_mtral.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk+1, mtrurecv);

					time3 = clock ();

					timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					_mtral.SplitMatrix (nlistloc, listloc,
										mtrlrecv, 
										gmtrl2.mtrarr);
					_mtral.SplitMatrix (nlistloc, listloc,
										mtrurecv, 
										gmtru2.mtrarr);

					mtrlrecv = mtrdummy;
					mtrurecv = mtrdummy;

				};

			} else {

				for (ilist=0;ilist<nlistloc;ilist++) {

					iblk = listloc[ilist];

					time2 = clock ();

					_mtral.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk,   mtrlrecv);
					_mtral.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk+1, mtrurecv);

					time3 = clock ();

					timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					gmtrl2.mtrarr[iblk] = mtrlrecv;
					gmtru2.mtrarr[iblk] = mtrurecv;

					mtrlrecv = mtrdummy;
					mtrurecv = mtrdummy;

					gmtrl2.WriteBlock (iblk);
					gmtrl2.FreeBlock  (iblk);

					gmtru2.WriteBlock (iblk);
					gmtru2.FreeBlock  (iblk);

				};

			};

		};

// Load diagonal data 

//		cout << " Current node:" << _tree.nodes[inodecurr] << endl;
//		_fout << " Current node:" << _tree.nodes[inodecurr] << endl;

		time2 = clock ();

		fctdiag.Ilu2LoadDiagonalData (inodecurr, _tree, mvm, iabext, jabext,
										_mtral, fct);

		time3 = clock ();

		timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Perform computations with blocks of the current node taking into account Schur data

		if (inodetype != exchangenode) {

			int iblkendlocal = _tree.nodes[inodecurr].indend;
			int isupbrd = _mtral.blksc[iblkendlocal+1];

//			if (_param.fcttyp != 2 && _param.fcttyp != 3) isupbrd = _mtral.nsupr;
			if (_param.fcttyp != 2) isupbrd = _mtral.nsupr;

			int nfrees_per_node = 10;
			int istepfree = (_tree.nodes[inodecurr].indend-_tree.nodes[inodecurr].indbeg+1) / nfrees_per_node;

			if (istepfree <= 0) istepfree = 1;

			int icount = -1;

			for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {

				icount++;

				if (_param.msglev > 1) {
					cout << " Myid = " << _tree.myid << " Iblk = " << iblk << " Isupbrd = " << isupbrd << endl;
				};
//				_fout << " Myid = " << _tree.myid << " Iblk = " << iblk << endl;

// Load diagonal data into fct structure

				int nlstblk = iabext[iblk+1]-iabext[iblk];
				int ibs = iabext[iblk];

				time2 = clock ();

				fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtral, fct);

				time3 = clock ();

				timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Init current block row of the matrix

				_mtral.ReadBlock (iblk);
				_mtrau.ReadBlock (iblk);

				time2 = clock ();

				fct.Ilu2InitBlkRow (iblk, _param.dshift,
									_mtral, _mtral.mtrarr[iblk], _mtrau.mtrarr[iblk], 
									mtrltemp, mtrutemp);

				time3 = clock ();

				timeini += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				_mtral.FreeBlock (iblk);
				_mtrau.FreeBlock (iblk);

// Update current block row by the Schur data

				if (inodecurr != inodesave) {

					gmtrl2.ReadBlock (iblk);
					gmtru2.ReadBlock (iblk);

					time2 = clock ();

					mtrl2temp = fct.AddBlocks (_mtral, mtrltemp, gmtrl2.mtrarr[iblk]);
					mtru2temp = fct.AddBlocks (_mtral, mtrutemp, gmtru2.mtrarr[iblk]);

					time3 = clock ();

					timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					mtrltemp = mtrl2temp;
					mtrutemp = mtru2temp;

					mtrl2temp = mtrdummy;
					mtru2temp = mtrdummy;

					gmtrl2.mtrarr[iblk] = mtrdummy;
					gmtru2.mtrarr[iblk] = mtrdummy;

				};

// Update current block row by the previous ones

				if (_tree.nodes[inodecurr].nchilds == 1 && inodetype != compnode2) {

					for (int j=iabextt[iblk];j<iabextt[iblk+1]-1;j++) {

						int jblk = jabextt[j];

//						cout << " update Jblk = " << jblk << endl;
//						_fout << " update Jblk = " << jblk << endl;

						gmtrl2.ReadBlock (jblk);
						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						if (icase == 1) {

							fct.Ilu2UpdateBlkRow (iblk, _param,
													_mtral, mtrltemp, mtrutemp,
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						} else {

							fct.Ilu2PointUpdateBlkRow (iblk, _param,
													_mtral, mtrltemp, mtrutemp,
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						};

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtrl2.FreeBlock (jblk);
						gmtru2.FreeBlock (jblk);

						mtrltemp = mtrl2temp;
						mtrutemp = mtru2temp;

						mtrl2temp = mtrdummy;
						mtru2temp = mtrdummy;

					};

				} else {

					for (int jblk=_tree.nodes[inodecurr].indbeg;jblk<iblk;jblk++) {

						gmtrl2.ReadBlock (jblk);
						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						if (icase == 1) {

							fct.Ilu2UpdateBlkRow (iblk, _param,
													_mtral, mtrltemp, mtrutemp,
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						} else {

							fct.Ilu2PointUpdateBlkRow (iblk, _param,
													_mtral, mtrltemp, mtrutemp,
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						};

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtrl2.FreeBlock (jblk);
						gmtru2.FreeBlock (jblk);

						mtrltemp = mtrl2temp;
						mtrutemp = mtru2temp;

						mtrl2temp = mtrdummy;
						mtru2temp = mtrdummy;

					};

				};

// Print Schur sparsity if necessary

				if (false && inodecurr != inodesave) {

					int isupbeg = _mtral.blksc[iblk];
					int isupend = _mtral.blksc[iblk+1]-1;

					CSMatrix aschur;

					aschur = mtrltemp.CSMatrix::FilterMatrix(isupbeg,isupend);

					char strbuff[256];

					sprintf(strbuff, "%s%s%d%s%d%s",_param.path,"DSchur_",_tree.myid,"_",iblk,".ps");

					aschur.A2Ps (3,strbuff,0,&isupbeg);

					int nproc = 2;
					int nchilds = 1;

					CTree tree (nproc,nchilds);

					CSMatrix aschursymm;

					aschursymm = aschur.SymmMtr ();

					int nblks;
					int *blks, *order;

					int n = aschur.GetN ();

					order = new int [n];

					aschursymm.PartOrdMtr (tree, order);

					tree.Partition (nblks,blks);

					int nsupmax=1000;

					tree.TreeInBlocks (nsupmax,nblks,blks);

					CSMatrix ao;

					ao = aschursymm.OrdMtr (order);

					sprintf(strbuff, "%s%s%d%s%d%s",_param.path,"DSchurOrd_",_tree.myid,"_",iblk,".ps");

					ao.A2Ps (3,strbuff,nblks,blks);

				};

// Factorize current block row

				time2 = clock ();

				if (icase == 1) {

					fct.Ilu2FctBlkRow (iblk, _param, isupbrd,
										_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
										_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
										gmtrl2.mtrarr[iblk], gmtru2.mtrarr[iblk]);

				} else {

					fct.Ilu2PointFctBlkRow (iblk, _param, isupbrd,
										_mtral, mtrltemp, mtrutemp, _mtral.mtrarr[iblk],
										_mtrl.mtrarr[iblk], _mtru.mtrarr[iblk], 
										gmtrl2.mtrarr[iblk], gmtru2.mtrarr[iblk]);

				};

				time3 = clock ();

				timefct += (double) (time3-time2) / (double) CLOCKS_PER_SEC;
				if (nodeend == inodecurr) {
					timefctd += (double) (time3-time2) / (double) CLOCKS_PER_SEC;
				};

				_mtrl.WriteBlock (iblk);
				_mtrl.FreeBlock  (iblk);

				_mtru.WriteBlock (iblk);
				_mtru.FreeBlock  (iblk);

				gmtrl2.WriteBlock (iblk);
				gmtrl2.FreeBlock  (iblk);

				gmtru2.WriteBlock (iblk);
				gmtru2.FreeBlock  (iblk);

				mtrltemp = mtrdummy;
				mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

				nlstblk = iabext[iblk+1]-iabext[iblk];
				ibs = iabext[iblk];

				time2 = clock ();

				fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtral, fct);

				time3 = clock ();

				timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Free data if necessary

				if (icount%istepfree == 0) {

					int isupend = _mtral.blksc[iblk+1]-1;

					for (int jblk=_tree.nodes[inodecurr].indbeg;jblk<iblk;jblk++) {

						gmtrl2.ReadBlock    (jblk);
						gmtrl2.FilterBlock  (jblk,isupend);
						gmtrl2.ReWriteBlock (jblk);
						gmtrl2.FreeBlock    (jblk);

						gmtru2.ReadBlock    (jblk);
						gmtru2.FilterBlock  (jblk,isupend);
						gmtru2.ReWriteBlock (jblk);
						gmtru2.FreeBlock    (jblk);

					};
				};

			};

// Scan second order matrices and filter unnecessary data

			int iblkbeg = _tree.nodes[inodecurr].indbeg;
			int iblkend = _tree.nodes[inodecurr].indend;
			int isupend = _mtral.blksc[iblkend+1]-1;

			for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {

				gmtrl2.ReadBlock    (iblk);
				gmtrl2.FilterBlock  (iblk,isupend);
				gmtrl2.ReWriteBlock (iblk);
				gmtrl2.FreeBlock    (iblk);

				gmtru2.ReadBlock    (iblk);
				gmtru2.FilterBlock  (iblk,isupend);
				gmtru2.ReWriteBlock (iblk);
				gmtru2.FreeBlock    (iblk);

			};

// Compute the block sparsity structure for the whole bordering

			int *iabbord, *jabbord;

			pmtral = &_mtral;

			int iblk;
			for (iblk=0;iblk<_mtral.nblksr;iblk++) {
				pmtrarr[iblk] = &gmtrl2.mtrarr[iblk];
			};

			pmtral->MtrBlockSparsity (pmtrarr, iblkbeg, iblkend, // Compute block sparsity of the matrix
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

// Compute the list of Schur blocks to be updated

			ilev = _tree.nodes[inodecurr].nodelv;

			nlistschur = mvm.nlistfctarr[ilev];
			listschur = mvm.listfctarr[ilev];

// Update Schur blocks according to the transposed block sparsity

			for (ilist=0;ilist<nlistschur;ilist++) {

				iblk = listschur[ilist];

				if (iabbordtr[iblk+1] > iabbordtr[iblk]) {

// Load diagonal data

					int nlstblk = iabext[iblk+1]-iabext[iblk];
					int ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
													_mtral, fct);

					time3 = clock ();

					timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Read data if necessary

					gmtrl2.ReadBlock (iblk);
					gmtru2.ReadBlock (iblk);

// Update

					for (int j=iabbordtr[iblk];j<iabbordtr[iblk+1];j++) {

						int jblk = jabbordtr[j];

						gmtrl2.ReadBlock (jblk);
						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						if (icase == 1) {

							fct.Ilu2UpdateBlkRow (iblk, _param,
													_mtral, gmtrl2.mtrarr[iblk], gmtru2.mtrarr[iblk],
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						} else {

							fct.Ilu2PointUpdateBlkRow (iblk, _param,
													_mtral, gmtrl2.mtrarr[iblk], gmtru2.mtrarr[iblk],
													gmtrl2.mtrarr[jblk], gmtru2.mtrarr[jblk],
													mtrl2temp, mtru2temp);

						};

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtrl2.FreeBlock (jblk);
						gmtru2.FreeBlock (jblk);

						gmtrl2.mtrarr[iblk] = mtrl2temp;
						gmtru2.mtrarr[iblk] = mtru2temp;

						mtrl2temp = mtrdummy;
						mtru2temp = mtrdummy;

					};

// Store data to the disk if necessary

					gmtrl2.WriteBlock (iblk);
					gmtrl2.FreeBlock  (iblk);

					gmtru2.WriteBlock (iblk);
					gmtru2.FreeBlock  (iblk);

// Unload diagonal data

					nlstblk = iabext[iblk+1]-iabext[iblk];
					ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtral, fct);

					time3 = clock ();

					timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				};

			};

// Free work arrays

			delete [] iabbord;
			delete [] jabbord;
			delete [] iabbordtr;
			delete [] jabbordtr;

// Free second order matrices of the current node

			for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {
				gmtrl2.mtrarr[iblk] = mtrdummy;
				gmtru2.mtrarr[iblk] = mtrdummy;
			};

		} else {

// Perform filtering computations in the exchange node

			int iblkbegfather = _tree.nodes[fathernode].indbeg;

			for (ichild=0;ichild<_tree.nodes[inodecurr].nchilds2;ichild++) {

				childnodecurr = _tree.nodes[inodecurr].childs2[ichild];

				int iblkbegloc = _tree.nodes[childnodecurr].indbeg;
				int iblkendloc = _tree.nodes[childnodecurr].indend;

				int isupbeg = _mtral.blksr[iblkbegloc];
				int isupend = _mtral.blksr[iblkendloc+1]-1;
				int isuprest = _mtral.blksr[iblkbegfather];

				double theta2 = _param.theta / 2.0e0;

				for (iblk=iblkbegloc;iblk<=iblkendloc;iblk++) {

// Load diagonal data into fct structure

					int nlstblk = iabext[iblk+1]-iabext[iblk];
					int ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
													_mtral, fct);

					time3 = clock ();

					timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Filter block data

					gmtrl2.ReadBlock    (iblk);
					gmtru2.ReadBlock    (iblk);

					int nzaini = gmtrl2.mtrarr[iblk].nzatot;

					time2 = clock ();

					gmtrl2.FilterSchurBlock (iblk, isupbeg, isupend, isuprest, theta2, fct);
					gmtru2.FilterSchurBlock (iblk, isupbeg, isupend, isuprest, theta2, fct);

					time3 = clock ();

					timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					int nzafin = gmtrl2.mtrarr[iblk].nzatot;
					if (_param.msglev > 1) {
						cout << " Filtering Iblk = " << iblk << " nzaini = " << nzaini << " nzafin = " << nzafin << endl;
					};

					nzstot += nzaini-nzafin;

					gmtrl2.ReWriteBlock (iblk);
					gmtrl2.FreeBlock    (iblk);

					gmtru2.ReWriteBlock (iblk);
					gmtru2.FreeBlock    (iblk);

// Unload diagonal data from fct structure

					nlstblk = iabext[iblk+1]-iabext[iblk];
					ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
													_mtral, fct);

					time3 = clock ();

					timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				};

			};

		};

// Exchange the data if necessary and perform updates

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = mvm.nlistupsndarr[ilev];
				int *listloc = mvm.listupsndarr[ilev];

//				cout << " Send matrices to " << fathercpuid << endl;
//				OutArr(cout," List of send ",nlistloc,listloc);

				if (_param.ooctyp == 0) {

					if (nlistloc > 0) {

						_mtral.CombineMatrices (nlistloc, listloc, 
												gmtrl2.mtrarr, mtrlsend);
						_mtral.CombineMatrices (nlistloc, listloc, 
												gmtru2.mtrarr, mtrusend);

						iblk = listloc[0];

						time2 = clock ();

						_mtral.SendMatrix (_tree.comm, fathercpuid, 2*iblk,   mtrlsend);
						_mtral.SendMatrix (_tree.comm, fathercpuid, 2*iblk+1, mtrusend);

						time3 = clock ();

						timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						mtrlsend = mtrdummy;
						mtrusend = mtrdummy;

						for (ilist=0;ilist<nlistloc;ilist++) {

							iblk = listloc[ilist];

							gmtrl2.mtrarr[iblk] = mtrdummy;
							gmtru2.mtrarr[iblk] = mtrdummy;

						};

					};

				} else {

					for (ilist=0;ilist<nlistloc;ilist++) {

						iblk = listloc[ilist];

						gmtrl2.ReadBlock (iblk);
						gmtru2.ReadBlock (iblk);

						time2 = clock ();

						_mtral.SendMatrix (_tree.comm, fathercpuid, 2*iblk,   gmtrl2.mtrarr[iblk]);
						_mtral.SendMatrix (_tree.comm, fathercpuid, 2*iblk+1, gmtru2.mtrarr[iblk]);

						time3 = clock ();

						timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtrl2.mtrarr[iblk] = mtrdummy;
						gmtru2.mtrarr[iblk] = mtrdummy;

					};

				};

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inodecurr].nchilds2;ichild++) {

				int childnode = _tree.nodes[inodecurr].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = mvm.nlistupsndarr[ilev];
				int *listloc = mvm.listupsndarr[ilev];

//				cout << " Send matrices to " << childcpuid << endl;
				if (childcpuid != _tree.myid) {

					if (_param.ooctyp == 0) {

// Compute sublist

						int nlistsub = 0;

						for (ilist=0;ilist<nlistloc;ilist++) {

							iblk = listloc[ilist];

							if (mvm.blk2cpu[iblk] == childcpuid) {
								listsub[nlistsub] = iblk;
								nlistsub++;
							};

						};

						if (nlistsub > 0) {

							_mtral.CombineMatrices (nlistsub, listsub, 
													gmtrl2.mtrarr, mtrlsend);
							_mtral.CombineMatrices (nlistsub, listsub, 
													gmtru2.mtrarr, mtrusend);

							iblk = listsub[0];

							time2 = clock ();

							_mtral.SendMatrix (_tree.comm, childcpuid, 2*iblk,   mtrlsend);
							_mtral.SendMatrix (_tree.comm, childcpuid, 2*iblk+1, mtrusend);

							time3 = clock ();

							timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

							mtrlsend = mtrdummy;
							mtrusend = mtrdummy;

							for (ilist=0;ilist<nlistsub;ilist++) {

								iblk = listsub[ilist];

								gmtrl2.mtrarr[iblk] = mtrdummy;
								gmtru2.mtrarr[iblk] = mtrdummy;

							};

						};

					} else {

						for (ilist=0;ilist<nlistloc;ilist++) {

							iblk = listloc[ilist];

							if (mvm.blk2cpu[iblk] == childcpuid) {

//								cout << " send iblk = " << iblk << endl;
								gmtrl2.ReadBlock (iblk);
								gmtru2.ReadBlock (iblk);

								time2 = clock ();

								_mtral.SendMatrix (_tree.comm, childcpuid, 2*iblk,   gmtrl2.mtrarr[iblk]);
								_mtral.SendMatrix (_tree.comm, childcpuid, 2*iblk+1, gmtru2.mtrarr[iblk]);

								time3 = clock ();

								timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								gmtrl2.mtrarr[iblk] = mtrdummy;
								gmtru2.mtrarr[iblk] = mtrdummy;

							};

						};

					};

				};

			};

		};

// Unload diagonal data from fct structure

		time2 = clock ();

		fctdiag.Ilu2UnloadDiagonalData (inodecurr, _tree, mvm, iabext, jabext,
										_mtral, fct);

		time3 = clock ();

		timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

NextNode:;

	};

// Store point/block scaling data in L and U structures

	if (_mtrl.ivect < 2) {
		fctdiag.StorePointScaling (_mtrl);
		fctdiag.StorePointScaling (_mtru);
	} else {
		fctdiag.StoreBlockScaling ('L', _mtrl);
		fctdiag.StoreBlockScaling ('U', _mtru);
	};

// Close files

	if (_param.ooctyp > 0) {

		gmtrl2.CloseFiles();
		gmtru2.CloseFiles();

		fctdiag.CloseFiles();

	};

// Finalize time measurement

	if (_tree.nproc != 1) {
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

// Output timing results

	if (_tree.myid == _tree.rootcpu) {
		if (_param.msglev > 1) {
			cout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
			cout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s.; TFctD = " << timefctd << endl;
			cout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
		};
		if (_param.msglev > 0) {
			_fout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
			_fout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s.; TFctD = " << timefctd << endl;
			_fout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
		};
	};

// Output Ilu2 factorization statistics

	int nzatotal = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzatotal += _mtral.mtrarr[iblk].nzatot;
	};

	int nzmema = 2 * nzatotal;

	for (ilist=0;ilist<_mtral.nlist;ilist++) {
		iblk = _mtral.listb[ilist];
		for (isup=_mtral.blksr[iblk];isup<_mtral.blksr[iblk+1];isup++) {
			blai = _mtral.sprndr[isup+1]-_mtral.sprndr[isup];
			nzmema -= blai*blai;
		};
	};

	int nzutot = 0;
	for (iblk=0;iblk<_mtral.nblksr;iblk++) {
		nzutot += _mtrl.mtrarr[iblk].nzatot;
	};

	int nzrtot = fct.nzr;

// Compute global data

	double arr0[7], arr1[7];

	arr0[0] = (double) nzutot;
	arr0[1] = (double) nzrtot;
	arr0[2] = (double) nzstot;
	arr0[3] = (double) nzatotal;
	arr0[4] = (double) nzmema;
	arr0[5] = fct.ops;
	arr0[6] = (double) fct.nmodsv;

	if (_tree.nproc != 1) {
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
	double dnza = arr1[3];
	double dnzamem = arr1[4];
	ops = arr1[5];
	int nmodtot = (int) arr1[6];

// Output resulting statistics

	densu = 1.0e2 * dnzu / dnza;
	densr = 1.0e2 * dnzr / dnza;
	denss = 1.0e2 * dnzs / dnza;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / dnzamem;

	if (_tree.myid == _tree.rootcpu) {

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

};

// Author: Kharchenko S.A.
// CFctDiagR: Exchange the scaling
//========================================================================================
void CFctDiagR::ExchangeScaling (bool _is2index, const CTree &_tree, const CMvm &_mvm, const CGSMatrix &_gmtra) { // Exchange the scaling

	const char *funcname = "ExchangeScaling";

//	double dzero = 0.0e0;

// Allocate working array

	int nlistarr=0;
	int *listarr;

	listarr = new int [_gmtra.nblksr];
	if (!listarr) MemoryFail (funcname);

// Search the tree

	int nlistloc;
	int *listloc;

	int inode, inodetype, ichild;
	int iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	int ilev;
	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode;

		inodetype = _tree.nodes[inode].nodetype;

		iblkbeg = _tree.nodes[inode].indbeg;
		iblkend = _tree.nodes[inode].indend;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Take the list

		nlistloc = _mvm.nlistdnrcvarr[ilev];
		listloc = _mvm.listdnrcvarr[ilev];

// Switch the cases

		if (inodetype != exchangenode) {

			if (fathercpuid != _tree.myid) {

// Allocate and register the memory for scaling

				AllocateScaling (_is2index, nlistloc, listloc, _gmtra);

// Receive scaling data if necessary

//				cout << " Receive scaling from " << fathercpuid << endl;
//				OutArr(cout," List scaling recv = ",nlistloc,listloc);
				ReceiveScaling (_is2index, _tree.comm, nlistloc, listloc, fathercpuid, inode, _gmtra);

// Store the data to the disk if necessary

				if (nfiles > 0) {
					WriteScaling (nlistloc, listloc);
				};

			};

		} else {

			for (ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

// Rebuild the list

					nlistarr = 0;

					for (int ilist=0;ilist<nlistloc;ilist++) {
						int iblk = listloc[ilist];
						if (_mvm.blk2cpu[iblk] == childcpuid) {
							listarr[nlistarr] = iblk;
							nlistarr++;
						};
					};

// Allocate and register the memory for scaling

					AllocateScaling (_is2index, nlistarr, listarr, _gmtra);

// Receive scaling data if necessary

//					cout << " Receive scaling from " << childcpuid << endl;
//					OutArr(cout," List scaling recv = ",nlistarr,listarr);
					ReceiveScaling (_is2index, _tree.comm, nlistarr, listarr, childcpuid, inode, _gmtra);

// Store the data to the disk if necessary

					if (nfiles > 0) {
						WriteScaling (nlistarr, listarr);
					};

				};
			};
		};

// Send the scaling data to the childs if necessary

		nlistloc = _mvm.nlistdnsndarr[ilev];
		listloc = _mvm.listdnsndarr[ilev];

// Switch the cases

		if (inodetype != compnode2) {

			for (ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

// Rebuild the list

					nlistarr = 0;

					for (int ilist=0;ilist<nlistloc;ilist++) {
						int iblk = listloc[ilist];
						if (_mvm.blk2cpu[iblk] != childcpuid) {
							listarr[nlistarr] = iblk;
							nlistarr++;
						};
					};

// Read the data from disk if necessary

					if (nfiles > 0) {
						ReadScaling (nlistarr, listarr);
					};

// Send

//					cout << " Send scaling to " << childcpuid << endl;
//					OutArr(cout," List scaling send = ",nlistarr,listarr);
					SendScaling (_is2index, _tree.comm, nlistarr, listarr, childcpuid, childnode, _gmtra);

// Free the scaling data from memory

					if (nfiles > 0) {
						FreeScaling (nlistarr, listarr);
					};

				};

			};

		} else {

// Read the data from disk if necessary

			if (nfiles > 0) {
				ReadScaling (nlistloc, listloc);
			};

// Send

//			cout << " Send scaling to " << fathercpuid2 << endl;
//			OutArr(cout," List scaling send = ",nlistloc,listloc);
			SendScaling (_is2index, _tree.comm, nlistloc, listloc, fathercpuid2, fathernode2, _gmtra);

// Free the scaling data from memory

			if (nfiles > 0) {
				FreeScaling (nlistloc, listloc);
			};

		};

NextNode:;

	};

// Free work array

	delete [] listarr;

};

// Author: Kharchenko S.A.
// CFctDiagR: Exchange the scaling in the general case
//========================================================================================
void CFctDiagR::ExchangeScalingGeneral (bool _is2index, const CTree &_tree, int *_blk2cpu, int _nlistschur, int *_listschur,
														int *_iab, int *_jab, const CGSMatrix &_gmtra) { // Exchange the scaling in the general case

	const char *funcname = "ExchangeScalingGeneral";

	int nproc = _tree.GetNproc ();
	int myid = _tree.GetMyid ();

// Allocate arrays for scaling for Schur data

//	AllocateScaling (_is2index, _nlistschur, _listschur, _gmtra);

// Create the local list of scaling data to be received

	int nblksloc = _gmtra.GetNblksr ();

	int nlistarr=0;
	int *listarr, *imask;

	listarr = new int [nblksloc];
	if (!listarr) MemoryFail (funcname);
	imask = new int [nblksloc];
	if (!imask) MemoryFail (funcname);

	int icycle = -1;

	int i, j, jj;

	for (i=0;i<nblksloc;i++) imask[i] = icycle;

	icycle++;

	for (i=0;i<nblksloc;i++) {
		if (_blk2cpu[i] == myid) {
			for (j=_iab[i];j<_iab[i+1];j++) {
				jj = _jab[j];
				if (imask[jj] != icycle && _blk2cpu[jj] != myid) {
					listarr[nlistarr] = jj;
					nlistarr++;
					imask[jj] = icycle;
				};
			};
		};
	};

	for (i=0;i<_nlistschur;i++) {
		jj = _listschur[i];
		if (imask[jj] != icycle && _blk2cpu[jj] != myid) {
			listarr[nlistarr] = jj;
			nlistarr++;
			imask[jj] = icycle;
		};
	};

	qsort (listarr, nlistarr, sizeof(int), compint);

// Allocate arrays for scaling

	AllocateScaling (_is2index, nlistarr, listarr, _gmtra);

// Split the total list over the processors

	int *iacpu;
	int *iptr;
	int *jacpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);
	jacpu = new int [nlistarr];
	if (!jacpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iblk, iproc;

	for (i=0;i<nlistarr;i++) {
		iblk = listarr[i];
		iproc = _blk2cpu[iblk];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int k;

	for (i=0;i<nlistarr;i++) {
		iblk = listarr[i];
		iproc = _blk2cpu[iblk];
		k = iptr[iproc];
		jacpu[k] = iblk;
		iptr[iproc]++;
	};

// Create local send arrays which contain lists of hypercells to be obtained

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nproc;

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

	int* pIntLoc;

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) (jacpu+iacpu[i]);
	};

// Exchange the data

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info;

	info = CMPIExchange::DataExchangeMPI (pcomm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if (info) throw " CFctDiagR::ExchangeScalingGeneral: Error in DataExchangeMPI";

// Free send data

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare sends of matrix data (in one send structure)

	int NObjSend2;
	int* ObjTypeSend2;
	int* ObjIDSend2;
	int* CpuIDSend2;
	int* ObjSizeSend2;
	char** ObjSend2;

	NObjSend2 = NObjRecv;

	ObjTypeSend2 = new int [NObjSend2];
	if(!ObjTypeSend2) MemoryFail(funcname);
	ObjIDSend2 = new int [NObjSend2];
	if(!ObjIDSend2) MemoryFail(funcname);
	CpuIDSend2 = new int [NObjSend2];
	if(!CpuIDSend2) MemoryFail(funcname);
	ObjSizeSend2 = new int [NObjSend2];
	if(!ObjSizeSend2) MemoryFail(funcname);
	ObjSend2 = new char* [NObjSend2];
	if(!ObjSend2) MemoryFail(funcname);

	int nLoc;

	for(i = 0; i < NObjSend2; i++)
	{
		ObjTypeSend2[i] = 1;
		ObjIDSend2[i] = i;
		CpuIDSend2[i] = CpuIDRecv[i];
		pIntLoc = (int *) ObjRecv[i];
		nLoc = (int) (ObjSizeRecv[i] / sizeof(int));
		CreateSetOfScalings (_is2index, nLoc, pIntLoc,
									ObjSend2[i], ObjSizeSend2[i], _gmtra);
	};

// Exchange the data

	int NObjRecv2;
	int* ObjTypeRecv2;
	int* ObjIDRecv2;
	int* CpuIDRecv2;
	int* ObjSizeRecv2;
	char** ObjRecv2;

	info = CMPIExchange::DataExchangeMPI (pcomm,
														NObjSend2, ObjTypeSend2, ObjIDSend2, CpuIDSend2,
														ObjSizeSend2, ObjSend2,
														NObjRecv2, ObjTypeRecv2, ObjIDRecv2, CpuIDRecv2,
														ObjSizeRecv2, ObjRecv2);
	if (info) throw " CFctDiagR::ExchangeScalingGeneral: Error in DataExchangeMPI";

// Free send data

	delete [] ObjTypeSend2;
	delete [] ObjIDSend2;
	delete [] CpuIDSend2;
	delete [] ObjSizeSend2;
	for(i = 0; i < NObjSend2; i++)
	{
		delete [] ObjSend2[i];
	};
	delete [] ObjSend2;

// Store received data

	for(i = 0; i < NObjRecv2; i++)
	{
		iproc = CpuIDRecv2[i];
		RestoreSetOfScalings (_is2index, iacpu[iproc+1]-iacpu[iproc], jacpu+iacpu[iproc],
										ObjRecv2[i], ObjSizeRecv2[i], _gmtra);
	};

// Free received data (second and first)

	delete [] ObjTypeRecv2;
	delete [] ObjIDRecv2;
	delete [] CpuIDRecv2;
	delete [] ObjSizeRecv2;
	for(i = 0; i < NObjRecv2; i++)
	{
		delete [] ObjRecv2[i];
	};
	delete [] ObjRecv2;

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++)
	{
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] listarr;
	delete [] imask;
	delete [] iacpu;
	delete [] iptr;
	delete [] jacpu;

};

// Author: Kharchenko S.A.
// CFctDiagR: Create scaling data for set of rows
//========================================================================================
void CFctDiagR::CreateSetOfScalings (bool _is2index, int _nlistloc, int *_listloc, // Create scaling data for set of rows
													char *&_pdata, int &_length, const CGSMatrix &_gmtra) {

	const char *funcname = "CreateSetOfScalings";

	double *sendarr;

	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		int nzloc = diasize[iblk];

		nztot += nlocblk+4*nzloc;

	};

	sendarr = new double [nztot];
	if (!sendarr) MemoryFail (funcname);

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc  = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sendarr[nztot+kii] = sclptloc[kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scllloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scluloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvlloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvuloc[kii];
		nztot += nzloc;

	};

	_pdata = (char *) sendarr;
	_length = nztot * sizeof(double);

};

// Author: Kharchenko S.A.
// CFctDiagR: Restore received scaling data
//========================================================================================
void CFctDiagR::RestoreSetOfScalings (bool _is2index, int _nlistloc, int *_listloc, // Restore received scaling data
													char *_pdata, int _length,
													const CGSMatrix &_gmtra) {

//	const char *funcname = "RestoreSetOfScalings";

//	double dzero = 0.0e0;

	double *recvarr;
	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		nztot += nlocblk+4*nzloc;

	};

	recvarr = (double *)_pdata;

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = recvarr[nztot+kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) scllloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) scluloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvlloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvuloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

	};

	if (nztot * sizeof(double) != _length) throw " CFctDiagR::RestoreSetOfScalings: wrong size of received block ";

};

// Author: Kharchenko S.A.
// CFctDiagR: Allocate the memory for scaling according to the prescribed list
//========================================================================================
void CFctDiagR::AllocateScaling (bool _is2index, int _nlistloc, int *_listloc, // Allocate the memory for scaling according to the prescribed list
									const CGSMatrix &_gmtra) {

	const char *funcname = "AllocateScaling";

	double dzero = 0.0e0;

	for (int ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int isupbeg = _gmtra.blksr[iblk];
		int isupend = _gmtra.blksr[iblk+1]-1;

		int nz=0;
		int nlocblk=0;

		if (_is2index) {
			nz = isupend-isupbeg+1;
			nlocblk = isupend-isupbeg+1;
		} else {
			for (int isup=isupbeg;isup<=isupend;isup++) {
				int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
				nz += blai*blai;
				nlocblk += blai;
			};
		};

		double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc;
		double *dpivloc;

		sclptloc = new double [nlocblk];
		if (!sclptloc) MemoryFail (funcname);
		scllloc = new double [nz];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new double [nz];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new double [nz];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new double [nz];
		if (!sclinvuloc) MemoryFail (funcname);
		diaggloc = new double [nz];
		if (!diaggloc) MemoryFail (funcname);
		dpivloc = new double [nlocblk];
		if (!dpivloc) MemoryFail (funcname);

		for (int kii=0;kii<nz;kii++) diaggloc[kii] = dzero;

		delete [] sclpt  [iblk];
		delete [] scll   [iblk];
		delete [] sclu   [iblk];
		delete [] sclinvl[iblk];
		delete [] sclinvu[iblk];
		delete [] diagg  [iblk];
		delete [] dpiv   [iblk];

		diasize[iblk] = nz;
		sclpt  [iblk] = sclptloc;
		scll   [iblk] = scllloc;
		sclu   [iblk] = scluloc;
		sclinvl[iblk] = sclinvlloc;
		sclinvu[iblk] = sclinvuloc;
		diagg  [iblk] = diaggloc;
		dpiv   [iblk] = dpivloc;

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Free the scaling
//========================================================================================
void CFctDiagR::FreeScaling (int _nlistloc, int *_listloc) { // Free the scaling

	const char *funcname = "FreeScaling";

	double *sclptloc;
	double *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	for (int ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		delete [] sclptloc;
		delete [] scllloc;
		delete [] scluloc;
		delete [] sclinvlloc;
		delete [] sclinvuloc;

		sclptloc = new double [0];
		if (!sclptloc) MemoryFail (funcname);
		scllloc = new double [0];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new double [0];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new double [0];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new double [0];
		if (!sclinvuloc) MemoryFail (funcname);

		sclpt  [iblk] = sclptloc;
		scll   [iblk] = scllloc;
		sclu   [iblk] = scluloc;
		sclinvl[iblk] = sclinvlloc;
		sclinvu[iblk] = sclinvuloc;

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Write scaling data to the disk
//========================================================================================
void CFctDiagR::WriteScaling (int _nlistloc, int *_listloc) { // Write scaling data to the disk

	const char *funcname = "WriteScaling";

	for (int ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		double *scllloc, *scluloc, *sclinvlloc, *sclinvuloc, *diaggloc;
		int nz;

		nz         = diasize[iblk];

		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];
		diaggloc   = diagg  [iblk];

		int ibsloc = ibsfile[ifile];

		blk2file[5*iblk] = ifile;
		blk2bs  [5*iblk] = ibsloc;

		FPut (files[ifile], nz, scllloc, ibsloc);

		ibsloc += nz;

		blk2file[5*iblk+1] = ifile;
		blk2bs  [5*iblk+1] = ibsloc;

		FPut (files[ifile], nz, scluloc, ibsloc);

		ibsloc += nz;

		blk2file[5*iblk+2] = ifile;
		blk2bs  [5*iblk+2] = ibsloc;

		FPut (files[ifile], nz, sclinvlloc, ibsloc);

		ibsloc += nz;

		blk2file[5*iblk+3] = ifile;
		blk2bs  [5*iblk+3] = ibsloc;

		FPut (files[ifile], nz, sclinvuloc, ibsloc);

		ibsloc += nz;

		blk2file[5*iblk+4] = ifile;
		blk2bs  [5*iblk+4] = ibsloc;

		FPut (files[ifile], nz, diaggloc, ibsloc);

		ibsloc += nz;

		ibsfile[ifile] = ibsloc;

		ifile++;

		if (ifile >= nfiles) ifile = 0;

		delete [] scllloc;
		delete [] scluloc;
		delete [] sclinvlloc;
		delete [] sclinvuloc;
		delete [] diaggloc;

		scllloc = new double [0];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new double [0];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new double [0];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new double [0];
		if (!sclinvuloc) MemoryFail (funcname);
		diaggloc = new double [0];
		if (!diaggloc) MemoryFail (funcname);

		scll   [iblk] = scllloc;
		sclu   [iblk] = scluloc;
		sclinvl[iblk] = sclinvlloc;
		sclinvu[iblk] = sclinvuloc;
		diagg  [iblk] = diaggloc;

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Read scaling data from the disk
//========================================================================================
void CFctDiagR::ReadScaling (int _nlistloc, int *_listloc) { // Read scaling data from the disk

	const char *funcname = "ReadScaling";

	for (int ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		double *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

		int nz  = diasize[iblk];

		delete [] scll   [iblk];
		delete [] sclu   [iblk];
		delete [] sclinvl[iblk];
		delete [] sclinvu[iblk];

		scllloc = new double [nz];
		if (!scllloc) MemoryFail (funcname);
		scluloc = new double [nz];
		if (!scluloc) MemoryFail (funcname);
		sclinvlloc = new double [nz];
		if (!sclinvlloc) MemoryFail (funcname);
		sclinvuloc = new double [nz];
		if (!sclinvuloc) MemoryFail (funcname);

		scll   [iblk] = scllloc;
		sclu   [iblk] = scluloc;
		sclinvl[iblk] = sclinvlloc;
		sclinvu[iblk] = sclinvuloc;

		int ifileloc  = blk2file[iblk*5];
		int ibsloc    = blk2bs  [iblk*5];

		FGet (files[ifileloc],nz,scllloc,ibsloc);

		ifileloc  = blk2file[iblk*5+1];
		ibsloc    = blk2bs  [iblk*5+1];

		FGet (files[ifileloc],nz,scluloc,ibsloc);

		ifileloc  = blk2file[iblk*5+2];
		ibsloc    = blk2bs  [iblk*5+2];

		FGet (files[ifileloc],nz,sclinvlloc,ibsloc);

		ifileloc  = blk2file[iblk*5+3];
		ibsloc    = blk2bs  [iblk*5+3];

		FGet (files[ifileloc],nz,sclinvuloc,ibsloc);

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Load diagonal data for the computations with the set of block rows in parallel mode
//========================================================================================
void CFctDiagR::Ilu2LoadDiagonalData (int _inode, // Load diagonal data for the computations with the set of block rows in parallel mode
										const CTree &_tree, const CMvm &_mvm, 
										const int *_iabext, const int *_jabext,
										const CGSMatrix &_gmtra, CFctR &_fct) {

//	const char *funcname = "Ilu2LoadDiagonalData";

// Take the list of blocks

	int ilev = _tree.nodes[_inode].nodelv;

	int nlistloc;
	int *listloc;

	nlistloc = _mvm.nlistuprcvarr[ilev];
	listloc  = _mvm.listuprcvarr[ilev];

// Take control data

	int inodetype = _tree.nodes[_inode].nodetype;

//	int fathernode = _tree.nodes[_inode].fatherid;
//	int fathercpuid = _tree.nodes[_inode].fathercpu;
//	int fathernode2 = _tree.nodes[_inode].fatherid2;
	int fathercpuid2 = _tree.nodes[_inode].fathercpu2;

// Switch the cases

	if (inodetype != compnode2) {

// Receive diagg arrays from child nodes if necessary and update local array

		int nchilds = _tree.nodes[_inode].nchilds;

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[_inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			if (childcpuid != _tree.myid) {

//				cout << " Recv diag data from " << childcpuid << endl;
//				OutArr(cout," List recv diag",nlistloc,listloc);
				ReceiveDiag (_tree.comm, '+', nlistloc, listloc, childcpuid, _inode);

			};

		};

	} else {

		if (fathercpuid2 != _tree.myid) {

//			cout << " Recv diag data from " << fathercpuid2 << endl;
//			OutArr(cout," List recv diag",nlistloc,listloc);
			ReceiveDiag (_tree.comm, '=', nlistloc, listloc, fathercpuid2, _inode);

		};

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Unload diagonal data for the computations with the set of block rows in parallel mode
//========================================================================================
void CFctDiagR::Ilu2UnloadDiagonalData (int _inode, // Unload diagonal data for the computations with the set of block rows in parallel mode
										const CTree &_tree, const CMvm &_mvm,
										const int *_iabext, const int *_jabext,
										const CGSMatrix &_gmtra, CFctR &_fct) {

	const char *funcname = "Ilu2UnloadDiagonalData";

// Take the list of blocks

	int ilev = _tree.nodes[_inode].nodelv;

	int nlistloc;
	int *listloc;

	nlistloc = _mvm.nlistupsndarr[ilev];
	listloc  = _mvm.listupsndarr [ilev];

// Take control data

	int inodetype = _tree.nodes[_inode].nodetype;

	int fathernode = _tree.nodes[_inode].fatherid;
	int fathercpuid = _tree.nodes[_inode].fathercpu;
//	int fathernode2 = _tree.nodes[_inode].fatherid2;
//	int fathercpuid2 = _tree.nodes[_inode].fathercpu2;

// Switch the cases

	if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

		if (fathercpuid != _tree.myid) {

//			cout << " Send diag data to " << fathercpuid << endl;
//			OutArr(cout," List send diag",nlistloc,listloc);
			SendDiag (_tree.comm, nlistloc, listloc, fathercpuid, fathernode, _gmtra);

		};

	} else {

		for (int ichild=0;ichild<_tree.nodes[_inode].nchilds2;ichild++) {

			int childnode = _tree.nodes[_inode].childs2[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			if (childcpuid != _tree.myid) {

// Rebuild the list

				int nlistarr;
				int *listarr;

				nlistarr = 0;

				listarr = new int [nlistloc];
				if (!listarr) MemoryFail (funcname);

				for (int ilist=0;ilist<nlistloc;ilist++) {
					int iblk = listloc[ilist];
					if (_mvm.blk2cpu[iblk] == childcpuid) {
						listarr[nlistarr] = iblk;
						nlistarr++;
					};
				};

//				cout << " Send diag data to " << childcpuid << endl;
//				OutArr(cout," List send diag",nlistarr,listarr);
				SendDiag (_tree.comm, nlistarr, listarr, childcpuid, childnode, _gmtra);

				delete [] listarr;

			};

		};

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Find supernodes with small singular values
//========================================================================================
void CGSMatrixRS::FindSups (ofstream &_fout, int _nthresh, double *_thresh, int *_iorder) { // Find supernodes with small singular values

	const char *funcname = "FindSups";

// Init work structures

	int iblk = listb[0];

	int blamxloc = mtrarr[iblk].blamx;

	CFctR fct (nblksr, nsupr, blamxloc);

	CFctDiagR fctdiag (nblksr);

// Allocate working arrays

	int *nlisttot;
	int *listtot;
	int *listtotnew;
	int *nlistarr;
	int **listarr;

	nlisttot = new int [_nthresh];
	if (!nlisttot) MemoryFail (funcname);
	listtot = new int [_nthresh*nsupr];
	if (!listtot) MemoryFail (funcname);
	listtotnew = new int [_nthresh*nsupr];
	if (!listtotnew) MemoryFail (funcname);

	nlistarr = new int [_nthresh];
	if (!nlistarr) MemoryFail (funcname);
	listarr = new int * [_nthresh];
	if (!listarr) MemoryFail (funcname);

// Main cycle over the blocks

	int i;
	for (i=0;i<_nthresh;i++) nlisttot[i] = 0;

	for (int ilist=0;ilist<nlist;ilist++) {

// Find lists of supernodes

		iblk = listb[ilist];

		ReadBlock (iblk);

		fctdiag.FindSups (iblk, _nthresh, _thresh, *this, mtrarr[iblk], nlistarr, listarr, fct);

		FreeBlock (iblk);

// Store data

		for (i=0;i<_nthresh;i++) {
			int j = nlistarr[i];
			int ibs = nlisttot[i];
			for (int kii=0;kii<j;kii++) {
				listtot[i*nsupr+ibs+kii] = listarr[i][kii];
			};
			nlisttot[i] += j;
		};

// Free work arrays

		for (i=0;i<_nthresh;i++) {
			delete [] listarr[i];
		};

	};

// Exchange data among processors

	int nproc, myid;

	nproc = CMPIExchange::GetNprocMPIGlobal ();
	myid = CMPIExchange::GetMyidMPIGlobal ();

	for (i=0;i<nsupr*_nthresh;i++) listtotnew[i] = 0;

	int ithresh;
	for (ithresh=0;ithresh<_nthresh;ithresh++) {
		int j = nlisttot[ithresh];
		for (i=0;i<j;i++) {
			int jj=listtot[ithresh*nsupr+i];
			listtotnew[ithresh*nsupr+jj] = 1;
		};
	};

	AddSparseIArray (myid, nproc, nsupr*_nthresh, listtotnew);

	for (ithresh=0;ithresh<_nthresh;ithresh++) {
		int icount = 0;
		for (i=0;i<nsupr;i++) {
			if (listtotnew[ithresh*nsupr+i] == 1) {
				int iord = _iorder[i];
				listtot[ithresh*nsupr+icount] = iord;
				icount++;
			};
		};
		nlisttot[ithresh] = icount;
		qsort (listtot+ithresh*nsupr, icount, sizeof(int), compint);
	};

// Output sorted reordered lists

	if (myid == 0) {
		cout  << " The list of supernodes for each threshold value: " << endl;
		_fout << " The list of supernodes for each threshold value: " << endl;
		for (ithresh=0;ithresh<_nthresh;ithresh++) {
			cout  << " Ithr = " << ithresh << " Value = " << _thresh[ithresh] << endl;
			_fout << " Ithr = " << ithresh << " Value = " << _thresh[ithresh] << endl;
			OutArr (cout, " SupList = ",nlisttot[ithresh],listtot+ithresh*nsupr);
			OutArr (_fout," SupList = ",nlisttot[ithresh],listtot+ithresh*nsupr);
		};
	};

// Free work arrays

	delete [] nlisttot;
	delete [] listtot;
	delete [] listtotnew;
	delete [] nlistarr;
	delete [] listarr;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif// 
