//------------------------------------------------------------------------------------------------
// File: ich2gcs.cpp
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

#include "tree.h"
#include "fct.h"
#include "fctdiag.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "slvparam.h"
#include "mvm.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CGSMatrixCS: Compute second order Incomplete Cholessky decomposition in parallel mode
//========================================================================================
void CGSMatrixCS::Ich2Schur (ofstream &_fout, const CTree &_tree, const CSlvParam &_param, // Compute second order Incomplete Cholessky decomposition in parallel mode
								CGSMatrixCS &_mtrau,
								CGSMatrixCS &_mtru) { 

	const char *funcname = "Ich2Schur";

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::SynchronizeMPI (pcomm);
	};

// Check the data on entry

	int numprocs;

	CMPIComm pcomm = _tree.GetComm ();
	numprocs = pcomm.GetNproc ();

	if (_tree.nproc != numprocs) {
		throw " Parallel Ilu2 function called with wrong number of processors ";
	};

// Statistics data

	double ops, densu, densr, denss, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	double timescl=0.0e0, timesclex=0.0e0;
	double timeload=0.0e0, timeunload=0.0e0;
	double timeini=0.0e0, timeupdlocal=0.0e0, timeupdnonlocal=0.0e0, timefct=0.0e0;
	double timeprep=0.0e0, timesend=0.0e0, timerecv=0.0e0;
	int nzstot=0;

// Init time measurement

	time0 = clock ();

// Compute bl2cpu

	int *bl2cpu;
	int *listsub;

	bl2cpu = new int [_mtrau.nblksr];
	if (!bl2cpu) MemoryFail (funcname);
	listsub = new int [_mtrau.nblksr];
	if (!listsub) MemoryFail (funcname);

	_tree.Block2Cpu (bl2cpu);

// Compute block sparsity structure of the matrix

	CSMatrix **pmtrarr;

	pmtrarr = new CSMatrix * [_mtrau.nblksr];
	if (!pmtrarr) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<_mtrau.nblksr;iblk++) {
		pmtrarr[iblk] = &_mtrau.mtrarr[iblk];
	};

	int *iab, *jab;

	const CGSMatrix *pmtrau;

	pmtrau = &_mtrau;

	pmtrau->MtrBlockSparsity (_tree, bl2cpu, pmtrarr, iab, jab);

// Compute extended block sparsity structure

	int *iabext, *jabext;

	pmtrau->ExtendBlockSparsity (_mtrau.nblksr, 
									iab, jab, iabext, jabext);

// Compute transposed extended block sparsity structure

	int *iabextt, *jabextt;

	int nzjext = iabext[_mtrau.nblksr];

	iabextt = new int [_mtrau.nblksr+1];
	if (!iabextt) MemoryFail (funcname);
	jabextt = new int [nzjext];
	if (!jabextt) MemoryFail (funcname);

	TransposeMatrix (_mtrau.nblksr, iabext, jabext,
						iabextt, jabextt);

// Determine blamx

	int blamxl=0;
	int blai;

	int isup;
	for (isup=0;isup<_mtrau.nsupr;isup++) {
		blai = _mtrau.sprndr[isup+1]-_mtrau.sprndr[isup];
		if (blai > blamxl) blamxl = blai;
	};

// Allocate Fct data

	CFctC fct (_mtrau.nblksr, _mtrau.nsupr, blamxl); 

// Allocate FctDiag data

	CFctDiagC fctdiag (_mtrau.nblksr);

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"DiagLLt_",_tree.myid,"_");

		fctdiag.SetupFiles (_param.ooctyp, strbuff);

	};

// Create mvm structure

	CMvm mvm (false, _tree, _mtrau);

// Compute the local scaling of the matrix

	for (iblk=0;iblk<_mtrau.nblksr;iblk++) {

		if (bl2cpu[iblk] == _tree.myid) {

			_mtrau.ReadBlock (iblk);

			time2 = clock ();

			fctdiag.Ich2ScaleBlkRow (iblk, _param.sclpttype, _param.sclmin,
										_mtrau, _mtrau.mtrarr[iblk], fct);

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
		cout  << " PointScaling" << fctdiag.hyst2;
		_fout << " PointScaling" << fctdiag.hyst2;
		cout  << " Scaling" << fctdiag.hyst;
		_fout << " Scaling" << fctdiag.hyst;
	};

// Exchange the scaling data from the bordering

	time2 = clock ();

	fctdiag.ExchangeScaling (_tree, mvm, _mtrau);

	time3 = clock ();

	timesclex += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Create dummy matrix

	CSMatrixCS mtrdummy;

// Init mtrl and mtru

	_mtru = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

		sprintf(strbuff, "%s%s%d%s",_param.path,"UMatr_",_tree.myid,"_");

		_mtru.SetupFiles (_param.ooctyp, strbuff);

	};

	_mtru.ivect = _param.sclpttype;

// Init gmtrl2 and gmtru2

	CGSMatrixCS gmtru2;
	CGSMatrixCS gmtru3;

	gmtru2 = _mtrau;

	gmtru3 = _mtrau;

	if (_param.ooctyp > 0) {

		char strbuff[256];

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
		int isupbeg = _mtrau.blksc[iblk];
		int isupend = _mtrau.blksc[iblk+1]-1;

		gmtru2.mtrarr[iblk] = fct.CreateBlock (isupbeg, isupend, _mtrau);

		gmtru2.WriteBlock (iblk);
		gmtru2.FreeBlock  (iblk);

	};

// Main cycle over the block rows according to the tree

	CSMatrixCS mtrutemp, mtru2temp, mtrurecv, mtrusend;

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

								_mtrau.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk+1, mtrurecv);

								time3 = clock ();

								timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								_mtrau.SplitMatrix (nlistloc, listloc,
													mtrurecv, 
													gmtru3.mtrarr);

								mtrurecv = mtrdummy;

// Update and free

								for (ilist=0;ilist<nlistloc;ilist++) {

									iblk = listloc[ilist];

									gmtru2.ReadBlock (iblk);

									time2 = clock ();

									mtru2temp = fct.AddBlocks (_mtrau, gmtru2.mtrarr[iblk], gmtru3.mtrarr[iblk]);

									time3 = clock ();

									timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

									gmtru2.mtrarr[iblk] = mtru2temp;

									gmtru3.mtrarr[iblk] = mtrdummy;

									mtru2temp = mtrdummy;

									mtrurecv = mtrdummy;

									gmtru2.WriteBlock (iblk);
									gmtru2.FreeBlock  (iblk);

								};

							};

						} else {

// Receive and update

							for (ilist=0;ilist<nlistloc;ilist++) {

								iblk = listloc[ilist];

								time2 = clock ();

								_mtrau.ReceiveMatrix (_tree.comm, childcpucurr, 2*iblk+1, mtrurecv);

								time3 = clock ();

								timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								gmtru2.ReadBlock (iblk);

								time2 = clock ();

								mtru2temp = fct.AddBlocks (_mtrau, gmtru2.mtrarr[iblk], mtrurecv);

								time3 = clock ();

								timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

								gmtru2.mtrarr[iblk] = mtru2temp;

								mtru2temp = mtrdummy;

								mtrurecv = mtrdummy;

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

					_mtrau.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk+1, mtrurecv);

					time3 = clock ();

					timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					_mtrau.SplitMatrix (nlistloc, listloc,
										mtrurecv, 
										gmtru2.mtrarr);

					mtrurecv = mtrdummy;

				};

			} else {

				for (ilist=0;ilist<nlistloc;ilist++) {

					iblk = listloc[ilist];

					time2 = clock ();

					_mtrau.ReceiveMatrix (_tree.comm, fathercpuid2, 2*iblk+1, mtrurecv);

					time3 = clock ();

					timerecv += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					gmtru2.mtrarr[iblk] = mtrurecv;

					mtrurecv = mtrdummy;

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
										_mtrau, fct);

		time3 = clock ();

		timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Perform computations with blocks of the current node taking into account Schur data

		if (inodetype != exchangenode) {

			for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {

				cout << " Myid = " << _tree.myid << " Iblk = " << iblk << endl;
//				_fout << " Myid = " << _tree.myid << " Iblk = " << iblk << endl;

// Load diagonal data into fct structure

				int nlstblk = iabext[iblk+1]-iabext[iblk];
				int ibs = iabext[iblk];

				time2 = clock ();

				fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtrau, fct);

				time3 = clock ();

				timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Init current block row of the matrix

				_mtrau.ReadBlock (iblk);

				time2 = clock ();

				fct.Ich2InitBlkRow (iblk, _param.dshift, _param.dshiftim,
									_mtrau, _mtrau.mtrarr[iblk], 
									mtrutemp);
//				_fout << " Iblk = " << iblk << " Data after init " << mtrutemp << endl;

				time3 = clock ();

				timeini += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				_mtrau.FreeBlock (iblk);

// Update current block row by the Schur data

				if (inodecurr != inodesave) {

					gmtru2.ReadBlock (iblk);

					time2 = clock ();

					mtru2temp = fct.AddBlocks (_mtrau, mtrutemp, gmtru2.mtrarr[iblk]);

					time3 = clock ();

					timeupdnonlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					mtrutemp = mtru2temp;

					mtru2temp = mtrdummy;

					gmtru2.mtrarr[iblk] = mtrdummy;

				};

// Update current block row by the previous ones

				if (_tree.nodes[inodecurr].nchilds == 1 && inodetype != compnode2) {

					for (int j=iabextt[iblk];j<iabextt[iblk+1]-1;j++) {

						int jblk = jabextt[j];

//						cout << " update Jblk = " << jblk << endl;
//						_fout << " update Jblk = " << jblk << endl;

						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						fct.Ich2UpdateBlkRow (iblk, _param,
												_mtrau, mtrutemp,
												gmtru2.mtrarr[jblk],
												mtru2temp);

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtru2.FreeBlock (jblk);

						mtrutemp = mtru2temp;

						mtru2temp = mtrdummy;

					};

				} else {

					for (int jblk=_tree.nodes[inodecurr].indbeg;jblk<iblk;jblk++) {

						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						fct.Ich2UpdateBlkRow (iblk, _param,
												_mtrau, mtrutemp,
												gmtru2.mtrarr[jblk],
												mtru2temp);

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtru2.FreeBlock (jblk);

						mtrutemp = mtru2temp;

						mtru2temp = mtrdummy;

					};

				};

// Print Schur sparsity if necessary

				if (false && inodecurr != inodesave) {

					int isupbeg = _mtrau.blksc[iblk];
					int isupend = _mtrau.blksc[iblk+1]-1;

					CSMatrix aschur;

					aschur = mtrutemp.CSMatrix::FilterMatrix(isupbeg,isupend);

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

				fct.Ich2FctBlkRow (iblk, _param,
									_mtrau, mtrutemp,
									_mtru.mtrarr[iblk], 
									gmtru2.mtrarr[iblk]);

				time3 = clock ();

				timefct += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

				_mtru.WriteBlock (iblk);
				_mtru.FreeBlock  (iblk);

				gmtru2.WriteBlock (iblk);
				gmtru2.FreeBlock  (iblk);

				mtrutemp = mtrdummy;

// Unload diagonal data from fct structure

				nlstblk = iabext[iblk+1]-iabext[iblk];
				ibs = iabext[iblk];

				time2 = clock ();

				fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtrau, fct);

				time3 = clock ();

				timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			};

// Scan second order matrices and filter unnecessary data

			int iblkbeg = _tree.nodes[inodecurr].indbeg;
			int iblkend = _tree.nodes[inodecurr].indend;
			int isupend = _mtrau.blksc[iblkend+1]-1;

			for (iblk=_tree.nodes[inodecurr].indbeg;iblk<=_tree.nodes[inodecurr].indend;iblk++) {

				gmtru2.ReadBlock    (iblk);
				gmtru2.FilterBlock  (iblk,isupend);
				gmtru2.ReWriteBlock (iblk);
				gmtru2.FreeBlock    (iblk);

			};

// Compute the block sparsity structure for the whole bordering

			int *iabbord, *jabbord;

			pmtrau = &_mtrau;

			int iblk;
			for (iblk=0;iblk<_mtrau.nblksr;iblk++) {
				pmtrarr[iblk] = &gmtru2.mtrarr[iblk];
			};

			pmtrau->MtrBlockSparsity (pmtrarr, iblkbeg, iblkend, // Compute block sparsity of the matrix
										iabbord, jabbord);

// Compute transposed block sparsity structure

			int *iabbordtr, *jabbordtr;

			int nzjbord = iabbord[_mtrau.nblksr];

			iabbordtr = new int [_mtrau.nblksr+1];
			if (!iabbordtr) MemoryFail (funcname);
			jabbordtr = new int [nzjbord];
			if (!jabbordtr) MemoryFail (funcname);

			TransposeMatrix (_mtrau.nblksr, iabbord, jabbord,
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
													_mtrau, fct);

					time3 = clock ();

					timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Read data if necessary

					gmtru2.ReadBlock (iblk);

// Update

					for (int j=iabbordtr[iblk];j<iabbordtr[iblk+1];j++) {

						int jblk = jabbordtr[j];

						gmtru2.ReadBlock (jblk);

						time2 = clock ();

						fct.Ich2UpdateBlkRow (iblk, _param,
												_mtrau, gmtru2.mtrarr[iblk],
												gmtru2.mtrarr[jblk],
												mtru2temp);

						time3 = clock ();

						timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						gmtru2.FreeBlock (jblk);

						gmtru2.mtrarr[iblk] = mtru2temp;

						mtru2temp = mtrdummy;

					};

// Store data to the disk if necessary

					gmtru2.WriteBlock (iblk);
					gmtru2.FreeBlock  (iblk);

// Unload diagonal data

					nlstblk = iabext[iblk+1]-iabext[iblk];
					ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
												_mtrau, fct);

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
				gmtru2.mtrarr[iblk] = mtrdummy;
			};

		} else {

// Perform filtering computations in the exchange node

			int iblkbegfather = _tree.nodes[fathernode].indbeg;

			for (ichild=0;ichild<_tree.nodes[inodecurr].nchilds2;ichild++) {

				childnodecurr = _tree.nodes[inodecurr].childs2[ichild];

				int iblkbegloc = _tree.nodes[childnodecurr].indbeg;
				int iblkendloc = _tree.nodes[childnodecurr].indend;

				int isupbeg = _mtrau.blksr[iblkbegloc];
				int isupend = _mtrau.blksr[iblkendloc+1]-1;
				int isuprest = _mtrau.blksr[iblkbegfather];

				double theta2 = _param.theta / 2.0e0;

				for (iblk=iblkbegloc;iblk<=iblkendloc;iblk++) {

// Load diagonal data into fct structure

					int nlstblk = iabext[iblk+1]-iabext[iblk];
					int ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2LoadDiagonalData (iblk, nlstblk, jabext+ibs,
													_mtrau, fct);

					time3 = clock ();

					timeload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

// Filter block data

					gmtru2.ReadBlock    (iblk);

					int nzaini = gmtru2.mtrarr[iblk].nzatot;

					time2 = clock ();

					gmtru2.FilterSchurBlock (iblk, isupbeg, isupend, isuprest, theta2, fct);

					time3 = clock ();

					timeupdlocal += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

					int nzafin = gmtru2.mtrarr[iblk].nzatot;
					cout << " Filtering Iblk = " << iblk << " nzaini = " << nzaini << " nzafin = " << nzafin << endl;

					nzstot += nzaini-nzafin;

					gmtru2.ReWriteBlock (iblk);
					gmtru2.FreeBlock    (iblk);

// Unload diagonal data from fct structure

					nlstblk = iabext[iblk+1]-iabext[iblk];
					ibs = iabext[iblk];

					time2 = clock ();

					fctdiag.Ilu2UnloadDiagonalData (iblk, nlstblk, jabext+ibs,
													_mtrau, fct);

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

						_mtrau.CombineMatrices (nlistloc, listloc, 
												gmtru2.mtrarr, mtrusend);

						iblk = listloc[0];

						time2 = clock ();

						_mtrau.SendMatrix (_tree.comm, fathercpuid, 2*iblk+1, mtrusend);

						time3 = clock ();

						timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

						mtrusend = mtrdummy;

						for (ilist=0;ilist<nlistloc;ilist++) {

							iblk = listloc[ilist];

							gmtru2.mtrarr[iblk] = mtrdummy;

						};

					};

				} else {

					for (ilist=0;ilist<nlistloc;ilist++) {

						iblk = listloc[ilist];

						gmtru2.ReadBlock (iblk);

						time2 = clock ();

						_mtrau.SendMatrix (_tree.comm, fathercpuid, 2*iblk+1, gmtru2.mtrarr[iblk]);

						time3 = clock ();

						timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

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

							_mtrau.CombineMatrices (nlistsub, listsub, 
													gmtru2.mtrarr, mtrusend);

							iblk = listsub[0];

							time2 = clock ();

							_mtrau.SendMatrix (_tree.comm, childcpuid, 2*iblk+1, mtrusend);

							time3 = clock ();

							timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

							mtrusend = mtrdummy;

							for (ilist=0;ilist<nlistsub;ilist++) {

								iblk = listsub[ilist];

								gmtru2.mtrarr[iblk] = mtrdummy;

							};

						};

					} else {

						for (ilist=0;ilist<nlistloc;ilist++) {

							iblk = listloc[ilist];

							if (mvm.blk2cpu[iblk] == childcpuid) {

//								cout << " send iblk = " << iblk << endl;
								gmtru2.ReadBlock (iblk);

								time2 = clock ();

								_mtrau.SendMatrix (_tree.comm, childcpuid, 2*iblk+1, gmtru2.mtrarr[iblk]);

								time3 = clock ();

								timesend += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

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
										_mtrau, fct);

		time3 = clock ();

		timeunload += (double) (time3-time2) / (double) CLOCKS_PER_SEC;

NextNode:;

	};

// Store point/block scaling data in U structure

	if (_mtru.ivect < 2) {
		fctdiag.StorePointScaling (_mtru);
	} else {
		fctdiag.StoreBlockScaling ('U', _mtru);
	};

// Close files

	if (_param.ooctyp > 0) {

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
		cout  << " Factorization" << fct.hyst;
		_fout << " Factorization" << fct.hyst;
	};

// Output timing results

	if (_tree.myid == _tree.rootcpu) {
		cout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
		cout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s." << endl;
		cout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
		_fout << " TScl = " << timescl << " s.; TSclEx = " << timesclex << " s.; TLoad = " << timeload << " s.; TUnload = " << timeunload << " s." << endl;
		_fout << " TIni = " << timeini << " s.; TUpdLoc = " << timeupdlocal << " s.; TUpdNonLoc = " << timeupdnonlocal << " s.; TFct = " << timefct << " s." << endl;
		_fout << " TPrep = " << timeprep << " s.; TSend = " << timesend << " s.; TRecv = " << timerecv << " s." << endl;
	};

// Output Ich2 factorization statistics

	int nzatotal = 0;
	for (iblk=0;iblk<_mtrau.nblksr;iblk++) {
		nzatotal += _mtrau.mtrarr[iblk].nzatot;
	};

	int nzmema = 2 * nzatotal;

	for (ilist=0;ilist<_mtrau.nlist;ilist++) {
		iblk = _mtrau.listb[ilist];
		for (isup=_mtrau.blksr[iblk];isup<_mtrau.blksr[iblk+1];isup++) {
			blai = _mtrau.sprndr[isup+1]-_mtrau.sprndr[isup];
			nzmema -= blai*blai;
		};
	};

	int nzutot = 0;
	for (iblk=0;iblk<_mtrau.nblksr;iblk++) {
		nzutot += _mtru.mtrarr[iblk].nzatot;
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
		for (int i=0;i<7;i++) arr1[i] = arr0[i];
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

		cout  << " Ich2 preconditioner generation statistics: " << endl;
		if (nmodtot != 0) {
			cout  << "     NmodS = " << nmodtot << endl;
		};
		cout  << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		_fout << " Ich2 preconditioner generation statistics: " << endl;
		if (nmodtot != 0) {
			_fout << "     NmodS = " << nmodtot << endl;
		};
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " %  DensS = " << denss << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

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

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
