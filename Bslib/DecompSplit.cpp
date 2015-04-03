//------------------------------------------------------------------------------------------------
// File: DecompSplit.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "globals.h"
#include "ExchangeMPI.h"
#include "Decomp.h"
#include "CpuTree.h"
#include "tree.h"
#include "smatrix.h"
#include "SolverPar.h"

using namespace std;

#ifdef __FV_2004__
//using namespace flowvision::EqnSolver;
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel partitioning of the matrix into many parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::SplitBlockMatrixRecursive (CMPIComm &_comm, int _nparts, // Compute in parallel partitioning of the matrix into many parts taking into account blocks partitioning
														int _nblks, int *_blks, 
														int *_blk2cpu,
														int *_ordernd, int *_orderblk) const {

	const char *funcname = "SplitBlockMatrixRecursive_00";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that the number of parts is a power of 2

	int nparts_pwr2 = 1;
	while (nparts_pwr2*2 <= _nparts) nparts_pwr2 *= 2;

	if (nparts_pwr2 != _nparts) {
		cout << " throw CSMatrix::SplitBlockMatrixRecursive: the number of nodes is not a power of 2 " << endl;
		throw " CSMatrix::SplitBlockMatrixRecursive: the number of nodes is not a power of 2 ";
	};

// Check that the number of nproc is a power of 2

	int nproc_pwr2 = 1;
	while (nproc_pwr2*2 <= nproc) nproc_pwr2 *= 2;

	if (nproc_pwr2 != nproc) {
		cout << " throw CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is not a power of 2 " << endl;
		throw " CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is not a power of 2 ";
	};

	if (nproc > _nparts) {
		cout << " throw CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is larger than the number of parts " << endl;
		throw " CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is larger than the number of parts ";
	};

	int nlev = 0;

	int np = 1;

	while (np*2 <= nproc) {
		np *= 2;
		nlev++;
	};

// Allocate the memory that supports the levels computations

	int nblkmax = 3*nlev+7;

	int *blkstot;
	int *blkstotblk;

	blkstot = new int [nblkmax+1];
	if (!blkstot) MemoryFail (funcname);
	blkstotblk = new int [nblkmax+1];
	if (!blkstotblk) MemoryFail (funcname);

	int i;

	for (i=0;i<=nblkmax;i++) blkstot[i] = 0;
	for (i=0;i<=nblkmax;i++) blkstotblk[i] = 0;

	int **pporderlev, *nblkslev, **ppblkslev;
	int **pporderlevblk;
	int *nblkslevsups, **ppblkslevsups;
	int *nsupslev, **ppsupslev;

	pporderlev = new int * [nlev+2];
	if (!pporderlev) MemoryFail (funcname);
	nblkslev = new int [nlev+2];
	if (!nblkslev) MemoryFail (funcname);
	ppblkslev = new int * [nlev+2];
	if (!ppblkslev) MemoryFail (funcname);
	pporderlevblk = new int * [nlev+2];
	if (!pporderlevblk) MemoryFail (funcname);
	nblkslevsups = new int [nlev+2];
	if (!nblkslevsups) MemoryFail (funcname);
	ppblkslevsups = new int * [nlev+2];
	if (!ppblkslevsups) MemoryFail (funcname);
	nsupslev = new int [nlev+2];
	if (!nsupslev) MemoryFail (funcname);
	ppsupslev = new int * [nlev+2];
	if (!ppsupslev) MemoryFail (funcname);

	CMPIComm *commlev;

	commlev = new CMPIComm [nlev+1];
	if (!commlev) MemoryFail (funcname);

// Create the whole set of communicators

	int *listcpu;

	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);

	commlev[0] = _comm;

	int ilev, nprocloc, myidloc, nproc1;
	CMPIComm comm1, comm2;

	for (ilev=1;ilev<=nlev;ilev++) {
		nprocloc = commlev[ilev-1].GetNproc ();
		myidloc = commlev[ilev-1].GetMyid ();
		nproc1 = nprocloc / 2;
		for (i=0;i<nproc1;i++) listcpu[i] = i;
		comm1 = commlev[ilev-1].CreateComm (nproc1, listcpu);
		for (i=nproc1;i<nproc;i++) listcpu[i-nproc1] = i;
		comm2 = commlev[ilev-1].CreateComm (nproc1, listcpu);
		if (myidloc < nproc1) {
			commlev[ilev] = comm1;
		} else {
			commlev[ilev] = comm2;
		};
	};

// Init zero level data

	int *blksloc;

	blksloc = new int [nproc+1];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) blksloc[i] = 0;

	blksloc[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksloc, blksloc);

	for (i=0;i<nproc;i++) blksloc[i+1] = blksloc[i]+blksloc[i+1];

	if (nlist > 0) {
		int indbeg = list[0];
		int indend = list[nlist-1];
		if (indbeg != blksloc[myid] || indend != blksloc[myid+1]-1) {
			cout << " throw CSMatrix::SplitBlockMatrixRecursive: wrong part of the matrix is considered" << endl;
			throw " CSMatrix::SplitBlockMatrixRecursive: wrong part of the matrix is considered";
		};
	};

	nblkslev[0] = nproc;
	ppblkslev[0] = blksloc;

	blksloc = new int [nproc+1];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) blksloc[i] = 0;

	blksloc[myid+1] = _nblks;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksloc, blksloc);

	for (i=0;i<nproc;i++) blksloc[i+1] = blksloc[i]+blksloc[i+1];

	nblkslevsups[0] = nproc;
	ppblkslevsups[0] = blksloc;

	blksloc = new int [_nblks+1];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<_nblks+1;i++) blksloc[i] = _blks[i];

	nsupslev[0] = _nblks;
	ppsupslev[0] = blksloc;

	int ni = nlist;

// Main cycle over the levels of the tree

	CSMatrix mtrlev;
	CSMatrix mtrdummy;

	mtrlev = *this;

	int igroup, myid1, ibegloc, iendloc, ibegsub, iendsub;
	int nlistloc, nz, j, jj, nloc;
	int *plist, *pia, *pja, *ialoc;

	for (ilev=0;ilev<nlev;ilev++) {

// Detect the group number for the current processor on the next level

		nprocloc = commlev[ilev].GetNproc ();
		myidloc = commlev[ilev].GetMyid ();

		nproc1 = nprocloc / 2;

		if (myidloc < nproc1) {
			igroup = 0;
		} else {
			igroup = 1;
		};

		myid1 = commlev[ilev+1].GetMyid ();

// Compute local ordering and partitioning, and reorder the matrix indices

		mtrlev.SplitBlockMatrix (commlev[ilev], 
											nsupslev[ilev], ppsupslev[ilev],
											pporderlevblk[ilev], blkstotblk+ilev*3,
											pporderlev[ilev], blkstot+ilev*3,
											nproc1, igroup, myid1,
											nblkslev[ilev+1], ppblkslev[ilev+1],
											nblkslevsups[ilev+1], ppblkslevsups[ilev+1],
											nsupslev[ilev+1], ppsupslev[ilev+1]);

// Reorder the matrix in-place

		mtrlev.ChangeIndices (commlev[ilev], pporderlev[ilev]);

// Exchange matrix data according to the partitioning

		ibegloc = blkstot[ilev*3+igroup];
		iendloc = blkstot[ilev*3+igroup+1]-1;

		ibegsub = ppblkslev[ilev+1][myid1];
		iendsub = ppblkslev[ilev+1][myid1+1]-1;

		ibegsub += ibegloc;
		iendsub += ibegloc;

		CSMatrix mtrnew;

		mtrnew = mtrlev.GetSubmatrixViaIntervals (commlev[ilev],
																1, &ibegsub, &iendsub);

// Filter and shift the rows/columns indices

		nlistloc = mtrnew.GetNlist ();
		plist = mtrnew.GetList ();
		pia = mtrnew.GetIa ();
		pja = mtrnew.GetJa ();

		ialoc = new int [nlistloc+1];
		if (!ialoc) MemoryFail (funcname);

		ialoc[0] = 0;
		nz = 0;

		for (i=0;i<nlistloc;i++) {
			plist[i] -= ibegloc;
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				if (jj >= ibegloc && jj <= iendloc) {
					pja[nz] = jj-ibegloc;
					nz++;
				};
			};
			ialoc[i+1] = nz;
		};

		for (i=0;i<=nlistloc;i++) pia[i] = ialoc[i];

		nloc = ppblkslev[ilev+1][nblkslev[ilev+1]];

		mtrnew.SetM (nloc);
		mtrnew.SetN (nloc);
		mtrnew.SetNsupr (nloc);
		mtrnew.SetNsupcBase (nloc);
		mtrnew.SetNzja (nz);

		delete [] ialoc;

		mtrlev = mtrnew;

		mtrnew = mtrdummy;

	};

// Reorder the last submatrix

	int npartsloc = _nparts / nproc;
	int nsupsloc = nsupslev[nlev];

	nloc = mtrlev.GetN ();

	pporderlev[nlev] = new int [nloc];
	if (!pporderlev[nlev]) MemoryFail (funcname);
	pporderlevblk[nlev] = new int [nsupsloc];
	if (!pporderlevblk[nlev]) MemoryFail (funcname);
	ppblkslev[nlev+1] = new int [npartsloc+1];
	if (!ppblkslev[nlev+1]) MemoryFail (funcname);
	ppblkslevsups[nlev+1] = new int [npartsloc+1];
	if (!ppblkslevsups[nlev+1]) MemoryFail (funcname);

	int *blk2cpuloc;

	blk2cpuloc = new int [nsupsloc];
	if (!blk2cpuloc) MemoryFail (funcname);
/*
	if (true) {

		int *pialoc = mtrlev.GetIa ();
		int *pjaloc = mtrlev.GetJa ();

		int icheck = 1;

		for (int ilist=0;ilist<nsupsloc;ilist++) {
			int ielem = ilist;
			int icheckloc = -1;
			for (int j=pialoc[ielem];j<pialoc[ielem+1];j++) {
				int jj = pjaloc[j];
				if (jj != ielem) {
					icheckloc = 1;
				};
			};
			if (icheckloc == -1) {
				cout << " Wrong matrix elem = " << ielem << endl;
				icheck = -1;
			};
		};

	};
*/
	mtrlev.SplitBlockMatrixRecursive (npartsloc,
													nsupslev[nlev], ppsupslev[nlev], 
													blk2cpuloc);
/*
	if (false) {

		int *pialoc = mtrlev.GetIa ();
		int *pjaloc = mtrlev.GetJa ();

		for (int ipart=0;ipart<npartsloc;ipart++) {

			int icheck = 1;

			for (int ilist=0;ilist<nsupsloc;ilist++) {
				int ielem = ilist;
				if (blk2cpuloc[ielem] == ipart) {
					int icheckloc = -1;
					for (int j=pialoc[ielem];j<pialoc[ielem+1];j++) {
						int jj = pjaloc[j];
						if (jj != ielem) {
							if (blk2cpuloc[jj] == ipart) {
								icheckloc = 1;
							};
						};
					};
					if (icheckloc == -1) {
						cout << " Wrong part = " << ipart << " Wrong elem = " << ielem << endl;
						icheck = -1;
					};
				};
			};
		};
	};
*/
	nblkslev[nlev+1] = npartsloc;
	nblkslevsups[nlev+1] = npartsloc;

	for (i=0;i<=npartsloc;i++) ppblkslevsups[nlev+1][i] = 0;

	int ipart;

	for (i=0;i<nsupslev[nlev];i++) {
		ipart = blk2cpuloc[i];
		ppblkslevsups[nlev+1][ipart+1]++;
	};

	for (i=0;i<npartsloc;i++) ppblkslevsups[nlev+1][i+1] = ppblkslevsups[nlev+1][i]+ppblkslevsups[nlev+1][i+1];
	for (i=0;i<npartsloc;i++) ppblkslev[nlev+1][i] = ppblkslevsups[nlev+1][i];

	int k;

	for (i=0;i<nsupslev[nlev];i++) {
		ipart = blk2cpuloc[i];
		k = ppblkslev[nlev+1][ipart];
		pporderlevblk[nlev][i] = k;
		ppblkslev[nlev+1][ipart]++;
	};

	int *sprndso;

	sprndso = new int [nsupsloc+1];
	if (!sprndso) MemoryFail (funcname);

	sprndso[0] = 0;

	int iord;

	for (i=0;i<nsupslev[nlev];i++) {
		iord = pporderlevblk[nlev][i];
		ni = ppsupslev[nlev][i+1]-ppsupslev[nlev][i];
		sprndso[iord+1] = ni;
	};

	for (i=0;i<nsupslev[nlev];i++) sprndso[i+1] = sprndso[i]+sprndso[i+1];

	int jold, jnew;

	for (i=0;i<nsupslev[nlev];i++) {
		iord = pporderlevblk[nlev][i];
		ni = ppsupslev[nlev][i+1]-ppsupslev[nlev][i];
		for (j=0;j<ni;j++) {
			jold = ppsupslev[nlev][i]+j;
			jnew = sprndso[iord]+j;
			pporderlev[nlev][jold] = jnew;
		};
	};

// Compute extended blocks partitionings

	int *blkstotext;
	int *blkstotblkext;

	blkstotext = new int [nlev*4];
	if (!blkstotext) MemoryFail (funcname);
	blkstotblkext = new int [nlev*4];
	if (!blkstotblkext) MemoryFail (funcname);

	for (i=0;i<nlev;i++) {
		for (j=0;j<3;j++) {
			blkstotext[i*4+j] = blkstot[i*3+j];;
			blkstotblkext[i*4+j] = blkstotblk[i*3+j];;
		};
		blkstotext[i*4+3] = blkstotext[i*4+2];;
		blkstotblkext[i*4+3] = blkstotblkext[i*4+2];;
	};

// Backward cycle with the ordering restoration

	if (nproc > 1) {

		for (ilev=nlev-1;ilev>=0;ilev--) {

			CSMatrix::RestoreOrderNd (ilev, nlev,
												commlev, 
												blkstotext, pporderlev, nblkslev, ppblkslev);

		};

		for (ilev=nlev-1;ilev>=0;ilev--) {

			CSMatrix::RestoreOrderNd (ilev, nlev,
												commlev, 
												blkstotblkext, pporderlevblk, nblkslevsups, ppblkslevsups);

		};

	};

	for (i=0;i<nlist;i++) _ordernd[i] = pporderlev[0][i];

	for (i=0;i<_nblks;i++) _orderblk[i] = pporderlevblk[0][i];

// Compute the final blksblk array

	int nblksblk = _nparts;

	int *blksblk;

	blksblk = new int [nblksblk+1];
	if (!blksblk) MemoryFail (funcname);

	for (i=0;i<nblksblk+1;i++) blksblk[i] = 0;

	for (i=0;i<npartsloc;i++) {
		blksblk[myid*npartsloc+i+1] = ppblkslevsups[nlev+1][i+1]-ppblkslevsups[nlev+1][i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nparts+1, blksblk, blksblk);

	for (i=0;i<nblksblk;i++) blksblk[i+1] = blksblk[i]+blksblk[i+1];

// Compute distribution array

	int icpuprev = 0;

	int iblk, iblknew, icpu, icpubeg, icpuend;

	for (iblk=0;iblk<_nblks;iblk++) {
		iblknew = _orderblk[iblk];
		if (iblknew >= blksblk[icpuprev] && iblknew < blksblk[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			icpubeg = 0;
			icpuend = _nparts-1;
			while (icpubeg != icpuend) {
				icpu = (icpubeg+icpuend)/2;
				if (iblknew >= blksblk[icpu] && iblknew < blksblk[icpu+1]) {
					icpubeg = icpu;
					icpuend = icpu;
				} else if (iblknew < blksblk[icpu]) {
					icpuend = icpu-1;
				} else if (iblknew >= blksblk[icpu+1]) {
					icpubeg = icpu+1;
				};
			};
			icpu = icpubeg;
			icpuprev = icpu;
		};
		_blk2cpu[iblk] = icpu;
	};

// Free work arrays

	delete [] blkstot;
	delete [] blkstotblk;

	for (i=0;i<=nlev;i++) {
		delete [] pporderlev[i];
		delete [] ppblkslev[i];
		delete [] pporderlevblk[i];
		delete [] ppblkslevsups[i];
		delete [] ppsupslev[i];
	};

	delete [] ppblkslev[nlev+1];
	delete [] ppblkslevsups[nlev+1];

	delete [] pporderlev;
	delete [] pporderlevblk;
	delete [] nblkslev;
	delete [] nblkslevsups;
	delete [] ppblkslev;
	delete [] ppblkslevsups;
	delete [] nsupslev;
	delete [] ppsupslev;
	delete [] commlev;
	delete [] listcpu;
	delete [] blkstotext;
	delete [] blkstotblkext;
	delete [] blksblk;
	delete [] sprndso;
	delete [] blk2cpuloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute partitioning of the matrix into many parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::SplitBlockMatrixRecursive (int _nparts, // Compute partitioning of the matrix into many parts taking into account blocks partitioning
														int _nblks, int *_blks, 
														int *_blk2cpu) const {

	const char *funcname = "SplitBlockMatrixRecursive_01";

// Check that the number of parts is a power of 2

	int nparts_pwr2 = 1;

	while (nparts_pwr2*2 <= _nparts) nparts_pwr2 *= 2;

	if (nparts_pwr2 != _nparts) {
		cout << " throw CSMatrix::SplitBlockMatrixRecursive: the number of nodes is not a power of 2 " << endl;
		throw " CSMatrix::SplitBlockMatrixRecursive: the number of nodes is not a power of 2 ";
	};

// Allocate work memory

	int *listnd;
	int *imasknd;
	int *irow2ind;
	int *listblk_old;
	int *listblk_new;
	int *hblks_old;
	int *hblks_new;
	int *blks_new;
	int *blk2cpu;
	int *nd2blk;

	listnd = new int [n];
	if(!listnd) MemoryFail(funcname);
	imasknd = new int [n];
	if(!imasknd) MemoryFail(funcname);
	irow2ind = new int [n];
	if(!irow2ind) MemoryFail(funcname);
	listblk_old = new int [_nblks];
	if(!listblk_old) MemoryFail(funcname);
	listblk_new = new int [_nblks];
	if(!listblk_new) MemoryFail(funcname);
	hblks_old = new int [_nparts*2+1];
	if(!hblks_old) MemoryFail(funcname);
	hblks_new = new int [_nparts*2+1];
	if(!hblks_new) MemoryFail(funcname);
	blks_new = new int [_nblks*2+1];
	if(!blks_new) MemoryFail(funcname);
	blk2cpu = new int [_nblks];
	if(!blk2cpu) MemoryFail(funcname);
	nd2blk = new int [n];
	if(!nd2blk) MemoryFail(funcname);

	int i, j;

	int icycle = -1;

	for (i=0;i<n;i++) imasknd[i] = icycle;

// Implement main cycle over the levels

	nparts_pwr2 = 1;

	hblks_old[0] = 0;
	hblks_old[1] = _nblks;

	for (i=0;i<_nblks;i++) listblk_old[i] = i;

	int nparts_pwr2_new, nlistblk, ibs, iblk, jblk, n1blk, n2blk, nlistnd, iblkold;

	CSMatrix mtrsub;
	CSMatrix mtrdummy;

	hblks_new[0] = 0;

	while (nparts_pwr2 != _nparts) {
		nparts_pwr2_new = nparts_pwr2 * 2;

		for (iblk=0;iblk<nparts_pwr2;iblk++) {

// Get the submatrix corresponding to the sublist

			nlistblk = hblks_old[iblk+1]-hblks_old[iblk];
			ibs = hblks_old[iblk];

			nlistnd = 0;
			blks_new[0] = 0;

			for (jblk=0;jblk<nlistblk;jblk++) {
				iblkold = listblk_old[ibs+jblk];
				for (j=_blks[iblkold];j<_blks[iblkold+1];j++) {
					listnd[nlistnd] = j;
					nd2blk[nlistnd] = jblk;
					nlistnd++;
				};
				blks_new[jblk+1] = blks_new[jblk] + _blks[iblkold+1]-_blks[iblkold];
			};

			SubmatrixList (nlistnd, listnd,
								icycle, imasknd, irow2ind,
								mtrsub);

// Compute new partitioning
/*
			if (true) {

				int *pialoc = mtrsub.GetIa ();
				int *pjaloc = mtrsub.GetJa ();

				int icheck = 1;

				for (int ilist=0;ilist<nlistnd;ilist++) {
					int ielem = ilist;
					int icheckloc = -1;
					for (int j=pialoc[ielem];j<pialoc[ielem+1];j++) {
						int jj = pjaloc[j];
						if (jj != ielem) {
							icheckloc = 1;
						};
					};
					if (icheckloc == -1) {
						cout << " Wrong matrix elem = " << ielem << endl;
						icheck = -1;
					};
				};

			};
*/
			mtrsub.SplitBlockMatrix (nlistblk, blks_new,
												blk2cpu);
/*
			if (true) {

				int *pialoc = mtrsub.GetIa ();
				int *pjaloc = mtrsub.GetJa ();

				for (int ipart=0;ipart<2;ipart++) {

					int icheck = 1;

					for (int ilist=0;ilist<nlistnd;ilist++) {
						int ielem = ilist;
						int iblk = nd2blk[ilist];
						if (blk2cpu[iblk] == ipart) {
							int icheckloc = -1;
							for (int j=pialoc[ielem];j<pialoc[ielem+1];j++) {
								int jj = pjaloc[j];
								if (jj != ielem) {
									jblk = nd2blk[jj];
									if (blk2cpu[jblk] == ipart) {
										icheckloc = 1;
									};
								};
							};
							if (icheckloc == -1) {
								cout << " After local split: Wrong part = " << ipart << " Wrong elem = " << ielem << endl;
								icheck = -1;
							};
						};
					};
				};
			};
*/
			mtrsub = mtrdummy;

// Store computed local partitioning

			n1blk = 0;
			n2blk = 0;

			for (i=0;i<nlistblk;i++) {
				if (blk2cpu[i] == 0) n1blk++;
				if (blk2cpu[i] == 1) n2blk++;
			};

			hblks_new[iblk*2+1] = hblks_new[iblk*2]+n1blk;
			hblks_new[iblk*2+2] = hblks_new[iblk*2+1]+n2blk;

			n2blk = n1blk;
			n1blk = 0;

			for (i=0;i<nlistblk;i++) {
				if (blk2cpu[i] == 0) {
					listblk_new[ibs+n1blk] = listblk_old[ibs+i];
					n1blk++;
				};
				if (blk2cpu[i] == 1) {
					listblk_new[ibs+n2blk] = listblk_old[ibs+i];
					n2blk++;
				};
			};

		};

// Prepare new partioning

		for (i=0;i<nparts_pwr2_new+1;i++) hblks_old[i] = hblks_new[i];
		for (i=0;i<_nblks;i++) listblk_old[i] = listblk_new[i];

		nparts_pwr2 = nparts_pwr2_new;

	};

// Return the result

	int jj;

	for (i=0;i<_nparts;i++) {
		for (j=hblks_old[i];j<hblks_old[i+1];j++) {
			jj = listblk_old[j];
			_blk2cpu[jj] = i;
		};
	};

// Free work arrays

	delete [] listnd;
	delete [] imasknd;
	delete [] irow2ind;
	delete [] listblk_old;
	delete [] listblk_new;
	delete [] hblks_old;
	delete [] hblks_new;
	delete [] blks_new;
	delete [] blk2cpu;
	delete [] nd2blk;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::SplitBlockMatrix (CMPIComm &_comm, // Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
											int _nblks, int *_blks,
											int *&_orderblk, int *_blkspartblk,
											int *&_ordernd, int *_blkspartnd,
											int _nproc1, int _igroup, int _myid1,
											int &_nblkslev, int *&_blkslev,
											int &_nblkslevblk, int *&_blkslevblk,
											int &_nblks1, int *&_blks1) const {

	const char *funcname = "SplitBlockMatrix_00";

// Collect sizes statistics for the whole set of submatrices

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	int *blkscpu;
	int *blk2cpu;

	blkscpu = new int [nproc+1];
	if (!blkscpu) MemoryFail (funcname);
	blk2cpu = new int [nproc];
	if (!blk2cpu) MemoryFail (funcname);

	int i;

	for (i=0;i<nproc+1;i++) blkscpu[i] = 0;
	for (i=0;i<nproc;i++) blk2cpu[i] = i;

	blkscpu[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blkscpu, blkscpu);

	for (i=0;i<nproc;i++) blkscpu[i+1] = blkscpu[i]+blkscpu[i+1];

// Eliminate diagonal

	CSMatrix asfltdini = *this;

	int j, jj, irow;

	int nzjaloc = 0;

	asfltdini.ia[0] = 0;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (irow != jj) {
				asfltdini.ja[nzjaloc] = jj;
				nzjaloc++;
			};
		};
		asfltdini.ia[i+1] = nzjaloc;
	};

	asfltdini.nzja = nzjaloc;
	asfltdini.nsupc = asfltdini.n;
	asfltdini.nsupr = asfltdini.n;

// Allocate the memory

	int *partition;
	int *partition_blk;
	int *nd2blk;

	partition = new int [nlist];
	if (!partition) MemoryFail (funcname);
	partition_blk = new int [_nblks];
	if (!partition_blk) MemoryFail (funcname);
	nd2blk = new int [nlist];
	if (!nd2blk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			nd2blk[j] = i;
		};
	};

// Connectivity test on entry
/*
	if (false) {

		int *partition_nd_new;

		partition_nd_new = new int [nlist];
		if (!partition_nd_new) MemoryFail (funcname);

		for (i=0;i<_nblks;i++) {
			for (j=_blks[i];j<_blks[i+1];j++) {
				partition_nd_new[j] = 0;
			};
		};

		int nsets;
		int *iasets, *jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 0,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || initial data connectivity test = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partition_nd_new;

	};
*/
// Call ParMETIS

	int options[10];
	options[0] = 0;

	int numflag = 0;

	int ncon = 1;
	int nparts = 2;

	float *tpwgts, *ubvec;

	tpwgts = new float [ncon*nparts];
	if (!tpwgts) MemoryFail (funcname);
	ubvec = new float [ncon];
	if (!ubvec) MemoryFail (funcname);

	for (i=0;i<nparts;i++) tpwgts[i] = (float)1.0 / ((float)nparts);

	ubvec[0] = (float) 1.05;

	int edgecut;
	int iopt = 0;
/*
	if (true) {

		int *blk2cpuloc;

		blk2cpuloc = new int [nproc];
		if (!blk2cpuloc) MemoryFail (funcname);

		for (i=0;i<nproc;i++) blk2cpuloc[i] = i;

		CSMatrix::CheckSymmMtr (_comm,
										nproc, blkscpu, blk2cpuloc,
										asfltdini.ia, asfltdini.ja,
										false, NULL);

		delete [] blk2cpuloc;

	};
*/
//	if (myid == 0) {
//		cout << " CSMatrix::SplitBlockMatrix: before PartKway " << endl;
//	};

	if (false) {
		CParMETIS::PartKway (_comm, blkscpu, asfltdini.ia, asfltdini.ja, 
									NULL, NULL, &iopt,
									&numflag, &ncon, &nparts,
									tpwgts, ubvec,
									options, &edgecut, partition);
	} else {

		int *listcpu;
		int *blkscpunew;

		listcpu = new int [nproc];
		if (!listcpu) MemoryFail (funcname);
		blkscpunew = new int [nproc+1];
		if (!blkscpunew) MemoryFail (funcname);

		int nlistcpu = 0;

		blkscpunew[0] = 0;

		for (i=0;i<nproc;i++) {
			if (blkscpu[i+1] > blkscpu[i]) {
				listcpu[nlistcpu] = i;
				blkscpunew[nlistcpu+1] = blkscpu[i+1];
				nlistcpu++;
			};
		};

		if (nlistcpu == nproc) {
			CParMETIS::PartKway (_comm, blkscpu, asfltdini.ia, asfltdini.ja, 
										NULL, NULL, &iopt,
										&numflag, &ncon, &nparts,
										tpwgts, ubvec,
										options, &edgecut, partition);
		} else {

			if (nlistcpu > 0) {

				CMPIComm commnew;

				commnew = _comm.CreateComm (nlistcpu, listcpu);

				if (blkscpu[myid+1] > blkscpu[myid]) {
					if (nlistcpu > 1) {
						CParMETIS::PartKway (commnew, blkscpunew, asfltdini.ia, asfltdini.ja, 
													NULL, NULL, &iopt,
													&numflag, &ncon, &nparts,
													tpwgts, ubvec,
													options, &edgecut, partition);
					} else if (nlistcpu == 1) {

						int n1tot = blkscpunew[nlistcpu];

						if (n1tot > 0) {
							METIS_WPartGraphRecursive (&n1tot, asfltdini.ia, asfltdini.ja, 
																NULL, NULL, &iopt,
																&numflag, &nparts, tpwgts,
																options, &edgecut, partition);
						};

					};
				};

			};

		};

		delete [] listcpu;
		delete [] blkscpunew;

	};

//	if (myid == 0) {
//		cout << " CSMatrix::SplitBlockMatrix: after PartKway " << endl;
//	};
/*
	int *wvertloc;

	wvertloc = new int [nlist];
	if (!wvertloc) MemoryFail (funcname);

	for (i=0;i<nlist;i++) wvertloc[i] = 0;

	ModifyDisjointParts (_comm,
								blkscpu,
								asfltdini.ia, asfltdini.ja, wvertloc, partition);

	delete [] wvertloc;
*/
	CSMatrix mtrdummy;

// Check connectivity of the resulting point partitioning
/*
	if (false) {

		int *partition_nd_new;

		partition_nd_new = new int [nlist];
		if (!partition_nd_new) MemoryFail (funcname);

		for (i=0;i<nlist;i++) {
			partition_nd_new[i] = partition[i];
		};

		int nsets;
		int *iasets, *jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 0,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || point connectivity test, index 0 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 1,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || point connectivity test, index 1 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partition_nd_new;

	};
*/
//	asfltdini = mtrdummy;

	int iblk, ilist;

	for (iblk=0;iblk<_nblks;iblk++) {
		ilist = _blks[iblk];
		if (ilist < _blks[iblk+1]) {
			partition_blk[iblk] = partition[ilist];
		};
	};

// Find local point two-sided boundary

	int nlistbnd;
	int *listbnd;

	ComputeTwoSidedBoundary (_comm, blkscpu, partition,
										nlistbnd, listbnd);

// Check connectivity after splitting the boundary
/*
	if (false) {

		int *partition_nd_new;

		partition_nd_new = new int [nlist];
		if (!partition_nd_new) MemoryFail (funcname);

		for (i=0;i<nlist;i++) {
			partition_nd_new[i] = partition[i];
		};

		int ishift = blkscpu[myid];

		for (i=0;i<nlistbnd;i++) {
			ilist = listbnd[i]-ishift;
			partition_nd_new[ilist] = -1;
		};

		int nsets;
		int *iasets, *jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 0,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || point connectivity test after splitting, index 0 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 1,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || point connectivity test after splitting, index 1 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partition_nd_new;

	};
*/
// Main extension cycle

	int icycle;
	int ishift = blkscpu[myid];

	int ncycle_ext = 1;

	for (icycle=0;icycle<ncycle_ext;icycle++) {

// Mark some blocks as the boundary ones

		for (i=0;i<nlistbnd;i++) {
			ilist = listbnd[i]-ishift;
			iblk = nd2blk[ilist];
			partition_blk[iblk] = -1;
		};

		delete [] listbnd;

// Mark as boundary the blocks with different nodes types inside

		int ipart, jpart;

		for (iblk=0;iblk<_nblks;iblk++) {
			ipart = partition_blk[iblk];
			for (j=_blks[iblk];j<_blks[iblk+1];j++) {
				jpart = partition[j];
				if (jpart != ipart) {
					partition_blk[iblk] = -1;
					break;
				};
			};
		};

// Create the list including the whole blocks

		nlistbnd = 0;

		for (iblk=0;iblk<_nblks;iblk++) {
			if (partition_blk[iblk] == -1) {
				nlistbnd += (_blks[iblk+1]-_blks[iblk]);
			};
		};

		listbnd = new int [nlistbnd];
		if (!listbnd) MemoryFail (funcname);

		nlistbnd = 0;

		for (iblk=0;iblk<_nblks;iblk++) {
			if (partition_blk[iblk] == -1) {
				for (j=_blks[iblk];j<_blks[iblk+1];j++) {
					listbnd[nlistbnd] = j+ishift;
					nlistbnd++;
				};
			};
		};

// Extend the boundary nodes in parallel

		ExtendNodes (_comm, 
							blkscpu, nlistbnd, listbnd);

// Mark again some blocks as the boundary ones

		for (i=0;i<nlistbnd;i++) {
			ilist = listbnd[i]-ishift;
			iblk = nd2blk[ilist];
			partition_blk[iblk] = -1;
		};

	};

// Check connectivity of the resulting partitioning
/*
	if (false) {

		int *partition_nd_new;

		partition_nd_new = new int [nlist];
		if (!partition_nd_new) MemoryFail (funcname);

		for (i=0;i<_nblks;i++) {
			for (j=_blks[i];j<_blks[i+1];j++) {
				partition_nd_new[j] = partition_blk[i];
			};
		};

		int nsets;
		int *iasets, *jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 0,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || connectivity test after extend, index 0 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 1,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || connectivity test after extend, index 1 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partition_nd_new;

	};
*/

// Compute block sparsity

	CSMatrix ablk;

	ablk = ComputeBlockSparsityLocal (_comm, _nblks, _blks);

// Compute partitioning in blocks

	int *gblks;

	gblks = new int [nproc+1];
	if (!gblks) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) gblks[i] = 0;

	gblks[myid+1] = _nblks;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, gblks, gblks);

	for (i=0;i<nproc;i++) gblks[i+1] = gblks[i]+gblks[i+1];

// Perform computation of the sets and mark the blocks

	int nsets0;
	int *iasets0, *jasets0;

	ComputeDisjointSets (_comm,
								gblks, 0,
								ablk.ia, ablk.ja, partition_blk,
								nsets0, iasets0, jasets0);

	int nsets1;
	int *iasets1, *jasets1;

	ComputeDisjointSets (_comm,
								gblks, 1,
								ablk.ia, ablk.ja, partition_blk,
								nsets1, iasets1, jasets1);

// Reassign block numbers

	int *partition_blk_sets;

	partition_blk_sets = new int [_nblks];
	if (!partition_blk_sets) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) partition_blk_sets[i] = partition_blk[i];

	for (i=0;i<nsets0;i++) {
		for (j=iasets0[i];j<iasets0[i+1];j++) {
			jj = jasets0[j];
			partition_blk_sets[jj] = i;
		};
	};
	for (i=0;i<nsets1;i++) {
		for (j=iasets1[i];j<iasets1[i+1];j++) {
			jj = jasets1[j];
			partition_blk_sets[jj] = nsets0+i;
		};
	};

	delete [] gblks;
	delete [] iasets0;
	delete [] jasets0;
	delete [] iasets1;
	delete [] jasets1;

	int nsets12 = nsets0+nsets1;

// Count the total number of blocks for each part

	int *nblks_sets;

	nblks_sets = new int [nsets12];
	if (!nblks_sets) MemoryFail (funcname);

	for (i=0;i<nsets12;i++) nblks_sets[i] = 0;

	int nblksbnd=0;

	int iset;

	for (i=0;i<_nblks;i++) {
		iset = partition_blk_sets[i];
		if (iset >= 0) {
			nblks_sets[iset]++;
		} else {
			nblksbnd++;
		};
	};

	int nblksbndloc = nblksbnd;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsets12, nblks_sets, nblks_sets);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &nblksbnd, &nblksbnd);

	int *gblkscpu;

	gblkscpu = new int [nproc+1];
	if (!gblkscpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) gblkscpu[i] = 0;

	gblkscpu[myid+1] = nblksbndloc;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, gblkscpu, gblkscpu);

	gblkscpu[1] += nsets12;

	for (i=0;i<nproc;i++) gblkscpu[i+1] = gblkscpu[i]+gblkscpu[i+1];

// Perform computation of the sparsity pattern and weights data in terms of blocks

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkWei_",myid,".dat");
//	ofstream fout (strbuff);

	int *wvertices, *wedges;

	CSMatrix asfltd;

//	asfltd = ComputeBlockSparsityWeights (_comm, _nblks, _blks, 
//														blkscpu, gblkscpu, partition_blk,
//														wvertices, wedges);
	asfltd = ComputeBlockSparsityWeightsSets (_comm, _nblks, _blks, nsets12,
															blkscpu, gblkscpu, partition_blk_sets,
															wvertices, wedges);

//	fout << " Asp = " << asfltd << endl;
//	OutArr (fout," Wvert = ",asfltd.GetNlist(),wvertices);
//	OutArr (fout," Wedges = ",asfltd.GetNzja(),wedges);
//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","Afiltr_",myid,".ps");
//	asfltd.A2Ps (strbuff,0,&i);

// Perform decomposition with weights

	int *partition_blk_new;

	partition_blk_new = new int [nblksbndloc+nsets12];
	if (!partition_blk_new) MemoryFail (funcname);

	iopt = 3;
/*
	if (false) {

		int *blk2cpuloc;

		blk2cpuloc = new int [nproc];
		if (!blk2cpuloc) MemoryFail (funcname);

		for (i=0;i<nproc;i++) blk2cpuloc[i] = i;

		CSMatrix::CheckSymmMtr (_comm,
										nproc, gblkscpu, blk2cpuloc,
										asfltd.ia, asfltd.ja,
										true, wedges);

		delete [] blk2cpuloc;

	};
*/
//	if (myid == 0) {
//		cout << " CSMatrix::SplitBlockMatrix weights: before PartKway " << endl;
//		OutArr (cout, " Gblkscpu =",nproc+1,gblkscpu);
//	};

	if (false) {

		CParMETIS::PartKway (_comm, gblkscpu, asfltd.ia, asfltd.ja, 
									wvertices, wedges, &iopt,
									&numflag, &ncon, &nparts,
									tpwgts, ubvec,
									options, &edgecut, partition_blk_new);

	} else {

		int nploc = 0;
		int *ibegarrloc;
		int *iendarrloc;

		ibegarrloc = new int [nproc];
		if (!ibegarrloc) MemoryFail (funcname);
		iendarrloc = new int [nproc];
		if (!iendarrloc) MemoryFail (funcname);

		if (myid == 0) {
			for (i=0;i<nproc;i++) {
				if (gblkscpu[i] < gblkscpu[i+1]) {
					ibegarrloc[nploc] = gblkscpu[i];
					iendarrloc[nploc] = gblkscpu[i+1]-1;
					nploc++;
				};
			};
		};

		CSMatrix a1cpu;
		a1cpu = asfltd.GetSubmatrixViaIntervals (_comm,
																nploc, ibegarrloc, iendarrloc);

		delete [] ibegarrloc;
		delete [] iendarrloc;

		int nnntot = gblkscpu[nproc];

		int *partitionblk_1cpu;
		int *wvertices_tot;

		partitionblk_1cpu = new int [nnntot];
		if (!partitionblk_1cpu) MemoryFail (funcname);
		wvertices_tot = new int [nnntot];
		if (!wvertices_tot) MemoryFail (funcname);

		for (i=0;i<nnntot;i++) partitionblk_1cpu[i] = 0;
		for (i=0;i<nnntot;i++) wvertices_tot[i] = 0;

		for (i=gblkscpu[myid];i<gblkscpu[myid+1];i++) {
			wvertices_tot[i] = wvertices[i-gblkscpu[myid]];
		};

		CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nnntot, wvertices_tot, wvertices_tot);

		int *ianzcpu;

		ianzcpu = new int [nproc+1];
		if (!ianzcpu) MemoryFail (funcname);

		for (i=0;i<=nproc;i++) ianzcpu[i] = 0;

		ianzcpu[myid+1] = asfltd.nzja;

		CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, ianzcpu, ianzcpu);

		for (i=0;i<nproc;i++) ianzcpu[i+1] = ianzcpu[i]+ianzcpu[i+1];

		int nnztot = ianzcpu[nproc];

		int *wedges_tot;

		wedges_tot = new int [nnztot];
		if (!wedges_tot) MemoryFail (funcname);

		for (i=0;i<nnztot;i++) wedges_tot[i] = 0;

		for (i=ianzcpu[myid];i<ianzcpu[myid+1];i++) {
			wedges_tot[i] = wedges[i-ianzcpu[myid]];
		};

		CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nnztot, wedges_tot, wedges_tot);

		if (myid == 0) {

			if (nnntot > 0) {
				METIS_WPartGraphRecursive (&nnntot, a1cpu.ia, a1cpu.ja, wvertices_tot, wedges_tot, &iopt,
													&numflag, &nparts, tpwgts,
													options, &edgecut, partitionblk_1cpu);
			};

// Compute local disjoint sets and locally modify the partitioning

			int nsetsloc0, nsetsloc1;
			int *iasetsloc0, *jasetsloc0, *iasetsloc1, *jasetsloc1;

			CSMatrix::ComputeDisjointSets (0, nnntot-1, 0,
														a1cpu.ia, a1cpu.ja, partitionblk_1cpu,
														nsetsloc0, iasetsloc0, jasetsloc0);
			CSMatrix::ComputeDisjointSets (0, nnntot-1, 1,
														a1cpu.ia, a1cpu.ja, partitionblk_1cpu,
														nsetsloc1, iasetsloc1, jasetsloc1);

			int *wsets0, *wsets1;

			wsets0 = new int [nsetsloc0];
			if (!wsets0) MemoryFail (funcname);
			wsets1 = new int [nsetsloc1];
			if (!wsets1) MemoryFail (funcname);

			int iaux;

			for (i=0;i<nsetsloc0;i++) {
				iaux = 0;
				for (j=iasetsloc0[i];j<iasetsloc0[i+1];j++) {
					jj = jasetsloc0[j];
					iaux += wvertices_tot[jj];
				};
				wsets0[i] = iaux;
			};
			for (i=0;i<nsetsloc1;i++) {
				iaux = 0;
				for (j=iasetsloc1[i];j<iasetsloc1[i+1];j++) {
					jj = jasetsloc1[j];
					iaux += wvertices_tot[jj];
				};
				wsets1[i] = iaux;
			};

			int wtot = 0;

			for (i=0;i<nsetsloc0;i++) wtot += wsets0[i];
			for (i=0;i<nsetsloc1;i++) wtot += wsets1[i];

			int wthresh = wtot / 100;

			for (i=0;i<nsetsloc0;i++) {
				if (wsets0[i] < wthresh) {
					for (j=iasetsloc0[i];j<iasetsloc0[i+1];j++) {
						jj = jasetsloc0[j];
						partitionblk_1cpu[jj] = 1;
					};
				};
			};
			for (i=0;i<nsetsloc1;i++) {
				if (wsets1[i] < wthresh) {
					for (j=iasetsloc1[i];j<iasetsloc1[i+1];j++) {
						jj = jasetsloc1[j];
						partitionblk_1cpu[jj] = 0;
					};
				};
			};

			delete [] iasetsloc0;
			delete [] jasetsloc0;
			delete [] iasetsloc1;
			delete [] jasetsloc1;
			delete [] wsets0;
			delete [] wsets1;

		};

		CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nnntot, partitionblk_1cpu, partitionblk_1cpu);

		for (i=gblkscpu[myid];i<gblkscpu[myid+1];i++) {
			partition_blk_new[i-gblkscpu[myid]] = partitionblk_1cpu[i];
		};

		delete [] partitionblk_1cpu;
		delete [] wvertices_tot;
		delete [] ianzcpu;
		delete [] wedges_tot;

	};

//	if (myid == 0) {
//		cout << " CSMatrix::SplitBlockMatrix weights: after PartKway " << endl;
//	};

// Transform the result

	int *ipartarr;

	ipartarr = new int [nsets12];
	if (!ipartarr) MemoryFail (funcname);

	for (i=0;i<nsets12;i++) ipartarr[i] = 0;

	if (myid == 0) {
		for (i=0;i<nsets12;i++) ipartarr[i] = partition_blk_new[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsets12, ipartarr, ipartarr);

	int nz = 0;

	if (myid == 0) nz = nsets12;

	int *blk2cpuloc;

	blk2cpuloc = new int [_nblks];
	if (!blk2cpuloc) MemoryFail (funcname);

	int ipart;

	for (i=0;i<_nblks;i++) {
		if (partition_blk_sets[i] >= 0) {
			ipart = partition_blk_sets[i];
			blk2cpuloc[i] = ipartarr[ipart];
		} else {
			blk2cpuloc[i] = partition_blk_new[nz];
			nz++;
		};
	};

// Check connectivity of the resulting partitioning
/*
	if (false) {

		int *partition_nd_new;

		partition_nd_new = new int [nlist];
		if (!partition_nd_new) MemoryFail (funcname);

		for (i=0;i<_nblks;i++) {
			for (j=_blks[i];j<_blks[i+1];j++) {
				partition_nd_new[j] = blk2cpuloc[i];
			};
		};

		int nsets;
		int *iasets, *jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 0,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || connectivity test, index 0 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		ComputeDisjointSets (_comm,
									blkscpu, 1,
									asfltdini.ia, asfltdini.ja, partition_nd_new,
									nsets, iasets, jasets);

		cout << " || connectivity test, index 1 = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partition_nd_new;

	};
*/
// Create final output data

	int *nicpuarr;

	nicpuarr = new int [nproc*2+1];
	if (!nicpuarr) MemoryFail (funcname);

	for (i=0;i<=nproc*2;i++) nicpuarr[i] = 0;

	int ni;

	for (i=0;i<_nblks;i++) {
		ni = _blks[i+1]-_blks[i];
		if (blk2cpuloc[i] == 0) {
			nicpuarr[myid+1] += ni;
		} else {
			nicpuarr[nproc+myid+1] += ni;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 2*nproc+1, nicpuarr, nicpuarr);

	for (i=0;i<nproc*2;i++) nicpuarr[i+1] = nicpuarr[i]+nicpuarr[i+1];

	_blkspartnd[0] = 0;
	_blkspartnd[1] = nicpuarr[nproc];
	_blkspartnd[2] = nicpuarr[2*nproc];

	int ibeg0 = nicpuarr[myid];
	int ibeg1 = nicpuarr[nproc+myid];

	_ordernd = new int [nlist];
	if (!_ordernd) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) {
		if (blk2cpuloc[i] == 0) {
			for (j=_blks[i];j<_blks[i+1];j++) {
				_ordernd[j] = ibeg0;
				ibeg0++;
			};
		} else {
			for (j=_blks[i];j<_blks[i+1];j++) {
				_ordernd[j] = ibeg1;
				ibeg1++;
			};
		};
	};

// Create final output block data

	for (i=0;i<=nproc*2;i++) nicpuarr[i] = 0;

	for (i=0;i<_nblks;i++) {
		ni = _blks[i+1]-_blks[i];
		if (blk2cpuloc[i] == 0) {
			nicpuarr[myid+1]++;
		} else {
			nicpuarr[nproc+myid+1]++;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 2*nproc+1, nicpuarr, nicpuarr);

	for (i=0;i<nproc*2;i++) nicpuarr[i+1] = nicpuarr[i]+nicpuarr[i+1];

	_blkspartblk[0] = 0;
	_blkspartblk[1] = nicpuarr[nproc];
	_blkspartblk[2] = nicpuarr[2*nproc];

	ibeg0 = nicpuarr[myid];
	ibeg1 = nicpuarr[nproc+myid];

	_orderblk = new int [_nblks];
	if (!_orderblk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) {
		if (blk2cpuloc[i] == 0) {
			_orderblk[i] = ibeg0;
			ibeg0++;
		} else {
			_orderblk[i] = ibeg1;
			ibeg1++;
		};
	};

// Allocate partioning data

	_nblkslev = _nproc1;
	_nblkslevblk = _nproc1;

	_blkslev = new int [_nblkslev+1];
	if (!_blkslev) MemoryFail (funcname);

	_blkslevblk = new int [_nblkslev+1];
	if (!_blkslevblk) MemoryFail (funcname);

// Compute optimal partitioning of the remaining data over the set of processors

	int ldloc = _nproc1+1;
	int isize = 4*ldloc+3;

	int *ilevel2cpu;

	ilevel2cpu = new int [isize];
	if (!ilevel2cpu) MemoryFail (funcname);

	for (i=0;i<isize;i++) ilevel2cpu[i] = -1;

	int *blk2cpusend;

	blk2cpusend = new int [_nblks];
	if (!blk2cpusend) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) blk2cpusend[i] = -1;

	int ni1tot = _blkspartnd[1]-_blkspartnd[0];
	int ni2tot = _blkspartnd[2]-_blkspartnd[1];

	int ni1 = ni1tot / _nproc1;

	ilevel2cpu[0] = 0;
	ilevel2cpu[_nproc1] = nproc-1;

	ilevel2cpu[ldloc+0] = 0;
	ilevel2cpu[ldloc+_nproc1] = nicpuarr[nproc]-nicpuarr[nproc-1];

	ilevel2cpu[2*ldloc+0] = 0;
	for (i=1;i<=_nproc1;i++) ilevel2cpu[2*ldloc+i] = _blkspartnd[1]-_blkspartnd[0];
	ilevel2cpu[3*ldloc+0] = 0;
	for (i=1;i<=_nproc1;i++) ilevel2cpu[3*ldloc+i] = _blkspartblk[1]-_blkspartblk[0];

	ilevel2cpu[isize-3] = 0;
	ilevel2cpu[isize-2] = 0;
	ilevel2cpu[isize-1] = 0;

	int icpu, nicurr, niblkcurr, niloc;

	CMPIStatus stat;

	for (icpu=0;icpu<nproc;icpu++) {

// Scan current data and prepare the sends

		if (myid == icpu) {
			ipart = ilevel2cpu[isize-3];
			niblkcurr = ilevel2cpu[isize-2];
			nicurr = ilevel2cpu[isize-1];
			iblk = -1;
			while (iblk<_nblks-1) {
				iblk++;
				if (blk2cpuloc[iblk] == 0) {
					niloc = _blks[iblk+1]-_blks[iblk];
					if (nicurr+niloc >= ni1 && ipart < _nproc1-1) {
						ilevel2cpu[ipart+1] = myid;
						ilevel2cpu[ldloc+ipart+1] = iblk+1;
						ilevel2cpu[2*ldloc+ipart+1] = ilevel2cpu[2*ldloc+ipart]+nicurr+niloc;
						ilevel2cpu[3*ldloc+ipart+1] = ilevel2cpu[3*ldloc+ipart]+niblkcurr+1;
						blk2cpusend[iblk] = ipart;
						nicurr = 0;
						niblkcurr = 0;
						ipart++;
					} else {
						nicurr += niloc;
						niblkcurr++;
						blk2cpusend[iblk] = ipart;
					};
				};
			};
			ilevel2cpu[isize-3] = ipart;
			ilevel2cpu[isize-2] = niblkcurr;
			ilevel2cpu[isize-1] = nicurr;
		};

// Send the data

		if (icpu < nproc-1 && myid == icpu) {
			CMPIExchange::Send (_comm, icpu+1, icpu+1,
										isize*sizeof(int), (char *)ilevel2cpu);
		};

// Receive the data

		if (myid == icpu+1) {
			CMPIExchange::Recv (_comm, icpu, icpu+1,
										isize*sizeof(int), (char *)ilevel2cpu, stat);
		};

	};

	if (myid < nproc-1) {
		for (i=0;i<isize;i++) ilevel2cpu[i] = 0;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, isize, ilevel2cpu, ilevel2cpu);

	if (_igroup == 0) {
		for (i=0;i<=_nproc1;i++) _blkslev[i] = ilevel2cpu[2*ldloc+i];
		for (i=0;i<=_nproc1;i++) _blkslevblk[i] = ilevel2cpu[3*ldloc+i];
	};

	int ni2 = ni2tot / _nproc1;

	for (i=0;i<isize;i++) ilevel2cpu[i] = -1;

	ilevel2cpu[0] = 0;
	ilevel2cpu[_nproc1] = nproc-1;

	ilevel2cpu[ldloc+0] = 0;
	ilevel2cpu[ldloc+_nproc1] = nicpuarr[nproc]-nicpuarr[nproc-1];

	ilevel2cpu[2*ldloc+0] = 0;
	for (i=1;i<=_nproc1;i++) ilevel2cpu[2*ldloc+i] = _blkspartnd[2]-_blkspartnd[1];
	ilevel2cpu[3*ldloc+0] = 0;
	for (i=1;i<=_nproc1;i++) ilevel2cpu[3*ldloc+i] = _blkspartblk[2]-_blkspartblk[1];

	ilevel2cpu[isize-3] = 0;
	ilevel2cpu[isize-2] = 0;
	ilevel2cpu[isize-1] = 0;

	for (icpu=0;icpu<nproc;icpu++) {

// Scan current data and prepare the sends

		if (myid == icpu) {
			ipart = ilevel2cpu[isize-3];
			niblkcurr = ilevel2cpu[isize-2];
			nicurr = ilevel2cpu[isize-1];
			iblk = -1;
			while (iblk<_nblks-1) {
				iblk++;
				if (blk2cpuloc[iblk] == 1) {
					niloc = _blks[iblk+1]-_blks[iblk];
					if (nicurr+niloc >= ni2 && ipart < _nproc1-1) {
						ilevel2cpu[ipart+1] = myid;
						ilevel2cpu[ldloc+ipart+1] = iblk+1;
						ilevel2cpu[2*ldloc+ipart+1] = ilevel2cpu[2*ldloc+ipart]+nicurr+niloc;
						ilevel2cpu[3*ldloc+ipart+1] = ilevel2cpu[3*ldloc+ipart]+niblkcurr+1;
						blk2cpusend[iblk] = _nproc1+ipart;
						nicurr = 0;
						niblkcurr = 0;
						ipart++;
					} else {
						nicurr += niloc;
						niblkcurr++;
						blk2cpusend[iblk] = _nproc1+ipart;
					};
				};
			};
			ilevel2cpu[isize-3] = ipart;
			ilevel2cpu[isize-2] = niblkcurr;
			ilevel2cpu[isize-1] = nicurr;
		};

// Send the data

		if (icpu < nproc-1 && myid == icpu) {
			CMPIExchange::Send (_comm, icpu+1, icpu+1,
										isize*sizeof(int), (char *)ilevel2cpu);
		};

// Receive the data

		if (myid == icpu+1) {
			CMPIExchange::Recv (_comm, icpu, icpu+1,
										isize*sizeof(int), (char *)ilevel2cpu, stat);
		};

	};

	if (myid < nproc-1) {
		for (i=0;i<isize;i++) ilevel2cpu[i] = 0;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, isize, ilevel2cpu, ilevel2cpu);

	if (_igroup == 1) {
		for (i=0;i<=_nproc1;i++) _blkslev[i] = ilevel2cpu[2*ldloc+i];
		for (i=0;i<=_nproc1;i++) _blkslevblk[i] = ilevel2cpu[3*ldloc+i];
	};

// Prepare the send arrays

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iproc;

	for (i=0;i<_nblks;i++) {
		iproc = blk2cpusend[i];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int *isendarr;

	isendarr = new int [_nblks];
	if (!isendarr) MemoryFail (funcname);

	int iobj, k;

	for (i=0;i<_nblks;i++) {
		iproc = blk2cpusend[i];
		k = iptr[iproc];
		isendarr[k] = _blks[i+1]-_blks[i];
		iptr[iproc]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) (isendarr+iacpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::SplitBlockMatrix: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::SplitBlockMatrix: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Store the result

	CIntInt *iiarr;

	iiarr = new CIntInt [nproc];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nproc;i++) {
		iiarr[i].intvalue = CpuIDRecv[i];
		iiarr[i].int2value = i;
	};

	//qsort (iiarr, nproc, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + nproc);

	_nblks1 = 0;

	for (i=0;i<nproc;i++) {
		iobj = iiarr[i].int2value;
		ni = ObjSizeRecv[iobj] / (sizeof(int));
		_nblks1 += ni;
	};

	_blks1 = new int [_nblks1+1];
	if(!_blks1) MemoryFail(funcname);

	_blks1[0] = 0;
	_nblks1 = 0;

	int *piarr;

	for (i=0;i<nproc;i++) {
		iobj = iiarr[i].int2value;
		ni = ObjSizeRecv[iobj] / (sizeof(int));
		piarr = (int *) ObjRecv[iobj];
		for (j=0;j<ni;j++) {
			niloc = piarr[j];
			_blks1[_nblks1+1] = _blks1[_nblks1]+niloc;
			_nblks1++;
		};
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] blkscpu;
	delete [] blk2cpu;
	delete [] partition;
	delete [] partition_blk;
	delete [] nd2blk;
	delete [] tpwgts;
	delete [] ubvec;
	delete [] listbnd;
	delete [] gblkscpu;
	delete [] wvertices;
	delete [] wedges;
	delete [] partition_blk_new;
	delete [] blk2cpuloc;
	delete [] nicpuarr;
	delete [] ilevel2cpu;
	delete [] blk2cpusend;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;
	delete [] iiarr;
	delete [] ipartarr;
	delete [] nblks_sets;
	delete [] partition_blk_sets;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute partitioning of the matrix into two parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::SplitBlockMatrix (int _nblks, int *_blks, // Compute partitioning of the matrix into two parts taking into account blocks partitioning
											int *_blk2cpu) const {

	const char *funcname = "SplitBlockMatrix_01";

// Compute initial Nd ordering

//	A2Ps (10,"MtrIni.ps",0,&_nblks);

	int *order;

	order = new int [n];
	if (!order) MemoryFail (funcname);

	int ordtype = -1;

	OrderPrfMtr (ordtype, order);

// Reorder and find initial separators

	CSMatrix mtro;

	mtro = this->OrdMtr (order);

//	mtro.A2Ps (10,"MtrOrd.ps",0,&_nblks);

	int ind1, ind2;

	mtro.FindSeparator (ind1, ind2);

	if (ind1 < 0 || ind2 < 0) {
		ind1 = 0;
		ind2 = 0;
	};

// Compute inverse ordering and ind2blk

	int *iorder;
	int *ind2blk;

	iorder = new int [n];
	if(!iorder) MemoryFail(funcname);
	ind2blk = new int [n];
	if(!ind2blk) MemoryFail(funcname);

	int i, j;

	for (i=0;i<n;i++) iorder[order[i]] = i;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			ind2blk[j] = i;
		};
	};

// Mark all blocks according to the separators

	int *markblk;

	markblk = new int [_nblks];
	if(!markblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) markblk[i] = -2;

	int inew, iold, iblkold;

	for (inew=0;inew<ind1;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		if (markblk[iblkold] == -2) {
			markblk[iblkold] = 0;
		} else if (markblk[iblkold] == 1) {
			markblk[iblkold] = -1;
		};
	};

	for (inew=ind1;inew<ind2;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		if (markblk[iblkold] == -2) {
			markblk[iblkold] = 1;
		} else if (markblk[iblkold] == 0) {
			markblk[iblkold] = -1;
		};
	};

	for (inew=ind2;inew<n;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		markblk[iblkold] = -1;
	};

// Compute block sparsity of the initial matrix

	CSMatrix aspblk;
	CSMatrix aspblk_save;

	aspblk = this->BlockSparsity (_nblks, ind2blk);
	aspblk_save = aspblk;

// Extend the boundary according to the parameter

	int ncycleext = 1;

	int *listblk;

	listblk = new int [_nblks];
	if(!listblk) MemoryFail(funcname);

	int nlistblk = 0;

	for (i=0;i<_nblks;i++) {
		if (markblk[i] < 0) {
			listblk[nlistblk] = i;
			nlistblk++;
		};
	};

	int *piablk = aspblk.GetIa ();
	int *pjablk = aspblk.GetJa ();

	int icycleext, nlistblk1, iblk, jblk;

	for (icycleext=0;icycleext<ncycleext;icycleext++) {
		nlistblk1 = nlistblk;
		for (i=0;i<nlistblk;i++) {
			iblk = listblk[i];
			for (j=piablk[iblk];j<piablk[iblk+1];j++) {
				jblk = pjablk[j];
				if (markblk[jblk] >= 0) {
					listblk[nlistblk1] = jblk;
					nlistblk1++;
					markblk[jblk] = -1;
				};
			};
		};
		nlistblk = nlistblk1;
	};

	delete [] listblk;

// Search for non-connected small subparts and add them into the boundary

	int *iacomp;
	int *iacompnd;
	int *listcomp;
	int *maskcomp;

	iacomp = new int [_nblks+1];
	if(!iacomp) MemoryFail(funcname);
	iacompnd = new int [_nblks+1];
	if(!iacompnd) MemoryFail(funcname);
	listcomp = new int [_nblks];
	if(!listcomp) MemoryFail(funcname);
	maskcomp = new int [_nblks];
	if(!maskcomp) MemoryFail(funcname);

	int ipart;

	for (ipart=0;ipart<2;ipart++) {
		int ncomp = 0;
		iacomp[0] = 0;
		iacompnd[0] = 0;
		int nz = 0;
		int nznd = 0;
		for (i=0;i<_nblks;i++) maskcomp[i] = -1;
		while (true) {
			int icheck = -1;
			for (i=0;i<_nblks;i++) {
				if (maskcomp[i] == -1 && markblk[i] == ipart) {
					maskcomp[i] = ncomp;
					listcomp[nz] = i;
					nz++;
					nznd += _blks[i+1]-_blks[i];
					ncomp++;
					icheck = 1;
					break;
				};
			};
			if (icheck == -1) break;
			for (i=iacomp[ncomp-1];i<nz;i++) {
				iblk = listcomp[i];
				for (j=piablk[iblk];j<piablk[iblk+1];j++) {
					int jj = pjablk[j];
					if (maskcomp[jj] == -1 && markblk[jj] == ipart) {
						maskcomp[jj] = ncomp-1;
						listcomp[nz] = jj;
						nz++;
						nznd += _blks[jj+1]-_blks[jj];
					};
				};
			};
			iacomp[ncomp] = nz;
			iacompnd[ncomp] = nznd;
		};
		if (ncomp > 1) {
			int icompmax = 0;
			int nzmax = iacompnd[1]-iacompnd[0];
			for (i=1;i<ncomp;i++) {
				if (iacompnd[i+1]-iacompnd[i] > nzmax) {
					icompmax = i;
					nzmax = iacompnd[i+1]-iacompnd[i];
				};
			};
			for (i=0;i<ncomp;i++) {
				if (i != icompmax) {
					for (j=iacomp[i];j<iacomp[i+1];j++) {
						iblk = listcomp[j];
						markblk[iblk] = -1;
					};
				};
			};
		};
	};

	delete [] iacomp;
	delete [] iacompnd;
	delete [] listcomp;
	delete [] maskcomp;

// Compute block ordering

	int nblks1=0;
	int nblks2=0;
	int nblksbnd=0;

	for (i=0;i<_nblks;i++) {
		if (markblk[i] == 0) nblks1++;
		if (markblk[i] == 1) nblks2++;
		if (markblk[i] == -1) nblksbnd++;
	};

	int ibs1 = 0;
	int ibs2 = ibs1+nblks1;
	int ibsbnd = ibs2+nblks2;

	int *orderblk;
	int *iorderblk;

	orderblk = new int [_nblks];
	if(!orderblk) MemoryFail(funcname);
	iorderblk = new int [_nblks];
	if(!iorderblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) {
		if (markblk[i] == 0) {
			orderblk[i] = ibs1;
			ibs1++;
		};
		if (markblk[i] == 1) {
			orderblk[i] = ibs2;
			ibs2++;
		};
		if (markblk[i] == -1) {
			orderblk[i] = ibsbnd;
			ibsbnd++;
		};
	};

	for (i=0;i<_nblks;i++) iorderblk[orderblk[i]] = i;

// Recompute point ordering and new blocks partitioning

	int *blkso;
	int *orderloc;

	blkso = new int [_nblks+1];
	if(!blkso) MemoryFail(funcname);
	orderloc = new int [n];
	if(!orderloc) MemoryFail(funcname);

	blkso[0] = 0;

	for (i=0;i<_nblks;i++) {
		iold = iorderblk[i];
		blkso[i+1] = (_blks[iold+1]-_blks[iold]);
	};

	for (i=0;i<_nblks;i++) blkso[i+1] = blkso[i]+blkso[i+1];

	int ibs, ibsn, ni;

	for (i=0;i<_nblks;i++) {
		inew = orderblk[i];
		ibs = _blks[i];
		ibsn = blkso[inew];
		ni = _blks[i+1]-_blks[i];
		for (j=0;j<ni;j++) {
			orderloc[ibs+j] = ibsn+j;
		};
	};

// Reorder the matrix

	mtro = this->OrdMtr (orderloc);

// Recompute blocks partitioning

	int nblksloc = nblksbnd+2;

	int *blksobnd;

	blksobnd = new int [nblksloc+1];
	if(!blksobnd) MemoryFail(funcname);

	blksobnd[0] = 0;
	blksobnd[1] = blkso[nblks1];
	blksobnd[2] = blkso[nblks1+nblks2];

	for (i=0;i<nblksbnd;i++) blksobnd[i+3] = blksobnd[i+2]+(blkso[nblks1+nblks2+i+1]-blkso[nblks1+nblks2+i]);

// Compute weights data

	int *sp2blkloc;

	sp2blkloc = new int [n];
	if(!sp2blkloc) MemoryFail(funcname);

	for (i=0;i<nblksloc;i++) {
		for (j=blksobnd[i];j<blksobnd[i+1];j++) {
			sp2blkloc[j] = i;
		};
	};

	aspblk = mtro.BlockSparsityDiag (nblksloc, sp2blkloc);

	int nzjablk = aspblk.GetNzja ();

	int *jaweights;

	jaweights = new int [nzjablk];
	if(!jaweights) MemoryFail(funcname);

	for (i=0;i<nzjablk;i++) jaweights[i] = 0;

	int *nzblkarr;

	nzblkarr = new int [nblksloc];
	if(!nzblkarr) MemoryFail(funcname);

	piablk = aspblk.GetIa ();
	pjablk = aspblk.GetJa ();
	int *piamtr = mtro.GetIa ();
	int *pjamtr = mtro.GetJa ();

	int jj;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			nzblkarr[jblk] = 0;
		};
		for (i=blksobnd[iblk];i<blksobnd[iblk+1];i++) {
			for (j=piamtr[i];j<piamtr[i+1];j++) {
				jj = pjamtr[j];
				jblk = sp2blkloc[jj];
				if (iblk != jblk) nzblkarr[jblk]++;
			};
		};
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			jaweights[j] = nzblkarr[jblk];
		};
	};

	int *iablkflt;
	int *jablkflt;
	int *ndweights;
	int *adjweights;
	int *partition;

	iablkflt = new int [nblksloc+1];
	if(!iablkflt) MemoryFail(funcname);
	jablkflt = new int [nzjablk];
	if(!jablkflt) MemoryFail(funcname);
	ndweights = new int [nblksloc];
	if(!ndweights) MemoryFail(funcname);
	adjweights = new int [nzjablk];
	if(!adjweights) MemoryFail(funcname);
	partition = new int [nblksloc];
	if(!partition) MemoryFail(funcname);

	iablkflt[0] = 0;
	nzjablk = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			if (iblk != jblk) {
				jablkflt[nzjablk] = pjablk[j];
				adjweights[nzjablk] = jaweights[j];
				nzjablk++;
			};
		};
		iablkflt[iblk+1] = nzjablk;
		ndweights[iblk] = blksobnd[iblk+1]-blksobnd[iblk];
	};

// Partition

	int pzero = 0;
	int ptwo = 2;
	int pthree = 3;
	int pvol;

	int options[5];

	options[0] = 0;

	METIS_PartGraphRecursive (&nblksloc, iablkflt, jablkflt, 
										ndweights, adjweights, 
										&pthree, &pzero, &ptwo, 
										options, &pvol, partition);

// Modify disjoint parts

//	ModifyDisjointParts (nblksloc, iablkflt, jablkflt, ndweights, partition);

// Store partitioning

	int *partitionext;

	partitionext = new int [_nblks];
	if(!partitionext) MemoryFail(funcname);

	for (i=0;i<nblks1;i++) partitionext[i] = partition[0];
	for (i=0;i<nblks2;i++) partitionext[nblks1+i] = partition[1];
	for (i=0;i<nblksbnd;i++) partitionext[nblks1+nblks2+i] = partition[i+2];

	for (i=0;i<_nblks;i++) {
		iblkold = iorderblk[i];
		_blk2cpu[iblkold] = partitionext[i];
	};

// Search again for non-connected small subparts and add them into the appropriate part

	piablk = aspblk_save.GetIa ();
	pjablk = aspblk_save.GetJa ();

	iacomp = new int [_nblks+1];
	if(!iacomp) MemoryFail(funcname);
	iacompnd = new int [_nblks+1];
	if(!iacompnd) MemoryFail(funcname);
	listcomp = new int [_nblks];
	if(!listcomp) MemoryFail(funcname);
	maskcomp = new int [_nblks];
	if(!maskcomp) MemoryFail(funcname);

	int ipart2;

	for (ipart2=0;ipart2<3;ipart2++) {
		ipart = ipart2%2;
		int ncomp = 0;
		iacomp[0] = 0;
		iacompnd[0] = 0;
		int nz = 0;
		int nznd = 0;
		for (i=0;i<_nblks;i++) maskcomp[i] = -1;
		while (true) {
			int icheck = -1;
			for (i=0;i<_nblks;i++) {
				if (maskcomp[i] == -1 && _blk2cpu[i] == ipart) {
					maskcomp[i] = ncomp;
					listcomp[nz] = i;
					nz++;
					nznd += _blks[i+1]-_blks[i];
					ncomp++;
					icheck = 1;
					break;
				};
			};
			if (icheck == -1) break;
			for (i=iacomp[ncomp-1];i<nz;i++) {
				iblk = listcomp[i];
				for (j=piablk[iblk];j<piablk[iblk+1];j++) {
					int jj = pjablk[j];
					if (maskcomp[jj] == -1 && _blk2cpu[jj] == ipart) {
						maskcomp[jj] = ncomp-1;
						listcomp[nz] = jj;
						nz++;
						nznd += _blks[jj+1]-_blks[jj];
					};
				};
			};
			iacomp[ncomp] = nz;
			iacompnd[ncomp] = nznd;
		};
		if (ncomp > 1) {
			int icompmax = 0;
			int nzmax = iacompnd[1]-iacompnd[0];
			for (i=1;i<ncomp;i++) {
				if (iacompnd[i+1]-iacompnd[i] > nzmax) {
					icompmax = i;
					nzmax = iacompnd[i+1]-iacompnd[i];
				};
			};
			int ipart1 = (ipart+1)%2;
			for (i=0;i<ncomp;i++) {
				if (i != icompmax) {
					for (j=iacomp[i];j<iacomp[i+1];j++) {
						iblk = listcomp[j];
						_blk2cpu[iblk] = ipart1;
					};
				};
			};
		};
	};

	delete [] iacomp;
	delete [] iacompnd;
	delete [] listcomp;
	delete [] maskcomp;

// Free work arrays

	delete [] order;
	delete [] iorder;
	delete [] ind2blk;
	delete [] markblk;
	delete [] orderblk;
	delete [] iorderblk;
	delete [] blkso;
	delete [] orderloc;
	delete [] blksobnd;
	delete [] sp2blkloc;
	delete [] jaweights;
	delete [] nzblkarr;
	delete [] iablkflt;
	delete [] jablkflt;
	delete [] ndweights;
	delete [] adjweights;
	delete [] partition;
	delete [] partitionext;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute two sided boundary in parallel mode
//========================================================================================
void CSMatrix::ComputeTwoSidedBoundary (CMPIComm &_comm, // Compute two sided boundary in parallel mode
														int *_blkscpu, int *_partition,
														int &_nlistbnd, int *&_listbnd) const {

	const char *funcname = "ComputeTwoSidedBoundary";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create local final partition

	int *partition;

	partition = new int [nlist];
	if (!partition) MemoryFail(funcname);

	int i;

	for (i=0;i<nlist;i++) partition[i] = -2;

// For each block row determine if there are adjacent nodes from other cpus

	int itype, j, jj, ishift, jjloc;

	ishift = _blkscpu[myid];

	for (i=0;i<nlist;i++) {
		itype = 1;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < _blkscpu[myid] || jj >= _blkscpu[myid+1]) {
				itype = -1;
			} else {
				jjloc = jj-ishift;
				if (_partition[jjloc] != _partition[i]) {
					partition[i] = -1;
				};
			};
		};
		if (itype == 1) {
			partition[i] = _partition[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjloc = jj-ishift;
				if (_partition[jjloc] != _partition[i]) {
					partition[i] = -1;
					break;
				};
			};
		};
	};

// Create the list of sends to each cpu

	int nz = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < _blkscpu[myid] || jj >= _blkscpu[myid+1]) {
				nz++;
			};
		};
	};

	int *isenddata;

	isenddata = new int [nz*2];
	if (!isenddata) MemoryFail(funcname);

	nz = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < _blkscpu[myid] || jj >= _blkscpu[myid+1]) {
				isenddata[nz*2] = i;
				isenddata[nz*2+1] = jj;
				nz++;
			};
		};
	};

	int *jj2cpu;

	jj2cpu = new int [nz];
	if (!jj2cpu) MemoryFail(funcname);

	int icpu, ibegcpu, iendcpu;

	int icpuprev = 0;

	for (i=0;i<nz;i++) {
		jj = isenddata[i*2+1];
		if (jj >= _blkscpu[icpuprev] && jj < _blkscpu[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			ibegcpu = 0;
			iendcpu = nproc-1;
			while (ibegcpu != iendcpu) {
				icpu = (ibegcpu+iendcpu) / 2;
				if (jj >= _blkscpu[icpu] && jj < _blkscpu[icpu+1]) {
					ibegcpu = icpu;
					iendcpu = icpu;
				} else if (jj < _blkscpu[icpu]) {
					iendcpu = icpu-1;
				} else if (jj >= _blkscpu[icpu+1]) {
					ibegcpu = icpu+1;
				};
			};
			icpu = ibegcpu;
		};
		jj2cpu[i] = icpu;
		icpuprev = icpu;
	};

	int *nelems2cpu;

	nelems2cpu = new int [nproc+1];
	if (!nelems2cpu) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) nelems2cpu[i] = 0;

	int iproc;

	for (i=0;i<nz;i++) {
		iproc = jj2cpu[i];
		nelems2cpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) nelems2cpu[i+1] = nelems2cpu[i]+nelems2cpu[i+1];

	int *ibscpu;

	ibscpu = new int [nproc];
	if (!ibscpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibscpu[i] = nelems2cpu[i];

	int *isendarr;

	isendarr = new int [nz*2];
	if (!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nz;i++) {
		iproc = jj2cpu[i];
		k = ibscpu[iproc];
		isendarr[k*2] = isenddata[i*2];
		isendarr[k*2+1] = isenddata[i*2+1];
		ibscpu[iproc]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (nelems2cpu[i+1]-nelems2cpu[i])*2*sizeof(int);
		ObjSend[i] = (char *) (isendarr+2*nelems2cpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeTwoSidedBoundary: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeTwoSidedBoundary: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int jjnew, ni;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[2*j+1];
			jjloc = jj-ishift;
			jjnew = _partition[jjloc];
			piarr[2*j+1] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) {
		cout << " throw CSMatrix::ComputeTwoSidedBoundary: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeTwoSidedBoundary: Error in DataExchangeMPI";
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect partitioning data

	int irow, ipart;

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (2*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*2];
			ipart = piarr[j*2+1];
			if (ipart != _partition[irow]) partition[irow] = -1;
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Correct not initialized data

	for (i=0;i<nlist;i++) {
		if (partition[i] == -2) partition[i] = _partition[i];
	};

// Create the local list of boundary nodes

	_nlistbnd = 0;

	for (i=0;i<nlist;i++) {
		if (partition[i] == -1) {
			_nlistbnd++;
		};
	};

	_listbnd = new int [_nlistbnd];
	if(!_listbnd) MemoryFail(funcname);


	_nlistbnd = 0;

	for (i=0;i<nlist;i++) {
		if (partition[i] == -1) {
			_listbnd[_nlistbnd] = ishift+i;
			_nlistbnd++;
		};
	};

// Free work arrays

	delete [] partition;
	delete [] isenddata;
	delete [] jj2cpu;
	delete [] nelems2cpu;
	delete [] ibscpu;
	delete [] isendarr;

};

// Author: Kharchenko S.A.
// CSMatrix: Extend the nodes in parallel mode
//========================================================================================
void CSMatrix::ExtendNodes (CMPIComm &_comm, // Extend the nodes in parallel mode
										int *_blkscpu, 
										int &_nlistbnd, int *&_listbnd) const {

	const char *funcname = "ExtendNodes";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create the mask

	int ishift = _blkscpu[myid];

	int *imask;

	imask = new int [nlist];
	if(!imask) MemoryFail(funcname);

	int i, irow;

	for (i=0;i<nlist;i++) imask[i] = 1;

	for (i=0;i<_nlistbnd;i++) {
		irow = _listbnd[i]-ishift;
		imask[irow] = -1;
	};

	int nz = 0;

	int j, jj, jjloc;

	for (i=0;i<_nlistbnd;i++) {
		irow = _listbnd[i]-ishift;
		for (j=ia[irow];j<ia[irow+1];j++) {
			jj = ja[j];
			if (jj < _blkscpu[myid] || jj >= _blkscpu[myid+1]) {
				nz++;
			} else {
				jjloc = jj-ishift;
				imask[jjloc] = -1;
			};
		};
	};

	int *isenddata;

	isenddata = new int [nz];
	if (!isenddata) MemoryFail(funcname);

	nz = 0;

	for (i=0;i<_nlistbnd;i++) {
		irow = _listbnd[i]-ishift;
		for (j=ia[irow];j<ia[irow+1];j++) {
			jj = ja[j];
			if (jj < _blkscpu[myid] || jj >= _blkscpu[myid+1]) {
				isenddata[nz] = jj;
				nz++;
			};
		};
	};

	int *jj2cpu;

	jj2cpu = new int [nz];
	if (!jj2cpu) MemoryFail(funcname);

	int icpu, ibegcpu, iendcpu;

	int icpuprev = 0;

	for (i=0;i<nz;i++) {
		jj = isenddata[i];
		if (jj >= _blkscpu[icpuprev] && jj < _blkscpu[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			ibegcpu = 0;
			iendcpu = nproc-1;
			while (ibegcpu != iendcpu) {
				icpu = (ibegcpu+iendcpu) / 2;
				if (jj >= _blkscpu[icpu] && jj < _blkscpu[icpu+1]) {
					ibegcpu = icpu;
					iendcpu = icpu;
				} else if (jj < _blkscpu[icpu]) {
					iendcpu = icpu-1;
				} else if (jj >= _blkscpu[icpu+1]) {
					ibegcpu = icpu+1;
				};
			};
			icpu = ibegcpu;
		};
		jj2cpu[i] = icpu;
		icpuprev = icpu;
	};

	int *nelems2cpu;

	nelems2cpu = new int [nproc+1];
	if (!nelems2cpu) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) nelems2cpu[i] = 0;

	int iproc;

	for (i=0;i<nz;i++) {
		iproc = jj2cpu[i];
		nelems2cpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) nelems2cpu[i+1] = nelems2cpu[i]+nelems2cpu[i+1];

	int *ibscpu;

	ibscpu = new int [nproc];
	if (!ibscpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibscpu[i] = nelems2cpu[i];

	int *isendarr;

	isendarr = new int [nz];
	if (!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nz;i++) {
		iproc = jj2cpu[i];
		k = ibscpu[iproc];
		isendarr[k] = isenddata[i];
		ibscpu[iproc]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (nelems2cpu[i+1]-nelems2cpu[i])*sizeof(int);
		ObjSend[i] = (char *) (isendarr+nelems2cpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ExtendNodes: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ExtendNodes: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Register the result

	int ni;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / sizeof(int);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[j];
			jjloc = jj-ishift;
			imask[jjloc] = -1;
		};
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Recreate the list array

	delete [] _listbnd;

	_nlistbnd = 0;

	for (i=0;i<nlist;i++) {
		if (imask[i] == -1) _nlistbnd++;
	};

	_listbnd = new int [_nlistbnd];
	if(!_listbnd) MemoryFail(funcname);

	_nlistbnd = 0;

	for (i=0;i<nlist;i++) {
		if (imask[i] == -1) {
			_listbnd[_nlistbnd] = i+ishift;
			_nlistbnd++;
		};
	};

// Free work arrays

	delete [] imask;
	delete [] isenddata;
	delete [] jj2cpu;
	delete [] nelems2cpu;
	delete [] ibscpu;
	delete [] isendarr;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity
//========================================================================================
CSMatrix CSMatrix::ComputeBlockSparsity (CMPIComm &_comm, int _nblks, int *_blks) const { // Compute block sparsity

	const char *funcname = "ComputeBlockSparsity";

// For each list index compute its block number

	int *list2blk;

	list2blk = new int [nlist];
	if(!list2blk) MemoryFail(funcname);

	int iblkprev = 0;

	int i, jj, iblk, ibegblk, iendblk;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
			iblk = iblkprev;
		} else {
			ibegblk = 0;
			iendblk = _nblks-1;
			while (ibegblk != iendblk) {
				iblk = (ibegblk+iendblk) / 2;
				if (jj >= _blks[iblk] && jj < _blks[iblk+1]) {
					ibegblk = iblk;
					iendblk = iblk;
				} else if (jj < _blks[iblk]) {
					iendblk = iblk-1;
				} else if (jj >= _blks[iblk+1]) {
					ibegblk = iblk+1;
				};
			};
			iblk = ibegblk;
		};
		list2blk[i] = iblk;
		iblkprev = iblk;
	};

// Create the sublists with the same block numbers

	int *ialist2blk;
	int *iptr;
	int *jalist2blk;

	ialist2blk = new int [_nblks+1];
	if(!ialist2blk) MemoryFail(funcname);
	iptr = new int [_nblks];
	if(!iptr) MemoryFail(funcname);
	jalist2blk = new int [nlist];
	if(!jalist2blk) MemoryFail(funcname);

	for (i=0;i<=_nblks;i++) ialist2blk[i] = 0;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		ialist2blk[iblk+1]++;
	};

	for (i=0;i<_nblks;i++) ialist2blk[i+1] = ialist2blk[i]+ialist2blk[i+1];

	for (i=0;i<_nblks;i++) iptr[i] = ialist2blk[i];

	int k;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		k = iptr[iblk];
		jalist2blk[k] = i;
		iptr[iblk]++;
	};

// For each column index compute its block number

	int *ja2blk;

	ja2blk = new int [nzja];
	if(!ja2blk) MemoryFail(funcname);

/*	iblkprev = 0;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
			iblk = iblkprev;
		} else {
			ibegblk = 0;
			iendblk = _nblks-1;
			while (ibegblk != iendblk) {
				iblk = (ibegblk+iendblk) / 2;
				if (jj >= _blks[iblk] && jj < _blks[iblk+1]) {
					ibegblk = iblk;
					iendblk = iblk;
				} else if (jj < _blks[iblk]) {
					iendblk = iblk-1;
				} else if (jj >= _blks[iblk+1]) {
					ibegblk = iblk+1;
				};
			};
			iblk = ibegblk;
		};
		ja2blk[i] = iblk;
		iblkprev = iblk;
	};
*/
	CSMatrix::ComputeJa2Blk (_nblks, _blks,
										nlist, ia, ja,
										ja2blk);

// Compute the local sparsity

	int *iablk;
	int *imaskblk;
	int *listblk;

	iablk = new int [_nblks+1];
	if(!iablk) MemoryFail(funcname);
	imaskblk = new int [_nblks];
	if(!imaskblk) MemoryFail(funcname);
	listblk = new int [_nblks];
	if(!listblk) MemoryFail(funcname);

	int icycle = -1;

	for (i=0;i<_nblks;i++) imaskblk[i] = icycle;

	int nz = 0;
	iablk[0] = 0;

	int nlistblk, ilistnew, ilist, j, jblk;

	for (iblk=0;iblk<_nblks;iblk++) {
		icycle++;
		nlistblk = 0;
		for (ilistnew=ialist2blk[iblk];ilistnew<ialist2blk[iblk+1];ilistnew++) {
			ilist = jalist2blk[ilistnew];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jblk = ja2blk[j];
				if (imaskblk[jblk] != icycle) {
					nlistblk++;
					imaskblk[jblk] = icycle;
				};
			};
		};
		nz += nlistblk;
		iablk[iblk+1] = nz;
	};

	int *jablk;

	jablk = new int [nz];
	if(!jablk) MemoryFail(funcname);

	nz = 0;
	for (iblk=0;iblk<_nblks;iblk++) {
		icycle++;
		nlistblk = 0;
		for (ilistnew=ialist2blk[iblk];ilistnew<ialist2blk[iblk+1];ilistnew++) {
			ilist = jalist2blk[ilistnew];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jblk = ja2blk[j];
				if (imaskblk[jblk] != icycle) {
					listblk[nlistblk] = jblk;
					nlistblk++;
					imaskblk[jblk] = icycle;
				};
			};
		};
		sort (listblk, listblk + nlistblk);
		for (j=0;j<nlistblk;j++) jablk[nz+j] = listblk[j];
		nz += nlistblk;
	};

// Add all block sparsities

	int *iablktot;

	iablktot = new int [_nblks+1];
	if(!iablktot) MemoryFail(funcname);

	iablktot[0] = 0;
	for (i=0;i<_nblks;i++) iablktot[i+1] = iablk[i+1]-iablk[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblks+1, iablktot, iablktot);

	for (i=0;i<_nblks;i++) iablktot[i+1] = iablktot[i]+iablktot[i+1];

	for (i=0;i<_nblks;i++) {
		if (iablk[i+1]-iablk[i] > 0) {
			if (iablk[i+1]-iablk[i] != iablktot[i+1]-iablktot[i]) {
				cout << " throw CSMatrix::ComputeBlockSparsity: some block rows are distributed among cpu's " << endl;
				throw " CSMatrix::ComputeBlockSparsity: some block rows are distributed among cpu's ";
			};
		};
	};

	int nzjablktot = iablktot[_nblks];

	int *jablktot;

	jablktot = new int [nzjablktot];
	if(!jablktot) MemoryFail(funcname);

	for (i=0;i<nzjablktot;i++) jablktot[i] = 0;

	int ibeg, jloc;

	for (i=0;i<_nblks;i++) {
		if (iablk[i+1]-iablk[i] > 0) {
			ibeg = iablktot[i];
			for (j=iablk[i];j<iablk[i+1];j++) {
				jloc = j-iablk[i];
				jablktot[ibeg+jloc] = jablk[j];
			};
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nzjablktot, jablktot, jablktot);

// Check diagonal values

	int ndiagadd = 0;

	int icheck;

	for (i=0;i<_nblks;i++) {
		icheck = -1;
		for (j=iablktot[i];j<iablktot[i+1];j++) {
			jj = jablktot[j];
			if (jj == i) icheck = 1;
		};
		if (icheck == -1) ndiagadd++;
	};

// Store as a sparse matrix with diagonal values

	CSMatrix aloc (_nblks, nzjablktot+ndiagadd);

	for (i=0;i<_nblks;i++) aloc.list[i] = i;

	nz = 0;
	aloc.ia[0] = 0;

	for (i=0;i<_nblks;i++) {
		icheck = -1;
		for (j=iablktot[i];j<iablktot[i+1];j++) {
			jj = jablktot[j];
			aloc.ja[nz] = jj;
			nz++;
			if (jj == i) icheck = 1;
		};
		if (icheck == -1) {
			aloc.ja[nz] = i;
			nz++;
		};
		aloc.ia[i+1] = nz;
		ibeg = aloc.ia[i];
		sort(aloc.ja+ibeg,aloc.ja + aloc.ia[i+1]);
	};
	aloc.m = _nblks;
	aloc.n = _nblks;
	aloc.nsupc = _nblks;
	aloc.nsupr = _nblks;
	aloc.nlist = _nblks;
	aloc.nzja = nz;

// Symmetrize on return

	CSMatrix asymm = aloc.SymmMtr ();

// Free work arrays

	delete [] list2blk;
	delete [] ialist2blk;
	delete [] iptr;
	delete [] jalist2blk;
	delete [] ja2blk;
	delete [] iablk;
	delete [] imaskblk;
	delete [] listblk;
	delete [] jablk;
	delete [] iablktot;
	delete [] jablktot;

	return asymm;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity with local blocks partitioning
//========================================================================================
CSMatrix CSMatrix::ComputeBlockSparsityLocal (CMPIComm &_comm, int _nblks, int *_blks) const { // Compute block sparsity with local blocks partitioning

	const char *funcname = "ComputeBlockSparsityLocal";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global blocks partitionings

	int i;

	int *gblks;
	int *blkscpu;

	gblks = new int [nproc+1];
	if (!gblks) MemoryFail(funcname);
	blkscpu = new int [nproc+1];
	if (!blkscpu) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) gblks[i] = 0;
	for (i=0;i<=nproc;i++) blkscpu[i] = 0;

	gblks[myid+1] = _nblks;
	blkscpu[myid+1] = _blks[_nblks];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, gblks, gblks);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blkscpu, blkscpu);

	for (i=0;i<nproc;i++) gblks[i+1] = gblks[i]+gblks[i+1];
	for (i=0;i<nproc;i++) blkscpu[i+1] = blkscpu[i]+blkscpu[i+1];

// Create support arrays

	int *nd2blk;
	int *blk2pack;

	nd2blk = new int [nlist];
	if (!nd2blk) MemoryFail(funcname);
	blk2pack = new int [_nblks];
	if (!blk2pack) MemoryFail(funcname);

	int j;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			nd2blk[j] = i;
		};
	};

	for (i=0;i<_nblks;i++) {
		blk2pack[i] = gblks[myid]+i;
	};

// For each index find its cpu number

	int *ja2cpu;

	ja2cpu = new int [nzja];
	if (!ja2cpu) MemoryFail(funcname);

/*	int icpuprev = 0;

	int jj, icpu, ibegcpu, iendcpu;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= blkscpu[icpuprev] && jj < blkscpu[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			ibegcpu = 0;
			iendcpu = nproc-1;
			while (ibegcpu != iendcpu) {
				icpu = (ibegcpu+iendcpu) / 2;
				if (jj >= blkscpu[icpu] && jj < blkscpu[icpu+1]) {
					ibegcpu = icpu;
					iendcpu = icpu;
				} else if (jj < blkscpu[icpu]) {
					iendcpu = icpu-1;
				} else if (jj >= blkscpu[icpu+1]) {
					ibegcpu = icpu+1;
				};
			};
			icpu = ibegcpu;
		};
		ja2cpu[i] = icpu;
		icpuprev = icpu;
	};
*/
	CSMatrix::ComputeJa2Blk (nproc, blkscpu,
										nlist, ia, ja,
										ja2cpu);

// Prepare the data for each cpu

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail(funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int icpu;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int *isendarr;

	isendarr = new int [nzja*4];
	if (!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		k = iptr[icpu];
		isendarr[k*4] = i;
		isendarr[k*4+1] = ja[i];
		isendarr[k*4+2] = 0;
		isendarr[k*4+3] = 0;
		iptr[icpu]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*4*sizeof(int);
		ObjSend[i] = (char *) (isendarr+4*iacpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int ishift = blkscpu[myid];
	int ni, jjloc, jblk, ipart;
	int *piarr;

	int jj;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (4*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[4*j+1];
			jjloc = jj-ishift;
			jblk = nd2blk[jjloc];
			ipart = 0;
			piarr[4*j+1] = jblk+gblks[myid];
			piarr[4*j+2] = ipart;
			piarr[4*j+3] = blk2pack[jblk];
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store partitioning data

	int irow, ind, iblk, iblkpack;

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (4*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*4];
			iblk = piarr[j*4+1];
			ipart = piarr[j*4+2];
			iblkpack = piarr[j*4+3];
			isendarr[ind*4+1] = iblk;
			isendarr[ind*4+2] = ipart;
			isendarr[ind*4+3] = iblkpack;
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Create new packed block structure

	int nzpair, nzpairmax = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		nzpair = 0;
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			nzpair += (ia[j+1]-ia[j]);
		};
		if (nzpair > nzpairmax) nzpairmax = nzpair;
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nzpairmax];
	if(!iiarr) MemoryFail(funcname);

// Allocate containers and pointers

	int isize_cont_max = _nblks;
	int n_pointers_max = _nblks;

	int **ja_cont;
	int **w_cont;

	ja_cont = new int * [n_pointers_max];
	if(!ja_cont) MemoryFail(funcname);
	w_cont = new int * [n_pointers_max];
	if(!w_cont) MemoryFail(funcname);

	int *ia_cont;
	int **pja_cont;
	int **pw_cont;

	ia_cont = new int [n_pointers_max+1];
	if(!ia_cont) MemoryFail(funcname);
	pja_cont = new int * [n_pointers_max];
	if(!pja_cont) MemoryFail(funcname);
	pw_cont = new int * [n_pointers_max];
	if(!pw_cont) MemoryFail(funcname);

// Compute the sparsity and weights for remaining blocks and store results in containers

	int isize_cont = isize_cont_max+10;
	int n_cont_alloc = 0;
	int n_pointers = 0;

	ia_cont[0] = 0;

	int nz, jblknew, isizeloc, iblkprev, icount, iblknew;

	for (iblk=0;iblk<_nblks;iblk++) {

// Create block row

		nz = 0;
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jblknew = isendarr[j*4+3];
				if (nz != 0) {
					iblkprev = iiarr[nz-1].intvalue;
					if (iblkprev == jblknew) {
						iiarr[nz-1].int2value++;
					} else {
						iiarr[nz].intvalue = jblknew;
						iiarr[nz].int2value = 1;
						nz++;
					};
				} else {
					iiarr[nz].intvalue = jblknew;
					iiarr[nz].int2value = 1;
					nz++;
				};
			};
		};

		//qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nz);

		int nznew = 1;
		if (nz == 0) nznew = 0;

		for (i=1;i<nz;i++) {
			iblknew = iiarr[i].intvalue;
			icount = iiarr[i].int2value;
			iblkprev = iiarr[nznew-1].intvalue;
			if (iblknew == iblkprev) {
				iiarr[nznew-1].int2value += icount;
			} else {
				iiarr[nznew].intvalue = iblknew;
				iiarr[nznew].int2value = icount;
				nznew++;
			};
		};

// Store block row in containers

		if (isize_cont+nznew >= isize_cont_max) {
			isize_cont = 0;
			isizeloc = isize_cont_max;
			if (nznew > isizeloc) isizeloc = nznew;
			ja_cont[n_cont_alloc] = new int [isizeloc];
			if(!ja_cont[n_cont_alloc]) MemoryFail(funcname);
			w_cont[n_cont_alloc] = new int [isizeloc];
			if(!w_cont[n_cont_alloc]) MemoryFail(funcname);
			n_cont_alloc++;
		};

		pja_cont[n_pointers] = ja_cont[n_cont_alloc-1]+isize_cont;
		pw_cont [n_pointers] = w_cont [n_cont_alloc-1]+isize_cont;
		isize_cont += nznew;

		for (i=0;i<nznew;i++) {
			pja_cont[n_pointers][i] = iiarr[i].intvalue;
			pw_cont [n_pointers][i] = iiarr[i].int2value;
		};

		ia_cont[n_pointers+1] = ia_cont[n_pointers] + nznew;

		n_pointers++;

	};

// Move containers data into the sparsity

	int nzjaloc = ia_cont[n_pointers];

	CSMatrix aloc (n_pointers,nzjaloc);

	n_pointers = 0;

	aloc.ia[0] = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		aloc.list[n_pointers] = blk2pack[iblk];
		for (j=ia_cont[n_pointers];j<ia_cont[n_pointers+1];j++) {
			aloc.ja[j] = pja_cont[n_pointers][j-ia_cont[n_pointers]];
		};
		aloc.ia[n_pointers+1] = ia_cont[n_pointers+1];
		n_pointers++;
	};

	int ntot = blkscpu[nproc];

	aloc.m = ntot;
	aloc.n = ntot;
	aloc.nsupr = ntot;
	aloc.nsupc = ntot;
	aloc.nlist = n_pointers;
	aloc.nzja = aloc.ia[n_pointers];

// Free containers data

	for (i=0;i<n_cont_alloc;i++) {
		delete [] ja_cont[i];
		delete [] w_cont[i];
	};

	delete [] ja_cont;
	delete [] w_cont;
	delete [] pja_cont;
	delete [] pw_cont;

// Filter diagonal values

	int nlistloc = aloc.nlist;
	int nzloc = 0;

	delete [] ia_cont;

	ia_cont = new int [nlistloc+1];
	if(!ia_cont) MemoryFail(funcname);

	ia_cont[0] = 0;

	for (i=0;i<nlistloc;i++) {
		irow = aloc.list[i];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jj = aloc.ja[j];
			if (jj != irow) {
				aloc.ja[nzloc] = aloc.ja[j];
				nzloc++;
			};
		};
		ia_cont[i+1] = nzloc;
	};

	for (i=0;i<=nlistloc;i++) aloc.ia[i] = ia_cont[i];

	aloc.nzja = nzjaloc;

// Free work arrays

	delete [] gblks;
	delete [] blkscpu;
	delete [] nd2blk;
	delete [] blk2pack;
	delete [] ja2cpu;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;
	delete [] iiarr;
	delete [] ia_cont;

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity with weights and partitioning indices
//========================================================================================
CSMatrix CSMatrix::ComputeBlockSparsityWeights (CMPIComm &_comm, int _nblks, int *_blks, // Compute block sparsity with weights and partitioning indices
																	int *_blkscpu, int *_gblkscpu, int *_partition_blk,
																	int *&_wvertices, int *&_wedges) const {

	const char *funcname = "ComputeBlockSparsityWeights";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global blocks partitionings

	int i;

	int *gblks;

	gblks = new int [nproc+1];
	if (!gblks) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) gblks[i] = 0;

	gblks[myid+1] = _nblks;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, gblks, gblks);

	for (i=0;i<nproc;i++) gblks[i+1] = gblks[i]+gblks[i+1];

// Create support arrays

	int *nd2blk;
	int *blk2pack;

	nd2blk = new int [nlist];
	if (!nd2blk) MemoryFail(funcname);
	blk2pack = new int [_nblks];
	if (!blk2pack) MemoryFail(funcname);

	int j;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			nd2blk[j] = i;
		};
	};

	int nz = 0;
	if (myid == 0) nz += 2;

	for (i=0;i<_nblks;i++) {
		if (_partition_blk[i] == 0) {
			blk2pack[i] = 0;
		} else if (_partition_blk[i] == 1) {
			blk2pack[i] = 1;
		} else {
			blk2pack[i] = _gblkscpu[myid]+nz;
			nz++;
		};
	};

// For each index find its cpu number

	int *ja2cpu;

	ja2cpu = new int [nzja];
	if (!ja2cpu) MemoryFail(funcname);

/*	int icpuprev = 0;

	int jj, icpu, ibegcpu, iendcpu;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= _blkscpu[icpuprev] && jj < _blkscpu[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			ibegcpu = 0;
			iendcpu = nproc-1;
			while (ibegcpu != iendcpu) {
				icpu = (ibegcpu+iendcpu) / 2;
				if (jj >= _blkscpu[icpu] && jj < _blkscpu[icpu+1]) {
					ibegcpu = icpu;
					iendcpu = icpu;
				} else if (jj < _blkscpu[icpu]) {
					iendcpu = icpu-1;
				} else if (jj >= _blkscpu[icpu+1]) {
					ibegcpu = icpu+1;
				};
			};
			icpu = ibegcpu;
		};
		ja2cpu[i] = icpu;
		icpuprev = icpu;
	};
*/
	CSMatrix::ComputeJa2Blk (nproc, _blkscpu,
										nlist, ia, ja,
										ja2cpu);

// Prepare the data for each cpu

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail(funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int icpu;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int *isendarr;

	isendarr = new int [nzja*4];
	if (!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		k = iptr[icpu];
		isendarr[k*4] = i;
		isendarr[k*4+1] = ja[i];
		isendarr[k*4+2] = 0;
		isendarr[k*4+3] = 0;
		iptr[icpu]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*4*sizeof(int);
		ObjSend[i] = (char *) (isendarr+4*iacpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int ishift = _blkscpu[myid];
	int ni, jjloc, jblk, ipart;
	int *piarr;

	int jj;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (4*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[4*j+1];
			jjloc = jj-ishift;
			jblk = nd2blk[jjloc];
			ipart = _partition_blk[jblk];
			piarr[4*j+1] = jblk+gblks[myid];
			piarr[4*j+2] = ipart;
			piarr[4*j+3] = blk2pack[jblk];
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store partitioning data

	int irow, ind, iblk, iblkpack;

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (4*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*4];
			iblk = piarr[j*4+1];
			ipart = piarr[j*4+2];
			iblkpack = piarr[j*4+3];
			isendarr[ind*4+1] = iblk;
			isendarr[ind*4+2] = ipart;
			isendarr[ind*4+3] = iblkpack;
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Create new packed block structure

// Create local sparsities of first and second block row

	int nlistbnd = 0;

	int *listbnd;
	int *listbnd2blk;

	listbnd = new int [nlist];
	if(!listbnd) MemoryFail(funcname);
	listbnd2blk = new int [nlist];
	if(!listbnd2blk) MemoryFail(funcname);

	int n1tot = 0;
	int n2tot = 0;
	int nzpair = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_partition_blk[iblk] == 0) {
			n1tot += (_blks[iblk+1]-_blks[iblk]);
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ia[i];j<ia[i+1];j++) {
					ipart = isendarr[j*4+2];
					if (ipart == 1) {
						cout << " throw CSMatrix::ComputeBlockSparsityWeights: error in block splitting " << endl;
						throw " CSMatrix::ComputeBlockSparsityWeights: error in block splitting ";
					};
				};
			};
		};
		if (_partition_blk[iblk] == 1) {
			n2tot += (_blks[iblk+1]-_blks[iblk]);
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ia[i];j<ia[i+1];j++) {
					ipart = isendarr[j*4+2];
					if (ipart == 0) {
						cout << " throw CSMatrix::ComputeBlockSparsityWeights: error in block splitting " << endl;
						throw " CSMatrix::ComputeBlockSparsityWeights: error in block splitting ";
					};
				};
			};
		};
		if (_partition_blk[iblk] == -1) {
			for (j=_blks[iblk];j<_blks[iblk+1];j++) {
				listbnd[nlistbnd] = j;
				listbnd2blk[nlistbnd] = iblk;
				nlistbnd++;
				nzpair += (ia[j+1]-ia[j]);
			};
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &n1tot, &n1tot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &n2tot, &n2tot);

	CIntInt *iiarr;

	iiarr = new CIntInt [nzpair];
	if(!iiarr) MemoryFail(funcname);

	nz = 0;

	int ilist, iblkprev, iblknew, icount;

	for (ilist=0;ilist<nlistbnd;ilist++) {
		i = listbnd[ilist];
		iblk = listbnd2blk[ilist];
		for (j=ia[i];j<ia[i+1];j++) {
			ipart = isendarr[j*4+2];
			if (ipart == 0) {
				if (nz != 0) {
					iblkprev = iiarr[nz-1].intvalue;
					if (iblkprev == iblk) {
						iiarr[nz-1].int2value++;
					} else {
						iiarr[nz].intvalue = iblk;
						iiarr[nz].int2value = 1;
						nz++;
					};
				} else {
					iiarr[nz].intvalue = iblk;
					iiarr[nz].int2value = 1;
					nz++;
				};
			};
		};
	};

	//qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + nz);

	int nznew = 1;
	if (nz == 0) nznew = 0;

	for (i=1;i<nz;i++) {
		iblknew = iiarr[i].intvalue;
		icount = iiarr[i].int2value;
		iblkprev = iiarr[nznew-1].intvalue;
		if (iblknew == iblkprev) {
			iiarr[nznew-1].int2value += icount;
		} else {
			iiarr[nznew].intvalue = iblknew;
			iiarr[nznew].int2value = icount;
			nznew++;
		};
	};

	int nzj1 = nznew;

	int *listblk1;
	int *nzblk1;

	listblk1 = new int [nzj1];
	if(!listblk1) MemoryFail(funcname);
	nzblk1 = new int [nzj1];
	if(!nzblk1) MemoryFail(funcname);

	for (i=0;i<nzj1;i++) {
		listblk1[i] = iiarr[i].intvalue;
		nzblk1[i] = iiarr[i].int2value;
	};

	nz = 0;

	for (ilist=0;ilist<nlistbnd;ilist++) {
		i = listbnd[ilist];
		iblk = listbnd2blk[ilist];
		for (j=ia[i];j<ia[i+1];j++) {
			ipart = isendarr[j*4+2];
			if (ipart == 1) {
				if (nz != 0) {
					iblkprev = iiarr[nz-1].intvalue;
					if (iblkprev == iblk) {
						iiarr[nz-1].int2value++;
					} else {
						iiarr[nz].intvalue = iblk;
						iiarr[nz].int2value = 1;
						nz++;
					};
				} else {
					iiarr[nz].intvalue = iblk;
					iiarr[nz].int2value = 1;
					nz++;
				};
			};
		};
	};

	//qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + nz);

	nznew = 1;
	if (nz == 0) nznew = 0;

	for (i=1;i<nz;i++) {
		iblknew = iiarr[i].intvalue;
		icount = iiarr[i].int2value;
		iblkprev = iiarr[nznew-1].intvalue;
		if (iblknew == iblkprev) {
			iiarr[nznew-1].int2value += icount;
		} else {
			iiarr[nznew].intvalue = iblknew;
			iiarr[nznew].int2value = icount;
			nznew++;
		};
	};

	int nzj2 = nznew;

	int *listblk2;
	int *nzblk2;

	listblk2 = new int [nzj2];
	if(!listblk2) MemoryFail(funcname);
	nzblk2 = new int [nzj2];
	if(!nzblk2) MemoryFail(funcname);

	for (i=0;i<nzj2;i++) {
		listblk2[i] = iiarr[i].intvalue;
		nzblk2[i] = iiarr[i].int2value;
	};

// Allocate containers and pointers

	int isize_cont_max = _nblks;
	int n_pointers_max = _nblks;

	int **ja_cont;
	int **w_cont;

	ja_cont = new int * [n_pointers_max];
	if(!ja_cont) MemoryFail(funcname);
	w_cont = new int * [n_pointers_max];
	if(!w_cont) MemoryFail(funcname);

	int *ia_cont;
	int **pja_cont;
	int **pw_cont;

	ia_cont = new int [n_pointers_max+1];
	if(!ia_cont) MemoryFail(funcname);
	pja_cont = new int * [n_pointers_max];
	if(!pja_cont) MemoryFail(funcname);
	pw_cont = new int * [n_pointers_max];
	if(!pw_cont) MemoryFail(funcname);

// Compute the sparsity and weights for remaining blocks and store results in containers

	int isize_cont = isize_cont_max+10;
	int n_cont_alloc = 0;
	int n_pointers = 0;

	ia_cont[0] = 0;

	int jblknew, isizeloc;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_partition_blk[iblk] == -1) {

// Create block row

			nz = 0;
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ia[i];j<ia[i+1];j++) {
					jblknew = isendarr[j*4+3];
					if (nz != 0) {
						iblkprev = iiarr[nz-1].intvalue;
						if (iblkprev == jblknew) {
							iiarr[nz-1].int2value++;
						} else {
							iiarr[nz].intvalue = jblknew;
							iiarr[nz].int2value = 1;
							nz++;
						};
					} else {
						iiarr[nz].intvalue = jblknew;
						iiarr[nz].int2value = 1;
						nz++;
					};
				};
			};

			//qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
			sort (iiarr, iiarr + nz);

			int nznew = 1;
			if (nz == 0) nznew = 0;

			for (i=1;i<nz;i++) {
				iblknew = iiarr[i].intvalue;
				icount = iiarr[i].int2value;
				iblkprev = iiarr[nznew-1].intvalue;
				if (iblknew == iblkprev) {
					iiarr[nznew-1].int2value += icount;
				} else {
					iiarr[nznew].intvalue = iblknew;
					iiarr[nznew].int2value = icount;
					nznew++;
				};
			};

// Store block row in containers

			if (isize_cont+nznew >= isize_cont_max) {
				isize_cont = 0;
				isizeloc = isize_cont_max;
				if (nznew > isizeloc) isizeloc = nznew;
				ja_cont[n_cont_alloc] = new int [isizeloc];
				if(!ja_cont[n_cont_alloc]) MemoryFail(funcname);
				w_cont[n_cont_alloc] = new int [isizeloc];
				if(!w_cont[n_cont_alloc]) MemoryFail(funcname);
				n_cont_alloc++;
			};

			pja_cont[n_pointers] = ja_cont[n_cont_alloc-1]+isize_cont;
			pw_cont [n_pointers] = w_cont [n_cont_alloc-1]+isize_cont;
			isize_cont += nznew;

			for (i=0;i<nznew;i++) {
				pja_cont[n_pointers][i] = iiarr[i].intvalue;
				pw_cont [n_pointers][i] = iiarr[i].int2value;
			};

			ia_cont[n_pointers+1] = ia_cont[n_pointers] + nznew;

			n_pointers++;

		};
	};

// Move containers data into the sparsity

	int nzjaloc = ia_cont[n_pointers];

	CSMatrix aloc (n_pointers,nzjaloc);

	_wvertices = new int [n_pointers];
	if(!_wvertices) MemoryFail(funcname);
	_wedges = new int [nzjaloc];
	if(!_wedges) MemoryFail(funcname);

	n_pointers = 0;

	aloc.ia[0] = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_partition_blk[iblk] == -1) {
			aloc.list[n_pointers] = blk2pack[iblk];
			for (j=ia_cont[n_pointers];j<ia_cont[n_pointers+1];j++) {
				aloc.ja[j] = pja_cont[n_pointers][j-ia_cont[n_pointers]];
				_wedges[j] = pw_cont[n_pointers][j-ia_cont[n_pointers]];
			};
			_wvertices[n_pointers] = _blks[iblk+1]-_blks[iblk];
			aloc.ia[n_pointers+1] = ia_cont[n_pointers+1];
			n_pointers++;
		};
	};

	int ntot = _gblkscpu[nproc];

	aloc.m = ntot;
	aloc.n = ntot;
	aloc.nsupr = ntot;
	aloc.nsupc = ntot;
	aloc.nlist = n_pointers;
	aloc.nzja = aloc.ia[n_pointers];

// Free containers data

	for (i=0;i<n_cont_alloc;i++) {
		delete [] ja_cont[i];
		delete [] w_cont[i];
	};

	delete [] ja_cont;
	delete [] w_cont;
	delete [] pja_cont;
	delete [] pw_cont;

// Send the first and the second row onto the first cpu and add

	int *nsize1row;
	int *nsize2row;

	nsize1row = new int [nproc+1];
	if(!nsize1row) MemoryFail(funcname);
	nsize2row = new int [nproc+1];
	if(!nsize2row) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) nsize1row[i] = 0;
	for (i=0;i<=nproc;i++) nsize2row[i] = 0;

	nsize1row[myid+1] = nzj1;
	nsize2row[myid+1] = nzj2;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, nsize1row, nsize1row);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, nsize2row, nsize2row);

	for (i=0;i<nproc;i++) nsize1row[i+1] = nsize1row[i]+nsize1row[i+1];
	for (i=0;i<nproc;i++) nsize2row[i+1] = nsize2row[i]+nsize2row[i+1];

	delete [] isendarr;

	isendarr = new int [2*(nzj1+nzj2)];
	if(!isendarr) MemoryFail(funcname);

	for (i=0;i<nzj1;i++) isendarr[i] = blk2pack[listblk1[i]];
	for (i=0;i<nzj1;i++) isendarr[nzj1+i] = nzblk1[i];
	for (i=0;i<nzj2;i++) isendarr[2*nzj1+i] = blk2pack[listblk2[i]];
	for (i=0;i<nzj2;i++) isendarr[2*nzj1+nzj2+i] = nzblk2[i];

	NObjSend = 1;

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

	ObjTypeSend[0] = 1;
	ObjIDSend[0] = 1;
	CpuIDSend[0] = 0;
	ObjSizeSend[0] = (2*(nzj1+nzj2))*sizeof(int);
	ObjSend[0] = (char *) isendarr;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Get received data

	if (myid == 0) {

		delete [] iiarr;

		int nz1tot = nsize1row[nproc];
		int nz2tot = nsize2row[nproc];

		int nzmax = nz1tot;
		if (nz2tot > nz1tot) nzmax = nz2tot;

		iiarr = new CIntInt [nzmax];
		if(!iiarr) MemoryFail(funcname);

		int *piarr;

		int iproc, ibeg, ni1loc, ni2loc, nzsum;

		for (i=0;i<NObjRecv;i++) {
			iproc = CpuIDRecv[i];
			ibeg = nsize1row[iproc];
			ni1loc = nsize1row[iproc+1]-nsize1row[iproc];
			piarr = (int *) ObjRecv[i];
			for (j=0;j<ni1loc;j++) {
				jblk = piarr[j];
				nzsum = piarr[ni1loc+j];
				iiarr[ibeg+j].intvalue = jblk;
				iiarr[ibeg+j].int2value = nzsum;
			};
		};

// 		qsort (iiarr, nz1tot, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nz1tot);

		int nznew = 1;
		if (nz1tot == 0) nznew = 0;

		for (i=1;i<nz1tot;i++) {
			iblknew = iiarr[i].intvalue;
			icount = iiarr[i].int2value;
			iblkprev = iiarr[nznew-1].intvalue;
			if (iblknew == iblkprev) {
				iiarr[nznew-1].int2value += icount;
			} else {
				iiarr[nznew].intvalue = iblknew;
				iiarr[nznew].int2value = icount;
				nznew++;
			};
		};

		delete [] listblk1;
		delete [] nzblk1;

		nzj1 = nznew;

		listblk1 = new int [nzj1];
		if(!listblk1) MemoryFail(funcname);
		nzblk1 = new int [nzj1];
		if(!nzblk1) MemoryFail(funcname);

		for (i=0;i<nzj1;i++) {
			listblk1[i] = iiarr[i].intvalue;
			nzblk1[i] = iiarr[i].int2value;
		};

		for (i=0;i<NObjRecv;i++) {
			iproc = CpuIDRecv[i];
			ibeg = nsize2row[iproc];
			ni1loc = nsize1row[iproc+1]-nsize1row[iproc];
			ni2loc = nsize2row[iproc+1]-nsize2row[iproc];
			piarr = (int *) ObjRecv[i];
			for (j=0;j<ni2loc;j++) {
				jblk = piarr[ni1loc*2+j];
				nzsum = piarr[ni1loc*2+ni2loc+j];
				iiarr[ibeg+j].intvalue = jblk;
				iiarr[ibeg+j].int2value = nzsum;
			};
		};

		//qsort (iiarr, nz2tot, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nz2tot);

		nznew = 1;
		if (nz2tot == 0) nznew = 0;

		for (i=1;i<nz2tot;i++) {
			iblknew = iiarr[i].intvalue;
			icount = iiarr[i].int2value;
			iblkprev = iiarr[nznew-1].intvalue;
			if (iblknew == iblkprev) {
				iiarr[nznew-1].int2value += icount;
			} else {
				iiarr[nznew].intvalue = iblknew;
				iiarr[nznew].int2value = icount;
				nznew++;
			};
		};

		delete [] listblk2;
		delete [] nzblk2;

		nzj2 = nznew;

		listblk2 = new int [nzj2];
		if(!listblk2) MemoryFail(funcname);
		nzblk2 = new int [nzj2];
		if(!nzblk2) MemoryFail(funcname);

		for (i=0;i<nzj2;i++) {
			listblk2[i] = iiarr[i].intvalue;
			nzblk2[i] = iiarr[i].int2value;
		};

	};

// Free receive data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Modify data on the first cpu

	if (myid == 0) {

		int nlistloc = aloc.nlist+2;
		nzjaloc = aloc.nzja + nzj1 + nzj2;

		CSMatrix anew (nlistloc,nzjaloc);

		int *wverticesloc;
		int *wedgesloc;

		wverticesloc = new int [nlistloc];
		if(!wverticesloc) MemoryFail(funcname);
		wedgesloc = new int [nzjaloc];
		if(!wedgesloc) MemoryFail(funcname);

		anew.list[0] = 0;
		anew.list[1] = 1;
		for (i=0;i<aloc.nlist;i++) anew.list[i+2] = aloc.list[i];

		for (i=0;i<nzj1;i++) anew.ja[i] = listblk1[i];
		for (i=0;i<nzj2;i++) anew.ja[nzj1+i] = listblk2[i];
		for (i=0;i<aloc.nzja;i++) anew.ja[nzj1+nzj2+i] = aloc.ja[i];

		for (i=0;i<nzj1;i++) wedgesloc[i] = nzblk1[i];
		for (i=0;i<nzj2;i++) wedgesloc[nzj1+i] = nzblk2[i];
		for (i=0;i<aloc.nzja;i++) wedgesloc[nzj1+nzj2+i] = _wedges[i];

		wverticesloc[0] = n1tot;
		wverticesloc[1] = n2tot;
		for (i=0;i<aloc.nlist;i++) wverticesloc[2+i] = _wvertices[i];

		anew.ia[0] = 0;
		anew.ia[1] = nzj1;
		anew.ia[2] = nzj1+nzj2;
		for (i=0;i<aloc.nlist;i++) anew.ia[i+3] = anew.ia[i+2]+(aloc.ia[i+1]-aloc.ia[i]);

		anew.m = ntot;
		anew.n = ntot;
		anew.nsupr = ntot;
		anew.nsupc = ntot;
		anew.nlist = nlistloc;
		anew.nzja = nzjaloc;

		aloc = anew;

		delete [] _wedges;
		_wedges = wedgesloc;

		delete [] _wvertices;
		_wvertices = wverticesloc;

	};

// Filter diagonal values

	int nlistloc = aloc.nlist;
	int nzloc = 0;

	delete [] ia_cont;

	ia_cont = new int [nlistloc+1];
	if(!ia_cont) MemoryFail(funcname);

	ia_cont[0] = 0;

	for (i=0;i<nlistloc;i++) {
		irow = aloc.list[i];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jj = aloc.ja[j];
			if (jj != irow) {
				aloc.ja[nzloc] = aloc.ja[j];
				_wedges[nzloc] = _wedges[j];
				nzloc++;
			};
		};
		ia_cont[i+1] = nzloc;
	};

	for (i=0;i<=nlistloc;i++) aloc.ia[i] = ia_cont[i];

	aloc.nzja = nzjaloc;

// Free work arrays

	delete [] gblks;
	delete [] nd2blk;
	delete [] blk2pack;
	delete [] ja2cpu;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;
	delete [] listbnd;
	delete [] listbnd2blk;
	delete [] iiarr;
	delete [] listblk1;
	delete [] nzblk1;
	delete [] listblk2;
	delete [] nzblk2;
	delete [] nsize1row;
	delete [] nsize2row;
	delete [] ia_cont;

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity with weights and partitioning indices
//========================================================================================
CSMatrix CSMatrix::ComputeBlockSparsityWeightsSets (CMPIComm &_comm, int _nblks, int *_blks, int _nsets, // Compute block sparsity with weights and partitioning indices
																	int *_blkscpu, int *_gblkscpu, int *_partition_blk,
																	int *&_wvertices, int *&_wedges) const {

	const char *funcname = "ComputeBlockSparsityWeightsSets";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global blocks partitionings

	int i;

	int *gblks;

	gblks = new int [nproc+1];
	if (!gblks) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) gblks[i] = 0;

	gblks[myid+1] = _nblks;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, gblks, gblks);

	for (i=0;i<nproc;i++) gblks[i+1] = gblks[i]+gblks[i+1];

// Create support arrays

	int *nd2blk;
	int *blk2pack;

	nd2blk = new int [nlist];
	if (!nd2blk) MemoryFail(funcname);
	blk2pack = new int [_nblks];
	if (!blk2pack) MemoryFail(funcname);

	int j;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			nd2blk[j] = i;
		};
	};

	int nz = 0;
	if (myid == 0) nz += _nsets;

	for (i=0;i<_nblks;i++) {
		if (_partition_blk[i] >= 0) {
			blk2pack[i] = _partition_blk[i];
		} else {
			blk2pack[i] = _gblkscpu[myid]+nz;
			nz++;
		};
	};

// For each index find its cpu number

	int *ja2cpu;

	ja2cpu = new int [nzja];
	if (!ja2cpu) MemoryFail(funcname);

/*	int icpuprev = 0;

	int jj, icpu, ibegcpu, iendcpu;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= _blkscpu[icpuprev] && jj < _blkscpu[icpuprev+1]) {
			icpu = icpuprev;
		} else {
			ibegcpu = 0;
			iendcpu = nproc-1;
			while (ibegcpu != iendcpu) {
				icpu = (ibegcpu+iendcpu) / 2;
				if (jj >= _blkscpu[icpu] && jj < _blkscpu[icpu+1]) {
					ibegcpu = icpu;
					iendcpu = icpu;
				} else if (jj < _blkscpu[icpu]) {
					iendcpu = icpu-1;
				} else if (jj >= _blkscpu[icpu+1]) {
					ibegcpu = icpu+1;
				};
			};
			icpu = ibegcpu;
		};
		ja2cpu[i] = icpu;
		icpuprev = icpu;
	};
*/
	CSMatrix::ComputeJa2Blk (nproc, _blkscpu,
										nlist, ia, ja,
										ja2cpu);

// Prepare the data for each cpu

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail(funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int icpu;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int *isendarr;

	isendarr = new int [nzja*4];
	if (!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		k = iptr[icpu];
		isendarr[k*4] = i;
		isendarr[k*4+1] = ja[i];
		isendarr[k*4+2] = 0;
		isendarr[k*4+3] = 0;
		iptr[icpu]++;
	};

// Exchange the data

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*4*sizeof(int);
		ObjSend[i] = (char *) (isendarr+4*iacpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int ishift = _blkscpu[myid];
	int ni, jjloc, jblk, ipart;
	int *piarr;

	int jj;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (4*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[4*j+1];
			jjloc = jj-ishift;
			jblk = nd2blk[jjloc];
			ipart = _partition_blk[jblk];
			piarr[4*j+1] = jblk+gblks[myid];
			piarr[4*j+2] = ipart;
			piarr[4*j+3] = blk2pack[jblk];
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store partitioning data

	int irow, ind, iblk, iblkpack;

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (4*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*4];
			iblk = piarr[j*4+1];
			ipart = piarr[j*4+2];
			iblkpack = piarr[j*4+3];
			isendarr[ind*4+1] = iblk;
			isendarr[ind*4+2] = ipart;
			isendarr[ind*4+3] = iblkpack;
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Create new packed block structure

// Create local sparsities of first and second block row

	int nlistbnd = 0;

	int *nitotarr;
	int *listbnd;
	int *listbnd2blk;

	nitotarr = new int [_nsets];
	if(!nitotarr) MemoryFail(funcname);
	listbnd = new int [nlist];
	if(!listbnd) MemoryFail(funcname);
	listbnd2blk = new int [nlist];
	if(!listbnd2blk) MemoryFail(funcname);

	for (i=0;i<_nsets;i++) nitotarr[i] = 0;

	int nzpair = 0;

	int jpart;

	for (iblk=0;iblk<_nblks;iblk++) {
		ipart = _partition_blk[iblk];
		if (ipart >= 0) {
			nitotarr[ipart] += (_blks[iblk+1]-_blks[iblk]);
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ia[i];j<ia[i+1];j++) {
					jpart = isendarr[j*4+2];
					if (jpart >= 0 && jpart != ipart) {
						cout << " throw CSMatrix::ComputeBlockSparsityWeights: error in block splitting " << endl;
						throw " CSMatrix::ComputeBlockSparsityWeights: error in block splitting ";
					};
				};
			};
		} else {
			for (j=_blks[iblk];j<_blks[iblk+1];j++) {
				listbnd[nlistbnd] = j;
				listbnd2blk[nlistbnd] = iblk;
				nlistbnd++;
				nzpair += (ia[j+1]-ia[j]);
			};
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nsets, nitotarr, nitotarr);

	CIntInt *iiarr;

	iiarr = new CIntInt [nzpair];
	if(!iiarr) MemoryFail(funcname);

	int ilist, iblkprev, iblknew, icount;

	int *nzjarr;
	int **plistblk;
	int **pnzblk;

	nzjarr = new int [_nsets];
	if(!nzjarr) MemoryFail(funcname);
	plistblk = new int * [_nsets];
	if(!plistblk) MemoryFail(funcname);
	pnzblk = new int * [_nsets];
	if(!pnzblk) MemoryFail(funcname);

	int iset, nznew;

	for (iset=0;iset<_nsets;iset++) {

		nz = 0;

		for (ilist=0;ilist<nlistbnd;ilist++) {
			i = listbnd[ilist];
			iblk = listbnd2blk[ilist];
			for (j=ia[i];j<ia[i+1];j++) {
				ipart = isendarr[j*4+2];
				if (ipart == iset) {
					if (nz != 0) {
						iblkprev = iiarr[nz-1].intvalue;
						if (iblkprev == iblk) {
							iiarr[nz-1].int2value++;
						} else {
							iiarr[nz].intvalue = iblk;
							iiarr[nz].int2value = 1;
							nz++;
						};
					} else {
						iiarr[nz].intvalue = iblk;
						iiarr[nz].int2value = 1;
						nz++;
					};
				};
			};
		};

		//qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nz);

		nznew = 1;
		if (nz == 0) nznew = 0;

		for (i=1;i<nz;i++) {
			iblknew = iiarr[i].intvalue;
			icount = iiarr[i].int2value;
			iblkprev = iiarr[nznew-1].intvalue;
			if (iblknew == iblkprev) {
				iiarr[nznew-1].int2value += icount;
			} else {
				iiarr[nznew].intvalue = iblknew;
				iiarr[nznew].int2value = icount;
				nznew++;
			};
		};

		nzjarr[iset] = nznew;

		plistblk[iset] = new int [nznew];
		if(!plistblk[iset]) MemoryFail(funcname);
		pnzblk[iset] = new int [nznew];
		if(!pnzblk[iset]) MemoryFail(funcname);

		for (i=0;i<nznew;i++) {
			plistblk[iset][i] = iiarr[i].intvalue;
			pnzblk[iset][i] = iiarr[i].int2value;
		};

	};

// Allocate containers and pointers

	int isize_cont_max = _nblks;
	int n_pointers_max = _nblks;

	int **ja_cont;
	int **w_cont;

	ja_cont = new int * [n_pointers_max];
	if(!ja_cont) MemoryFail(funcname);
	w_cont = new int * [n_pointers_max];
	if(!w_cont) MemoryFail(funcname);

	int *ia_cont;
	int **pja_cont;
	int **pw_cont;

	ia_cont = new int [n_pointers_max+1];
	if(!ia_cont) MemoryFail(funcname);
	pja_cont = new int * [n_pointers_max];
	if(!pja_cont) MemoryFail(funcname);
	pw_cont = new int * [n_pointers_max];
	if(!pw_cont) MemoryFail(funcname);

// Compute the sparsity and weights for remaining blocks and store results in containers

	int isize_cont = isize_cont_max+10;
	int n_cont_alloc = 0;
	int n_pointers = 0;

	ia_cont[0] = 0;

	int jblknew, isizeloc;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_partition_blk[iblk] == -1) {

// Create block row

			nz = 0;
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ia[i];j<ia[i+1];j++) {
					jblknew = isendarr[j*4+3];
					if (nz != 0) {
						iblkprev = iiarr[nz-1].intvalue;
						if (iblkprev == jblknew) {
							iiarr[nz-1].int2value++;
						} else {
							iiarr[nz].intvalue = jblknew;
							iiarr[nz].int2value = 1;
							nz++;
						};
					} else {
						iiarr[nz].intvalue = jblknew;
						iiarr[nz].int2value = 1;
						nz++;
					};
				};
			};

// 			qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
			sort (iiarr, iiarr + nz);

			int nznew = 1;
			if (nz == 0) nznew = 0;

			for (i=1;i<nz;i++) {
				iblknew = iiarr[i].intvalue;
				icount = iiarr[i].int2value;
				iblkprev = iiarr[nznew-1].intvalue;
				if (iblknew == iblkprev) {
					iiarr[nznew-1].int2value += icount;
				} else {
					iiarr[nznew].intvalue = iblknew;
					iiarr[nznew].int2value = icount;
					nznew++;
				};
			};

// Store block row in containers

			if (isize_cont+nznew >= isize_cont_max) {
				isize_cont = 0;
				isizeloc = isize_cont_max;
				if (nznew > isizeloc) isizeloc = nznew;
				ja_cont[n_cont_alloc] = new int [isizeloc];
				if(!ja_cont[n_cont_alloc]) MemoryFail(funcname);
				w_cont[n_cont_alloc] = new int [isizeloc];
				if(!w_cont[n_cont_alloc]) MemoryFail(funcname);
				n_cont_alloc++;
			};

			pja_cont[n_pointers] = ja_cont[n_cont_alloc-1]+isize_cont;
			pw_cont [n_pointers] = w_cont [n_cont_alloc-1]+isize_cont;
			isize_cont += nznew;

			for (i=0;i<nznew;i++) {
				pja_cont[n_pointers][i] = iiarr[i].intvalue;
				pw_cont [n_pointers][i] = iiarr[i].int2value;
			};

			ia_cont[n_pointers+1] = ia_cont[n_pointers] + nznew;

			n_pointers++;

		};
	};

// Move containers data into the sparsity

	int nzjaloc = ia_cont[n_pointers];

	CSMatrix aloc (n_pointers,nzjaloc);

	_wvertices = new int [n_pointers];
	if(!_wvertices) MemoryFail(funcname);
	_wedges = new int [nzjaloc];
	if(!_wedges) MemoryFail(funcname);

	n_pointers = 0;

	aloc.ia[0] = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_partition_blk[iblk] == -1) {
			aloc.list[n_pointers] = blk2pack[iblk];
			for (j=ia_cont[n_pointers];j<ia_cont[n_pointers+1];j++) {
				aloc.ja[j] = pja_cont[n_pointers][j-ia_cont[n_pointers]];
				_wedges[j] = pw_cont[n_pointers][j-ia_cont[n_pointers]];
			};
			_wvertices[n_pointers] = _blks[iblk+1]-_blks[iblk];
			aloc.ia[n_pointers+1] = ia_cont[n_pointers+1];
			n_pointers++;
		};
	};

	int ntot = _gblkscpu[nproc];

	aloc.m = ntot;
	aloc.n = ntot;
	aloc.nsupr = ntot;
	aloc.nsupc = ntot;
	aloc.nlist = n_pointers;
	aloc.nzja = aloc.ia[n_pointers];

// Free containers data

	for (i=0;i<n_cont_alloc;i++) {
		delete [] ja_cont[i];
		delete [] w_cont[i];
	};

	delete [] ja_cont;
	delete [] w_cont;
	delete [] pja_cont;
	delete [] pw_cont;

// Send the set of initial rows onto the first cpu and add

	int *nsizerowarr;

	nsizerowarr = new int [_nsets*(nproc+1)];
	if(!nsizerowarr) MemoryFail(funcname);

	for (i=0;i<(nproc+1)*_nsets;i++) nsizerowarr[i] = 0;

	for (i=0;i<_nsets;i++) {
		nsizerowarr[(nproc+1)*i+myid+1] = nzjarr[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, (nproc+1)*_nsets, nsizerowarr, nsizerowarr);

	for (iset=0;iset<_nsets;iset++) {
		for (i=0;i<nproc;i++) nsizerowarr[iset*(nproc+1)+i+1] = nsizerowarr[iset*(nproc+1)+i]+nsizerowarr[iset*(nproc+1)+i+1];
	};

	delete [] isendarr;

	int nzj = 0;

	for (i=0;i<_nsets;i++) nzj += nzjarr[i];

	isendarr = new int [2*nzj];
	if(!isendarr) MemoryFail(funcname);

	nzj = 0;

	for (i=0;i<_nsets;i++) {
		for (j=0;j<nzjarr[i];j++) {
			isendarr[nzj*2+j] = blk2pack[plistblk[i][j]];
			isendarr[nzj*2+nzjarr[i]+j] = pnzblk[i][j];
		};
		nzj += nzjarr[i];
	};

	NObjSend = 1;

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

	ObjTypeSend[0] = 1;
	ObjIDSend[0] = 1;
	CpuIDSend[0] = 0;
	ObjSizeSend[0] = (2*nzj)*sizeof(int);
	ObjSend[0] = (char *) isendarr;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ComputeBlockSparsityWeights: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Get received data

	int nzloc, niloc;

	if (myid == 0) {

		int *iptrcpu;

		iptrcpu = new int [nproc];
		if(!iptrcpu) MemoryFail(funcname);

		for (i=0;i<nproc;i++) iptrcpu[i] = 0;

		delete [] iiarr;

		int nzmax = 0;
		for (i=0;i<_nsets;i++) {
			nzloc = nsizerowarr[i*(nproc+1)+nproc];
			if (nzloc > nzmax) nzmax = nzloc;
		};

		iiarr = new CIntInt [nzmax];
		if(!iiarr) MemoryFail(funcname);

		int *piarr;

		int iproc, ibeg, nzsum;

		for (iset=0;iset<_nsets;iset++) {

			nzloc = 0;

			for (i=0;i<NObjRecv;i++) {
				iproc = CpuIDRecv[i];
				ibeg = iptrcpu[iproc];
				niloc = nsizerowarr[iset*(nproc+1)+iproc+1]-nsizerowarr[iset*(nproc+1)+iproc];
				piarr = (int *) ObjRecv[i];
				for (j=0;j<niloc;j++) {
					jblk = piarr[ibeg*2+j];
					nzsum = piarr[ibeg*2+niloc+j];
					iiarr[nzloc].intvalue = jblk;
					iiarr[nzloc].int2value = nzsum;
					nzloc++;
				};
				iptrcpu[iproc] += niloc;
			};

			//qsort (iiarr, nzloc, sizeof(CIntInt), CIntInt::CompareIntInt);
			sort (iiarr, iiarr + nzloc);

			nznew = 1;
			if (nzloc == 0) nznew = 0;

			for (i=1;i<nzloc;i++) {
				iblknew = iiarr[i].intvalue;
				icount = iiarr[i].int2value;
				iblkprev = iiarr[nznew-1].intvalue;
				if (iblknew == iblkprev) {
					iiarr[nznew-1].int2value += icount;
				} else {
					iiarr[nznew].intvalue = iblknew;
					iiarr[nznew].int2value = icount;
					nznew++;
				};
			};

			delete [] plistblk[iset];
			delete [] pnzblk[iset];

			nzjarr[iset] = nznew;

			plistblk[iset] = new int [nznew];
			if(!plistblk[iset]) MemoryFail(funcname);
			pnzblk[iset] = new int [nznew];
			if(!pnzblk[iset]) MemoryFail(funcname);

			for (i=0;i<nznew;i++) {
				plistblk[iset][i] = iiarr[i].intvalue;
				pnzblk[iset][i] = iiarr[i].int2value;
			};

		};

		delete [] iptrcpu;

	};

// Free receive data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Modify data on the first cpu

	if (myid == 0) {

		int nlistloc = aloc.nlist+_nsets;
		nzjaloc = aloc.nzja;
		int nzjaini = 0;

		for (i=0;i<_nsets;i++) {
			nzjaini += nzjarr[i];
		};

		nzjaloc += nzjaini;

		CSMatrix anew (nlistloc,nzjaloc);

		int *wverticesloc;
		int *wedgesloc;

		wverticesloc = new int [nlistloc];
		if(!wverticesloc) MemoryFail(funcname);
		wedgesloc = new int [nzjaloc];
		if(!wedgesloc) MemoryFail(funcname);

		for (i=0;i<_nsets;i++) anew.list[i] = i;
		for (i=0;i<aloc.nlist;i++) anew.list[i+_nsets] = aloc.list[i];

		nzjaini = 0;
		for (i=0;i<_nsets;i++) {
			for (j=0;j<nzjarr[i];j++) {
				anew.ja[nzjaini] = plistblk[i][j];
				nzjaini++;
			};
		};

		for (i=0;i<aloc.nzja;i++) anew.ja[nzjaini+i] = aloc.ja[i];

		nzjaini = 0;
		for (i=0;i<_nsets;i++) {
			for (j=0;j<nzjarr[i];j++) {
				wedgesloc[nzjaini] = pnzblk[i][j];
				nzjaini++;
			};
		};

		for (i=0;i<aloc.nzja;i++) wedgesloc[nzjaini+i] = _wedges[i];

		for (i=0;i<_nsets;i++) wverticesloc[i] = nitotarr[i];
		for (i=0;i<aloc.nlist;i++) wverticesloc[_nsets+i] = _wvertices[i];

		nzjaini = 0;
		anew.ia[0] = 0;
		for (i=0;i<_nsets;i++) {
			nzjaini += nzjarr[i];
			anew.ia[i+1] = nzjaini;
		};
		for (i=0;i<aloc.nlist;i++) anew.ia[i+_nsets+1] = anew.ia[i+_nsets]+(aloc.ia[i+1]-aloc.ia[i]);

		anew.m = ntot;
		anew.n = ntot;
		anew.nsupr = ntot;
		anew.nsupc = ntot;
		anew.nlist = nlistloc;
		anew.nzja = nzjaloc;

		aloc = anew;

		delete [] _wedges;
		_wedges = wedgesloc;

		delete [] _wvertices;
		_wvertices = wverticesloc;

	};

// Filter diagonal values

	int nlistloc = aloc.nlist;
	nzloc = 0;

	delete [] ia_cont;

	ia_cont = new int [nlistloc+1];
	if(!ia_cont) MemoryFail(funcname);

	ia_cont[0] = 0;

	for (i=0;i<nlistloc;i++) {
		irow = aloc.list[i];
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jj = aloc.ja[j];
			if (jj != irow) {
				aloc.ja[nzloc] = aloc.ja[j];
				_wedges[nzloc] = _wedges[j];
				nzloc++;
			};
		};
		ia_cont[i+1] = nzloc;
	};

	for (i=0;i<=nlistloc;i++) aloc.ia[i] = ia_cont[i];

//	aloc.nzja = nzjaloc;
	aloc.nzja = nzloc;

// Free work arrays

	for (i=0;i<_nsets;i++) {
		delete [] plistblk[i];
		delete [] pnzblk[i];
	};

	delete [] gblks;
	delete [] nd2blk;
	delete [] blk2pack;
	delete [] ja2cpu;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;
	delete [] listbnd;
	delete [] listbnd2blk;
	delete [] iiarr;
	delete [] nzjarr;
	delete [] plistblk;
	delete [] pnzblk;
	delete [] nsizerowarr;
	delete [] nitotarr;
	delete [] ia_cont;

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Change rows/column indices
//========================================================================================
void CSMatrix::ChangeIndices (CMPIComm &_comm, // Change rows/column indices
										int *_listord) {

	const char *funcname = "ChangeIndices";

	int nproc = _comm.GetNproc ();

// Compute the total size

	int ntot = -1;

	int i, jj;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj > ntot) ntot = jj;
	};

	CMPIExchange::ExchangeArrayMPI (_comm,
												INTEGERVALUE, MAXIMUM, 1, &ntot, &ntot);

	ntot++;

	if (ntot <= 0) return;

// Create the initial matrix partitioning

	int *ibegarr;
	int *iendarr;

	ibegarr = new int [nlist];
	if(!ibegarr) MemoryFail(funcname);
	iendarr = new int [nlist];
	if(!iendarr) MemoryFail(funcname);

	int npartsloc = 0;
	int nz;

	if (nlist > 0) {
		npartsloc = 1;
		ibegarr[npartsloc-1] = list[0];
		iendarr[npartsloc-1] = list[0];
		nz = 0;
		while (nz < nlist-1) {
			nz++;
			if (list[nz] == iendarr[npartsloc-1]+1) {
				iendarr[npartsloc-1] = list[nz];
			} else {
				npartsloc++;
				ibegarr[npartsloc-1] = list[nz];
				iendarr[npartsloc-1] = list[nz];
			};
		};
	};

// Compute the global partitioning

	int nptot;
	int *ibegtot, *iendtot, *icputot;

	CSMatrix::ComputeIntervals (_comm,
											npartsloc, ibegarr, iendarr,
											nptot, ibegtot, iendtot, icputot);

	int *blkstot;

	blkstot = new int [nptot+1];
	if(!blkstot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		blkstot[i] = ibegtot[i];
	};
	blkstot[nptot] = ntot;

	int *ibstot;

	ibstot = new int [nptot];
	if(!ibstot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibstot[i] = -1;

	int iblk = 0;
	int niloc = 0;

	int ibeg;

	for (i=0;i<npartsloc;i++) {
		ibeg = ibegarr[i];
		while (ibegtot[iblk] < ibeg) iblk++;
		ibstot[iblk] = niloc;
		niloc += blkstot[iblk+1]-blkstot[iblk];
	};

	delete [] ibegarr;
	delete [] iendarr;

// Check that the intervals are correct

	if (ibegtot[0] != 0) {
		cout << " throw CSMatrix::ChangeIndices: incorrect intervals " << endl;
		throw " CSMatrix::ChangeIndices: incorrect intervals ";
	};

	for (i=0;i<nptot-1;i++) {
		if (iendtot[i] != ibegtot[i+1]-1) {
			cout << " throw CSMatrix::ChangeIndices: incorrect intervals " << endl;
			throw " CSMatrix::ChangeIndices: incorrect intervals ";
		};
	};

	if (iendtot[nptot-1] != ntot-1) {
		cout << " throw CSMatrix::ChangeIndices: incorrect intervals " << endl;
		throw " CSMatrix::ChangeIndices: incorrect intervals ";
	};

// For each column index get its cpu number

	int *ja2blk;
	int *ja2cpu;

	ja2blk = new int [nzja];
	if(!ja2blk) MemoryFail(funcname);
	ja2cpu = new int [nzja];
	if(!ja2cpu) MemoryFail(funcname);

/*	int iblkbeg, iblkend;

	int iblkprev = 0;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= ibegtot[iblkprev] && jj <= iendtot[iblkprev]) {
			ja2blk[i] = iblkprev;
			ja2cpu[i] = icputot[iblkprev];
		} else {
			iblkbeg = 0;
			iblkend = nptot-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= ibegtot[iblk] && jj <= iendtot[iblk]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < ibegtot[iblk]) {
					iblkend = iblk-1;
				} else if (jj > iendtot[iblk]) {
					iblkbeg = iblk+1;
				};
			};
			iblkprev = iblkbeg;
			ja2blk[i] = iblkprev;
			ja2cpu[i] = icputot[iblkprev];
		};
	};
*/
	int *blkstot1;

	blkstot1 = new int [nptot+1];
	if(!blkstot1) MemoryFail(funcname);

	for (i=0;i<nptot;i++) blkstot1[i] = ibegtot[i];
	blkstot1[nptot] = ntot;

	CSMatrix::ComputeJa2Blk (nptot, blkstot1,
										nlist, ia, ja,
										ja2blk);

	delete [] blkstot1;

	for (i=0;i<nzja;i++) {
		iblk = ja2blk[i];
		ja2cpu[i] = icputot[iblk];
	};

// Prepare send list for each cpu

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if(!iacpu) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int icpu;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int *isendarr;

	isendarr = new int [nzja*3];
	if(!isendarr) MemoryFail(funcname);

	int k;

	for (i=0;i<nzja;i++) {
		icpu = ja2cpu[i];
		k = iptr[icpu];
		isendarr[k*3] = i;
		isendarr[k*3+1] = ja2blk[i];
		isendarr[k*3+2] = ja[i];
		iptr[icpu]++;
	};

// Exchange indices information

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = (iacpu[i+1]-iacpu[i])*3*sizeof(int);
		ObjSend[i] = (char *) (isendarr+3*iacpu[i]);
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) {
		cout << " throw CSMatrix::ChangeIndices: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ChangeIndices: Error in DataExchangeMPI";
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int j, ni, jjloc, jjnew, jjblk, ishift;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (3*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jjblk = piarr[3*j+1];
			jj = piarr[3*j+2];
			ibeg = blkstot[jjblk];
			ishift = ibstot[jjblk];
			jjloc = jj-ibeg+ishift;
			jjnew = _listord[jjloc];
			piarr[3*j+2] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) {
		cout << " throw CSMatrix::ChangeIndices: Error in DataExchangeMPI" << endl;
		throw " CSMatrix::ChangeIndices: Error in DataExchangeMPI";
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Modify the sparsity data

	for (i=0;i<nlist;i++) list[i] = _listord[i];

	int ind;

	for (i=0;i<NObjSend;i++) {
		ni = ObjSizeSend[i] / (3*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ind = piarr[3*j];
			jjnew = piarr[3*j+2];
			ja[ind] = jjnew;
		};
	};

// Free receive arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Sort

	SortListAndColumns ();

// Free work arrays

	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] blkstot;
	delete [] ibstot;
	delete [] ja2blk;
	delete [] ja2cpu;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;

};

// Author: Kharchenko S.A.
// CSMatrix: Modify in parallel the disjoint sets
//========================================================================================
void CSMatrix::ModifyDisjointParts (CMPIComm &_comm, // Modify in parallel the disjoint sets
												int *_blkscpu, 
												int *_ia, int *_ja, int *_wvertices, int *_partition) {

	const char *funcname = "ModifyDisjointParts";

// Compute disjoint sets

	int nsets;
	int *iasets;
	int *jasets;

	ComputeDisjointSets (_comm,
								_blkscpu, 0,
								_ia, _ja, _partition,
								nsets, iasets, jasets);

	cout << " ||cpu, index 0 = " << nsets << endl;

// For all sets compute their weights and modify

	int *isetsize;

	isetsize = new int [nsets];
	if(!isetsize) MemoryFail(funcname);

	int i, j, jj;

	for (i=0;i<nsets;i++) isetsize[i] = 0;

	for (i=0;i<nsets;i++) {
		for (j=iasets[i];j<iasets[i+1];j++) {
			jj = jasets[j];
			isetsize[i] += _wvertices[jj];
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsets, isetsize, isetsize);

	int indmax = 0;
	int isizemax = isetsize[0];

	for (i=1;i<nsets;i++) {
		if (isetsize[i] > isizemax) {
			indmax = i;
			isizemax = isetsize[i];
		};
	};

	for (i=0;i<nsets;i++) {
		if (i != indmax && (double) isetsize[i]*1.0e2 < isizemax) {
			for (j=iasets[i];j<iasets[i+1];j++) {
				jj = jasets[j];
				_partition[jj] = 1;
			};
		};
	};

	delete [] isetsize;
	delete [] iasets;
	delete [] jasets;

// Compute disjoint sets

	ComputeDisjointSets (_comm,
								_blkscpu, 1,
								_ia, _ja, _partition,
								nsets, iasets, jasets);

	cout << " ||cpu, index 1 = " << nsets << endl;

// For all sets compute their weights and modify

	isetsize = new int [nsets];
	if(!isetsize) MemoryFail(funcname);

	for (i=0;i<nsets;i++) isetsize[i] = 0;

	for (i=0;i<nsets;i++) {
		for (j=iasets[i];j<iasets[i+1];j++) {
			jj = jasets[j];
			isetsize[i] += _wvertices[jj];
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsets, isetsize, isetsize);

	indmax = 0;
	isizemax = isetsize[0];

	for (i=1;i<nsets;i++) {
		if (isetsize[i] > isizemax) {
			indmax = i;
			isizemax = isetsize[i];
		};
	};

	for (i=0;i<nsets;i++) {
		if (i != indmax && (double) isetsize[i]*1.0e2 < isizemax) {
			for (j=iasets[i];j<iasets[i+1];j++) {
				jj = jasets[j];
				_partition[jj] = 0;
			};
		};
	};

	delete [] isetsize;
	delete [] iasets;
	delete [] jasets;

};

// Author: Kharchenko S.A.
// CSMatrix: Modify the disjoint sets
//========================================================================================
void CSMatrix::ModifyDisjointParts ( // Modify the disjoint sets
												int _ni, int *_ia, int *_ja, int *_wvertices, int *_partition) {

	const char *funcname = "ModifyDisjointParts";

// Init the sizes

	int ibeg = 0;
	int iend = _ni-1;

	int nsets;
	int *iasets;
	int *jasets;

// Check connectivity on entry
/*
	if (false) {

		int *partloc;

		partloc = new int [_ni];
		if(!partloc) MemoryFail(funcname);

		int i;

		for (i=0;i<_ni;i++) partloc[i] = 0;

		ComputeDisjointSets (ibeg, iend, 0,
									_ia, _ja, partloc,
									nsets, iasets, jasets);

		cout << " Connectivity check: = " << nsets << endl;

		delete [] iasets;
		delete [] jasets;

		delete [] partloc;

	};
*/
// Compute disjoint sets

	ComputeDisjointSets (ibeg, iend, 0,
								_ia, _ja, _partition,
								nsets, iasets, jasets);

	cout << " 1cpu, index 0 = " << nsets << endl;

// For all sets compute their weights and modify

	int *isetsize;

	isetsize = new int [nsets];
	if(!isetsize) MemoryFail(funcname);

	int i, j, jj;

	for (i=0;i<nsets;i++) isetsize[i] = 0;

	for (i=0;i<nsets;i++) {
		for (j=iasets[i];j<iasets[i+1];j++) {
			jj = jasets[j];
			isetsize[i] += _wvertices[jj];
		};
	};

	int indmax = 0;
	int isizemax = isetsize[0];

	for (i=1;i<nsets;i++) {
		if (isetsize[i] > isizemax) {
			indmax = i;
			isizemax = isetsize[i];
		};
	};

	for (i=0;i<nsets;i++) {
		if (i != indmax && (double) isetsize[i]*1.0e2 < isizemax) {
			for (j=iasets[i];j<iasets[i+1];j++) {
				jj = jasets[j];
				_partition[jj] = 1;
			};
		};
	};

	delete [] isetsize;
	delete [] iasets;
	delete [] jasets;

// Compute disjoint sets

	ComputeDisjointSets (ibeg, iend, 1,
								_ia, _ja, _partition,
								nsets, iasets, jasets);

	cout << " 1cpu, index 1 = " << nsets << endl;

// For all sets compute their weights and modify

	isetsize = new int [nsets];
	if(!isetsize) MemoryFail(funcname);

	for (i=0;i<nsets;i++) isetsize[i] = 0;

	for (i=0;i<nsets;i++) {
		for (j=iasets[i];j<iasets[i+1];j++) {
			jj = jasets[j];
			isetsize[i] += _wvertices[jj];
		};
	};

	indmax = 0;
	isizemax = isetsize[0];

	for (i=1;i<nsets;i++) {
		if (isetsize[i] > isizemax) {
			indmax = i;
			isizemax = isetsize[i];
		};
	};

	for (i=0;i<nsets;i++) {
		if (i != indmax && (double) isetsize[i]*1.0e2 < isizemax) {
			for (j=iasets[i];j<iasets[i+1];j++) {
				jj = jasets[j];
				_partition[jj] = 0;
			};
		};
	};

	delete [] isetsize;
	delete [] iasets;
	delete [] jasets;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel the set of disjoint sets for the interval and ipart
//========================================================================================
void CSMatrix::ComputeDisjointSets (CMPIComm &_comm, // Compute in parallel the set of disjoint sets for the interval and ipart
												int *_blkscpu, int _ipart,
												int *_ia, int *_ja, int *_partition,
												int &_nsets, int *&_iasets, int *&_jasets) {

	const char *funcname = "ComputeDisjointSets";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Compute initial disjoint sets

	int ibeg = _blkscpu[myid];
	int iend = _blkscpu[myid+1]-1;
	int niloc = iend-ibeg+1;

	ComputeDisjointSets (ibeg, iend, _ipart,
								_ia, _ja, _partition,
								_nsets, _iasets, _jasets);

// Register each index

	int *irow2set;

	irow2set = new int [niloc];
	if(!irow2set) MemoryFail(funcname);

	int i, j, jj;

	for (i=0;i<niloc;i++) irow2set[i] = -1;

	for (i=0;i<_nsets;i++) {
		for (j=_iasets[i];j<_iasets[i+1];j++) {
			jj = _jasets[j];
			irow2set[jj] = i;
		};
	};

// Count the whole number of sets

	int *iasetscpu;

	iasetscpu = new int [nproc+1];
	if(!iasetscpu) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iasetscpu[i] = 0;

	iasetscpu[myid+1] = _nsets;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, iasetscpu, iasetscpu);

	for (i=0;i<nproc;i++) iasetscpu[i+1] = iasetscpu[i]+iasetscpu[i+1];

	int nsetstot = iasetscpu[nproc];
	int ishift = iasetscpu[myid];

	for (i=0;i<niloc;i++) {
		if (irow2set[i] >= 0) irow2set[i] += ishift;
	};

// Mark each set by its own value

	int *jasetscpu;

	jasetscpu = new int [nsetstot];
	if(!jasetscpu) MemoryFail(funcname);

	for (i=0;i<nsetstot;i++) jasetscpu[i] = i;

// For each index find cpu number to be sent

	int nzjaloc = _ia[niloc];

	int *ja2blk;

	ja2blk = new int [nzjaloc];
	if(!ja2blk) MemoryFail(funcname);

/*	int iblk, ibegblk, iendblk;

	int iblkprev = 0;

	for (i=0;i<nzjaloc;i++) {
		jj = _ja[i];
		if (jj >= _blkscpu[iblkprev] && jj < _blkscpu[iblkprev+1]) {
			iblk = iblkprev;
		} else {
			ibegblk = 0;
			iendblk = nproc-1;
			while (ibegblk != iendblk) {
				iblk = (ibegblk+iendblk) / 2;
				if (jj >= _blkscpu[iblk] && jj < _blkscpu[iblk+1]) {
					ibegblk = iblk;
					iendblk = iblk;
				} else if (jj < _blkscpu[iblk]) {
					iendblk = iblk-1;
				} else if (jj >= _blkscpu[iblk+1]) {
					ibegblk = iblk+1;
				};
			};
			iblk = ibegblk;
		};
		ja2blk[i] = iblk;
		iblkprev = iblk;
	};
*/
	CSMatrix::ComputeJa2Blk (nproc, _blkscpu,
										niloc, _ia, _ja,
										ja2blk);

	int *iacpu;
	int *iptrcpu;

	iacpu = new int [nproc+1];
	if(!iacpu) MemoryFail(funcname);
	iptrcpu = new int [nproc];
	if(!iptrcpu) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iproc;

	for (i=0;i<nzjaloc;i++) {
		iproc = ja2blk[i];
		iacpu[iproc+1]++;
	};
	iacpu[myid+1] = 0;

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	int nzjasend = iacpu[nproc];

	int ldsend = 3;

	int *sendarr;

	sendarr = new int [nzjasend*ldsend];
	if(!sendarr) MemoryFail(funcname);

	int *ialists;
	int *jalists;
	int *ialists1;
	int *jasetscpu1;

	ialists = new int [nsetstot+1];
	if(!ialists) MemoryFail(funcname);
	jalists = new int [nsetstot];
	if(!jalists) MemoryFail(funcname);
	ialists1 = new int [nsetstot+1];
	if(!ialists1) MemoryFail(funcname);
	jasetscpu1 = new int [nsetstot];
	if(!jasetscpu1) MemoryFail(funcname);

// Perform all computations in parallel

	int ivalue, jvalue, k, iset, ni, jblk, jjloc, jset, jvalmax, kk;

	while (true) {

// Create current registration lists

		for (i=0;i<=nsetstot;i++) ialists[i] = 0;

		for (i=0;i<nsetstot;i++) {
			ivalue = jasetscpu[i];
			ialists[ivalue+1]++;
		};

		for (i=0;i<nsetstot;i++) ialists[i+1] = ialists[i]+ialists[i+1];

		for (i=0;i<nsetstot;i++) ialists1[i] = ialists[i];

		for (i=0;i<nsetstot;i++) {
			ivalue = jasetscpu[i];
			k = ialists1[ivalue];
			jalists[k] = i;
			ialists1[ivalue]++;
		};

		for (i=0;i<nsetstot;i++) jasetscpu1[i] = jasetscpu[i];

// Prepare current bordering data

		for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

		for (i=0;i<niloc;i++) {
			iset = irow2set[i];
			if (iset >= 0) {
				ivalue = jasetscpu[iset];
			} else {
				ivalue = -1;
			};
			for (j=_ia[i];j<_ia[i+1];j++) {
				jj = _ja[j];
				iproc = ja2blk[j];
				if (iproc != myid) {
					k = iptrcpu[iproc];
					sendarr[k*ldsend] = _ja[j];
					sendarr[k*ldsend+1] = iproc;
					sendarr[k*ldsend+2] = ivalue;
					iptrcpu[iproc]++;
				};
			};
		};

// Exchange current bordering data

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

		for(i = 0; i < NObjSend; i++) {
			ObjTypeSend[i] = 1;
			ObjIDSend[i] = i;
			CpuIDSend[i] = i;
			ObjSizeSend[i] = ldsend*(iacpu[i+1]-iacpu[i])*sizeof(int);
			ObjSend[i] = (char *) (sendarr+ldsend*iacpu[i]);
		};

		int NObjRecv;
		int* ObjTypeRecv;
		int* ObjIDRecv;
		int* CpuIDRecv;
		int* ObjSizeRecv;
		char** ObjRecv;

		int info;

		info = CMPIExchange::DataExchangeMPI (_comm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
		if(info) {
			cout << " throw CSMatrix::SplitBlockMatrix: Error in DataExchangeMPI" << endl;
			throw " CSMatrix::SplitBlockMatrix: Error in DataExchangeMPI";
		};

		delete [] ObjTypeSend;
		delete [] ObjIDSend;
		delete [] CpuIDSend;
		delete [] ObjSizeSend;
		delete [] ObjSend;

// Perform reregistration

		int *piarr;

		for (i=0;i<NObjRecv;i++) {
			ni = ObjSizeRecv[i] / (ldsend*sizeof(int));
			piarr = (int *) ObjRecv[i];
			for (j=0;j<ni;j++) {
				jj = piarr[ldsend*j];
				jblk = piarr[ldsend*j+1];
				ivalue = piarr[ldsend*j+2];
				jjloc = jj-ibeg;
				jset = irow2set[jjloc];
				if (jset >= 0) {
					jvalue = jasetscpu1[jset];
				} else {
					jvalue = -1;
				};
				if (ivalue >= 0 && jvalue >= 0) {
					if (ivalue != jvalue) {
						jvalmax = ivalue;
						if (jvalmax < jvalue) jvalmax = jvalue;
						for (k=ialists[ivalue];k<ialists[ivalue+1];k++) {
							kk = jalists[k];
							jasetscpu1[kk] = jvalmax;
						};
						for (k=ialists[jvalue];k<ialists[jvalue+1];k++) {
							kk = jalists[k];
							jasetscpu1[kk] = jvalmax;
						};
					};
				};
			};
		};

		delete [] ObjTypeRecv;
		delete [] ObjIDRecv;
		delete [] CpuIDRecv;
		delete [] ObjSizeRecv;
		for(i = 0; i < NObjRecv; i++) {
			delete [] ObjRecv[i];
		};
		delete [] ObjRecv;

// Compare

		CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, nsetstot, jasetscpu1, jasetscpu1);

		int icheck = 1;

		for (i=0;i<nsetstot;i++) {
			if (jasetscpu1[i] != jasetscpu[i]) icheck = 0;
		};

		if (icheck == 1) break;

		for (i=0;i<nsetstot;i++) jasetscpu[i] = jasetscpu1[i];

	};

// Restructure the sets according to global information

	for (i=0;i<=nsetstot;i++) ialists[i] = 0;

	for (i=0;i<nsetstot;i++) {
		ivalue = jasetscpu[i];
		ialists[ivalue+1]++;
	};

	for (i=0;i<nsetstot;i++) ialists[i+1] = ialists[i]+ialists[i+1];

	for (i=0;i<nsetstot;i++) ialists1[i] = ialists[i];

	for (i=0;i<nsetstot;i++) {
		ivalue = jasetscpu[i];
		k = ialists1[ivalue];
		jalists[k] = i;
		ialists1[ivalue]++;
	};

	int nsetsnew = 0;

	for (i=0;i<nsetstot;i++) {
		if (ialists[i+1] != ialists[i]) nsetsnew++;
	};

	int *iasetsnew;

	iasetsnew = new int [nsetsnew+1];
	if(!iasetsnew) MemoryFail(funcname);

	int *jasetsnew;

	jasetsnew = new int [niloc];
	if(!jasetsnew) MemoryFail(funcname);

	int nz = 0;
	nsetsnew = 0;
	iasetsnew[0] = 0;

	int jsetloc;

	for (i=0;i<nsetstot;i++) {
		if (ialists[i+1] != ialists[i]) {
			for (j=ialists[i];j<ialists[i+1];j++) {
				jset = jalists[j];
				if (jset >= iasetscpu[myid] && jset < iasetscpu[myid+1]) {
					jsetloc = jset-iasetscpu[myid];
					for (k=_iasets[jsetloc];k<_iasets[jsetloc+1];k++) {
						jasetsnew[nz] = _jasets[k];
						nz++;
					};
				};
			};
			iasetsnew[nsetsnew+1] = nz;
			nsetsnew++;
		};
	};

	for (i=0;i<niloc;i++) _jasets[i] = jasetsnew[i];

	delete [] _iasets;

	_nsets = nsetsnew;
	_iasets = iasetsnew;

// Free work arrays

	delete [] irow2set;
	delete [] iasetscpu;
	delete [] jasetscpu;
	delete [] ja2blk;
	delete [] iacpu;
	delete [] iptrcpu;
	delete [] sendarr;
	delete [] ialists;
	delete [] jalists;
	delete [] ialists1;
	delete [] jasetscpu1;
	delete [] jasetsnew;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute the set of dijoint sets for the interval and ipart
//========================================================================================
void CSMatrix::ComputeDisjointSets (int _ibeg, int _iend, int _ipart, // Compute the set of dijoint sets for the interval and ipart
												int *_ia, int *_ja, int *_partition,
												int &_nsets, int *&_iasets, int *&_jasets) {

	const char *funcname = "ComputeDisjointSets";

// Init the mask array

	int niloc = _iend-_ibeg+1;

	int *imask;

	imask = new int [niloc];
	if(!imask) MemoryFail(funcname);

	int i;

	for (i=0;i<niloc;i++) imask[i] = 0;

// Main cycle over the nodes

	_iasets = new int [niloc+1];
	if(!_iasets) MemoryFail(funcname);
	_jasets = new int [niloc];
	if(!_jasets) MemoryFail(funcname);

	_nsets = 0;
	_iasets[0] = 0;

	int nz0, nz, nz1, ilist, j, jj, jlist;

	nz = 0;

	while (true) {

// Find initial search node

		nz0 = nz;

		for (i=0;i<niloc;i++) {
			if (imask[i] == 0 && _partition[i] == _ipart) {
				_nsets++;
				imask[i] = _nsets;
				_jasets[nz] = i;
				nz++;
				break;
			};
		};

		if (nz == nz0) break;

		while (nz != nz0) {
			nz1 = nz;
			for (i=nz0;i<nz;i++) {
				ilist = _jasets[i];
				for (j=_ia[ilist];j<_ia[ilist+1];j++) {
					jj = _ja[j];
					if (jj >= _ibeg && jj <= _iend) {
						jlist = jj-_ibeg;
						if (imask[jlist] == 0 && _partition[jlist] == _ipart) {
							_jasets[nz1] = jlist;
							nz1++;
							imask[jlist] = _nsets;
						};
					};
				};
			};
			nz0 = nz;
			nz = nz1;
		};

		_iasets[_nsets] = nz;

	};

// Free work arrays

	delete [] imask;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
