//------------------------------------------------------------------------------------------------
// File: DecompBnd.cpp
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
#include "gsmatrix.h"
#include "smatrix.h"
#include "SolverPar.h"

using namespace std;


// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel decomposition into the set of volume and boundary blocks
//========================================================================================
void CSMatrix::DecompSchurBndRecursive (CMPIComm &_comm, int _decomptype, int _ncycle, int _ncolors, // Compute in parallel decomposition into the set of volume and boundary blocks
														int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
														int *_order,
														int &_nblksnew, int *&_blknew2ext, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const {

	const char *funcname = "DecompSchurBndRecursive_00";

//	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Perform major surface partitioning

	DecompSurfaceRecursive (_comm, _decomptype, _ncolors,
										_nblks, _blks, _blk2cpu, _blk2color,
										_order,
										_nblksnew, _blksnew, _blk2cpunew, _blk2colornew);

	_blknew2ext = new int [_nblksnew+1];
	if (!_blknew2ext) MemoryFail (funcname);

	int i;

	for (i=0;i<=_nblksnew;i++) _blknew2ext[i] = i;

// Perform further partitioning if necessary

	if (_decomptype > 0) {

// Reassing previous data

		int nblksnew, *blknew2ext, *blksnew, *blk2cpunew, *blk2colornew;

		nblksnew = _nblksnew;

		blknew2ext = _blknew2ext;
		blksnew = _blksnew;
		blk2cpunew = _blk2cpunew;
		blk2colornew = _blk2colornew;

// Transform the matrix

		CSMatrix mtrloc;

		mtrloc = *this;

		mtrloc.ChangeIndices (_comm, _order);

		int *ibegarrloc;
		int *iendarrloc;

		ibegarrloc = new int [nblksnew];
		if (!ibegarrloc) MemoryFail (funcname);
		iendarrloc = new int [nblksnew];
		if (!iendarrloc) MemoryFail (funcname);

		int nploc = 0;

		for (i=0;i<nblksnew;i++) {
			if (blk2cpunew[i] == myid && blksnew[i] < blksnew[i+1]) {
				ibegarrloc[nploc] = blksnew[i];
				iendarrloc[nploc] = blksnew[i+1]-1;
				nploc++;
			};
		};

		CSMatrix mtrnew;

		mtrnew = mtrloc.GetSubmatrixViaIntervals (_comm,
																nploc, ibegarrloc, iendarrloc);

		CSMatrix mtrdummy;

		mtrloc = mtrdummy;

// Compute new ordering and partitioning

		int nlist1loc = mtrnew.GetNlist ();

		int *orderloc;

		orderloc = new int [nlist1loc];
		if (!orderloc) MemoryFail (funcname);

		mtrnew.DecompSchurSurface (_comm, _decomptype, _ncycle, _ncolors,
											nblksnew, blksnew, blk2cpunew, blk2colornew,
											orderloc,
											_nblksnew, _blknew2ext, _blksnew, _blk2cpunew, _blk2colornew);

// Split major first blocks for each color

		int nsplit = 1;

		if (nsplit > 1) {

			int nblkstot = _blknew2ext[_nblksnew];

			int *iblkfirst;
			int *imaskblk;

			iblkfirst = new int [_ncolors];
			if (!iblkfirst) MemoryFail (funcname);
			imaskblk = new int [nblkstot];
			if (!imaskblk) MemoryFail (funcname);

			for (i=0;i<_ncolors;i++) iblkfirst[i] = -1;
			for (i=0;i<nblkstot;i++) imaskblk[i] = -1;

			int icolor;

			for (i=0;i<nblkstot;i++) {
				icolor = _blk2colornew[i];
				if (iblkfirst[icolor] == -1) {
					iblkfirst[icolor] = i;
					imaskblk[i] = 0;
				};
			};

			for (i=0;i<_ncolors;i++) {
				if (iblkfirst[i] == -1) throw " CSMatrix::DecompSchurBndRecursive: some colors have no block ";
			};

			int nblkstotext = nblkstot + (nsplit-1)*_ncolors;

			int *gblknew2ext;
			int *blksnewext;
			int *blk2cpunewext;
			int *blk2colornewext;

			gblknew2ext = new int [_nblksnew+1];
			if (!gblknew2ext) MemoryFail (funcname);
			blksnewext = new int [nblkstotext+1];
			if (!blksnewext) MemoryFail (funcname);
			blk2cpunewext = new int [nblkstotext];
			if (!blk2cpunewext) MemoryFail (funcname);
			blk2colornewext = new int [nblkstotext];
			if (!blk2colornewext) MemoryFail (funcname);

			int nblkloc = 0;

			gblknew2ext[0] = 0;
			blksnewext[0] = 0;

			int j, niloc, nilocsplit, k;

			for (i=0;i<_nblksnew;i++) {
				for (j=_blknew2ext[i];j<_blknew2ext[i+1];j++) {
					if (imaskblk[j] == 0) {
						niloc = _blksnew[j+1]-_blksnew[j];
						nilocsplit = niloc / nsplit;
						for (k=0;k<nsplit;k++) {
							blksnewext[nblkloc+1] = nilocsplit;
							blk2cpunewext[nblkloc] = _blk2cpunew[j];
							blk2colornewext[nblkloc] = _blk2colornew[j];
							nblkloc++;
						};
						blksnewext[nblkloc] = niloc-(nsplit-1)*nilocsplit;
					} else {
						blksnewext[nblkloc+1] = _blksnew[j+1]-_blksnew[j];
						blk2cpunewext[nblkloc] = _blk2cpunew[j];
						blk2colornewext[nblkloc] = _blk2colornew[j];
						nblkloc++;
					};
				};
				gblknew2ext[i+1] = nblkloc;
			};

			for (i=0;i<nblkloc;i++) blksnewext[i+1] = blksnewext[i]+blksnewext[i+1];

			delete [] _blknew2ext;
			delete [] _blksnew;
			delete [] _blk2cpunew;
			delete [] _blk2colornew;

			_blknew2ext = gblknew2ext;
			_blksnew = blksnewext;
			_blk2cpunew = blk2cpunewext;
			_blk2colornew = blk2colornewext;

			delete [] iblkfirst;
			delete [] imaskblk;

		};

// Transform the ordering array

// Compute the direct ordering for the local partitioning

		int nlistord = 0;

		for (i=0;i<nblksnew;i++) {
			if (blk2cpunew[i] == myid) {
				nlistord += (blksnew[i+1]-blksnew[i]);
			};
		};

		int *listord;

		listord = new int [nlistord];
		if (!listord) MemoryFail (funcname);

		nlistord = 0;

		int j;

		for (i=0;i<nblksnew;i++) {
			if (blk2cpunew[i] == myid) {
				for (j=blksnew[i];j<blksnew[i+1];j++) {
					listord[nlistord] = j;
					nlistord++;
				};
			};
		};

		int nlistpart;
		int *orderpart;

		CSMatrix::DirOrderViaIntervals (_comm,
													nlistord, listord, orderloc,
													nploc, ibegarrloc, iendarrloc, 
													nlistpart, orderpart);

		int *order;

		order = new int [nlist];
		if (!order) MemoryFail (funcname);

		for (i=0;i<nlist;i++) order[i] = _order[i];

		CSMatrix::OrderViaIntervals (_comm,
												nlist, order,
												nploc, ibegarrloc, iendarrloc, 
												orderpart,
												_order);

		delete [] ibegarrloc;
		delete [] iendarrloc;

		delete [] listord;
		delete [] orderpart;

// Free work arrays

		delete [] blknew2ext;
		delete [] blksnew;
		delete [] blk2cpunew;
		delete [] blk2colornew;
		delete [] orderloc;
		delete [] order;

	};

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel the surfaces that separate processors data
//========================================================================================
void CSMatrix::DecompSurfaceRecursive (CMPIComm &_comm, int _decomptype, int _ncolors, // Compute in parallel the surfaces that separate processors data
														int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
														int *_order,
														int &_nblksnew, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const {

	const char *funcname = "DecompSurfaceRecursive";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that the number of colors on entry is the number of blocks

	if (_ncolors != _nblks) throw " CSMatrix::DecompSchurBndRecursive: the number of colors is not the number of blocks ";

// Check that the number of colors and cpu's are powers of 2

	int ncolors_pwr2 = 1;
	while (ncolors_pwr2*2 <= _ncolors) ncolors_pwr2 *= 2;

	if (ncolors_pwr2 != _ncolors) throw " CSMatrix::DecompSchurBndRecursive: the number of colors is not a power of 2 ";

	int nproc_pwr2 = 1;
	while (nproc_pwr2*2 <= nproc) nproc_pwr2 *= 2;

	if (nproc_pwr2 != nproc) throw " CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is not a power of 2 ";

	if (nproc > _ncolors) throw " CSMatrix::SplitBlockMatrixRecursive: the number of cpu's is larger than the number of colors ";

	int nlev = 0;

	int np = 1;

	while (np*2 <= nproc) {
		np *= 2;
		nlev++;
	};

// Create the color2cpu and back correspondence

	int *color2cpu;
	int *iacpu2color, *jacpu2color;

	color2cpu = new int [_ncolors];
	if (!color2cpu) MemoryFail (funcname);
	iacpu2color = new int [nproc+1];
	if (!iacpu2color) MemoryFail (funcname);
	jacpu2color = new int [_ncolors];
	if (!jacpu2color) MemoryFail (funcname);

	int ni = _ncolors / nproc;

	int icolor = 0;
	iacpu2color[0] = 0;

	int i, j;

	for (i=0;i<nproc;i++) {
		for (j=0;j<ni;j++) {
			color2cpu[icolor] = i;
			jacpu2color[icolor] = icolor;
			icolor++;
		};
		iacpu2color[i+1] = icolor;
	};

// Check the matrix data on entry

	int nlistloc = 0;

	for (i=0;i<_nblks;i++) {
		if (_blk2cpu[i] == myid) {
			if (nlist < nlistloc) throw " CSMatrix::DecompSchurBndRecursive: wrong part of the matrix is considered";
			if (nlist < nlistloc+(_blks[i+1]-_blks[i])) throw " CSMatrix::DecompSchurBndRecursive: wrong part of the matrix is considered";
			if (_blks[i] != _blks[i+1]) {
				int indbeg = nlistloc;
				int indend = nlistloc+(_blks[i+1]-_blks[i])-1;
				int irowbeg = list[indbeg];
				int irowend = list[indend];
				if (irowbeg != _blks[i] || irowend != _blks[i+1]-1) throw " CSMatrix::DecompSchurBndRecursive: wrong part of the matrix is considered";
			};
			nlistloc += (_blks[i+1]-_blks[i]);
		};
	};

	if (nlist != nlistloc) throw " CSMatrix::DecompSchurBndRecursive: wrong part of the matrix is considered";

// Allocate the memory that supports the levels computations

	int **pporderlev;
	int *nblkslev, **ppblkslev, **ppblk2cpu, **pppartblk;
	int *nblkslevprep, **ppblkslevprep, **ppblk2cpuprep, **pppartblkprep;

	pporderlev = new int * [nlev+2];
	if (!pporderlev) MemoryFail (funcname);
	nblkslev = new int [nlev+2];
	if (!nblkslev) MemoryFail (funcname);
	ppblkslev = new int * [nlev+2];
	if (!ppblkslev) MemoryFail (funcname);
	ppblk2cpu = new int * [nlev+2];
	if (!ppblk2cpu) MemoryFail (funcname);
	pppartblk = new int * [nlev+2];
	if (!pppartblk) MemoryFail (funcname);
	nblkslevprep = new int [nlev+2];
	if (!nblkslevprep) MemoryFail (funcname);
	ppblkslevprep = new int * [nlev+2];
	if (!ppblkslevprep) MemoryFail (funcname);
	ppblk2cpuprep = new int * [nlev+2];
	if (!ppblk2cpuprep) MemoryFail (funcname);
	pppartblkprep = new int * [nlev+2];
	if (!pppartblkprep) MemoryFail (funcname);

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

	int nblkmax = 4*nlev+9;

	int *blkstot;

	blkstot = new int [nblkmax+1];
	if (!blkstot) MemoryFail (funcname);

	for (i=0;i<=nblkmax;i++) blkstot[i] = 0;

	int *blksloc;

	blksloc = new int [_nblks+1];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<_nblks+1;i++) blksloc[i] = _blks[i];

	nblkslevprep[0] = _nblks;
	ppblkslevprep[0] = blksloc;

	blksloc = new int [_nblks];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) blksloc[i] = _blk2cpu[i];

	ppblk2cpuprep[0] = blksloc;

	blksloc = new int [_nblks];
	if (!blksloc) MemoryFail (funcname);

	int nblks1 = _nblks / 2;
	if (nblks1 > 0) {
		for (i=0;i<_nblks;i++) blksloc[i] = i / nblks1;
	} else {
		for (i=0;i<_nblks;i++) blksloc[i] = 0;
	};

	pppartblkprep[0] = blksloc;

// Find first and last block indices

	int iblkbeg = -1;
	int iblkend = -1;

	for (i=0;i<_nblks;i++) {
		if (_blk2cpu[i] == myid && iblkbeg == -1) iblkbeg = i;
	};
	for (i=_nblks-1;i>=0;i--) {
		if (_blk2cpu[i] == myid && iblkend == -1) iblkend = i;
	};
	for (i=iblkbeg;i<=iblkend;i++) {
		if (_blk2cpu[i] != myid) throw " CSMatrix::DecompSchurBndRecursive: wrong part of the matrix is considered";
	};

	int ibegblkloc = iblkbeg;
	int iendblkloc = iblkend;

	int nblksloc = _nblks;

// Main cycle over the levels of the tree

	CSMatrix mtrlev;
	CSMatrix mtrdummy;

	mtrlev = *this;

	int igroup, myid1;
	int nz, jj, nloc;
	int *plist, *pia, *pja, *ialoc;

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkDecompBnd_",myid,".dat");
//	ofstream ffffout (strbuff);

	for (ilev=0;ilev<nlev;ilev++) {

// Compute local ordering and partitioning, and reorder the matrix indices

//		ffffout << " Ilev = " << ilev << " Matr = " << mtrlev << endl;
//		OutArr (ffffout," blks ",nblkslevprep[ilev]+1,ppblkslevprep[ilev]);
//		OutArr (ffffout," blk2cpu ",nblkslevprep[ilev],ppblk2cpuprep[ilev]);
//		OutArr (ffffout," partblk ",nblkslevprep[ilev],pppartblkprep[ilev]);

		mtrlev.ComputeGroupsBoundary (commlev[ilev],
												nblkslevprep[ilev], ppblkslevprep[ilev], ppblk2cpuprep[ilev], pppartblkprep[ilev],
												pporderlev[ilev], 
												nblkslev[ilev+1], ppblkslev[ilev+1], ppblk2cpu[ilev+1], pppartblk[ilev+1]);

// Reorder the matrix in-place

		mtrlev.ChangeIndices (commlev[ilev], pporderlev[ilev]);

// Create modified block partitioning

		nprocloc = commlev[ilev].GetNproc ();
		myidloc = commlev[ilev].GetMyid ();

		nproc1 = nprocloc / 2;

		if (myidloc < nproc1) {
			igroup = 0;
		} else {
			igroup = 1;
		};

		myid1 = commlev[ilev+1].GetMyid ();

		int nblks1loc = nblksloc / 2;

		blkstot[ilev*4+0] = 0;
		blkstot[ilev*4+1] = ppblkslev[ilev+1][nblks1loc];
		blkstot[ilev*4+2] = ppblkslev[ilev+1][2*nblks1loc];
		blkstot[ilev*4+3] = ppblkslev[ilev+1][nblkslev[ilev+1]];

		int ishiftblk = 0;
		if (igroup == 1) ishiftblk = nblks1loc;

		int ibegblk1loc = ibegblkloc - ishiftblk;
		int iendblk1loc = iendblkloc - ishiftblk;

		nblkslevprep[ilev+1] = nblks1loc;

		ppblkslevprep[ilev+1] = new int [nblks1loc+1];
		if (!ppblkslevprep[ilev+1]) MemoryFail (funcname);

		ppblkslevprep[ilev+1][0] = 0;

		int iold;

		for (i=0;i<nblks1loc;i++) {
			iold = i+ishiftblk;
			ppblkslevprep[ilev+1][i+1] = ppblkslevprep[ilev+1][i] + (ppblkslev[ilev+1][iold+1]-ppblkslev[ilev+1][iold]);
		};

		ppblk2cpuprep[ilev+1] = new int [nblks1loc];
		if (!ppblk2cpuprep[ilev+1]) MemoryFail (funcname);

		int niblk = iblkend-iblkbeg+1;

		int icpu;

		for (icpu=0;icpu<nproc1;icpu++) {
			for (j=0;j<niblk;j++) {
				ppblk2cpuprep[ilev+1][icpu*niblk+j] = icpu;
			};
		};

		int nblks2loc = nblks1loc / 2;
		if (nblks2loc == 0) nblks2loc = 1;

		pppartblkprep[ilev+1] = new int [nblks1loc];
		if (!pppartblkprep[ilev+1]) MemoryFail (funcname);

		for (j=0;j<nblks1loc;j++) pppartblkprep[ilev+1][j] = j / nblks2loc;

// Get the new submatrix according to the intervals

		CSMatrix mtrnew;

		int *ibegarrloc;
		int *iendarrloc;

		ibegarrloc = new int [nblkslev[ilev+1]];
		if (!ibegarrloc) MemoryFail (funcname);
		iendarrloc = new int [nblkslev[ilev+1]];
		if (!iendarrloc) MemoryFail (funcname);

		int nploc = 0;

		for (i=0;i<nblkslev[ilev+1];i++) {
			if (ppblk2cpu[ilev+1][i] == myidloc && ppblkslev[ilev+1][i] < ppblkslev[ilev+1][i+1])  {
				ibegarrloc[nploc] = ppblkslev[ilev+1][i];
				iendarrloc[nploc] = ppblkslev[ilev+1][i+1]-1;
				nploc++;
			};
		};

		mtrnew = mtrlev.GetSubmatrixViaIntervals (commlev[ilev],
																nploc, ibegarrloc, iendarrloc);

//		OutArr (ffffout," Blocks parti ",nblkslev[ilev+1]+1,ppblkslev[ilev+1]);
//		ffffout << " Mtrnew after get submatrix " << mtrnew << endl;

//		sprintf (strbuff,"%s%i%s","OrdMtr_",myid,".ps");
//		mtrnew.A2Ps (strbuff,nblkslev[ilev+1],ppblkslev[ilev+1]);

		delete [] ibegarrloc;
		delete [] iendarrloc;

// Filter and shift the rows/columns indices

		int ibegtot = ppblkslev[ilev+1][ishiftblk];
		int iendtot = ppblkslev[ilev+1][ishiftblk+nblks1loc]-1;

		int iendloc = ppblkslev[ilev+1][iendblkloc+1]-1;

		nlistloc = mtrnew.GetNlist ();
		plist = mtrnew.GetList ();
		pia = mtrnew.GetIa ();
		pja = mtrnew.GetJa ();

		ialoc = new int [nlistloc+1];
		if (!ialoc) MemoryFail (funcname);

		int nlistnew = 0;
		ialoc[0] = 0;
		nz = 0;

		for (i=0;i<nlistloc;i++) {
			if (plist[i] <= iendloc) {
				plist[nlistnew] = plist[i] - ibegtot;
				for (j=pia[i];j<pia[i+1];j++) {
					jj = pja[j];
					if (jj >= ibegtot && jj <= iendtot) {
						pja[nz] = jj-ibegtot;
						nz++;
					};
				};
				ialoc[nlistnew+1] = nz;
				nlistnew++;
			};
		};

		for (i=0;i<=nlistnew;i++) pia[i] = ialoc[i];

		nloc = ppblkslevprep[ilev+1][nblkslevprep[ilev+1]];

		mtrnew.SetM (nloc);
		mtrnew.SetN (nloc);
		mtrnew.SetNsupr (nloc);
		mtrnew.SetNsupcBase (nloc);
		mtrnew.SetNlist (nlistnew);
		mtrnew.SetNzja (nz);

		delete [] ialoc;

		mtrlev = mtrnew;

		mtrnew = mtrdummy;

// Reassign local block index data

		nblksloc = nblks1loc;

		ibegblkloc = ibegblk1loc;
		iendblkloc = iendblk1loc;

	};

// Reorder and partition the last submatrix

	mtrlev.DecompSurfaceRecursive1Cpu (commlev[nlev], _decomptype,
														nblkslevprep[nlev], ppblkslevprep[nlev], 
														pporderlev[nlev],
														nblkslev[nlev+1], ppblkslev[nlev+1], ppblk2cpu[nlev+1]);

// Modify blocks partitionings for ordering computations

	for (ilev=0;ilev<=nlev;ilev++) {
		nblksloc = 0;
		ppblkslevprep[ilev][0] = 0;
		int iproc = -1;
		for (j=0;j<nblkslevprep[ilev];j++) {
			int jproc = ppblk2cpuprep[ilev][j];
			if (jproc != iproc) {
				ppblkslevprep[ilev][nblksloc+1] = ppblkslevprep[ilev][j+1];
				nblksloc++;
				iproc = jproc;
			} else {
				ppblkslevprep[ilev][nblksloc] = ppblkslevprep[ilev][j+1];
			};
		};
		nblkslevprep[ilev] = nblksloc;
	};

// Backward cycle with the ordering restoration

	if (nproc > 1) {

		for (ilev=nlev-1;ilev>=0;ilev--) {

			CSMatrix::RestoreOrderNd (ilev, nlev,
												commlev, 
												blkstot, pporderlev, nblkslevprep, ppblkslevprep);

		};

	};

	for (i=0;i<nlist;i++) _order[i] = pporderlev[0][i];

// Compute the final arrays

	_nblksnew = 2*_ncolors-1;

	_blksnew = new int [_nblksnew+1];
	if (!_blksnew) MemoryFail (funcname);
	_blk2cpunew = new int [_nblksnew];
	if (!_blk2cpunew) MemoryFail (funcname);
	_blk2colornew = new int [_nblksnew];
	if (!_blk2colornew) MemoryFail (funcname);

	int *imaskblk;
	int *blk2lev;

	imaskblk = new int [_nblksnew];
	if (!imaskblk) MemoryFail (funcname);
	blk2lev = new int [_nblksnew];
	if (!blk2lev) MemoryFail (funcname);

	for (i=0;i<_nblksnew;i++) imaskblk[i] = -1;
	for (i=0;i<_nblksnew;i++) blk2lev[i] = -1;

	double *perf;
	double *memory;

	perf = new double [_ncolors];
	if(!perf) MemoryFail(funcname);
	memory = new double [_ncolors];
	if(!memory) MemoryFail(funcname);

	for (i=0;i<_ncolors;i++) perf[i] = 1.0;
	for (i=0;i<_ncolors;i++) memory[i] = 1.0;

	CTree treebinary (_ncolors, 2, perf, memory);

	treebinary.SortNodes ();

	delete [] perf;
	delete [] memory;

	int ntails = _nblksnew-_ncolors;
	int nlevloc = treebinary.GetNlev ();
	int nnodes = treebinary.GetNnodes ();
	CNode *pnodes = treebinary.GetNode ();

	nz = 0;

	for (ilev=0;ilev<nlevloc;ilev++) {
		if (nz < ntails) {
			for (j=0;j<nnodes;j++) {
				if (pnodes[j].GetNodelv() == ilev) {
					imaskblk[j] = nz;
					_blk2colornew[j] = pnodes[j].GetNodecpu();
					blk2lev[j] = ilev;
					nz++;
				};
			};
		};
	};

	int iblk = 0;

	for (i=0;i<_nblksnew;i++) {
		if (imaskblk[i] == -1) {
			if (iblk < _ncolors) {
				_blk2colornew[i] = iblk;
				blk2lev[i] = nlev;
				iblk++;
			};
		};
	};

	for (i=0;i<_nblksnew;i++) imaskblk[i] = -1;

	int niproc = _ncolors / nproc;

	for (i=0;i<_nblksnew;i++) {
		int icpu = _blk2colornew[i];
		int jcpu = icpu / niproc;
		_blk2cpunew[i] = jcpu;
	};

	for (i=0;i<=_nblksnew;i++) _blksnew[i] = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		for (i=0;i<_nblksnew;i++) {
			if (blk2lev[i] == ilev && _blk2cpunew[i] == myid) {
				_blksnew[i+1] = blkstot[ilev*4+3]-blkstot[ilev*4+2];
				break;
			};
		};
	};

	for (i=0;i<_nblksnew;i++) {
		if (blk2lev[i] == nlev && _blk2cpunew[i] == myid) {
			for (j=0;j<nblkslev[nlev+1];j++) {
				_blksnew[i+j+1] = ppblkslev[nlev+1][j+1]-ppblkslev[nlev+1][j];
			};
			break;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblksnew+1, _blksnew, _blksnew);

	for (i=0;i<_nblksnew;i++) _blksnew[i+1] = _blksnew[i]+_blksnew[i+1];

// Free work arrays

	for (i=0;i<=nlev;i++) {
		delete [] pporderlev[i];
		delete [] ppblkslev[i+1];
		delete [] ppblk2cpu[i+1];
		delete [] ppblkslevprep[i];
		delete [] ppblk2cpuprep[i];
		delete [] pppartblkprep[i];
	};
	for (i=0;i<nlev;i++) {
		delete [] pppartblk[i+1];
	};

	delete [] color2cpu;
	delete [] iacpu2color;
	delete [] jacpu2color;
	delete [] pporderlev;
	delete [] nblkslev;
	delete [] ppblkslev;
	delete [] ppblk2cpu;
	delete [] pppartblk;
	delete [] nblkslevprep;
	delete [] ppblkslevprep;
	delete [] ppblk2cpuprep;
	delete [] pppartblkprep;
	delete [] commlev;
	delete [] listcpu;
	delete [] blkstot;
	delete [] imaskblk;
	delete [] blk2lev;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel decomposition into the set of boundary blocks
//========================================================================================
void CSMatrix::DecompSchurSurface (CMPIComm &_comm, int _decomptype, int _ncycle, int _ncolors, // Compute in parallel decomposition into the set of boundary blocks
												int _nblks, int *_blks, int *_blk2cpu, int *_blk2color,
												int *_order,
												int &_nblksnew, int *&_blknew2ext, int *&_blksnew, int *&_blk2cpunew, int *&_blk2colornew) const {

	const char *funcname = "DecompSchurSurface";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Detect the major valume blocks

	int *imaskblk;
	int *imaskcolor;

	imaskblk = new int [_nblks];
	if(!imaskblk) MemoryFail(funcname);
	imaskcolor = new int [_ncolors];
	if(!imaskcolor) MemoryFail(funcname);

	int i, icolor, iblk;

	for (i=0;i<_nblks;i++) imaskblk[i] = 0;
	for (i=0;i<_ncolors;i++) imaskcolor[i] = -1;

	for (i=0;i<_nblks;i++) {
		icolor = _blk2color[i];
		if (imaskcolor[icolor] == -1) {
			imaskcolor[icolor] = 0;
			imaskblk[i] = -1;
		};
	};

// Compute the set of Schur submatrices

	CSMatrix *mtrarr;
	CSMatrix *mtrblkarr;

	mtrarr = new CSMatrix [_nblks];
	if(!mtrarr) MemoryFail(funcname);
	mtrblkarr = new CSMatrix [_nblks];
	if(!mtrblkarr) MemoryFail(funcname);

	ComputeSetOfSchurComplements (_comm, _ncycle, 
												_nblks, _blks, _blk2cpu, imaskblk,
												mtrarr, mtrblkarr);

// Compute new partitionings

	int **porderblk;
	int *nblksnewarr;
	int **pblksnew;
	int **pblk2colornew;

	porderblk = new int * [_nblks];
	if(!porderblk) MemoryFail(funcname);
	nblksnewarr = new int [_nblks];
	if(!nblksnewarr) MemoryFail(funcname);
	pblksnew = new int * [_nblks];
	if(!pblksnew) MemoryFail(funcname);
	pblk2colornew = new int * [_nblks];
	if(!pblk2colornew) MemoryFail(funcname);

	for (iblk=0;iblk<_nblks;iblk++) nblksnewarr[iblk] = 0;

	int decomptypeloc;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (imaskblk[iblk] == 0 && _blk2cpu[iblk] == myid) {

			icolor = _blk2color[iblk];

			int niloc = _blks[iblk+1]-_blks[iblk];

			porderblk[iblk] = new int [niloc];
			if(!porderblk[iblk]) MemoryFail(funcname);

			decomptypeloc = 1;

			mtrarr[iblk].PartSchurSubmatrix (decomptypeloc, icolor, 
														_blk2color, imaskblk, mtrblkarr[iblk],
														porderblk[iblk],
														nblksnewarr[iblk], pblksnew[iblk], pblk2colornew[iblk]);

		};
		if (imaskblk[iblk] == -1 && _blk2cpu[iblk] == myid) nblksnewarr[iblk] = 1;
	};

// Create global arrays

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblks, nblksnewarr, nblksnewarr);

	_blknew2ext = new int [_nblks+1];
	if(!_blknew2ext) MemoryFail(funcname);

	_blknew2ext[0] = 0;
	for (iblk=0;iblk<_nblks;iblk++) _blknew2ext[iblk+1] = nblksnewarr[iblk];

	for (iblk=0;iblk<_nblks;iblk++) _blknew2ext[iblk+1] = _blknew2ext[iblk]+_blknew2ext[iblk+1];

	_nblksnew = _blknew2ext[_nblks];

	_blksnew = new int [_nblksnew+1];
	if(!_blksnew) MemoryFail(funcname);
	_blk2cpunew = new int [_nblksnew];
	if(!_blk2cpunew) MemoryFail(funcname);
	_blk2colornew = new int [_nblksnew];
	if(!_blk2colornew) MemoryFail(funcname);

	for (iblk=0;iblk<=_nblksnew;iblk++) _blksnew[iblk] = 0;
	for (iblk=0;iblk<_nblksnew;iblk++) _blk2colornew[iblk] = 0;

	int j, icpu, ibeg;
	int ibs = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (imaskblk[iblk] == 0 && _blk2cpu[iblk] == myid) {
			for (j=_blknew2ext[iblk];j<_blknew2ext[iblk+1];j++) {
				_blksnew[j+1] = pblksnew[iblk][j-_blknew2ext[iblk]+1]-pblksnew[iblk][j-_blknew2ext[iblk]];
				_blk2colornew[j] = pblk2colornew[iblk][j-_blknew2ext[iblk]];
			};
			ibeg = _blks[iblk];
			for (j=_blks[iblk];j<_blks[iblk+1];j++) {
				_order[ibs] = porderblk[iblk][j-ibeg]+ibeg;
				ibs++;
			};
		};
		if (imaskblk[iblk] == -1 && _blk2cpu[iblk] == myid) {
			j = _blknew2ext[iblk];
			_blksnew[j+1] = _blks[iblk+1]-_blks[iblk];
			_blk2colornew[j] = _blk2color[iblk];
			for (j=_blks[iblk];j<_blks[iblk+1];j++) {
				_order[ibs] = j;
				ibs++;
			};
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblksnew+1, _blksnew, _blksnew);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblksnew, _blk2colornew, _blk2colornew);

	for (iblk=0;iblk<_nblksnew;iblk++) _blksnew[iblk+1] = _blksnew[iblk]+_blksnew[iblk+1];

	int ncolor2cpu = _ncolors / nproc;

	for (iblk=0;iblk<_nblksnew;iblk++) {
		icolor = _blk2colornew[iblk];
		icpu = icolor / ncolor2cpu;
		_blk2cpunew[iblk] = icpu;
	};

// Free work arrays for set of blocks

	for (iblk=0;iblk<_nblks;iblk++) {
		if (imaskblk[iblk] == 0 && _blk2cpu[iblk] == myid) {
			delete [] porderblk[iblk];
			delete [] pblksnew[iblk];
			delete [] pblk2colornew[iblk];
		};
	};

// Restore the number of blocks

	_nblksnew = _nblks;

// Free work arrays

	delete [] imaskblk;
	delete [] imaskcolor;
	delete [] mtrarr;
	delete [] mtrblkarr;
	delete [] porderblk;
	delete [] nblksnewarr;
	delete [] pblksnew;
	delete [] pblk2colornew;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute on one cpu partitioning of the matrix into many parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::DecompSurfaceRecursive1Cpu (CMPIComm &_comm1cpu, int _decomptype, // Compute on one cpu partitioning of the matrix into many parts taking into account blocks partitioning
														int _nblks, int *_blks, 
														int *&_order,
														int &_nblksnew, int *&_blksnew, int *&_blk2colornew) const {

	const char *funcname = "DecompSurfaceRecursive1Cpu";

// Check that the number of parts is a power of 2

	int nparts_pwr2 = 1;

	while (nparts_pwr2*2 <= _nblks) nparts_pwr2 *= 2;

	if (nparts_pwr2 != _nblks) throw " CSMatrix::DecompSchurBndRecursive1Cpu: the number of parts is not a power of 2 ";

// Allocate work memory

	int *listnd;
	int *listnd_new;
	int *iorderloc;
	int *imasknd;
	int *irow2ind;
	int *nsizes_tails;
	int *list_tails;
	int *blks_old;
	int *blks_new;
	int *blksloc;
	int *blk2colorloc;
	int *partblkloc;

	listnd = new int [n];
	if(!listnd) MemoryFail(funcname);
	listnd_new = new int [n];
	if(!listnd_new) MemoryFail(funcname);
	iorderloc = new int [n];
	if(!iorderloc) MemoryFail(funcname);
	imasknd = new int [n];
	if(!imasknd) MemoryFail(funcname);
	irow2ind = new int [n];
	if(!irow2ind) MemoryFail(funcname);
	nsizes_tails = new int [_nblks];
	if(!nsizes_tails) MemoryFail(funcname);
	list_tails = new int [n];
	if(!list_tails) MemoryFail(funcname);
	blks_old = new int [_nblks+1];
	if(!blks_old) MemoryFail(funcname);
	blks_new = new int [_nblks+1];
	if(!blks_new) MemoryFail(funcname);
	blksloc = new int [_nblks+1];
	if(!blksloc) MemoryFail(funcname);
	blk2colorloc = new int [_nblks+1];
	if(!blk2colorloc) MemoryFail(funcname);
	partblkloc = new int [_nblks+1];
	if(!partblkloc) MemoryFail(funcname);

	int i, j;

	int icycle = -1;

	for (i=0;i<n;i++) imasknd[i] = icycle;
	for (i=0;i<n;i++) listnd[i] = i;
	for (i=0;i<n;i++) listnd_new[i] = i;

// Implement main cycle over the levels

	nparts_pwr2 = 1;

	for (i=0;i<=_nblks;i++) blks_old[i] = _blks[i];

	CSMatrix mtrsub;
	CSMatrix mtrdummy;

	int ntails = 0;
	int ntot_tails = 0;

	int nparts_pwr2_new, ibs, iblk, nlistnd;
	int niblk, ihblk, ibegblk, iendblk, nblksloc, nblks1loc, nblksnew, iold, iblkloc, nz;
	int *orderloc, *blksnew, *blk2cpunew, *partblknew;

	while (nparts_pwr2 != _nblks) {
		nparts_pwr2_new = nparts_pwr2 * 2;

		niblk = _nblks / nparts_pwr2;

		for (ihblk=0;ihblk<nparts_pwr2;ihblk++) {

// Get the submatrix corresponding to the sublist

			ibegblk = ihblk*niblk;
			iendblk = (ihblk+1)*niblk-1;

			nlistnd = blks_old[iendblk+1]-blks_old[ibegblk];
			ibs = blks_old[ibegblk];

			SubmatrixList (nlistnd, listnd+ibs,
								icycle, imasknd, irow2ind,
								mtrsub);

// Create local partitionings

			nblksloc = iendblk-ibegblk+1;

			blksloc[0] = 0;
			for (iblk=ibegblk;iblk<=iendblk;iblk++) {
				blksloc[iblk-ibegblk+1] = blksloc[iblk-ibegblk]+(blks_old[iblk+1]-blks_old[iblk]);
			};
			nblks1loc = nblksloc / 2;
			if (nblks1loc == 0) nblks1loc = 1;
			for (iblk=0;iblk<nblksloc;iblk++) {
				blk2colorloc[iblk] = 0;
				partblkloc[iblk] = iblk / nblks1loc;
			};

// Compute new partitioning

			mtrsub.ComputeGroupsBoundary (_comm1cpu, 
													nblksloc, blksloc, blk2colorloc, partblkloc,
													orderloc, 
													nblksnew, blksnew, blk2cpunew, partblknew);

			mtrsub = mtrdummy;

// Inverse order

			for (i=0;i<nlistnd;i++) iorderloc[orderloc[i]] = i;;

// Store computed local partitioning and ordering

			if (nblksnew != nblksloc+1) throw " CSMatrix::DecompSchurBndRecursive1Cpu: the number of parts is not correct ";

			for (i=blksnew[nblksloc];i<blksnew[nblksloc+1];i++) {
				iold = iorderloc[i];
				list_tails[ntot_tails] = listnd[ibs+iold];
				ntot_tails++;
			};
			nsizes_tails[ntails] = blksnew[nblksloc+1]-blksnew[nblksloc];
			ntails++;

			blks_new[0] = 0;
			for (iblk=ibegblk;iblk<=iendblk;iblk++) {
				iblkloc = iblk-ibegblk;
				nz = blks_new[iblk];
				for (i=blksnew[iblkloc];i<blksnew[iblkloc+1];i++) {
					iold = iorderloc[i];
					listnd_new[nz] = listnd[ibs+iold];
					nz++;
				};
				blks_new[iblk+1] = nz;
			};

// Free computed arrays

			delete [] orderloc;
			delete [] blksnew;
			delete [] blk2cpunew;
			delete [] partblknew;

		};

// Prepare new partioning

		for (i=0;i<n;i++) listnd[i] = listnd_new[i];

		for (i=0;i<=_nblks;i++) blks_old[i] = blks_new[i];

		nparts_pwr2 = nparts_pwr2_new;

	};

// Perform final ND reordering

	for (iblk=0;iblk<_nblks;iblk++) {

		nlistnd = blks_old[iblk+1]-blks_old[iblk];
		ibs = blks_old[iblk];

		SubmatrixList (nlistnd, listnd+ibs,
							icycle, imasknd, irow2ind,
							mtrsub);

		if (nlistnd > 0) mtrsub.OrderPrfMtr (-1, listnd_new);

		for (i=0;i<nlistnd;i++) iorderloc[listnd_new[i]] = i;

		for (i=0;i<nlistnd;i++) {
			iold = iorderloc[i];
			listnd_new[i] = listnd[ibs+iold];
		};

		for (i=0;i<nlistnd;i++) listnd[ibs+i] = listnd_new[i];

	};

// Compute the ordering of blocks and tails

	int *blk2ind;
	int *tailblk2ind;
	int *imaskblk;

	blk2ind = new int [_nblks];
	if(!blk2ind) MemoryFail(funcname);
	tailblk2ind = new int [ntails];
	if(!tailblk2ind) MemoryFail(funcname);
	imaskblk = new int [2*_nblks];
	if(!imaskblk) MemoryFail(funcname);

	for (i=0;i<2*_nblks;i++) imaskblk[i] = -1;

	double *perf;
	double *memory;

	perf = new double [_nblks];
	if(!perf) MemoryFail(funcname);
	memory = new double [_nblks];
	if(!memory) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) perf[i] = 1.0;
	for (i=0;i<_nblks;i++) memory[i] = 1.0;

	CTree treebinary (_nblks, 2, perf, memory);

	treebinary.SortNodes ();

	delete [] perf;
	delete [] memory;

	int nlev = treebinary.GetNlev ();
	int nnodes = treebinary.GetNnodes ();
	CNode *pnodes = treebinary.GetNode ();

	nz = 0;

	int ilev;

	for (ilev=0;ilev<nlev;ilev++) {
		if (nz < ntails) {
			for (j=0;j<nnodes;j++) {
				if (pnodes[j].GetNodelv() == ilev) {
					tailblk2ind[nz] = j;
					nz++;
				};
			};
		};
	};

	int ind;

	for (i=0;i<ntails;i++) {
		ind = tailblk2ind[i];
		imaskblk[ind] = 1;
	};

	iblk = 0;

	for (i=0;i<2*_nblks;i++) {
		if (imaskblk[i] == -1) {
			if (iblk < _nblks) {
				blk2ind[iblk] = i;
				iblk++;
			};
		};
	};

// Prepare the final output

	_nblksnew = _nblks+ntails;

	_order = new int [n];
	if(!_order) MemoryFail(funcname);
	_blksnew = new int [_nblksnew+1];
	if(!_blksnew) MemoryFail(funcname);
	_blk2colornew = new int [_nblksnew];
	if(!_blk2colornew) MemoryFail(funcname);

	for (i=0;i<=_nblksnew;i++) _blksnew[i] = 0;

	for (i=0;i<_nblks;i++) {
		ind = blk2ind[i];
		_blksnew[ind+1] = blks_old[i+1]-blks_old[i];
	};
	for (i=0;i<ntails;i++) {
		ind = tailblk2ind[i];
		_blksnew[ind+1] = nsizes_tails[i];
	};

	for (i=0;i<_nblksnew;i++) _blksnew[i+1] = _blksnew[i]+_blksnew[i+1];

	int ibsnew;

	for (i=0;i<_nblks;i++) {
		ibs = blks_old[i];
		nlistnd = blks_old[i+1]-blks_old[i];
		ind = blk2ind[i];
		ibsnew = _blksnew[ind];
		for (j=0;j<nlistnd;j++) {
			listnd_new[ibsnew+j] = listnd[ibs+j];
		};
	};
	nz = 0;
	for (i=0;i<ntails;i++) {
		nlistnd = nsizes_tails[i];
		ind = tailblk2ind[i];
		ibsnew = _blksnew[ind];
		for (j=0;j<nlistnd;j++) {
			listnd_new[ibsnew+j] = list_tails[nz];
			nz++;
		};
	};

	for (i=0;i<n;i++) _order[listnd_new[i]] = i;

	for (i=0;i<_nblks;i++) {
		ind = blk2ind[i];
		_blk2colornew[ind] = i;
	};

	int ncolorsloc = 1;
	nz = 0;

	while (nz < ntails) {
		int ncolorsshft = _nblks / ncolorsloc;
		for (i=0;i<ncolorsloc;i++) {
			ind = tailblk2ind[nz];
			_blk2colornew[ind] = i*ncolorsshft;
			nz++;
		};
		ncolorsloc *= 2;
	};

// Free work arrays

	delete [] listnd;
	delete [] listnd_new;
	delete [] iorderloc;
	delete [] imasknd;
	delete [] irow2ind;
	delete [] nsizes_tails;
	delete [] list_tails;
	delete [] blks_old;
	delete [] blks_new;
	delete [] blksloc;
	delete [] blk2colorloc;
	delete [] partblkloc;
	delete [] blk2ind;
	delete [] tailblk2ind;
	delete [] imaskblk;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
//========================================================================================
void CSMatrix::ComputeGroupsBoundary (CMPIComm &_comm, // Compute in parallel partitioning of the matrix into two parts taking into account blocks partitioning
													int _nblks, int *_blks, int *_blk2cpu, int *_partitionblk,
													int *&_order, 
													int &_nblksnew, int *&_blksnew, int *&_blk2cpunew, int *&_partitionblknew) const {

	const char *funcname = "ComputeGroupsBoundary";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkBnd_",myid,".dat");
//	ofstream fout (strbuff);

// Compute statistics on entry

	int statpart[] = {0,0};

	int i;

	for (i=0;i<_nblks;i++) {
		if (_partitionblk[i] == 0) {
			statpart[0] += _blks[i+1]-_blks[i];
		} else {
			statpart[1] += _blks[i+1]-_blks[i];
		};
	};

// Compute local base

	int *ibsblk;

	ibsblk = new int [_nblks];
	if (!ibsblk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) ibsblk[i] = -1;

	int nz = 0;

	for (i=0;i<_nblks;i++) {
		if (_blk2cpu[i] == myid) {
			ibsblk[i] = nz;
			nz += (_blks[i+1]-_blks[i]);
		};
	};

// Check partitioning index for own data

	int *list2blk;
	int *imaskndini;
	int *imasknd;

	list2blk = new int [nlist];
	if (!list2blk) MemoryFail (funcname);
	imaskndini = new int [nlist];
	if (!imaskndini) MemoryFail (funcname);
	imasknd = new int [nlist];
	if (!imasknd) MemoryFail (funcname);

	int iblkprev = 0;

	int jj, iblk, ibegblk, iendblk, icpu, j;

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

// For each element find its block number

	int *ja2blk;

	ja2blk = new int [nzja];
	if (!ja2blk) MemoryFail (funcname);
/*
	iblkprev = 0;

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

// Create the list of local boundary elements

	int nlistbnd = 0;

	int *listbnd;
	int *listbndind;

	listbnd = new int [nlist];
	if (!listbnd) MemoryFail (funcname);
	listbndind = new int [nlist];
	if (!listbndind) MemoryFail (funcname);

	int ipart, icheck, jblk, jpart;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		ipart = _partitionblk[iblk];
		imaskndini[i] = ipart;
		imasknd[i] = ipart;
		icheck = 1;
		for (j=ia[i];j<ia[i+1];j++) {
			jblk = ja2blk[j];
			jpart = _partitionblk[jblk];
			if (jpart != ipart) icheck = -1;
		};
		if (icheck == -1) {
			listbnd[nlistbnd] = list[i];
			listbndind[nlistbnd] = i;
			nlistbnd++;
		};
	};

// Count the number of boundary nodes on each cpu

	int *blksbnd;

	blksbnd = new int [nproc+1];
	if (!blksbnd) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) blksbnd[i] = 0;

	blksbnd[myid+1] = nlistbnd;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksbnd, blksbnd);

	for (i=0;i<nproc;i++) blksbnd[i+1] = blksbnd[i]+blksbnd[i+1];

// Create the new ordering

	int *orderbnd;

	orderbnd = new int [nlist];
	if (!orderbnd) MemoryFail (funcname);

	nz = blksbnd[myid];

	for (i=0;i<nlist;i++) orderbnd[i] = -1;

	int jind;

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		orderbnd[jind] = nz;
		imasknd[jind] = -1;
		nz++;
	};

// Count statistics of the remaining volume nodes

	statpart[0] = 0;
	statpart[1] = 0;

	for (i=0;i<nlist;i++) {
		if (imasknd[i] == 0) {
			statpart[0]++;
		} else if (imasknd[i] == 1) {
			statpart[1]++;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 2, statpart, statpart);

// Create the local submatrix with the old rows indices

	int nzjaloc = 0;

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		nzjaloc += (ia[jind+1]-ia[jind]);
	};

	CSMatrix aloc (nlistbnd,nzjaloc);

	int *ja2blkloc;

	ja2blkloc = new int [nzjaloc];
	if (!ja2blkloc) MemoryFail (funcname);

	nzjaloc = 0;
	aloc.ia[0] = 0;

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		aloc.list[i] = orderbnd[jind];
		for (j=ia[jind];j<ia[jind+1];j++) {
			aloc.ja[nzjaloc] = ja[j];
			ja2blkloc[nzjaloc] = ja2blk[j];
			nzjaloc++;
		};
		aloc.ia[i+1] = nzjaloc;
	};

	aloc.m = blksbnd[nproc];
	aloc.n = blksbnd[nproc];
	aloc.nlist = nlistbnd;
	aloc.nzja = nzjaloc;

// For all indices in the submatrix find its position

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	for (i=0;i<nlistbnd;i++) {
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jblk = ja2blkloc[j];
			icpu = _blk2cpu[jblk];
			iacpu[icpu+1]++;
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	nz = iacpu[nproc];

	int *isendarr;

	isendarr = new int [nz*3];
	if (!isendarr) MemoryFail (funcname);

	int k;

	for (i=0;i<nlistbnd;i++) {
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jblk = ja2blkloc[j];
			icpu = _blk2cpu[jblk];
			k = iptr[icpu];
			isendarr[k*3] = j;
			isendarr[k*3+1] = aloc.ja[j];
			isendarr[k*3+2] = jblk;
			iptr[icpu]++;
		};
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
	if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int jjnew, ni, ibs, jjloc;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (3*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[3*j+1];
			jblk = piarr[3*j+2];
			ibs = ibsblk[jblk];
			jjloc = jj-_blks[jblk];
			jjnew = orderbnd[jjloc+ibs];
			piarr[3*j+1] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Renumber the column indices

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (3*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			jind = piarr[j*3];
			jjnew = piarr[j*3+1];
			aloc.ja[jind] = jjnew;
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

// Filter undefined column indices

	int *ialoc;

	ialoc = new int [nlistbnd+1];
	if(!ialoc) MemoryFail(funcname);

	ialoc[0] = 0;
	nz = 0;

	for (i=0;i<nlistbnd;i++) {
		for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
			jj = aloc.ja[j];
			if (jj >= 0 && jj != aloc.list[i]) {
				aloc.ja[nz] = jj;
				nz++;
			};
		};
		ialoc[i+1] = nz;
	};

	for (i=0;i<=nlistbnd;i++) aloc.ia[i] = ialoc[i];

// Partition in parallel the boundary into two parts

	int *partitionbnd;

	partitionbnd = new int [nlistbnd];
	if(!partitionbnd) MemoryFail(funcname);

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

// Compute weights

//	for (i=0;i<nparts;i++) tpwgts[i] = (float)1.0 / (float (nparts));

	int nlistbnd_2 = nlistbnd / 2;

	int deltav = statpart[1]-statpart[0];

	if (deltav >= nlistbnd_2) {
		tpwgts[0] = 1.0;
		tpwgts[1] = 0.0;
	} else if (deltav <= -nlistbnd_2) {
		tpwgts[0] = 0.0;
		tpwgts[1] = 1.0;
	} else {
		double aux;
		float faux;
		aux = (double) (deltav+nlistbnd_2) / (double) nlistbnd;
		faux = (float) aux;
		if (faux > 1.0) faux = 1.0;
		if (faux < 0.0) faux = 0.0;
		float faux1 = (float) (1.0)-faux;
		if (faux1 > 1.0) faux1 = 1.0;
		if (faux1 < 0.0) faux1 = 0.0;
		tpwgts[0] = faux;
		tpwgts[1] = faux1;
	};

	double dmin = 0.005;
	double dmax = 0.995;

	if (tpwgts[0] < dmin) tpwgts[0] = (float)dmin;
	if (tpwgts[0] > dmax) tpwgts[0] = (float)dmax;

	tpwgts[1] = (float) (1.0) - tpwgts[0];

	ubvec[0] = (float) 1.05;

	int edgecut;
	int iopt = 0;

	int nlistcpu = 0;
	int *listcpu;
	int *blksbndsm;

	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	blksbndsm = new int [nproc+1];
	if (!blksbndsm) MemoryFail (funcname);

	blksbndsm[0] = 0;

	int icheckcpu = -1;

	for (i=0;i<nproc;i++) {
		if (blksbnd[i+1] > blksbnd[i]) {
			listcpu[nlistcpu] = i;
			blksbndsm[nlistcpu+1] = blksbnd[i+1];
			nlistcpu++;
			if (i == myid) icheckcpu = 1;
		};
	};

	CMPIComm commloc = _comm.CreateComm (nlistcpu, listcpu);

	if (icheckcpu == 1) {
/*
		if (true) {
			int myidloc = commloc.GetMyid ();
			CSMatrix atot;
			int nploc = 0;
			int ibegloc = 0;
			int iendloc = -1;
			if (myidloc == 0) {
				nploc = 1;
				ibegloc = 0;
				iendloc = aloc.GetN ();
			};
			atot = aloc.GetSubmatrixViaIntervals (commloc, nploc, &ibegloc, &iendloc);
			int nzloc = atot.GetNzja ();
			if (myidloc == 0) {
				char strbuff[256];
				sprintf (strbuff,"%s%i%s","Matr_",nzloc,".dat");
				ofstream fout (strbuff);
				fout << " Matr = " << atot << endl;
				sprintf (strbuff,"%s%i%s","Matr_",nzloc,".ps");
				atot.A2Ps (strbuff,0,&i);
			};
		};
*/
		if (nlistcpu > 1) {
/*
			if (false) {

				int *blk2cpuloc;

				blk2cpuloc = new int [nproc];
				if (!blk2cpuloc) MemoryFail (funcname);

				for (i=0;i<nproc;i++) blk2cpuloc[i] = i;

				CSMatrix::CheckSymmMtr (commloc,
												nlistcpu, blksbndsm, blk2cpuloc,
												aloc.ia, aloc.ja,
												false, NULL);

				delete [] blk2cpuloc;

			};
*/
//			if (myid == 0) {
//				cout << " CSMatrix::ComputeGroupsBoundary: before PartKway " << endl;
//			};

			if (false) {

				CParMETIS::PartKway (commloc, blksbndsm, aloc.ia, aloc.ja, 
											NULL, NULL, &iopt,
											&numflag, &ncon, &nparts,
											tpwgts, ubvec,
											options, &edgecut, partitionbnd);

			} else {
				int nploc = 0;
				int *ibegarrloc;
				int *iendarrloc;

				ibegarrloc = new int [nlistcpu];
				if (!ibegarrloc) MemoryFail (funcname);
				iendarrloc = new int [nlistcpu];
				if (!iendarrloc) MemoryFail (funcname);

				int myidloc = commloc.GetMyid ();

				if (myidloc == 0) {
					for (i=0;i<nlistcpu;i++) {
						if (blksbndsm[i] < blksbndsm[i+1]) {
							ibegarrloc[nploc] = blksbndsm[i];
							iendarrloc[nploc] = blksbndsm[i+1]-1;
							nploc++;
						};
					};
				};

				CSMatrix a1cpu;
				a1cpu = aloc.GetSubmatrixViaIntervals (commloc,
																		nploc, ibegarrloc, iendarrloc);

				delete [] ibegarrloc;
				delete [] iendarrloc;

				int nnntot = blksbndsm[nlistcpu];

				int *partitionbnd_1cpu;

				partitionbnd_1cpu = new int [nnntot];
				if (!partitionbnd_1cpu) MemoryFail (funcname);

				for (i=0;i<nnntot;i++) partitionbnd_1cpu[i] = 0;

				if (myidloc == 0) {

					METIS_WPartGraphRecursive (&nnntot, a1cpu.ia, a1cpu.ja, NULL, NULL, &iopt,
														&numflag, &nparts, tpwgts,
														options, &edgecut, partitionbnd_1cpu);

				};

				CMPIExchange::ExchangeArrayMPI (commloc, INTEGERVALUE, ADD, nnntot, partitionbnd_1cpu, partitionbnd_1cpu);

				for (i=blksbndsm[myidloc];i<blksbndsm[myidloc+1];i++) {
					partitionbnd[i-blksbndsm[myidloc]] = partitionbnd_1cpu[i];
				};

				delete [] partitionbnd_1cpu;

			};

//			if (myid == 0) {
//				cout << " CSMatrix::ComputeGroupsBoundary: after PartKway " << endl;
//			};

		} else {
			int nloc = blksbndsm[nlistcpu];
//			METIS_PartGraphKway (&nloc, aloc.ia, aloc.ja, NULL, NULL, &iopt,
//										&numflag, &nparts, 
//										options, &edgecut, partitionbnd);
			METIS_WPartGraphRecursive (&nloc, aloc.ia, aloc.ja, NULL, NULL, &iopt,
										&numflag, &nparts, tpwgts,
										options, &edgecut, partitionbnd);
		};
	};

	delete [] listcpu;
	delete [] blksbndsm;

	CSMatrix mtrdummy;

	aloc = mtrdummy;

// Mark nodes

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		ipart = partitionbnd[i];
		ipart += 2;
		imasknd[jind] = ipart;
	};

// Exchange boundary information once again

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		for (j=ia[jind];j<ia[jind+1];j++) {
			jblk = ja2blk[j];
			icpu = _blk2cpu[jblk];
			iacpu[icpu+1]++;
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	nz = iacpu[nproc];

	delete [] isendarr;

	isendarr = new int [nz*3];
	if (!isendarr) MemoryFail (funcname);

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		for (j=ia[jind];j<ia[jind+1];j++) {
			jblk = ja2blk[j];
			icpu = _blk2cpu[jblk];
			k = iptr[icpu];
			isendarr[k*3] = j;
			isendarr[k*3+1] = ja[j];
			isendarr[k*3+2] = jblk;
			iptr[icpu]++;
		};
	};

// Exchange the data

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

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (3*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[3*j+1];
			jblk = piarr[3*j+2];
			ibs = ibsblk[jblk];
			jjloc = jj-_blks[jblk];
			jjnew = imasknd[jjloc+ibs];
			piarr[3*j+1] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store the result

	int *ja2ipart;

	ja2ipart = new int [nzja];
	if (!ja2ipart) MemoryFail (funcname);

	for (i=0;i<nzja;i++) ja2ipart[i] = -1;

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / (3*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			jind = piarr[j*3];
			jjnew = piarr[j*3+1];
			ja2ipart[jind] = jjnew;
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

// Modify the partitioning and create the final boundary

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		ipart = imaskndini[jind];
		for (j=ia[jind];j<ia[jind+1];j++) {
			jpart = ja2ipart[j];
			if (jpart != ipart && jpart != ipart+2) {
				imasknd[jind] = -1;
			};
		};
	};

	for (i=0;i<nlistbnd;i++) {
		jind = listbndind[i];
		ipart = imaskndini[jind];
		if (imasknd[jind] != -1) imasknd[jind] = ipart;
	};

// Make the boundary smaller if possible

	if (true) {

		int ncycles = 1;

		int ic, iipart, ipart1, nmodmax;

		for (ic=0;ic<ncycles;ic++) {

			for (ipart=0;ipart<2;ipart++) {

				ipart1 = (ipart+1) % 2;

// Exchange the data

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

				info = CMPIExchange::DataExchangeMPI (_comm,
																	NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
																	ObjSizeSend, ObjSend,
																	NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
																	ObjSizeRecv, ObjRecv);
				if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

				delete [] ObjTypeSend;
				delete [] ObjIDSend;
				delete [] CpuIDSend;
				delete [] ObjSizeSend;
				delete [] ObjSend;

// Prepare reply

				for (i=0;i<NObjRecv;i++) {
					ni = ObjSizeRecv[i] / (3*sizeof(int));
					piarr = (int *) ObjRecv[i];
					for (j=0;j<ni;j++) {
						jj = piarr[3*j+1];
						jblk = piarr[3*j+2];
						ibs = ibsblk[jblk];
						jjloc = jj-_blks[jblk];
						jjnew = imasknd[jjloc+ibs];
						piarr[3*j+1] = jjnew;
					};
				};

// Send indices back

				info = CMPIExchange::DataExchangeMPI (_comm,
																	NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
																	ObjSizeRecv, ObjRecv,
																	NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
																	ObjSizeSend, ObjSend);
				if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

				delete [] ObjTypeRecv;
				delete [] ObjIDRecv;
				delete [] CpuIDRecv;
				delete [] ObjSizeRecv;
				for(i = 0; i < NObjRecv; i++) {
					delete [] ObjRecv[i];
				};
				delete [] ObjRecv;

// Store the result

				for (i=0;i<NObjSend;i++) {
					icpu = CpuIDSend[i];
					ni = ObjSizeSend[i] / (3*sizeof(int));
					piarr = (int *) ObjSend[i];
					for (j=0;j<ni;j++) {
						jind = piarr[j*3];
						jjnew = piarr[j*3+1];
						ja2ipart[jind] = jjnew;
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

// Compute the maximal allowable data to be modified

				nmodmax = nlistbnd;

				if (ipart == 0) {
					nmodmax = 0;

					for (i=0;i<nlistbnd;i++) {
						jind = listbndind[i];
						iipart = imaskndini[jind];
//						if (iipart == ipart && imasknd[jind] == -1) {
						if (imasknd[jind] == -1) {
							icheck = 1;
							for (j=ia[jind];j<ia[jind+1];j++) {
								jpart = ja2ipart[j];
								if (jpart == ipart1) {
									icheck = -1;
								};
							};
							if (icheck == 1) {
								nmodmax++;
//								fout << " Modification: jind = " << jind << " irow = " << list[jind] << endl;
//								OutArr (fout," Indices ",ia[jind+1]-ia[jind],ja+ia[jind]);
//								OutArr (fout," Parts ",ia[jind+1]-ia[jind],ja2ipart+ia[jind]);
							};
						};
					};

					nmodmax = nmodmax / 2;

				};

// Modify the partitioning and create the final boundary

				int nmod = 0;

				for (i=0;i<nlistbnd;i++) {
					jind = listbndind[i];
					iipart = imaskndini[jind];
					if (iipart == ipart && imasknd[jind] == -1) {
//					if (imasknd[jind] == -1) {
						icheck = 1;
						for (j=ia[jind];j<ia[jind+1];j++) {
							jpart = ja2ipart[j];
							if (jpart == ipart1) {
								icheck = -1;
							};
						};
						if (icheck == 1 && nmod < nmodmax) {
							imasknd[jind] = ipart;
							nmod++;
//							fout << " Modification: jind = " << jind << " irow = " << list[jind] << endl;
//							OutArr (fout," Indices ",ia[jind+1]-ia[jind],ja+ia[jind]);
//							OutArr (fout," Parts ",ia[jind+1]-ia[jind],ja2ipart+ia[jind]);
						};
					};
				};

				int icheck1;

				for (i=0;i<nlistbnd;i++) {
					jind = listbndind[i];
					iipart = imaskndini[jind];
//					if (iipart == ipart && imasknd[jind] == -1) {
					if (imasknd[jind] == -1) {
						icheck = 1;
						icheck1 = -1;
						for (j=ia[jind];j<ia[jind+1];j++) {
							jpart = ja2ipart[j];
							if (jpart == ipart1) {
								icheck = -1;
							};
							if (jpart >= 0) icheck1 = 1;
						};
						if (icheck == 1 && icheck1 == 1 && nmod < nmodmax) {
							imasknd[jind] = ipart;
							nmod++;
//							fout << " Modification: jind = " << jind << " irow = " << list[jind] << endl;
//							OutArr (fout," Indices ",ia[jind+1]-ia[jind],ja+ia[jind]);
//							OutArr (fout," Parts ",ia[jind+1]-ia[jind],ja2ipart+ia[jind]);
						};
					};
				};

			};
		};
	};

// Final check of the partitioning
/*
	if (true) {

		for (i=0;i<nproc+1;i++) iacpu[i] = 0;

		for (i=0;i<nlist;i++) {
			jind = i;
			for (j=ia[jind];j<ia[jind+1];j++) {
				jblk = ja2blk[j];
				icpu = _blk2cpu[jblk];
				iacpu[icpu+1]++;
			};
		};

		for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
		for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

		nz = iacpu[nproc];

		delete [] isendarr;

		isendarr = new int [nz*3];
		if (!isendarr) MemoryFail (funcname);

		for (i=0;i<nlist;i++) {
			jind = i;
			for (j=ia[jind];j<ia[jind+1];j++) {
				jblk = ja2blk[j];
				icpu = _blk2cpu[jblk];
				k = iptr[icpu];
				isendarr[k*3] = j;
				isendarr[k*3+1] = ja[j];
				isendarr[k*3+2] = jblk;
				iptr[icpu]++;
			};
		};

// Exchange the data

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

		info = CMPIExchange::DataExchangeMPI (_comm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
		if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

		delete [] ObjTypeSend;
		delete [] ObjIDSend;
		delete [] CpuIDSend;
		delete [] ObjSizeSend;
		delete [] ObjSend;

// Prepare reply

		for (i=0;i<NObjRecv;i++) {
			ni = ObjSizeRecv[i] / (3*sizeof(int));
			piarr = (int *) ObjRecv[i];
			for (j=0;j<ni;j++) {
				jj = piarr[3*j+1];
				jblk = piarr[3*j+2];
				ibs = ibsblk[jblk];
				jjloc = jj-_blks[jblk];
				jjnew = imasknd[jjloc+ibs];
				piarr[3*j+1] = jjnew;
			};
		};

// Send indices back

		info = CMPIExchange::DataExchangeMPI (_comm,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend);
		if(info) throw " CSMatrix::ComputeGroupsBoundary: Error in DataExchangeMPI";

		delete [] ObjTypeRecv;
		delete [] ObjIDRecv;
		delete [] CpuIDRecv;
		delete [] ObjSizeRecv;
		for(i = 0; i < NObjRecv; i++) {
			delete [] ObjRecv[i];
		};
		delete [] ObjRecv;

// Store the result

		for (i=0;i<NObjSend;i++) {
			icpu = CpuIDSend[i];
			ni = ObjSizeSend[i] / (3*sizeof(int));
			piarr = (int *) ObjSend[i];
			for (j=0;j<ni;j++) {
				jind = piarr[j*3];
				jjnew = piarr[j*3+1];
				ja2ipart[jind] = jjnew;
			};
		};

		OutArr(fout," Listbnd ",nlistbnd,listbnd);

		for (i=0;i<nlist;i++) {
			ipart = imasknd[i];
			jind = i;
			if (ipart >= 0) {
				int icheck = 1;
				for (j=ia[jind];j<ia[jind+1];j++) {
					if (ja2ipart[j] >= 0 && ja2ipart[j] != ipart) {
						icheck = -1;
					};
				};
				if (icheck < 0) {
					cout << " Error in partitioning: i = " << i << " irow = " << list[i] << endl;
					fout << " Error in partitioning: i = " << i << " irow = " << list[i] << endl;
					fout << "     ini = " << imaskndini[i] << " fin = " << imasknd[i] << " Bndtype = " << orderbnd[i] << endl;
					OutArr (fout," Indices ",ia[i+1]-ia[i],ja+ia[i]);
					OutArr (fout," Parts ",ia[i+1]-ia[i],ja2ipart+ia[i]);
					for (j=ia[jind];j<ia[jind+1];j++) {
						jj = ja[j];
						jblk = ja2blk[j];
						jjloc = jj-_blks[jblk]+ibsblk[jblk];
						fout << "       jj = " << jj << " imaskini " << imaskndini[jjloc] << " imaskfin " << imasknd[jjloc] << endl;
					};
				};
			};
		};

	};
*/
// Determine the modified block number for each index

	int *list2blknew;

	list2blknew = new int [nlist];
	if (!list2blknew) MemoryFail (funcname);

	int iipart, iblknew;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		ipart = _partitionblk[iblk];
		iipart = imasknd[i];
		if (iipart > 1) iipart -= 2;
		if (iipart == -1) {
			iblknew = _nblks;
		} else {
			if (iipart == ipart) {
				iblknew = iblk;
			} else {
				icheck = -1;
				for (j=ia[i];j<ia[i+1];j++) {
					jblk = ja2blk[j];
					jpart = ja2ipart[j];
					if (jpart == iipart && _partitionblk[jblk] == iipart) {
						iblknew = jblk;
						icheck = 1;
						break;
					};
				};
				if (icheck == -1) {
					throw " CSMatrix::ComputeGroupsBoundary: Error when searching the block number";
				};
			};
		};
		list2blknew[i] = iblknew;
	};

// Create the final ordering of the matrix

	_nblksnew = _nblks + nproc;

	_blksnew = new int [_nblksnew+1];
	if (!_blksnew) MemoryFail (funcname);
	_blk2cpunew = new int [_nblksnew];
	if (!_blk2cpunew) MemoryFail (funcname);
	_partitionblknew = new int [_nblksnew];
	if (!_partitionblknew) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) _blk2cpunew[i] = _blk2cpu[i];
	for (i=0;i<nproc;i++) _blk2cpunew[_nblks+i] = i;
	for (i=0;i<_nblks;i++) _partitionblknew[i] = _partitionblk[i];
	for (i=0;i<nproc;i++) _partitionblknew[_nblks+i] = -1;

	for (i=0;i<=_nblksnew;i++) _blksnew[i] = 0;

	int nblksnewloc = (_nblks+1)*nproc;

	int *blksnewloc;

	blksnewloc = new int [nblksnewloc+1];
	if (!blksnewloc) MemoryFail (funcname);

	for (i=0;i<=nblksnewloc;i++) blksnewloc[i] = 0;

	for (i=0;i<nlist;i++) {
		iblk = list2blknew[i];
		blksnewloc[iblk*nproc+myid+1]++;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblksnewloc+1, blksnewloc, blksnewloc);

	for (i=0;i<nblksnewloc;i++) blksnewloc[i+1] = blksnewloc[i]+blksnewloc[i+1];

	for (i=0;i<_nblks;i++) {
		_blksnew[i+1] = blksnewloc[(i+1)*nproc];
	};
	for (i=0;i<nproc;i++) {
		_blksnew[_nblks+i+1] = blksnewloc[_nblks*nproc+i+1];
	};

	_order = new int [nlist];
	if (!_order) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		iblk = list2blknew[i];
		k = blksnewloc[iblk*nproc+myid];
		_order[i] = k;
		blksnewloc[iblk*nproc+myid]++;
	};

// Free work arrays

	delete [] ibsblk;
	delete [] list2blk;
	delete [] imaskndini;
	delete [] imasknd;
	delete [] ja2blk;
	delete [] listbnd;
	delete [] listbndind;
	delete [] blksbnd;
	delete [] orderbnd;
	delete [] ja2blkloc;
	delete [] iacpu;
	delete [] iptr;
	delete [] isendarr;
	delete [] ialoc;
	delete [] partitionbnd;
	delete [] tpwgts;
	delete [] ubvec;
	delete [] ja2ipart;
	delete [] list2blknew;
	delete [] blksnewloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel the set of Schur complements
//========================================================================================
void CSMatrix::ComputeSetOfSchurComplements (CMPIComm &_comm, int _ncycle, // Compute in parallel the sparsity of the extended boundary
															int _nblks, int *_blks, int *_blk2cpu, int *_imaskblk,
															CSMatrix *_mtrarr, CSMatrix *_mtrblkarr) const {

	const char *funcname = "ComputeSetOfSchurComplements";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Find partitioning for own data

	int *list2blk;

	list2blk = new int [nlist];
	if (!list2blk) MemoryFail (funcname);

	int jj, iblk, ibegblk, iendblk, icpu, j, i;

	int iblkprev = 0;

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

// For each element find its block number

	int *ja2blk;

	ja2blk = new int [nzja];
	if (!ja2blk) MemoryFail (funcname);

//	iblkprev = 0;
/*
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

// Determine the local block sparsity

	int *imaskblk;
	int *listblk;

	imaskblk = new int [_nblks];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [_nblks];
	if (!listblk) MemoryFail (funcname);

	int icycleblk = -1;

	for (i=0;i<_nblks;i++) imaskblk[i] = icycleblk;

	int nlistblk = 0;

	icycleblk++;
/*
	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		if (imaskblk[iblk] != icycleblk) {
			listblk[nlistblk] = iblk;
			nlistblk++;
			imaskblk[iblk] = icycleblk;
		};
	};
*/
	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			listblk[nlistblk] = iblk;
			nlistblk++;
			imaskblk[iblk] = icycleblk;
		};
	};

	int *ialist2blk;
	int *iablk;
	int *iptrblk;

	ialist2blk = new int [nlistblk+1];
	if (!ialist2blk) MemoryFail (funcname);
	iablk = new int [nlistblk+1];
	if (!iablk) MemoryFail (funcname);
	iptrblk = new int [_nblks];
	if (!iptrblk) MemoryFail (funcname);

	for (i=0;i<=nlistblk;i++) ialist2blk[i] = 0;
	for (i=0;i<=nlistblk;i++) iablk[i] = 0;

	for (i=0;i<_nblks;i++) iptrblk[i] = -1;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		iptrblk[iblk] = i;
	};

	int ilist, ind;

	for (ilist=0;ilist<nlist;ilist++) {
		iblk = list2blk[ilist];
		ind = iptrblk[iblk];
		if (ind >= 0) {
			ialist2blk[ind+1]++;
		} else {
			throw " CSMatrix::ComputeSetOfSchurComplements: the block element not found ";
		};
	};

	for (i=0;i<nlistblk;i++) ialist2blk[i+1] = ialist2blk[i]+ialist2blk[i+1];

	int ilistblk, nlistloc, jblk;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		icycleblk++;
		nlistloc = 0;
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jblk = ja2blk[j];
				if (imaskblk[jblk] != icycleblk) {
					nlistloc++;
					imaskblk[jblk] = icycleblk;
				};
			};
		};
		iablk[ilistblk+1] = nlistloc;
	};

	for (i=0;i<nlistblk;i++) iablk[i+1] = iablk[i]+iablk[i+1];

	for (i=0;i<nlistblk;i++) iptrblk[i] = iablk[i];

	int nzjablk = iablk[nlistblk];

	int *jablk;

	jablk = new int [nzjablk];
	if (!jablk) MemoryFail (funcname);

	int k, ibeg, ni, nj;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		icycleblk++;
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jblk = ja2blk[j];
				if (imaskblk[jblk] != icycleblk) {
					k = iptrblk[ilistblk];
					jablk[k] = jblk;
					imaskblk[jblk] = icycleblk;
					iptrblk[ilistblk]++;
				};
			};
		};
		ibeg = iablk[ilistblk];
		ni = iablk[ilistblk+1]-iablk[ilistblk];
		sort (jablk+ibeg, jablk+ibeg + ni);
	};

// Compute global block sparsity

	int *iablkglb;

	iablkglb = new int [_nblks+1];
	if (!iablkglb) MemoryFail (funcname);

	for (i=0;i<=_nblks;i++) iablkglb[i] = 0;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		iablkglb[iblk+1] = iablk[ilistblk+1]-iablk[ilistblk];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblks+1, iablkglb, iablkglb);

	for (i=0;i<_nblks;i++) iablkglb[i+1] = iablkglb[i]+iablkglb[i+1];

	int nzjablkglb = iablkglb[_nblks];

	int *jablkglb;

	jablkglb = new int [nzjablkglb];
	if (!jablkglb) MemoryFail (funcname);

	for (i=0;i<nzjablkglb;i++) jablkglb[i] = 0;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		for (j=iablkglb[iblk];j<iablkglb[iblk+1];j++) {
			jablkglb[j] = jablk[j-iablkglb[iblk]+iablk[ilistblk]];
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nzjablkglb, jablkglb, jablkglb);

// Compute the block sparsity of the A^TA

	CSMatrix ablkmtr (_nblks, nzjablkglb);

	int *plistblk = ablkmtr.GetList ();
	for (i=0;i<_nblks;i++) plistblk[i] = i;

	int *piablk = ablkmtr.GetIa ();
	delete [] piablk;
	ablkmtr.SetIa (iablkglb);

	int *pjablk = ablkmtr.GetJa ();
	delete [] pjablk;
	ablkmtr.SetJa (jablkglb);

	ablkmtr.m = _nblks;
	ablkmtr.n = _nblks;
	ablkmtr.nsupc = _nblks;
	ablkmtr.nsupr = _nblks;
	ablkmtr.nlist = _nblks;
	ablkmtr.nzja = nzjablkglb;

	CSMatrix atablkmtr;

	atablkmtr = ablkmtr.AtAMatrSpars ();

// Allocate the mask

	int nimax = _nblks;
	int niblkmax = _nblks;

	for (i=0;i<_nblks;i++) {
		ni = _blks[i+1]-_blks[i];
		if (niblkmax < ni) niblkmax = ni;
	};

	if (nimax < niblkmax) nimax = niblkmax;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ni = _blks[iblk+1]-_blks[iblk];
		if (nimax < ni) nimax = ni;
		ni = 0;
		for (j=iablk[i];j<iablk[i+1];j++) {
			jblk = jablk[j];
			nj = _blks[jblk+1]-_blks[jblk];
			ni += nj;
		};
		if (nimax < ni) nimax = ni;
	};

	int *ibsmask;
	int *imask;
	int *imask1;
	int *listarr;
	int *iptrarr;

	ibsmask = new int [_nblks];
	if (!ibsmask) MemoryFail (funcname);
	imask = new int [nimax];
	if (!imask) MemoryFail (funcname);
	imask1 = new int [niblkmax];
	if (!imask1) MemoryFail (funcname);
	listarr = new int [nimax];
	if (!listarr) MemoryFail (funcname);
	iptrarr = new int [nimax];
	if (!iptrarr) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) ibsmask[i] = -1;

	int icycle = -1;

	for (i=0;i<nimax;i++) imask[i] = icycle;

	int icycle1 = -1;

	for (i=0;i<niblkmax;i++) imask1[i] = icycle;

// Store the subblocks in the supersparse format

	int *niarr;
	int *nzjaarr;

	niarr = new int [_nblks];
	if (!niarr) MemoryFail (funcname);
	nzjaarr = new int [_nblks];
	if (!nzjaarr) MemoryFail (funcname);

	CSMatrix *amatrarr;
	CSMatrix *amatrrowarr;
	CSMatrix *amatrblkrowarr;

	amatrarr = new CSMatrix [nzjablk];
	if (!amatrarr) MemoryFail (funcname);
	amatrrowarr = new CSMatrix [_nblks];
	if (!amatrrowarr) MemoryFail (funcname);
	amatrblkrowarr = new CSMatrix [_nblks];
	if (!amatrblkrowarr) MemoryFail (funcname);

	CSMatrix mtrdummy;

	int ibs, jjloc, niloc, nzjloc, kkloc, irow, irowloc, kind, jbeg, njloc;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		ni = _blks[iblk+1]-_blks[iblk];
		ibs = 0;
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			ibsmask[jblk] = ibs;
			ibs += (_blks[jblk+1]-_blks[jblk]);
			niarr[jblk] = 0;
			nzjaarr[jblk] = 0;
		};
		icycle++;
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jblk = ja2blk[j];
				jjloc = jj-_blks[jblk];
				ind = jjloc + ibsmask[jblk];
				if (imask[ind] != icycle) {
					niarr[jblk]++;
					imask[ind] = icycle;
				};
				nzjaarr[jblk]++;
			};
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			niloc = niarr[jblk];
			nzjloc = nzjaarr[jblk];
			CSMatrix temp (niloc,nzjloc);
			amatrrowarr[jblk] = temp;
			niarr[jblk] = 0;
		};
		icycle++;
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jblk = ja2blk[j];
				jjloc = jj-_blks[jblk];
				ind = jjloc + ibsmask[jblk];
				if (imask[ind] != icycle) {
					amatrrowarr[jblk].list[niarr[jblk]] = jjloc;
					niarr[jblk]++;
					imask[ind] = icycle;
				};
			};
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			niloc = niarr[jblk];
			sort (amatrrowarr[jblk].list, amatrrowarr[jblk].list + niloc);
			for (k=0;k<=niloc;k++) amatrrowarr[jblk].ia[k] = 0;
			ibs = ibsmask[jblk];
			for (k=0;k<niloc;k++) {
				kkloc = amatrrowarr[jblk].list[k];
				listarr[ibs+kkloc] = k;
			};
		};
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jblk = ja2blk[j];
				jjloc = jj-_blks[jblk];
				ind = jjloc + ibsmask[jblk];
				k = listarr[ind];
				amatrrowarr[jblk].ia[k+1]++;
			};
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			niloc = niarr[jblk];
			for (k=0;k<niloc;k++) amatrrowarr[jblk].ia[k+1] = amatrrowarr[jblk].ia[k]+amatrrowarr[jblk].ia[k+1];
			ibs = ibsmask[jblk];
			for (k=0;k<niloc;k++) iptrarr[ibs+k] = amatrrowarr[jblk].ia[k];
		};
		for (i=ialist2blk[ilistblk];i<ialist2blk[ilistblk+1];i++) {
			irow = list[i];
			irowloc = irow-_blks[iblk];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jblk = ja2blk[j];
				jjloc = jj-_blks[jblk];
				ind = jjloc + ibsmask[jblk];
				ibs = ibsmask[jblk];
				k = listarr[ind];
				kind = iptrarr[ibs+k];
				amatrrowarr[jblk].ja[kind] = irowloc;
				iptrarr[ibs+k]++;
			};
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			niloc = niarr[jblk];
			for (k=0;k<niloc;k++) {
				jbeg = amatrrowarr[jblk].ia[k];
				njloc = amatrrowarr[jblk].ia[k+1]-amatrrowarr[jblk].ia[k];
				sort (amatrrowarr[jblk].ja+jbeg, amatrrowarr[jblk].ja+jbeg + njloc);
			};
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			nj = _blks[jblk+1]-_blks[jblk];
//			amatrrowarr[jblk].m = ni;
//			amatrrowarr[jblk].n = nj;
//			amatrrowarr[jblk].nsupr = ni;
//			amatrrowarr[jblk].nsupc = nj;
			amatrrowarr[jblk].m = niblkmax;
			amatrrowarr[jblk].n = niblkmax;
			amatrrowarr[jblk].nsupr = niblkmax;
			amatrrowarr[jblk].nsupc = niblkmax;
			amatrrowarr[jblk].nlist = niarr[jblk];
			amatrrowarr[jblk].nzja = nzjaarr[jblk];
			niarr[jblk] = 0;
			nzjaarr[jblk] = 0;
		};
		for (j=iablk[ilistblk];j<iablk[ilistblk+1];j++) {
			jblk = jablk[j];
			amatrarr[j] = amatrrowarr[jblk].TranspMtrMask (icycle, imask);
			amatrrowarr[jblk] = mtrdummy;
		};
	};

// Create the sends list

	int *iacpu;
	int *iptrcpu;
	int *imaskcpu;

	iacpu = new int [nproc+1];
	if(!iacpu) MemoryFail(funcname);
	iptrcpu = new int [nproc];
	if(!iptrcpu) MemoryFail(funcname);
	imaskcpu = new int [nproc];
	if(!imaskcpu) MemoryFail(funcname);

	int icyclecpu = -1;
	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int kblk;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=iablk[i];j<iablk[i+1];j++) {
			jblk = jablk[j];
			if (jblk > iblk) {
				icyclecpu++;
				for (k=iablk[i];k<j;k++) {
					kblk = jablk[k];
					icpu = _blk2cpu[kblk];
					if (imaskcpu[icpu] != icyclecpu) {
						iacpu[icpu+1]++;
						imaskcpu[icpu] = icyclecpu;
					};
				};
			};
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	int nzjablksend = iacpu[nproc];

	int *jindcpu;

	jindcpu = new int [3*nzjablksend];
	if(!jindcpu) MemoryFail(funcname);

	int kkk;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=iablk[i];j<iablk[i+1];j++) {
			jblk = jablk[j];
			if (jblk > iblk) {
				icyclecpu++;
				for (k=iablk[i];k<j;k++) {
					kblk = jablk[k];
					icpu = _blk2cpu[kblk];
					if (imaskcpu[icpu] != icyclecpu) {
						imaskcpu[icpu] = icyclecpu;
						kkk = iptrcpu[icpu];
						jindcpu[kkk*3] = iblk;
						jindcpu[kkk*3+1] = jblk;
						jindcpu[kkk*3+2] = j;
						iptrcpu[icpu]++;
					};
				};
			};
		};
	};

// Exchange block data for future multiply

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nzjablksend;

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

	for (i=0;i<nproc;i++) {
		for (j=iacpu[i];j<iacpu[i+1];j++) {
			iblk = jindcpu[3*j];
			jblk = jindcpu[3*j+1];
			ind = jindcpu[3*j+2];
			ObjTypeSend[j] = iblk;
			ObjIDSend[j] = jblk;
			CpuIDSend[j] = i;
			amatrarr[ind].PackMatrix (ObjSizeSend[j], ObjSend[j]);
		};
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
	if(info) throw " CSMatrix::ComputeSetOfSchurComplements: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for (i=0;i<NObjSend;i++) delete [] ObjSend[i];
	delete [] ObjSend;

// Store received data as submatrices

	int *imaskblkt;

	imaskblkt = new int [_nblks];
	if(!imaskblkt) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) imaskblkt[i] = 0;

	for (i=0;i<NObjRecv;i++) {
		jblk = ObjIDRecv[i];
		imaskblkt[jblk]++;
	};

	int *listblkt;
	int *iablkt;

	listblkt = new int [_nblks];
	if(!listblkt) MemoryFail(funcname);
	iablkt = new int [_nblks+1];
	if(!iablkt) MemoryFail(funcname);

	int nlistblkt = 0;

	iablkt[0] = 0;

	for (i=0;i<_nblks;i++) {
		if (imaskblkt[i] > 0) {
			listblkt[nlistblkt] = i;
			iablkt[nlistblkt+1] = iablkt[nlistblkt]+imaskblkt[i];
			nlistblkt++;
		};
	};

	for (i=0;i<_nblks;i++) imaskblkt[i] = -1;

	for (i=0;i<nlistblkt;i++) {
		jj = listblkt[i];
		imaskblkt[jj] = i;
	};

	int nzjablkt = iablkt[nlistblkt];

	int *iptrblkt;
	int *jablkt;
	CSMatrix *amtrblkt;

	iptrblkt = new int [nlistblkt];
	if(!iptrblkt) MemoryFail(funcname);
	jablkt = new int [nzjablkt];
	if(!jablkt) MemoryFail(funcname);
	amtrblkt = new CSMatrix [nzjablkt];
	if(!amtrblkt) MemoryFail(funcname);

	for (i=0;i<nlistblkt;i++) iptrblkt[i] = iablkt[i];

	for (i=0;i<NObjRecv;i++) {
		iblk = ObjTypeRecv[i];
		jblk = ObjIDRecv[i];
		ind = imaskblkt[jblk];
		k = iptrblkt[ind];
		jablkt[k] = iblk;
		iptrblkt[ind]++;
	};

	for (i=0;i<nlistblkt;i++) {
		ibeg = iablkt[i];
		ni = iablkt[i+1]-iablkt[i];
		sort (jablkt+ibeg, jablkt+ibeg + ni);
	};

	for (i=0;i<NObjRecv;i++) {
		iblk = ObjTypeRecv[i];
		jblk = ObjIDRecv[i];
		ind = imaskblkt[jblk];
		k = -1;
		for (j=iablkt[ind];j<iablkt[ind+1];j++) {
			if (jablkt[j] == iblk) {
				k = j;
				break;
			};
		};
		if (k < 0) throw " CSMatrix::ComputeSetOfSchurComplements: the element not found ";
		amtrblkt[k].UnPackMatrix (ObjTypeRecv[i], ObjRecv[i]);
	};

// Free packed received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for (i=0;i<NObjRecv;i++) delete [] ObjRecv[i];
	delete [] ObjRecv;

// Prepare the multiplications block sparsity

	int *piatablk = atablkmtr.GetIa ();
	int *pjatablk = atablkmtr.GetJa ();

	int *iamult;

	iamult = new int [nlistblk+1];
	if (!iamult) MemoryFail (funcname);

	iamult[0] = 0;
	int nz = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piatablk[iblk];j<piatablk[iblk+1];j++) {
			jblk = pjatablk[j];
			if (jblk > iblk) nz++;
		};
		iamult[i+1] = nz;
	};

	int *jamult;
	CSMatrix *amultarr;
	CSMatrix *amultblkarr;

	int nzjamult = nz;

	jamult = new int [nzjamult];
	if (!jamult) MemoryFail (funcname);
	amultarr = new CSMatrix [nzjamult];
	if (!amultarr) MemoryFail (funcname);
	amultblkarr = new CSMatrix [nzjamult];
	if (!amultblkarr) MemoryFail (funcname);

	nz = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piatablk[iblk];j<piatablk[iblk+1];j++) {
			jblk = pjatablk[j];
			if (jblk > iblk) {
				jamult[nz] = jblk;
				nz++;
			};
		};
	};

// Perform all necessary local multiplications

	CSMatrix mtrmult, mtrblkmult, aprod;

	int ibeg1, iend1, ibeg2, iend2, ip1, ip2, kblk1, kblk2;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		for (j=iamult[ilistblk];j<iamult[ilistblk+1];j++) {
			jblk = jamult[j];
			ind = imaskblkt[jblk];
			if (ind >= 0) {
				ibeg1 = iablk[ilistblk];
				iend1 = iablk[ilistblk+1]-1;
				ibeg2 = iablkt[ind];
				iend2 = iablkt[ind+1]-1;
				ip1 = ibeg1;
				ip2 = ibeg2;
				mtrmult = mtrdummy;
				mtrblkmult = mtrdummy;
				while (ip1 <= iend1 && ip2 <= iend2) {
					kblk1 = jablk[ip1];
					kblk2 = jablkt[ip2];
					if (kblk1 == kblk2) {
						if (kblk1 < jblk) {
							aprod = amatrarr[ip1].SuperSparseMultiplyByRows (icycle, imask, 
																								icycle1, imask1, 
																								amtrblkt[ip2]);
							int nzjaprod = aprod.GetNzja ();
							CSMatrix aprodstr (nzjaprod, nzjaprod);
							int *piaprodstr = aprodstr.GetIa ();
							int *pjaprodstr = aprodstr.GetJa ();
							piaprodstr[0] = 0;
							for (k=0;k<nzjaprod;k++) {
								piaprodstr[k+1] = piaprodstr[k]+1;
								pjaprodstr[k] = iblk;
							};
							aprod.SuperSparseAddWithBlockList (icycle, imask, icycle1, imask1,
																			aprodstr,
																			mtrmult, mtrblkmult);
							aprod = mtrdummy;
						};
						ip1++;
						ip2++;
					} else if (kblk1 < kblk2) {
						ip1++;
					} else if (kblk1 > kblk2) {
						ip2++;
					};
				};
				amultarr[j] = mtrmult;
				amultblkarr[j] = mtrblkmult;
			};
		};
	};

// Free transposed matrices

	for (i=0;i<nzjablkt;i++) amtrblkt[i] = mtrdummy;

// Exchange the results of multiplications

	NObjSend = nzjamult;

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

	char *p1matr, *p2matr, *ptot;
	int *pint;
	int isize1, isize2;

	for (ilistblk=0;ilistblk<nlistblk;ilistblk++) {
		iblk = listblk[ilistblk];
		for (j=iamult[ilistblk];j<iamult[ilistblk+1];j++) {
			jblk = jamult[j];
			icpu = _blk2cpu[jblk];
			ObjTypeSend[j] = iblk;
			ObjIDSend[j] = jblk;
			CpuIDSend[j] = icpu;
			amultarr[j].PackMatrix (isize1, p1matr);
			amultblkarr[j].PackMatrix (isize2, p2matr);
			amultarr[j] = mtrdummy;
			amultblkarr[j] = mtrdummy;
			ObjSizeSend[j] = isize1+isize2+2*sizeof(int);
			ObjSend[j] = new char [ObjSizeSend[j]];
			if(!ObjSend[j]) MemoryFail(funcname);
			ptot = ObjSend[j];
			pint = (int *)ptot;
			pint[0] = isize1;
			pint[1] = isize2;
			ptot += 2*sizeof(int);
			memcpy (ptot, p1matr, isize1);
			ptot += isize1;
			memcpy (ptot, p2matr, isize2);
			delete [] p1matr;
			delete [] p2matr;
		};
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::ComputeSetOfSchurComplements: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for (i=0;i<NObjSend;i++) delete [] ObjSend[i];
	delete [] ObjSend;

// Store received data as submatrices

	for (i=0;i<_nblks;i++) imaskblkt[i] = 0;

	for (i=0;i<NObjRecv;i++) {
		jblk = ObjIDRecv[i];
		imaskblkt[jblk]++;
	};

	nlistblkt = 0;

	iablkt[0] = 0;

	for (i=0;i<_nblks;i++) {
		if (imaskblkt[i] > 0) {
			listblkt[nlistblkt] = i;
			iablkt[nlistblkt+1] = iablkt[nlistblkt]+imaskblkt[i];
			nlistblkt++;
		};
	};

	for (i=0;i<_nblks;i++) imaskblkt[i] = -1;

	for (i=0;i<nlistblkt;i++) {
		jj = listblkt[i];
		imaskblkt[jj] = i;
	};

	nzjablkt = iablkt[nlistblkt];

	delete [] iptrblkt;
	delete [] jablkt;
	delete [] amtrblkt;

	CSMatrix *amtrblkstrt;

	iptrblkt = new int [nlistblkt];
	if(!iptrblkt) MemoryFail(funcname);
	jablkt = new int [nzjablkt];
	if(!jablkt) MemoryFail(funcname);
	amtrblkt = new CSMatrix [nzjablkt];
	if(!amtrblkt) MemoryFail(funcname);
	amtrblkstrt = new CSMatrix [nzjablkt];
	if(!amtrblkstrt) MemoryFail(funcname);

	for (i=0;i<nlistblkt;i++) iptrblkt[i] = iablkt[i];

	for (i=0;i<NObjRecv;i++) {
		iblk = ObjTypeRecv[i];
		jblk = ObjIDRecv[i];
		ind = imaskblkt[jblk];
		k = iptrblkt[ind];
		jablkt[k] = iblk;
		iptrblkt[ind]++;
	};

	for (i=0;i<nlistblkt;i++) {
		ibeg = iablkt[i];
		ni = iablkt[i+1]-iablkt[i];
		sort (jablkt+ibeg, jablkt+ibeg + ni);
	};

	for (i=0;i<NObjRecv;i++) {
		iblk = ObjTypeRecv[i];
		jblk = ObjIDRecv[i];
		ind = imaskblkt[jblk];
		k = -1;
		for (j=iablkt[ind];j<iablkt[ind+1];j++) {
			if (jablkt[j] == iblk) {
				k = j;
				break;
			};
		};
		if (k < 0) throw " CSMatrix::ComputeSetOfSchurComplements: the element not found ";
		ptot = ObjRecv[i];
		pint = (int *)ptot;
		isize1 = pint[0];
		isize2 = pint[1];
		ptot += 2*sizeof(int);
		amtrblkt[k].UnPackMatrix (isize1, ptot);
		ptot += isize1;
		amtrblkstrt[k].UnPackMatrix (isize2, ptot);
	};

// Free packed received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for (i=0;i<NObjRecv;i++) delete [] ObjRecv[i];
	delete [] ObjRecv;

// Create sparsity pointers to the transposed block structure up to the diagonal

	int *iablktini;
	int *iptrblktini;

	iablktini = new int [_nblks+1];
	if(!iablktini) MemoryFail(funcname);
	iptrblktini = new int [_nblks];
	if(!iptrblktini) MemoryFail(funcname);

	for (i=0;i<=_nblks;i++) iablktini[i] = 0;

	int icheckdiag;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		icheckdiag = 0;
		for (j=iablk[i];j<iablk[i+1];j++) {
			jblk = jablk[j];
			if (jblk <= iblk) iablktini[iblk+1]++;
			if (jblk == iblk) icheckdiag = 1;
		};
		if (icheckdiag == 0) iablktini[iblk+1]++;
	};

	for (i=0;i<_nblks;i++) iablktini[i+1] = iablktini[i]+iablktini[i+1];

	for (i=0;i<_nblks;i++) iptrblktini[i] = iablktini[i];

	int nzjablktini = iablktini[_nblks];

	int *jablktini;
	CSMatrix **pablktini;

	jablktini = new int [nzjablktini];
	if(!jablktini) MemoryFail(funcname);
	pablktini = new CSMatrix * [nzjablktini];
	if(!pablktini) MemoryFail(funcname);

	CSMatrix mtrdiag;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		icheckdiag = 0;
		for (j=iablk[i];j<iablk[i+1];j++) {
			jblk = jablk[j];
			if (jblk <= iblk) {
				k = iptrblktini[iblk];
				jablktini[k] = jblk;
				pablktini[k] = amatrarr+j;
				iptrblktini[iblk]++;
			};
			if (jblk == iblk) icheckdiag = 1;
		};
		if (icheckdiag == 0) {
			k = iptrblktini[iblk];
			jablktini[k] = iblk;
			pablktini[k] = &mtrdiag;
			iptrblktini[iblk]++;
		};
	};

// Add parts of the Schur complement

	CSMatrix aprodstr;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_imaskblk[iblk] >= 0 && _blk2cpu[iblk] == myid) {

// Init by the diagonal block data

			ind = iablktini[iblk+1]-1;

			if (jablktini[ind] != iblk) {
				throw " CSMatrix::ComputeSetOfSchurComplements: the diagonal element not found ";
			};

			mtrmult = *pablktini[ind];

			int nzjadia = mtrmult.GetNzja ();

			CSMatrix mtrstr (nzjadia,0);

			int *piastr = mtrstr.GetIa ();

			for (i=0;i<=nzjadia;i++) piastr[i] = 0;

// Multiply and add

			ilistblk = imaskblkt[iblk];

			if (ilistblk >= 0) {

				ibeg1 = iablktini[iblk];
				iend1 = iablktini[iblk+1]-1;
				ibeg2 = iablkt[ilistblk];
				iend2 = iablkt[ilistblk+1]-1;
				ip1 = ibeg1;
				ip2 = ibeg2;
				while (ip1 <= iend1 && ip2 <= iend2) {
					kblk1 = jablktini[ip1];
					kblk2 = jablkt[ip2];
					if (kblk1 == kblk2) {
						if (kblk1 < iblk) {

							pablktini[ip1]->SuperSparseMultiplyByRowsWithBlockList (
																		icycle, imask, 
																		icycle1, imask1, 
																		amtrblkt[ip2], 
																		kblk1, amtrblkstrt[ip2],
																		aprod, aprodstr);

							aprod.SuperSparseAddWithBlockList (icycle, imask, icycle1, imask1,
																			aprodstr,
																			mtrmult, mtrstr);
							aprod = mtrdummy;
							aprodstr = mtrdummy;
						};
						ip1++;
						ip2++;
					} else if (kblk1 < kblk2) {
						ip1++;
					} else if (kblk1 > kblk2) {
						ip2++;
					};
				};

			};

// Perform final filtering of the local block description

			int nlistmlt = mtrmult.GetNlist ();
			int *plistmlt = mtrmult.GetList ();
			int *piamlt = mtrmult.GetIa ();
			int *pjamlt = mtrmult.GetJa ();
			int *piamltstr = mtrstr.GetIa ();
			int *pjamltstr = mtrstr.GetJa ();

			int nzjanew = 0;

			for (i=0;i<nlistmlt;i++) {
				irow = plistmlt[i];
				for (j=piamlt[i];j<piamlt[i+1];j++) {
					jj = pjamlt[j];
					if (irow == jj) {
						nzjanew += piamltstr[j+1]-piamltstr[j];
					};
				};
			};

			CSMatrix mtrstrnew (nlistmlt,nzjanew);

			int *pianew = mtrstrnew.GetIa ();
			int *pjanew = mtrstrnew.GetJa ();

			pianew[0] = 0;
			nzjanew = 0;

			for (i=0;i<nlistmlt;i++) {
				irow = plistmlt[i];
				for (j=piamlt[i];j<piamlt[i+1];j++) {
					jj = pjamlt[j];
					if (irow == jj) {
						for (k=piamltstr[j];k<piamltstr[j+1];k++) {
							pjanew[nzjanew] = pjamltstr[k];
							nzjanew++;
						};
						pianew[i+1] = nzjanew;
					};
				};
			};

			mtrstrnew.m = nlistmlt;
			mtrstrnew.n = nlistmlt;
			mtrstrnew.nsupc = nlistmlt;
			mtrstrnew.nsupr = nlistmlt;
			mtrstrnew.nlist = nlistmlt;
			mtrstrnew.nzja = nzjanew;

// Store the result

			_mtrarr[iblk] = mtrmult;
			_mtrblkarr[iblk] = mtrstrnew;

		};
	};

// Free work arrays

	delete [] list2blk;
	delete [] ja2blk;
	delete [] imaskblk;
	delete [] listblk;
	delete [] ialist2blk;
	delete [] iablk;
	delete [] iptrblk;
	delete [] jablk;
	delete [] ibsmask;
	delete [] imask;
	delete [] imask1;
	delete [] listarr;
	delete [] iptrarr;
	delete [] niarr;
	delete [] nzjaarr;
	delete [] amatrarr;
	delete [] amatrrowarr;
	delete [] amatrblkrowarr;
	delete [] iacpu;
	delete [] iptrcpu;
	delete [] imaskcpu;
	delete [] jindcpu;
	delete [] imaskblkt;
	delete [] listblkt;
	delete [] iablkt;
	delete [] iptrblkt;
	delete [] jablkt;
	delete [] amtrblkt;
	delete [] iamult;
	delete [] jamult;
	delete [] amultarr;
	delete [] amultblkarr;
	delete [] iablktini;
	delete [] iptrblktini;
	delete [] jablktini;
	delete [] pablktini;
	delete [] amtrblkstrt;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute in parallel the sparsity of the extended boundary
//========================================================================================
CSMatrix CSMatrix::ExtendGroupsBoundary (CMPIComm &_comm, int _ncycle, // Compute in parallel the sparsity of the extended boundary
														int _nblks, int *_blks, int *_blk2cpu, int *_partitionblk) const {

	const char *funcname = "ExtendGroupsBoundary";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	if (_ncycle != 1 && _ncycle != 2) throw " CSMatrix::ExtendGroupsBoundary: incorrect parameter ncycle ";

// Compute local base

	int *ibsblklist;

	ibsblklist = new int [_nblks];
	if (!ibsblklist) MemoryFail (funcname);

	int i;

	for (i=0;i<_nblks;i++) ibsblklist[i] = -1;

	int nz = 0;

	for (i=0;i<_nblks;i++) {
		if (_blk2cpu[i] == myid) {
			ibsblklist[i] = nz;
			nz += (_blks[i+1]-_blks[i]);
		};
	};

// Check partitioning index for own data

	int *list2blk;

	list2blk = new int [nlist];
	if (!list2blk) MemoryFail (funcname);

	int jj, iblk, ibegblk, iendblk, icpu, j;

	int iblkprev = 0;

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

// For each element find its block number

	int *ja2blk;

	ja2blk = new int [nzja];
	if (!ja2blk) MemoryFail (funcname);

/*
	iblkprev = 0;

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

// Count the number of local elements to be sent for each cpu

	int *iacpu;
	int *iptr;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int ni, ibs, jblk, ilist;

	for (iblk=0;iblk<_nblks;iblk++) {
		ni = _blks[iblk+1]-_blks[iblk];
		if (_blk2cpu[iblk] == myid && _partitionblk[iblk] == -1) {
			ibs = ibsblklist[iblk];
			for (i=0;i<ni;i++) {
				ilist = ibs+i;
				for (j=ia[ilist];j<ia[ilist+1];j++) {
					jblk = ja2blk[j];
					icpu = _blk2cpu[jblk];
					iacpu[icpu+1]++;
				};
			};
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	nz = iacpu[nproc];

	int *isendarr;

	isendarr = new int [3*nz];
	if (!isendarr) MemoryFail (funcname);

	int k;

	for (iblk=0;iblk<_nblks;iblk++) {
		ni = _blks[iblk+1]-_blks[iblk];
		if (_blk2cpu[iblk] == myid && _partitionblk[iblk] == -1) {
			ibs = ibsblklist[iblk];
			for (i=0;i<ni;i++) {
				ilist = ibs+i;
				for (j=ia[ilist];j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2blk[j];
					icpu = _blk2cpu[jblk];
					k = iptr[icpu];
					isendarr[k*3] = ilist;
					isendarr[k*3+1] = jj;
					isendarr[k*3+2] = jblk;
					iptr[icpu]++;
				};
			};
		};
	};

// Compute the set of submatrices

	CSMatrix *mtrarr;

	mtrarr = new CSMatrix [nproc];
	if (!mtrarr) MemoryFail (funcname);

	int *imask;
	int *listloc;
	int *list2loc;
	int *col2list;
	int *iptrloc;

	imask = new int [nlist];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlist];
	if (!list2loc) MemoryFail (funcname);
	col2list = new int [nlist];
	if (!col2list) MemoryFail (funcname);
	iptrloc = new int [nlist];
	if (!iptrloc) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlist;i++) imask[i] = icycle;

	int nlistloc, nzjaloc, ind;

	for (icpu=0;icpu<nproc;icpu++) {
		icycle++;
		nlistloc = 0;
		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ilist = isendarr[j*3];
			if (imask[ilist] != icycle) {
				imask[ilist] = icycle;
				listloc[nlistloc] = ilist;
				nlistloc++;
			};
		};
		sort (listloc, listloc + nlistloc);
		for (i=0;i<nlistloc;i++) {
			ilist = listloc[i];
			col2list[ilist] = i;
		};
		nzjaloc = iacpu[icpu+1]-iacpu[icpu];
		CSMatrix temp (nlistloc,nlistloc,nzjaloc,nzjaloc);
		for (i=0;i<=nlistloc;i++) temp.ia[i] = 0;
		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ilist = isendarr[j*3];
			ind = col2list[ilist];
			temp.ia[ind+1]++;
		};
		for (i=0;i<nlistloc;i++) temp.ia[i+1] = temp.ia[i]+temp.ia[i+1];
		for (i=0;i<nlistloc;i++) iptrloc[i] = temp.ia[i];
		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ilist = isendarr[j*3];
			jj = isendarr[j*3+1];
			jblk = isendarr[j*3+2];
			ind = col2list[ilist];
			k = iptrloc[ind];
			temp.ja[k] = jj;
			temp.ja2[k] = jblk;
			iptrloc[ind]++;
		};
		for (i=0;i<nlistloc;i++) {
			ilist = listloc[i];
			iblk = list2blk[ilist];
			temp.list[i] = list[ilist];
			temp.list2[i] = iblk;
		};
		temp.m = m;
		temp.n = n;
		temp.nsupc = nsupc;
		temp.nsupr = nsupr;
		temp.nlist = nlistloc;
		temp.nlist2 = nlistloc;
		temp.nzja = nzjaloc;
		temp.nzja2 = nzjaloc;
	};

	delete [] imask;
	delete [] listloc;
	delete [] list2loc;
	delete [] col2list;
	delete [] iptrloc;

	delete [] isendarr;

// Exchange the matrices

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

	CSMatrix mtrdummy;

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		mtrarr[i].PackMatrix (ObjSizeSend[i], ObjSend[i]);
		mtrarr[i] = mtrdummy;
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
	if(info) throw " CSMatrix::ExtendGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for (i=0;i<nproc;i++) delete [] ObjSend[i];
	delete [] ObjSend;

// Unpack received data

	for(i = 0; i < NObjRecv; i++) {
		icpu = CpuIDRecv[i];
		mtrarr[icpu].UnPackMatrix (ObjSizeSend[i], ObjSend[i]);
	};

// Free receive structures

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Modify indices to create the set of two index matrices

	int nlistcpu;
	int *plist, *plist2;
	int jjnew;

	for (icpu=0;icpu<nproc;icpu++) {
		nlistcpu = mtrarr[icpu].nlist;
		plist = mtrarr[icpu].list;
		plist2 = mtrarr[icpu].list2;
		for (i=0;i<nlistcpu;i++) {
			jj = plist[i];
			jblk = plist2[i];
			jjnew = jj-_blks[jblk];
			plist[i] = jjnew;
		};
		nlistcpu = mtrarr[icpu].nzja;
		plist = mtrarr[icpu].ja;
		plist2 = mtrarr[icpu].ja2;
		for (i=0;i<nlistcpu;i++) {
			jj = plist[i];
			jblk = plist2[i];
			jjnew = jj-_blks[jblk];
			plist[i] = jjnew;
		};
	};

// Combine received matrices into one

	int *listcpu;

	listcpu = new int [nproc];
	if(!listcpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) listcpu[i] = i;

	CSMatrix mtratot;

	CGSMatrix::CombineMatricesSort2Index (n, nproc, listcpu,
														mtrarr,
														mtratot);

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

// Modify indices back in the result

	nlistcpu = mtratot.nlist;
	plist = mtratot.list;
	plist2 = mtratot.list2;
	for (i=0;i<nlistcpu;i++) {
		jj = plist[i];
		jblk = plist2[i];
		jjnew = jj+_blks[jblk];
		plist[i] = jjnew;
	};
	nlistcpu = mtratot.nzja;
	plist = mtratot.ja;
	plist2 = mtratot.ja2;
	for (i=0;i<nlistcpu;i++) {
		jj = plist[i];
		jblk = plist2[i];
		jjnew = jj+_blks[jblk];
		plist[i] = jjnew;
	};

// Prepare multiplication and Schur computations block mask

	int *imaskblk;

	imaskblk = new int [_nblks];
	if (!imaskblk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) imaskblk[i] = -1;

	for (i=0;i<nzja;i++) {
		jblk = ja2blk[i];
		if (_partitionblk[jblk] != -1) imaskblk[jblk] = 0;
	};

	for (i=0;i<_nblks;i++) {
		if (_partitionblk[jblk] == -1) imaskblk[i] = 0;
	};

	int *ibsblk;

	ibsblk = new int [_nblks];
	if (!ibsblk) MemoryFail (funcname);

	nz = 0;

	for (i=0;i<_nblks;i++) {
		if (imaskblk[i] == 0) {
			ibsblk[i] = nz;
			nz += _blks[i+1]-_blks[i];
		};
	};

	int *imasknd;

	imasknd = new int [nz];
	if (!imasknd) MemoryFail (funcname);

	icycle = -1;

	for (i=0;i<nz;i++) imasknd[i] = icycle;

// Perform multiplication

	nlistloc = mtratot.nlist;
	plist = mtratot.list;
	plist2 = mtratot.list2;

	int *pia = mtratot.ia;
	int *pja = mtratot.ja;
	int *pja2 = mtratot.ja2;

	int nztot = 0;

	int jlist, kk, kblk, kind, kbs;

	for (i=0;i<nlistloc;i++) {
		icycle++;
		for (j=pia[i];j<pia[i+1];j++) {
			jj = pja[j];
			jblk = pja2[j];
			ibs = ibsblk[jblk];
			jlist = ibs+jj-_blks[jblk];
			for (k=ia[jlist];k<ia[jlist+1];k++) {
				kk = ja[k];
				kblk = ja2blk[k];
				kbs = ibsblk[kblk];
				if (kbs >= 0) {
					kind = kbs+kk-_blks[kblk];
					if (imasknd[kind] != icycle) {
						nztot++;
						imasknd[kind] = icycle;
					};
				};
			};
		};
	};

	CSMatrix aloc (nlistloc,nlistloc,nztot,nztot);

	nztot = 0;
	aloc.ia[0] = 0;

	for (i=0;i<nlistloc;i++) {
		aloc.list[i] = plist[i];
		aloc.list2[i] = plist2[i];
		icycle++;
		for (j=pia[i];j<pia[i+1];j++) {
			jj = pja[j];
			jblk = pja2[j];
			ibs = ibsblk[jblk];
			jlist = ibs+jj-_blks[jblk];
			for (k=ia[jlist];k<ia[jlist+1];k++) {
				kk = ja[k];
				kblk = ja2blk[k];
				kbs = ibsblk[kblk];
				if (kbs >= 0) {
					kind = kbs+kk-_blks[kblk];
					if (imasknd[kind] != icycle) {
						aloc.ja[nztot] = kk;
						aloc.ja2[nztot] = kblk;
						nztot++;
						imasknd[kind] = icycle;
					};
				};
			};
		};
		aloc.ia[nlist+1] = nztot;
	};

	aloc.m = m;
	aloc.n = n;
	aloc.nsupc = nsupc;
	aloc.nsupr = nsupr;
	aloc.nlist = nlistloc;
	aloc.nlist2 = nlistloc;
	aloc.nzja = nztot;
	aloc.nzja2 = nztot;

	CSMatrix mtramult = aloc;

	aloc = mtrdummy;

// Split the matrix into the set of submatrices

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	for (i=0;i<mtramult.nlist;i++) {
		for (j=mtramult.ia[i];j<mtramult.ia[i+1];j++) {
			jblk = mtramult.ja2[j];
			icpu = _blk2cpu[jblk];
			iacpu[icpu+1]++;
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	nz = iacpu[nproc];

	isendarr = new int [nz*4];
	if (!isendarr) MemoryFail (funcname);

	int irow;

	for (i=0;i<mtramult.nlist;i++) {
		irow = mtramult.list[i];
		iblk = mtramult.list2[i];
		for (j=mtramult.ia[i];j<mtramult.ia[i+1];j++) {
			jj = mtramult.ja[j];
			jblk = mtramult.ja2[j];
			icpu = _blk2cpu[jblk];
			k = iptr[icpu];
			isendarr[k*4  ] = irow;
			isendarr[k*4+1] = iblk;
			isendarr[k*4+2] = jj;
			isendarr[k*4+3] = jblk;
			iptr[icpu]++;
		};
	};

	mtramult = mtrdummy;

// Exchange the data

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
		ObjSend[i] = (char *)(isendarr+iacpu[i]*4);
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::ExtendGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	delete [] isendarr;

// Combine received data into one matrix

	nztot = 0;

	for (i=0;i<NObjRecv;i++) {
		nztot += ObjSizeRecv[i] / (4*sizeof(int));
	};

	isendarr = new int [nztot*4];
	if (!isendarr) MemoryFail (funcname);

	nztot = 0;

	int *piarr;
	int nzloc;

	for (i=0;i<NObjRecv;i++) {
		piarr = (int *)ObjRecv[i];
		nzloc = ObjSizeRecv[i] / (4*sizeof(int));
		for (j=0;j<nzloc;j++) {
			isendarr[nztot*4] = piarr[j*4];
			isendarr[nztot*4+1] = piarr[j*4+1];
			isendarr[nztot*4+2] = piarr[j*4+2];
			isendarr[nztot*4+3] = piarr[j*4+3];
			nztot++;
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

	CIntInt *iiarr;

	iiarr = new CIntInt [nztot];
	if(!iiarr) MemoryFail(funcname);

	for (j=0;j<nztot;j++) {
		iiarr[j].intvalue = isendarr[j*4];
		iiarr[j].int2value = j;
	};

	//qsort (iiarr, nztot, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + nztot);

	int *isendarr1;

	isendarr1 = new int [nztot*4];
	if (!isendarr1) MemoryFail (funcname);

	for (i=0;i<nztot;i++) {
		ind = iiarr[i].int2value;
		isendarr1[i*4] = isendarr[ind*4];
		isendarr1[i*4+1] = isendarr[ind*4+1];
		isendarr1[i*4+2] = isendarr[ind*4+2];
		isendarr1[i*4+3] = isendarr[ind*4+3];
	};

	delete [] iiarr;
	delete [] isendarr;

	isendarr = isendarr1;

	int *ialist;

	ialist = new int [nztot+1];
	if (!ialist) MemoryFail (funcname);

	ialist[0] = 0;

	nlistloc = 0;

	j = 0;
	while (j<nztot) {
		if (j>0) {
			if (isendarr[(j-1)*4] != isendarr[j*4]) {
				nlistloc++;
			};
		};
		j++;
		ialist[nlistloc+1] = j;
	};
	nlistloc++;

	listloc = new int [nlistloc];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistloc];
	if (!list2loc) MemoryFail (funcname);

	for (i=0;i<nlistloc;i++) {
		ind = ialist[i];
		listloc[i] = isendarr[ind*4];
		list2loc[i] = isendarr[ind*4+1];
	};

	int nzmax = 0;

	for (i=0;i<nlistloc;i++) {
		nzloc = ialist[i+1]-ialist[i];
		if (nzloc > nzmax) nzmax = nzloc;
	};

	iiarr = new CIntInt [nzmax];
	if(!iiarr) MemoryFail(funcname);

	int *ialoc;
	int *jaloc;
	int *ja2loc;

	ialoc = new int [nlistloc+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nztot];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nztot];
	if (!ja2loc) MemoryFail (funcname);

	ialoc[0] = 0;
	nz = 0;

	for (i=0;i<nlistloc;i++) {
		nzloc = 0;
		for (j=ialist[i];j<ialist[i+1];j++) {
			iiarr[nzloc].intvalue = isendarr[j*4+2];
			iiarr[nzloc].int2value = nzloc;
			nzloc++;
		};
		//qsort (iiarr, nzloc, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nzloc);
		for (j=0;j<nzloc;j++) {
			ind = ialist[i]+iiarr[j].int2value;
			jj = isendarr[ind*4+2];
			jblk = isendarr[ind*4+3];
			if (j == 0) {
				jaloc[nz] = jj;
				ja2loc[nz] = jblk;
				nz++;
			} else if (jj != jaloc[nz-1]) {
				jaloc[nz] = jj;
				ja2loc[nz] = jblk;
				nz++;
			};
		};
		ialoc[i+1] = nz;
	};

	delete [] iiarr;

	delete [] mtramult.list;
	delete [] mtramult.list2;
	delete [] mtramult.ia;
	delete [] mtramult.ja;
	delete [] mtramult.ja2;

	mtramult.list = listloc;
	mtramult.list2 = list2loc;
	mtramult.ia = ialoc;
	mtramult.ja = jaloc;
	mtramult.ja2 = ja2loc;

	mtramult.m = m;
	mtramult.n = n;
	mtramult.nsupc = nsupc;
	mtramult.nsupr = nsupr;
	mtramult.nlist = nlistloc;
	mtramult.nlist2 = nlistloc;
	mtramult.nzja = nz;
	mtramult.nzja2 = nz;

// Compute local Schur complement

// Open submatrices

	int nlistL, *plistL, *plist2L, *piaL, *pjaL, *pja2L;
	int nlistR, *plistR, *plist2R, *piaR, *pjaR, *pja2R;

	if (_ncycle == 1) {
		nlistL = mtratot.nlist;
		plistL = mtratot.list;
		plist2L = mtratot.list2;
		piaL = mtratot.ia;
		pjaL = mtratot.ja;
		pja2L = mtratot.ja2;
	} else {
		nlistL = mtramult.nlist;
		plistL = mtramult.list;
		plist2L = mtramult.list2;
		piaL = mtramult.ia;
		pjaL = mtramult.ja;
		pja2L = mtramult.ja2;
	};
	nlistR = mtramult.nlist;
	plistR = mtramult.list;
	plist2R = mtramult.list2;
	piaR = mtramult.ia;
	pjaR = mtramult.ja;
	pja2R = mtramult.ja2;

// Transpose right matrix

	int nlistrmax = nlistloc;

	int *listRT;
	int *list2RT;

	listRT = new int [nlistrmax];
	if (!listRT) MemoryFail (funcname);
	list2RT = new int [nlistrmax];
	if (!list2RT) MemoryFail (funcname);

	icycle++;

	int nlistRT = 0;

	int jind;

	for (i=0;i<nlistR;i++) {
		irow = plistR[i];
		iblk = plist2R[i];
		for (j=piaR[i];j<piaR[i+1];j++) {
			jj = pjaR[j];
			jblk = pja2R[j];
			ibs = ibsblk[jblk];
			jind = jj-_blks[jblk]+ibs;
			if (imasknd[jind] != icycle) {
				listRT[nlistRT] = jj;
				list2RT[nlistRT] = jblk;
				imasknd[jind] = icycle;
			};
		};
	};

	iiarr = new CIntInt [nlistRT];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nlistRT;i++) {
		iiarr[i].intvalue = listRT[i];
		iiarr[i].int2value = i;
	};

	//qsort (iiarr, nlistRT, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + nlistRT);

	int *listRTn, *list2RTn;

	listRTn = new int [nlistRT];
	if (!listRTn) MemoryFail (funcname);
	list2RTn = new int [nlistRT];
	if (!list2RTn) MemoryFail (funcname);

	for (i=0;i<nlistRT;i++) {
		ind = iiarr[i].int2value;
		listRTn[i] = listRT[ind];
		list2RTn[i] = list2RT[ind];
	};

	delete [] iiarr;
	delete [] listRT;
	delete [] list2RT;

	listRT = listRTn;
	list2RT = list2RTn;

	int *nd2ind;

	nd2ind = new int [nlistrmax];
	if (!nd2ind) MemoryFail (funcname);

	for (i=0;i<nlistrmax;i++) nd2ind[i] = -1;

	for (i=0;i<nlistRT;i++) {
		jj = listRT[i];
		jblk = list2RT[i];
		ibs = ibsblk[jblk];
		jind = jj-_blks[jblk]+ibs;
		nd2ind[jind] = i;
	};

	int *iaRT;
	int *iptrRT;

	iaRT = new int [nlistRT+1];
	if (!iaRT) MemoryFail (funcname);
	iptrRT = new int [nlistRT];
	if (!iptrRT) MemoryFail (funcname);

	for (i=0;i<=nlistRT;i++) iaRT[i] = 0;

	int indt;

	for (i=0;i<nlistR;i++) {
		irow = plistR[i];
		iblk = plist2R[i];
		for (j=piaR[i];j<piaR[i+1];j++) {
			jj = pjaR[j];
			jblk = pja2R[j];
			ibs = ibsblk[jblk];
			jind = jj-_blks[jblk]+ibs;
			indt = nd2ind[jind];
			iaRT[indt+1]++;
		};
	};

	for (i=0;i<nlistRT;i++) iaRT[i+1] = iaRT[i]+iaRT[i+1];

	for (i=0;i<nlistRT;i++) iptrRT[i] = iaRT[i];

	int nzjaRT = iaRT[nlistRT];

	int *jaRT;
	int *ja2RT;

	jaRT = new int [nzjaRT];
	if (!jaRT) MemoryFail (funcname);
	ja2RT = new int [nzjaRT];
	if (!ja2RT) MemoryFail (funcname);

	for (i=0;i<nlistR;i++) {
		irow = plistR[i];
		iblk = plist2R[i];
		for (j=piaR[i];j<piaR[i+1];j++) {
			jj = pjaR[j];
			jblk = pja2R[j];
			ibs = ibsblk[jblk];
			jind = jj-_blks[jblk]+ibs;
			indt = nd2ind[jind];
			k = iptrRT[indt];
			jaRT[k] = irow;
			ja2RT[k] = iblk;
			iptrRT[indt]++;
		};
	};

// Multiply

	listloc = new int [nlistL];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistL];
	if (!list2loc) MemoryFail (funcname);
	iiarr = new CIntInt [nlistL];
	if (!iiarr) MemoryFail (funcname);

	int *iamult;

	iamult = new int [nlistL+1];
	if (!iamult) MemoryFail (funcname);

	iamult[0] = 0;
	int nzmult = 0;

	for (i=0;i<nlistL;i++) {
		icycle++;
		nlistloc = 0;
		for (j=piaL[i];j<piaL[i+1];j++) {
			jj = pjaL[j];
			jblk = pja2L[j];
			ibs = ibsblk[jblk];
			jind = jj-_blks[jblk]+ibs;
			indt = nd2ind[jind];
			if (indt >= 0) {
				for (k=iaRT[indt];k<iaRT[indt+1];k++) {
					kk = jaRT[k];
					kblk = ja2RT[k];
					kbs = ibsblk[jblk];
					kind = kk-_blks[kblk]+kbs;
					if (imasknd[kind] != icycle) {
						nlistloc++;
						imasknd[kind] = icycle;
					};
				};
			};
		};
		nzmult += nlistloc;
		iamult[i+1] = nzmult;
	};

	int *jamult;
	int *ja2mult;

	jamult = new int [nzmult];
	if (!jamult) MemoryFail (funcname);
	ja2mult = new int [nzmult];
	if (!ja2mult) MemoryFail (funcname);

	nzmult = 0;

	for (i=0;i<nlistL;i++) {
		icycle++;
		nlistloc = 0;
		for (j=piaL[i];j<piaL[i+1];j++) {
			jj = pjaL[j];
			jblk = pja2L[j];
			ibs = ibsblk[jblk];
			jind = jj-_blks[jblk]+ibs;
			indt = nd2ind[jind];
			if (indt >= 0) {
				for (k=iaRT[indt];k<iaRT[indt+1];k++) {
					kk = jaRT[k];
					kblk = ja2RT[k];
					kbs = ibsblk[jblk];
					kind = kk-_blks[kblk]+kbs;
					if (imasknd[kind] != icycle) {
						listloc[nlistloc] = kk;
						list2loc[nlistloc] = kblk;
						nlistloc++;
						imasknd[kind] = icycle;
					};
				};
			};
		};
		for (i=0;i<nlistloc;i++) {
			iiarr[i].intvalue = listloc[i];
			iiarr[i].int2value = i;
		};
		//qsort (iiarr, nlistloc, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nlistloc);
		for (i=0;i<nlistloc;i++) {
			ind = iiarr[i].int2value;
			jamult[nzmult+i] = listloc[ind];
			ja2mult[nzmult+i] = list2loc[ind];
		};
		nzmult += nlistloc;
	};

	delete [] iiarr;
	delete [] listloc;
	delete [] list2loc;

	int nlistmult = nlistL;

	int *listmult;
	int *list2mult;

	listmult = new int [nlistmult];
	if (!listmult) MemoryFail (funcname);
	list2mult = new int [nlistmult];
	if (!list2mult) MemoryFail (funcname);

	for (i=0;i<nlistmult;i++) listmult[i] = plistL[i];
	for (i=0;i<nlistmult;i++) list2mult[i] = plist2L[i];

	mtratot = mtrdummy;
	mtramult = mtrdummy;

// Split computed Schur complement into the set of submatrices

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	for (i=0;i<nlistmult;i++) {
		jblk = list2mult[i];
		icpu = _blk2cpu[jblk];
		iacpu[icpu+1] += iamult[i+1]-iamult[i];
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	nz = iacpu[nproc];

	isendarr = new int [nz*4];
	if (!isendarr) MemoryFail (funcname);

	for (i=0;i<nlistmult;i++) {
		irow = listmult[i];
		iblk = list2mult[i];
		icpu = _blk2cpu[iblk];
		for (j=iamult[i];j<iamult[i+1];j++) {
			jj = jamult[j];
			jblk = ja2mult[j];
			k = iptr[icpu];
			isendarr[k*4] = irow;
			isendarr[k*4+1] = iblk;
			isendarr[k*4+2] = jj;
			isendarr[k*4+3] = jblk;
		};
	};

// Exchange submatrices

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
		ObjSend[i] = (char *)(isendarr+iacpu[i]*4);
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::ExtendGroupsBoundary: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	delete [] isendarr;

// Combine received data into one array

	nztot = 0;

	for (i=0;i<NObjRecv;i++) {
		nztot += ObjSizeRecv[i] / (4*sizeof(int));
	};

	isendarr = new int [nztot*4];
	if (!isendarr) MemoryFail (funcname);

	nztot = 0;

	for (i=0;i<NObjRecv;i++) {
		piarr = (int *)ObjRecv[i];
		nzloc = ObjSizeRecv[i] / (4*sizeof(int));
		for (j=0;j<nzloc;j++) {
			isendarr[nztot*4] = piarr[j*4];
			isendarr[nztot*4+1] = piarr[j*4+1];
			isendarr[nztot*4+2] = piarr[j*4+2];
			isendarr[nztot*4+3] = piarr[j*4+3];
			nztot++;
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

// Add received data with the local Schur submatrix

	int nschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid && _partitionblk[iblk] == -1) {
			nschur += (_blks[iblk+1]-_blks[iblk]);
		};
	};

	int *listschur;
	int *list2schur;
	int *iaschur;
	int *iptrschur;

	listschur = new int [nschur];
	if (!listschur) MemoryFail (funcname);
	list2schur = new int [nschur];
	if (!list2schur) MemoryFail (funcname);
	iaschur = new int [nschur+1];
	if (!iaschur) MemoryFail (funcname);
	iptrschur = new int [nschur];
	if (!iptrschur) MemoryFail (funcname);

	nschur = 0;

	int iind;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid && _partitionblk[iblk] == -1) {
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				listschur[nschur] = i;
				list2schur[nschur] = iblk;
				ibs = ibsblk[iblk];
				iind = i-_blks[iblk]+ibs;
				nd2ind[iind] = nschur;
				nschur++;
			};
		};
	};

	for (i=0;i<=nschur;i++) iaschur[i] = 0;

	int ipos;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		if (_partitionblk[iblk] == -1) {
			ibs = ibsblk[iblk];
			iind = irow-_blks[iblk]+ibs;
			ipos = nd2ind[iind];
			iaschur[ipos+1] += (ia[i+1]-ia[i]);
		};
	};

	for (i=0;i<nztot;i++) {
		irow = isendarr[i*4];
		iblk = isendarr[i*4+1];
		ibs = ibsblk[iblk];
		iind = irow-_blks[iblk]+ibs;
		ipos = nd2ind[iind];
		iaschur[ipos+1]++;
	};

	for (i=0;i<nschur;i++) iaschur[i+1] = iaschur[i]+iaschur[i+1];

	for (i=0;i<nschur;i++) iptrschur[i] = iaschur[i];

	int nzschur = iaschur[nschur];

	int *jaschur;
	int *ja2schur;

	jaschur = new int [nzschur];
	if (!jaschur) MemoryFail (funcname);
	ja2schur = new int [nzschur];
	if (!ja2schur) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		if (_partitionblk[iblk] == -1) {
			ibs = ibsblk[iblk];
			iind = irow-_blks[iblk]+ibs;
			ipos = nd2ind[iind];
			for (j=ia[i];j<ia[i+1];j++) {
				k = iptrschur[ipos];
				jaschur[k] = ja[j];
				ja2schur[k] = ja2blk[j];
				iptrschur[ipos]++;
			};
		};
	};

	for (i=0;i<nztot;i++) {
		irow = isendarr[i*4];
		iblk = isendarr[i*4+1];
		jj = isendarr[i*4+2];
		jblk = isendarr[i*4+3];
		ibs = ibsblk[iblk];
		iind = irow-_blks[iblk]+ibs;
		ipos = nd2ind[iind];
		k = iptrschur[ipos];
		jaschur[k] = jj;
		ja2schur[k] = jblk;
		iptrschur[ipos]++;
	};

	nzmax = 0;

	for (i=0;i<nschur;i++) {
		nzloc = iaschur[i+1]-iaschur[i];
		if (nzloc > nzmax) nzmax = nzloc;
	};

	iiarr = new CIntInt [nzmax];
	if (!iiarr) MemoryFail (funcname);
	listloc = new int [nzmax];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nzmax];
	if (!list2loc) MemoryFail (funcname);

	int jloc;

	for (i=0;i<nschur;i++) {
		nzloc = 0;
		for (j=iaschur[i];j<iaschur[i+1];j++) {
			iiarr[nzloc].intvalue = jaschur[j];
			iiarr[nzloc].int2value = nzloc;
		};
		//qsort (iiarr, nzloc, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + nzloc);
		for (j=iaschur[i];j<iaschur[i+1];j++) {
			jloc = j-iaschur[i];
			ind = iiarr[jloc].int2value;
			jind = ind+iaschur[i];
			listloc[jloc] = jaschur[jind];
			list2loc[jloc] = ja2schur[jind];
		};
		for (j=0;j<nzloc;j++) {
			jind = j+iaschur[i];
			jaschur[jind] = listloc[j];
			ja2schur[jind] = list2loc[j];
		};
	};

	int *iaschur1;

	iaschur1 = new int [nschur+1];
	if (!iaschur1) MemoryFail (funcname);

	iaschur1[0] = 0;
	nz = 0;

	for (i=0;i<nschur;i++) {
		j = iaschur[i];
		jaschur[nz] = jaschur[j];
		ja2schur[nz] = ja2schur[j];
		nz++;
		for (j=iaschur[i]+1;j<iaschur[i+1];j++) {
			if (jaschur[j] != jaschur[nz-1]) {
				jaschur[nz] = jaschur[j];
				ja2schur[nz] = ja2schur[j];
				nz++;
			};
		};
		iaschur1[i+1] = nz;
	};

	for (i=0;i<=nschur;i++) iaschur[i] = iaschur1[i];

// Store schur data in the matrix format

	nzschur = iaschur[nschur];

	CSMatrix aschur;

	delete [] aschur.list;
	delete [] aschur.list2;
	delete [] aschur.ia;
	delete [] aschur.ja;
	delete [] aschur.ja2;

	aschur.list = listschur;
	aschur.list2 = list2schur;
	aschur.ia = iaschur;
	aschur.ja = jaschur;
	aschur.ja2 = ja2schur;

	aschur.m = m;
	aschur.n = n;
	aschur.nsupc = nsupc;
	aschur.nsupr = nsupr;
	aschur.nlist = nschur;
	aschur.nlist2 = nschur;
	aschur.nzja = nzschur;
	aschur.nzja2 = nzschur;

// Free work arrays

	delete [] ibsblklist;
	delete [] list2blk;
	delete [] ja2blk;
	delete [] iacpu;
	delete [] iptr;
	delete [] mtrarr;
	delete [] listcpu;

	return aschur;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute 3D partitioning of the surface boundary according to the partitioning type
//========================================================================================
void CSMatrix::Complete3DReduction (int _decomp_type, // Compute 3D partitioning of the surface boundary according to the partitioning type
													CSMatrix &_amtrschur, CSMatrix &_amtrblk,
													int _nblks, int *_blktype, int *_blk2cpu,
													int *_order, 
													int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew) {

	if (_decomp_type == 0) {
		Complete3DReduction_3Dtype0 (_amtrschur, _amtrblk,
											_nblks, _blktype, _blk2cpu,
											_order, 
											_nblksnew, _blksnew, _blktypenew, _blk2cpunew);
	} else if (_decomp_type == 1 || _decomp_type == 2) {
		Complete3DReduction_3Dtype1_2 (_decomp_type,
													_amtrschur, _amtrblk,
													_nblks, _blktype, _blk2cpu,
													_order, 
													_nblksnew, _blksnew, _blktypenew, _blk2cpunew);
/*
	} else if (_decomp_type == 2) {
		Complete3DReduction_3Dtype2 (_amtrschur, _amtrblk,
											_nblks, _blktype, _blk2cpu,
											_order, 
											_nblksnew, _blksnew, _blktypenew, _blk2cpunew);
	} else if (_decomp_type == 3) {
		Complete3DReduction_3Dtype3 (_amtrschur, _amtrblk,
											_nblks, _blktype, _blk2cpu,
											_order, 
											_nblksnew, _blksnew, _blktypenew, _blk2cpunew);
	} else if (_decomp_type == 4) {
		Complete3DReduction_3Dtype4 (_amtrschur, _amtrblk,
											_nblks, _blktype, _blk2cpu,
											_order, 
											_nblksnew, _blksnew, _blktypenew, _blk2cpunew);
*/
	};

};

// Author: Kharchenko S.A.
// CSMatrix: 3D partitioning of the surface boundary according to the partitioning type==0
//========================================================================================
void CSMatrix::Complete3DReduction_3Dtype0 ( // 3D partitioning of the surface boundary according to the partitioning type==0
															CSMatrix &_amtrschur, CSMatrix &_amtrblk,
															int _nblks, int *_blktype, int *_blk2cpu,
															int *_order, 
															int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew) {

	const char *funcname = "Complete3DReduction_3Dtype0";

// Fill the structures

	int ntot = _amtrschur.n;

	_nblksnew = 1;

	_blksnew = new int [_nblksnew+1];
	if (!_blksnew) MemoryFail (funcname);

	_blksnew[0] = 0;
	_blksnew[1] = ntot;

	_blktypenew = new int [_nblksnew];
	if (!_blktypenew) MemoryFail (funcname);

	_blktypenew[0] = 2;

	_blk2cpunew = new int [_nblksnew];
	if (!_blk2cpunew) MemoryFail (funcname);

	_blk2cpunew[0] = 0;

};

// Author: Kharchenko S.A.
// CSMatrix: 3D partitioning of the surface boundary according to the partitioning type==1 || type==2
//========================================================================================
void CSMatrix::Complete3DReduction_3Dtype1_2 (int _decomp_type, // 3D partitioning of the surface boundary according to the partitioning type==1 || type==2
																CSMatrix &_amtrschur, CSMatrix &_amtrblk,
																int _nblks, int *_blktype, int *_blk2cpu,
																int *_order, 
																int &_nblksnew, int *&_blksnew, int *&_blktypenew, int *&_blk2cpunew) {

	const char *funcname = "Complete3DReduction_3Dtype1_2";

// Find all surface pairs and their parititioning

	int npairs;
	int *iapairs, *listpairs, *cpupairs;
	int *orderloc, *decomp_pairs;

	Complete3DReduction_FindAndPartitionSurfacePairs (_amtrschur, _amtrblk,
																		_nblks, _blktype, _blk2cpu,
																		npairs, iapairs, listpairs, cpupairs,
																		orderloc, decomp_pairs);

// Sort pair blocks according to the cpu numbers

	CIntInt *iiarr;

	iiarr = new CIntInt [2*npairs];
	if (!iiarr) MemoryFail (funcname);

	int i;

	for (i=0;i<npairs;i++) {
		iiarr[i*2].intvalue = cpupairs[i*2];
		iiarr[i*2].int2value = i*2;
		iiarr[i*2+1].intvalue = cpupairs[i*2+1];
		iiarr[i*2+1].int2value = i*2+1;
	};

	//qsort (iiarr, 2*npairs, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + 2*npairs);

	int *orderblk;

	orderblk = new int [2*npairs];
	if (!orderblk) MemoryFail (funcname);

	for (i=0;i<2*npairs;i++) orderblk[iiarr[i].int2value] = 0;

// Fill return structures

	int ntot = _amtrschur.n;

	_nblksnew = npairs*3+1;

	_blksnew = new int [_nblksnew+1];
	if (!_blksnew) MemoryFail (funcname);
	_blktypenew = new int [_nblksnew];
	if (!_blktypenew) MemoryFail (funcname);
	_blk2cpunew = new int [_nblksnew];
	if (!_blk2cpunew) MemoryFail (funcname);

	_blksnew[0] = 0;

	int iold, ipair, inum, ni;

	for (i=0;i<2*npairs;i++) {
		iold = iiarr[i].int2value;
		ipair = iold / 2;
		inum = iold - ipair*2;
		if (inum == 0) {
			ni = decomp_pairs[ipair*4+1]-decomp_pairs[ipair*4];
		} else {
			ni = decomp_pairs[ipair*4+2]-decomp_pairs[ipair*4+1];
		};
		_blksnew[i+1] = _blksnew[i] + ni;
		_blk2cpunew[i] = cpupairs[iold];
		_blktypenew[i] = 2;
	};

	for (i=0;i<npairs;i++) {
		ni = decomp_pairs[ipair*4+3]-decomp_pairs[ipair*4+2];
		_blksnew[2*npairs+i+1] = _blksnew[2*npairs+i] + ni;
		_blk2cpunew[2*npairs+i] = cpupairs[i*2];
		_blktypenew[2*npairs+i] = 1;
	};

	_blksnew[3*npairs+2] = ntot;
	_blk2cpunew[3*npairs+1] = 0;
	_blktypenew[3*npairs+1] = 1;

// Compute order

	int *imask;
	int *itypearr;
	int *ishiftarr;

	imask = new int [ntot];
	if (!imask) MemoryFail (funcname);
	itypearr = new int [ntot];
	if (!itypearr) MemoryFail (funcname);
	ishiftarr = new int [ntot];
	if (!ishiftarr) MemoryFail (funcname);

	for (i=0;i<ntot;i++) imask[i] = -1;

	int j, ind, jloc;

	for (ipair=0;ipair<npairs;ipair++) {
		for (j=iapairs[ipair];j<iapairs[ipair+1];j++) {
			ind = listpairs[j];
			imask[ind] = ipair;
			jloc = orderloc[j];
			if (jloc < decomp_pairs[4*ipair+1]) {
				itypearr[ind] = 0;
				ishiftarr[ind] = jloc;
			} else if (jloc >= decomp_pairs[4*ipair+1] && jloc < decomp_pairs[4*ipair+2]) {
				itypearr[ind] = 1;
				ishiftarr[ind] = jloc-decomp_pairs[4*ipair+1];
			} else {
				itypearr[ind] = 2;
				ishiftarr[ind] = jloc-decomp_pairs[4*ipair+2];
			};
		};
	};

	int *iptr;

	iptr = new int [_nblksnew];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<_nblksnew;i++) iptr[i] = _blksnew[i];

	int itype, ishift, iblk, iblknew;

	for (i=0;i<ntot;i++) {
		ipair = imask[ind];
		itype = itypearr[ind];
		ishift = ishiftarr[ind];
		if (ipair != -1) {
			if (itype == 2) {
				_order[i] = _blksnew[2*npairs+ipair]+ishift;
			} else {
				iblk = ipair*2+itype;
				iblknew = orderblk[iblk];
				_order[i] = _blksnew[iblknew]+ishift;
			};
		};
	};

	int k;

	for (i=0;i<ntot;i++) {
		ipair = imask[ind];
		itype = itypearr[ind];
		if (ipair == -1) {
			k = iptr[_nblksnew];
			_order[i] = k;
			iptr[_nblksnew]++;
		};
	};

// Reassign if necessary

	if (_decomp_type == 1) {
		_nblksnew = 2*npairs+1;
		_blksnew[2*npairs+2] = ntot;
		_blk2cpunew[2*npairs+1] = 0;
		_blktypenew[2*npairs+1] = 1;
	};

// Free work arrays

	delete [] iapairs;
	delete [] listpairs;
	delete [] cpupairs;
	delete [] orderloc;
	delete [] decomp_pairs;
	delete [] iiarr;
	delete [] orderblk;
	delete [] imask;
	delete [] itypearr;
	delete [] ishiftarr;
	delete [] iptr;

};

// Author: Kharchenko S.A.
// CSMatrix: 3D partitioning of the surface boundary according to the partitioning type==1 || type==2
//========================================================================================
void CSMatrix::Complete3DReduction_FindAndPartitionSurfacePairs ( // Compute the set of surface pairs and partition them
												CSMatrix &_amtrschur, CSMatrix &_amtrblk,
												int _nblks, int *_blktype, int *_blk2cpu,
												int &_npairs, int *&_iapairs, int *&_listpairs, int *&_cpupairs,
												int *&_order, int *&_decomp_pairs) {

	const char *funcname = "Complete3DReduction_FindAndPartitionSurfacePairs";

// Compute the maximal number of cpus

	int nproc = 0;

	int iblk, icpu;

	for (iblk=0;iblk<_nblks;iblk++) {
		icpu = _blk2cpu[iblk];
		if (icpu > nproc) nproc = icpu;
	};

	nproc++;

// Mark all nodes that have only volume neibours

	int ntot = _amtrschur.n;

	int *imask;

	imask = new int [ntot];
	if (!imask) MemoryFail (funcname);

	int i, j, jblk;

	for (i=0;i<ntot;i++) imask[i] = -1;

	for (i=0;i<ntot;i++) {
		int icheck = 1;
		for (j=_amtrblk.ia[i];j<_amtrblk.ia[i+1];j++) {
			jblk = _amtrblk.ja[j];
			if (_blktype[jblk] != 3) {
				icheck = -1;
			};
		};
		if (_amtrblk.ia[i+1]-_amtrblk.ia[i] != 2) icheck = -1;
		if (icheck == 1) imask[i] = 1;
	};

	int *listnd;

	listnd = new int [ntot];
	if (!listnd) MemoryFail (funcname);

	int nlistnd = 0;

	for (i=0;i<ntot;i++) {
		if (imask[i] == 1) {
			listnd[nlistnd] = i;
			nlistnd++;
		};
	};

// Create per cpu lists

	int *iandcpu;
	int *iptrnd;

	iandcpu = new int [nproc+1];
	if (!iandcpu) MemoryFail (funcname);
	iptrnd = new int [nproc];
	if (!iptrnd) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iandcpu[i] = 0;

	int irow, ibs, iblk1, iblk2, icpu1, icpu2;

	for (i=0;i<nlistnd;i++) {
		irow = listnd[i];
		ibs = _amtrblk.ia[irow];
		iblk1 = _amtrblk.ja[ibs];
		iblk2 = _amtrblk.ja[ibs+1];
		icpu1 = _blk2cpu[iblk1];
		icpu2 = _blk2cpu[iblk2];
		if (icpu1 == icpu2) throw " CSMatrix::Complete3DReduction_FindAndPartitionSurfacePairs: threre coincident cpus in block sparsity ";
		iandcpu[icpu1+1]++;
		iandcpu[icpu2+1]++;
	};

	for (i=0;i<nproc;i++) iandcpu[i+1] = iandcpu[i]+iandcpu[i+1];
	for (i=0;i<nproc;i++) iptrnd[i] = iandcpu[i];

	int nzjand = iandcpu[nproc];

	int *jandcpu;

	jandcpu = new int [nzjand];
	if (!jandcpu) MemoryFail (funcname);

	int k;

	for (i=0;i<nlistnd;i++) {
		irow = listnd[i];
		ibs = _amtrblk.ia[irow];
		iblk1 = _amtrblk.ja[ibs];
		iblk2 = _amtrblk.ja[ibs+1];
		icpu1 = _blk2cpu[iblk1];
		icpu2 = _blk2cpu[iblk2];
		k = iptrnd[icpu1];
		jandcpu[k] = irow;
		iptrnd[icpu1]++;
		k = iptrnd[icpu2];
		jandcpu[k] = irow;
		iptrnd[icpu2]++;
	};

// Sort data in the lists of each cpu

	int ni;
	int nzmax = 0;

	for (i=0;i<nproc;i++) {
		ni = iandcpu[i+1]-iandcpu[i];
		if (ni > nzmax) nzmax = ni;
	};

	CIntInt *iiarr;
	int *listloc;

	iiarr = new CIntInt [nzmax];
	if (!iiarr) MemoryFail (funcname);
	listloc = new int [nzmax];
	if (!listloc) MemoryFail (funcname);

	int jloc, jlocold, jold;

	for (icpu=0;icpu<nproc;icpu++) {
		ni = iandcpu[i+1]-iandcpu[i];
		for (j=iandcpu[icpu];j<iandcpu[icpu+1];j++) {
			jloc = j-iandcpu[icpu];
			irow = jandcpu[j];
			ibs = _amtrblk.ia[irow];
			iblk1 = _amtrblk.ja[ibs];
			iblk2 = _amtrblk.ja[ibs+1];
			icpu1 = _blk2cpu[iblk1];
			icpu2 = _blk2cpu[iblk2];
			if (icpu1 != icpu) {
				iiarr[jloc].intvalue = icpu1;
			} else {
				iiarr[jloc].intvalue = icpu2;
			};
			iiarr[jloc].int2value = jloc;
		};
		//qsort (iiarr, ni, sizeof(CIntInt), CIntInt::CompareIntInt);
		sort (iiarr, iiarr + ni);
		for (j=iandcpu[icpu];j<iandcpu[icpu+1];j++) {
			jloc = j-iandcpu[icpu];
			jlocold = iiarr[jloc].int2value;
			jold = jlocold+iandcpu[icpu];
			irow = jandcpu[jold];
			listloc[jloc] = irow;
		};
		for (j=iandcpu[icpu];j<iandcpu[icpu+1];j++) {
			jloc = j-iandcpu[icpu];
			jandcpu[j] = listloc[jloc];
		};
	};

// Create block description of the sublists

	int *iablkndcpu;

	iablkndcpu = new int [nproc+1];
	if (!iablkndcpu) MemoryFail (funcname);

	iablkndcpu[0] = 0;
	int nzblk = 0;

	int jcpu, jcpuprev;

	for (icpu=0;icpu<nproc;icpu++) {
		for (j=iandcpu[icpu];j<iandcpu[icpu+1];j++) {
			irow = jandcpu[j];
			ibs = _amtrblk.ia[irow];
			iblk1 = _amtrblk.ja[ibs];
			iblk2 = _amtrblk.ja[ibs+1];
			icpu1 = _blk2cpu[iblk1];
			icpu2 = _blk2cpu[iblk2];
			jcpu = icpu1;
			if (jcpu == icpu) jcpu = icpu2;
			if (j==iandcpu[icpu]) {
				jcpuprev = jcpu;
				nzblk++;
			} else {
				if (jcpu != jcpuprev) {
					jcpuprev = jcpu;
					nzblk++;
				};
			};
		};
		iablkndcpu[icpu+1] = nzblk;
	};

	int *jablkndcpu;
	int *ibsblkndcpu;

	jablkndcpu = new int [nzblk];
	if (!jablkndcpu) MemoryFail (funcname);
	ibsblkndcpu = new int [nzblk+1];
	if (!ibsblkndcpu) MemoryFail (funcname);

	nzblk = 0;

	for (icpu=0;icpu<nproc;icpu++) {
		for (j=iandcpu[icpu];j<iandcpu[icpu+1];j++) {
			irow = jandcpu[j];
			ibs = _amtrblk.ia[irow];
			iblk1 = _amtrblk.ja[ibs];
			iblk2 = _amtrblk.ja[ibs+1];
			icpu1 = _blk2cpu[iblk1];
			icpu2 = _blk2cpu[iblk2];
			jcpu = icpu1;
			if (jcpu == icpu) jcpu = icpu2;
			if (j==iandcpu[icpu]) {
				ibsblkndcpu[nzblk] = j;
				jcpuprev = jcpu;
				jablkndcpu[nzblk] = jcpu;
				nzblk++;
			} else {
				if (jcpu != jcpuprev) {
					ibsblkndcpu[nzblk] = j;
					jcpuprev = jcpu;
					nzblk++;
				};
			};
		};
	};
	ibsblkndcpu[nzblk] = nzjand;

// Split all nodes into pairs

	int nlistpairs = 0;
	_npairs = 0;

	for (icpu=0;icpu<nproc;icpu++) {
		for (j=iablkndcpu[icpu];j<iablkndcpu[icpu+1];j++) {
			jcpu = jablkndcpu[j];
			if (jcpu > icpu) {
				_npairs++;
				nlistpairs += ibsblkndcpu[j+1]-ibsblkndcpu[j];
			};
		};
	};

	_iapairs = new int [_npairs+1];
	if (!_iapairs) MemoryFail (funcname);
	_listpairs = new int [nlistpairs];
	if (!_listpairs) MemoryFail (funcname);
	_cpupairs = new int [_npairs*2];
	if (!_cpupairs) MemoryFail (funcname);

	_npairs = 0;

	_iapairs[0] = 0;
	nlistpairs = 0;

	for (icpu=0;icpu<nproc;icpu++) {
		for (j=iablkndcpu[icpu];j<iablkndcpu[icpu+1];j++) {
			jcpu = jablkndcpu[j];
			if (jcpu > icpu) {
				_cpupairs[_npairs*2] = icpu;
				_cpupairs[_npairs*2+1] = jcpu;
				for (k=ibsblkndcpu[j];k<ibsblkndcpu[j+1];k++) {
					_listpairs[nlistpairs] = jandcpu[k];
					nlistpairs++;
				};
				_npairs++;
				_iapairs[_npairs] = nlistpairs;
			};
		};
	};

	for (i=0;i<_npairs;i++) {
		ni = _iapairs[i+1]-_iapairs[i];
		ibs = _iapairs[i];
		sort (_listpairs+ibs, _listpairs+ibs + ni);
	};

// Free current work arrays

	delete [] imask;
	delete [] listnd;
	delete [] iandcpu;
	delete [] iptrnd;
	delete [] jandcpu;
	delete [] iiarr;
	delete [] listloc;
	delete [] iablkndcpu;
	delete [] jablkndcpu;
	delete [] ibsblkndcpu;

// Partition pairs data

	_order = new int [nlistpairs];
	if (!_order) MemoryFail (funcname);
	_decomp_pairs = new int [_npairs*4];
	if (!_decomp_pairs) MemoryFail (funcname);

	int *irow2ind;

	imask = new int [ntot];
	if (!imask) MemoryFail (funcname);
	irow2ind = new int [ntot];
	if (!irow2ind) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<ntot;i++) imask[i] = icycle;

	CSMatrix asub;
	CSMatrix aord;

	int ipair, ind1, ind2;

	for (ipair=0;ipair<_npairs;ipair++) {

		icycle++;

		ni = _iapairs[ipair+1]-_iapairs[ipair];
		ibs = _iapairs[ipair];

		_amtrschur.SubmatrixList (ni, _listpairs+ibs,
											icycle, imask, irow2ind,
											asub);

		asub.OrderPrfMtr (-1, _order+ibs);

		aord = asub.OrdMtr (_order+ibs);

		aord.FindSeparator (ind1, ind2);

		if (ind1 < 0 || ind2 < 0) {
			ind1 = 0;
			ind2 = 0;
		};

		_decomp_pairs[ipair*4] = 0;
		_decomp_pairs[ipair*4+1] = ind1;
		_decomp_pairs[ipair*4+2] = ind2;
		_decomp_pairs[ipair*4+3] = ni;

	};

// Free work arrays

	delete [] imask;
	delete [] irow2ind;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition Schur complement submatrix by ND METIS with boundary local assignments
//========================================================================================
void CSMatrix::PartSchurSubmatrix (int _decomptype, int _icolor, // Partition Schur complement submatrix by ND METIS with boundary local assignments
												int *_blk2color, int *_imaskblk, CSMatrix &_ablkref,
												int *_order,
												int &_nblksnew, int *&_blksnew, int *&_blk2colornew) {

	const char *funcname = "PartSchurSubmatrix";

	if (nlist != _ablkref.GetNlist ()) {
		throw " CSMatrix::PartSchurSubmatrix: incorrect block submatrix ";
	};

// Scan the reference block matrix and find all valid pairs

	CInd2Int *ind2arr;

	ind2arr = new CInd2Int [nlist];
	if (!ind2arr) MemoryFail (funcname);

	int *piablkref = _ablkref.GetIa ();
	int *pjablkref = _ablkref.GetJa ();

	int nz = 0;

	int i, nzj, j, iblk1, iblk2;

	for (i=0;i<nlist;i++) {
		nzj = piablkref[i+1]-piablkref[i];
		if (nzj == 2) {
			j = piablkref[i];
			iblk1 = pjablkref[j];
			iblk2 = pjablkref[j+1];
			if (_imaskblk[iblk1] == -1 && _imaskblk[iblk2] == -1) {
				ind2arr[nz].indx = iblk1;
				ind2arr[nz].indy = iblk2;
				ind2arr[nz].intvalue = i;
				nz++;
			};
		};
	};

// Perform pairs sort

	//qsort (ind2arr, nz, sizeof(CInd2Int), CInd2Int::CompareInd2Int);
	sort(ind2arr, ind2arr + nz);

// Detect all colors pairs subsets

	int *iapairs;

	iapairs = new int [nz+1];
	if (!iapairs) MemoryFail (funcname);

	int npairs = 0;

	iapairs[0] = 0;

	int ip = -1;

	if (nz > 0) {
		ip = 0;
		npairs++;
		iapairs[npairs] = ip+1;
	};

	while (ip < nz-1) {
		if (ind2arr[ip+1].indx == ind2arr[ip].indx && ind2arr[ip+1].indy == ind2arr[ip].indy) {
			ip++;
			iapairs[npairs] = ip+1;
		} else {
			ip++;
			npairs++;
			iapairs[npairs] = ip+1;
		};
	};

// Filter too small set sizes

	int isizemax = 0;

	int isize;

	for (i=0;i<npairs;i++) {
		isize = iapairs[i+1]-iapairs[i];
		if (isize > isizemax) isizemax = isize;
	};

	isizemax = isizemax / 30;

	int *ialoc;

	ialoc = new int [npairs+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	int npairsnew = 0;

	nz = 0;

	for (i=0;i<npairs;i++) {
		isize = iapairs[i+1]-iapairs[i];
		if (isize >= isizemax) {
			for (j=iapairs[i];j<iapairs[i+1];j++) {
				ind2arr[nz] = ind2arr[j];
				nz++;
			};
			ialoc[npairsnew+1] = nz;
			npairsnew++;
		};
	};

	npairs = npairsnew;

	delete [] iapairs;

	iapairs = ialoc;

// Create initial direct and inverse orderings

	int *order;
	int *iorder;

	order = new int [nlist];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlist];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlist;i++) order[i] = -1;

	nz = iapairs[npairsnew];

	for (i=0;i<nz;i++) {
		j = ind2arr[i].intvalue;
		order[j] = i;
		iorder[i] = j;
	};

	for (i=0;i<nlist;i++) {
		if (order[i] == -1) {
			order[i] = nz;
			iorder[nz] = i;
			nz++;
		};
	};

// Reorder the matrix

	CSMatrix aini = *this;

	aini.m = nlist;
	aini.n = nlist;
	aini.nsupc = nlist;
	aini.nsupr = nlist;

	CSMatrix ao;

	ao = aini.OrdMtr (order);

// For each submatrix perform its partitioning

	int *orderpairs;
	int *iorderpairs;
	int *separr;

	orderpairs = new int [nlist];
	if (!orderpairs) MemoryFail (funcname);
	iorderpairs = new int [nlist];
	if (!iorderpairs) MemoryFail (funcname);
	separr = new int [npairs*2];
	if (!separr) MemoryFail (funcname);

	CSMatrix asub, asubo;

	int ni, ibeg, iend, indsep1, indsep2;

	int nz0 = 0;

	for (i=0;i<npairs;i++) {

		ni = iapairs[i+1]-iapairs[i];
		ibeg = iapairs[i];
		iend = iapairs[i+1]-1;

		ao.Submatrix (ibeg, iend, asub);

		asub.OrderPrfMtr (-1, orderpairs+ibeg);

		asubo = asub.OrdMtr (orderpairs+ibeg);

		asubo.FindSeparator (indsep1, indsep2);

		if (indsep1 < 0 || indsep2 < 0) {
			indsep1 = 0;
			indsep2 = 0;
		};

		separr[i*2] = indsep1;
		separr[i*2+1] = indsep2;

		nz0 += indsep2;

		for (j=0;j<ni;j++) orderpairs[ibeg+j] = orderpairs[ibeg+j]+ibeg;

	};

	nz = iapairs[npairs];

	for (i=0;i<nz;i++) iorderpairs[orderpairs[i]] = i;

// Sort each part according to the color number

	CIntInt *iiarr;

	iiarr = new CIntInt [npairs*2];
	if (!iiarr) MemoryFail (funcname);

	for (i=0;i<npairs;i++) {
//		iiarr[i*2].intvalue = ind2arr[i].indx;
		iiarr[i*2].intvalue = _blk2color[ind2arr[i].indx];
		iiarr[i*2].int2value = i*2;
//		iiarr[i*2+1].intvalue = ind2arr[i].indy;
		iiarr[i*2+1].intvalue = _blk2color[ind2arr[i].indy];
		iiarr[i*2+1].int2value = i*2+1;
	};

	//qsort (iiarr, npairs*2, sizeof(CIntInt), CIntInt::CompareIntInt);
	sort (iiarr, iiarr + npairs*2);

	int *ordblk;

	ordblk = new int [npairs*2];
	if (!ordblk) MemoryFail (funcname);

	for (i=0;i<npairs*2;i++) ordblk[iiarr[i].int2value] = i;

// Initial partitioning

	int nblkspairs = 2*npairs;
	if (_decomptype > 0) nblkspairs = 3*npairs;

	int *blksini;
	int *blk2colorini;

	blksini = new int [nblkspairs+1];
	if (!blksini) MemoryFail (funcname);
	blk2colorini = new int [nblkspairs];
	if (!blk2colorini) MemoryFail (funcname);

	for (i=0;i<=nblkspairs;i++) blksini[i] = 0;

	for (i=0;i<npairs;i++) {
		blksini[i*2+1] = separr[i*2];
		blksini[i*2+2] = separr[i*2+1]-separr[i*2];
		blk2colorini[i*2] = _blk2color[ind2arr[i].indx];
		blk2colorini[i*2+1] = _blk2color[ind2arr[i].indy];
	};
	if (_decomptype > 0) {
		for (i=0;i<npairs;i++) {
			blksini[2*npairs+i+1] = (iapairs[i+1]-iapairs[i])-separr[i*2+1];
			blk2colorini[2*npairs+i] = _blk2color[ind2arr[i].indx];
		};
	};

	for (i=0;i<nblkspairs;i++) blksini[i+1] = blksini[i]+blksini[i+1];

// Create new matrix ordering and partitioning

	_nblksnew = nblkspairs+1;

	_blksnew = new int [_nblksnew+1];
	if (!_blksnew) MemoryFail (funcname);
	_blk2colornew = new int [_nblksnew];
	if (!_blk2colornew) MemoryFail (funcname);

	_blksnew[0] = 0;

	for (i=0;i<2*npairs;i++) {
		j = iiarr[i].int2value;
		_blksnew[i+1] = blksini[j+1]-blksini[j];
		_blk2colornew[i] = blk2colorini[j];
	};

	if (_decomptype > 0) {
		for (i=0;i<npairs;i++) {
			_blksnew[2*npairs+i+1] = blksini[2*npairs+i+1]-blksini[2*npairs+i];
			_blk2colornew[2*npairs+i] = blk2colorini[2*npairs+i];
		};
	};

	nz = iapairs[npairs];

	_blksnew[_nblksnew] = nlist-nz0;
	if (_decomptype > 0) {
		_blksnew[_nblksnew] = nlist-nz;
	};
	_blk2colornew[_nblksnew-1] = _icolor;

	for (i=0;i<_nblksnew;i++) _blksnew[i+1] = _blksnew[i]+_blksnew[i+1];

	for (i=0;i<nz;i++) {
		j = iorderpairs[i];
		_order[i] = iorder[j];
	};
	for (i=nz;i<nlist;i++) {
		_order[i] = iorder[i];
	};

	for (i=0;i<nlist;i++) {
		iorder[i] = _order[i];
		order[_order[i]] = i;
	};

	int ibs = _blksnew[_nblksnew-1];

	int ip0, ip1, ip2, ip3, iblk1o, iblk2o;

	if (_decomptype == 0) {
		for (i=0;i<npairs;i++) {
			ip0 = iapairs[i];
			ip1 = ip0+(blksini[i*2+1]-blksini[i*2]);
			ip2 = ip1+(blksini[i*2+2]-blksini[i*2+1]);
			ip3 = iapairs[i+1];
			iblk1 = i*2;
			iblk2 = i*2+1;
			iblk1o = ordblk[iblk1];
			iblk2o = ordblk[iblk2];
			for (j=ip0;j<ip1;j++) {
				orderpairs[j] = _blksnew[iblk1o]+j-ip0;
			};
			for (j=ip1;j<ip2;j++) {
				orderpairs[j] = _blksnew[iblk2o]+j-ip1;
			};
			for (j=ip2;j<ip3;j++) {
				orderpairs[j] = ibs;
				ibs++;
			};
		};
	} else {
		for (i=0;i<npairs;i++) {
			ip0 = iapairs[i];
			ip1 = ip0+(blksini[i*2+1]-blksini[i*2]);
			ip2 = ip1+(blksini[i*2+2]-blksini[i*2+1]);
			ip3 = iapairs[i+1];
			iblk1 = i*2;
			iblk2 = i*2+1;
			iblk1o = ordblk[iblk1];
			iblk2o = ordblk[iblk2];
			for (j=ip0;j<ip1;j++) {
				orderpairs[j] = _blksnew[iblk1o]+j-ip0;
			};
			for (j=ip1;j<ip2;j++) {
				orderpairs[j] = _blksnew[iblk2o]+j-ip1;
			};
			for (j=ip2;j<ip3;j++) {
				orderpairs[j] = _blksnew[2*npairs+i]+j-ip2;
			};
		};
	};

	for (i=nz;i<nlist;i++) {
		orderpairs[i] = ibs;
		ibs++;
	};

	for (i=0;i<nlist;i++) iorderpairs[orderpairs[i]] = i;

	for (i=0;i<nlist;i++) {
		j = iorderpairs[i];
		_order[i] = iorder[j];
	};

	for (i=0;i<nlist;i++) iorder[_order[i]] = i;

	for (i=0;i<nlist;i++) _order[i] = iorder[i];

// Remove zero length data in partitioning

	if (_nblksnew > 1) {

		int nblksloc = 0;

		for (i=0;i<_nblksnew;i++) {
			if (_blksnew[i+1]-_blksnew[i] != 0) {
				_blksnew[nblksloc+1] = _blksnew[i+1];
				_blk2colornew[nblksloc] = _blk2colornew[i];
				nblksloc++;
			};
		};

		_nblksnew = nblksloc;

	};

// Free work arrays

	delete [] ind2arr;
	delete [] iapairs;
	delete [] order;
	delete [] iorder;
	delete [] orderpairs;
	delete [] iorderpairs;
	delete [] separr;
	delete [] iiarr;
	delete [] blksini;
	delete [] blk2colorini;
	delete [] ordblk;

};
