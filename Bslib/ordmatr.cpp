//------------------------------------------------------------------------------------------------
// File: ordmatr.cpp
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
#include <ctime>
#include <iomanip>

#include "globals.h"
#include "smatrix.h"
#include "gsmatrix.h"
#include "tree.h"
#include "slvparam.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSMatrix: Compute Nested Dissection or profile optimization ordering of the matrix
//========================================================================================
void CSMatrix::OrderPrfMtr (int _ordtyp, int *_order) const { // Compute Nested Dissection or profile optimization ordering of the matrix

	const char *funcname = "OrderPrfMtr";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Call METIS ordering routine

	int options[8], i;

	options[0] = 0;

	int nloc = nsupr;

	int *pn = &nloc;
	int *pzero = &i;

	i=0;

	int *iorder;
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);

	if (_ordtyp == -1) {

		int *ialoc;
		int *jaloc;

		ialoc = new int [nsupr+1];
		if (!ialoc) MemoryFail (funcname);
		jaloc = new int [as.ia[nsupr]];
		if (!jaloc) MemoryFail (funcname);

		int i, j, jj;

		ialoc[0] = 0;
		int nz = 0;

		for (i=0;i<nsupr;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj != i) {
					jaloc[nz] = jj;
					nz++;
				};
			};
			ialoc[i+1] = nz;
		};

		METIS_NodeND (pn, ialoc, jaloc, pzero, options, iorder, _order);

		delete [] ialoc;
		delete [] jaloc;

	} else if (_ordtyp == 0) {
		for (i=0;i<n;i++) {
			_order[i] = i;
		};
	} else {

// Reorder by GIBBS-POOLE-STOCKMEYER and GIBBS-KING algorithms

		if (false) {

// Allocate all necessary local arrays

			int wrklen = 6*nsupr+10;
			int *ian, *degree, *work;

			ian = new int [nsupr+1];
			if (!ian) MemoryFail (funcname);
			degree = new int [nsupr];
			if (!degree) MemoryFail (funcname);
			work = new int [wrklen];
			if (!work) MemoryFail (funcname);

// Recompute ja and compute new ia

			ian[0] = 0;
			int nz = 0;

			for (i=0;i<nsupr;i++) {
				for (int j=as.ia[i];j<as.ia[i+1];j++) {
					int jj = as.ja[j];
					if (jj != i) {
						as.ja[nz] = as.ja[j];
						nz++;
					};
				};
				ian[i+1] = nz;
			};

// Set order and degree arrays

			for (i=0;i<nsupr;i++) _order[i] = i+1;
			for (i=0;i<nsupr;i++) {
				degree[i] = ian[i+1]-ian[i];
			};

// Shift ia and ja arrays by 1 (Fortran programming style)

			for (i=0;i<ian[nsupr];i++) as.ja[i]++;
			for (i=0;i<nsupr;i++) ian[i]++;

// Optimize the profile and check

			int nloc = nsupr;
			int optpro = 1;

			int bandwd=0, profil=0, ergps=0, spaceloc=0;

			GPSKCA (nloc, degree, ian, as.ja, 
					optpro, wrklen, _order, work, 
					bandwd, profil, ergps, spaceloc);

//			if (ergps != 0) {
//				cout << " Error in GPSKCA" << endl;
//				throw " Error in GPSKCA";
//			};

//			if (spaceloc == -1) {
//				cout << " Error in GPSKCA add" << endl;
//				throw " Error in GPSKCA add";
//			};

// Shift order array back by 1

			for (i=0;i<nsupr;i++) _order[i]--;

// Delete wokr arrays

			delete [] ian;
			delete [] degree;
			delete [] work;

		} else {

/*
// Shift matrix data

			int nzas = as.ia[nsupr];

			for (i=0;i<nzas;i++) as.ja[i] += 1;
		
			int out = 66;
			int length = 2*nzas+10*nsupr;
			int *intbuf;

			intbuf = new int [length];
			if (!intbuf) MemoryFail (funcname);

			int nloc = n;

//			NEWOPT (&_ordtyp, &nloc, as.ia, &nzas, as.ja, _order, 
//					&out, &length, intbuf); 

			for (i=0;i<nsupr;i++) {
				_order[i] -= 1;
//				cout << " i = " << i << " order = " << _order[i] << endl;
			};

			delete [] intbuf;
*/
			int *xls;
			int *mask;

			int nzas = as.ia[nsupr];
			for (i=0;i<nzas;i++) as.ja[i] += 1;

			for (i=0;i<=nloc;i++) as.ia[i] += 1;

			xls = new int [nloc+1];
			if (!xls) MemoryFail (funcname);
			mask = new int [nloc+1];
			if (!mask) MemoryFail (funcname);

			genrcm (nloc, as.ia, as.ja, _order, xls, mask);

			for (i=0;i<nloc;i++) {
				_order[i] -= 1;
			};

			for (i=0;i<nloc;i++) {
				mask[_order[i]] = i;
			};

			for (i=0;i<nloc;i++) {
				_order[i] = mask[i];
			};

			delete [] xls;
			delete [] mask;

		};

	};

	delete [] iorder;

};

// Author: Kharchenko S.A.
// CSMatrix: Reorder the matrix according to the bordering parameter and compute profile optimization ordering for the remaining subblock
//========================================================================================
void CSMatrix::OrderBrdMtr (int _nbord, int *_order) const { // Reorder the matrix according to the bordering parameter and compute profile optimization ordering for the remaining subblock

	const char *funcname = "OrderBrdMtr";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Compute bordering nodes

// Allocate work arrays

	int *iorder;
	int *imask;
	int *lstloc;

	iorder = new int [nsupr+1];
	if (!iorder) MemoryFail (funcname);
	imask = new int [nsupr+1];
	if (!imask) MemoryFail (funcname);
	lstloc = new int [nsupr+1];
	if (!lstloc) MemoryFail (funcname);

// Init arrays

	int i;
	for (i=0;i<nsupr;i++) {
		imask[i] = as.ia[i+1]-as.ia[i];
		iorder[i] = 0;
	};
	iorder[nsupr] = 0;

// Count the number of nodes for each number of elements

	for (i=0;i<nsupr;i++) {
		int ii = imask[i];
		iorder[ii]++;
	};

// Compute the list of count numbers

	int icount = 0;

	int nlstloc = 0;

	for (i=nsupr;i>=0;i--) {
		if (iorder[i]>0) {
			icount += iorder[i];
			lstloc[nlstloc] = i;
			nlstloc++;
			if (icount >= _nbord) goto nextstage;
		};
	};

nextstage:;

// Create the list of bordering nodes

	for (i=0;i<nsupr;i++) _order[i] = -1;
	for (i=0;i<nsupr;i++) iorder[i] = -1;

	int nb = 0;

	for (int ilist=0;ilist<nlstloc;ilist++) {
		int nelem = lstloc[ilist];
		for (i=0;i<nsupr;i++) {
			if (imask[i] == nelem) {
				_order[i] = nsupr-nb-1;
				nb++;
			};
		};
	};

// Fill the ordering of the remaining nodes

	int nloc = 0;

	for (i=0;i<nsupr;i++) {
		if (_order[i] == -1) {
			_order[i] = nloc;
			nloc++;
		};
	};

	for (i=0;i<nsupr;i++) iorder[_order[i]] = i;

	for (i=0;i<nsupr;i++) {
		if (iorder[i] == -1) {
			cout << " Error in order array " << endl;
			throw "Error in order array";
		};
	};

// Reorder the matrix

	CSMatrix ao;

	ao = as.OrdMtr (_order);

// Compute the remaining submatrix

	int nzjas = ao.nzja;

	CSMatrix aloc (nsupr,nzjas);

// Compute local submatrix

	nloc = 0;
	int nzloc = 0;
	aloc.ia[0] = 0;
	for (i=0;i<nsupr-nb;i++) {
		for (int j=ao.ia[i];j<ao.ia[i+1];j++) {
			int jj = ao.ja[j];
			if (jj < nsupr-nb) {
				aloc.ja[nzloc] = jj;
				nzloc++;
			};
		};
		nloc++;
		aloc.ia[nloc] = nzloc;
	};

	for (i=0;i<nloc;i++) aloc.list[i] = i;

	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupc = nloc;
	aloc.nsupr = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

// Compute local ordering 

	int ordtyp = 1;

	aloc.OrderPrfMtr (ordtyp, lstloc); 

// Recompute global ordering

	for (i=0;i<nsupr;i++) iorder[_order[i]] = i;

	for (i=0;i<nsupr-nb;i++) {
		int iold = iorder[i];
		int inew = lstloc[i];
		_order[iold] = inew;
	};

// Delete working arrays

	delete [] iorder;
	delete [] lstloc;
	delete [] imask;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS into not equal parts with weights of each node/edge and compute its ordering
//========================================================================================
void CSMatrix::PartMtr (int _nparts, float *_weights, // Partition matrix by METIS into not equal parts with weights of each node/edge and compute its ordering
						int *_ndweights, int *_adjweights, int *_blks, int *_partition, int *_order) const {

	const char *funcname = "PartMtr_00";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local adjacency weights

	int *adjwloc;

	adjwloc = new int [as.nzja];
	if (!adjwloc) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<as.nzja;i++) adjwloc[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwloc[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwloc[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwloc[j] = _adjweights[k];
					};
				};
			};
			if (adjwloc[j] < 0) throw " PartMtr_00: index was not found or incorrect weights array";
		};
	};

// Call METIS 

	int options[5], volume, i1;

	options[0] = 0;

	int nloc = nsupr;

	int *pn = &nloc;
	int *pnp = &_nparts;
	int *pvol = &volume;
	int *pzero = &i;
	int *pone = &i1;

	i=0;
	i1=3;

	if (nloc > 0) {
		METIS_WPartGraphRecursive (pn, as.ia, as.ja, _ndweights, adjwloc, pone, pzero, pnp, _weights, options, pvol, _partition);
	};
//	METIS_PartGraphRecursive (pn, as.ia, as.ja, NULL, NULL, pzero, pzero, pnp, options, pvol, _partition);

// Search independent nodes for each processor

	int *pblks, *mask, *listlc, *lstmax, *iorder;

	pblks = new int [_nparts+2];
	if (!pblks) MemoryFail (funcname);
	mask = new int [nsupr];
	if (!mask) MemoryFail (funcname);
	listlc = new int [nsupr];
	if (!listlc) MemoryFail (funcname);
	lstmax = new int [nsupr];
	if (!lstmax) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nsupr;i++) mask[i] = -1;

	int icpu, jcpu;

	int nn=0;
	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			jcpu = _partition[jj];
			if (jcpu != icpu) {
				listlc[nn] = i;
				nn++;
				goto label1;
			};
		};
		mask[i] = icpu;
label1:;
	};

// Modify independent data (move small unrelated subblocks back to the bordering)

// Compute preliminary blocks partitioning and preliminary ordering

	for (i=0;i<_nparts+2;i++) _blks[i] = 0;

	for (i=0;i<nsupr;i++) {
		icpu = mask[i];
		if (icpu >= 0) {
			_blks[icpu+1]++;
		};
	};

	for (i=0;i<_nparts;i++) _blks[i+1] += _blks[i];

	_blks[_nparts+1] = nsupr;

// Set ordering array

	for (i=0;i<nsupr;i++) _order[i] = -1;
	for (i=0;i<_nparts;i++) pblks[i] = _blks[i];

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		if (mask[i] >= 0) {
			_order[i] = pblks[icpu];
			pblks[icpu]++;
		};
	};

	int nz = pblks[_nparts-1];
	for (i=0;i<nsupr;i++) {
		if (mask[i] == -1) {
			_order[i] = nz;
			nz++;
		};
	};

	for (i=0;i<nsupr;i++) {
		iorder[_order[i]] = i;
	};

// Compute reordered matrix

	CSMatrix ao;

	ao = as.OrdMtr (_order);

// Main cycle over the blocks

	int iblk;

	for (iblk=0;iblk<_nparts;iblk++) {

		int isupbeg = _blks[iblk];
		int isupend = _blks[iblk+1]-1;

		if (isupbeg <= isupend) {

			int *orderindep, *iorderindep;

			orderindep = new int [isupend-isupbeg+1];
			if (!orderindep) MemoryFail (funcname);
			iorderindep = new int [isupend-isupbeg+1];
			if (!iorderindep) MemoryFail (funcname);

// Compute submatrix for current block

			CSMatrix asub;

			ao.Submatrix (isupbeg,isupend,asub);

// Split each major block into unrelated subblocks

			int nblksindep, *blksindep;

			if (isupend-isupbeg+1 > 0) {
				asub.IndependentSubsets (nblksindep,blksindep,orderindep);
			} else {
				nblksindep = 0;
				blksindep = new int [nblksindep+1];
				if (!blksindep) MemoryFail (funcname);
				blksindep[0] = 0;
			};

// Compute inverse ordering

			for (i=0;i<isupend-isupbeg+1;i++) iorderindep[orderindep[i]] = i;

// Find block index with the largest related data

			int iblkmax = 0;
			int nlocblkmax = blksindep[1]-blksindep[0];

			int iblkindep;

			for (iblkindep=1;iblkindep<nblksindep;iblkindep++) {
				int ni = blksindep[iblkindep+1]-blksindep[iblkindep];
				if (ni > nlocblkmax) {
					iblkmax = iblkindep;
					nlocblkmax = ni;
				};
			};

// Move back small unrelated subblocks into the bordering

			for (iblkindep=0;iblkindep<nblksindep;iblkindep++) {
				if (iblkindep != iblkmax) {
					for (i=blksindep[iblkindep];i<blksindep[iblkindep+1];i++) {
						int iold = iorderindep[i]+isupbeg;
						int iini = iorder[iold];
						mask[iini] = -1;
						listlc[nn] = iini;
						nn++;
					};
				};
			};

// Free work arrays

			delete [] orderindep;
			delete [] iorderindep;
			delete [] blksindep;

		};

	};

// Free work memory

	CSMatrix mtrdummy;

	ao = mtrdummy;

// Compute again preliminary blocks partitioning and preliminary ordering

	for (i=0;i<_nparts+2;i++) _blks[i] = 0;

	for (i=0;i<nsupr;i++) {
		icpu = mask[i];
		if (icpu >= 0) {
			_blks[icpu+1]++;
		};
	};

	for (i=0;i<_nparts;i++) _blks[i+1] += _blks[i];

	_blks[_nparts+1] = nsupr;

// Set ordering array

	for (i=0;i<nsupr;i++) _order[i] = -1;
	for (i=0;i<_nparts;i++) pblks[i] = _blks[i];

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		if (mask[i] >= 0) {
			_order[i] = pblks[icpu];
			pblks[icpu]++;
		};
	};

	nz = pblks[_nparts-1];
	for (i=0;i<nsupr;i++) {
		if (mask[i] == -1) {
			_order[i] = nz;
			nz++;
		};
	};

	for (i=0;i<nsupr;i++) {
		iorder[_order[i]] = i;
	};

// Compute reordered matrix

	ao = as.OrdMtr (_order);

// Scan the last block

	int isupbeg = _blks[_nparts];
	int isupend = _blks[_nparts+1]-1;

	int *orderindep, *iorderindep;

	orderindep = new int [isupend-isupbeg+1];
	if (!orderindep) MemoryFail (funcname);
	iorderindep = new int [isupend-isupbeg+1];
	if (!iorderindep) MemoryFail (funcname);

// Compute submatrix for current block

	CSMatrix asub;

	ao.Submatrix (isupbeg,isupend,asub);

// Split each major block into unrelated subblocks

	int nblksindep, *blksindep;

	if (isupend-isupbeg+1 > 0) {
		asub.IndependentSubsets (nblksindep,blksindep,orderindep);
	} else {
		nblksindep = 0;
		blksindep = new int [nblksindep+1];
		if (!blksindep) MemoryFail (funcname);
		blksindep[0] = 0;
	};

// Compute inverse ordering

	for (i=0;i<isupend-isupbeg+1;i++) iorderindep[orderindep[i]] = i;

// Scan all independent blocks and move blocks related to just one processor

	int iblkindep;

	int icpusave, *imaskcpu, nprocloc;

	imaskcpu = new int [_nparts];
	if (!imaskcpu) MemoryFail (funcname);

	int icyclecpu = -1;

	for (i=0;i<_nparts;i++) imaskcpu[i] = icyclecpu;

	for (iblkindep=0;iblkindep<nblksindep;iblkindep++) {
		icyclecpu++;
		nprocloc = 0;
		icpusave = -1;
		for (i=blksindep[iblkindep];i<blksindep[iblkindep+1];i++) {
			int iold = iorderindep[i]+isupbeg;
			int iini = iorder[iold];
			if (mask[iini] >= 0) throw " PartMtr_00: wrong imask value";
			for (j=as.ia[iini];j<as.ia[iini+1];j++) {
				jj = as.ja[j];
				if (mask[jj] >= 0) {
					icpu = _partition[jj];
					if (imaskcpu[icpu] != icyclecpu) {
						imaskcpu[icpu] = icyclecpu;
						nprocloc++;
						icpusave = icpu;
					};
				};
			};
		};
		if (nprocloc == 1) {
			for (i=blksindep[iblkindep];i<blksindep[iblkindep+1];i++) {
				int iold = iorderindep[i]+isupbeg;
				int iini = iorder[iold];
				mask[iini] = icpusave;
				_partition[iini] = icpusave;
			};
		} else if (icpusave != -1) {
			for (i=blksindep[iblkindep];i<blksindep[iblkindep+1];i++) {
				int iold = iorderindep[i]+isupbeg;
				int iini = iorder[iold];
				mask[iini] = -1;
				_partition[iini] = icpusave;
			};
		};
	};

// Free work arrays

	delete [] orderindep;
	delete [] iorderindep;
	delete [] blksindep;
	delete [] imaskcpu;

// Free work memory

	ao = mtrdummy;

// Recompute listlc array

	nn=0;
	for (i=0;i<nsupr;i++) {
		if (mask[i] < 0) {
			listlc[nn] = i;
			nn++;
		};
	};

// Restore mask array

	int nnini = nn;

	for (i=0;i<nsupr;i++) {
		if (mask[i] >= 0) mask[i] = 0;
	};

// Go out of the cycle

	int ilist;

	for (ilist=0;ilist<nn;ilist++) {
		i = listlc[ilist];
		mask[i] = 1;
	};

// Search the list of remaining nodes

// Note: current implementation is O(nn^2) algorithm in the worst case

	int iwmax, nnodes, indmax, iwloc, nlstmax;

label2:;

	nnodes = 0;
	iwmax = -1;
	indmax = -1;
	nlstmax = 0;

// Compute the largest weight and count the number of nodes

	for (ilist=0;ilist<nn;ilist++) {
		i = listlc[ilist];
		if (mask[i] == -1) {
			nnodes++;
			iwloc = 0;
			icpu = _partition[i];
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				jcpu = _partition[jj];
				if (jcpu != icpu) {
					if (mask[jj] != 1) iwloc++;
				};
			};
			if (iwloc > iwmax) {
				indmax = i;
				iwmax = iwloc;
			};
		};
	};

// Return if the number of nodes is zero

	if (nnodes == 0) goto label3;

// Process all nodes if the weight is zero

	if (iwmax == 0) {
		for (ilist=0;ilist<nn;ilist++) {
			i = listlc[ilist];
			if (mask[i] == -1) mask[i] = 0;
		};
		goto label3;
	};

// Find the list of nodes with the largest weight

	for (ilist=0;ilist<nn;ilist++) {
		i = listlc[ilist];
		if (mask[i] == -1) {
			iwloc = 0;
			icpu = _partition[i];
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				jcpu = _partition[jj];
				if (jcpu != icpu) {
					if (mask[jj] != 1) iwloc++;
				};
			};
			if (iwloc == iwmax) {
				lstmax[nlstmax] = i;
				nlstmax++;
			};
		};
	};

// Move nodes with the largest weight to the bordering

	if (nlstmax >= 0) {
//		cout << " Indmax = " << indmax << " iwmax = " << iwmax << endl;
		for (ilist=0;ilist<nlstmax;ilist++) {
			indmax = lstmax[ilist];
			mask[indmax] = 1;
		};
		goto label2;
	};

// Out of the cycle

label3:;

// Compute new ordering

// Count the number of elements in the independent blocks

	for (i=0;i<=_nparts;i++) _blks[i] = 0;

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		if (mask[i] == 0) _blks[icpu+1]++;
	};

	for (i=0;i<_nparts;i++) _blks[i+1] += _blks[i];

	_blks[_nparts+1] = nsupr;

// Set ordering array

	for (i=0;i<nsupr;i++) _order[i] = -1;
	for (i=0;i<_nparts;i++) pblks[i] = _blks[i];

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		if (mask[i] == 0) {
			_order[i] = pblks[icpu];
			pblks[icpu]++;
		};
	};

	nz = pblks[_nparts-1];
	for (i=0;i<nsupr;i++) {
		if (mask[i] == 1) {
			_order[i] = nz;
			nz++;
		};
	};

	int myid;
	myid = CMPIExchange::GetMyidMPIGlobal ();

	if (myid == 0) {
		cout << " Partition: N = " << nsupr << " Np = " << _nparts << " NBini = " << nnini << " NBfin = " << nsupr-_blks[_nparts] << endl;
	};

	delete [] adjwloc;
	delete [] iorder;
	delete [] lstmax;
	delete [] listlc;
	delete [] mask;
	delete [] pblks;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS into not equal parts and compute its ordering
//========================================================================================
void CSMatrix::PartMtr (int _nparts, float *_weights, int *_blks, int *_partition, int *_order) const { // Partition matrix by METIS into not equal parts and compute its ordering

	const char *funcname = "PartMtr_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartMtr (_nparts, _weights, ndweights, adjweights, _blks, _partition, _order);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS and compute its ordering
//========================================================================================
void CSMatrix::PartMtr (int _nparts, int *_blks, int *_partition, int *_order) const { // Partition matrix by METIS and compute its ordering

	const char *funcname = "PartMtr_02";

// Compute weights array

	float *weights;

	weights = new float [_nparts];
	if (!weights) MemoryFail (funcname);

	double aux = 1.0e0 / (double) _nparts;

	for (int i=0;i<_nparts;i++) weights[i] = (float) aux;

// Partition

	PartMtr (_nparts, weights, _blks, _partition, _order);

	delete [] weights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks
//========================================================================================
void CSMatrix::PartOrdMtr (int _nparts, float *_weights, int *_ndweights, int *_adjweights, int *_blks, int *_order) const { // Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks

	const char *funcname = "PartOrdMtr_00";

// Compute profile optimization ordering if _nparts == 1

	if (_nparts == 1) {

		_blks[0] = 0;
		_blks[1] = nsupr;
		_blks[2] = nsupr;

		int ordtyp = 1;

		OrderPrfMtr (ordtyp, _order); 

	} else {

// Allocate working array

		int *partition;

		partition = new int [nsupr];
		if (!partition) MemoryFail (funcname);

// Compute new partition

		PartMtr (_nparts, _weights, _ndweights, _adjweights, _blks, partition, _order); 

// Compute symmetrized matrix

		CSMatrix as;

		as = SymmMtr ();

// Reorder the matrix

		CSMatrix ao;

		ao = as.OrdMtr (_order);

// Compute the number of nonzeroes in each diagonal subblock

		int ndiag, iblk, i, j, jj;
		int *nzblk;

		ndiag = _blks[_nparts];

		nzblk = new int [_nparts];
		if (!nzblk) MemoryFail (funcname);
		for (i=0;i<_nparts;i++) nzblk[i] = 0;

		for (iblk=0;iblk<_nparts;iblk++) {
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ao.ia[i];j<ao.ia[i+1];j++) {
					jj = ao.ja[j];
					if (jj < ndiag) nzblk[iblk]++;
				};
			};
		};

// Compute the maximal number of elements in the diagonal block

		int nmax=0, nzmax=0, nloc, nzloc;

		for (iblk=0;iblk<_nparts;iblk++) {
			nloc = _blks[iblk+1]-_blks[iblk];
			if (nloc > nmax) nmax = nloc;
			if (nzblk[iblk] > nzmax) nzmax = nzblk[iblk];
		};

// For each diagonal subblock create the corresponding main submatrix 
// and compute the corresponding local ordering

		int ishft, ordtyp;

		CSMatrix aloc (nmax,nzmax);

//	ofstream fout ("Submat.dat");

		for (iblk=0;iblk<_nparts;iblk++) {

// Compute local submatrix

			nloc = 0;
			nzloc = 0;
			ishft = _blks[iblk];
			aloc.ia[0] = 0;
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				for (j=ao.ia[i];j<ao.ia[i+1];j++) {
					jj = ao.ja[j];
					if (jj < ndiag) {
						aloc.ja[nzloc] = jj-ishft;
						nzloc++;
					};
				};
				nloc++;
				aloc.ia[nloc] = nzloc;
			};
			aloc.m = nloc;
			aloc.n = nloc;
			aloc.nsupr = nloc;
			aloc.nsupc = nloc;
			aloc.nzja = nzloc;
			aloc.nlist = nloc;

//			fout << " Iblk = " << iblk << aloc;
//			aloc.A2Ps ("Subm.ps",0,_blks);

// Compute local ordering 

			ordtyp = 1;

			aloc.OrderPrfMtr (ordtyp, partition+ishft); 

		};

// Recompute global ordering

		int *iorder;

		iorder = new int [nsupr];
		if (!iorder) MemoryFail (funcname);

		for (i=0;i<nsupr;i++) iorder[_order[i]] = i;

		for (iblk=0;iblk<_nparts;iblk++) {
			ishft = _blks[iblk];
			for (i=_blks[iblk];i<_blks[iblk+1];i++) {
				int iold = iorder[i];
				int inew = ishft+partition[i];
				_order[iold] = inew;
			};
		};

// Delete working arrays

		delete [] iorder;
		delete [] nzblk;
		delete [] partition;

	};

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks
//========================================================================================
void CSMatrix::PartOrdMtr (int _nparts, float *_weights, int *_blks, int *_order) const { // Partition matrix by METIS into not equal parts and compute profile optimization ordering for subblocks

	const char *funcname = "PartOrdMtr_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartOrdMtr (_nparts, _weights, ndweights, adjweights, _blks, _order);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS and compute profile optimization ordering for subblocks
//========================================================================================
void CSMatrix::PartOrdMtr (int _nparts, int *_blks, int *_order) const { // Partition matrix by METIS and compute profile optimization ordering for subblocks

	const char *funcname = "PartOrdMtr_02";

// Compute weights array

	float *weights;

	weights = new float [_nparts];
	if (!weights) MemoryFail (funcname);

	double aux = 1.0e0 / (double) _nparts;

	for (int i=0;i<_nparts;i++) weights[i] = (float) aux;

// Partition

	PartOrdMtr (_nparts, weights, _blks, _order);

	delete [] weights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix according to the prescribed binary tree, matrix ordering is assumed to be nested dessection
//========================================================================================
void CSMatrix::PartBinaryTreeND (CTree &_tree) const { // Partition matrix according to the prescribed binary tree, matrix ordering is assumed to be nested dessection

	const char *funcname = "PartBinaryTreeND";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	_tree.nodes[rootnode].indbeg = 0;
	_tree.nodes[rootnode].indend = nsupr-1;
	_tree.nodes[rootnode].indbegtot = 0;
	_tree.nodes[rootnode].indendtot = nsupr-1;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	CSMatrix mtrdummy;

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichild, isupbeg, isupend, nloc, nzloc;

	int i, j, jj;

	int nbordtot = 0;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

// Take current submatrix

		nloc = isupend-isupbeg+1;
		nzloc = 0;

		for (i=isupbeg;i<=isupend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= isupbeg && jj <= isupend) nzloc++;
			};
		};

		CSMatrix aloc (nloc,nzloc);

		nloc = 0;
		nzloc = 0;
		aloc.ia[0] = 0;
		for (i=isupbeg;i<=isupend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= isupbeg && jj <= isupend) {
					aloc.ja[nzloc] = jj-isupbeg;
					nzloc++;
				};
			};
			aloc.ia[nloc+1] = nzloc;
			nloc++;
		};
		aloc.m = nloc;
		aloc.n = nloc;
		aloc.nsupr = nloc;
		aloc.nsupc = nloc;
		aloc.nzja = nzloc;
		aloc.nlist = nloc;

// Compute separators

		int isep1, isep2;

		aloc.FindSeparator (isep1, isep2);

		if (isep1 == -1) {
//			throw " CSMatrix::PartBinaryTreeND error: separator was not found ";
			isep1 = 0;
			isep2 = 0;
		};

//		cout << " Ncurr = " << nloc << " Isep1 = " << isep1 << " Isep2 = " << isep2 << " Nbord = " << nloc - isep2 << endl;

		nbordtot += (nloc-isep2);

// Free submatrix data

		aloc = mtrdummy;

// Register the partitioning in the child nodes if necessary

		int blksloc[] = {0,isep1,isep2,nloc};

		if (nchildscurr > 1) {
			for (ichild=0;ichild<nchildscurr;ichild++) {
				ichildnode = _tree.nodes[inodecurr].childs[ichild];
				_tree.nodes[ichildnode].indbeg = isupbeg+blksloc[ichild];
				_tree.nodes[ichildnode].indend = isupbeg+blksloc[ichild+1]-1;
				_tree.nodes[ichildnode].indbegtot = isupbeg+blksloc[ichild];
				_tree.nodes[ichildnode].indendtot = isupbeg+blksloc[ichild+1]-1;
			};
			_tree.nodes[inodecurr].indbeg = isupbeg+blksloc[nchildscurr];
			_tree.nodes[inodecurr].indend = isupend;
		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in CSMatrix::PartBinaryTreeND routine";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

//	cout << " NbordTot = " << nbordtot << endl;

// Free work arrays

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS including Schur complement
//========================================================================================
void CSMatrix::PartBinaryTreeNDSchur (CTree &_tree, int *_order, const CSlvParam &_param) const { // Partition matrix by METIS including Schur complement

	const char *funcname = "PartBinaryTreeNDSchur";

// Compute Nested Dissection global ordering

//	ofstream fout ("ChkPart.dat");

	OrderPrfMtr (-1, _order);

//	OutArr (fout," Initial Order = ",n,_order);

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Reorder

	CSMatrix ao;

	ao = as.OrdMtr(_order);

	CSMatrix mtrdummy;

	as = mtrdummy;

// Compute initial partitioning

	ao.PartBinaryTreeND (_tree);

// Compute number of procs for each node

	CTree treetest;

	treetest = _tree;

	treetest.CpuListsForNodes ();

// Allocate the subtree array

	CTree *psubtreearr;

	psubtreearr = _tree.GetSubtreearr ();

	if (psubtreearr != 0) delete [] psubtreearr;

	int nnodesloc = _tree.GetNnodes ();

	psubtreearr = new CTree [nnodesloc];
	if (!psubtreearr) MemoryFail (funcname);

	_tree.SetNsubtree (nnodesloc);
	_tree.SetSubtreearr (psubtreearr);

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;
	int *order;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);
	order = new int [nsupr];
	if (!order) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupr;i++) order[_order[i]] = i;

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int isupbeg, isupend;

	int nsubtree = 0;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

		if (nchildscurr != 1) {

// Take submatrix corresponding to the whole set of indices

			CSMatrix aloc;

			int inode2 = _tree.nodes[inodecurr].childs[0];
			int nchilds2loc = _tree.nodes[inodecurr].nchilds;
			int isupbeg0 = _tree.nodes[inode2].indbeg;

			while (nchilds2loc != 1) {
				int ichild2 = _tree.nodes[inode2].childs[0];
				inode2 = ichild2;
				nchilds2loc = _tree.nodes[inode2].nchilds;
				isupbeg0 = _tree.nodes[inode2].indbeg;
			};

			cout << " Inodecurr = " << inodecurr << " isupbeg0 = " << isupbeg0 << " isupbeg = " << isupbeg << " isupend = " << isupend << endl;
//			fout << " Inodecurr = " << inodecurr << " isupbeg0 = " << isupbeg0 << " isupbeg = " << isupbeg << " isupend = " << isupend << endl;

			ao.Submatrix (isupbeg0, isupend, aloc);

//			char strbuff[256];
//			sprintf(strbuff, "%s%d%s","SubmIni_",inodecurr,".ps");
//			aloc.A2Ps (1,strbuff,0,&i);

// Compute extended Schur complement

			int nbord = isupend+1-isupbeg;

			CSMatrix shloc;

			if (_param.ncycle == 0) throw " CSMatrix::PartBinaryTreeNDSchur: wrong value of parameter ncycle ";

			aloc.ComputeSchur (nbord, _param.ncycle, 
										shloc);

// Compute partitioning of the Schur complement

			int nloc = nbord;

			int *orderloc, *iorderloc;

			orderloc = new int [nloc];
			if (!orderloc) MemoryFail (funcname);
			iorderloc = new int [nloc];
			if (!iorderloc) MemoryFail (funcname);

			if (nloc != 0) {
				shloc.OrderPrfMtr (-1, orderloc);
			};

//			OutArr (fout, " Order loc =",nloc,orderloc);

			for (i=0;i<nloc;i++) iorderloc[orderloc[i]] = i;

			CSMatrix aoloc;

			aoloc = shloc.OrdMtrSymm (orderloc);

//			char strbufff[256];
//			sprintf(strbufff, "%s%d%s","ShlocOrd_",inodecurr,".ps");
//			aoloc.A2Ps (1,strbufff,0,&i);

// Compute the number of processors in the subtree

			int nprocloc = treetest.nodes[inodecurr].nprocnode;
//			int nprocloc = _param.ncpuextsch;

// Create and store local binary tree

			int nchilds = 2;

			double *memory2cpu;
			double *procweight;

			memory2cpu = new double [nprocloc+1];
			if (!memory2cpu) MemoryFail (funcname);
			procweight = new double [nprocloc+1];
			if (!procweight) MemoryFail (funcname);

			for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
			for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

			CTree treeloc (nprocloc, nchilds, procweight, memory2cpu);

			delete [] memory2cpu;
			delete [] procweight;

			treeloc.SortNodes ();

			aoloc.PartBinaryTreeND (treeloc);

			psubtreearr[nsubtree] = treeloc;

			_tree.nodes[inodecurr].indtree = nsubtree;

			nsubtree++;

// Update the ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderloc[i-isupbeg];
				orderloc[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderloc[i-isupbeg];
			};

// Free work arrays

			delete [] orderloc;
			delete [] iorderloc;

		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " CSMatrix::PartBinaryTreeNDSchur: Internal error";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

	_tree.SetNsubtree (nsubtree);

// Delete working arrays

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] order;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by ND METIS
//========================================================================================
void CSMatrix::PartBinaryTreeND (CTree &_tree, int *_order, const CSlvParam &_param, // Partition matrix by ND METIS
														int &_nblks, int *&_blks, int *&_blk2cpu) const {

	const char *funcname = "PartBinaryTreeND";

// Compute Nested Dissection global ordering

	OrderPrfMtr (-1, _order);

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Reorder

	CSMatrix ao;

	ao = as.OrdMtr(_order);

	CSMatrix mtrdummy;

	as = mtrdummy;

// Compute initial partitioning

	ao.PartBinaryTreeND (_tree);

// Compute blocks partitioning

	if (!_tree.IsTreeSorted ()) throw " CSMatrix::PartBinaryTreeND: tree on entry is not sorted ";

	_tree.PartitionTree (_nblks, _blks, _blk2cpu);

// Store blocks partitionings

	int *pblkstree = _tree.GetBlkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();

	delete [] pblkstree;
	delete [] pblk2cputree;

	pblkstree = new int [_nblks+1];
	if (!pblkstree) MemoryFail (funcname);
	pblk2cputree = new int [_nblks];
	if (!pblk2cputree) MemoryFail (funcname);

	int i, j;

	for (i=0;i<=_nblks;i++) pblkstree[i] = _blks[i];
	for (i=0;i<_nblks;i++) pblk2cputree[i] = _blk2cpu[i];

	_tree.SetNblkstree (_nblks);
	_tree.SetBlkstree (pblkstree);
	_tree.SetBlk2cputree (pblk2cputree);

// Compute block sparsity structure

	int *sp2blk;

	sp2blk = new int [nsupr];
	if (!sp2blk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			sp2blk[j] = i;
		};
	};

	CSMatrix ablkstr;

	ablkstr = ao.BlockSparsityDiag (_nblks, sp2blk);

	delete [] sp2blk;

	CSMatrix * pablktree = _tree.GetAblktree ();

	*pablktree = ablkstr.SymmMtr ();

// Extend the block sparsity

	int *piablk = ablkstr.GetIa ();
	int *pjablk = ablkstr.GetJa ();

	int *iabext, *jabext;

	CGSMatrix gmtrdummy;

	gmtrdummy.ExtendBlockSparsityOpt (_nblks,
													piablk, pjablk, iabext, jabext);

	CSMatrix ablkext = ablkstr;

	piablk = ablkext.GetIa ();
	pjablk = ablkext.GetJa ();

	delete [] piablk;
	delete [] pjablk;

	ablkext.SetIa (iabext);
	ablkext.SetJa (jabext);

	ablkext.SetNzja (iabext[_nblks]);

	*pablktree = ablkext.SymmMtr ();

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by ND METIS including Schur complement with boundary local assignments
//========================================================================================
void CSMatrix::PartBinaryTreeNDSchurBnd (CTree &_tree, int *_order, const CSlvParam &_param,
														int &_nblks, int *&_blks, int *&_blk2cpu) const { // Partition matrix by ND METIS including Schur complement with boundary local assignments

//	const char *funcname = "PartBinaryTreeNDSchurBnd";

// Compute Nested Dissection global ordering

//	ofstream fout ("ChkPart.dat");

	OrderPrfMtr (-1, _order);

//	OutArr (fout," Initial Order = ",n,_order);

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Reorder

	CSMatrix ao;

	ao = as.OrdMtr(_order);

	CSMatrix mtrdummy;

	as = mtrdummy;

// Compute initial partitioning

	ao.PartBinaryTreeND (_tree);

// Perform all other computations

	ao.PartBinaryTreeNDSchurBnd_FixedInitialOrdering (_tree, _order, _param,
																		_nblks, _blks, _blk2cpu);

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by ND METIS including Schur complement with boundary local assignments (fixed initial ordering on entry)
//========================================================================================
void CSMatrix::PartBinaryTreeNDSchurBnd_FixedInitialOrdering ( // Partition matrix by ND METIS including Schur complement with boundary local assignments (fixed initial ordering on entry)
																					CTree &_tree, int *_order, const CSlvParam &_param,
																					int &_nblks, int *&_blks, int *&_blk2cpu) const {

	const char *funcname = "PartBinaryTreeNDSchurBnd_FixedInitialOrdering";

// Compute symmetrized matrix

	CSMatrix ao;

	ao = SymmMtr ();

// Compute number of procs for each node

	CTree treetest;

	treetest = _tree;

	treetest.CpuListsForNodes ();

// Allocate the subtree array

	int nnodesloc = _tree.GetNnodes ();

/*
	CTree *psubtreearr;

	psubtreearr = _tree.GetSubtreearr ();

	if (psubtreearr != 0) delete [] psubtreearr;

	psubtreearr = new CTree [nnodesloc];
	if (!psubtreearr) MemoryFail (funcname);

	_tree.SetNsubtree (nnodesloc);
	_tree.SetSubtreearr (psubtreearr);
*/
// Allocate work arrays

	int nprocmax = _tree.nproc;
	int nblksmax = _tree.nnodes + 3*nprocmax*(2*nprocmax-2) + 10;

	int *order;
	int *orderloc;
	int *iorderloc;
	int *blksloc;
	int *blksloc1;
	int *blk2cpu;
	int *imaskcpu2loc;

	order = new int [nsupr];
	if (!order) MemoryFail (funcname);
	orderloc = new int [nsupr];
	if (!orderloc) MemoryFail (funcname);
	iorderloc = new int [nsupr];
	if (!iorderloc) MemoryFail (funcname);
	blksloc = new int [nblksmax+1];
	if (!blksloc) MemoryFail (funcname);
	blksloc1 = new int [nblksmax+1];
	if (!blksloc1) MemoryFail (funcname);
	blk2cpu = new int [nblksmax];
	if (!blk2cpu) MemoryFail (funcname);
	imaskcpu2loc = new int [nprocmax];
	if (!imaskcpu2loc) MemoryFail (funcname);

	int *npartscpu;
	int *ibspartcpu;
	int *nparts1cpu;
	int *ibspart1cpu;
	int *imaskcpu;
	int *listcpu;
	int *imaskind;
	int *listind;
	int *ind2list;
	int *list1cpu;
	int *list2cpu;
	int *nlistcpupairs;
	int *ibscpupairs;

	npartscpu = new int [nprocmax];
	if (!npartscpu) MemoryFail (funcname);
	ibspartcpu = new int [nprocmax];
	if (!ibspartcpu) MemoryFail (funcname);
	nparts1cpu = new int [nprocmax];
	if (!npartscpu) MemoryFail (funcname);
	ibspart1cpu = new int [nprocmax];
	if (!ibspartcpu) MemoryFail (funcname);
	imaskcpu = new int [nprocmax];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nprocmax];
	if (!listcpu) MemoryFail (funcname);
	imaskind = new int [nprocmax*nprocmax];
	if (!imaskind) MemoryFail (funcname);
	listind = new int [nprocmax*nprocmax];
	if (!listind) MemoryFail (funcname);
	ind2list = new int [nprocmax*nprocmax];
	if (!ind2list) MemoryFail (funcname);
	list1cpu = new int [nprocmax*nprocmax];
	if (!list1cpu) MemoryFail (funcname);
	list2cpu = new int [nprocmax*nprocmax];
	if (!list2cpu) MemoryFail (funcname);
	nlistcpupairs = new int [nprocmax*nprocmax];
	if (!nlistcpupairs) MemoryFail (funcname);
	ibscpupairs = new int [nprocmax*nprocmax];
	if (!ibscpupairs) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupr;i++) order[_order[i]] = i;

	for (i=0;i<nsupr;i++) orderloc[i] = i;
	for (i=0;i<nsupr;i++) iorderloc[i] = i;

	int icyclecpu = -1;

	for (i=0;i<nprocmax;i++) imaskcpu[i] = icyclecpu;

	int icycleind = -1;

	for (i=0;i<nprocmax*nprocmax;i++) imaskind[i] = icycleind;

// Search the tree

	int inodecurr, nchildscurr;
	int isupbeg, isupend, iproc, k, jj, jjord, iblk, iblkbeg, iblkend, ncpuloc;
	int j, jproc, iiproc, jjproc, ind, ind1, nz, ilistind, ibs, nmatrloc, indsep1, indsep2;
	int niloc, ni1, ni2, ni3, j1, j2, j3, ishft, jjold, jjoldind, nblksfirst, ibs0;
	int iiprocn, jjprocn;

//	int nsubtree = 0;
	int nblksloc = 0;
	int nblksloc0, nblksloc1;

	CSMatrix mtrdummy;

	blksloc[0] = 0;

	for (inodecurr=0;inodecurr<nnodesloc;inodecurr++) {

// Take data of the node

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;
		iproc = _tree.nodes[inodecurr].nodecpu;

		nblksloc0 = nblksloc;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

		if (nchildscurr == 1) {

			blk2cpu[nblksloc] = iproc;
			blksloc[nblksloc+1] = isupend+1;
			nblksloc++;

		} else {

// Take submatrix corresponding to the whole set of indices

			CSMatrix aloc;

			int inode2 = _tree.nodes[inodecurr].childs[0];
			int nchilds2loc = _tree.nodes[inodecurr].nchilds;
			int isupbeg0 = _tree.nodes[inode2].indbeg;

			while (nchilds2loc != 1) {
				int ichild2 = _tree.nodes[inode2].childs[0];
				inode2 = ichild2;
				nchilds2loc = _tree.nodes[inode2].nchilds;
				isupbeg0 = _tree.nodes[inode2].indbeg;
			};

//			cout << " Inodecurr = " << inodecurr << " isupbeg0 = " << isupbeg0 << " isupbeg = " << isupbeg << " isupend = " << isupend << endl;
//			cout << "     niloc = " << isupend-isupbeg+1 << endl;
//			fout << " Inodecurr = " << inodecurr << " isupbeg0 = " << isupbeg0 << " isupbeg = " << isupbeg << " isupend = " << isupend << endl;

			ao.Submatrix (isupbeg0, isupend, aloc);

//			char strbuff[256];
//			sprintf(strbuff, "%s%d%s","SubmIni_",inodecurr,".ps");
//			aloc.A2Ps (1,strbuff,0,&i);

// Compute extended Schur complement

			int nbord = isupend+1-isupbeg;

			CSMatrix shloc;

			if (_param.ncycle < 1) throw " CSMatrix::PartBinaryTreeNDSchurBnd: wrong value of parameter ncycle ";

			aloc.ComputeSchur (nbord, _param.ncycle, 
										shloc);
//			cout << " nbord*7 = " << nbord*7 << " nzja schur " << shloc.GetNzja () << endl;

// Find local rows submatrix

			CSMatrix arowsloc;

			ao.SubmatrixByRows (isupbeg, isupend, arowsloc);

// Create for each row index the list of processors that references to it

			int nzloc = arowsloc.GetNzja ();

			int *jarowsord;
			int *jarows2cpu;

			jarowsord = new int [nzloc];
			if (!jarowsord) MemoryFail (funcname);
			jarows2cpu = new int [nzloc];
			if (!jarows2cpu) MemoryFail (funcname);

			int nzjarows = arowsloc.GetNzja ();
			int *pjarows = arowsloc.GetJa ();

			int iblkprev = 0;

			for (k=0;k<nzjarows;k++) {
				jj = pjarows[k];
				jjord = orderloc[jj];
				jarowsord[k] = jjord;
				if (jjord >= isupbeg) {
					jarows2cpu[k] = iproc;
				} else {
					if (jjord >= blksloc[iblkprev] && jjord < blksloc[iblkprev+1]) {
						iblk = iblkprev;
					} else {
						iblkbeg = 0;
						iblkend = nblksloc-1;
						while (iblkbeg != iblkend) {
							iblk = (iblkbeg+iblkend) / 2;
							if (jjord >= blksloc[iblk] && jjord < blksloc[iblk+1]) {
								iblkbeg = iblk;
								iblkend = iblk;
							} else if (jjord < blksloc[iblk]) {
								iblkend = iblk-1;
							} else if (jjord >= blksloc[iblk+1]) {
								iblkbeg = iblk+1;
							};
						};
						iblk = iblkbeg;
						iblkprev = iblk;
					};
					jarows2cpu[k] = blk2cpu[iblk];
				};
			};

// Separate all row indices that reference to a pair of processors

			int *row2ind;

			row2ind = new int [nbord];
			if (!row2ind) MemoryFail (funcname);

			int *piarows = arowsloc.GetIa ();

			icycleind++;
			int nlistind = 0;

			for (i=0;i<nbord;i++) {
				icyclecpu++;
				ncpuloc = 0;
				row2ind[i] = -1;
				for (j=piarows[i];j<piarows[i+1];j++) {
					jjord = jarowsord[j];
					jproc = jarows2cpu[j];
					if (jjord < isupbeg) {
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[ncpuloc] = jproc;
							ncpuloc++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
				if (ncpuloc < 1) {
//					cout << " Before throw: CSMatrix::PartBinaryTreeNDSchurBnd: there are rows with no reference " << endl;
//					throw " CSMatrix::PartBinaryTreeNDSchurBnd: there are rows with no reference ";
				};
				if (ncpuloc == 1) {
					iiproc = listcpu[0];
					ind = iiproc*nprocmax+iiproc;
					if (imaskind[ind] != icycleind) {
						nlistcpupairs[ind] = 0;
						listind[nlistind] = ind;
						list1cpu[nlistind] = iiproc;
						list2cpu[nlistind] = iiproc;
						ind2list[ind] = nlistind;
						nlistind++;
						imaskind[ind] = icycleind;
					};
					nlistcpupairs[ind]++;
					row2ind[i] = ind;
				};
				if (ncpuloc == 2) {
					iiproc = listcpu[0];
					jjproc = listcpu[1];
					ind = iiproc*nprocmax+jjproc;
					ind1 = jjproc*nprocmax+iiproc;
					if (imaskind[ind] != icycleind) {
						nlistcpupairs[ind] = 0;
						nlistcpupairs[ind1] = 0;
						listind[nlistind] = ind;
						list1cpu[nlistind] = iiproc;
						list2cpu[nlistind] = jjproc;
						ind2list[ind] = nlistind;
						ind2list[ind1] = nlistind;
						nlistind++;
						imaskind[ind] = icycleind;
						imaskind[ind1] = icycleind;
					};
					nlistcpupairs[ind]++;
					nlistcpupairs[ind1]++;
					row2ind[i] = ind;
				};
			};

// Init the counters of the cpu pairs lists

			nz = 0;

			for (i=0;i<nlistind;i++) {
				ind = listind[i];
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				ind = iiproc*nprocmax+jjproc;
				ind1 = jjproc*nprocmax+iiproc;
				ibscpupairs[ind] = nz;
				ibscpupairs[ind1] = nz;
				nz += nlistcpupairs[ind];
				nlistcpupairs[ind] = 0;
				nlistcpupairs[ind1] = 0;
			};

// Store the lists

			int *listpairs;

			listpairs = new int [nbord];
			if (!listpairs) MemoryFail (funcname);

			for (i=0;i<nbord;i++) {
				if (row2ind[i] >= 0) {
					ind = row2ind[i];
					ilistind = ind2list[ind];
					iiproc = list1cpu[ilistind];
					jjproc = list2cpu[ilistind];
					ind = iiproc*nprocmax+jjproc;
					ind1 = jjproc*nprocmax+iiproc;
					ibs = ibscpupairs[ind]+nlistcpupairs[ind];
					listpairs[ibs] = i;
					nlistcpupairs[ind]++;
					if (iiproc != jjproc) nlistcpupairs[ind1]++;
				};
			};

// Sort data in the lists

			for (i=0;i<nlistind;i++) {
				ind = listind[i];
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				ind = iiproc*nprocmax+jjproc;
				ind1 = jjproc*nprocmax+iiproc;
				ibs = ibscpupairs[ind];
				qsort (listpairs+ibs, nlistcpupairs[ind], sizeof(int), compint);
			};

// Perform filtering of the data in the lists according to the sparsity of the Schur complement

			int *imasksubm;
			int *irow2indsubm;

			imasksubm = new int [nbord];
			if (!imasksubm) MemoryFail (funcname);
			irow2indsubm = new int [nbord];
			if (!irow2indsubm) MemoryFail (funcname);

			int icyclesubm = -1;

			for (i=0;i<nbord;i++) imasksubm[i] = icyclesubm;

			int *piash = shloc.GetIa ();
			int *pjash = shloc.GetJa ();

			int nn, icheck, kk;

			for (i=0;i<nlistind;i++) {
				icyclesubm++;
				ind = listind[i];
				ibs = ibscpupairs[ind];
				nmatrloc = nlistcpupairs[ind];
				for (j=0;j<nmatrloc;j++) {
					jj = listpairs[ibs+j];
					imasksubm[jj] = icyclesubm;
					irow2indsubm[j] = jj;
				};
				nn = 0;
				for (j=0;j<nmatrloc;j++) {
					jj = irow2indsubm[j];
					icheck = 1;
					for (k=piash[jj];k<piash[jj+1];k++) {
						kk = pjash[k];
//						if (imasksubm[kk] != icyclesubm) icheck = -1;
						if (imasksubm[kk] != icyclesubm && row2ind[kk] > 0) icheck = -1;
					};
					if (icheck == 1) {
						listpairs[ibs+nn] = jj;
						nn++;
					} else {
						row2ind[jj] = -1;
					};
				};
//				cout << "     Nlist bef flt = " << nlistcpupairs[ind] << " After = " << nn << endl;
				nlistcpupairs[ind] = nn;
			};

// Scan all the lists and try to split the lists into three parts and 
// store the corresponding ordering and partitioning

			int nprocloc = _tree.nodes[inodecurr].nprocnode;

			int *blkscpuarr;

			blkscpuarr = new int [nprocloc+2];
			if (!blkscpuarr) MemoryFail (funcname);

			for (i=0;i<nprocloc+2;i++) blkscpuarr[i] = 0;

			int *indsep1arr;
			int *indsep2arr;
			int *orderpairs;
			int *iorderpairs;

			indsep1arr = new int [nlistind];
			if (!indsep1arr) MemoryFail (funcname);
			indsep2arr = new int [nlistind];
			if (!indsep2arr) MemoryFail (funcname);
			orderpairs = new int [nbord];
			if (!orderpairs) MemoryFail (funcname);
			iorderpairs = new int [nbord];
			if (!iorderpairs) MemoryFail (funcname);

			CSMatrix shsubm, sho, shsubmsymm;

			int nmatrmax = 0;

			for (i=0;i<nlistind;i++) {
				ind = listind[i];
				nmatrloc = nlistcpupairs[ind];
				if (nmatrloc > nmatrmax) nmatrmax = nmatrloc;
			};

			nmatrmax = nmatrmax / 10;
			if (nmatrmax < 1) nmatrmax = 1;

			for (i=0;i<nlistind;i++) {

// Get submatrix

				ind = listind[i];
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				ind = iiproc*nprocmax+jjproc;
				ind1 = jjproc*nprocmax+iiproc;
				ibs = ibscpupairs[ind];
				nmatrloc = nlistcpupairs[ind];
				shloc.SubmatrixList (nmatrloc, listpairs+ibs, 
											icyclesubm, imasksubm, irow2indsubm,
											shsubm);

// Compute new local ordering

				if (nmatrloc != 0 && nmatrloc >= nmatrmax && iiproc != jjproc) {

//					cout << "    Iproc = " << iiproc << " Jproc = " << jjproc << endl;
					shsubmsymm = shsubm.SymmMtr ();

					shsubmsymm.OrderPrfMtr (-1, orderpairs+ibs);

// Compute and store new partitioning

					sho = shsubmsymm.OrdMtr (orderpairs+ibs);

					sho.FindSeparator (indsep1, indsep2);

					if (indsep1 < 0 || indsep2 < 0) {
						indsep1 = 0;
						indsep2 = 0;
					};

//					if (true) {
//						char strbuff[256];
//						int blkloc[] = {0,indsep1,indsep2,nmatrloc};
//						sprintf (strbuff, "%s%i%s%i%s","Subm_",inodecurr,"_nd_",i,"_lst.ps");
//						sho.A2Ps (1,strbuff,3,blkloc);
//					};

				} else {
					indsep1 = 0;
					indsep2 = 0;
					for (j=0;j<nmatrloc;j++) orderpairs[ibs+j] = j;
				};

				if (indsep1 != 0 && indsep2 != 0) {
//					cout << "       Ni1 = " << indsep1 << " Ni2 = " << indsep2-indsep1 << " Nbnd = " << nmatrloc-indsep2 << endl;
				};

				indsep1arr[i] = indsep1;
				indsep2arr[i] = indsep2;

// Inverse local order

				for (j=0;j<nmatrloc;j++) {
					iorderpairs[ibs+orderpairs[ibs+j]] = j;
				};

			};

			shsubm = mtrdummy;
			sho = mtrdummy;

// For each subpart mark its place in the global Schur ordering and partitioning

			for (i=0;i<_tree.nodes[inodecurr].nprocnode;i++) {
				iiproc = _tree.nodes[inodecurr].cpuidarrnode[i];
				imaskcpu2loc[iiproc] = i;
				npartscpu[iiproc] = 0;
				ibspartcpu[iiproc] = 0;
				nparts1cpu[iiproc] = 0;
				ibspart1cpu[iiproc] = 0;
			};

			for (i=0;i<nlistind;i++) {
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				npartscpu[iiproc] += 1;
				npartscpu[jjproc] += 1;
				nparts1cpu[iiproc] += 1;
			};

			nz = 0;
			for (i=0;i<_tree.nodes[inodecurr].nprocnode;i++) {
				iiproc = _tree.nodes[inodecurr].cpuidarrnode[i];
				ibspartcpu[iiproc] = nz;
				nz += npartscpu[iiproc];
				npartscpu[iiproc] = 0;
			};
			nblksfirst = nz;
			for (i=0;i<_tree.nodes[inodecurr].nprocnode;i++) {
				iiproc = _tree.nodes[inodecurr].cpuidarrnode[i];
				ibspart1cpu[iiproc] = nz;
				nz += nparts1cpu[iiproc];
				nparts1cpu[iiproc] = 0;
			};

			int *orderpart;

			orderpart = new int [nlistind*3];
			if (!orderpart) MemoryFail (funcname);

			for (i=0;i<nlistind;i++) {
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				ibs = ibspartcpu[iiproc]+npartscpu[iiproc];
				orderpart[i*3] = ibs;
				npartscpu[iiproc]++;
				ibs = ibspartcpu[jjproc]+npartscpu[jjproc];
				orderpart[i*3+1] = ibs;
				npartscpu[jjproc]++;
				ibs = ibspart1cpu[iiproc]+nparts1cpu[iiproc];
				orderpart[i*3+2] = ibs;
				nparts1cpu[iiproc]++;
			};

			int *sprndspart;

			sprndspart = new int [nlistind*3+1];
			if (!sprndspart) MemoryFail (funcname);

			for (i=0;i<=nlistind*3;i++) sprndspart[i] = 0;

			for (i=0;i<nlistind;i++) {
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				iiprocn = imaskcpu2loc[iiproc];
				jjprocn = imaskcpu2loc[jjproc];
				ind = listind[i];
				niloc = nlistcpupairs[ind];
				indsep1 = indsep1arr[i];
				indsep2 = indsep2arr[i];
				ni1 = indsep1;
				ni2 = indsep2-indsep1;
				ni3 = niloc - (ni1+ni2);
				j1 = orderpart[i*3];
				j2 = orderpart[i*3+1];
				j3 = orderpart[i*3+2];
				sprndspart[j1+1] = ni1;
				sprndspart[j2+1] = ni2;
				sprndspart[j3+1] = ni3;
				blkscpuarr[iiprocn+1] += ni1;
				blkscpuarr[jjprocn+1] += ni2;
			};

			for (i=0;i<nlistind*3;i++) sprndspart[i+1] = sprndspart[i]+sprndspart[i+1];

// Compute global partitioning and distribution

			if (nblksloc+3*nlistind+1 > nblksmax) throw " CSMatrix::PartBinaryTreeNDSchurBnd: insufficient number of blocks reserved ";

			for (i=0;i<nlistind;i++) {
				iiproc = list1cpu[i];
				jjproc = list2cpu[i];
				j1 = orderpart[i*3];
				j2 = orderpart[i*3+1];
				j3 = orderpart[i*3+2];
				blk2cpu[nblksloc+j1] = iiproc;
				blk2cpu[nblksloc+j2] = jjproc;
//				blk2cpu[nblksloc+j3] = iiproc;
				blk2cpu[nblksloc+j3] = iproc;
			};

			for (i=0;i<3*nlistind;i++) {
				blksloc[nblksloc+1] = isupbeg+sprndspart[i+1];
				nblksloc++;
			};

			if (blksloc[nblksloc] != isupend+1) {

				blk2cpu[nblksloc] = iproc;
				blksloc[nblksloc+1] = isupend+1;
				nblksloc++;

			};

// Pack blocks partitioning

			nblksloc = nblksloc0+nblksfirst;
			blk2cpu[nblksloc] = iproc;
			blksloc[nblksloc+1] = isupend+1;

			nblksloc++;

// Compute global inverse ordering of the Schur complement

			for (i=0;i<nlistind;i++) {
				ind = listind[i];
				niloc = nlistcpupairs[ind];
				indsep1 = indsep1arr[i];
				indsep2 = indsep2arr[i];
				ni1 = indsep1;
				ni2 = indsep2-indsep1;
				ni3 = niloc - (ni1+ni2);
				ibs = ibscpupairs[ind];
				ibs0 = ibs;
				j1 = orderpart[i*3];
				j2 = orderpart[i*3+1];
				j3 = orderpart[i*3+2];
				ishft = sprndspart[j1];
				for (j=0;j<ni1;j++) {
					jjold = iorderpairs[ibs+j];
					jjoldind = listpairs[ibs0+jjold];
					iorderloc[isupbeg+ishft+j] = jjoldind;
				};
				ibs += ni1;
				ishft = sprndspart[j2];
				for (j=0;j<ni2;j++) {
					jjold = iorderpairs[ibs+j];
					jjoldind = listpairs[ibs0+jjold];
					iorderloc[isupbeg+ishft+j] = jjoldind;
				};
				ibs += ni2;
				ishft = sprndspart[j3];
				for (j=0;j<ni3;j++) {
					jjold = iorderpairs[ibs+j];
					jjoldind = listpairs[ibs0+jjold];
					iorderloc[isupbeg+ishft+j] = jjoldind;
				};
			};

			ishft = sprndspart[nlistind*3];

			for (i=0;i<nbord;i++) {
				if (row2ind[i] < 0) {
					iorderloc[isupbeg+ishft] = i;
					ishft++;
				};
			};

// Inverse the local ordering in the complement

			for (i=0;i<nbord;i++) orderloc[isupbeg+iorderloc[isupbeg+i]] = isupbeg+i;

			for (i=0;i<nbord;i++) iorderloc[isupbeg+i] += isupbeg;

// Compute local blocks partitioning

//			for (i=0;i<nprocloc;i++) {
//				iiproc = _tree.nodes[inodecurr].cpuidarrnode[i];
//				ibs = ibspartcpu[iiproc]+npartscpu[iiproc];
//				blkscpuarr[i+1] = isupbeg+sprndspart[ibs+1];
//			};

			for (i=0;i<nprocloc;i++) blkscpuarr[i+1] = blkscpuarr[i]+blkscpuarr[i+1];
			for (i=0;i<nprocloc+1;i++) blkscpuarr[i] += isupbeg;

			blkscpuarr[nprocloc+1] = isupbeg+nbord;

// Perform final filtering of zero size blocks

			nblksloc1 = nblksloc0;
			blksloc1[nblksloc1] = blksloc[nblksloc1];

			for (i=nblksloc0;i<nblksloc;i++) {
				niloc = blksloc[i+1]-blksloc[i];
				if (niloc > 0) {
					blksloc1[nblksloc1+1] = blksloc[i+1];
					blk2cpu[nblksloc1] = blk2cpu[i];
					nblksloc1++;
				};
			};

			for (i=nblksloc0;i<=nblksloc1;i++) blksloc[i] = blksloc1[i];

			nblksloc = nblksloc1;

// Add zero size block if necessary

			if (nblksloc == nblksloc0) {
				blk2cpu[nblksloc] = iproc;
				blksloc[nblksloc+1] = isupend+1;
				nblksloc++;
			};

// Replace blkscpuarr by the block numbers
/*
			int irow;

			ind = nblksloc0;

			for (i=0;i<=nprocloc+1;i++) {
				irow = blkscpuarr[i];
				while (ind<nblksloc && blksloc[ind] != irow) {
					ind++;
				};
				if (irow != blksloc[ind]) throw " CSMatrix::PartBinaryTreeNDSchurBnd: block number not found";
				blkscpuarr[i] = ind;
			};

// Create and store local binary tree

			int nchilds = nprocloc;

			double *memory2cpu;
			double *procweight;

			memory2cpu = new double [nprocloc+1];
			if (!memory2cpu) MemoryFail (funcname);
			procweight = new double [nprocloc+1];
			if (!procweight) MemoryFail (funcname);

			for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
			for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

			CTree treeloc (nprocloc, nchilds, procweight, memory2cpu);

			delete [] memory2cpu;
			delete [] procweight;

			treeloc.SortNodes ();

			int nnodessub = treeloc.GetNnodes ();
			CNode *pnodessub = treeloc.GetNode ();

			for (i=0;i<nnodessub;i++) {
				pnodessub[i].indbeg = blkscpuarr[i];
				pnodessub[i].indbegtot = blkscpuarr[i];
				pnodessub[i].indend = blkscpuarr[i+1]-1;
				pnodessub[i].indendtot = blkscpuarr[i+1]-1;
			};

			pnodessub[nnodessub-1].indbegtot = nblksloc0;

			psubtreearr[nsubtree] = treeloc;

			_tree.nodes[inodecurr].indtree = nsubtree;

			nsubtree++;
*/
// Free work arrays

			delete [] jarowsord;
			delete [] jarows2cpu;
			delete [] row2ind;
			delete [] listpairs;
			delete [] indsep1arr;
			delete [] indsep2arr;
			delete [] imasksubm;
			delete [] irow2indsubm;
			delete [] orderpairs;
			delete [] iorderpairs;
			delete [] blkscpuarr;
			delete [] sprndspart;
			delete [] orderpart;
		};

// Replace index numbers by block numbers

		_tree.nodes[inodecurr].indbeg = nblksloc0;
		_tree.nodes[inodecurr].indend = nblksloc-1;

		int inode2 = _tree.nodes[inodecurr].childs[0];
		int nchilds2loc = _tree.nodes[inodecurr].nchilds;
		int iblkbeg0 = _tree.nodes[inode2].indbeg;

		while (nchilds2loc != 1) {
			int ichild2 = _tree.nodes[inode2].childs[0];
			inode2 = ichild2;
			nchilds2loc = _tree.nodes[inode2].nchilds;
			iblkbeg0 = _tree.nodes[inode2].indbeg;
		};

		_tree.nodes[inodecurr].indbegtot = iblkbeg0;
		_tree.nodes[inodecurr].indendtot = nblksloc-1;

	};

// Reorder the matrix once again

	ao = OrdMtr (orderloc);

// Compute inverse global order

	for (i=0;i<nsupr;i++) order[_order[i]] = i;

	int iold;

	for (i=0;i<nsupr;i++) {
		iold = iorderloc[i];
		orderloc[i] = order[iold];
	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[orderloc[i]] = i;

//	_tree.SetNsubtree (nsubtree);

// Store blocks partitionings

	_nblks = nblksloc;

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);
	_blk2cpu = new int [_nblks];
	if (!_blk2cpu) MemoryFail (funcname);

	for (i=0;i<=_nblks;i++) _blks[i] = blksloc[i];
	for (i=0;i<_nblks;i++) _blk2cpu[i] = blk2cpu[i];

	int *pblkstree = _tree.GetBlkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();

	delete [] pblkstree;
	delete [] pblk2cputree;

	pblkstree = new int [_nblks+1];
	if (!pblkstree) MemoryFail (funcname);
	pblk2cputree = new int [_nblks];
	if (!pblk2cputree) MemoryFail (funcname);

	for (i=0;i<=_nblks;i++) pblkstree[i] = _blks[i];
	for (i=0;i<_nblks;i++) pblk2cputree[i] = _blk2cpu[i];

	_tree.SetNblkstree (_nblks);
	_tree.SetBlkstree (pblkstree);
	_tree.SetBlk2cputree (pblk2cputree);

// Compute block sparsity structure

	int *sp2blk;

	sp2blk = new int [nsupr];
	if (!sp2blk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			sp2blk[j] = i;
		};
	};

	CSMatrix ablkstr;

	ablkstr = ao.BlockSparsityDiag (_nblks, sp2blk);

//	ablkstr.A2Ps (1, "BlkStr.ps",0, &i);

	delete [] sp2blk;

	CSMatrix * pablktree = _tree.GetAblktree ();

	*pablktree = ablkstr.SymmMtr ();

// Extend the block sparsity

	int *piablk = ablkstr.GetIa ();
	int *pjablk = ablkstr.GetJa ();

	int *iabext, *jabext;

	CGSMatrix gmtrdummy;

	gmtrdummy.ExtendBlockSparsityOpt (_nblks,
													piablk, pjablk, iabext, jabext);

	CSMatrix ablkext = ablkstr;

	piablk = ablkext.GetIa ();
	pjablk = ablkext.GetJa ();

	delete [] piablk;
	delete [] pjablk;

	ablkext.SetIa (iabext);
	ablkext.SetJa (jabext);

	ablkext.SetNzja (iabext[_nblks]);

//	ablkext.A2Ps (1, "BlkStrExt.ps",0, &i);

	if (_param.schurtype > 1) {
		*pablktree = ablkext.SymmMtr ();
	};

// Free working arrays

	delete [] order;
	delete [] orderloc;
	delete [] iorderloc;
	delete [] blksloc;
	delete [] blksloc1;
	delete [] blk2cpu;
	delete [] imaskcpu2loc;
	delete [] npartscpu;
	delete [] ibspartcpu;
	delete [] nparts1cpu;
	delete [] ibspart1cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] imaskind;
	delete [] listind;
	delete [] ind2list;
	delete [] list1cpu;
	delete [] list2cpu;
	delete [] nlistcpupairs;
	delete [] ibscpupairs;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks
//========================================================================================
void CSMatrix::PartOrdMtr (CTree &_tree, int *_ndweights, int *_adjweights, int *_order) const { // Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks

	const char *funcname = "PartOrdMtr_03";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local nodes/adjacency weights

	int *ndwarr;
	int *adjwarr;

	ndwarr = new int [nsupr];
	if (!ndwarr) MemoryFail (funcname);
	adjwarr = new int [as.nzja];
	if (!adjwarr) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<nsupr;i++) ndwarr[i] = _ndweights[i];

	for (i=0;i<as.nzja;i++) adjwarr[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwarr[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwarr[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwarr[j] = _adjweights[k];
					};
				};
			};
			if (adjwarr[j] < 0) throw " PartMtr_00: index was not found or incorrect weights array";
		};
	};

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;
	int *iorder;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	_tree.nodes[rootnode].indbeg = 0;
	_tree.nodes[rootnode].indend = nsupr-1;
	_tree.nodes[rootnode].indbegtot = 0;
	_tree.nodes[rootnode].indendtot = nsupr-1;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	CSMatrix mtrdummy;

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichild, isupbeg, isupend, nloc, nzloc;

	for (i=0;i<nsupr;i++) _order[i] = i;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

// Take current submatrix

		nloc = isupend-isupbeg+1;
		nzloc = 0;

		for (i=isupbeg;i<=isupend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= isupbeg && jj <= isupend) nzloc++;
			};
		};

		CSMatrix aloc (nloc,nzloc);

		int *ndwloc;
		int *adjwloc;

		ndwloc = new int [nloc];
		if (!ndwloc) MemoryFail (funcname);
		adjwloc = new int [nzloc];
		if (!adjwloc) MemoryFail (funcname);

		nloc = 0;
		nzloc = 0;
		aloc.ia[0] = 0;
		for (i=isupbeg;i<=isupend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= isupbeg && jj <= isupend) {
					aloc.ja[nzloc] = jj-isupbeg;
					adjwloc[nzloc] = adjwarr[j];
					nzloc++;
				};
			};
			aloc.ia[nloc+1] = nzloc;
			ndwloc[nloc] = ndwarr[i];
			nloc++;
		};
		aloc.m = nloc;
		aloc.n = nloc;
		aloc.nsupr = nloc;
		aloc.nsupc = nloc;
		aloc.nzja = nzloc;
		aloc.nlist = nloc;

// Compute array weights

		float *weights;

		weights = new float [nchildscurr+2];
		if (!weights) MemoryFail (funcname);

		for (i=0;i<nchildscurr;i++) weights[i] = 0.0e0;

		for (i=0;i<nchildscurr;i++) {
			ichildnode = _tree.nodes[inodecurr].childs[i];
//			int ichildcpu = _tree.nodes[ichildnode].nodecpu;
			weights[i] = (float) _tree.nodes[ichildnode].nodeperf;
		};

		float faux = 0.0e0;

		for (i=0;i<nchildscurr;i++) faux += weights[i];

		for (i=0;i<nchildscurr;i++) weights[i] /= faux;

// Compute new ordering

		int *blksloc, *partition, *orderloc, *iorderloc;

		blksloc = new int [nchildscurr+2];
		if (!blksloc) MemoryFail (funcname);
		partition = new int [nloc];
		if (!partition) MemoryFail (funcname);
		orderloc = new int [nloc];
		if (!orderloc) MemoryFail (funcname);
		iorderloc = new int [nloc];
		if (!iorderloc) MemoryFail (funcname);

		if (nchildscurr == 1) {
			int ordtyp=1;
			aloc.OrderPrfMtr (ordtyp, orderloc);
		} else {
			if (nloc >= nchildscurr+1) {
				aloc.PartMtr (nchildscurr, weights, ndwloc, adjwloc, blksloc, partition, orderloc);
			} else {
				for (i=0;i<=nchildscurr;i++) blksloc[i] = 0;
				blksloc[nchildscurr+1] = nloc;
				for (i=0;i<nloc;i++) partition[i] = 0;
				for (i=0;i<nloc;i++) orderloc[i] = i;
			};
		};

		for (i=0;i<nloc;i++) iorderloc[orderloc[i]] = i;

		delete [] ndwloc;
		delete [] adjwloc;
		delete [] weights;

// Free submatrix data

		aloc = mtrdummy;

// Update the ordering

		for (i=isupbeg;i<=isupend;i++) {
			int iold = iorderloc[i-isupbeg];
//			orderloc[iold] = _order[i];
			orderloc[i-isupbeg] = _order[iold+isupbeg];
		};

		for (i=isupbeg;i<=isupend;i++) {
			_order[i] = orderloc[i-isupbeg];
		};

		for (i=0;i<nloc;i++) orderloc[iorderloc[i]] = i;

// Reorder current part of the matrix in place

		int *ialoc, *jaloc;

		nzloc = 0;
		for (i=isupbeg;i<=isupend;i++) {
			nzloc += as.ia[i+1]-as.ia[i];
		};

		ialoc = new int [nloc+1];
		if (!ialoc) MemoryFail (funcname);
		jaloc = new int [nzloc];
		if (!jaloc) MemoryFail (funcname);

		ndwloc = new int [nloc];
		if (!ndwloc) MemoryFail (funcname);
		adjwloc = new int [nzloc];
		if (!adjwloc) MemoryFail (funcname);

		nzloc = 0;

		ialoc[0] = 0;

		for (i=isupbeg;i<=isupend;i++) {
			int nzloc0 = nzloc;
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj>=isupbeg && jj<=isupend) {
					int jjloc = jj-isupbeg;
					int jjnew = orderloc[jjloc];
					jj = jjnew + isupbeg;
				};
				jaloc[nzloc] = jj;
				adjwloc[nzloc] = adjwarr[j];
				nzloc++;
			};
//			qsort (jaloc+nzloc0, nzloc-nzloc0, sizeof(int), compint);
			SortElm (nzloc-nzloc0, jaloc+nzloc0, adjwloc+nzloc0);
			ialoc[i-isupbeg+1] = nzloc;
			ndwloc[i-isupbeg] = ndwarr[i];
		};

		nzloc = as.ia[isupbeg];

		for (i=isupbeg;i<=isupend;i++) {
			int iold = iorderloc[i-isupbeg];
			for (j=ialoc[iold];j<ialoc[iold+1];j++) {
				jj = jaloc[j];
				as.ja[nzloc] = jj;
				adjwarr[nzloc] = adjwloc[j];
				nzloc++;
			};
			as.ia[i+1] = nzloc;
			ndwarr[i] = ndwloc[iold];
		};

		delete [] ialoc;
		delete [] jaloc;
		delete [] ndwloc;
		delete [] adjwloc;

//		as.A2Ps (100,"LocalOrd.ps",0,&i);

// Register the partitioning in the child nodes if necessary

		if (nchildscurr > 1) {
			for (ichild=0;ichild<nchildscurr;ichild++) {
				ichildnode = _tree.nodes[inodecurr].childs[ichild];
				_tree.nodes[ichildnode].indbeg = isupbeg+blksloc[ichild];
				_tree.nodes[ichildnode].indend = isupbeg+blksloc[ichild+1]-1;
				_tree.nodes[ichildnode].indbegtot = isupbeg+blksloc[ichild];
				_tree.nodes[ichildnode].indendtot = isupbeg+blksloc[ichild+1]-1;
			};
			_tree.nodes[inodecurr].indbeg = isupbeg+blksloc[nchildscurr];
			_tree.nodes[inodecurr].indend = isupend;
		};

// Free local arrays

		delete [] blksloc;
		delete [] partition;
		delete [] orderloc;
		delete [] iorderloc;

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in PartOrdMtr routine";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) iorder[_order[i]] = i;
	for (i=0;i<nsupr;i++) _order[i] = iorder[i];

// Free work arrays

	delete [] adjwarr;
	delete [] ndwarr;

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] iorder;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks
//========================================================================================
void CSMatrix::PartOrdMtr (CTree &_tree, int *_order) const { // Partition matrix by METIS according to the prescribed matrix tree and compute profile optimization ordering for diagonal subblocks

	const char *funcname = "PartOrdMtr_04";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartOrdMtr (_tree, ndweights, adjweights, _order);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS including Schur complement
//========================================================================================
void CSMatrix::PartOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_ndweights, int *_adjweights, int *_order, const CSlvParam &_param) const { // Partition matrix by METIS including Schur complement

	const char *funcname = "PartOrdMtrSchur_00";

// Compute initial ordering

	PartOrdMtr (_tree, _ndweights, _adjweights, _order);

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Reorder

	CSMatrix ao;

	ao = as.OrdMtr(_order);

	CSMatrix mtrdummy;

	as = mtrdummy;

// Create extended tree

	_treeext = _tree.ExtendTree ();

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;
	int *order;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);
	order = new int [nsupr];
	if (!order) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupr;i++) order[_order[i]] = i;

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichild, isupbeg, isupend, isupbegtot, isupendtot, nloc, nzloc;
	int j, jj;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;
		isupbegtot = _tree.nodes[inodecurr].indbegtot;
		isupendtot = _tree.nodes[inodecurr].indendtot;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

		if (nchildscurr != 1) {

// Compute local blocks partitioning

			int *blksloc;

			blksloc = new int [nchildscurr+2];
			if (!blksloc) MemoryFail (funcname);

			blksloc[0] = 0;
			blksloc[nchildscurr+1] = isupendtot+1-isupbegtot;

			for (ichild=0;ichild<nchildscurr;ichild++) {
				int inode2 = _tree.nodes[inodecurr].childs[ichild];
				blksloc[ichild+1] = _tree.nodes[inode2].indendtot+1-isupbegtot;
			};

// Take submatrix corresponding to the whole set of indices

			nloc = isupendtot-isupbegtot+1;
			nzloc = 0;

			for (i=isupbegtot;i<=isupendtot;i++) {
				for (j=ao.ia[i];j<ao.ia[i+1];j++) {
					jj = ao.ja[j];
					if (jj >= isupbegtot && jj <= isupendtot) nzloc++;
				};
			};

			CSMatrix aloc (nloc,nzloc);

			nloc = 0;
			nzloc = 0;
			aloc.ia[0] = 0;
			for (i=isupbegtot;i<=isupendtot;i++) {
				for (j=ao.ia[i];j<ao.ia[i+1];j++) {
					jj = ao.ja[j];
					if (jj >= isupbegtot && jj <= isupendtot) {
						aloc.ja[nzloc] = jj-isupbegtot;
						nzloc++;
					};
				};
				aloc.ia[nloc+1] = nzloc;
				nloc++;
			};

			aloc.m = nloc;
			aloc.n = nloc;
			aloc.nsupr = nloc;
			aloc.nsupc = nloc;
			aloc.nzja = nzloc;
			aloc.nlist = nloc;

//			char strbuff[256];
//			sprintf(strbuff, "%s%s%d%s",_param.path,"Aloc_",inodecurr,".ps");
//			aloc.A2Ps (3,strbuff,0,&isupbeg);

// Compute the sparsity of the bordering

			int *imask, *listloc;
			int *ibord, *jbord;

			imask = new int [nloc+1];
			if (!imask) MemoryFail (funcname);
			listloc = new int [nloc+1];
			if (!listloc) MemoryFail (funcname);
			ibord = new int [nloc+1];
			if (!ibord) MemoryFail (funcname);

			int icycle = -1;
			int kii;
			for (kii=0;kii<nloc;kii++) imask[kii] = icycle;

			int ibegbrd=blksloc[nchildscurr];

			ibord[0] = 0;

			nzloc = 0;

			int i;
			for (i=0;i<blksloc[nchildscurr];i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					if (jj >= ibegbrd) {
						if (imask[jj] != icycle) {
							nzloc++;
							imask[jj] = icycle;
						};
					};
				};
				ibord[i+1] = nzloc;
			};

			jbord = new int [nzloc];
			if (!jbord) MemoryFail (funcname);

			nzloc = 0;

			for (i=0;i<blksloc[nchildscurr];i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					if (jj >= ibegbrd) {
						if (imask[jj] != icycle) {
							jbord[nzloc] = jj;
							nzloc++;
							imask[jj] = icycle;
						};
					};
				};
			};

			for (i=0;i<blksloc[nchildscurr];i++) {
				nzloc = ibord[i+1]-ibord[i];
				int ibs = ibord[i];
				qsort (jbord+ibs, nzloc, sizeof(int), compint);
			};

			for (kii=blksloc[nchildscurr];kii<nloc;kii++) ibord[kii+1] = ibord[kii];

// Compute transposed bordering matrix

			nzloc = ibord[nloc];

			int *ibordt, *jbordt;

			ibordt = new int [nloc+1];
			if (!ibordt) MemoryFail (funcname);
			jbordt = new int [nzloc];
			if (!jbordt) MemoryFail (funcname);

			TransposeMatrix (nloc, ibord, jbord, ibordt, jbordt);

// Recompute aloc without bordering

			nzloc = 0;

			for (i=0;i<blksloc[nchildscurr];i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					if (jj < ibegbrd) {
						aloc.ja[nzloc] = jj;
						nzloc++;
					};
				};
				ibord[i+1] = nzloc;
			};

			for (i=blksloc[nchildscurr];i<nloc;i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					aloc.ja[nzloc] = jj;
					nzloc++;
				};
				ibord[i+1] = nzloc;
			};

			for (i=0;i<=nloc;i++) aloc.ia[i] = ibord[i];

// Compute extended bordering via ncycle parameter

			int nnloc = nloc-ibegbrd;

			for (i=0;i<=nnloc;i++) ibordt[i] = ibordt[i+ibegbrd];

			int **jbordarr;

			jbordarr = new int * [nnloc];
			if (!jbordarr) MemoryFail (funcname);

			for (int icyc=0;icyc<_param.ncycle;icyc++) {

// Extend

				nzloc = 0;

				ibord[i] = 0;

				for (i=0;i<nnloc;i++) {

					icycle++;

					int nlistloc = 0;

					for (j=ibordt[i];j<ibordt[i+1];j++) {
						jj = jbordt[j];
						for (int k=aloc.ia[jj];k<aloc.ia[jj+1];k++) {
							int kk = aloc.ja[k];
							if (imask[kk] != icycle) {
								listloc[nlistloc] = kk;
								nlistloc++;
								imask[kk] = icycle;
							};
						};
					};

					qsort (listloc, nlistloc, sizeof(int), compint);

					int *jloc;

					jloc = new int [nlistloc];
					if (!jloc) MemoryFail (funcname);

					for (kii=0;kii<nlistloc;kii++) jloc[kii] = listloc[kii];

					jbordarr[i] = jloc;

					nzloc += nlistloc;
					ibord[i+1] = nzloc;

				};

// Store

				delete [] jbord;

				jbord = new int [nzloc];
				if (!jbord) MemoryFail (funcname);

				for (i=0;i<nnloc;i++) {
					int nlistloc = ibord[i+1]-ibord[i];
					for (kii=0;kii<nlistloc;kii++) jbord[ibord[i]+kii] = jbordarr[i][kii];
				};

// Free memory

				for (i=0;i<nnloc;i++) delete [] jbordarr[i];

// Reassign for the next cycle

				int *iptr;

				iptr = ibordt;
				ibordt = ibord;
				ibord = iptr;

				iptr = jbordt;
				jbordt = jbord;
				jbord = iptr;

			};

			for (i=nnloc;i>=0;i--) ibordt[i+ibegbrd] = ibordt[i];
			for (i=0;i<=blksloc[nchildscurr];i++) ibordt[i] = 0;

// Recompute ibord and jbord arrays

			nzloc = ibordt[nloc];

			delete [] jbord;

			jbord = new int [nzloc];
			if (!jbord) MemoryFail (funcname);

			TransposeMatrix (nloc, ibordt, jbordt, ibord, jbord);

// Compute the memory required to store Schur complement

			nzloc = 0;

			for (i=ibegbrd;i<nloc;i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					if (jj >= i) {
						imask[jj] = icycle;
						nzloc++;
					};
				};
				for (j=ibordt[i];j<ibordt[i+1];j++) {
					jj = jbordt[j];
					for (int k=ibord[jj];k<ibord[jj+1];k++) {
						int kk = jbord[k];
						if (kk >= i) {
							if (imask[kk] != icycle) {
								imask[kk] = icycle;
								nzloc++;
							};
						};
					};
				};
			};

// Create Schur complement submatrix

			nnloc = nloc-ibegbrd;

			CSMatrix shloc (nnloc,nzloc);

			nnloc = 0;
			nzloc = 0;
			shloc.ia[0] = 0;
			for (i=ibegbrd;i<nloc;i++) {
				icycle++;
				for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
					jj = aloc.ja[j];
					if (jj >= i) {
						imask[jj] = icycle;
						shloc.ja[nzloc] = jj-ibegbrd;
						nzloc++;
					};
				};
				for (j=ibordt[i];j<ibordt[i+1];j++) {
					jj = jbordt[j];
					for (int k=ibord[jj];k<ibord[jj+1];k++) {
						int kk = jbord[k];
						if (kk >= i) {
							if (imask[kk] != icycle) {
								imask[kk] = icycle;
								shloc.ja[nzloc] = kk-ibegbrd;
								nzloc++;
							};
						};
					};
				};
				shloc.list[nnloc] = nnloc;
				nnloc++;
				shloc.ia[nnloc] = nzloc;
			};

			shloc.m = nnloc;
			shloc.n = nnloc;
			shloc.nsupr = nnloc;
			shloc.nsupc = nnloc;
			shloc.nzja = nzloc;
			shloc.nlist = nnloc;

			for (i=0;i<nnloc;i++) {
				nzloc = shloc.ia[i+1]-shloc.ia[i];
				int ibs = shloc.ia[i];
				qsort (shloc.ja+ibs, nzloc, sizeof(int), compint);
			};

//			sprintf(strbuff, "%s%s%d%s",_param.path,"Shloc_",inodecurr,".ps");
//			shloc.A2Ps (3,strbuff,0,&isupbeg);

// Compute array weights

			float *weights;

			weights = new float [nchildscurr+2];
			if (!weights) MemoryFail (funcname);

			for (i=0;i<nchildscurr;i++) weights[i] = 0.0e0;

			for (i=0;i<nchildscurr;i++) {
				ichildnode = _tree.nodes[inodecurr].childs[i];
				int ichildcpu = _tree.nodes[ichildnode].nodecpu;
				weights[i] = (float) _tree.perf2cpu[ichildcpu];
			};

			float faux = 0.0e0;

			for (i=0;i<nchildscurr;i++) faux += weights[i];

			for (i=0;i<nchildscurr;i++) weights[i] /= faux;

// Compute partitioning of the Schur complement

			int *orderloc, *iorderloc;

			orderloc = new int [nloc];
			if (!orderloc) MemoryFail (funcname);
			iorderloc = new int [nloc];
			if (!iorderloc) MemoryFail (funcname);

			if (nnloc >= nchildscurr+1) {
				shloc.PartOrdMtr (nchildscurr, weights, blksloc, orderloc);
			} else {
				for (i=0;i<=nchildscurr;i++) blksloc[i] = 0;
				blksloc[nchildscurr+1] = nnloc;
				for (i=0;i<nnloc;i++) orderloc[i] = i;
			};

			for (i=0;i<nnloc;i++) iorderloc[orderloc[i]] = i;

			CSMatrix aoloc;

			aoloc = shloc.OrdMtrSymm (orderloc);

//			sprintf(strbuff, "%s%s%d%s",_param.path,"ShlocOrd_",inodecurr,".ps");
//			aoloc.A2Ps (3,strbuff,nchildscurr+1,blksloc);

			aoloc = mtrdummy;

			delete [] weights;

// Update the ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderloc[i-isupbeg];
				orderloc[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderloc[i-isupbeg];
			};

// Store computed indices in the extended tree

			for (ichild=0;ichild<nchildscurr;ichild++) {
				ichildnode = _treeext.nodes[inodecurr].childs[ichild];

				_treeext.nodes[ichildnode].indbeg = isupbeg+blksloc[ichild];
				_treeext.nodes[ichildnode].indbegtot = isupbeg;

				_treeext.nodes[ichildnode].indend = isupbeg+blksloc[ichild+1]-1;
				_treeext.nodes[ichildnode].indendtot = isupend;

			};

			_treeext.nodes[inodecurr].indbeg = isupbeg+blksloc[nchildscurr];
			_treeext.nodes[inodecurr].indbegtot = isupbeg;

			_treeext.nodes[inodecurr].indend = isupend;
			_treeext.nodes[inodecurr].indendtot = isupend;

			ichildnode = _treeext.nodes[inodecurr].childs[0];

			int fathernodeid2 = _treeext.nodes[ichildnode].fatherid2;

			_treeext.nodes[fathernodeid2].indbegtot = isupbeg;
			_treeext.nodes[fathernodeid2].indendtot = isupend;

// Free work data

			aloc = mtrdummy;
			shloc = mtrdummy;

			delete [] blksloc;
			delete [] imask;
			delete [] listloc;
			delete [] ibord;
			delete [] jbord;
			delete [] ibordt;
			delete [] jbordt;
			delete [] orderloc;
			delete [] iorderloc;
			delete [] jbordarr;

		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in PartOrdMtr routine";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

// Delete working arrays

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] order;

// Condense tree back

	_tree = _treeext.CollapseTree ();

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS including Schur complement
//========================================================================================
void CSMatrix::PartOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_order, const CSlvParam &_param) const { // Partition matrix by METIS including Schur complement

	const char *funcname = "PartOrdMtrSchur_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartOrdMtrSchur (_tree, _treeext, ndweights, adjweights, _order, _param);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition AtA matrix by METIS including Schur complement
//========================================================================================
void CSMatrix::PartAtAOrdMtrSchur (CTree &_tree, CTree &_treeext, int *_order, const CSlvParam &_param) const { // Partition AtA matrix by METIS including Schur complement

//	const char *funcname = "PartAtAOrdMtrSchur";

// Compute AtA matrix

	CSMatrix ata;

	ata = AtAMatrSpars ();

// Compute new partitioning and ordering

	ata.PartOrdMtrSchur (_tree, _treeext, _order, _param);

};

// Author: Kharchenko S.A.
// CSMatrix: Compute the independent subsets of the matrix nodes
//========================================================================================
void CSMatrix::IndependentSubsets (int &_nblks, int *&_blks, int *_order) const { // Compute the independent subsets of the matrix nodes

	const char *funcname = "IndependentSubsets";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Allocate work arrays

	int *listloc, *imask;

	listloc = new int [nsupr+1];
	if (!listloc) MemoryFail (funcname);
	imask = new int [nsupr+1];
	if (!imask) MemoryFail (funcname);

	int i;

// Init work arrays

	for (i=0;i<nsupr;i++) imask[i] = -1;

	int nlistloc = 1;
	listloc[0] = 0;

// Main cycle

	int iblk = 0;

	while (nlistloc != 0) {
		nlistloc = 1;
		i = 0;
		while (i<nlistloc) {
			int isup = listloc[i];
			if (i==0) imask[isup] = iblk;
			if (imask[isup] != iblk) {
				throw " IndependentSubsets: wrong imask value";
			};
			for (int j=as.ia[isup];j<as.ia[isup+1];j++) {
				int jsup = as.ja[j];
				if (imask[jsup] != -1 && imask[jsup] != iblk) {
					throw " IndependentSubsets: wrong imask value for j index";
				};
				if (imask[jsup] == -1) {
					imask[jsup] = iblk;
					listloc[nlistloc] = jsup;
					nlistloc++;
				};
			};
			i++;
		};
		nlistloc = 0;
		for (i=0;i<nsupr;i++) {
			if (imask[i] == -1) {
				nlistloc = 1;
				listloc[0] = i;
				break;
			};
		};
		iblk++;
	};

// Return ordering and blocks partitioning

	_nblks = iblk;

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	for (i=0;i<=_nblks;i++) _blks[i] = 0;

	for (i=0;i<nsupr;i++) {
		iblk = imask[i];
		_blks[iblk+1]++;
	};

	for (i=0;i<_nblks;i++) _blks[i+1] = _blks[i]+_blks[i+1];

	for (i=0;i<_nblks;i++) listloc[i] = _blks[i];

	for (i=0;i<nsupr;i++) {
		iblk = imask[i];
		int k = listloc[iblk];
		_order[i] = k;
		listloc[iblk]++;
	};

// Delete work arrays

	delete [] listloc;
	delete [] imask;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
//========================================================================================
void CSMatrix::PartitionMtr (int _nparts, // Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
								int *_ndweights, int *_adjweights, int *_blks, int *_partition, int *_order) const {

	const char *funcname = "PartitionMtr_00";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local adjacency weights

	int *adjwloc;

	adjwloc = new int [as.nzja];
	if (!adjwloc) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<as.nzja;i++) adjwloc[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwloc[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwloc[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwloc[j] = _adjweights[k];
					};
				};
			};
			if (adjwloc[j] < 0) throw " PartitionMtr_00: index was not found or incorrect weights array";
		};
	};

// Call METIS 

	int options[5], volume, i1;

	options[0] = 0;

	int nloc = nsupr;

	int *pn = &nloc;
	int *pnp = &_nparts;
	int *pvol = &volume;
	int *pzero = &i;
	int *pone = &i1;

	i=0;
//	i1=3;
	i1=2;

//	METIS_PartGraphKway (pn, as.ia, as.ja, _ndweights, adjwloc, pone, pzero, pnp, options, pvol, _partition);
	METIS_PartGraphKway (pn, as.ia, as.ja, _ndweights, NULL, pone, pzero, pnp, options, pvol, _partition);

// Search independent nodes for each processor

	int *pblks;

	pblks = new int [_nparts];
	if (!pblks) MemoryFail (funcname);

	int icpu;

// Compute new ordering

// Count the number of elements in the independent blocks

	for (i=0;i<=_nparts;i++) _blks[i] = 0;

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		_blks[icpu+1]++;
	};

	for (i=0;i<_nparts;i++) _blks[i+1] += _blks[i];

// Set ordering array

	for (i=0;i<nsupr;i++) _order[i] = -1;
	for (i=0;i<_nparts;i++) pblks[i] = _blks[i];

	for (i=0;i<nsupr;i++) {
		icpu = _partition[i];
		_order[i] = pblks[icpu];
		pblks[icpu]++;
	};

	delete [] adjwloc;
	delete [] pblks;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
//========================================================================================
void CSMatrix::PartitionMtr (int _nparts, // Partition matrix by METIS without bordering splitting with weights of each node/edge and compute its ordering
								int *_blks, int *_partition, int *_order) const {

	const char *funcname = "PartitionMtr_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartitionMtr (_nparts, ndweights, adjweights, _blks, _partition, _order);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Reorder matrix and weights data in place
//========================================================================================
void CSMatrix::ReordMtrWeights (int *_ndweights, int *_adjweights, int *_order) { // Reorder matrix and weights data in place

	const char *funcname = "ReordMtrWeights";

	int *ialoc, *jaloc, *iorder;
	int *ndwloc, *adjwloc;

	ialoc = new int [nsupr+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);

	ndwloc = new int [nsupr];
	if (!ndwloc) MemoryFail (funcname);
	adjwloc = new int [nzja];
	if (!adjwloc) MemoryFail (funcname);

	int nzloc = 0;

	int i, j, jj;

	for (i=0;i<nsupr; i++) iorder[_order[i]] = i;

	ialoc[0] = 0;

	for (i=0;i<nsupr;i++) {
		int nzloc0 = nzloc;
		int iold = iorder[i];
		for (j=ia[iold];j<ia[iold+1];j++) {
			jj = ja[j];
			jaloc[nzloc] = _order[jj];
			adjwloc[nzloc] = _adjweights[j];
			nzloc++;
		};
//		qsort (jaloc+nzloc0, nzloc-nzloc0, sizeof(int), compint);
		SortElm (nzloc-nzloc0, jaloc+nzloc0, adjwloc+nzloc0);
		ialoc[i+1] = nzloc;
		ndwloc[i] = _ndweights[iold];
	};

	for (i=0;i<=nsupr; i++) ia[i] = ialoc[i];
	for (i=0;i<nzja;i++) ja[i] = jaloc[i];
	for (i=0;i<nsupr;  i++) _ndweights[i] = ndwloc[i];
	for (i=0;i<nzja;i++) _adjweights[i] = adjwloc[i];

	delete [] iorder;
	delete [] ialoc;
	delete [] jaloc;
	delete [] ndwloc;
	delete [] adjwloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix
//========================================================================================
void CSMatrix::Submatrix (int _ilistbeg, int _ilistend, int _isupbeg, int _isupend, // Take current submatrix
							CSMatrix &_asub) const {

//	const char *funcname = "Submatrix";

// Create local submatrix

	int nloc = _isupend-_isupbeg+1;
	int nzloc = 0;

	int i, j, jj;

	for (i=_ilistbeg;i<=_ilistend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) nzloc++;
		};
	};

	CSMatrix aloc (nloc,nzloc);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=_ilistbeg;i<=_ilistend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) {
				aloc.ja[nzloc] = jj-_isupbeg;
				nzloc++;
			};
		};
		aloc.ia[nloc+1] = nzloc;
		nloc++;
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix
//========================================================================================
void CSMatrix::Submatrix (int _isupbeg, int _isupend, // Take current submatrix
							CSMatrix &_asub) const {

//	const char *funcname = "Submatrix";

// Create local submatrix

	int nloc = _isupend-_isupbeg+1;
	int nzloc = 0;

	int i, j, jj;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) nzloc++;
		};
	};

	CSMatrix aloc (nloc,nzloc);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) {
				aloc.ja[nzloc] = jj-_isupbeg;
				nzloc++;
			};
		};
		aloc.ia[nloc+1] = nzloc;
		nloc++;
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix according to the list
//========================================================================================
void CSMatrix::SubmatrixList (int _nlist, int *_list, // Take current submatrix according to the list
										int &_icycle, int *_imask, int *_irow2ind,
										CSMatrix &_asub) const {

//	const char *funcname = "SubmatrixList";

// Mark the list

	int i, j;

	_icycle++;

	for (i=0;i<_nlist;i++) {
		j = _list[i];
		_imask[j] = _icycle;
		_irow2ind[j] = i;
	};

// Create local submatrix

	int nloc = _nlist;
	int nzloc = 0;

	int ii, jj;

	for (i=0;i<_nlist;i++) {
		ii = _list[i];
		for (j=ia[ii];j<ia[ii+1];j++) {
			jj = ja[j];
			if (_imask[jj] == _icycle) nzloc++;
		};
	};

	CSMatrix aloc (nloc,nzloc);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=0;i<_nlist;i++) {
		ii = _list[i];
		for (j=ia[ii];j<ia[ii+1];j++) {
			jj = ja[j];
			if (_imask[jj] == _icycle) {
				aloc.ja[nzloc] = _irow2ind[jj];
				nzloc++;
			};
		};
		qsort (aloc.ja+aloc.ia[nloc],nzloc-aloc.ia[nloc],sizeof(int),compint);
		aloc.list[nloc] = i;
		aloc.ia[nloc+1] = nzloc;
		nloc++;
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix by rows
//========================================================================================
void CSMatrix::SubmatrixByRows (int _isupbeg, int _isupend, // Take current submatrix by rows
							CSMatrix &_asub) const {

//	const char *funcname = "SubmatrixByRows";

// Create local submatrix

	int nloc = _isupend-_isupbeg+1;
	int nzloc = 0;

	int i, j, jj;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			nzloc++;
		};
	};

	CSMatrix aloc (nloc,nzloc);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
//			aloc.ja[nzloc] = jj-_isupbeg;
			aloc.ja[nzloc] = jj;
			nzloc++;
		};
		aloc.list[nloc] = i;
		aloc.ia[nloc+1] = nzloc;
		nloc++;
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix
//========================================================================================
void CSMatrix::Submatrix2Index (int _iblkbeg, int _iblkend, // Take current submatrix
											CSMatrix &_asub) const {

//	const char *funcname = "Submatrix";

// Create local submatrix

	int nloc = 0;
	int nzloc = 0;

	int i, j, jj, iblk, inode, jblk;

	for (i=0;i<nlist;i++) {
		iblk = list2[i];
		if (iblk >= _iblkbeg && iblk <= _iblkend) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja2[j];
				if (jj >= _iblkbeg && jj <= _iblkend) nzloc++;
			};
			nloc++;
		};
	};

	CSMatrix aloc (nloc,nloc,nzloc,nzloc);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=0;i<nlist;i++) {
		inode = list[i];
		iblk = list2[i];
		if (iblk >= _iblkbeg && iblk <= _iblkend) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jblk = ja2[j];
				if (jblk >= _iblkbeg && jblk <= _iblkend) {
					aloc.ja[nzloc] = jj;
					aloc.ja2[nzloc] = jblk;
					nzloc++;
				};
			};
			aloc.list[nloc] = inode;
			aloc.list2[nloc] = iblk;
			aloc.ia[nloc+1] = nzloc;
			nloc++;
		};
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nzja2 = nzloc;
	aloc.nlist = nloc;
	aloc.nlist2 = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Take current submatrix with weights
//========================================================================================
void CSMatrix::SubmatrixWeights (int _isupbeg, int _isupend, // Take current submatrix with weights
									int *_ndweights, int *_adjweights, 
									CSMatrix &_asub, int *&_ndweightssub, int *&_adjweightssub) const {

	const char *funcname = "SubmatrixWeights";

// Create local submatrix

	int nloc = _isupend-_isupbeg+1;
	int nzloc = 0;

	int i, j, jj;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) nzloc++;
		};
	};

	CSMatrix aloc (nloc,nzloc);

	_ndweightssub = new int [nloc];
	if (!_ndweightssub) MemoryFail (funcname);
	_adjweightssub = new int [nzloc];
	if (!_adjweightssub) MemoryFail (funcname);

	nloc = 0;
	nzloc = 0;
	aloc.ia[0] = 0;

	for (i=_isupbeg;i<=_isupend;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= _isupbeg && jj <= _isupend) {
				aloc.ja[nzloc] = jj-_isupbeg;
				_adjweightssub[nzloc] = _adjweights[j];
				nzloc++;
			};
		};
		aloc.ia[nloc+1] = nzloc;
		_ndweightssub[nloc] = _ndweights[i];
		nloc++;
	};
	aloc.m = nloc;
	aloc.n = nloc;
	aloc.nsupr = nloc;
	aloc.nsupc = nloc;
	aloc.nzja = nzloc;
	aloc.nlist = nloc;

	_asub = aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute the sparsity of the AtA matrix in prescribed rows and the resulting weights assuming that matrix is symmetric
//========================================================================================
void CSMatrix::A2MatrWeights (int _nrows, int *_listrows, // Compute the sparsity of the AtA matrix in prescribed rows and the resulting weights assuming that matrix is symmetric
									int *_ndweights, int *_adjweights, 
									CSMatrix &_asub, int *&_ndweightssub, int *&_adjweightssub) const {

	const char *funcname = "A2MatrWeights";

// Modify diagonal values of adjweights

	int i, j, jj, k, kk, ilist, isup;

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj == i) _adjweights[j] = _ndweights[i];
		};
	};

// Compute the upper estimate of the number of elements in the reduced AtA matrix

	int *imask, *imaskrow;
	double *fmask;
	int *lstloc;

	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	fmask = new double [nsupr];
	if (!fmask) MemoryFail (funcname);
	imaskrow = new int [nsupr];
	if (!imaskrow) MemoryFail (funcname);
	lstloc = new int [nsupr];
	if (!lstloc) MemoryFail (funcname);

	int icycle = -1;
	for (i=0;i<nsupr;i++) {
		imask[i] = icycle;
		fmask[i] = 0;
		imaskrow[i] = 0;
	};

	for (ilist=0;ilist<_nrows;ilist++) {
		isup = _listrows[ilist];
		imaskrow[isup] = 1;
	};

	int nz = 0;

	for (ilist=0;ilist<_nrows;ilist++) {
		isup = _listrows[ilist];
		icycle++;
		for (j=ia[isup];j<ia[isup+1];j++) {
			jj = ja[j];
			for (k=ia[jj];k<ia[jj+1];k++) {
				kk = ja[k];
				if (imask[kk] != icycle && imaskrow[kk]) {
					nz++;
					imask[kk] = icycle;
				};
			};
		};
	};

	nz += nsupr;

// Allocate local arrays

	int *ialoc;
	int *jaloc;

	_ndweightssub = new int [nsupr];
	if (!_ndweightssub) MemoryFail (funcname);
	_adjweightssub = new int [nz];
	if (!_adjweightssub) MemoryFail (funcname);

	ialoc = new int [nsupr+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nz];
	if (!jaloc) MemoryFail (funcname);

// Fill arrays

	int ijmax;

	ialoc[0] = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		icycle++;
		int nzloc = 0;
		if (imaskrow[isup]) {
			for (j=ia[isup];j<ia[isup+1];j++) {
				jj = ja[j];
				for (k=ia[jj];k<ia[jj+1];k++) {
					kk = ja[k];
					if (imask[kk] != icycle && imaskrow[kk]) {
						imask[kk] = icycle;
						fmask[kk] = 0.0e0;
						lstloc[nzloc] = kk;
						nzloc++;
					};
					if (imaskrow[kk]) {
						fmask[kk] += (double)_adjweights[j] * (double)_adjweights[k];
					};
				};
			};
			if (nzloc != 0) qsort (lstloc, nzloc, sizeof(int), compint);
			for (j=0;j<nzloc;j++) jaloc[nz+j] = lstloc[j];
			for (j=0;j<nzloc;j++) {
				jj = lstloc[j];
				ijmax = _ndweights[isup] * _ndweights[jj];
				if (fmask[jj] < ijmax) {
					_adjweightssub[nz+j] = (int)fmask[lstloc[j]];
				} else {
					_adjweightssub[nz+j] = ijmax;
				};
			};
		} else {
			nzloc = 1;
			jaloc[nz] = isup;
			_adjweightssub[nz] = _ndweights[isup];
		};
		nz += nzloc;
		ialoc[isup+1] = nz;
		_ndweightssub[isup] = _ndweights[isup];
	};

// Create ata structure

	CSMatrix ata (nsupc,nz);

	for (i=0;i<nsupc;i++) ata.list[i] = i;
	for (i=0;i<=nsupc;i++) ata.ia[i] = ialoc[i];
	for (i=0;i<nz;i++) ata.ja[i] = jaloc[i];

	ata.m = nsupc;
	ata.n = nsupc;
	ata.nsupc = nsupc;
	ata.nsupr = nsupc;
	ata.nlist = nsupc;
	ata.nzja = nz;

// Free work arrays

	delete [] imask;
	delete [] fmask;
	delete [] imaskrow;
	delete [] lstloc;
	delete [] ialoc;
	delete [] jaloc;

	_asub = ata;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrForSortedTree (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtrForSortedTree_02";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local nodes/adjacency weights

	int *order;
	int *iorder;
	int *imask;
	int *ndwarr;
	int *adjwarr;

	order = new int [nsupr];
	if (!order) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);
	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	ndwarr = new int [nsupr];
	if (!ndwarr) MemoryFail (funcname);
	adjwarr = new int [as.nzja];
	if (!adjwarr) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<nsupr;i++) ndwarr[i] = _ndweights[i];

	for (i=0;i<as.nzja;i++) adjwarr[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwarr[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwarr[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwarr[j] = _adjweights[k];
					};
				};
			};
			if (adjwarr[j] < 0) throw " PartitionOrdMtr_02: index was not found or incorrect weights array";
		};
	};

// Reorder matrix and weights data in place

	as.ReordMtrWeights (ndwarr, adjwarr, _order);

// Create the list of bordering supernodes

	int nrows = 0;
	int *listrows;

	listrows = new int [nsupr];
	if (!listrows) MemoryFail (funcname);

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int isupbeg, isupend, nloc, iblk, inode;

	for (inode=0;inode<_tree.nnodes;inode++) {

		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nchildscurr = _tree.nodes[inode].nchilds;

		if (nchildscurr > 1) {
			for (i=isupbeg;i<=isupend;i++) {
				listrows[nrows] = i;
				nrows++;
			};
		};

	};

// Add into the list all adjacent rows/columns

	int *imaskloc;

	imaskloc = new int [nsupr];
	if (!imaskloc) MemoryFail (funcname);

	int icycleloc = -1;

	for (i=0;i<nsupr;i++) {
		imaskloc[i] = icycleloc;
	};

	icycleloc++;

	for (i=0;i<nrows;i++) {
		j = listrows[i];
		imaskloc[j] = icycleloc;
	};

	int nrowsini = nrows;

	int irow, ii;

	for (i=0;i<nrowsini;i++) {
		irow = listrows[i];
		for (ii=as.ia[irow];ii<as.ia[irow+1];ii++) {
			jj = as.ja[ii];
			if (imaskloc[jj] != icycleloc) {
				listrows[nrows] = jj;
				nrows++;
				imaskloc[jj] = icycleloc;
			};
		};
	};

	delete [] imaskloc;

	qsort (listrows, nrows, sizeof(int), compint);

// Compute second degree bordering extended structure with weights

	CSMatrix asext;
	CSMatrix asext1;

	int *ndwarrext;
	int *adjwarrext;
	int *ndwarrext1;
	int *adjwarrext1;

	as.A2MatrWeights (nrows, listrows, ndwarr, adjwarr, 
						asext1, ndwarrext1, adjwarrext1);

	asext1.A2MatrWeights (nrows, listrows, ndwarrext1, adjwarrext1, 
						asext, ndwarrext, adjwarrext);

	delete [] listrows;
	delete [] ndwarrext1;
	delete [] adjwarrext1;

// Allocate work arrays

	int nlev = _tree.nlev;
	int nnodes = _tree.nnodes;

	int *inodelv, *ichildlv, *nchildslv;
	int *pnblks;
	int **pblks;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);
	pnblks = new int [nnodes];
	if (!pnblks) MemoryFail (funcname);
	pblks = new int* [nnodes];
	if (!pblks) MemoryFail (funcname);

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	CSMatrix mtrdummy;

	int *ndwloc, *adjwloc;
	int *ndwlocblk, *adjwlocblk;

	for (i=0;i<nsupr;i++) order[_order[i]] = i;
	for (i=0;i<nsupr;i++) imask[i] = -1;

	int icycle = 0;

//	ofstream fout ("ChkBlkLoc.dat");

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

//		fout << " Inodecurr = " << inodecurr << endl;
//		fout << " Isupbeg = " << isupbeg << endl;
//		fout << " Isupend = " << isupend << endl;
//		fout << " Nchildscurr = " << nchildscurr << endl;

// Take current submatrix if necessary

		if (true) {
//		if (nchildscurr == 1) {

// Take current submatrix

			nloc = isupend-isupbeg+1;

			CSMatrix aloc;

			if (nchildscurr == 1) {

				as.SubmatrixWeights (isupbeg, isupend, ndwarr, adjwarr, 
										aloc, ndwloc, adjwloc);

			} else {

				asext.SubmatrixWeights (isupbeg, isupend, ndwarrext, adjwarrext,
										aloc, ndwloc, adjwloc);

			};

//			char buff[256];
//			sprintf(buff, "%s%i%s","IniMtr_",inodecurr,".ps");
//			aloc.A2Ps (buff,0,&nloc);

// Compute the independent blocks

			int *blksindep, *orderindep, *iorderindep;

			orderindep = new int [nloc];
			if (!orderindep) MemoryFail (funcname);
			iorderindep = new int [nloc];
			if (!iorderindep) MemoryFail (funcname);

			int nblksindep;

			if (nloc > 0) {
				aloc.IndependentSubsets (nblksindep,blksindep,orderindep);
			} else {
				nblksindep = 0;
				blksindep = new int [nblksindep+1];
				if (!blksindep) MemoryFail (funcname);
				blksindep[0] = 0;
			};

			for (i=0;i<nloc;i++) iorderindep[orderindep[i]] = i;

//			OutArr (fout," Blksindep = ",nblksindep+1,blksindep);

// Update the global ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderindep[i-isupbeg];
//					orderloc[iold] = _order[i];
				orderindep[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderindep[i-isupbeg];
			};

			for (i=0;i<nloc;i++) orderindep[iorderindep[i]] = i;

// Reorder the matrix and weights arrays locally

			aloc.ReordMtrWeights (ndwloc, adjwloc, orderindep);

//			sprintf(buff, "%s%i%s","OrdMtr_",inodecurr,".ps");
//			aloc.A2Ps (buff,0,blksindep);
//			sprintf(buff, "%s%i%s","OrdMtrBlk_",inodecurr,".ps");
//			aloc.A2Ps (buff,nblksindep,blksindep);

// Allocate work arrays to store independent data

			int *pnblksindep;
			int **pblksindep;

			pnblksindep = new int [nblksindep];
			if (!pnblksindep) MemoryFail (funcname);
			pblksindep = new int* [nblksindep];
			if (!pblksindep) MemoryFail (funcname);

// Cycle over the independent subblocks

			for (int iblkindep=0;iblkindep<nblksindep;iblkindep++) {

				int isupbegloc = blksindep[iblkindep];
				int isupendloc = blksindep[iblkindep+1]-1;

				int isupbegglob = isupbeg+isupbegloc;
				int isupendglob = isupbeg+isupendloc;

				int nlocblk = isupendloc-isupbegloc+1;

				CSMatrix alocblk;

				aloc.SubmatrixWeights (isupbegloc, isupendloc, ndwarr, adjwarr, 
										alocblk, ndwlocblk, adjwlocblk);

// Determine the number of partitions in the current block

				int wtotal = 0;

				for (i=0;i<nlocblk;i++) wtotal += ndwlocblk[i];

				int npartsloc;

				if (nchildscurr == 1) {
					npartsloc = (wtotal+_weightcellvolume-1) / _weightcellvolume;
				} else {
					npartsloc = (wtotal+_weightcellboundary-1) / _weightcellboundary;
				};

				if (npartsloc < 1) npartsloc = 1;

// Compute new ordering and partitioning

				if (npartsloc > 1) {

					int *blksloc, *partitionloc, *orderloc, *iorderloc;
					int *blks1loc;

					blksloc = new int [npartsloc+2];
					if (!blksloc) MemoryFail (funcname);
					blks1loc = new int [npartsloc+2];
					if (!blks1loc) MemoryFail (funcname);
					partitionloc = new int [nlocblk];
					if (!partitionloc) MemoryFail (funcname);
					orderloc = new int [nlocblk];
					if (!orderloc) MemoryFail (funcname);
					iorderloc = new int [nlocblk];
					if (!iorderloc) MemoryFail (funcname);

					alocblk.PartitionMtr (npartsloc, ndwlocblk, adjwlocblk, blksloc, partitionloc, orderloc);
//					OutArr (fout," Partition after PartitionMtr = ",npartsloc+1,blksloc);

					for (i=0;i<nlocblk;i++) iorderloc[orderloc[i]] = i;

// Free submatrix data

					alocblk = mtrdummy;

// Update the ordering

					for (i=isupbegglob;i<=isupendglob;i++) {
						int iold = iorderloc[i-isupbegglob];
//						orderloc[iold] = _order[i];
						orderloc[i-isupbegglob] = order[iold+isupbegglob];
					};

					for (i=isupbegglob;i<=isupendglob;i++) {
						order[i] = orderloc[i-isupbegglob];
					};

// Register partitioning data

					for (iblk=0;iblk<npartsloc;iblk++) {
						icycle++;
//						fout << " Imask arr ini = " << isupbegglob+blksloc[iblk] << " fin = " << isupbegglob+blksloc[iblk+1]-1 << endl;
						for (i=blksloc[iblk];i<blksloc[iblk+1];i++) {
							imask[isupbegglob+i] = icycle;
						};
					};

// Free local arrays

					delete [] blksloc;
					delete [] partitionloc;
					delete [] orderloc;
					delete [] iorderloc;

					for (i=0;i<npartsloc+1;i++) blks1loc[i] = blksloc[i];

					pnblksindep[iblkindep] = npartsloc;
					pblksindep[iblkindep] = blks1loc;

				} else {

					icycle++;
//					fout << " Imask arr ini = " << isupbegglob << " fin = " << isupendglob << endl;
					for (i=isupbegglob;i<=isupendglob;i++) {
						imask[i] = icycle;
					};

					int *blks1loc;

					blks1loc = new int [2];
					if (!blks1loc) MemoryFail (funcname);

					blks1loc[0] = 0;
					blks1loc[1] = nlocblk;

					pnblksindep[iblkindep] = 1;
					pblksindep[iblkindep] = blks1loc;

				};

				delete [] ndwlocblk;
				delete [] adjwlocblk;

			};

// Free work arrays

			delete [] ndwloc;
			delete [] adjwloc;

			delete [] blksindep;
			delete [] orderindep;
			delete [] iorderindep;

// Combine partitionings of the independent subblocks into one and free work arrays

			int nblksnode = 0;
			for (i=0;i<nblksindep;i++) nblksnode += pnblksindep[i];

			pnblks[inodecurr] = nblksnode;

			int *blksnode;

			blksnode = new int [nblksnode+1];
			if (!blksnode) MemoryFail (funcname);

			blksnode[0] = 0;

			int ip = 0;
			for (i=0;i<nblksindep;i++) {
				for (j=0;j<pnblksindep[i];j++) {
					blksnode[ip+1] = blksnode[ip]+(pblksindep[i][j+1]-pblksindep[i][j]);
					ip++;
				};
			};

			pblks[inodecurr] = blksnode;

			delete [] pnblksindep;
			for (i=0;i<nblksindep;i++) {
				delete [] pblksindep[i];
			};
			delete [] pblksindep;

		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in PartOrdMtr routine";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

// Check blocks partitioning

	for (i=0;i<nsupr;i++) {
		if (imask[i] == -1) {
			throw " PartitionOrdMtr_02: not all nodes were processed";
		};
	};

// Transform local partitioning data into the global partitioning data

	_nblks = 0;
	for (i=0;i<nnodes;i++) _nblks += pnblks[i];

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	_blks[0] = 0;

	int ip = 0;
	for (i=0;i<nnodes;i++) {
		_tree.nodes[i].indbeg = ip;
		_tree.nodes[i].indbegtot = ip;
		for (j=0;j<pnblks[i];j++) {
			_blks[ip+1] = _blks[ip]+(pblks[i][j+1]-pblks[i][j]);
			ip++;
		};
		_tree.nodes[i].indend = ip-1;
		_tree.nodes[i].indendtot = ip-1;
	};

	delete [] pnblks;
	for (i=0;i<nnodes;i++) {
		delete [] pblks[i];
	};
	delete [] pblks;

//	OutArr (fout," Final blocks partitioning = ",_nblks+1,_blks);

// Perform final modification of the extended tree data

	int nchildsloc, iblkbeg, iblkend;

	for (inode=0;inode<nnodes;inode++) {
		nchildsloc = _tree.nodes[inode].nchilds;
		if (nchildsloc > 1) {
			ifathernode = inode;
			ichildnode = _tree.nodes[inode].childs[0];
			iblkbeg = _tree.nodes[ichildnode].indbegtot;
			iblkend = _tree.nodes[ifathernode].indend;
			_tree.nodes[inode].indbegtot = iblkbeg;
			_tree.nodes[inode].indendtot = iblkend;
		};
	};

// Free work arrays

	delete [] order;
	delete [] iorder;
	delete [] imask;
	delete [] adjwarr;
	delete [] ndwarr;
	delete [] adjwarrext;
	delete [] ndwarrext;

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrForSortedTree (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtr_03";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartitionOrdMtrForSortedTree (_tree, ndweights, adjweights, _order,
						_weightcellvolume, _weightcellboundary, _nblks, _blks);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtr (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtr_02";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local nodes/adjacency weights

	int *order;
	int *iorder;
	int *imask;
	int *ndwarr;
	int *adjwarr;

	order = new int [nsupr];
	if (!order) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);
	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	ndwarr = new int [nsupr];
	if (!ndwarr) MemoryFail (funcname);
	adjwarr = new int [as.nzja];
	if (!adjwarr) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<nsupr;i++) ndwarr[i] = _ndweights[i];

	for (i=0;i<as.nzja;i++) adjwarr[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwarr[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwarr[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwarr[j] = _adjweights[k];
					};
				};
			};
			if (adjwarr[j] < 0) throw " PartitionOrdMtr_02: index was not found or incorrect weights array";
		};
	};

// Reorder matrix and weights data in place

	as.ReordMtrWeights (ndwarr, adjwarr, _order);

// Create the list of bordering supernodes

	int nrows = 0;
	int *listrows;

	listrows = new int [nsupr];
	if (!listrows) MemoryFail (funcname);

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int isupbeg, isupend, nloc, iblk, iblkloc, inode;

	for (inode=0;inode<_tree.nnodes;inode++) {

		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nchildscurr = _tree.nodes[inode].nchilds;

		if (nchildscurr > 1) {
			for (i=isupbeg;i<=isupend;i++) {
				listrows[nrows] = i;
				nrows++;
			};
		};

	};

// Add into the list all adjacent rows/columns

	int *imaskloc;

	imaskloc = new int [nsupr];
	if (!imaskloc) MemoryFail (funcname);

	int icycleloc = -1;

	for (i=0;i<nsupr;i++) {
		imaskloc[i] = icycleloc;
	};

	icycleloc++;

	for (i=0;i<nrows;i++) {
		j = listrows[i];
		imaskloc[j] = icycleloc;
	};

	int nrowsini = nrows;

	int irow, ii;

	for (i=0;i<nrowsini;i++) {
		irow = listrows[i];
		for (ii=as.ia[irow];ii<as.ia[irow+1];ii++) {
			jj = as.ja[ii];
			if (imaskloc[jj] != icycleloc) {
				listrows[nrows] = jj;
				nrows++;
				imaskloc[jj] = icycleloc;
			};
		};
	};

	delete [] imaskloc;

	qsort (listrows, nrows, sizeof(int), compint);

// Compute second degree bordering extended structure with weights

	CSMatrix asext;
	CSMatrix asext1;

	int *ndwarrext;
	int *adjwarrext;
	int *ndwarrext1;
	int *adjwarrext1;

	as.A2MatrWeights (nrows, listrows, ndwarr, adjwarr, 
						asext1, ndwarrext1, adjwarrext1);

	asext1.A2MatrWeights (nrows, listrows, ndwarrext1, adjwarrext1, 
						asext, ndwarrext, adjwarrext);

	delete [] listrows;
	delete [] ndwarrext1;
	delete [] adjwarrext1;

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	CSMatrix mtrdummy;

	int *ndwloc, *adjwloc;
	int *ndwlocblk, *adjwlocblk;

	for (i=0;i<nsupr;i++) order[_order[i]] = i;
	for (i=0;i<nsupr;i++) imask[i] = -1;

	int icycle = 0;

//	ofstream fout ("ChkBlkLoc.dat");

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

//		fout << " Inodecurr = " << inodecurr << endl;
//		fout << " Isupbeg = " << isupbeg << endl;
//		fout << " Isupend = " << isupend << endl;
//		fout << " Nchildscurr = " << nchildscurr << endl;

// Take current submatrix if necessary

		if (true) {
//		if (nchildscurr == 1) {

// Take current submatrix

			nloc = isupend-isupbeg+1;

			CSMatrix aloc;

			if (nchildscurr == 1) {

				as.SubmatrixWeights (isupbeg, isupend, ndwarr, adjwarr, 
										aloc, ndwloc, adjwloc);

			} else {

				asext.SubmatrixWeights (isupbeg, isupend, ndwarrext, adjwarrext,
										aloc, ndwloc, adjwloc);

			};

//			char buff[256];
//			sprintf(buff, "%s%i%s","IniMtr_",inodecurr,".ps");
//			aloc.A2Ps (buff,0,&nloc);

// Compute the independent blocks

			int *blksindep, *orderindep, *iorderindep;

			orderindep = new int [nloc];
			if (!orderindep) MemoryFail (funcname);
			iorderindep = new int [nloc];
			if (!iorderindep) MemoryFail (funcname);

			int nblksindep;

			aloc.IndependentSubsets (nblksindep,blksindep,orderindep);

			for (i=0;i<nloc;i++) iorderindep[orderindep[i]] = i;

//			OutArr (fout," Blksindep = ",nblksindep+1,blksindep);

// Update the global ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderindep[i-isupbeg];
//					orderloc[iold] = _order[i];
				orderindep[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderindep[i-isupbeg];
			};

			for (i=0;i<nloc;i++) orderindep[iorderindep[i]] = i;

// Reorder the matrix and weights arrays locally

			aloc.ReordMtrWeights (ndwloc, adjwloc, orderindep);

//			sprintf(buff, "%s%i%s","OrdMtr_",inodecurr,".ps");
//			aloc.A2Ps (buff,0,blksindep);
//			sprintf(buff, "%s%i%s","OrdMtrBlk_",inodecurr,".ps");
//			aloc.A2Ps (buff,nblksindep,blksindep);

// Cycle over the independent subblocks

			for (int iblkindep=0;iblkindep<nblksindep;iblkindep++) {

				int isupbegloc = blksindep[iblkindep];
				int isupendloc = blksindep[iblkindep+1]-1;

				int isupbegglob = isupbeg+isupbegloc;
				int isupendglob = isupbeg+isupendloc;

				int nlocblk = isupendloc-isupbegloc+1;

				CSMatrix alocblk;

				aloc.SubmatrixWeights (isupbegloc, isupendloc, ndwarr, adjwarr, 
										alocblk, ndwlocblk, adjwlocblk);

// Determine the number of partitions in the current block

				int wtotal = 0;

				for (i=0;i<nlocblk;i++) wtotal += ndwlocblk[i];

				int npartsloc;

				if (nchildscurr == 1) {
					npartsloc = (wtotal+_weightcellvolume-1) / _weightcellvolume;
				} else {
					npartsloc = (wtotal+_weightcellboundary-1) / _weightcellboundary;
				};

				if (npartsloc < 1) npartsloc = 1;

// Compute new ordering and partitioning

				if (npartsloc > 1) {

					int *blksloc, *partitionloc, *orderloc, *iorderloc;

					blksloc = new int [npartsloc+2];
					if (!blksloc) MemoryFail (funcname);
					partitionloc = new int [nlocblk];
					if (!partitionloc) MemoryFail (funcname);
					orderloc = new int [nlocblk];
					if (!orderloc) MemoryFail (funcname);
					iorderloc = new int [nlocblk];
					if (!iorderloc) MemoryFail (funcname);

					alocblk.PartitionMtr (npartsloc, ndwlocblk, adjwlocblk, blksloc, partitionloc, orderloc);
//					OutArr (fout," Partition after PartitionMtr = ",npartsloc+1,blksloc);

					for (i=0;i<nlocblk;i++) iorderloc[orderloc[i]] = i;

// Free submatrix data

					alocblk = mtrdummy;

// Update the ordering

					for (i=isupbegglob;i<=isupendglob;i++) {
						int iold = iorderloc[i-isupbegglob];
//						orderloc[iold] = _order[i];
						orderloc[i-isupbegglob] = order[iold+isupbegglob];
					};

					for (i=isupbegglob;i<=isupendglob;i++) {
						order[i] = orderloc[i-isupbegglob];
					};

// Register partitioning data

					for (iblk=0;iblk<npartsloc;iblk++) {
						icycle++;
//						fout << " Imask arr ini = " << isupbegglob+blksloc[iblk] << " fin = " << isupbegglob+blksloc[iblk+1]-1 << endl;
						for (i=blksloc[iblk];i<blksloc[iblk+1];i++) {
							imask[isupbegglob+i] = icycle;
						};
					};

// Free local arrays

					delete [] blksloc;
					delete [] partitionloc;
					delete [] orderloc;
					delete [] iorderloc;

				} else {

					icycle++;
//					fout << " Imask arr ini = " << isupbegglob << " fin = " << isupendglob << endl;
					for (i=isupbegglob;i<=isupendglob;i++) {
						imask[i] = icycle;
					};

				};

				delete [] ndwlocblk;
				delete [] adjwlocblk;

			};

// Free work arrays

			delete [] ndwloc;
			delete [] adjwloc;

			delete [] blksindep;
			delete [] orderindep;
			delete [] iorderindep;

		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in PartOrdMtr routine";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

// Check blocks partitioning

	for (i=0;i<nsupr;i++) {
		if (imask[i] == -1) {
			throw " PartitionOrdMtr_02: not all nodes were processed";
		};
	};

// Transform local partitioning data into the global partitioning data

	_blks = new int [icycle+2];
	if (!_blks) MemoryFail (funcname);

	_blks[0] = 0;

	iblk = 0;

	i = 0;
	iblkloc = imask[0];

	while (i<nsupr) {
		if (imask[i] != iblkloc) {
			_blks[iblk+1] = i;
			if (i != nsupr-1) iblk++;
			if (i<nsupr-1) {
				iblkloc = imask[i];
			};
		};
		i++;
		if (i == nsupr-1 && nsupr > 1) {
			_blks[iblk+1] = nsupr;
			iblk++;
		} else if (nsupr == 1) {
			_blks[iblk+1] = 1;
			iblk++;
		};
	};

	_nblks = iblk;

//	OutArr (fout," Final blocks partitioning = ",_nblks+1,_blks);

// Free work arrays

	delete [] order;
	delete [] iorder;
	delete [] imask;
	delete [] adjwarr;
	delete [] ndwarr;
	delete [] adjwarrext;
	delete [] ndwarrext;

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtr (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtr_03";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartitionOrdMtr (_tree, ndweights, adjweights, _order,
						_weightcellvolume, _weightcellboundary, _nblks, _blks);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrSchurForSortedTree (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtrSchurForSortedTree_00";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local nodes/adjacency weights

	int *order;
	int *iorder;
	int *imask;
	int *ndwarr;
	int *adjwarr;

	order = new int [nsupr];
	if (!order) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);
	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	ndwarr = new int [nsupr];
	if (!ndwarr) MemoryFail (funcname);
	adjwarr = new int [as.nzja];
	if (!adjwarr) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<nsupr;i++) ndwarr[i] = _ndweights[i];

	for (i=0;i<as.nzja;i++) adjwarr[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwarr[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwarr[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwarr[j] = _adjweights[k];
					};
				};
			};
			if (adjwarr[j] < 0) throw " PartitionOrdMtrSchur_00: index was not found or incorrect weights array";
		};
	};

// Reorder matrix and weights data in place

	as.ReordMtrWeights (ndwarr, adjwarr, _order);

// Create the list of bordering supernodes

	int nrows = 0;
	int *listrows;

	listrows = new int [nsupr];
	if (!listrows) MemoryFail (funcname);

	int isupbeg, isupend, nloc, iblk, inode, itype, nchildscurr;

	for (inode=0;inode<_tree.nnodes;inode++) {

		itype = _tree.nodes[inode].nodetype;
		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nchildscurr = _tree.nodes[inode].nchilds;

		if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {
			for (i=isupbeg;i<=isupend;i++) {
				listrows[nrows] = i;
				nrows++;
			};
		};

	};

// Add into the list all adjacent rows/columns

	int *imaskloc;

	imaskloc = new int [nsupr];
	if (!imaskloc) MemoryFail (funcname);

	int icycleloc = -1;

	for (i=0;i<nsupr;i++) {
		imaskloc[i] = icycleloc;
	};

	icycleloc++;

	for (i=0;i<nrows;i++) {
		j = listrows[i];
		imaskloc[j] = icycleloc;
	};

	int nrowsini = nrows;

	int irow, ii;

	for (i=0;i<nrowsini;i++) {
		irow = listrows[i];
		for (ii=as.ia[irow];ii<as.ia[irow+1];ii++) {
			jj = as.ja[ii];
			if (imaskloc[jj] != icycleloc) {
				listrows[nrows] = jj;
				nrows++;
				imaskloc[jj] = icycleloc;
			};
		};
	};

	delete [] imaskloc;

	qsort (listrows, nrows, sizeof(int), compint);

// Compute second degree bordering extended structure with weights

	CSMatrix asext;
	CSMatrix asext1;

	int *ndwarrext;
	int *adjwarrext;
	int *ndwarrext1;
	int *adjwarrext1;

	as.A2MatrWeights (nrows, listrows, ndwarr, adjwarr, 
						asext1, ndwarrext1, adjwarrext1);

	asext1.A2MatrWeights (nrows, listrows, ndwarrext1, adjwarrext1, 
						asext, ndwarrext, adjwarrext);

	delete [] listrows;
	delete [] ndwarrext1;
	delete [] adjwarrext1;

// Search the tree

	CSMatrix mtrdummy;

	int *ndwloc, *adjwloc;
	int *ndwlocblk, *adjwlocblk;

	int nnodes = _tree.nnodes;

	int *pnblks;
	int **pblks;

	pnblks = new int [nnodes];
	if (!pnblks) MemoryFail (funcname);
	pblks = new int* [nnodes];
	if (!pblks) MemoryFail (funcname);


	for (i=0;i<nsupr;i++) order[_order[i]] = i;
	for (i=0;i<nsupr;i++) imask[i] = -1;

	int icycle = 0;

	for (inode=0;inode<_tree.nnodes;inode++) {

		itype = _tree.nodes[inode].nodetype;
		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nloc = isupend-isupbeg+1;

		nchildscurr = _tree.nodes[inode].nchilds;

// Take current submatrix if necessary

		if (itype != exchangenode && nloc > 0) {

// Take current submatrix

			CSMatrix aloc;

			if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {

				asext.SubmatrixWeights (isupbeg, isupend, ndwarrext, adjwarrext,
										aloc, ndwloc, adjwloc);

			} else {

				as.SubmatrixWeights (isupbeg, isupend, ndwarr, adjwarr, 
										aloc, ndwloc, adjwloc);

			};

// Compute the independent blocks

			int *blksindep, *orderindep, *iorderindep;

			orderindep = new int [nloc];
			if (!orderindep) MemoryFail (funcname);
			iorderindep = new int [nloc];
			if (!iorderindep) MemoryFail (funcname);

			int nblksindep;

			aloc.IndependentSubsets (nblksindep,blksindep,orderindep);

			for (i=0;i<nloc;i++) iorderindep[orderindep[i]] = i;

// Update the global ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderindep[i-isupbeg];
				orderindep[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderindep[i-isupbeg];
			};

			for (i=0;i<nloc;i++) orderindep[iorderindep[i]] = i;

// Reorder the matrix and weights arrays locally

			aloc.ReordMtrWeights (ndwloc, adjwloc, orderindep);

// Allocate work arrays to store independent data

			int *pnblksindep;
			int **pblksindep;

			pnblksindep = new int [nblksindep];
			if (!pnblksindep) MemoryFail (funcname);
			pblksindep = new int* [nblksindep];
			if (!pblksindep) MemoryFail (funcname);

// Cycle over the independent subblocks

			for (int iblkindep=0;iblkindep<nblksindep;iblkindep++) {

				int isupbegloc = blksindep[iblkindep];
				int isupendloc = blksindep[iblkindep+1]-1;

				int isupbegglob = isupbeg+isupbegloc;
				int isupendglob = isupbeg+isupendloc;

				int nlocblk = isupendloc-isupbegloc+1;

				CSMatrix alocblk;

				aloc.SubmatrixWeights (isupbegloc, isupendloc, ndwloc, adjwloc, 
										alocblk, ndwlocblk, adjwlocblk);

// Determine the number of partitions in the current block

				int wtotal = 0;

				for (i=0;i<nlocblk;i++) wtotal += ndwlocblk[i];

				int npartsloc;

				if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {
					npartsloc = (wtotal+_weightcellboundary-1) / _weightcellboundary;
				} else {
					npartsloc = (wtotal+_weightcellvolume-1) / _weightcellvolume;
				};

				if (npartsloc < 1) npartsloc = 1;

// Compute new ordering and partitioning

				if (npartsloc > 1) {

					int *blksloc, *partitionloc, *orderloc, *iorderloc;
					int *blks1loc;

					blksloc = new int [npartsloc+2];
					if (!blksloc) MemoryFail (funcname);
					blks1loc = new int [npartsloc+2];
					if (!blks1loc) MemoryFail (funcname);
					partitionloc = new int [nlocblk];
					if (!partitionloc) MemoryFail (funcname);
					orderloc = new int [nlocblk];
					if (!orderloc) MemoryFail (funcname);
					iorderloc = new int [nlocblk];
					if (!iorderloc) MemoryFail (funcname);

					alocblk.PartitionMtr (npartsloc, ndwlocblk, adjwlocblk, blksloc, partitionloc, orderloc);

					for (i=0;i<nlocblk;i++) iorderloc[orderloc[i]] = i;

// Free submatrix data

					alocblk = mtrdummy;

// Update the ordering

					for (i=isupbegglob;i<=isupendglob;i++) {
						int iold = iorderloc[i-isupbegglob];
						orderloc[i-isupbegglob] = order[iold+isupbegglob];
					};

					for (i=isupbegglob;i<=isupendglob;i++) {
						order[i] = orderloc[i-isupbegglob];
					};

// Register partitioning data

					for (iblk=0;iblk<npartsloc;iblk++) {
						icycle++;
						for (i=blksloc[iblk];i<blksloc[iblk+1];i++) {
							imask[isupbegglob+i] = icycle;
						};
					};

					for (i=0;i<npartsloc+1;i++) blks1loc[i] = blksloc[i];

					pnblksindep[iblkindep] = npartsloc;
					pblksindep[iblkindep] = blks1loc;

// Free local arrays

					delete [] blksloc;
					delete [] partitionloc;
					delete [] orderloc;
					delete [] iorderloc;

				} else {

					icycle++;
					for (i=isupbegglob;i<=isupendglob;i++) {
						imask[i] = icycle;
					};

					int *blks1loc;

					blks1loc = new int [2];
					if (!blks1loc) MemoryFail (funcname);

					blks1loc[0] = 0;
					blks1loc[1] = nlocblk;

					pnblksindep[iblkindep] = 1;
					pblksindep[iblkindep] = blks1loc;

				};

				delete [] ndwlocblk;
				delete [] adjwlocblk;

			};

// Free work arrays

			delete [] ndwloc;
			delete [] adjwloc;

			delete [] blksindep;
			delete [] orderindep;
			delete [] iorderindep;

// Combine partitionings of the independent subblocks into one and free work arrays

			int nblksnode = 0;
			for (i=0;i<nblksindep;i++) nblksnode += pnblksindep[i];

			pnblks[inode] = nblksnode;

			int *blksnode;

			blksnode = new int [nblksnode+1];
			if (!blksnode) MemoryFail (funcname);

			blksnode[0] = 0;

			int ip = 0;
			for (i=0;i<nblksindep;i++) {
				for (j=0;j<pnblksindep[i];j++) {
					blksnode[ip+1] = blksnode[ip]+(pblksindep[i][j+1]-pblksindep[i][j]);
					ip++;
				};
			};

			pblks[inode] = blksnode;

			delete [] pnblksindep;
			for (i=0;i<nblksindep;i++) {
				delete [] pblksindep[i];
			};
			delete [] pblksindep;

		} else {

			int *blks1loc;

			blks1loc = new int [2];
			if (!blks1loc) MemoryFail (funcname);

			blks1loc[0] = 0;
			blks1loc[1] = nloc;

			pnblks[inode] = 0;
			pblks[inode] = blks1loc;

		};

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

// Check blocks partitioning

	for (i=0;i<nsupr;i++) {
		if (imask[i] == -1) {
			throw " PartitionOrdMtrSchur_00: not all nodes were processed";
		};
	};

// Transform local partitioning data into the global partitioning data

	_nblks = 0;
	for (i=0;i<nnodes;i++) _nblks += pnblks[i];

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	_blks[0] = 0;

	int ip = 0;
	for (i=0;i<nnodes;i++) {
		_tree.nodes[i].indbeg = ip;
		_tree.nodes[i].indbegtot = ip;
		for (j=0;j<pnblks[i];j++) {
			_blks[ip+1] = _blks[ip]+(pblks[i][j+1]-pblks[i][j]);
			ip++;
		};
		_tree.nodes[i].indend = ip-1;
		_tree.nodes[i].indendtot = ip-1;
	};

	delete [] pnblks;
	for (i=0;i<nnodes;i++) {
		delete [] pblks[i];
	};
	delete [] pblks;

// Perform final modification of the extended tree data

	int ifathernode, nchilds2loc, ichildnode, iblkbeg, iblkend;

	for (inode=0;inode<nnodes;inode++) {
		itype = _tree.nodes[inode].nodetype;
		if (itype == exchangenode) {
			ifathernode = _tree.nodes[inode].fatherid;
			nchilds2loc = _tree.nodes[inode].nchilds2;
			ichildnode = _tree.nodes[inode].childs2[0];
			iblkbeg = _tree.nodes[ichildnode].indbeg;
			iblkend = _tree.nodes[ifathernode].indend;
			_tree.nodes[inode].indbeg = -1;
			_tree.nodes[inode].indend = -1;
			_tree.nodes[inode].indbegtot = iblkbeg;
			_tree.nodes[inode].indendtot = iblkend;
			_tree.nodes[ifathernode].indbegtot = iblkbeg;
			_tree.nodes[ifathernode].indendtot = iblkend;
			for (i=0;i<nchilds2loc;i++) {
				ichildnode = _tree.nodes[inode].childs2[i];
				_tree.nodes[ichildnode].indbegtot = iblkbeg;
				_tree.nodes[ichildnode].indendtot = iblkend;
			};
		};
	};

// Free work arrays

	delete [] order;
	delete [] iorder;
	delete [] imask;
	delete [] adjwarr;
	delete [] ndwarr;
	delete [] adjwarrext;
	delete [] ndwarrext;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrSchurForSortedTree (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtrSchur_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartitionOrdMtrSchurForSortedTree (_tree, ndweights, adjweights, _order,
						_weightcellvolume, _weightcellboundary, _nblks, _blks);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrSchur (CTree &_tree, int *_ndweights, int *_adjweights, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtrSchur_00";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Assign local nodes/adjacency weights

	int *order;
	int *iorder;
	int *imask;
	int *ndwarr;
	int *adjwarr;

	order = new int [nsupr];
	if (!order) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);
	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	ndwarr = new int [nsupr];
	if (!ndwarr) MemoryFail (funcname);
	adjwarr = new int [as.nzja];
	if (!adjwarr) MemoryFail (funcname);

	int i, j, j1, j2, jj, jj1, jj2, iend1, iend2, k;

	for (i=0;i<nsupr;i++) ndwarr[i] = _ndweights[i];

	for (i=0;i<as.nzja;i++) adjwarr[i] = -1;

	for (i=0;i<nsupr;i++) {
		j1 = ia[i];
		j2 = as.ia[i];
		iend1 = ia[i+1];
		iend2 = as.ia[i+1];
		while (j1 < iend1 && j2 < iend2) {
			jj1 = ja[j1];
			jj2 = as.ja[j2];
			if (jj1 == jj2) {
				adjwarr[j2] = _adjweights[j1];
				j1++;
				j2++;
			} else if (jj1 < jj2) {
				j1++;
			} else if (jj1 > jj2) {
				j2++;
			};
		};
	};

	for (i=0;i<nsupr;i++) {
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (adjwarr[j] < 0) {
				for (k=ia[jj];k<ia[jj+1];k++) {
					if (ja[k] == i) {
						adjwarr[j] = _adjweights[k];
					};
				};
			};
			if (adjwarr[j] < 0) throw " PartitionOrdMtrSchur_00: index was not found or incorrect weights array";
		};
	};

// Reorder matrix and weights data in place

	as.ReordMtrWeights (ndwarr, adjwarr, _order);

// Create the list of bordering supernodes

	int nrows = 0;
	int *listrows;

	listrows = new int [nsupr];
	if (!listrows) MemoryFail (funcname);

	int isupbeg, isupend, nloc, iblk, iblkloc, inode, itype, nchildscurr;

	for (inode=0;inode<_tree.nnodes;inode++) {

		itype = _tree.nodes[inode].nodetype;
		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nchildscurr = _tree.nodes[inode].nchilds;

		if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {
			for (i=isupbeg;i<=isupend;i++) {
				listrows[nrows] = i;
				nrows++;
			};
		};

	};

// Add into the list all adjacent rows/columns

	int *imaskloc;

	imaskloc = new int [nsupr];
	if (!imaskloc) MemoryFail (funcname);

	int icycleloc = -1;

	for (i=0;i<nsupr;i++) {
		imaskloc[i] = icycleloc;
	};

	icycleloc++;

	for (i=0;i<nrows;i++) {
		j = listrows[i];
		imaskloc[j] = icycleloc;
	};

	int nrowsini = nrows;

	int irow, ii;

	for (i=0;i<nrowsini;i++) {
		irow = listrows[i];
		for (ii=as.ia[irow];ii<as.ia[irow+1];ii++) {
			jj = as.ja[ii];
			if (imaskloc[jj] != icycleloc) {
				listrows[nrows] = jj;
				nrows++;
				imaskloc[jj] = icycleloc;
			};
		};
	};

	delete [] imaskloc;

	qsort (listrows, nrows, sizeof(int), compint);

// Compute second degree bordering extended structure with weights

	CSMatrix asext;
	CSMatrix asext1;

	int *ndwarrext;
	int *adjwarrext;
	int *ndwarrext1;
	int *adjwarrext1;

	as.A2MatrWeights (nrows, listrows, ndwarr, adjwarr, 
						asext1, ndwarrext1, adjwarrext1);

	asext1.A2MatrWeights (nrows, listrows, ndwarrext1, adjwarrext1, 
						asext, ndwarrext, adjwarrext);

	delete [] listrows;
	delete [] ndwarrext1;
	delete [] adjwarrext1;

// Search the tree

	CSMatrix mtrdummy;

	int *ndwloc, *adjwloc;
	int *ndwlocblk, *adjwlocblk;

	for (i=0;i<nsupr;i++) order[_order[i]] = i;
	for (i=0;i<nsupr;i++) imask[i] = -1;

	int icycle = 0;

	for (inode=0;inode<_tree.nnodes;inode++) {

		itype = _tree.nodes[inode].nodetype;
		isupbeg = _tree.nodes[inode].indbeg;
		isupend = _tree.nodes[inode].indend;

		nloc = isupend-isupbeg+1;

		nchildscurr = _tree.nodes[inode].nchilds;

// Take current submatrix if necessary

		if (itype != exchangenode && nloc > 0) {

// Take current submatrix

			CSMatrix aloc;

			if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {

				asext.SubmatrixWeights (isupbeg, isupend, ndwarrext, adjwarrext,
										aloc, ndwloc, adjwloc);

			} else {

				as.SubmatrixWeights (isupbeg, isupend, ndwarr, adjwarr, 
										aloc, ndwloc, adjwloc);

			};

// Compute the independent blocks

			int *blksindep, *orderindep, *iorderindep;

			orderindep = new int [nloc];
			if (!orderindep) MemoryFail (funcname);
			iorderindep = new int [nloc];
			if (!iorderindep) MemoryFail (funcname);

			int nblksindep;

			aloc.IndependentSubsets (nblksindep,blksindep,orderindep);

			for (i=0;i<nloc;i++) iorderindep[orderindep[i]] = i;

// Update the global ordering

			for (i=isupbeg;i<=isupend;i++) {
				int iold = iorderindep[i-isupbeg];
				orderindep[i-isupbeg] = order[iold+isupbeg];
			};

			for (i=isupbeg;i<=isupend;i++) {
				order[i] = orderindep[i-isupbeg];
			};

			for (i=0;i<nloc;i++) orderindep[iorderindep[i]] = i;

// Reorder the matrix and weights arrays locally

			aloc.ReordMtrWeights (ndwloc, adjwloc, orderindep);

// Cycle over the independent subblocks

			for (int iblkindep=0;iblkindep<nblksindep;iblkindep++) {

				int isupbegloc = blksindep[iblkindep];
				int isupendloc = blksindep[iblkindep+1]-1;

				int isupbegglob = isupbeg+isupbegloc;
				int isupendglob = isupbeg+isupendloc;

				int nlocblk = isupendloc-isupbegloc+1;

				CSMatrix alocblk;

				aloc.SubmatrixWeights (isupbegloc, isupendloc, ndwarr, adjwarr, 
										alocblk, ndwlocblk, adjwlocblk);

// Determine the number of partitions in the current block

				int wtotal = 0;

				for (i=0;i<nlocblk;i++) wtotal += ndwlocblk[i];

				int npartsloc;

				if ((nchildscurr > 1 && itype == compnode) || itype == compnode2) {
					npartsloc = (wtotal+_weightcellboundary-1) / _weightcellboundary;
				} else {
					npartsloc = (wtotal+_weightcellvolume-1) / _weightcellvolume;
				};

				if (npartsloc < 1) npartsloc = 1;

// Compute new ordering and partitioning

				if (npartsloc > 1) {

					int *blksloc, *partitionloc, *orderloc, *iorderloc;

					blksloc = new int [npartsloc+2];
					if (!blksloc) MemoryFail (funcname);
					partitionloc = new int [nlocblk];
					if (!partitionloc) MemoryFail (funcname);
					orderloc = new int [nlocblk];
					if (!orderloc) MemoryFail (funcname);
					iorderloc = new int [nlocblk];
					if (!iorderloc) MemoryFail (funcname);

					alocblk.PartitionMtr (npartsloc, ndwlocblk, adjwlocblk, blksloc, partitionloc, orderloc);

					for (i=0;i<nlocblk;i++) iorderloc[orderloc[i]] = i;

// Free submatrix data

					alocblk = mtrdummy;

// Update the ordering

					for (i=isupbegglob;i<=isupendglob;i++) {
						int iold = iorderloc[i-isupbegglob];
						orderloc[i-isupbegglob] = order[iold+isupbegglob];
					};

					for (i=isupbegglob;i<=isupendglob;i++) {
						order[i] = orderloc[i-isupbegglob];
					};

// Register partitioning data

					for (iblk=0;iblk<npartsloc;iblk++) {
						icycle++;
						for (i=blksloc[iblk];i<blksloc[iblk+1];i++) {
							imask[isupbegglob+i] = icycle;
						};
					};

// Free local arrays

					delete [] blksloc;
					delete [] partitionloc;
					delete [] orderloc;
					delete [] iorderloc;

				} else {

					icycle++;
					for (i=isupbegglob;i<=isupendglob;i++) {
						imask[i] = icycle;
					};

				};

				delete [] ndwlocblk;
				delete [] adjwlocblk;

			};

// Free work arrays

			delete [] ndwloc;
			delete [] adjwloc;

			delete [] blksindep;
			delete [] orderindep;
			delete [] iorderindep;

		};

	};

// Inverse order

	for (i=0;i<nsupr;i++) _order[order[i]] = i;

// Check blocks partitioning

	for (i=0;i<nsupr;i++) {
		if (imask[i] == -1) {
			throw " PartitionOrdMtrSchur_00: not all nodes were processed";
		};
	};

// Transform local partitioning data into the global partitioning data

	_blks = new int [icycle+2];
	if (!_blks) MemoryFail (funcname);

	_blks[0] = 0;

	iblk = 0;

	i = 0;
	iblkloc = imask[0];

	while (i<nsupr) {
		if (imask[i] != iblkloc) {
			_blks[iblk+1] = i;
			if (i != nsupr-1) iblk++;
			if (i<nsupr-1) {
				iblkloc = imask[i];
			};
		};
		i++;
		if (i == nsupr-1 && nsupr > 1) {
			_blks[iblk+1] = nsupr;
			iblk++;
		} else if (nsupr == 1) {
			_blks[iblk+1] = 1;
			iblk++;
		};
	};

	_nblks = iblk;

// Free work arrays

	delete [] order;
	delete [] iorder;
	delete [] imask;
	delete [] adjwarr;
	delete [] ndwarr;
	delete [] adjwarrext;
	delete [] ndwarrext;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix by METIS without bordering splitting according to the prescribed matrix tree
//========================================================================================
void CSMatrix::PartitionOrdMtrSchur (CTree &_tree, int *_order, // Partition matrix by METIS according to the prescribed matrix tree
								int _weightcellvolume, int _weightcellboundary, int &_nblks, int *&_blks) const {

	const char *funcname = "PartitionOrdMtrSchur_01";

	int *ndweights;
	int *adjweights;

	ndweights = new int [nsupr];
	if (!ndweights) MemoryFail (funcname);
	adjweights = new int [nzja];
	if (!adjweights) MemoryFail (funcname);

	int i;

	for (i=0;i<nsupr;i++) ndweights [i] = 1;
	for (i=0;i<nzja; i++) adjweights[i] = 1;

	PartitionOrdMtrSchur (_tree, ndweights, adjweights, _order,
						_weightcellvolume, _weightcellboundary, _nblks, _blks);

	delete [] ndweights;
	delete [] adjweights;

};

// Author: Kharchenko S.A.
// CSMatrix: Sparsify the extended block sparsity assuming binary tree partitioning
//========================================================================================
void CSMatrix::SparsifyExtendedBinaryTree (const CTree &_tree, // Sparsify Schur part of the extended block sparsity assuming binary tree partitioning
															int _nblks, int *_iabext, int *_jabext,
															int *_iabextflt, int *_jabextflt) {

	const char *funcname = "SparsifyExtendedBinaryTree";

// Allocate mask array

	int nzjaloc = _iabext[_nblks];

	int *jamask;

	jamask = new int [nzjaloc];
	if (!jamask) MemoryFail (funcname);

	int i;

	for (i=0;i<nzjaloc;i++) jamask[i] = 1;

// Scan the tree and change mask values

	int ibeg1, iend1, ibeg2, iend2, inode, inode1_ch1, inode1_ch2;
	int j, jj;

	for (inode=0;inode<_tree.nnodes;inode++) {
		if (_tree.nodes[inode].nchilds == 2) {
			inode1_ch1 = _tree.nodes[inode].childs[0];
			inode1_ch2 = _tree.nodes[inode].childs[1];
//			ibeg1 = _tree.nodes[inode1_ch1].indbeg;
//			iend1 = _tree.nodes[inode1_ch1].indend;
//			ibeg2 = _tree.nodes[inode1_ch2].indbeg;
//			iend2 = _tree.nodes[inode1_ch2].indend;
			ibeg1 = _tree.nodes[inode1_ch1].indbegtot;
			iend1 = _tree.nodes[inode1_ch1].indendtot;
			ibeg2 = _tree.nodes[inode1_ch2].indbegtot;
			iend2 = _tree.nodes[inode1_ch2].indendtot;
			for (i=ibeg1;i<=iend1;i++) {
				for (j=_iabext[i];j<_iabext[i+1];j++) {
					jj = _jabext[j];
					if (jj >= ibeg2 && jj <= iend2) jamask[j] = 0;
				};
			};
		};
	};

// Compute final arrays

	for (i=0;i<=_nblks;i++) _iabextflt[i] = 0;

	for (i=0;i<_nblks;i++) {
		for (j=_iabext[i];j<_iabext[i+1];j++) {
			if (jamask[j] > 0) _iabextflt[i+1]++;
		};
	};

	for (i=0;i<_nblks;i++) _iabextflt[i+1] = _iabextflt[i]+_iabextflt[i+1];

	int nz = 0;

	for (i=0;i<_nblks;i++) {
		for (j=_iabext[i];j<_iabext[i+1];j++) {
			if (jamask[j] > 0) {
				_jabextflt[nz] = _jabext[j];
				nz++;
			};
		};
	};

// Free work arrays

	delete [] jamask;

};

// Author: Kharchenko S.A.
// CSMatrix: Sparsify Schur part of the extended block sparsity assuming binary tree partitioning
//========================================================================================
void CSMatrix::SparsifyExtendedSchurBinaryTree (const CTree &_tree, // Sparsify Schur part of the extended block sparsity assuming binary tree partitioning
																int _nblks, int *_iabext, int *_jabext,
																int *_iabextflt, int *_jabextflt) {

	const char *funcname = "SparsifyExtendedSchurBinaryTree";

// Allocate mask array

	int nzjaloc = _iabext[_nblks];

	int *jamask;

	jamask = new int [nzjaloc];
	if (!jamask) MemoryFail (funcname);

	int i;

	for (i=0;i<nzjaloc;i++) jamask[i] = 1;

// Scan the tree and change mask values

	int ibeg1, iend1, ibeg2, iend2, inode, indtreeloc, inode1, inode1_ch1, inode1_ch2;
	int j, jj;

	for (inode=0;inode<_tree.nnodes;inode++) {
		indtreeloc = _tree.nodes[inode].indtree;
		if (indtreeloc >= 0) {
			for (inode1=0;inode1<_tree.subtreearr[indtreeloc].nnodes;inode1++) {
				if (_tree.subtreearr[indtreeloc].nodes[inode1].nchilds == 2) {
					inode1_ch1 = _tree.subtreearr[indtreeloc].nodes[inode1].childs[0];
					inode1_ch2 = _tree.subtreearr[indtreeloc].nodes[inode1].childs[1];
					ibeg1 = _tree.subtreearr[indtreeloc].nodes[inode1_ch1].indbeg;
					iend1 = _tree.subtreearr[indtreeloc].nodes[inode1_ch1].indend;
					ibeg2 = _tree.subtreearr[indtreeloc].nodes[inode1_ch2].indbeg;
					iend2 = _tree.subtreearr[indtreeloc].nodes[inode1_ch2].indend;
					for (i=ibeg1;i<=iend1;i++) {
						for (j=_iabext[i];j<_iabext[i+1];j++) {
							jj = _jabext[j];
							if (jj >= ibeg2 && jj <= iend2) jamask[j] = 0;
						};
					};
				};
			};
		};
	};

// Compute final arrays

	for (i=0;i<=_nblks;i++) _iabextflt[i] = 0;

	for (i=0;i<_nblks;i++) {
		for (j=_iabext[i];j<_iabext[i+1];j++) {
			if (jamask[j] > 0) _iabextflt[i+1]++;
		};
	};

	for (i=0;i<_nblks;i++) _iabextflt[i+1] = _iabextflt[i]+_iabextflt[i+1];

	int nz = 0;

	for (i=0;i<_nblks;i++) {
		for (j=_iabext[i];j<_iabext[i+1];j++) {
			if (jamask[j] > 0) {
				_jabextflt[nz] = _jabext[j];
				nz++;
			};
		};
	};

// Free work arrays

	delete [] jamask;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
