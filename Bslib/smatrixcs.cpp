//------------------------------------------------------------------------------------------------
// File: smatrixcs.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "smatrix.h"
#include "gsmatrix.h"
#include "globals.h"
#include "slvparam.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixCS: Memory allocation zero data constructor
//========================================================================================
CSMatrixCS::CSMatrixCS (): CSMatrixC () { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixCS_00";

	blamx = 0;
	sprndr = new int [1];
	if (!sprndr) MemoryFail (funcname);
	sprndr[0] = 0;
	sprndc = new int [1];
	if (!sprndc) MemoryFail (funcname);
	sprndc[0] = 0;
	bsa = new int [1];
	if (!bsa) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixCS: Memory allocation zero data constructor
//========================================================================================
CSMatrixCS::CSMatrixCS (int _nsupa, int _nzja, int _nza): CSMatrixC (_nsupa,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixCS_01";

	dcmplx czero (0.0e0,0.0e0);

	delete [] a;

	a = new dcmplx [_nza];
	if (!a) MemoryFail (funcname);

	int i;

	for (i=0;i<_nza;i++) {
		a[i] = czero;
	};
	nza = _nza;
	nzatot = nza;
	blamx = 0;
	sprndr = new int [_nsupa+1];
	if (!sprndr) MemoryFail (funcname);
	for (i=0;i<=nsupr;i++) {
		sprndr[i] = 0;
	};
	sprndc = new int [_nsupa+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<=nsupc;i++) {
		sprndc[i] = 0;
	};
	bsa = new int [nzja];
	if (!bsa) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		bsa[i] = 0;
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: // Laplace^2 structure unsymmetric complex matrix constructor with diagonal shift
//========================================================================================
CSMatrixCS::CSMatrixCS (int _nn, double _dshft, // Laplace^2 structure unsymmetric complex matrix constructor with diagonal shift
						double _dunsyr, double _dunsyi,
						double _alphar, double _alphai,
						int _nsupamx, const int *_sprnda): CSMatrixC () {

	const char *funcname = "CSMatrixCS_02";

	blamx = 0;
	sprndr = new int [1];
	if (!sprndr) MemoryFail (funcname);
	sprndc = new int [1];
	if (!sprndc) MemoryFail (funcname);
	bsa = new int [1];
	if (!bsa) MemoryFail (funcname);

// Create real and imaginary parts of the matrix

	CSMatrixRS mtrdummy;

	CSMatrixRS ar (_nn, _dshft, _dunsyr, _nsupamx, _sprnda); 
	CSMatrixRS ai (_nn, _dshft, _dunsyi, _nsupamx, _sprnda); 

// Allocate temporary data

	CSMatrixCS temp (ar.nsupr,ar.nzja,ar.nza);

// Assign the sparsity

	temp.CSMatrix::operator=((const CSMatrix &) ar);

// Assign supernodes and bsa arrays

	int i;

	for (i=0;i<=ar.nsupr;i++) temp.sprndc[i] = ar.sprndr[i];
	for (i=0;i<=ar.nsupr;i++) temp.sprndr[i] = ar.sprndr[i];

	for (i=0;i<ar.nzja;i++) temp.bsa[i] = ar.bsa[i];

	temp.blamx = ar.blamx;

// Assign the elements

	for (i=0;i<ar.nza;i++) {
		temp.a[i].real(_alphar*ar.a[i]);
		temp.a[i].imag(_alphai*ai.a[i]);
	};

// Free the data

	ar = mtrdummy;
	ai = mtrdummy;

// Return the result

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Copy constructor
//========================================================================================
CSMatrixCS::CSMatrixCS (const CSMatrixCS &_aa): CSMatrixC (_aa) { // Copy constructor

	const char *funcname = "CSMatrixCS_copy";

	blamx = _aa.blamx;
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	int i;
	for (i=0;i<=nsupr;i++) {
		sprndr[i] = _aa.sprndr[i];
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<=nsupc;i++) {
		sprndc[i] = _aa.sprndc[i];
	};
	bsa = new int [nzja];
	if (!bsa) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		bsa[i] = _aa.bsa[i];
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Equality operator
//========================================================================================
CSMatrixCS &CSMatrixCS::operator= (const CSMatrixCS &_aa) { // Equality operator

	const char *funcname = "CSMatrixCS_=";

	delete [] sprndr;
	delete [] sprndc;
	delete [] bsa;

	this->CSMatrixC::operator= ((const CSMatrixC &) _aa);

	blamx = _aa.blamx;
	sprndr = new int [nsupr+1];
	if (!sprndr) MemoryFail (funcname);
	int i;
	for (i=0;i<=nsupr;i++) {
		sprndr[i] = _aa.sprndr[i];
	};
	sprndc = new int [nsupc+1];
	if (!sprndc) MemoryFail (funcname);
	for (i=0;i<=nsupc;i++) {
		sprndc[i] = _aa.sprndc[i];
	};
	bsa = new int [nzja];
	if (!bsa) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		bsa[i] = _aa.bsa[i];
	};

	return *this;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Construct complex matrix from two real matrices with the same sparsity
//========================================================================================
CSMatrixCS CSMatrixCS::RS2CS (CSMatrixRS &_ar, CSMatrixRS &_ai, // Construct complex matrix from two real matrices with the same sparsity
								double _alphar, double _alphai) {

//	const char *funcname = "RS2CS";

// Allocate temporary data

	CSMatrixCS temp (_ar.nsupr,_ar.nzja,_ar.nza);

// Assign the sparsity

	temp.CSMatrix::operator=((const CSMatrix &) _ar);

// Assign supernodes and bsa arrays

	int i;

	for (i=0;i<=_ar.nsupr;i++) temp.sprndc[i] = _ar.sprndr[i];
	for (i=0;i<=_ar.nsupr;i++) temp.sprndr[i] = _ar.sprndr[i];

	for (i=0;i<_ar.nzja;i++) temp.bsa[i] = _ar.bsa[i];

	temp.blamx = _ar.blamx;

// Assign the elements

	for (i=0;i<_ar.nza;i++) {
		temp.a[i].real(_alphar*_ar.a[i]);
		temp.a[i].imag(_alphai*_ai.a[i]);
	};

// Return the result

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Transform matrix from point to sprnda format
//========================================================================================
CSMatrixCS CSMatrixCS::C2CS (const CSMatrixC &_ac, int _nsupamx, const int *_sprnda) { // Transform matrix from point to sprnda format

	const char *funcname = "C2CS";

	dcmplx czero (0.0e0,0.0e0);

// Check input sprnda array

	int nloc = _ac.n;

	if (_sprnda[_nsupamx] < nloc) throw " C2CS: Error in input array sprnda ";

// Compute nsupa

	int nsupa = 0;

	int isup, jsup, jj;

	for (isup=1;isup<=_nsupamx;isup++) {
		if (_sprnda[isup] >= nloc) {
			nsupa = isup;
			break;
		};
	};

	if (nsupa == 0) throw " C2CS: Error in input array sprnda ";

// Compute nd2sp and sp2nd arrays

	int *nd2sp, *sp2nd;
	int i, j;

	nd2sp = new int [nloc];
	if (!nd2sp) MemoryFail (funcname);
	sp2nd = new int [nsupa+1];
	if (!sp2nd) MemoryFail (funcname);

	sp2nd[0] = 0;
	for (i=0;i<nsupa-1;i++) {
		sp2nd[i+1] = _sprnda[i+1];
	};
	sp2nd[nsupa] = nloc;

	for (i=0;i<nsupa;i++) {
		for (j=sp2nd[i];j<sp2nd[i+1] && j<nloc;j++) {
			nd2sp[j] = i;
		};
	};

// Count the total number of elements in the collapsed matrix

	int icycle=0;

	int *imask, *lstloc;

	imask = new int [nsupa];
	if (!imask) MemoryFail (funcname);
	lstloc = new int [nsupa];
	if (!lstloc) MemoryFail (funcname);

	for (i=0;i<nsupa;i++) imask[i] = icycle;

	int nzj = 0;
	int nz = 0;
	int nlstloc;

	int blai, blaj;

	for (isup=0;isup<nsupa;isup++) {
		blai = sp2nd[isup+1]-sp2nd[isup];
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = nd2sp[jj];
				if (imask[jsup] != icycle) {
					lstloc[nlstloc] = jsup;
					imask[jsup] = icycle;
					nlstloc++;
				};
			};
		};
		nzj += nlstloc;
		for (j=0;j<nlstloc;j++) {
			jsup = lstloc[j];
			blaj = sp2nd[jsup+1]-sp2nd[jsup];
			nz += blai*blaj;
		};
	};

// Compute sparsity structure of the collapsed matrix

	CSMatrixCS temp (nsupa, nzj, nz);

	temp.ia[0] = 0;

	for (isup=0;isup<nsupa;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = nd2sp[jj];
				if (imask[jsup] != icycle) {
					lstloc[nlstloc] = jsup;
					imask[jsup] = icycle;
					nlstloc++;
				};
			};
		};
		if (nlstloc != 0) qsort (lstloc, nlstloc, sizeof(int), compint);
		for (j=0;j<nlstloc;j++) temp.ja[temp.ia[isup]+j] = lstloc[j];
		temp.ia[isup+1] = temp.ia[isup]+nlstloc;
	};
	for (i=0;i<nsupa;i++) temp.list[i] = i;

// Init arrays sprndr, sprndc and bsa

	for (isup=0;isup<=nsupa;isup++) {
		temp.sprndr[isup] = sp2nd[isup];
		temp.sprndc[isup] = sp2nd[isup];
	};

	nz = 0;

	temp.blamx = 0;

	for (isup=0;isup<nsupa;isup++) {
		blai = sp2nd[isup+1]-sp2nd[isup];
		if (blai > temp.blamx) temp.blamx = blai;
		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			blaj = sp2nd[jsup+1]-sp2nd[jsup];
			temp.bsa[j] = nz;
			nz += blai*blaj;
		};
	};

// Fill elements of the collapsed matrix

// Init elements by zeroes

	for (i=0;i<nz;i++) temp.a[i] = czero;

// Fill all remaining elements

	for (isup=0;isup<nsupa;isup++) {
		blai = sp2nd[isup+1]-sp2nd[isup];

// Init addresses of the nonzero supernodes

		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			imask[jsup] = temp.bsa[j];
		};

// Fill elements of the supernode row

		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			int kii = i-sp2nd[isup];
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = nd2sp[jj];
				blaj = sp2nd[jsup+1]-sp2nd[jsup];
				int kjj = jj-sp2nd[jsup];
				int ibs = imask[jsup];
				temp.a[ibs+kjj*blai+kii] = _ac.a[j];
			};
		};

	};

// Reassign m and n

	temp.m = sp2nd[nsupa];
	temp.n = sp2nd[nsupa];

// Delete work arrays

	delete [] sp2nd;
	delete [] nd2sp;
	delete [] imask;
	delete [] lstloc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Transform matrix stored by columns from point to sprnda format
//========================================================================================
CSMatrixCS CSMatrixCS::C2CSRectCol (const CSMatrixC &_ac, int _nsupamx, const int *_sprnda) { // Transform matrix stored by columns from point to sprnda format

	const char *funcname = "C2CSRectCol";

	dcmplx czero (0.0e0,0.0e0);

// Check input sprnda array

	int nloc = _ac.n;
	int mloc = _ac.m;

	if (_sprnda[_nsupamx] < nloc) throw " C2CSRectCol: Error in input array sprnda ";
	if (_sprnda[_nsupamx] < mloc) throw " C2CSRectCol: Error in input array sprnda ";

// Compute nsupa

	int isup, jsup, jj;

	int nsupcloc = 0;

	for (isup=1;isup<=_nsupamx;isup++) {
		if (_sprnda[isup] >= nloc) {
			nsupcloc = isup;
			break;
		};
	};

	if (nsupcloc == 0) throw " C2CSRectCol: Error in input array sprnda ";

	int nsuprloc = 0;

	for (isup=1;isup<=_nsupamx;isup++) {
		if (_sprnda[isup] >= mloc) {
			nsuprloc = isup;
			break;
		};
	};

	if (nsuprloc == 0) throw " C2CSRectCol: Error in input array sprnda ";

// Compute nd2sp and sp2nd arrays

	int *ndc2sp, *spc2nd;
	int *ndr2sp, *spr2nd;
	int i, j;

	ndc2sp = new int [nloc];
	if (!ndc2sp) MemoryFail (funcname);
	spc2nd = new int [nsupcloc+1];
	if (!spc2nd) MemoryFail (funcname);
	ndr2sp = new int [mloc];
	if (!ndr2sp) MemoryFail (funcname);
	spr2nd = new int [nsuprloc+1];
	if (!spr2nd) MemoryFail (funcname);

	spc2nd[0] = 0;
	for (i=0;i<nsupcloc-1;i++) {
		spc2nd[i+1] = _sprnda[i+1];
	};
	spc2nd[nsupcloc] = nloc;

	for (i=0;i<nsupcloc;i++) {
		for (j=spc2nd[i];j<spc2nd[i+1] && j<nloc;j++) {
			ndc2sp[j] = i;
		};
	};

	spr2nd[0] = 0;
	for (i=0;i<nsuprloc-1;i++) {
		spr2nd[i+1] = _sprnda[i+1];
	};
	spr2nd[nsuprloc] = mloc;

	for (i=0;i<nsuprloc;i++) {
		for (j=spr2nd[i];j<spr2nd[i+1] && j<mloc;j++) {
			ndr2sp[j] = i;
		};
	};

// Count the total number of elements in the collapsed matrix

	int icycle=0;

	int *imaskr, *lstrloc;

	imaskr = new int [nsuprloc];
	if (!imaskr) MemoryFail (funcname);
	lstrloc = new int [nsuprloc];
	if (!lstrloc) MemoryFail (funcname);

	for (i=0;i<nsuprloc;i++) imaskr[i] = icycle;

	int nzj = 0;
	int nz = 0;
	int nlstloc;

	int blai, blaj;

	for (isup=0;isup<nsupcloc;isup++) {
		blai = spc2nd[isup+1]-spc2nd[isup];
		icycle++;
		nlstloc = 0;
		for (i=spc2nd[isup];i<spc2nd[isup+1] && i<nloc;i++) {
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = ndr2sp[jj];
				if (imaskr[jsup] != icycle) {
					lstrloc[nlstloc] = jsup;
					imaskr[jsup] = icycle;
					nlstloc++;
				};
			};
		};
		nzj += nlstloc;
		for (j=0;j<nlstloc;j++) {
			jsup = lstrloc[j];
			blaj = spr2nd[jsup+1]-spr2nd[jsup];
			nz += blai*blaj;
		};
	};

// Compute sparsity structure of the collapsed matrix

	int nsupmx = nsupcloc;
	if (nsuprloc > nsupmx) nsupmx = nsuprloc;

	CSMatrixCS temp (nsupmx, nzj, nz);

	temp.ia[0] = 0;

	for (isup=0;isup<nsupcloc;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=spc2nd[isup];i<spc2nd[isup+1] && i<nloc;i++) {
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = ndr2sp[jj];
				if (imaskr[jsup] != icycle) {
					lstrloc[nlstloc] = jsup;
					imaskr[jsup] = icycle;
					nlstloc++;
				};
			};
		};
		if (nlstloc != 0) qsort (lstrloc, nlstloc, sizeof(int), compint);
		for (j=0;j<nlstloc;j++) temp.ja[temp.ia[isup]+j] = lstrloc[j];
		temp.ia[isup+1] = temp.ia[isup]+nlstloc;
	};
	for (i=0;i<nsupcloc;i++) temp.list[i] = i;

// Init arrays sprndr, sprndc and bsa

	for (isup=0;isup<=nsuprloc;isup++) {
		temp.sprndr[isup] = spr2nd[isup];
	};
	for (isup=0;isup<=nsupcloc;isup++) {
		temp.sprndc[isup] = spc2nd[isup];
	};

	nz = 0;

	temp.blamx = 0;

	for (isup=0;isup<nsupcloc;isup++) {
		blai = spc2nd[isup+1]-spc2nd[isup];
		if (blai > temp.blamx) temp.blamx = blai;
		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			blaj = spr2nd[jsup+1]-spr2nd[jsup];
			temp.bsa[j] = nz;
			nz += blai*blaj;
		};
	};

// Fill elements of the collapsed matrix

// Init elements by zeroes

	for (i=0;i<nz;i++) temp.a[i] = czero;

// Fill all remaining elements

	for (isup=0;isup<nsupcloc;isup++) {
		blai = spc2nd[isup+1]-spc2nd[isup];

// Init addresses of the nonzero supernodes

		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			imaskr[jsup] = temp.bsa[j];
		};

// Fill elements of the supernode row

		for (i=spc2nd[isup];i<spc2nd[isup+1] && i<nloc;i++) {
			int kii = i-spc2nd[isup];
			for (j=_ac.ia[i];j<_ac.ia[i+1];j++) {
				jj = _ac.ja[j];
				jsup = ndr2sp[jj];
				blaj = spr2nd[jsup+1]-spr2nd[jsup];
				int kjj = jj-spr2nd[jsup];
				int ibs = imaskr[jsup];
				temp.a[ibs+kjj*blai+kii] = _ac.a[j];
			};
		};

	};

// Reassign m and n

	temp.m = spr2nd[nsuprloc];
	temp.n = spc2nd[nsupcloc];
	temp.nsupc = nsupcloc;
	temp.nsupr = nsuprloc;
	temp.nlist = nsupcloc;

// Delete work arrays

	delete [] spc2nd;
	delete [] ndc2sp;
	delete [] spr2nd;
	delete [] ndr2sp;
	delete [] imaskr;
	delete [] lstrloc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Copy matrix into disk
//========================================================================================
CSMatrixCS CSMatrixCS::CS2Disk (int _nblks, int *_blks, FILE **_mtrfiles) const { // Copy matrix into disk

	const char *funcname = "CS2Disk";

	CSMatrixCS temp;

	temp = *this;

// Remove a array from the memory

	delete [] temp.a;

	temp.nza = 0;

	temp.a = new dcmplx [temp.nza];
	if (!temp.a) MemoryFail (funcname);

// Copy a array to the disk and reassign bsa

	int nztot = 0;

	for (int iblk=0;iblk<_nblks;iblk++) {

// Reassign bsa for the current block

		int nz = 0;
		for (int isup=_blks[iblk];isup<_blks[iblk+1];isup++) {
			int blai=sprndc[isup+1]-sprndc[isup];
			for (int j=ia[isup];j<ia[isup+1];j++) {
				int jsup = ja[j];
				int blaj=sprndr[jsup+1]-sprndr[jsup];
				temp.bsa[j] = nz;
				nz += blai*blaj;
			};
		};

// Store current block

		int i=0;

		FPut (_mtrfiles[iblk], nz, a+nztot, i);

// Register

		nztot += nz;

	};

	temp.nza = 0;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Matrix ordering routine
//========================================================================================
CSMatrixCS CSMatrixCS::OrdMtr (const int *_order) const { // Matrix ordering routine

	const char *funcname = "OrdMtr_00";

// Compute inverse order

	int *iord;
	dcmplx **pdarr;

	iord = new int [nsupr];
	if (!iord) MemoryFail (funcname);
	pdarr = new dcmplx * [nsupr];
	if (!pdarr) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupr;i++) iord[_order[i]] = i;

// Reorder the matrix

	CSMatrixCS an(nsupr,nzja,nza);

// Compute new arrays sprndr and sprndc

	for (i=0;i<=nsupr;i++) an.sprndr[i] = 0;

	int j;

	for (i=0;i<nsupr;i++) {
		j = _order[i];
		an.sprndr[j+1] = sprndr[i+1]-sprndr[i];
	};

	for (i=0;i<nsupr;i++) an.sprndr[i+1] += an.sprndr[i];
	for (i=0;i<=nsupr;i++) an.sprndc[i] = an.sprndr[i];

// Compute reordered sparsity structure

	for (i=0;i<nsupr;i++) an.list[i] = i;

	for (i=0;i<=nsupr;i++) an.ia[i] = 0;

	for (i=0;i<nsupr;i++) {
		j = _order[i];
		an.ia[j+1] = ia[i+1]-ia[i];
	};

	for (i=0;i<nsupr;i++) an.ia[i+1] += an.ia[i];

	int inew, nzloc, ibs, jold, jnew, kii;
	int *pja;
	dcmplx *paloc;

	int blai, blaj;

	int nzj = 0;
	int nz = 0;

	for (inew=0;inew<nsupr;inew++) {
		i = iord[inew];
		blai = sprndr[i+1]-sprndr[i];
		nzloc = ia[i+1]-ia[i];
		ibs = nzj;
		for (j=ia[i];j<ia[i+1];j++) {
			jold = ja[j];
			pdarr[j-ia[i]] = a+bsa[j];
			jnew = _order[jold];
			an.ja[nzj] = jnew;
			nzj++;
		};

// Sort elements

		pja = an.ja;

		if (nzloc != 0) SortPtr (nzloc, pja+ibs, pdarr);

		for (j=an.ia[inew];j<an.ia[inew+1];j++) {
			jnew = an.ja[j];
			blaj = an.sprndr[jnew+1]-an.sprndr[jnew];
			an.bsa[j] = nz;
			paloc = pdarr[j-an.ia[inew]];
			for (kii=0;kii<blai*blaj;kii++) an.a[nz+kii] = paloc[kii];
			nz += blai*blaj;
		};

	};

// Init control data

	an.m      = m;
	an.n      = n;
	an.nsupc  = nsupc;
	an.nsupr  = nsupr;
	an.nlist  = nlist;
	an.nzja   = nzja;
	an.nza    = nza;
	an.nzatot = nzatot;
	an.blamx  = blamx;

// Delete ordering array

	delete [] iord;
	delete [] pdarr;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Matrix ordering routine
//========================================================================================
CSMatrixCS CSMatrixCS::OrdMtr (int _nblks, int *_blks, FILE **_mtrfiles, // Matrix ordering routine
								const int *_order, 
								FILE **_mtrofiles) const { 

	const char *funcname = "OrdMtr_01";

// Compute inverse order

	int *iord;
	int *iarr;
	int *sp2bl;
	dcmplx *celem;

	iord = new int [nsupr];
	if (!iord) MemoryFail (funcname);
	iarr = new int [nsupr];
	if (!iarr) MemoryFail (funcname);
	sp2bl = new int [nsupr];
	if (!sp2bl) MemoryFail (funcname);
	celem = new dcmplx [blamx*blamx];
	if (!celem) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupr;i++) iord[_order[i]] = i;

// Assign sp2bl

	for (int iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			sp2bl[i] = iblk;
		};
	};

// Reorder the matrix

	int nzaloc=0;

	CSMatrixCS an(nsupr,nzja,nzaloc);

// Compute new arrays sprndr and sprndc

	for (i=0;i<=nsupr;i++) an.sprndr[i] = 0;

	int j;

	for (i=0;i<nsupr;i++) {
		j = _order[i];
		an.sprndr[j+1] = sprndr[i+1]-sprndr[i];
	};

	for (i=0;i<nsupr;i++) an.sprndr[i+1] += an.sprndr[i];
	for (i=0;i<=nsupr;i++) an.sprndc[i] = an.sprndr[i];

// Compute reordered sparsity structure

	for (i=0;i<nsupr;i++) an.list[i] = i;

	for (i=0;i<=nsupr;i++) an.ia[i] = 0;

	for (i=0;i<nsupr;i++) {
		j = _order[i];
		an.ia[j+1] = ia[i+1]-ia[i];
	};

	for (i=0;i<nsupr;i++) an.ia[i+1] += an.ia[i];

	int inew, nzloc, ibs, jold, jnew;
	int iblknw, iblkol;
	int *pja;

	int blai, blaj;

	int nzj = 0;
	int nz;

	for (iblknw=0;iblknw<_nblks;iblknw++) {
		nz = 0;
		for (inew=_blks[iblknw];inew<_blks[iblknw+1];inew++) {
			i = iord[inew];
			iblkol = sp2bl[i];
			blai = sprndr[i+1]-sprndr[i];
			nzloc = ia[i+1]-ia[i];
			ibs = nzj;
			for (j=ia[i];j<ia[i+1];j++) {
				jold = ja[j];
				iarr[j-ia[i]] = bsa[j];
				jnew = _order[jold];
				an.ja[nzj] = jnew;
				nzj++;
			};

// Sort elements

			pja = an.ja;

			if (nzloc != 0) SortElm (nzloc, pja+ibs, iarr);

// Store the data

			for (j=an.ia[inew];j<an.ia[inew+1];j++) {
				jnew = an.ja[j];
				blaj = an.sprndr[jnew+1]-an.sprndr[jnew];
				an.bsa[j] = nz;
				ibs = iarr[j-an.ia[inew]];

				FGet (_mtrfiles[iblkol],blai*blaj,celem,ibs);

				FPut (_mtrofiles[iblknw],blai*blaj,celem,nz);

				nz += blai*blaj;
			};
		};
	};

// Init control data

	an.m      = m;
	an.n      = n;
	an.nsupc  = nsupc;
	an.nsupr  = nsupr;
	an.nlist  = nlist;
	an.nzja   = nzja;
	an.nza    = nza;
	an.nzatot = nzatot;
	an.blamx  = blamx;

// Delete ordering array

	delete [] iord;
	delete [] iarr;
	delete [] sp2bl;
	delete [] celem;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute the elements of the lower triangular part of the AtA matrix
//========================================================================================
CSMatrixCS CSMatrixCS::AtAMatrElem (ofstream &_fout, CSlvParam &_param) const { // Compute the elements of the lower triangular part of the AtA matrix

	const char *funcname = "AtAMatrElem";

	dcmplx czero (0.0e0,0.0e0);

	double ops, perf, tottim;
	clock_t time0, time1;
//
// Init statistics data
//
	time0 = clock ();

	int nzaloc = GetNza();

// Compute the sparsity of the AtA matrix

	CSMatrix atastr;

	atastr = AtAMatrSpars();

// Compute the number of elements in the AtA matrix

	int isup, j, jsup, blai, blaj;

	int nz = 0;

	for (isup=0;isup<nsupc;isup++) {
		blai = sprndc[isup+1]-sprndc[isup];
		for (j=atastr.ia[isup];j<atastr.ia[isup+1];j++) {
			jsup = atastr.ja[j];
			blaj = sprndc[jsup+1]-sprndc[jsup];
			nz += blai*blaj;
		};
	};

// Allocate AtA matrix

	CSMatrixCS ata (nsupc,atastr.nzja,nz);

// Copy the sparsity

	ata.CSMatrix::operator=(atastr);

// Free the sparsity

	CSMatrix mtrdummy;

	atastr = mtrdummy;

// Allocate temporary data

	dcmplx *elem;

	elem = new dcmplx [blamx*blamx];
	if (!elem) MemoryFail (funcname);

// Compute the elements

	nz = 0;

	int blak, kii, kjj, kkk, ibeg1, iend1, ibeg2, iend2;
	int jind1, jind2, jj1, jj2, ibs, jbs, i, blaij;
	dcmplx *pelem, *pai, *paj;

	ops = 0.0e0;

	for (isup=0;isup<nsupc;isup++) {

		blai = sprndc[isup+1]-sprndc[isup];
		for (j=ata.ia[isup];j<ata.ia[isup+1];j++) {

			ata.bsa[j] = nz;

			jsup = ata.ja[j];
			blaj = sprndc[jsup+1]-sprndc[jsup];

			blaij = blai*blaj;

			pelem = elem;
			for (kii=0;kii<blaij;kii++) *pelem++ = czero;

// Search current columns

			ibeg1 = ia[isup]; 
			iend1 = ia[isup+1]-1; 
			ibeg2 = ia[jsup]; 
			iend2 = ia[jsup+1]-1; 
			jind1 = ibeg1;
			jind2 = ibeg2;

			while (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = ja[jind2];
				if (jj1 < jj2) {
					jind1++;
				} else if (jj1 > jj2) {
					jind2++;
				} else if (jj1 == jj2) {

// Update the product

					blak = sprndr[jj1+1]-sprndr[jj1];

					ibs = bsa[jind1];
					jbs = bsa[jind2];

					for (kii=0;kii<blai;kii++) {
						pelem = elem+kii;
						for (kjj=0;kjj<blaj;kjj++) {
							pai = a+ibs+kii;
							paj = a+jbs+kjj;
							for (kkk=0;kkk<blak;kkk++) {
//								elem[kjj*blai+kii] += a[ibs+kkk*blai+kii]*a[jbs+kkk*blaj+kjj];
								*pelem += *pai * conj(*paj);
								pai += blai;
								paj += blaj;
							};
							pelem += blai;
						};
					};

					ops += blaij*blak;

					jind1++;
					jind2++;

				};
			};

// Store computed data

			pelem = elem;
			pai = ata.a+nz;
			for (kii=0;kii<blaij;kii++) *pai++ = *pelem++;

			nz += blaij;

		};
	};

// Store head data

	ata.m = n;
	ata.n = n;
	ata.nsupr = nsupc;
	ata.nsupc = nsupc;
	ata.nlist = nsupc;
	ata.nza = nz;
	ata.nzatot = nz;

	for (i=0;i<nsupc;i++) ata.list[i] = i;
	for (i=0;i<=nsupc;i++) ata.sprndr[i] = sprndc[i];
	for (i=0;i<=nsupc;i++) ata.sprndc[i] = sprndc[i];

	ata.blamx=0;

	for (i=0;i<nsupc;i++) {
		blai = sprndc[i+1]-sprndc[i];
		if (blai > ata.blamx) ata.blamx = blai;
	};

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output AtA computations statistics

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzaloc;

	double dens = (double) nz / (double) nzaloc * 1.0e2;

	if (_param.msglev > 0) {

		cout  << " AtA statistics: " << endl;
		cout  << "     Densty = " << dens << " % " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

		_fout << " AtA statistics: " << endl;
		_fout << "     Densty = " << dens << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	};

// Free work memory

	delete [] elem;

	return ata;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute L and U data for given matrix
//========================================================================================
void CSMatrixCS::CombineLU (CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const { // Compute L and U data for given matrix

//	const char *funcname = "CombineLU_00";

// Init dummy matrix

	CSMatrixCS mtrdummy;

// Prepare U data

// Count the number of elements

	int isup, jsup, nzj, nz, ilist, j, blai, blaj;

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			if (jsup >= isup) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary U

	CSMatrixCS mtru (nsupr,nzj,nz);

	mtru.m = n;
	mtru.n = n;
	mtru.nsupr = nsupr;
	mtru.nsupc = nsupr;
	mtru.nlist = nsupr;
	mtru.nzja = nzj;
	mtru.nza = nz;
	mtru.nzatot = nz;
	mtru.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) mtru.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) mtru.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) mtru.sprndc[isup] = sprndr[isup];

	mtru.ia[0] = 0;

	nzj = 0;
	nz = 0;

	int kii, blaij, ibs;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			blaij = blai*blaj;
			if (jsup >= isup) {
				mtru.ja[nzj] = jsup;
				mtru.bsa[nzj] = nz;
				ibs = bsa[j];
				for (kii=0;kii<blaij;kii++) mtru.a[nz+kii] = a[ibs+kii];
				nzj++;
				nz += blaij;
			};
		};
		mtru.ia[isup+1] = nzj;
	};

// Store U data

	_mtru = mtru;

	mtru = mtrdummy;

// Prepare L data

// Count the number of elements

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			if (jsup <= isup) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary L

	CSMatrixCS mtrl (nsupr,nzj,nz);

	mtrl.m = n;
	mtrl.n = n;
	mtrl.nsupr = nsupr;
	mtrl.nsupc = nsupr;
	mtrl.nlist = nsupr;
	mtrl.nzja = nzj;
	mtrl.nza = nz;
	mtrl.nzatot = nz;
	mtrl.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) mtrl.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) mtrl.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) mtrl.sprndc[isup] = sprndr[isup];

	mtrl.ia[0] = 0;

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			blaij = blai*blaj;
			if (jsup <= isup) {
				mtrl.ja[nzj] = jsup;
				mtrl.bsa[nzj] = nz;
				ibs = bsa[j];
				for (kii=0;kii<blaij;kii++) mtrl.a[nz+kii] = a[ibs+kii];
				nzj++;
				nz += blaij;
			};
		};
		mtrl.ia[isup+1] = nzj;
	};

// Transpose L

	CSMatrixCS mtrlt;

	mtrlt = mtrl.TranspMtrNoConj ();
	mtrl = mtrdummy;

// Copy the result

	_mtrl = mtrlt;

	mtrlt = mtrdummy;

// Extend the sparsity structures of L and U 

	_mtrl.ExtendSparsity (_mtru);
	_mtru.ExtendSparsity (_mtrl);

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute L and U data for given matrix
//========================================================================================
void CSMatrixCS::CombineLU (int _nblks, int *_blks, FILE **_mtrfiles, // Compute L and U data for given matrix
							FILE **_mtrlfiles, FILE **_mtrufiles,
							CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const {

	const char *funcname = "CombineLU_01";

	int *sp2bl;
	dcmplx *celem, *celem1;

	sp2bl = new int [nsupr];
	if (!sp2bl) MemoryFail (funcname);
	celem = new dcmplx [blamx*blamx];
	if (!celem) MemoryFail (funcname);
	celem1 = new dcmplx [blamx*blamx];
	if (!celem1) MemoryFail (funcname);

// Assign sp2bl

	int iblk;
	for (iblk=0;iblk<_nblks;iblk++) {
		for (int i=_blks[iblk];i<_blks[iblk+1];i++) {
			sp2bl[i] = iblk;
		};
	};

// Init dummy matrix

	CSMatrixCS mtrdummy;

// Prepare U data

// Count the number of elements

	int isup, jsup, nzj, nz, ilist, j, blai, blaj;

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			if (jsup >= isup) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary U

	int nzloc = 0;

	CSMatrixCS mtru (nsupr,nzj,nzloc);

	mtru.m = n;
	mtru.n = n;
	mtru.nsupr = nsupr;
	mtru.nsupc = nsupr;
	mtru.nlist = nsupr;
	mtru.nzja = nzj;
	mtru.nza = nzloc;
	mtru.nzatot = nz;
	mtru.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) mtru.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) mtru.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) mtru.sprndc[isup] = sprndr[isup];

	mtru.ia[0] = 0;

	nzj = 0;

	int kii, blaij, ibs;

	for (iblk=0;iblk<_nblks;iblk++) {
		nz = 0;
		for (isup=_blks[iblk];isup<_blks[iblk+1];isup++) {
			blai = sprndr[isup+1]-sprndr[isup];
			for (j=ia[isup];j<ia[isup+1];j++) {
				jsup = ja[j];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				blaij = blai*blaj;
				if (jsup >= isup) {
					mtru.ja[nzj] = jsup;
					mtru.bsa[nzj] = nz;
					ibs = bsa[j];

					FGet (_mtrfiles[iblk],blaij,celem,ibs);

					FPut (_mtrufiles[iblk],blaij,celem,nz);

					nzj++;
					nz += blaij;
				};
			};
			mtru.ia[isup+1] = nzj;
		};
	};

// Store U data

	_mtru = mtru;

	mtru = mtrdummy;

// Prepare L data

// Count the number of elements

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			if (jsup <= isup) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary L

	nzloc = 0;

	CSMatrixCS mtrl (nsupr,nzj,nzloc);

	mtrl.m = n;
	mtrl.n = n;
	mtrl.nsupr = nsupr;
	mtrl.nsupc = nsupr;
	mtrl.nlist = nsupr;
	mtrl.nzja = nzj;
	mtrl.nza = nzloc;
	mtrl.nzatot = nz;
	mtrl.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) mtrl.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) mtrl.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) mtrl.sprndc[isup] = sprndr[isup];

	mtrl.ia[0] = 0;

	nzj = 0;
	nz = 0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			blaij = blai*blaj;
			if (jsup <= isup) {
				mtrl.ja[nzj] = jsup;
				mtrl.bsa[nzj] = bsa[j];
				nzj++;
				nz += blaij;
			};
		};
		mtrl.ia[isup+1] = nzj;
	};

// Transpose the sparsity of L

	int *ilt, *iptr;
	int *jlt, *bslt, *bslt2l;

	ilt = new int [nsupr+1];
	if (!ilt) MemoryFail (funcname);
	iptr = new int [nsupr];
	if (!iptr) MemoryFail (funcname);
	jlt = new int [nzj];
	if (!jlt) MemoryFail (funcname);
	bslt = new int [nzj];
	if (!bslt) MemoryFail (funcname);
	bslt2l = new int [nzj];
	if (!bslt2l) MemoryFail (funcname);

	int i;
	for (i=0;i<=nsupr;i++) ilt[i] = 0;

	for (i=0;i<nsupr;i++) {
		for (j=mtrl.ia[i];j<mtrl.ia[i+1];j++) {
			int jj = mtrl.ja[j];
			ilt[jj+1]++;
		};
	};

	for (i=0;i<nsupr;i++) ilt[i+1] += ilt[i];

	for (i=0;i<nsupr;i++) iptr[i] = ilt[i];

	int k;
	for (i=0;i<nsupr;i++) {
		for (j=mtrl.ia[i];j<mtrl.ia[i+1];j++) {
			int jj = mtrl.ja[j];
			k = iptr[jj];
			jlt[k] = i;
			bslt2l[k] = mtrl.bsa[j];
			iptr[jj]++;
		};
	};

	for (iblk=0;iblk<_nblks;iblk++) {
		nz = 0;
		for (isup=_blks[iblk];isup<_blks[iblk+1];isup++) {
			blai = sprndr[isup+1]-sprndr[isup];
			for (j=ilt[isup];j<ilt[isup+1];j++) {
				jsup = jlt[j];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				blaij = blai*blaj;
				bslt[j] = nz;

				int jblk=sp2bl[jsup];
				ibs = bslt2l[j];

				FGet (_mtrfiles[jblk],blaij,celem,ibs);

				for (kii=0;kii<blai;kii++) {
					for (int kjj=0;kjj<blaj;kjj++) {
						celem1[kjj*blai+kii] = celem[kii*blaj+kjj];
					};
				};

				FPut (_mtrlfiles[iblk],blaij,celem1,nz);

				nz += blaij;
			};
		};
	};

// Copy the result

	for (i=0;i<=nsupr;i++) mtrl.ia[i] = ilt[i];
	for (i=0;i<nzj;i++) mtrl.ja[i] = jlt[i];
	for (i=0;i<nzj;i++) mtrl.bsa[i] = bslt[i];

	_mtrl = mtrl;

	mtrl = mtrdummy;

// Extend the sparsity structures of L and U 

//	_mtrl.ExtendSparsity (_mtru);
//	_mtru.ExtendSparsity (_mtrl);

// Check that sparsity structures of L and U coincide

	for (i=0;i<=nsupr;i++) {
		if (_mtrl.ia[i] != _mtru.ia[i]) {
			cout << " Error in the sparsity of L and U " << endl;
			throw " Error in the sparsity of L and U ";
		};
	};

	for (i=0;i<_mtrl.nzja;i++) {
		if (_mtrl.ja[i] != _mtru.ja[i]) {
			cout << " Error in the sparsity of L and U " << endl;
			throw " Error in the sparsity of L and U ";
		};
	};

// Free work arrays

	delete [] sp2bl;
	delete [] celem;
	delete [] celem1;
	delete [] ilt;
	delete [] iptr;
	delete [] jlt;
	delete [] bslt;
	delete [] bslt2l;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Extend sparsity of the matrix according to the second matrix
//========================================================================================
void CSMatrixCS::ExtendSparsity (const CSMatrixCS &_mtr2) { // Extend sparsity of the matrix according to the second matrix

	const char *funcname = "ExtendSparsity";

	dcmplx czero (0.0e0,0.0e0);

// Check the data on entry

	if (nlist != _mtr2.nlist) throw " Error in sparsity extend: the numbers of nonzero rows do not coincide";
//	if (nlist != nsupr) throw " Error in sparsity extend: wrong number of nonzero rows";

	int ilist;
	for (ilist=0;ilist<nlist;ilist++) {
		if (list[ilist] != _mtr2.list[ilist]) throw " Error in sparsity extend: the lists of nonzero rows do not coincide";
	};

// Count the new extended number of elements in each row

	int nzj, nz, ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;
	int blai, blaj, isup, jsup;

	nzj = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {

		isup = list[ilist];

		blai = sprndr[isup+1]-sprndr[isup];

// Scan both rows at the same time

		ibeg1 = ia[ilist]; 
		iend1 = ia[ilist+1]-1; 
		ibeg2 = _mtr2.ia[ilist]; 
		iend2 = _mtr2.ia[ilist+1]-1; 
		jind1 = ibeg1;
		jind2 = ibeg2;

		while (jind1 <= iend1 || jind2 <= iend2) {
			if (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = _mtr2.ja[jind2];
				if (jj1 < jj2) {
					jsup = jj1;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					nzj++;
					nz += blai*blaj;
					jind1++;
				} else if (jj1 > jj2) {
					jsup = jj2;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					nzj++;
					nz += blai*blaj;
					jind2++;
				} else if (jj1 == jj2) {
					jsup = jj1;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					nzj++;
					nz += blai*blaj;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				jsup = ja[jind1];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				nzj++;
				nz += blai*blaj;
				jind1++;
			} else if (jind2 <= iend2) {
				jsup = _mtr2.ja[jind2];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				nzj++;
				nz += blai*blaj;
				jind2++;
			};
		};

	};

// Allocate work memory for extended matrix

	CSMatrixCS temp (nlist,nzj,nz);

	temp.m = n;
	temp.n = n;
	temp.nsupr = nsupr;
	temp.nsupc = nsupr;
	temp.nlist = nlist;
	temp.nzja = nzj;
	temp.nza = nz;
	temp.nzatot = nz;
	temp.blamx = blamx;

	for (ilist=0;ilist<nlist;ilist++) temp.list[ilist] = list[ilist];

	delete [] temp.sprndr;
	delete [] temp.sprndc;

	temp.sprndr = new int [nsupr+1];
	if(!temp.sprndr) MemoryFail(funcname);
	temp.sprndc = new int [nsupr+1];
	if(!temp.sprndc) MemoryFail(funcname);

	for (isup=0;isup<=nsupr;isup++) temp.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) temp.sprndc[isup] = sprndr[isup];

// Fill ia, ja, bsa and a arrays

	int kii, ibs;

	temp.ia[0] = 0;

	nzj = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {

		isup = list[ilist];

		blai = sprndr[isup+1]-sprndr[isup];

// Scan both rows at the same time

		ibeg1 = ia[ilist]; 
		iend1 = ia[ilist+1]-1; 
		ibeg2 = _mtr2.ia[ilist]; 
		iend2 = _mtr2.ia[ilist+1]-1; 
		jind1 = ibeg1;
		jind2 = ibeg2;

		while (jind1 <= iend1 || jind2 <= iend2) {
			if (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = _mtr2.ja[jind2];
				if (jj1 < jj2) {
					jsup = jj1;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					temp.ja[nzj] = jsup;
					ibs = bsa[jind1];
					temp.bsa[nzj] = nz;
					for (kii=0;kii<blai*blaj;kii++) temp.a [nz+kii] = a[ibs+kii];
					nzj++;
					nz += blai*blaj;
					jind1++;
				} else if (jj1 > jj2) {
					jsup = jj2;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					temp.ja[nzj] = jsup;
					temp.bsa[nzj] = nz;
					for (kii=0;kii<blai*blaj;kii++) temp.a [nz+kii] = czero;
					nzj++;
					nz += blai*blaj;
					jind2++;
				} else if (jj1 == jj2) {
					jsup = jj1;
					blaj = sprndr[jsup+1]-sprndr[jsup];
					temp.ja[nzj] = jsup;
					ibs = bsa[jind1];
					temp.bsa[nzj] = nz;
					for (kii=0;kii<blai*blaj;kii++) temp.a [nz+kii] = a[ibs+kii];
					nzj++;
					nz += blai*blaj;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				jsup = ja[jind1];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				temp.ja[nzj] = jsup;
				ibs = bsa[jind1];
				temp.bsa[nzj] = nz;
				for (kii=0;kii<blai*blaj;kii++) temp.a [nz+kii] = a[ibs+kii];
				nzj++;
				nz += blai*blaj;
				jind1++;
			} else if (jind2 <= iend2) {
				jsup = _mtr2.ja[jind2];
				blaj = sprndr[jsup+1]-sprndr[jsup];
				temp.ja[nzj] = jsup;
				temp.bsa[nzj] = nz;
				for (kii=0;kii<blai*blaj;kii++) temp.a [nz+kii] = czero;
				nzj++;
				nz += blai*blaj;
				jind2++;
			};
		};
		temp.ia[ilist+1] = nzj;
	};

// Return obtained matrix

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute transposed matrix
//========================================================================================
CSMatrixCS CSMatrixCS::TranspMtr () const { // Compute transposed matrix 

	const char *funcname = "TranspMtr";

	CSMatrixCS at(nsupr,nzja,nza);

	int ilist, isup, blai, blaj;
	int i, j, jj;

	for (i=0;i<=nsupr;i++) at.ia[i] = 0;

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<nsupr;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [nsupr];

	if (!pat) MemoryFail (funcname);

	int k, kii, kjj, ibs, kbs;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			at.ja[k] = i;
			pat[jj]++;
		};
	};

	int nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=at.ia[i];j<at.ia[i+1];j++) {
			jj = at.ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			at.bsa[j] = nz;
			nz += blai*blaj;
		};
	};

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			k = pat[jj];
			ibs = bsa[j];
			kbs = at.bsa[k];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					at.a[kbs+kii*blaj+kjj] = conj(a[ibs+kjj*blai+kii]);
				};
			};
			pat[jj]++;
		};
	};

	delete [] pat;

// Assign control data

	at.m = n;
	at.n = n;
	at.nsupr = nsupr;
	at.nsupc = nsupr;
	at.nlist = nsupr;
	at.nzja = nzja;
	at.nza = nza;
	at.nzatot = nza;
	at.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) at.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) at.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) at.sprndc[isup] = sprndr[isup];

	return at;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute transposed matrix without conjugation
//========================================================================================
CSMatrixCS CSMatrixCS::TranspMtrNoConj () const { // Compute transposed matrix without conjugation

	const char *funcname = "TranspMtrNoConj";

	CSMatrixCS at(nsupr,nzja,nza);

	int ilist, isup, blai, blaj;
	int i, j, jj;

	for (i=0;i<=nsupr;i++) at.ia[i] = 0;

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<nsupr;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [nsupr];

	if (!pat) MemoryFail (funcname);

	int k, kii, kjj, ibs, kbs;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			at.ja[k] = i;
			pat[jj]++;
		};
	};

	int nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=at.ia[i];j<at.ia[i+1];j++) {
			jj = at.ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			at.bsa[j] = nz;
			nz += blai*blaj;
		};
	};

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			k = pat[jj];
			ibs = bsa[j];
			kbs = at.bsa[k];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					at.a[kbs+kii*blaj+kjj] = a[ibs+kjj*blai+kii];
				};
			};
			pat[jj]++;
		};
	};

	delete [] pat;

// Assign control data

	at.m = n;
	at.n = n;
	at.nsupr = nsupr;
	at.nsupc = nsupr;
	at.nlist = nsupr;
	at.nzja = nzja;
	at.nza = nza;
	at.nzatot = nza;
	at.blamx = blamx;

	for (ilist=0;ilist<nsupr;ilist++) at.list[ilist] = ilist;

	for (isup=0;isup<=nsupr;isup++) at.sprndr[isup] = sprndr[isup];
	for (isup=0;isup<=nsupr;isup++) at.sprndc[isup] = sprndr[isup];

	return at;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute transposed for the matrix stored by rows without conjugation
//========================================================================================
CSMatrixCS CSMatrixCS::TranspMtrRowsNoConj () const { // Compute transposed for the matrix stored by rows without conjugation

	const char *funcname = "TranspMtrRowsNoConj";

	int nsupmx = nsupc;
	if (nsupr > nsupmx) nsupmx = nsupr;

	CSMatrixCS at(nsupmx,nzja,nza);

	int ilist, isup, blai, blaj;
	int i, j, jj;

	for (i=0;i<=nsupc;i++) at.ia[i] = 0;

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<nsupc;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [nsupc];

	if (!pat) MemoryFail (funcname);

	int k, kii, kjj, ibs, kbs;

	for (i=0;i<nsupc;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			at.ja[k] = i;
			pat[jj]++;
		};
	};

	int nz = 0;

	for (i=0;i<nsupc;i++) {
		blai = sprndc[i+1]-sprndc[i];
		for (j=at.ia[i];j<at.ia[i+1];j++) {
			jj = at.ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			at.bsa[j] = nz;
			nz += blai*blaj;
		};
	};

	for (i=0;i<nsupc;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			blaj = sprndc[jj+1]-sprndc[jj];
			k = pat[jj];
			ibs = bsa[j];
			kbs = at.bsa[k];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					at.a[kbs+kii*blaj+kjj] = a[ibs+kjj*blai+kii];
				};
			};
			pat[jj]++;
		};
	};

	delete [] pat;

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupr = nsupr;
	at.nsupc = nsupc;
	at.nlist = nsupc;
	at.nzja = nzja;
	at.nza = nza;
	at.nzatot = nza;
	at.blamx = blamx;

	for (ilist=0;ilist<nsupc;ilist++) at.list[ilist] = ilist;

	for (isup=0;isup<=nsupc;isup++) at.sprndc[isup] = sprndc[isup];
	for (isup=0;isup<=nsupr;isup++) at.sprndr[isup] = sprndr[isup];

	return at;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute L and U data for given matrix in parallel
//========================================================================================
void CSMatrixCS::CombineLU (CMPIComm &_comm, int _nblks, int *_blks, int *_blks2cpu, // Compute L and U data for given matrix in parallel
									CSMatrixCS &_mtrl, CSMatrixCS &_mtru) const {

	const char *funcname = "CombineLU";

	int nproc = _comm.GetNproc ();
//	int myid = _comm.GetMyid ();

// Init dummy matrix

	CSMatrixCS mtrdummy;

// Prepare U data

// Count the number of elements

	int nz, nzj, ilist, irow, j, jj, blai, blaj;

	nzj = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		blai = sprndr[irow+1]-sprndr[irow];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			if (jj >= irow) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary U

	CSMatrixCS mtru (nlist,nzj,nz);

	mtru.m = n;
	mtru.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtru.list[ilist] = list[ilist];

	mtru.ia[0] = 0;
	nzj = 0;
	nz = 0;

	int ibs, k;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		blai = sprndr[irow+1]-sprndr[irow];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			if (jj >= irow) {
				mtru.ja[nzj] = jj;
				ibs = bsa[j];
				mtru.bsa[nzj] = nz;
				for (k=0;k<blai*blaj;k++) {
					mtru.a[nz] = a[ibs+k];
					nz++;
				};
				nzj++;
			};
		};
		mtru.ia[ilist+1] = nzj;
	};

	mtru.nlist = nlist;
	mtru.nsupr = nsupr;
	mtru.nsupc = nsupc;
	mtru.nzja = nzj;
	mtru.nza = nz;
	mtru.nzatot = nz;
	mtru.blamx = blamx;

	delete [] mtru.sprndr;
	delete [] mtru.sprndc;

	mtru.sprndr = new int [nsupr+1];
	if (!mtru.sprndr) MemoryFail (funcname);
	mtru.sprndc = new int [nsupc+1];
	if (!mtru.sprndc) MemoryFail (funcname);

	int i;

	for (i=0;i<=nsupr;i++) mtru.sprndr[i] = sprndr[i];
	for (i=0;i<=nsupc;i++) mtru.sprndc[i] = sprndc[i];

// Store U data

	_mtru = mtru;

	mtru = mtrdummy;

// Prepare local L data

	nzj = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		blai = sprndr[irow+1]-sprndr[irow];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			if (jj <= irow) {
				nzj++;
				nz += blai*blaj;
			};
		};
	};

// Allocate and init temporary L

	CSMatrixCS mtrl (nlist,nzj,nz);

	mtrl.m = n;
	mtrl.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtrl.list[ilist] = list[ilist];

	mtrl.ia[0] = 0;
	nzj = 0;
	nz = 0;

	int kii, kjj;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		blai = sprndr[irow+1]-sprndr[irow];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			if (jj <= irow) {
				mtrl.ja[nzj] = jj;
				mtrl.bsa[nzj] = nz;
				ibs = bsa[j];
				for (kii=0;kii<blai;kii++) {
					for (kjj=0;kjj<blaj;kjj++) {
						mtrl.a [nz+kjj*blai+kii] = a[ibs+kjj*blai+kii];
					};
				};
				nzj++;
				nz += blai*blaj;
			};
		};
		mtrl.ia[ilist+1] = nzj;
	};

	mtrl.nlist = nlist;
	mtrl.nsupr = nsupr;
	mtrl.nsupc = nsupc;
	mtrl.nzja = nzj;
	mtrl.nza = nz;
	mtrl.nzatot = nz;
	mtrl.blamx = blamx;

// Mark block number for each sparsity element

	int *jablk;

	jablk = new int [nzj];
	if (!jablk) MemoryFail (funcname);

	int iblkini = 0;
	int iblkprev = 0;

	int ip, iblkbeg, iblkend;

	for (i=0;i<nlist;i++) {
		for (j=mtrl.ia[i];j<mtrl.ia[i+1];j++) {
			jj = mtrl.ja[j];
			if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
				jablk[j] = iblkprev;
				iblkbeg = iblkprev;
				iblkend = iblkprev;
			} else {
				if (j>mtrl.ia[i]) {
					ip = iblkprev+1;
					if (jj >= _blks[ip] && jj < _blks[ip+1]) {
						jablk[j] = ip;
						iblkprev = ip;
						iblkbeg = ip;
						iblkend = ip;
						goto EndCycle;
					};
					iblkbeg = iblkprev+1;
					iblkend = _nblks-1;
				} else {
					ip = iblkini;
					if (jj >= _blks[ip] && jj < _blks[ip+1]) {
						jablk[j] = ip;
						iblkprev = ip;
						iblkbeg = ip;
						iblkend = ip;
						goto EndCycle;
					};
					iblkbeg = 0;
					iblkend = _nblks-1;
				};
				while (iblkbeg != iblkend) {
					ip = (iblkbeg+iblkend)/2;
					if (ip<iblkbeg) ip = iblkbeg;
					if (ip>iblkend) ip = iblkend;
					if (jj >= _blks[ip] && jj < _blks[ip+1]) {
						jablk[j] = ip;
						iblkbeg = ip;
						iblkprev = ip;
						break;
					};
					if (jj < _blks[ip]) {
						iblkend = ip-1;
						if (iblkend < iblkbeg) iblkend = iblkbeg;
					} else if (jj >= _blks[ip+1]) {
						iblkbeg = ip+1;
						if (iblkend < iblkbeg) iblkbeg = iblkend;
					};
				};
				jablk[j] = iblkbeg;
			};
EndCycle:;
			iblkprev = jablk[j];
			if (j == mtrl.ia[i]) iblkini = jablk[j];
		};
	};

// Create sublists for each cpu

	int *iacpu;
	int *nzcpu;
	int *iptrcpu;
	int *jacpuind;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	nzcpu = new int [nproc+1];
	if (!nzcpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);
	jacpuind = new int [nzj];
	if (!jacpuind) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;
	for (i=0;i<nproc;i++) nzcpu[i] = 0;

	int iblk, iproc;

	for (j=0;j<nzj;j++) {
		iblk = jablk[j];
		iproc = _blks2cpu[iblk];
		iacpu[iproc+1]++;
	};

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		blai = sprndr[irow+1]-sprndr[irow];
		for (j=mtrl.ia[ilist];j<mtrl.ia[ilist+1];j++) {
			jj = mtrl.ja[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			iblk = jablk[j];
			iproc = _blks2cpu[iblk];
			nzcpu[iproc] += blai*blaj;
		};
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (j=0;j<nzj;j++) {
		iblk = jablk[j];
		iproc = _blks2cpu[iblk];
		k = iptrcpu[iproc];
		jacpuind[k] = j;
		iptrcpu[iproc]++;
	};

	int *ja2irow;

	ja2irow = new int [nzj];
	if (!ja2irow) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		irow = mtrl.list[i];
		for (j=mtrl.ia[i];j<mtrl.ia[i+1];j++) {
			ja2irow[j] = irow;
		};
	};

// Create the lists of blocks per cpu

	int *imaskblk;
	int *iablkcpu;
	int *jablkcpu;
	int *ibsblk;

	imaskblk = new int [_nblks];
	if (!imaskblk) MemoryFail (funcname);
	iablkcpu = new int [nproc+1];
	if (!iablkcpu) MemoryFail (funcname);
	jablkcpu = new int [_nblks*nproc];
	if (!jablkcpu) MemoryFail (funcname);
	ibsblk = new int [_nblks];
	if (!ibsblk) MemoryFail (funcname);

	int icycleblk = -1;

	for (i=0;i<_nblks;i++) imaskblk[i] = icycleblk;
	for (i=0;i<nproc+1;i++) iablkcpu[i] = 0;

	int nlistloc, jind;

	for (iproc=0;iproc<nproc;iproc++) {
		icycleblk++;
		nlistloc = 0;
		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			iblk = jablk[jind];
			if (imaskblk[iblk] != icycleblk) {
				nlistloc++;
				imaskblk[iblk] = icycleblk;
			};
		};
		iablkcpu[iproc+1] = nlistloc;
	};

	for (i=0;i<nproc;i++) iablkcpu[i+1] = iablkcpu[i]+iablkcpu[i+1];

	int ishft;

	for (iproc=0;iproc<nproc;iproc++) {
		ishft = iablkcpu[iproc];
		icycleblk++;
		nlistloc = 0;
		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			iblk = jablk[jind];
			if (imaskblk[iblk] != icycleblk) {
				jablkcpu[ishft+nlistloc] = iblk;
				nlistloc++;
				imaskblk[iblk] = icycleblk;
			};
		};
		qsort (jablkcpu+ishft,nlistloc,sizeof(int),compint);
	};

	int nlistmax = 0, ni;

	for (iproc=0;iproc<nproc;iproc++) {
		nlistloc = 0;
		for (j=iablkcpu[iproc];j<iablkcpu[iproc+1];j++) {
			iblk = jablkcpu[j];
			ibsblk[j] = nlistloc;
			ni = _blks[iblk+1]-_blks[iblk];
			nlistloc += ni;
		};
		if (nlistloc > nlistmax) nlistmax = nlistloc;
	};

// Create transposed L structures for each cpu

	CSMatrixCS *mtrarr;

	mtrarr = new CSMatrixCS [nproc];
	if (!mtrarr) MemoryFail (funcname);

	int *imask;
	int *listloc;
	int *ibslist;
	int *iptr;

	imask = new int [nlistmax];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nlistmax];
	if (!listloc) MemoryFail (funcname);
	ibslist = new int [nlistmax];
	if (!ibslist) MemoryFail (funcname);
	iptr = new int [nlistmax];
	if (!iptr) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlistmax;i++) imask[i] = icycle;

	int jjloc, nzloc, ipos, nzjloc, ibstemp, ibsl;

	for (iproc=0;iproc<nproc;iproc++) {

		icycle++;

		for (j=iablkcpu[iproc];j<iablkcpu[iproc+1];j++) {
			iblk = jablkcpu[j];
			imaskblk[iblk] = j;
		};

		nlistloc = 0;
		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			jj = mtrl.ja[jind];
			iblk = jablk[jind];
			ipos = imaskblk[iblk];
			jjloc = jj-_blks[iblk];
			ibs = ibsblk[ipos];
			if (imask[ibs+jjloc] != icycle) {
				listloc[nlistloc] = jj;
				ibslist[ibs+jjloc] = nlistloc;
				nlistloc++;
				imask[ibs+jjloc] = icycle;
			};
		};

		nzjloc = iacpu[iproc+1]-iacpu[iproc];
		nzloc = nzcpu[iproc];

		CSMatrixCS temp (nlistloc,nzjloc, nzloc);

		for (i=0;i<nlistloc;i++) temp.list[i] = listloc[i];

		for (i=0;i<=nlistloc;i++) temp.ia[i] = 0;

		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			jj = mtrl.ja[jind];
			irow = ja2irow[jind];
			blai = sprndr[irow+1]-sprndr[irow];
			blaj = sprndr[jj+1]-sprndr[jj];
			iblk = jablk[jind];
			ipos = imaskblk[iblk];
			jjloc = jj-_blks[iblk];
			ibs = ibsblk[ipos];
			ilist = ibslist[ibs+jjloc];
			temp.ia[ilist+1]++;
		};

		for (i=0;i<nlistloc;i++) temp.ia[i+1] = temp.ia[i]+temp.ia[i+1];

		for (i=0;i<nlistloc;i++) iptr[i] = temp.ia[i];

		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			jj = mtrl.ja[jind];
			iblk = jablk[jind];
			ipos = imaskblk[iblk];
			jjloc = jj-_blks[iblk];
			ibs = ibsblk[ipos];
			ilist = ibslist[ibs+jjloc];
			k = iptr[ilist];
			temp.ja[k] = ja2irow[jind];
			iptr[ilist]++;
		};

		nzloc = 0;

		for (i=0;i<nlistloc;i++) {
			irow = temp.list[i];
			for (j=temp.ia[i];j<temp.ia[i+1];j++) {
				jj = temp.ja[j];
				blai = sprndr[irow+1]-sprndr[irow];
				blaj = sprndr[jj+1]-sprndr[jj];
				temp.bsa[j] = nzloc;
				nzloc += blai*blaj;
			};
		};

		for (i=0;i<nlistloc;i++) iptr[i] = temp.ia[i];

		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			jj = mtrl.ja[jind];
			iblk = jablk[jind];
			ipos = imaskblk[iblk];
			jjloc = jj-_blks[iblk];
			ibs = ibsblk[ipos];
			ilist = ibslist[ibs+jjloc];
			k = iptr[ilist];
			irow = ja2irow[jind];
			ibstemp = temp.bsa[k];
			ibsl = mtrl.bsa[jind];
			blai = sprndr[irow+1]-sprndr[irow];
			blaj = sprndr[jj+1]-sprndr[jj];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					temp.a[ibstemp+kii*blaj+kjj] = mtrl.a[ibsl+kjj*blai+kii];
				};
			};
			iptr[ilist]++;
		};

		delete [] temp.sprndr;
		delete [] temp.sprndc;

		temp.sprndr = new int [nsupr+1];
		if (!temp.sprndr) MemoryFail (funcname);
		temp.sprndc = new int [nsupc+1];
		if (!temp.sprndc) MemoryFail (funcname);

		for (i=0;i<=nsupr;i++) temp.sprndr[i] = sprndr[i];
		for (i=0;i<=nsupc;i++) temp.sprndc[i] = sprndc[i];

		temp.m = n;
		temp.n = n;
		temp.nsupr = nsupr;
		temp.nsupc = nsupc;
		temp.nlist = nlistloc;
		temp.nzja = nzjloc;
		temp.nza = nzloc;
		temp.nzatot = nzloc;
		temp.blamx = blamx;

		temp.SortListAndColumns ();

		mtrarr[iproc] = temp;

		temp = mtrdummy;

	};

// Echange the matrices

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

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		mtrarr[i].PackMatrix (ObjSizeSend[i],ObjSend[i]);
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
	if(info) throw " CSMatrixCS::CombineLU: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrixRS::CombineLU: wrong number of the received objects ";

	for(i = 0; i < NObjRecv; i++) mtrarr[i].UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect submatrices of the current cpu

	int *listblkloc;

	listblkloc = new int [nproc];
	if (!listblkloc) MemoryFail (funcname);

	for (i=0;i<nproc;i++) listblkloc[i] = i;

	CGSMatrixCS::CombineMatricesSort (n, nproc, listblkloc,
												nsupr, sprndr,
												mtrarr,
												_mtrl);

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

	_mtrl.RemoveEquivalentRowsNumbers ();

// Free work arrays

	delete [] jablk;
	delete [] iacpu;
	delete [] nzcpu;
	delete [] iptrcpu;
	delete [] jacpuind;
	delete [] ja2irow;
	delete [] imaskblk;
	delete [] iablkcpu;
	delete [] jablkcpu;
	delete [] ibsblk;
	delete [] mtrarr;
	delete [] imask;
	delete [] listloc;
	delete [] ibslist;
	delete [] iptr;
	delete [] listblkloc;

// Extend the sparsity structures of L and U 

	_mtrl.ExtendSparsity (_mtru);
	_mtru.ExtendSparsity (_mtrl);

};

// Author: Kharchenko S.A.
// CSMatrixCS: Remove from the list of rows equivalent rows numbers
//========================================================================================
void CSMatrixCS::RemoveEquivalentRowsNumbers () { // Remove from the list of rows equivalent rows numbers

	const char *funcname = "RemoveEquivalentRowsNumbers";

	if (nlist == 0) return;

	CSMatrixCS mtrnew = *this;

// Modify the rows list and ia array

	int ilistnew = 0;
	int irownew = list[0];

	mtrnew.list[ilistnew] = irownew;

	int i, irow;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		if (irow != irownew) {
			ilistnew++;
			irownew = irow;
			mtrnew.list[ilistnew] = irownew;
		};
		mtrnew.ia[ilistnew+1] = ia[i+1];
	};

	int nlistnew = ilistnew+1;

// Sort data in each row

	int nzmax = 0;
	int nz;
	int nzrowmax = 0;
	int nzrow;

	int blai, j, jj, blaj;

	for (i=0;i<nlistnew;i++) {
		irow = mtrnew.list[i];
		blai = mtrnew.sprndr[irow+1]-mtrnew.sprndr[irow];
		nz = mtrnew.ia[i+1]-mtrnew.ia[i];
		if (nz > nzmax) nzmax = nz;
		nzrow = 0;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			jj = mtrnew.ja[j];
			blaj = mtrnew.sprndr[jj+1]-mtrnew.sprndr[jj];
			nz++;
			nzrow += blai*blaj;
		};
		if (nzrowmax < nzrow) nzrowmax = nzrow;
	};

	CIntInt *iiarr;
	int *iarrloc;
	int *iarr2loc;
	dcmplx *darrloc;

	iiarr = new CIntInt [nzmax];
	if(!iiarr) MemoryFail(funcname);
	iarrloc = new int [nzmax];
	if(!iarrloc) MemoryFail(funcname);
	iarr2loc = new int [nzmax];
	if(!iarr2loc) MemoryFail(funcname);
	darrloc = new dcmplx [nzrowmax];
	if(!darrloc) MemoryFail(funcname);

	int nztot = 0;

	int ibs, ind, kii;

	for (i=0;i<nlistnew;i++) {
		nz = 0;
		irow = mtrnew.list[i];
		blai = mtrnew.sprndr[irow+1]-mtrnew.sprndr[irow];
		nzrow = 0;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			jj = mtrnew.ja[j];
			blaj = mtrnew.sprndr[jj+1]-mtrnew.sprndr[jj];
			iiarr[nz].intvalue = mtrnew.ja[j];
			iiarr[nz].int2value = j-mtrnew.ia[i];
			nz++;
			nzrow += blai*blaj;
		};
		qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
		nz = 0;
		nzrow = 0;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			ind = iiarr[nz].int2value;
			jj = mtrnew.ja[mtrnew.ia[i]+ind];
			ibs = mtrnew.bsa[mtrnew.ia[i]+ind];
			blaj = mtrnew.sprndr[jj+1]-mtrnew.sprndr[jj];
			iarrloc[j-mtrnew.ia[i]] = jj;
			iarr2loc[j-mtrnew.ia[i]] = nztot+nzrow;
			for (kii=0;kii<blai*blaj;kii++) darrloc[nzrow+kii] = mtrnew.a[ibs+kii];
			nz++;
			nzrow += blai*blaj;
		};
		for (j=0;j<nz;j++) mtrnew.ja[mtrnew.ia[i]+j] = iarrloc[j];
		for (j=0;j<nz;j++) mtrnew.bsa[mtrnew.ia[i]+j] = iarr2loc[j];
		for (j=0;j<nzrow;j++) mtrnew.a[nztot+j] = darrloc[j];
		nztot += nzrow;
	};

	delete [] iiarr;
	delete [] iarrloc;
	delete [] iarr2loc;
	delete [] darrloc;

	mtrnew.SetNlist (nlistnew);

// Return the result

	*this = mtrnew;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute symmetrized matrix
//========================================================================================
CSMatrixCS CSMatrixCS::SymmMtr () const { // Compute symmetrized matrix 

//	const char *funcname = "SymmMtr";

// Compute transposed matrix

	CSMatrixCS at;

	at = TranspMtr();

// Compute the symmetric matrix

	int nzjas, nzas;
	nzjas = 2*nzja;
	nzas = 2*nza;

	CSMatrixCS as(nsupr,nzjas,nzas);

	int i, nz, nzj;
	int ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;
	int blai, blaj, ibs, kii, kjj;

	as.ia[0] = 0;

	nzj = 0;
	nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		ibeg1 = ia[i];
		iend1 = ia[i+1]-1;
		ibeg2 = at.ia[i];
		iend2 = at.ia[i+1]-1;
		jind1 = ibeg1;
		jind2 = ibeg2;
		while (jind1 <= iend1 || jind2 <= iend2) {
			if (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = at.ja[jind2];
				if (jj1 < jj2) {
					blaj = sprndr[jj1+1]-sprndr[jj1];
					as.ja[nzj] = ja[jind1];
					as.bsa[nzj] = nz;
					ibs = bsa[jind1];
					for (kii=0;kii<blai*blaj;kii++) as.a[nz+kii] = a[ibs+kii];
					nzj++;
					nz += blai*blaj;
					jind1++;
				} else if (jj1 > jj2) {
					blaj = sprndr[jj2+1]-sprndr[jj2];
					as.ja[nzj] = at.ja[jind2];
					as.bsa[nzj] = nz;
					ibs = at.bsa[jind2];
					for (kii=0;kii<blai;kii++) {
						for (kjj=0;kjj<blaj;kjj++) {
							as.a[nz+kjj*blai+kii] = at.a[ibs+kjj*blai+kii];
						};
					};
					nzj++;
					nz += blai*blaj;
					jind2++;
				} else if (jj1 == jj2) {
					blaj = sprndr[jj1+1]-sprndr[jj1];
					as.ja[nzj] = ja[jind1];
					as.bsa[nzj] = nz;
					ibs = bsa[jind1];
					for (kii=0;kii<blai*blaj;kii++) as.a[nz+kii] = a[ibs+kii];
					nzj++;
					nz += blai*blaj;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				jj1 = ja[jind1];
				blaj = sprndr[jj1+1]-sprndr[jj1];
				as.ja[nzj] = ja[jind1];
				as.bsa[nzj] = nz;
				ibs = bsa[jind1];
				for (kii=0;kii<blai*blaj;kii++) as.a[nz+kii] = a[ibs+kii];
				nzj++;
				nz += blai*blaj;
				jind1++;
			} else if (jind2 <= iend2) {
				jj2 = at.ja[jind2];
				blaj = sprndr[jj2+1]-sprndr[jj2];
				as.ja[nzj] = at.ja[jind2];
				as.bsa[nzj] = nz;
				ibs = at.bsa[jind2];
				for (kii=0;kii<blai;kii++) {
					for (kjj=0;kjj<blaj;kjj++) {
						as.a[nz+kjj*blai+kii] = at.a[ibs+kjj*blai+kii];
					};
				};
				nzj++;
				nz += blai*blaj;
				jind2++;
			};
		};
		as.ia[i+1] = nzj;
	};

	for (i=0;i<nsupr;i++) as.list[i] = i;
	for (i=0;i<=nsupr;i++) as.sprndr[i] = sprndr[i];
	for (i=0;i<=nsupr;i++) as.sprndc[i] = sprndr[i];

	as.m = n;
	as.n = n;
	as.nsupr = nsupr;
	as.nsupc = nsupr;
	as.nlist = nsupr;
	as.nzja = nzj;
	as.nza = nz;
	as.nzatot = nz;
	as.blamx = blamx;

	return as;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Ordering of the matrix columns
//========================================================================================
CSMatrixCS CSMatrixCS::OrdMtrRectCol (const int *_order) const { // Ordering of the matrix columns

	const char *funcname = "OrdMtrRectCol";

// Compute inverse order

	int *iord;
	iord = new int [nsupc];
	if (!iord) MemoryFail (funcname);

	int i;
	for (i=0;i<nsupc;i++) iord[_order[i]] = i;

// Reorder the matrix

	int nsupmx = nsupr;

	if (nsupc > nsupmx) nsupmx = nsupc;

	CSMatrixCS an(nsupmx,nzja,nza);

// Copy list

	for (i=0;i<nsupc;i++) an.list[i] = list[i];

// Copy sprndr

	for (i=0;i<=nsupr;i++) an.sprndr[i] = sprndr[i];

// Create new sprndc

	for (i=0;i<=nsupc;i++) an.sprndc[i] = 0;

	int j;

	for (i=0;i<nsupc;i++) {
		j = _order[i];
		an.sprndc[j+1] = sprndc[i+1]-sprndc[i];
	};

	for (i=0;i<nsupc;i++) an.sprndc[i+1] += an.sprndc[i];

// Create new ia

	for (i=0;i<=nsupc;i++) an.ia[i] = 0;

	for (i=0;i<nsupc;i++) {
		j = _order[i];
		an.ia[j+1] = ia[i+1]-ia[i];
	};

	for (i=0;i<nsupc;i++) an.ia[i+1] += an.ia[i];

// Create new ja

	int nz = 0;
	int isup;

	int isupnew;
	for (isupnew=0;isupnew<nsupc;isupnew++) {
		isup = iord[isupnew];
		for (j=ia[isup];j<ia[isup+1];j++) {
			an.ja[nz] = ja[j];
			nz++;
		};
	};

// Compute new bsa

	nz = 0;

	int blai,blaj,jsup,kii,ibs;

	for (isup=0;isup<nsupc;isup++) {
		blai = an.sprndc[isup+1]-an.sprndc[isup];
		for (j=an.ia[isup];j<an.ia[isup+1];j++) {
			jsup = an.ja[j];
			blaj = an.sprndr[jsup+1]-an.sprndr[jsup];
			an.bsa[j] = nz;
			nz += blai * blaj;
		};
	};

// Compute new a

	nz = 0;

	for (isupnew=0;isupnew<nsupc;isupnew++) {
		blai = an.sprndc[isupnew+1]-an.sprndc[isupnew];
		isup = iord[isupnew];
		for (j=ia[isup];j<ia[isup+1];j++) {
			jsup = ja[j];
			blaj = sprndr[jsup+1]-sprndr[jsup];
			ibs = bsa[j];
			for (kii=0;kii<blai*blaj;kii++) an.a[nz+kii] = a[ibs+kii];
			nz += blai*blaj;
		};
	};

// Init head data

	an.m      = m;
	an.n      = n;
	an.nsupc  = nsupc;
	an.nsupr  = nsupr;
	an.nlist  = nlist;
	an.nzja   = nzja;
	an.nza    = nza;
	an.nzatot = nzatot;
	an.blamx  = blamx;

// Delete ordering array

	delete [] iord;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Symmetric matrix ordering routine
//========================================================================================
CSMatrixCS CSMatrixCS::OrdMtrSymm (const int *_order) const { // Symmetric matrix ordering routine

//	const char *funcname = "OrdMtrSymm";

// Compute the symmetric matrix

	CSMatrixCS as;

	as = SymmMtr ();

// Reorder the matrix

	CSMatrixCS an;

	an = as.OrdMtr (_order);

	CSMatrixCS mtrdummy;

	as = mtrdummy;

// Filter unnecessary elements

	CSMatrixCS af(nsupr,nzja,nza);

	int i, j, jj, blai, blaj, kii, ibs;

	af.ia[0] = 0;

	int nzj = 0;
	int nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = an.sprndr[i+1]-an.sprndr[i];
		for (j=an.ia[i];j<an.ia[i+1];j++) {
			jj = an.ja[j];
			blaj = an.sprndr[jj+1]-an.sprndr[jj];
			if (jj >= i) {
				af.ja[nzj] = an.ja[j];
				af.bsa[nzj] = nz;
				ibs = an.bsa[j];
				for (kii=0;kii<blai*blaj;kii++) af.a[nz+kii] = an.a[ibs+kii];
				nzj++;
				nz += blai*blaj;
			};
		};
		af.ia[i+1] = nzj;
	};

	af.m      = m;
	af.n      = n;
	af.nsupc  = nsupc;
	af.nsupr  = nsupr;
	af.nlist  = nlist;
	af.nzja   = nzja;
	af.nza    = nza;
	af.nzatot = nzatot;
	af.blamx  = blamx;

	for (i=0;i<=nsupr;i++) af.sprndr[i] = an.sprndr[i];
	for (i=0;i<=nsupr;i++) af.sprndc[i] = an.sprndr[i];

	return af;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Output matrix
//========================================================================================
void CSMatrixCS::OutputMatrix (ostream _stream, // Output matrix
								int _nblks, int *_blks, FILE **_mtrfiles) const {

// Output class CSMatrixCS

	_stream << *this << endl;

// Main read cycle

	dcmplx *aloc;

	aloc = new dcmplx [blamx*blamx];

	for (int iblk=0;iblk<_nblks;iblk++) {

// Read data of the current block

		for (int isup=_blks[iblk];isup<_blks[iblk+1];isup++) {
			int blai=sprndc[isup+1]-sprndc[isup];
			for (int j=ia[isup];j<ia[isup+1];j++) {
				int jsup = ja[j];
				int blaj=sprndr[jsup+1]-sprndr[jsup];
				int ibs = bsa[j];
				FGet (_mtrfiles[iblk], blai*blaj, aloc, ibs);
				OutArr(_stream," Local blk = ",blai*blaj,aloc);
			};
		};

	};

	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Pack matrix data
//========================================================================================
void CSMatrixCS::PackMatrix (int &_length, char *&_obj) { // Pack matrix data

	const char* funcname = "PackMatrix";

// Count the total size of the object

	_length = 12 * sizeof(int) + (2*nlist+1+nlist2+nzja+nzja2+nzja3) * sizeof(int) + (nza) * sizeof(dcmplx);
	if (blamx > 1) {
		_length += (nlist*2+nzja+2)*sizeof(int);
	};

// Allocate the structure with the head and init pointers

	_obj = new char [_length];
	if (!_obj) MemoryFail(funcname);

	char* pLoc;

	pLoc = _obj;

	int *pHead;
	int *plist;
	int *plist2;
	int *pia;
	int *pja;
	int *pja2;
	int *pja3;
	int *psprndr = 0;
	int *psprndc = 0;
	int *pbsa = 0;
	dcmplx *pa;

	pHead = (int *) pLoc;
	pLoc += 12 * sizeof(int);

	plist = (int *) pLoc;
	pLoc += nlist * sizeof(int);

	plist2 = (int *) pLoc;
	pLoc += nlist2 * sizeof(int);

	pia = (int *) pLoc;
	pLoc += (nlist+1) * sizeof(int);

	pja = (int *) pLoc;
	pLoc += (nzja) * sizeof(int);

	pja2 = (int *) pLoc;
	pLoc += (nzja2) * sizeof(int);

	pja3 = (int *) pLoc;
	pLoc += (nzja3) * sizeof(int);

	if (blamx > 1) {
		psprndr = (int *) pLoc;
		pLoc += (nlist+1) * sizeof(int);
		psprndc = (int *) pLoc;
		pLoc += (nlist+1) * sizeof(int);
		pbsa = (int *) pLoc;
		pLoc += (nzja) * sizeof(int);
	};

	pa = (dcmplx *) pLoc;
	pLoc += nza * sizeof(dcmplx);

	pHead[0]  = m;
	pHead[1]  = n;
	pHead[2]  = nsupr;
	pHead[3]  = nsupc;
	pHead[4]  = nlist;
	pHead[5]  = nlist2;
	pHead[6]  = nzja;
	pHead[7]  = nzja2;
	pHead[8]  = nzja3;
	pHead[9]  = nzatot;
	pHead[10] = nza;
	pHead[11] = blamx;

// Fill arrays

	int j;

	for(j = 0; j < nlist;   j++) plist[j] = list[j];
	for(j = 0; j < nlist2;  j++) plist2[j] = list2[j];
	for(j = 0; j < nlist+1; j++) pia[j] = ia[j];
	for(j = 0; j < nzja;    j++) pja[j] = ja[j];
	for(j = 0; j < nzja2;   j++) pja2[j] = ja2[j];
	for(j = 0; j < nzja3;   j++) pja3[j] = ja3[j];
	if (blamx > 1) {
		for(j = 0; j < nlist+1; j++) psprndr[j] = sprndr[j];
		for(j = 0; j < nlist+1; j++) psprndc[j] = sprndc[j];
		for(j = 0; j < nzja;    j++) pbsa[j] = bsa[j];
	};
	for(j = 0; j < nza;     j++) pa[j] = a[j];

};

// Author: Kharchenko S.A.
// CSMatrixCS: UnPack matrix data
//========================================================================================
void CSMatrixCS::UnPackMatrix (int _length, char *_obj) { // UnPack matrix data

	const char* funcname = "UnPackMatrix";

// Free previous arrays

	delete [] list;
	delete [] list2;
	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;
	delete [] sprndr;
	delete [] sprndc;
	delete [] bsa;
	delete [] a;

// Get head data

	char* pLoc;

	pLoc = _obj;

	int* pHead;
	int *plist;
	int *plist2;
	int *pia;
	int *pja;
	int *pja2;
	int *pja3;
	int *psprndr = 0;
	int *psprndc = 0;
	int *pbsa = 0;
	dcmplx *pa;

	pHead = (int *) pLoc;
	pLoc += 12 * sizeof(int);

	m       = pHead[0];
	n       = pHead[1];
	nsupr   = pHead[2];
	nsupc   = pHead[3];
	nlist   = pHead[4];
	nlist2  = pHead[5];
	nzja    = pHead[6];
	nzja2   = pHead[7];
	nzja3   = pHead[8];
	nzatot  = pHead[9];
	nza     = pHead[10];
	blamx   = pHead[11];

	plist = (int *) pLoc;
	pLoc += nlist * sizeof(int);

	plist2 = (int *) pLoc;
	pLoc += nlist2 * sizeof(int);

	pia = (int *) pLoc;
	pLoc += (nlist+1) * sizeof(int);

	pja = (int *) pLoc;
	pLoc += (nzja) * sizeof(int);

	pja2 = (int *) pLoc;
	pLoc += (nzja2) * sizeof(int);

	pja3 = (int *) pLoc;
	pLoc += (nzja3) * sizeof(int);

	if (blamx > 1) {
		psprndr = (int *) pLoc;
		pLoc += (nlist+1) * sizeof(int);
		psprndc = (int *) pLoc;
		pLoc += (nlist+1) * sizeof(int);
		pbsa = (int *) pLoc;
		pLoc += (nzja) * sizeof(int);
	};

	pa = (dcmplx *) pLoc;
	pLoc += nza * sizeof(dcmplx);

// Allocate and get the arrays

	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	if (blamx > 1) {
		sprndr = new int [nlist+1];
		if (!sprndr) MemoryFail (funcname);
		sprndc = new int [nlist+1];
		if (!sprndc) MemoryFail (funcname);
		bsa = new int [nzja];
		if (!bsa) MemoryFail (funcname);
	} else {
		sprndr = new int [0];
		if (!sprndr) MemoryFail (funcname);
		sprndc = new int [0];
		if (!sprndc) MemoryFail (funcname);
		bsa = new int [0];
		if (!bsa) MemoryFail (funcname);
	};
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);

	int j;

	for(j = 0; j < nlist;   j++) list[j] = plist[j];
	for(j = 0; j < nlist2;  j++) list2[j] = plist2[j];
	for(j = 0; j < nlist+1; j++) ia[j] = pia[j];
	for(j = 0; j < nzja;    j++) ja[j] = pja[j];
	for(j = 0; j < nzja2;   j++) ja2[j] = pja2[j];
	for(j = 0; j < nzja3;   j++) ja3[j] = pja3[j];
	if (blamx > 1) {
		for(j = 0; j < nlist+1; j++) sprndr[j] = psprndr[j];
		for(j = 0; j < nlist+1; j++) sprndc[j] = psprndc[j];
		for(j = 0; j < nzja;    j++) bsa[j] = pbsa[j];
	};
	for(j = 0; j < nza;     j++) a[j] = pa[j];

};

// Author: Kharchenko S.A.
// Description: Sort list indices and column indices
// CSMatrixCS::SortListAndColumns()
//========================================================================================
void CSMatrixCS::SortListAndColumns () { // Sort list indices and column indices

	const char *funcname = "SortListAndColumns";

// Reorder the rows

	CIntInt *iiarr;

	iiarr = new CIntInt [nlist];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlist;i++) {
		iiarr[i].intvalue = list[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr,nlist,sizeof(CIntInt),CIntInt::CompareIntInt);

	int *order, *iorder;

	order = new int [nlist];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlist];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		iorder[i] = iiarr[i].int2value;
	};

	for (i=0;i<nlist;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *ialoc;
	int *jaloc;
	int *bsaloc;
	dcmplx *aloc;

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);
	bsaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);
	aloc = new dcmplx [nza];
	if (!aloc) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		j = order[i];
		listloc[j] = list[i];
		ialoc[j+1] = ia[i+1]-ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlist;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	int ni, ibeg;

	for (i=0;i<nlist;i++) {
		j = order[i];
		ibeg = ialoc[j];
		ni = ialoc[j+1]-ialoc[j];
		for (k=ia[i];k<ia[i+1];k++) {
			jaloc[k-ia[i]+ialoc[j]] = ja[k];
		};
	};

	int nzloc = 0;
	int nzrowmax = 0;
	int nzrow;
	int irow, blai, blaj, ind, jj;

	for (i=0;i<nlist;i++) {
		irow = listloc[i];
		blai = sprndr[irow+1]-sprndr[irow];
		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			bsaloc[j] = nzloc;
			nzloc += blai*blaj;
			nzrow += blai*blaj;
		};
		if (nzrow > nzrowmax) nzrowmax = nzrow;
	};

	int ibsold, ibsnew, kii;

	for (i=0;i<nlist;i++) {
		j = order[i];
		irow = listloc[j];
		blai = sprndr[irow+1]-sprndr[irow];
		ibeg = ialoc[j];
		ni = ialoc[j+1]-ialoc[j];
		for (k=ia[i];k<ia[i+1];k++) {
			jj = ja[k];
			blaj = sprndr[jj+1]-sprndr[jj];
			ind = k-ia[i]+ialoc[j];
			ibsold = bsa[k];
			ibsnew = bsaloc[ind];
			for (kii=0;kii<blai*blaj;kii++) aloc[ibsnew+kii] = a[ibsold+kii];
		};
	};

	int nzmax = 0;
	int nz;

	for (i=0;i<nlist;i++) {
		nz = ialoc[i+1]-ialoc[i];
		if (nz > nzmax) nzmax = nz;
	};

	int *iarrloc;
	int *iarr2loc;
	dcmplx *darrloc;

	iiarr = new CIntInt [nzmax];
	if(!iiarr) MemoryFail(funcname);
	iarrloc = new int [nzmax];
	if(!iarrloc) MemoryFail(funcname);
	iarr2loc = new int [nzmax];
	if(!iarr2loc) MemoryFail(funcname);
	darrloc = new dcmplx [nzrowmax];
	if(!darrloc) MemoryFail(funcname);

	int nztot = 0;

	int ibs;

	for (i=0;i<nlist;i++) {
		nz = 0;
		irow = listloc[i];
		blai = sprndr[irow+1]-sprndr[irow];
		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			iiarr[nz].intvalue = jaloc[j];
			iiarr[nz].int2value = j-ialoc[i];
			nz++;
			nzrow += blai*blaj;
		};
		qsort (iiarr, nz, sizeof(CIntInt), CIntInt::CompareIntInt);
		nz = 0;
		nzrow = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			ind = iiarr[nz].int2value;
			jj = jaloc[ialoc[i]+ind];
			ibs = bsaloc[ialoc[i]+ind];
			blaj = sprndr[jj+1]-sprndr[jj];
			iarrloc[j-ialoc[i]] = jj;
			iarr2loc[j-ialoc[i]] = nztot+nzrow;
			for (kii=0;kii<blai*blaj;kii++) darrloc[nzrow+kii] = aloc[ibs+kii];
			nz++;
			nzrow += blai*blaj;
		};
		for (j=0;j<nz;j++) jaloc[ialoc[i]+j] = iarrloc[j];
		for (j=0;j<nz;j++) bsaloc[ialoc[i]+j] = iarr2loc[j];
		for (j=0;j<nzrow;j++) aloc[nztot+j] = darrloc[j];
		nztot += nzrow;
	};

	delete [] iiarr;
	delete [] iarrloc;
	delete [] iarr2loc;
	delete [] darrloc;

	for (i=0;i<nlist;i++) list[i] = listloc[i];
	for (i=0;i<=nlist;i++) ia[i] = ialoc[i];
	for (i=0;i<nzja;i++) ja[i] = jaloc[i];
	for (i=0;i<nzja;i++) bsa[i] = bsaloc[i];
	for (i=0;i<nza;i++) a[i] = aloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] bsaloc;
	delete [] aloc;

	delete [] order;
	delete [] iorder;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CSMatrixCS &_a) { // Output matrix

	_stream << " CSMatrixCS:" << endl;

	_stream << " Blamx = " << _a.blamx << endl;
	OutArr (_stream, " SprndR =", _a.nsupr+1, _a.sprndr);
	OutArr (_stream, " SprndC =", _a.nsupc+1, _a.sprndc);
	OutArr (_stream, " Bsa =", _a.nzja, _a.bsa);

	_stream << (const CSMatrixC &) _a << endl;

	return _stream;

};
