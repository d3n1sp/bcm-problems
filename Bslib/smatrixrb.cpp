#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "smatrix.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixRB: Memory allocation zero data constructor
//========================================================================================
CSMatrixRB::CSMatrixRB (): CSMatrixR () { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixRB_00";

	bla = 0;
	bla_2 = 0;
	bla_3 = 0;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Memory allocation zero data constructor
//========================================================================================
CSMatrixRB::CSMatrixRB (int _nsupa, int _nzja, int _bla): CSMatrixR (_nsupa,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixRB_01";

	delete [] a;

	nza = _nzja*_bla*_bla;
	nzatot = nza;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = 0.0e0;
	};
	bla = _bla;
	bla_2 = bla*bla;
	bla_3 = bla_2*bla;
	m = _nsupa*bla;
	n = _nsupa*bla;
};

// Author: Kharchenko S.A.
// CSMatrixRB: Constructor for Laplace
//========================================================================================
CSMatrixRB::CSMatrixRB (int _nn, int _bla): CSMatrixR () { // Constructor for Laplace

	const char *funcname = "CSMatrixRB_02";

// Create temporary object

	CSMatrixR tempr (_nn);

// Convert into bla format

	CSMatrixRB temp;

	temp = R2RBSymm (tempr, _bla);

// Reassign the data

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Constructor for Laplace^2 with diagonal shift
//========================================================================================
CSMatrixRB::CSMatrixRB (int _nn, double _dshft, int _bla): CSMatrixR () { // Constructor for Laplace^2 with diagonal shift

	const char *funcname = "CSMatrixRB_03";

// Create temporary object

	CSMatrixR tempr (_nn, _dshft);

// Convert into bla format

	CSMatrixRB temp;

	temp = R2RBSymm (tempr, _bla);

// Reassign the data

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Copy constructor
//========================================================================================
CSMatrixRB::CSMatrixRB (const CSMatrixRB &_aa): CSMatrixR (_aa) { // Copy constructor

	const char *funcname = "CSMatrixRB_copy";

	bla = _aa.bla;
	bla_2 = _aa.bla_2;
	bla_3 = _aa.bla_3;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Equality operator
//========================================================================================
CSMatrixRB &CSMatrixRB::operator= (const CSMatrixRB &_aa) { // Equality operator

	const char *funcname = "CSMatrixRB_=";

	this->CSMatrixR::operator= ((const CSMatrixR &) _aa);

	bla = _aa.bla;
	bla_2 = _aa.bla_2;
	bla_3 = _aa.bla_3;

	return *this;
};

// Author: Kharchenko S.A.
// CSMatrixRB: Transform matrix from point to bla format
//========================================================================================
CSMatrixRB CSMatrixRB::R2RBSymm (const CSMatrixR &_ar, int _bla) { // Transform matrix from point to bla format

	const char *funcname = "R2RBSymm";

// Compute nd2sp and sp2nd arrays

	int nloc = _ar.n;
	int nsupa = (nloc+_bla-1) / _bla;

	int *nd2sp, *sp2nd;
	int i, j;

	nd2sp = new int [nloc];
	if (!nd2sp) MemoryFail (funcname);
	sp2nd = new int [nsupa+1];
	if (!sp2nd) MemoryFail (funcname);

	sp2nd[0] = 0;
	for (i=0;i<nsupa-1;i++) {
		sp2nd[i+1] = sp2nd[i] + _bla;
	};
	sp2nd[nsupa] = nloc;

	for (i=0;i<nsupa;i++) {
		for (j=sp2nd[i];j<sp2nd[i+1] && j<nloc;j++) {
			nd2sp[j] = i;
		};
	};

// Count the total number of elements in the collapsed matrix

	int isup, jsup, jj;

	int icycle=0;

	int *imask, *lstloc;

	imask = new int [nsupa];
	if (!imask) MemoryFail (funcname);
	lstloc = new int [nsupa];
	if (!lstloc) MemoryFail (funcname);

	for (i=0;i<nsupa;i++) imask[i] = icycle;

	int nz = 0;
	int nlstloc;

	for (isup=0;isup<nsupa;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			for (j=_ar.ia[i];j<_ar.ia[i+1];j++) {
				jj = _ar.ja[j];
				jsup = nd2sp[jj];
				if (imask[jsup] != icycle) {
					lstloc[nlstloc] = jsup;
					imask[jsup] = icycle;
					nlstloc++;
				};
			};
		};
		nz += nlstloc;
	};

// Compute sparsity structure of the collapsed matrix

	CSMatrixRB temp (nsupa, nz, _bla);

	temp.ia[0] = 0;

	for (isup=0;isup<nsupa;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			for (j=_ar.ia[i];j<_ar.ia[i+1];j++) {
				jj = _ar.ja[j];
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

// Fill elements of the collapsed matrix

// Init elements by zeroes

	bla_2 = _bla*_bla;

	nz = temp.ia[nsupa] * bla_2;

	for (i=0;i<nz;i++) temp.a[i] = 0.0e0;

// Init diagonal elements by all ones

	for (isup=0;isup<nsupa;isup++) {
		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			if (jsup == isup) {
				for (int kii=0;kii<_bla;kii++) {
					temp.a[j*bla_2+kii*_bla+kii] = 1.0e0;
				};
			};
		};
	};

// Fill all remaining elements

	for (isup=0;isup<nsupa;isup++) {

// Init addresses of the nonzero supernodes

		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			imask[jsup] = j*bla_2;
		};

// Fill elements of the supernode row

		for (i=sp2nd[isup];i<sp2nd[isup+1] && i<nloc;i++) {
			int kii = i-sp2nd[isup];
			for (j=_ar.ia[i];j<_ar.ia[i+1];j++) {
				jj = _ar.ja[j];
				jsup = nd2sp[jj];
				int kjj = jj-sp2nd[jsup];
				int ibs = imask[jsup];
				temp.a[ibs+kjj*_bla+kii] = _ar.a[j];
			};
		};

	};

// Symmetrize diagonal supernodes

	for (isup=0;isup<nsupa;isup++) {
		for (j=temp.ia[isup];j<temp.ia[isup+1];j++) {
			jsup = temp.ja[j];
			if (jsup == isup) {
				for (int kii=0;kii<_bla;kii++) {
					for (int kjj=0;kjj<kii;kjj++) {
						temp.a[j*bla_2+kjj*_bla+kii] = temp.a[j*bla_2+kii*_bla+kjj];
					};
				};
			};
		};
	};

// Delete work arrays

	delete [] sp2nd;
	delete [] nd2sp;
	delete [] imask;
	delete [] lstloc;

// Reassign m and n

	temp.m = nsupa * _bla;
	temp.n = nsupa * _bla;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CSMatrixRB &_a) { // Output matrix

	_stream << " CSMatrixRB:" << endl;

	_stream << " Bla = " << _a.bla << endl;
	_stream << " Bla_2 = " << _a.bla_2 << endl;
	_stream << " Bla_3 = " << _a.bla_3 << endl;

	_stream << (const CSMatrixR &) _a << endl;

	return _stream;

};


